#ifndef parallel_rows_hpp
#define parallel_rows_hpp

#include <algorithm>
#include <condition_variable>
#include <deque>
#include <exception>
#include <functional>
#include <mutex>
#include <thread>
#include <vector>

namespace rowworkers
{
inline thread_local bool executing = false;

class Executor
{
public:
    explicit Executor(unsigned count)
    {
        try
        {
            for (unsigned i = 0; i < count; ++i)
                threads.emplace_back([this] {
                    for (;;)
                    {
                        std::function<void()> work;
                        {
                            std::unique_lock<std::mutex> lock(queueMutex);
                            ready.wait(lock, [&] { return stopping || !pendingWork.empty(); });
                            if (stopping && pendingWork.empty()) return;
                            work = std::move(pendingWork.front()); pendingWork.pop_front();
                        }
                        work();
                    }
                });
        }
        catch (...) { close(); throw; }
    }
    ~Executor() { close(); }
    void submit(std::function<void()> work)
    {
        { std::lock_guard<std::mutex> lock(queueMutex); pendingWork.push_back(std::move(work)); }
        ready.notify_one();
    }
private:
    void close()
    {
        { std::lock_guard<std::mutex> lock(queueMutex); stopping = true; }
        ready.notify_all();
        for (auto& thread : threads) if (thread.joinable()) thread.join();
    }
    std::mutex queueMutex;
    std::condition_variable ready;
    std::deque<std::function<void()>> pendingWork;
    std::vector<std::thread> threads;
    bool stopping = false;
};

inline unsigned concurrency()
{
    const unsigned available = std::thread::hardware_concurrency();
    return (std::min)(32u, available ? available : 4u);
}
inline Executor& executor()
{
    static Executor instance(concurrency() - 1);
    return instance;
}
}

// Partition an inclusive row range. Each callback owns disjoint complete rows;
// all workers finish before the caller begins the next simulation pass.
template <typename Fn>
inline void parallelforrows(int startrow, int endrow, Fn&& fn, int minrowsperworker = 32)
{
    if (endrow < startrow)
        return;

    const int totalrows = endrow - startrow + 1;
    minrowsperworker = (std::max)(1, minrowsperworker);
    unsigned int workerstouse = rowworkers::concurrency();

    if (rowworkers::executing || workerstouse <= 1 || totalrows <= minrowsperworker)
    {
        fn(startrow, endrow);
        return;
    }

    const int maxworkers = (std::max)(1, totalrows / minrowsperworker);
    workerstouse = (std::min)(workerstouse, static_cast<unsigned int>(maxworkers));

    if (workerstouse <= 1)
    {
        fn(startrow, endrow);
        return;
    }

    // Nested callbacks execute locally, so workers never wait for work queued
    // behind themselves. Independent callers may safely share the executor.
    auto& executor = rowworkers::executor();
    std::mutex completionMutex;
    std::condition_variable completed;
    unsigned remaining = workerstouse;
    std::exception_ptr failure;
    const auto run = [&](int first, int last)
    {
        rowworkers::executing = true;
        try { fn(first, last); }
        catch (...)
        {
            std::lock_guard<std::mutex> lock(completionMutex);
            if (!failure) failure = std::current_exception();
        }
        rowworkers::executing = false;
        {
            std::lock_guard<std::mutex> lock(completionMutex);
            --remaining;
            completed.notify_one();
        }
    };

    const int basechunk = totalrows / static_cast<int>(workerstouse);
    const int remainder = totalrows % static_cast<int>(workerstouse);

    int currentrow = startrow;

    for (unsigned int worker = 0; worker < workerstouse; worker++)
    {
        const int chunksize = basechunk + (worker < static_cast<unsigned int>(remainder) ? 1 : 0);
        const int chunkstart = currentrow;
        const int chunkend = currentrow + chunksize - 1;
        currentrow = chunkend + 1;

        if (worker + 1 == workerstouse)
            run(chunkstart, chunkend);
        else
        {
            try { executor.submit([&, chunkstart, chunkend] { run(chunkstart, chunkend); }); }
            catch (...) { run(chunkstart, chunkend); }
        }
    }

    std::unique_lock<std::mutex> lock(completionMutex);
    completed.wait(lock, [&] { return remaining == 0; });
    if (failure) std::rethrow_exception(failure);
}

#endif
