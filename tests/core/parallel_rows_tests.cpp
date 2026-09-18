#ifdef _WIN32
#include <Windows.h>
#endif
#include "parallel_rows.hpp"
#include <array>
#include <atomic>
#include <iostream>
#include <stdexcept>

int main()
{
    std::atomic<bool> correct{true};
    const auto exercise = [&] {
        for (int repetition = 0; repetition < 40; ++repetition)
        {
            std::array<int, 256> values{};
            parallelforrows(0, 63, [&](int first, int last) {
                for (int row = first; row <= last; ++row)
                    parallelforrows(0, 3, [&](int begin, int end) {
                        for (int x = begin; x <= end; ++x) ++values[row * 4 + x];
                    }, 1);
            }, 4);
            for (int value : values) if (value != 1) correct = false;
        }
    };
    std::thread other(exercise);
    exercise();
    other.join();
    std::atomic<int> finishedRows{0};
    bool caught = false;
    try
    {
        parallelforrows(0, 127, [&](int first, int last) {
            finishedRows += last - first + 1;
            if (first == 0) throw std::runtime_error("fixture");
        }, 4);
    }
    catch (const std::runtime_error&) { caught = true; }
    exercise(); // A failed callback must leave the executor reusable.
    parallelforrows(1, 0, [&](int, int) { correct = false; });
    if (!correct || !caught || finishedRows != 128)
    {
        std::cerr << "row workers must cover disjoint ranges, join before rethrow, and support nested/concurrent calls\n";
        return 1;
    }
    return 0;
}
