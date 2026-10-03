#pragma once

#include <memory>
#include <utility>

// Solver state carried between the steps of one run (for example the tectonics
// engine from the Geosphere initializer to the evolver). It is not world state:
// a copy is empty, so snapshots and copies never share or replay solver internals,
// and the consuming step re-creates the state when it finds none.
template<typename T>
class transientstate
{
public:
    transientstate() = default;
    transientstate(const transientstate&) noexcept {}
    transientstate& operator=(const transientstate&) noexcept
    {
        value_.reset();
        return *this;
    }
    transientstate(transientstate&&) noexcept = default;
    transientstate& operator=(transientstate&&) noexcept = default;

    void set(std::shared_ptr<T> value) { value_ = std::move(value); }
    std::shared_ptr<T> take() { return std::move(value_); }
    void reset() { value_.reset(); }
    bool empty() const { return value_ == nullptr; }
    T* get() const { return value_.get(); }
    // The state, created empty when there is none.
    T& ensure()
    {
        if (value_ == nullptr)
            value_ = std::make_shared<T>();
        return *value_;
    }

private:
    std::shared_ptr<T> value_;
};
