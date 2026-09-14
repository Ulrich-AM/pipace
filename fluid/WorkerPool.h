#pragma once

#include <algorithm>
#include <atomic>
#include <condition_variable>
#include <mutex>
#include <thread>
#include <vector>

// Persistent worker pool. Jobs are coarse index ranges, not per-cell tasks.
// The calling thread always participates. extraWorkers()==0 is a direct
// sequential call with no pool involvement.
class WorkerPool {
public:
    WorkerPool() = default;
    ~WorkerPool() { shutdown(); }

    WorkerPool(WorkerPool const &) = delete;
    WorkerPool &operator=(WorkerPool const &) = delete;

    // Extra worker threads beyond the calling thread. 0 => fully serial.
    void setExtraWorkers(int extra);
    int extraWorkers() const { return static_cast<int>(workers_.size()); }
    int totalWorkers() const { return extraWorkers() + 1; }

    // Parallel for over [begin, end). Direct sequential path when the pool is
    // empty or the span is below serialMin_ (no std::function allocation).
    template <typename Fn>
    void parallelFor(int begin, int end, Fn const &fn);

private:
    void shutdown();
    void workerLoop();

    using JobFn = void (*)(void *, int, int);

    std::vector<std::thread> workers_;
    std::mutex mutex_;
    std::condition_variable cvStart_;
    std::condition_variable cvDone_;
    JobFn jobFn_ = nullptr;
    void *jobCtx_ = nullptr;
    std::atomic<int> cursor_{0};
    int rangeEnd_ = 0;
    int grain_ = 64;
    int serialMin_ = 256;
    std::atomic<int> busy_{0};
    bool stop_ = false;
    unsigned generation_ = 0;
};

template <typename Fn>
void WorkerPool::parallelFor(int begin, int end, Fn const &fn) {
    if (end <= begin) return;
    int span = end - begin;
    if (workers_.empty() || span < serialMin_) {
        fn(begin, end);
        return;
    }

    int n = static_cast<int>(workers_.size()) + 1;
    grain_ = std::max(serialMin_, (span + n - 1) / n);
    rangeEnd_ = end;
    cursor_.store(begin, std::memory_order_relaxed);
    jobFn_ = [](void *p, int a, int b) { (*static_cast<Fn *>(p))(a, b); };
    jobCtx_ = const_cast<void *>(static_cast<void const *>(&fn));
    busy_.store(static_cast<int>(workers_.size()) + 1, std::memory_order_relaxed);
    {
        std::lock_guard<std::mutex> lock(mutex_);
        ++generation_;
    }
    cvStart_.notify_all();

    for (;;) {
        int a = cursor_.fetch_add(grain_, std::memory_order_relaxed);
        if (a >= end) break;
        int b = std::min(end, a + grain_);
        fn(a, b);
    }

    if (busy_.fetch_sub(1, std::memory_order_acq_rel) != 1) {
        std::unique_lock<std::mutex> lock(mutex_);
        cvDone_.wait(lock, [&] { return busy_.load(std::memory_order_acquire) == 0; });
    }
    jobFn_ = nullptr;
    jobCtx_ = nullptr;
}
