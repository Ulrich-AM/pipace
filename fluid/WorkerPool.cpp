#include "WorkerPool.h"

void WorkerPool::shutdown() {
    {
        std::lock_guard<std::mutex> lock(mutex_);
        stop_ = true;
        ++generation_;
    }
    cvStart_.notify_all();
    for (auto &t : workers_) {
        if (t.joinable()) t.join();
    }
    workers_.clear();
    stop_ = false;
    jobFn_ = nullptr;
    jobCtx_ = nullptr;
}

void WorkerPool::setExtraWorkers(int extra) {
    extra = std::max(0, extra);
    if (static_cast<int>(workers_.size()) == extra) return;
    shutdown();
    if (extra == 0) return;
    stop_ = false;
    workers_.reserve(static_cast<size_t>(extra));
    for (int i = 0; i < extra; ++i)
        workers_.emplace_back([this] { workerLoop(); });
}

void WorkerPool::workerLoop() {
    unsigned seen = 0;
    for (;;) {
        {
            std::unique_lock<std::mutex> lock(mutex_);
            cvStart_.wait(lock, [&] { return stop_ || generation_ != seen; });
            if (stop_) return;
            seen = generation_;
        }
        JobFn fn = jobFn_;
        void *ctx = jobCtx_;
        for (;;) {
            int a = cursor_.fetch_add(grain_, std::memory_order_relaxed);
            if (a >= rangeEnd_) break;
            int b = std::min(rangeEnd_, a + grain_);
            if (fn) fn(ctx, a, b);
        }
        if (busy_.fetch_sub(1, std::memory_order_acq_rel) == 1) {
            std::lock_guard<std::mutex> lock(mutex_);
            cvDone_.notify_one();
        }
    }
}
