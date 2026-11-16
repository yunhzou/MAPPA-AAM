// One cooperative deadline shared by native scientific operations.
#pragma once
#include <algorithm>
#include <chrono>
#include <cmath>
#include <optional>
#include <stdexcept>

namespace event_native {
struct SearchExpired : std::exception {};

class Deadline {
    using Clock = std::chrono::steady_clock;
    Clock::time_point end = Clock::time_point::max();
public:
    explicit Deadline(std::optional<double> seconds = std::nullopt) {
        if (seconds) {
            if (!std::isfinite(*seconds) || *seconds < 0)
                throw std::invalid_argument("allowance must be finite and nonnegative");
            end = Clock::now() + std::chrono::duration_cast<Clock::duration>(
                std::chrono::duration<double>(*seconds));
        }
    }
    void check() const {
        if (end != Clock::time_point::max() && Clock::now() >= end) throw SearchExpired();
    }
    std::optional<double> remaining() const {
        if (end == Clock::time_point::max()) return std::nullopt;
        return std::max(0., std::chrono::duration<double>(end - Clock::now()).count());
    }
};
}
