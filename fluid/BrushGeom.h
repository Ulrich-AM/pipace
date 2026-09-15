#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>

enum class BrushShape : uint8_t {
    Circle = 0,
    Square = 1,
    Triangle = 2
};

// Grid Y increases downward. Circle matches the existing r² disc.
inline bool brushContains(BrushShape shape, int cx, int cy, int x, int y, int radius) {
    int r = std::max(0, radius);
    int dx = x - cx;
    int dy = y - cy;
    switch (shape) {
        case BrushShape::Square:
            return std::abs(dx) <= r && std::abs(dy) <= r;
        case BrushShape::Triangle: {
            if (std::abs(dx) > r || dy < -r || dy > r) return false;
            if (r == 0) return dx == 0 && dy == 0;
            float t = static_cast<float>(dy + r) / static_cast<float>(2 * r);
            float halfW = t * static_cast<float>(r);
            return std::abs(dx) <= halfW + 0.51f;
        }
        case BrushShape::Circle:
        default:
            return dx * dx + dy * dy <= r * r;
    }
}

template <typename Fn>
inline void forEachBrushCell(BrushShape shape, int cx, int cy, int radius, Fn &&fn) {
    int r = std::max(0, radius);
    for (int y = cy - r; y <= cy + r; ++y)
        for (int x = cx - r; x <= cx + r; ++x)
            if (brushContains(shape, cx, cy, x, y, r)) fn(x, y);
}
