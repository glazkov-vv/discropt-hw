#pragma once

#include <vector>
#include <cmath>

using i64 = long long;

struct TPoint {
    i64 Id;
    double X, Y;
};

inline double Dist(const TPoint& a, const TPoint& b) {
    double dx = a.X - b.X, dy = a.Y - b.Y;
    return std::sqrt(dx * dx + dy * dy);
}
