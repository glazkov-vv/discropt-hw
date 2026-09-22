#pragma once

#include <cmath>
#include <vector>

using i64 = long long;

// Точка 0 — склад, его требование равно нулю
struct TPoint {
    i64 Id;
    double Demand;
    double X, Y;
};

inline double Dist(const TPoint& a, const TPoint& b) {
    double dx = a.X - b.X, dy = a.Y - b.Y;
    return std::sqrt(dx * dx + dy * dy);
}
