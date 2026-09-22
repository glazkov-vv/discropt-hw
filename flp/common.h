#pragma once

#include <cmath>
#include <vector>

using i64 = long long;
using ui64 = unsigned long long;

struct TFacility {
    i64 Id;
    double Cost;
    double Capacity;
    double X, Y;
};

struct TCustomer {
    i64 Id;
    double Demand;
    double X, Y;
};

inline double Dist(const TFacility& f, const TCustomer& c) {
    double dx = f.X - c.X, dy = f.Y - c.Y;
    return std::sqrt(dx * dx + dy * dy);
}
