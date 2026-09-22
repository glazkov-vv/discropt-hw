#pragma once

#include "common.h"

#include <vector>

// Returns routes: result[v] = points visited by vehicle v in order (without the depot)
std::vector<std::vector<i64>> Solve(const std::vector<TPoint>& points, i64 vehicles, double capacity);
