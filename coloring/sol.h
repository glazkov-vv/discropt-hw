#pragma once

#include "common.h"

#include <vector>

// Returns coloring: result[v] = color assigned to vertex v (0-indexed)
std::vector<i64> Solve(i64 n, const std::vector<TEdge>& edges);
