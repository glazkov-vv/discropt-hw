#include "check.h"

#include <unordered_set>

void Check(const std::vector<TPoint>& points, const std::vector<i64>& route) {
    i64 n = points.size();
    if ((i64)route.size() != n) {
        throw "Route length mismatch";
    }
    std::unordered_set<i64> visited;
    for (i64 id : route) {
        if (id < 0 || id >= n) {
            throw "Invalid point id";
        }
        visited.insert(id);
    }
    if ((i64)visited.size() != n) {
        throw "Route does not visit all points exactly once";
    }
}
