#include "check.h"

void Check(const std::vector<TPoint>& points, i64 vehicles, double capacity, const std::vector<std::vector<i64>>& routes)
{
    i64 n = points.size();
    if ((i64)routes.size() != vehicles) throw "Routes count mismatch";
    std::vector<bool> visited(n, false);
    i64 count = 0;
    for (const auto& route : routes) {
        double load = 0;
        for (i64 id : route) {
            if (id <= 0 || id >= n) {
                throw "Invalid point id";
            }
            if (visited[id]) {
                throw "Point visited twice";
            }
            visited[id] = true;
            count++;
            load += points[id].Demand;
        }
        if (load > capacity + 1e-6) {
            throw "Vehicle capacity exceeded";
        }
    }
    if (count != n - 1) throw "Not all customers are visited";
}
