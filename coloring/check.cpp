#include "check.h"

void Check(i64 n, const std::vector<TEdge>& edges, const std::vector<i64>& coloring)
{
    if ((i64)coloring.size() != n) throw "Coloring size mismatch";
    for (i64 c : coloring) {
        if (c < 0) {
            throw "Uncolored vertex";
        }
    }
    for (const auto& e : edges) {
        if (coloring[e.u] == coloring[e.v]) {
            throw "Adjacent vertices share a color";
        }
    }
}
