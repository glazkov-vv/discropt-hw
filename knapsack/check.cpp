#include "check.h"

#include <unordered_set>

void Check(i64 maxWeight, const std::vector<TObject>& objects, const std::vector<i64>& ids) {
    i64 totalWeight = 0;
    for (i64 id : ids) {
        totalWeight += objects[id].Weight;
    }
    if (totalWeight > maxWeight) {
        throw "Max size exceeded";
    }
}