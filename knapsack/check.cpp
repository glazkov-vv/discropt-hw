#include "check.h"

#include <unordered_set>

void Check(i64 maxWeight, const std::vector<TObject>& objects, const std::vector<i64>& ids) {
    i64 totalWeight = 0;
    std::unordered_set<i64> distinctIds;

    for (i64 id : ids) {
        if (id < 0 || id >= objects.size()) {
            throw "Invalid id";
        }

        totalWeight += objects[id].Weight;
        distinctIds.insert(id);
    }
    if (totalWeight > maxWeight) {
        throw "Max size exceeded";
    }
    if (ids.size() != distinctIds.size()) {
        throw "Repeated ids";
    }
}