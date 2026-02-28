#include "check.h"

#include <unordered_set>

void Check(i64 n, const std::vector<TObject>& objects, const std::vector<i64>& ids) {
    std::unordered_set<i64> allTaken;
    for (i64 id : ids) {
        for (i64 s : objects[id].ss) {
            allTaken.insert(s);
        }
    }
    if (allTaken.size() != n){
        throw "Not all objects are taken";
    }
}