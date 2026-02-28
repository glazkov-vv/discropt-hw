#include "sol.h"

#include <iostream>
#include <vector>
#include <sstream>
#include <unordered_set>
#include <numeric>
#include <ranges>
#include <fstream>

double GetObjectScore(const std::unordered_set<i64>& remaining, const TObject& object)
{
    i64 coversNew = 0;
    for (const auto& s : object.ss) {
        if (remaining.contains(s)) {
            coversNew++;
        }
    }
    return static_cast<double>(coversNew) / object.c;
}

std::vector<i64> Solve(i64 m, std::vector<TObject> objects)
{
    std::unordered_set<i64> remaining = std::views::iota(0, m) |
        std::ranges::to<std::unordered_set<i64>>();

    std::vector<i64> ans;
    while (!remaining.empty()) {
        if (objects.empty()) {
            return {};
        }
        auto bestObjectIt = std::max_element(objects.begin(), objects.end(), [&] (const auto& l, const auto& r) {
            return GetObjectScore(remaining, l) < GetObjectScore(remaining, r);
        });
        
        auto& object = *bestObjectIt;
        for (const auto& s : object.ss) {
            remaining.erase(s);
        }

        ans.push_back(object.Id);
        objects.erase(bestObjectIt);
    }

    return ans;
}