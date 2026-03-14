#include "sol.h"

#include <iostream>
#include <vector>
#include <sstream>
#include <unordered_set>
#include <set>
#include <numeric>
#include <ranges>
#include <fstream>
#include <random>
#include <time.h>
#include <cmath>

static constexpr i64 budget = 10'000'000'000;

i64 CalculateDepth(i64 n) {
    for (i64 depth = 1;; depth++) {
        if (std::pow(n, depth + 1) * std::log2(n) > budget) {
            return depth - 1;
        }
    }
}

double GetObjectScore(const TObject& object)
{
    return ((double)(object.Cost)) / object.Weight;
}

std::pair<i64, std::vector<i64>> Greedy(i64 maxWeight, std::set<std::pair<double, TObject>>& objects) {
    std::vector<i64> ans;
    i64 remainingWeight = maxWeight;
    i64 score = 0;
    for (const auto& [k, v] : objects) {
        if (v.Weight <= remainingWeight) {
            ans.push_back(v.Id);
            remainingWeight -= v.Weight;
            score += v.Cost;
        }
    }
    return {score, ans};
}

std::pair<i64, std::vector<i64>> FullSearch(i64 remainDepth, i64 maxWeight, std::set<std::pair<double, TObject>>& objects) {
    if (remainDepth == 0) {
        return Greedy(maxWeight, objects);
    }
    i64 bestScore = 0;
    std::vector<i64> bestAns;

    std::vector<TObject> allObjects;
    for (const auto& [k, v] : objects) {
        allObjects.push_back(v);
    }
    for (const auto& object : allObjects) {
        if (object.Weight > maxWeight) {
            continue;
        }
        auto objectScore = GetObjectScore(object);
        

        i64 myScore = object.Cost;
        objects.erase({-objectScore, object});
        auto [otherScore, otherAns] = FullSearch(remainDepth - 1, maxWeight - object.Weight, objects);
        objects.insert({-objectScore, object});
        if (myScore + otherScore > bestScore) {
            bestScore = myScore + otherScore;
            bestAns = otherAns;
            bestAns.push_back(object.Id);
        }
    }

    return {bestScore, bestAns};
}

std::vector<i64> Solve(i64 maxWeight, std::vector<TObject> objects)
{
    std::set<std::pair<double, TObject>> objectsByScore;
    for (const auto& object : objects) {
        objectsByScore.insert({-GetObjectScore(object), object});
    }

    auto [score, ans] = FullSearch(CalculateDepth(std::ssize(objects)), maxWeight, objectsByScore); 
    return ans;
}