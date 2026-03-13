#include "sol.h"

#include <iostream>
#include <vector>
#include <sstream>
#include <unordered_set>
#include <map>
#include <numeric>
#include <ranges>
#include <fstream>
#include <random>
#include <time.h>

double GetObjectScore(const TObject& object)
{
    return ((double)(object.Cost)) / object.Weight;
}

std::vector<i64> Solve(i64 maxWeight, std::vector<TObject> objects)
{
    srand(time(NULL));
    std::multimap<double, TObject> objByScore;
    for (const auto& object : objects) {
        objByScore.insert({-GetObjectScore(object), object});
    }

    std::unordered_set<i64> bestAns;
    i64 bestScore = 0;

    std::unordered_set<i64> ans;
    i64 score = 0;
    i64 remainingWeight = maxWeight;


    for (int iter = 0; iter < 100000; iter++) {
        if (iter % 1000 == 0) {
            std::cout << iter << std::endl;
        }
        if (rand() % 100 < 10) {
            if (ans.empty()) {
                continue;
            }

            i64 randomIndex = rand() % ans.size();
            auto it = std::next(ans.begin(), randomIndex);
            i64 id = *it;

            remainingWeight += objects[id].Weight;
            score -= objects[id].Cost;
            objByScore.insert({-GetObjectScore(objects[id]),objects[id]});
            ans.erase(id);
        } else {
            for (auto it = objByScore.begin(); it != objByScore.end(); it++) {
                auto obj = it->second;
                if (remainingWeight >= obj.Weight) {
                    remainingWeight -= obj.Weight;
                    score += obj.Cost;
                    objByScore.erase(it);
                    ans.insert(obj.Id);
                    break;
                }
            }
            if (score > bestScore) {
                bestScore = score;
                bestAns = ans;
                
            }
        }

    }
    

    return std::vector<i64>(bestAns.begin(), bestAns.end());
}