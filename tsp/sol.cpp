#include "sol.h"

#include <chrono>
#include <cmath>
#include <limits>
#include <random>

static constexpr double TimeBudget = 60.0;

static double RouteLength(const std::vector<TPoint>& points,
                          const std::vector<i64>& route) {
  i64 n = route.size();
  double total = 0;
  for (i64 i = 0; i < n; i++) {
    total += Dist(points[route[i]], points[route[(i + 1) % n]]);
  }
  return total;
}

static std::vector<i64> RandomizedGreedy(const std::vector<TPoint>& points,
                                         std::mt19937& rng) {
  i64 n = points.size();
  std::vector<bool> visited(n, false);
  std::vector<i64> route;

  std::uniform_int_distribution<i64> startDist(0, n - 1);
  i64 cur = startDist(rng);
  visited[cur] = true;
  route.push_back(cur);

  static constexpr double topProbs[] = {0.995, 0.005, 0.00};

  for (i64 step = 1; step < n; step++) {
    std::vector<std::pair<double, i64>> candidates;
    for (i64 i = 0; i < n; i++) {
      if (!visited[i]) {
        candidates.push_back({Dist(points[cur], points[i]), i});
      }
    }
    i64 k = std::min((i64)std::size(topProbs), (i64)candidates.size());
    std::partial_sort(candidates.begin(), candidates.begin() + k,
                      candidates.end());

    std::discrete_distribution<i64> dist(topProbs, topProbs + k);
    i64 next = candidates[dist(rng)].second;
    visited[next] = true;
    route.push_back(next);
    cur = next;
  }
  return route;
}

std::vector<i64> Solve(std::vector<TPoint> points) {
  std::mt19937 rng(std::random_device{}());
  auto start = std::chrono::steady_clock::now();

  std::vector<i64> bestRoute;
  double bestLength = std::numeric_limits<double>::max();

  while (true) {
    double elapsed =
        std::chrono::duration<double>(std::chrono::steady_clock::now() - start)
            .count();
    if (elapsed >= TimeBudget) break;

    auto route = RandomizedGreedy(points, rng);
    double length = RouteLength(points, route);
    if (length < bestLength) {
      bestLength = length;
      bestRoute = std::move(route);
    }
  }

  return bestRoute;
}
