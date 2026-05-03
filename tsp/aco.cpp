#include "aco.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>
#include <random>
#include <vector>

static constexpr double AcoTimeBudget = 60.0;
static constexpr i64 NumAnts = 50;
static constexpr i64 K = 5;
static constexpr double Alpha = 4.0;
static constexpr double Beta = 1.0;
static constexpr double EvaporationRate = 0.1;

static constexpr unsigned int Seed = 424242;

static double RouteLength(const std::vector<TPoint>& points,
                          const std::vector<i64>& route) {
  i64 n = route.size();
  double total = 0;
  for (i64 i = 0; i < n; i++) {
    total += Dist(points[route[i]], points[route[(i + 1) % n]]);
  }
  return total;
}

static std::vector<std::vector<i64>> BuildCandidateList(
    const std::vector<TPoint>& points, i64 k) {
  i64 n = points.size();
  k = std::min(k, n - 1);
  std::vector<std::vector<i64>> nearest(n, std::vector<i64>(k));
  std::vector<std::pair<double, i64>> dists(n);
  for (i64 i = 0; i < n; i++) {
    for (i64 j = 0; j < n; j++) {
      dists[j] = {Dist(points[i], points[j]), j};
    }
    dists[i].first = 1e18;
    std::partial_sort(dists.begin(), dists.begin() + k, dists.end());
    for (i64 ki = 0; ki < k; ki++) {
      nearest[i][ki] = dists[ki].second;
    }
  }
  return nearest;
}

static std::vector<i64> AntTour(
    const std::vector<TPoint>& points,
    const std::vector<std::vector<i64>>& nearest,
    const std::vector<std::vector<double>>& currentPheromones,
    std::mt19937& rng) {
  i64 n = points.size();
  i64 k = nearest[0].size();
  std::vector<bool> visited(n, false);
  std::vector<i64> route;
  route.reserve(n);

  std::uniform_int_distribution<i64> startDist(0, n - 1);
  i64 cur = startDist(rng);
  visited[cur] = true;
  route.push_back(cur);

  std::vector<double> probs;
  std::vector<i64> cands;
  probs.reserve(k);
  cands.reserve(k);

  for (i64 step = 1; step < n; step++) {
    probs.clear();
    cands.clear();

    for (i64 ki = 0; ki < k; ki++) {
      i64 j = nearest[cur][ki];
      if (!visited[j]) {
        double d = Dist(points[cur], points[j]);
        probs.push_back(std::pow(currentPheromones[cur][ki], Alpha) *
                        std::pow(1.0 / d, Beta));
        cands.push_back(j);
      }
    }

    i64 next;
    if (!cands.empty()) {
      std::discrete_distribution<i64> dist(probs.begin(), probs.end());
      next = cands[dist(rng)];
    } else {
      double best = std::numeric_limits<double>::max();
      next = -1;
      for (i64 j = 0; j < n; j++) {
        if (!visited[j]) {
          double d = Dist(points[cur], points[j]);
          if (d < best) {
            best = d;
            next = j;
          }
        }
      }
    }

    visited[next] = true;
    route.push_back(next);
    cur = next;
  }
  return route;
}

static void TwoOpt(const std::vector<TPoint>& points, std::vector<i64>& route) {
  i64 n = route.size();
  bool improved = true;
  while (improved) {
    improved = false;
    for (i64 i = 0; i < n - 1; i++) {
      for (i64 j = i + 2; j < n; j++) {
        if (j == n - 1 && i == 0) {
          continue;
        }
        double before = Dist(points[route[i]], points[route[i + 1]]) +
                        Dist(points[route[j]], points[route[(j + 1) % n]]);
        double after = Dist(points[route[i]], points[route[j]]) +
                       Dist(points[route[i + 1]], points[route[(j + 1) % n]]);
        if (after < before) {
          std::reverse(route.begin() + i + 1, route.begin() + j + 1);
          improved = true;
        }
      }
    }
  }
}

std::vector<i64> SolveACO(std::vector<TPoint> points) {
  i64 n = points.size();
  i64 k = std::min(K, n - 1);

  auto nearest = BuildCandidateList(points, k);

  std::vector<std::vector<double>> currentPheromones(
      n, std::vector<double>(k, 1.0));

  std::mt19937 rng(Seed);
  auto start = std::chrono::steady_clock::now();

  std::vector<i64> bestRoute;
  double bestLength = std::numeric_limits<double>::max();

  while (true) {
    double elapsed =
        std::chrono::duration<double>(std::chrono::steady_clock::now() - start)
            .count();
    if (elapsed >= AcoTimeBudget) break;

    std::vector<std::vector<i64>> routes(NumAnts);
    std::vector<double> lengths(NumAnts);

    for (i64 a = 0; a < NumAnts; a++) {
      routes[a] = AntTour(points, nearest, currentPheromones, rng);
      TwoOpt(points, routes[a]);
      lengths[a] = RouteLength(points, routes[a]);
      if (lengths[a] < bestLength) {
        bestLength = lengths[a];
        bestRoute = routes[a];
      }
    }

    for (auto& row : currentPheromones) {
      for (auto& t : row) {
        t *= (1.0 - EvaporationRate);
      }
    }

    for (i64 a = 0; a < NumAnts; a++) {
      double deposit = 1 / lengths[a];
      for (i64 i = 0; i < n; i++) {
        i64 from = routes[a][i];
        i64 to = routes[a][(i + 1) % n];

        for (i64 ki = 0; ki < k; ki++) {
          if (nearest[from][ki] == to) {
            currentPheromones[from][ki] += deposit;
            break;
          }
        }
        for (i64 ki = 0; ki < k; ki++) {
          if (nearest[to][ki] == from) {
            currentPheromones[to][ki] += deposit;
            break;
          }
        }
      }
    }
  }

  TwoOpt(points, bestRoute);
  return bestRoute;
}
