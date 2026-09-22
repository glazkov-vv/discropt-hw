#include "sa.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>
#include <random>
#include <vector>

static constexpr double SaTimeBudget = 60.0;
static constexpr double T0Factor = 1.0;
static constexpr double TFinalFactor = 1e-5;
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

static std::vector<i64> NearestNeighborGreedy(const std::vector<TPoint>& points,
                                              std::mt19937& rng) {
  i64 n = points.size();
  std::vector<bool> visited(n, false);
  std::vector<i64> route;
  route.reserve(n);

  std::uniform_int_distribution<i64> startDist(0, n - 1);
  i64 cur = startDist(rng);
  visited[cur] = true;
  route.push_back(cur);

  for (i64 step = 1; step < n; step++) {
    double best = std::numeric_limits<double>::max();
    i64 next = -1;
    for (i64 j = 0; j < n; j++) {
      if (!visited[j]) {
        double d = Dist(points[cur], points[j]);
        if (d < best) {
          best = d;
          next = j;
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

static double TwoOptDelta(const std::vector<TPoint>& points,
                          const std::vector<i64>& route, i64 i, i64 j) {
  i64 n = route.size();
  const TPoint& a = points[route[i]];
  const TPoint& b = points[route[i + 1]];
  const TPoint& c = points[route[j]];
  const TPoint& d = points[route[(j + 1) % n]];
  return (Dist(a, c) + Dist(b, d)) - (Dist(a, b) + Dist(c, d));
}

static double OrOptDelta(const std::vector<TPoint>& points,
                         const std::vector<i64>& route, i64 i, i64 l, i64 j,
                         bool reversed) {
  i64 n = route.size();
  const TPoint& a = points[route[(i - 1 + n) % n]];
  const TPoint& sf = points[route[i]];
  const TPoint& sl = points[route[i + l - 1]];
  const TPoint& b = points[route[(i + l) % n]];
  const TPoint& c = points[route[j]];
  const TPoint& d = points[route[(j + 1) % n]];
  double removed = Dist(a, sf) + Dist(sl, b) + Dist(c, d);
  double added = Dist(a, b) + (reversed ? Dist(c, sl) + Dist(sf, d)
                                        : Dist(c, sf) + Dist(sl, d));
  return added - removed;
}

static void ApplyOrOpt(std::vector<i64>& route, i64 i, i64 l, i64 j,
                       bool reversed) {
  std::vector<i64> seg(route.begin() + i, route.begin() + i + l);
  if (reversed) std::reverse(seg.begin(), seg.end());
  route.erase(route.begin() + i, route.begin() + i + l);
  i64 insertPos = (j < i) ? (j + 1) : (j - l + 1);
  route.insert(route.begin() + insertPos, seg.begin(), seg.end());
}

std::vector<i64> SolveSA(std::vector<TPoint> points) {
  i64 n = points.size();
  if (n < 4) {
    std::vector<i64> route(n);
    for (i64 i = 0; i < n; i++) route[i] = i;
    return route;
  }

  std::mt19937 rng(Seed);

  auto route = NearestNeighborGreedy(points, rng);
  double length = RouteLength(points, route);

  double avgEdge = length / n;
  double T0 = T0Factor * avgEdge;
  double Tf = TFinalFactor * avgEdge;
  double logRatio = std::log(Tf / T0);

  std::vector<i64> bestRoute = route;
  double bestLength = length;

  std::uniform_real_distribution<double> u01(0.0, 1.0);

  auto start = std::chrono::steady_clock::now();
  while (true) {
    double elapsed =
        std::chrono::duration<double>(std::chrono::steady_clock::now() - start)
            .count();
    if (elapsed >= SaTimeBudget) {
      break;
    }
    double frac = elapsed / SaTimeBudget;
    double T = T0 * std::exp(logRatio * frac);

    for (i64 step = 0; step < n; step++) {
      double delta;
      bool isTwoOpt = (u01(rng) < 0.5);
      i64 i = 0, j = 0, l = 0;
      bool reversed = false;

      if (isTwoOpt) {
        std::uniform_int_distribution<i64> distI(0, n - 3);
        i = distI(rng);
        std::uniform_int_distribution<i64> distJ(i + 2, n - 1);
        j = distJ(rng);
        if (i == 0 && j == n - 1) {
          continue;
        }
        delta = TwoOptDelta(points, route, i, j);
      } else {
        std::uniform_int_distribution<i64> lenDist(1, 3);
        l = std::min<i64>(lenDist(rng), n - 2);
        std::uniform_int_distribution<i64> iDist(0, n - l);
        i = iDist(rng);
        std::uniform_int_distribution<i64> jDist(0, n - 1);
        j = jDist(rng);
        i64 lo = (i - 1 + n) % n;
        bool invalid = false;
        for (i64 k = 0; k <= l; k++) {
          if (j == (lo + k) % n) {
            invalid = true;
            break;
          }
        }
        if (invalid) {
          continue;
        }
        reversed = (l > 1) && (u01(rng) < 0.5);
        delta = OrOptDelta(points, route, i, l, j, reversed);
      }

      if (delta < 0 || u01(rng) < std::exp(-delta / T)) {
        if (isTwoOpt) {
          std::reverse(route.begin() + i + 1, route.begin() + j + 1);
        } else {
          ApplyOrOpt(route, i, l, j, reversed);
        }
        length += delta;
        if (length < bestLength) {
          bestLength = length;
          bestRoute = route;
        }
      }
    }
  }

  TwoOpt(points, bestRoute);
  return bestRoute;
}
