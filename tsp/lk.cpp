#include "lk.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <deque>
#include <numeric>
#include <random>
#include <vector>

static constexpr double LkTimeBudget = 60.0;
static constexpr i64 K = 10;
static constexpr i64 MaxDepth = 50;
static constexpr i64 Breadth[] = {5, 3, 1};
static constexpr i64 MaxSegmentLength = 50;
static constexpr double Eps = 1e-9;

static constexpr unsigned int Seed = 424242;

// Grid with ~2 points per cell, cells are scanned in growing rings
static std::vector<std::vector<i64>> BuildCandidateList(
    const std::vector<TPoint>& points, i64 k) {
  i64 n = points.size();
  k = std::min(k, n - 1);

  double minX = 1e18, minY = 1e18, maxX = -1e18, maxY = -1e18;
  for (const auto& p : points) {
    minX = std::min(minX, p.X);
    minY = std::min(minY, p.Y);
    maxX = std::max(maxX, p.X);
    maxY = std::max(maxY, p.Y);
  }
  double width = std::max(maxX - minX, 1e-9);
  double height = std::max(maxY - minY, 1e-9);
  double cellSize = std::sqrt(width * height * 2 / n);
  cellSize = std::max(cellSize, std::max(width, height) / n);
  i64 gridX = width / cellSize + 1;
  i64 gridY = height / cellSize + 1;

  std::vector<std::vector<i64>> cells(gridX * gridY);
  std::vector<i64> cellX(n), cellY(n);
  for (i64 i = 0; i < n; i++) {
    cellX[i] = std::min(gridX - 1, (i64)((points[i].X - minX) / cellSize));
    cellY[i] = std::min(gridY - 1, (i64)((points[i].Y - minY) / cellSize));
    cells[cellY[i] * gridX + cellX[i]].push_back(i);
  }

  std::vector<std::vector<i64>> nearest(n, std::vector<i64>(k));
  std::vector<std::pair<double, i64>> dists;
  for (i64 i = 0; i < n; i++) {
    dists.clear();
    for (i64 r = 0;; r++) {
      for (i64 cy = cellY[i] - r; cy <= cellY[i] + r; cy++) {
        if (cy < 0 || cy >= gridY) {
          continue;
        }
        bool edgeRow = cy == cellY[i] - r || cy == cellY[i] + r;
        i64 step = edgeRow ? 1 : 2 * r;
        for (i64 cx = cellX[i] - r; cx <= cellX[i] + r; cx += step) {
          if (cx < 0 || cx >= gridX) {
            continue;
          }
          for (i64 j : cells[cy * gridX + cx]) {
            if (j != i) {
              dists.push_back({Dist(points[i], points[j]), j});
            }
          }
        }
      }

      if ((i64)dists.size() >= k) {
        std::partial_sort(dists.begin(), dists.begin() + k, dists.end());
        dists.resize(k);
        if (dists[k - 1].first <= r * cellSize) {
          break;
        }
      }
      if (r > gridX && r > gridY) {
        break;
      }
    }
    for (i64 ki = 0; ki < k; ki++) {
      nearest[i][ki] = dists[ki].second;
    }
  }
  return nearest;
}

static i64 FindRoot(std::vector<i64>& parent, i64 v) {
  while (parent[v] != v) {
    parent[v] = parent[parent[v]];
    v = parent[v];
  }
  return v;
}

// Takes shortest candidate edges keeping paths, then joins path ends greedily
static std::vector<i64> GreedyEdgeTour(
    const std::vector<TPoint>& points,
    const std::vector<std::vector<i64>>& nearest) {
  i64 n = points.size();

  std::vector<std::pair<double, std::pair<i64, i64>>> edges;
  for (i64 i = 0; i < n; i++) {
    for (i64 j : nearest[i]) {
      edges.push_back({Dist(points[i], points[j]), {i, j}});
    }
  }
  std::sort(edges.begin(), edges.end());

  std::vector<std::array<i64, 2>> adj(n, {-1, -1});
  std::vector<i64> degree(n, 0);
  std::vector<i64> parent(n);
  std::iota(parent.begin(), parent.end(), 0);
  for (const auto& [d, edge] : edges) {
    auto [u, v] = edge;
    if (degree[u] == 2 || degree[v] == 2) {
      continue;
    }
    i64 ru = FindRoot(parent, u);
    i64 rv = FindRoot(parent, v);
    if (ru == rv) {
      continue;
    }
    parent[ru] = rv;
    adj[u][degree[u]++] = v;
    adj[v][degree[v]++] = u;
  }

  std::vector<i64> ends;
  for (i64 i = 0; i < n; i++) {
    if (degree[i] < 2) {
      ends.push_back(i);
    }
  }

  std::vector<bool> visited(n, false);
  std::vector<i64> route;
  route.reserve(n);

  i64 cur = ends[0];
  while (true) {
    i64 prev = -1;
    while (true) {
      visited[cur] = true;
      route.push_back(cur);
      i64 next = -1;
      for (i64 a : adj[cur]) {
        if (a != -1 && a != prev) {
          next = a;
        }
      }
      if (next == -1) {
        break;
      }
      prev = cur;
      cur = next;
    }
    if ((i64)route.size() == n) {
      break;
    }

    i64 best = -1;
    double bestDist = 1e18;
    i64 kept = 0;
    for (i64 e : ends) {
      if (visited[e]) {
        continue;
      }
      ends[kept++] = e;
      double d = Dist(points[cur], points[e]);
      if (d < bestDist) {
        bestDist = d;
        best = e;
      }
    }
    ends.resize(kept);
    cur = best;
  }
  return route;
}

struct TFlip {
  i64 From;
  i64 Length;
  bool Toggled;
};

class TLinKernighan {
 public:
  TLinKernighan(const std::vector<TPoint>& points,
                const std::vector<std::vector<i64>>& nearest,
                std::vector<i64> route)
      : Points(points),
        Nearest(nearest),
        N(route.size()),
        Route(std::move(route)),
        Pos(N),
        InQueue(N, false) {
    for (i64 i = 0; i < N; i++) {
      Pos[Route[i]] = i;
      Length += D(Route[i], Route[(i + 1) % N]);
    }
  }

  double GetLength() const { return Length; }

  const std::vector<i64>& GetRoute() const { return Route; }

  void PushAll() {
    for (i64 v : Route) {
      Push(v);
    }
  }

  void Optimize() {
    while (!Queue.empty()) {
      i64 t1 = Queue.front();
      Queue.pop_front();
      InQueue[t1] = false;

      for (bool backward : {false, true}) {
        Reversed ^= backward;
        double gain = ImproveFrom(t1);
        Reversed ^= backward;
        if (gain > Eps) {
          Length -= gain;
          Push(t1);
          for (i64 v : Touched) {
            Push(v);
          }
          break;
        }
      }
    }
  }

  // Double bridge A B C D -> A C B D on short segments near a random vertex
  void Perturb(std::mt19937& rng) {
    std::uniform_int_distribution<i64> vertexDist(0, N - 1);
    std::uniform_int_distribution<i64> lengthDist(
        1, std::min(MaxSegmentLength, (N - 2) / 2));

    i64 a = vertexDist(rng);
    i64 b1 = Next(a);
    i64 b2 = Walk(b1, lengthDist(rng) - 1);
    i64 c1 = Next(b2);
    i64 c2 = Walk(c1, lengthDist(rng) - 1);
    i64 d = Next(c2);

    Length += D(a, c1) + D(c2, b1) + D(b2, d) - D(a, b1) - D(b2, c1) -
              D(c2, d);
    Flip(b1, c2);
    Flip(c2, c1);
    Flip(b2, b1);

    for (i64 v : {a, b1, b2, c1, c2, d}) {
      Push(v);
    }
  }

  void Commit() { Log.clear(); }

  void Rollback(double length) {
    while (!Log.empty()) {
      UndoFlip();
    }
    Length = length;
  }

 private:
  double D(i64 a, i64 b) const { return Dist(Points[a], Points[b]); }

  i64 Next(i64 v) const {
    return Route[(Pos[v] + (Reversed ? N - 1 : 1)) % N];
  }

  i64 Prev(i64 v) const {
    return Route[(Pos[v] + (Reversed ? 1 : N - 1)) % N];
  }

  i64 Walk(i64 v, i64 steps) const {
    for (i64 i = 0; i < steps; i++) {
      v = Next(v);
    }
    return v;
  }

  void Push(i64 v) {
    if (!InQueue[v]) {
      InQueue[v] = true;
      Queue.push_back(v);
    }
  }

  void Reverse(i64 from, i64 length) {
    for (i64 k = 0; k < length / 2; k++) {
      i64 i = (from + k) % N;
      i64 j = (from + length - 1 - k) % N;
      std::swap(Route[i], Route[j]);
      Pos[Route[i]] = i;
      Pos[Route[j]] = j;
    }
  }

  // Reverses path a..b, or the complement if it is shorter
  void Flip(i64 a, i64 b) {
    if (Reversed) {
      std::swap(a, b);
    }
    i64 from = Pos[a];
    i64 length = (Pos[b] - from + N) % N + 1;
    bool toggled = 2 * length > N;
    if (toggled) {
      from = (Pos[b] + 1) % N;
      length = N - length;
      Reversed = !Reversed;
    }
    Reverse(from, length);
    Log.push_back({from, length, toggled});
  }

  void UndoFlip() {
    TFlip flip = Log.back();
    Log.pop_back();
    Reverse(flip.From, flip.Length);
    if (flip.Toggled) {
      Reversed = !Reversed;
    }
  }

  bool IsAdded(i64 a, i64 b) const {
    for (auto [u, v] : Added) {
      if ((u == a && v == b) || (u == b && v == a)) {
        return true;
      }
    }
    return false;
  }

  double ImproveFrom(i64 t1) {
    Touched.clear();
    Added.clear();
    BestGain = 0;
    BestLogSize = Log.size();

    i64 t2 = Next(t1);
    Step(0, t1, t2, D(t1, t2));

    while ((i64)Log.size() > BestLogSize) {
      UndoFlip();
    }
    return BestGain;
  }

  // Tour is t1 t2 ... t4 t3 ..., replace (t1,t2), (t4,t3) with (t2,t3), (t1,t4)
  void Step(i64 depth, i64 t1, i64 t2, double gain) {
    std::vector<std::pair<double, i64>> cands;
    for (i64 t3 : Nearest[t2]) {
      if (gain - D(t2, t3) <= Eps) {
        break;
      }
      if (t3 == t1 || t3 == Next(t2)) {
        continue;
      }
      i64 t4 = Prev(t3);
      if (IsAdded(t3, t4)) {
        continue;
      }
      cands.push_back({D(t3, t4) - D(t2, t3), t3});
    }
    std::sort(cands.rbegin(), cands.rend());

    i64 breadth = Breadth[std::min(depth, (i64)std::size(Breadth) - 1)];
    breadth = std::min(breadth, (i64)cands.size());

    for (i64 ci = 0; ci < breadth; ci++) {
      i64 t3 = cands[ci].second;
      i64 t4 = Prev(t3);
      double newGain = gain - D(t2, t3) + D(t3, t4);

      Flip(t2, t4);
      Added.push_back({t2, t3});
      Touched.push_back(t2);
      Touched.push_back(t3);
      Touched.push_back(t4);

      double closedGain = newGain - D(t4, t1);
      if (closedGain > BestGain + Eps) {
        BestGain = closedGain;
        BestLogSize = Log.size();
      }

      if (depth + 1 < MaxDepth) {
        Step(depth + 1, t1, t4, newGain);
      }
      if (BestGain > Eps) {
        return;
      }

      UndoFlip();
      Added.pop_back();
    }
  }

  const std::vector<TPoint>& Points;
  const std::vector<std::vector<i64>>& Nearest;
  i64 N;
  std::vector<i64> Route;
  std::vector<i64> Pos;
  bool Reversed = false;
  double Length = 0;

  std::deque<i64> Queue;
  std::vector<bool> InQueue;

  std::vector<TFlip> Log;
  std::vector<std::pair<i64, i64>> Added;
  std::vector<i64> Touched;
  double BestGain = 0;
  i64 BestLogSize = 0;
};

std::vector<i64> SolveLK(std::vector<TPoint> points) {
  i64 n = points.size();
  if (n <= 3) {
    std::vector<i64> route(n);
    std::iota(route.begin(), route.end(), 0);
    return route;
  }

  auto nearest = BuildCandidateList(points, K);
  TLinKernighan lk(points, nearest, GreedyEdgeTour(points, nearest));
  lk.PushAll();
  lk.Optimize();
  lk.Commit();
  if (n < 8) {
    return lk.GetRoute();
  }

  std::mt19937 rng(Seed);
  auto start = std::chrono::steady_clock::now();

  while (true) {
    double elapsed =
        std::chrono::duration<double>(std::chrono::steady_clock::now() - start)
            .count();
    if (elapsed >= LkTimeBudget) break;

    double before = lk.GetLength();
    lk.Perturb(rng);
    lk.Optimize();
    if (lk.GetLength() > before + Eps) {
      lk.Rollback(before);
    }
    lk.Commit();
  }

  return lk.GetRoute();
}
