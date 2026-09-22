#include <algorithm>
#include <chrono>
#include <random>
#include <set>
#include <vector>

#include "sol.h"

static constexpr double TimeLimit = 55;

i64 NumColors(const std::vector<i64>& coloring) {
  if (coloring.empty()) {
    return 0;
  }
  return *std::max_element(coloring.begin(), coloring.end()) + 1;
}

std::vector<i64> Dsatur(const std::vector<std::vector<i64>>& adj) {
  i64 n = adj.size();
  std::vector<i64> color(n, -1);
  std::vector<std::set<i64>> satNumColors(n);
  std::vector<i64> uncoloredNumDeg(n);
  for (i64 i = 0; i < n; i++) {
    uncoloredNumDeg[i] = adj[i].size();
  }

  for (i64 step = 0; step < n; step++) {
    i64 best = -1;
    for (i64 v = 0; v < n; v++) {
      if (color[v] >= 0) {
        continue;
      }
      if (best == -1) {
        best = v;
        continue;
      }
      i64 sv = satNumColors[v].size();
      i64 sb = satNumColors[best].size();
      if (sv > sb || (sv == sb && uncoloredNumDeg[v] > uncoloredNumDeg[best])) {
        best = v;
      }
    }

    i64 c = 0;
    while (satNumColors[best].count(c)) {
      c++;
    }
    color[best] = c;
    for (i64 u : adj[best]) {
      if (color[u] < 0) {
        satNumColors[u].insert(c);
        uncoloredNumDeg[u]--;
      }
    }
  }
  return color;
}

bool Tabucol(std::vector<i64>& color,
                    const std::vector<std::vector<i64>>& adj, i64 k,
                    std::mt19937& rng, const auto& elapsed) {
  i64 n = adj.size();
  std::vector<std::vector<i64>> numAdjColored(n, std::vector<i64>(k, 0));
  std::vector<std::vector<i64>> tabuUntil(n, std::vector<i64>(k, 0));
  i64 conflicts = 0;
  for (i64 v = 0; v < n; v++) {
    for (i64 u : adj[v]) {
      numAdjColored[v][color[u]]++;
    }
    conflicts += numAdjColored[v][color[v]];
  }
  conflicts /= 2;
  i64 bestConflicts = conflicts;
  std::uniform_int_distribution<i64> tenureDist(0, 9);

  for (i64 iter = 0; conflicts > 0; iter++) {
    if (elapsed() >= TimeLimit) {
      return false;
    }

    i64 bestV = -1, bestC = -1, bestDelta = 0, numTies = 0;
    for (i64 v = 0; v < n; v++) {
      i64 cur = numAdjColored[v][color[v]];
      if (cur == 0) {
        continue;
      }
      for (i64 c = 0; c < k; c++) {
        if (c == color[v]) {
          continue;
        }
        i64 delta = numAdjColored[v][c] - cur;
        if (tabuUntil[v][c] > iter && conflicts + delta >= bestConflicts) {
          continue;
        }
        if (bestV == -1 || delta < bestDelta) {
          bestV = v;
          bestC = c;
          bestDelta = delta;
          numTies = 1;
        } else if (delta == bestDelta && rng() % ++numTies == 0) {
          bestV = v;
          bestC = c;
        }
      }
    }
    if (bestV == -1) {
      continue;
    }

    i64 oldColor = color[bestV];
    color[bestV] = bestC;
    for (i64 u : adj[bestV]) {
      numAdjColored[u][oldColor]--;
      numAdjColored[u][bestC]++;
    }
    conflicts += bestDelta;
    bestConflicts = std::min(bestConflicts, conflicts);
    tabuUntil[bestV][oldColor] =
        iter + static_cast<i64>(0.6 * conflicts) + tenureDist(rng);
  }
  return true;
}

std::vector<i64> Solve(i64 n, const std::vector<TEdge>& edges) {
  std::vector<std::vector<i64>> adj(n);
  for (const auto& e : edges) {
    adj[e.u].push_back(e.v);
    adj[e.v].push_back(e.u);
  }

  std::mt19937 rng(42);

  auto startTime = std::chrono::steady_clock::now();
  auto elapsed = [&]() {
    return std::chrono::duration<double>(std::chrono::steady_clock::now() -
                                         startTime)
        .count();
  };

  std::vector<i64> best = Dsatur(adj);

  while (elapsed() < TimeLimit) {
    i64 k = NumColors(best) - 1;
    auto color = best;
    std::uniform_int_distribution<i64> colorDist(0, k - 1);
    for (i64& c : color) {
      if (c >= k) {
        c = colorDist(rng);
      }
    }
    if (!Tabucol(color, adj, k, rng, elapsed)) {
      break;
    }
    best = std::move(color);
  }

  return best;
}
