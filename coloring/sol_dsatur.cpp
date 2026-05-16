#include <algorithm>
#include <chrono>
#include <iostream>
#include <numeric>
#include <queue>
#include <random>
#include <set>
#include <vector>

#include "sol.h"

static constexpr double TimeLimit = 55;

static i64 NumColors(const std::vector<i64>& coloring) {
  if (coloring.empty()) {
    return 0;
  }
  return *std::max_element(coloring.begin(), coloring.end()) + 1;
}

// DSATUR: at each step, pick the uncolored vertex whose colored neighbors use
// the largest number of distinct colors; break ties by uncolored-neighbor
// degree. Assign it the smallest color unused by its neighbors.
static std::vector<i64> Dsatur(const std::vector<std::vector<i64>>& adj) {
  i64 n = adj.size();
  std::vector<i64> color(n, -1);
  std::vector<std::set<i64>> satNbrColors(n);
  std::vector<i64> uncoloredNbrDeg(n);
  for (i64 i = 0; i < n; i++) {
    uncoloredNbrDeg[i] = adj[i].size();
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
      i64 sv = satNbrColors[v].size();
      i64 sb = satNbrColors[best].size();
      if (sv > sb || (sv == sb && uncoloredNbrDeg[v] > uncoloredNbrDeg[best])) {
        best = v;
      }
    }

    i64 c = 0;
    while (satNbrColors[best].count(c)) {
      c++;
    }
    color[best] = c;
    for (i64 u : adj[best]) {
      if (color[u] < 0) {
        satNbrColors[u].insert(c);
        uncoloredNbrDeg[u]--;
      }
    }
  }
  return color;
}

// Swap colors along the Kempe chain containing v in the bichromatic subgraph
// induced by colors {color[v], targetColor}. After the swap, v has color
// targetColor and the result is still a proper coloring.
static void KempeChainSwap(std::vector<i64>& color,
                           const std::vector<std::vector<i64>>& adj, i64 v,
                           i64 targetColor) {
  i64 oldColor = color[v];
  if (oldColor == targetColor) {
    return;
  }
  std::vector<i64> chain;
  std::vector<char> inChain(adj.size(), 0);
  std::queue<i64> q;
  q.push(v);
  inChain[v] = 1;
  chain.push_back(v);
  while (!q.empty()) {
    i64 u = q.front();
    q.pop();
    i64 look = (color[u] == oldColor) ? targetColor : oldColor;
    for (i64 w : adj[u]) {
      if (!inChain[w] && color[w] == look) {
        inChain[w] = 1;
        chain.push_back(w);
        q.push(w);
      }
    }
  }
  for (i64 u : chain) {
    color[u] = (color[u] == oldColor) ? targetColor : oldColor;
  }
}

// Try to empty the color class `targetColor` by relocating each of its
// vertices. First attempt a direct recolor; otherwise try Kempe chain swaps
// that shrink the targetColor class. Returns true iff the class is fully
// emptied.
static bool TryEliminateColor(std::vector<i64>& color,
                              const std::vector<std::vector<i64>>& adj,
                              i64 targetColor, std::mt19937& rng) {
  i64 n = color.size();
  auto backup = color;
  i64 maxColor = *std::max_element(color.begin(), color.end());

  auto classSize = [&](i64 c) {
    i64 s = 0;
    for (i64 x : color) {
      if (x == c) s++;
    }
    return s;
  };

  i64 stall = 0;
  while (true) {
    std::vector<i64> members;
    for (i64 v = 0; v < n; v++) {
      if (color[v] == targetColor) {
        members.push_back(v);
      }
    }
    if (members.empty()) {
      return true;
    }
    std::shuffle(members.begin(), members.end(), rng);

    bool progress = false;
    for (i64 v : members) {
      std::vector<char> nbrColor(maxColor + 1, 0);
      for (i64 u : adj[v]) {
        if (color[u] >= 0 && color[u] <= maxColor) {
          nbrColor[color[u]] = 1;
        }
      }
      i64 directColor = -1;
      for (i64 c = 0; c <= maxColor; c++) {
        if (c != targetColor && !nbrColor[c]) {
          directColor = c;
          break;
        }
      }
      if (directColor != -1) {
        color[v] = directColor;
        progress = true;
      }
    }
    if (progress) {
      stall = 0;
      continue;
    }

    // No direct move possible. Try Kempe chain swaps that strictly shrink
    // the target class.
    i64 curSize = std::ssize(members);
    i64 bestV = -1, bestC = -1, bestSize = curSize;
    for (i64 v : members) {
      for (i64 c = 0; c <= maxColor; c++) {
        if (c == targetColor) continue;
        auto trial = color;
        KempeChainSwap(trial, adj, v, c);
        i64 s = 0;
        for (i64 x : trial) {
          if (x == targetColor) s++;
        }
        if (s < bestSize) {
          bestSize = s;
          bestV = v;
          bestC = c;
        }
      }
    }
    if (bestV != -1) {
      KempeChainSwap(color, adj, bestV, bestC);
      stall = 0;
      continue;
    }

    // Random Kempe swap to perturb the configuration, hoping a direct move
    // opens up next iteration.
    if (stall < 4) {
      std::uniform_int_distribution<i64> mDist(0, std::ssize(members) - 1);
      std::uniform_int_distribution<i64> cDist(0, maxColor);
      i64 v = members[mDist(rng)];
      i64 c;
      do {
        c = cDist(rng);
      } while (c == targetColor);
      KempeChainSwap(color, adj, v, c);
      stall++;
      continue;
    }

    color = std::move(backup);
    return false;
  }
}

// After Kempe-chain shuffling, re-derive a vertex ordering from the current
// coloring (smaller color classes first) and re-run a greedy coloring.
static std::vector<i64> GreedyByOrdering(
    const std::vector<std::vector<i64>>& adj, const std::vector<i64>& order) {
  i64 n = adj.size();
  std::vector<i64> color(n, -1);
  for (i64 v : order) {
    i64 deg = std::ssize(adj[v]);
    std::vector<char> forbidden(deg + 2, 0);
    for (i64 u : adj[v]) {
      if (color[u] >= 0 && color[u] <= deg) {
        forbidden[color[u]] = 1;
      }
    }
    for (i64 c = 0;; c++) {
      if (!forbidden[c]) {
        color[v] = c;
        break;
      }
    }
  }
  return color;
}

static std::vector<i64> OrderingFromColoring(const std::vector<i64>& color,
                                             std::mt19937& rng) {
  i64 n = color.size();
  i64 k = NumColors(color);
  std::vector<std::vector<i64>> classes(k);
  for (i64 v = 0; v < n; v++) {
    classes[color[v]].push_back(v);
  }
  // Largest class first — the standard "largest-first" trick after a target
  // class has been eliminated. The smallest classes (often singletons from
  // failed Kempe attempts) are pushed to the end.
  std::stable_sort(
      classes.begin(), classes.end(),
      [](const auto& a, const auto& b) { return a.size() > b.size(); });
  std::vector<i64> order;
  order.reserve(n);
  for (auto& cls : classes) {
    std::shuffle(cls.begin(), cls.end(), rng);
    for (i64 v : cls) {
      order.push_back(v);
    }
  }
  return order;
}

std::vector<i64> Solve(i64 n, const std::vector<TEdge>& edges) {
  std::vector<std::vector<i64>> adj(n);
  for (const auto& e : edges) {
    adj[e.u].push_back(e.v);
    adj[e.v].push_back(e.u);
  }

  std::mt19937 rng(42);

  std::vector<i64> best = Dsatur(adj);
  i64 bestK = NumColors(best);

  auto startTime = std::chrono::steady_clock::now();
  auto elapsed = [&]() {
    return std::chrono::duration<double>(std::chrono::steady_clock::now() -
                                         startTime)
        .count();
  };

  std::vector<i64> current = best;

  while (elapsed() < TimeLimit) {
    std::cout << "iteration " << elapsed() << std::endl;
    i64 k = NumColors(current);

    // Identify the smallest color class and try to dissolve it via Kempe
    // chains.
    std::vector<i64> sizes(k, 0);
    for (i64 c : current) {
      sizes[c]++;
    }
    std::vector<i64> classOrder(k);
    std::iota(classOrder.begin(), classOrder.end(), 0);
    std::stable_sort(classOrder.begin(), classOrder.end(),
                     [&](i64 a, i64 b) { return sizes[a] < sizes[b]; });

    bool improved = false;
    for (i64 target : classOrder) {
      auto trial = current;
      if (TryEliminateColor(trial, adj, target, rng)) {
        // Compact colors so they stay 0..k-2.
        std::vector<i64> remap(k, -1);
        i64 next = 0;
        for (i64 v = 0; v < n; v++) {
          if (remap[trial[v]] == -1) {
            remap[trial[v]] = next++;
          }
        }
        for (i64 v = 0; v < n; v++) {
          trial[v] = remap[trial[v]];
        }
        current = std::move(trial);
        improved = true;
        if (NumColors(current) < bestK) {
          bestK = NumColors(current);
          best = current;
        }
        break;
      }
      if (elapsed() >= TimeLimit) {
        break;
      }
    }

    if (!improved) {
      // Stuck: perturb by greedy re-coloring with a randomized class-based
      // ordering, then restart improvement from there if it isn't worse.
      auto order = OrderingFromColoring(current, rng);
      auto reGreedy = GreedyByOrdering(adj, order);
      i64 rk = NumColors(reGreedy);
      if (rk < bestK) {
        bestK = rk;
        best = reGreedy;
      }
      if (rk <= NumColors(current)) {
        current = std::move(reGreedy);
      } else {
        // Random restart from a shuffled DSATUR-like baseline.
        current = best;
        std::vector<i64> shuf(n);
        std::iota(shuf.begin(), shuf.end(), 0);
        std::shuffle(shuf.begin(), shuf.end(), rng);
        for (i64 i = 0; i + 1 < n; i += 2) {
          std::uniform_int_distribution<i64> cd(0, NumColors(current) - 1);
          i64 v = shuf[i];
          i64 c = cd(rng);
          KempeChainSwap(current, adj, v, c);
        }
      }
    }
  }

  return best;
}
