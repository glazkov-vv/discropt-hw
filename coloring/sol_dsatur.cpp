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

static std::vector<i64> Dsatur(const std::vector<std::vector<i64>>& adj) {
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
    i64 look = ((color[u] == oldColor) ? targetColor : oldColor);
    for (i64 w : adj[u]) {
      if (!inChain[w] && (color[w] == look)) {
        inChain[w] = 1;
        chain.push_back(w);
        q.push(w);
      }
    }
  }
  for (i64 u : chain) {
    color[u] = ((color[u] == oldColor) ? targetColor : oldColor);
  }
}

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

  // std::cout << "Tryeliminatecolor " << std::endl;

  constexpr double InnerBudget = 1.0;
  auto innerStart = std::chrono::steady_clock::now();
  auto innerElapsed = [&] {
    return std::chrono::duration<double>(std::chrono::steady_clock::now() -
                                         innerStart)
        .count();
  };

  while (true) {
    if (innerElapsed() >= InnerBudget) {
      color = std::move(backup);
      return false;
    }

    std::vector<i64> members;
    for (i64 v = 0; v < n; v++) {
      if (color[v] == targetColor) {
        members.push_back(v);
      }
    }
    // std::cout << "Current smallest size " << members.size() << std::endl;
    if (members.empty()) {
      return true;
    }
    std::shuffle(members.begin(), members.end(), rng);

    bool progress = false;
    for (i64 v : members) {
      std::vector<char> numColor(maxColor + 1, 0);
      for (i64 u : adj[v]) {
        if (color[u] >= 0 && color[u] <= maxColor) {
          numColor[color[u]] = 1;
        }
      }
      i64 directColor = -1;
      for (i64 c = 0; c <= maxColor; c++) {
        if (c != targetColor && !numColor[c]) {
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
      continue;
    }

    const i64 curSize = std::ssize(members);
    const i64 slack = std::max<i64>(2, curSize / 2);
    const i64 sizeCap = curSize + slack;

    std::uniform_int_distribution<i64> mDist(0, curSize - 1);
    std::uniform_int_distribution<i64> cDist(0, maxColor);

    bool accepted = false;
    for (int attempt = 0; attempt < 20; attempt++) {
      i64 vCand = members[mDist(rng)];
      i64 cCand = cDist(rng);
      while (cCand == targetColor) {
        cCand = cDist(rng);
      }

      auto trial = color;
      KempeChainSwap(trial, adj, vCand, cCand);
      i64 newSize = 0;
      for (i64 x : trial) {
        if (x == targetColor) {
          newSize++;
        }
      }
      if (newSize <= sizeCap) {
        color = std::move(trial);
        accepted = true;
        break;
      }
    }
    if (!accepted) {
      color = std::move(backup);
      return false;
    }
  }
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
    // std::cout << "iteration " << elapsed() << std::endl;
    i64 k = NumColors(current);

    std::vector<i64> sizes(k, 0);
    for (i64 c : current) {
      sizes[c]++;
    }
    std::vector<i64> classOrder(k);
    std::iota(classOrder.begin(), classOrder.end(), 0);
    std::stable_sort(classOrder.begin(), classOrder.end(),
                     [&](i64 a, i64 b) { return sizes[a] < sizes[b]; });

    bool improved = false;
    for (i64 i = 0; i < std::min(k, 5ll); i++) {
      i64 target = classOrder[i];
      auto trial = current;
      if (TryEliminateColor(trial, adj, target, rng)) {
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
    else {
      break;
    }
  }

  return best;
}
