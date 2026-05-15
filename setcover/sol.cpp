#include "sol.h"

#include <algorithm>
#include <chrono>
#include <numeric>
#include <random>
#include <vector>

namespace {

struct TSolution {
  std::vector<bool> selected;
  std::vector<i64> coverCount;
  i64 totalCost = 0;
};

// Жадно достраиваем покрытие. На каждом шаге с вероятностью randomPickProb
// берём случайное множество, покрывающее хотя бы один непокрытый элемент,
// иначе — лучшее по отношению (новых элементов)/стоимость.
void GreedyExtend(const std::vector<TObject>& objects, TSolution& sol,
                  i64& uncovered, std::mt19937_64& rng, double randomPickProb) {
  const i64 N = static_cast<i64>(objects.size());
  std::uniform_real_distribution<double> uni01(0.0, 1.0);
  std::vector<i64> useful;
  useful.reserve(N);

  while (uncovered > 0) {
    bool pickRandom = uni01(rng) < randomPickProb;

    i64 chosen = -1;
    if (pickRandom) {
      useful.clear();
      for (i64 i = 0; i < N; i++) {
        if (sol.selected[i]) continue;
        for (i64 s : objects[i].ss) {
          if (sol.coverCount[s] == 0) {
            useful.push_back(i);
            break;
          }
        }
      }
      if (useful.empty()) return;
      std::uniform_int_distribution<i64> di(0, useful.size() - 1);
      chosen = useful[di(rng)];
    } else {
      double bestScore = -1.0;
      for (i64 i = 0; i < N; i++) {
        if (sol.selected[i]) continue;
        i64 newCov = 0;
        for (i64 s : objects[i].ss) {
          if (sol.coverCount[s] == 0) newCov++;
        }
        if (newCov == 0) continue;
        double sc = static_cast<double>(newCov) / objects[i].c;
        if (sc > bestScore) {
          bestScore = sc;
          chosen = i;
        }
      }
      if (chosen == -1) return;
    }

    sol.selected[chosen] = true;
    sol.totalCost += objects[chosen].c;
    for (i64 s : objects[chosen].ss) {
      if (sol.coverCount[s] == 0) uncovered--;
      sol.coverCount[s]++;
    }
  }
}

// Удаляем избыточные множества: те, у которых каждый элемент уже покрыт
// другим выбранным множеством. Пробуем выбрасывать сначала самые дорогие.
void DropRedundant(const std::vector<TObject>& objects, TSolution& sol) {
  const i64 N = static_cast<i64>(objects.size());
  std::vector<i64> order;
  order.reserve(N);
  for (i64 i = 0; i < N; i++) {
    if (sol.selected[i]) order.push_back(i);
  }
  std::sort(order.begin(), order.end(),
            [&](i64 a, i64 b) { return objects[a].c > objects[b].c; });
  for (i64 idx : order) {
    bool redundant = true;
    for (i64 s : objects[idx].ss) {
      if (sol.coverCount[s] <= 1) {
        redundant = false;
        break;
      }
    }
    if (redundant) {
      sol.selected[idx] = false;
      sol.totalCost -= objects[idx].c;
      for (i64 s : objects[idx].ss) sol.coverCount[s]--;
    }
  }
}

TSolution EmptySolution(i64 m, i64 N) {
  TSolution s;
  s.selected.assign(N, false);
  s.coverCount.assign(m, 0);
  s.totalCost = 0;
  return s;
}

// Полностью жадное стартовое решение.
TSolution MakeGreedy(i64 m, const std::vector<TObject>& objects,
                     std::mt19937_64& rng) {
  TSolution s = EmptySolution(m, objects.size());
  i64 uncovered = m;
  GreedyExtend(objects, s, uncovered, rng, 0.0);
  if (uncovered == 0) DropRedundant(objects, s);
  return s;
}

// Случайное стартовое решение: берём множества в случайном порядке, добавляя
// только те, что покрывают что-то новое; затем убираем избыточные.
TSolution MakeRandom(i64 m, const std::vector<TObject>& objects,
                     std::mt19937_64& rng) {
  const i64 N = static_cast<i64>(objects.size());
  TSolution sol = EmptySolution(m, N);
  i64 uncovered = m;

  std::vector<i64> order(N);
  std::iota(order.begin(), order.end(), 0);
  std::shuffle(order.begin(), order.end(), rng);

  for (i64 idx : order) {
    if (uncovered == 0) break;
    bool useful = false;
    for (i64 e : objects[idx].ss) {
      if (sol.coverCount[e] == 0) {
        useful = true;
        break;
      }
    }
    if (!useful) continue;
    sol.selected[idx] = true;
    sol.totalCost += objects[idx].c;
    for (i64 e : objects[idx].ss) {
      if (sol.coverCount[e] == 0) uncovered--;
      sol.coverCount[e]++;
    }
  }

  if (uncovered == 0) DropRedundant(objects, sol);
  return sol;
}

}  // namespace

std::vector<i64> Solve(i64 m, std::vector<TObject> objects) {
  const i64 N = static_cast<i64>(objects.size());
  if (N == 0 || m == 0) return {};

  std::mt19937_64 rng(123456789);

  // Пул решений: 1 жадное + 4 случайных.
  constexpr i64 PoolSize = 5;
  std::vector<TSolution> pool;
  pool.reserve(PoolSize);
  pool.push_back(MakeGreedy(m, objects, rng));
  for (i64 i = 1; i < PoolSize; i++) {
    pool.push_back(MakeRandom(m, objects, rng));
  }

  // Лучшее по пулу.
  i64 bestCost = pool[0].totalCost;
  i64 bestIdx = 0;
  for (i64 i = 1; i < PoolSize; i++) {
    if (pool[i].totalCost < bestCost) {
      bestCost = pool[i].totalCost;
      bestIdx = i;
    }
  }

  const auto startTime = std::chrono::steady_clock::now();
  const auto timeLimit = std::chrono::milliseconds(60000);
  const double randomPickProb = 0.13;

  std::vector<i64> sel;
  i64 turn = 0;
  while (true) {
    auto now = std::chrono::steady_clock::now();
    if (now - startTime > timeLimit) break;

    // Перебираем слоты пула циклически.
    i64 slot = turn++ % PoolSize;
    TSolution cand = pool[slot];
    i64 uncovered = 0;
    for (i64 c : cand.coverCount) {
      if (c == 0) uncovered++;
    }

    // Убираем случайную часть выбранных множеств.
    sel.clear();
    for (i64 i = 0; i < N; i++) {
      if (cand.selected[i]) sel.push_back(i);
    }
    if (sel.empty()) continue;
    std::shuffle(sel.begin(), sel.end(), rng);
    const i64 removeCnt = sel.size() / 4;
    for (i64 t = 0; t < removeCnt; t++) {
      i64 idx = sel[t];
      cand.selected[idx] = false;
      cand.totalCost -= objects[idx].c;
      for (i64 e : objects[idx].ss) {
        cand.coverCount[e]--;
        if (cand.coverCount[e] == 0) uncovered++;
      }
    }

    // Достраиваем жадно с малой долей случайных шагов.
    GreedyExtend(objects, cand, uncovered, rng, randomPickProb);
    if (uncovered != 0) continue;
    DropRedundant(objects, cand);

    if (cand.totalCost < pool[slot].totalCost) {
      pool[slot] = std::move(cand);
      if (pool[slot].totalCost < bestCost) {
        bestCost = pool[slot].totalCost;
        bestIdx = slot;
      }
    }
  }

  std::vector<i64> ans;
  ans.reserve(N);
  for (i64 i = 0; i < N; i++) {
    if (pool[bestIdx].selected[i]) ans.push_back(objects[i].Id);
  }
  return ans;
}
