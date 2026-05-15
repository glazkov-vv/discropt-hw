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
// берём случайное множество, покрывающее что-то новое, иначе — лучшее по
// отношению (новых элементов)/стоимость.
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
        for (i64 e : objects[i].ss) {
          if (sol.coverCount[e] == 0) {
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
        for (i64 e : objects[i].ss) {
          if (sol.coverCount[e] == 0) newCov++;
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
    for (i64 e : objects[chosen].ss) {
      if (sol.coverCount[e] == 0) uncovered--;
      sol.coverCount[e]++;
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
    for (i64 e : objects[idx].ss) {
      if (sol.coverCount[e] <= 1) {
        redundant = false;
        break;
      }
    }
    if (redundant) {
      sol.selected[idx] = false;
      sol.totalCost -= objects[idx].c;
      for (i64 e : objects[idx].ss) sol.coverCount[e]--;
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

TSolution MakeGreedy(i64 m, const std::vector<TObject>& objects,
                     std::mt19937_64& rng) {
  TSolution s = EmptySolution(m, objects.size());
  i64 uncovered = m;
  GreedyExtend(objects, s, uncovered, rng, 0.0);
  if (uncovered == 0) DropRedundant(objects, s);
  return s;
}

// Сходство по Жаккару между двумя выбранными подмножествами множеств.
double Similarity(const TSolution& a, const TSolution& b) {
  i64 inter = 0;
  i64 uni = 0;
  const i64 N = static_cast<i64>(a.selected.size());
  for (i64 i = 0; i < N; i++) {
    bool sa = a.selected[i];
    bool sb = b.selected[i];
    if (sa && sb) inter++;
    if (sa || sb) uni++;
  }
  return uni == 0 ? 1.0 : static_cast<double>(inter) / uni;
}

// Пытаемся улучшить слот: удаляем случайную часть и достраиваем жадно с шумом.
TSolution Improve(const std::vector<TObject>& objects, const TSolution& base,
                  std::mt19937_64& rng, double randomPickProb,
                  double removeFraction) {
  const i64 N = static_cast<i64>(objects.size());
  TSolution cand = base;

  std::vector<i64> sel;
  sel.reserve(N);
  for (i64 i = 0; i < N; i++) {
    if (cand.selected[i]) sel.push_back(i);
  }
  if (sel.empty()) return cand;

  std::shuffle(sel.begin(), sel.end(), rng);
  i64 removeCnt =
      std::max<i64>(1, static_cast<i64>(sel.size() * removeFraction));
  if (removeCnt > static_cast<i64>(sel.size())) removeCnt = sel.size();

  i64 uncovered = 0;
  for (i64 c : cand.coverCount) {
    if (c == 0) uncovered++;
  }

  for (i64 t = 0; t < removeCnt; t++) {
    i64 idx = sel[t];
    cand.selected[idx] = false;
    cand.totalCost -= objects[idx].c;
    for (i64 e : objects[idx].ss) {
      cand.coverCount[e]--;
      if (cand.coverCount[e] == 0) uncovered++;
    }
  }

  GreedyExtend(objects, cand, uncovered, rng, randomPickProb);
  if (uncovered == 0) DropRedundant(objects, cand);
  return cand;
}

}  // namespace

std::vector<i64> Solve(i64 m, std::vector<TObject> objects) {
  const i64 N = static_cast<i64>(objects.size());
  if (N == 0 || m == 0) return {};

  std::mt19937_64 rng(123456789);

  constexpr i64 PoolSize = 5;
  TSolution greedy = MakeGreedy(m, objects, rng);
  if (greedy.totalCost == 0) return {};
  std::vector<TSolution> pool(PoolSize, greedy);

  const auto startTime = std::chrono::steady_clock::now();
  const auto timeLimit = std::chrono::milliseconds(60000);

  // Со временем уменьшаем "разнообразие": в начале больше шума и более крупные
  // разрушения, к концу — точная локальная доводка.
  const double randomPickProbInit = 0.25;
  const double randomPickProbFinal = 0;
  const double removeFractionInit = 0.30;
  const double removeFractionFinal = 0.10;

  const double simThreshold =
      0.85;  // пара решений с большим сходством — "слишком похожие"

  while (true) {
    auto now = std::chrono::steady_clock::now();
    if (now - startTime > timeLimit) break;

    double frac = std::chrono::duration<double>(now - startTime).count() /
                  std::chrono::duration<double>(timeLimit).count();
    if (frac > 1.0) frac = 1.0;
    const double randomPickProb =
        randomPickProbInit + (randomPickProbFinal - randomPickProbInit) * frac;
    const double removeFraction =
        removeFractionInit + (removeFractionFinal - removeFractionInit) * frac;

    // 1) Пытаемся улучшить каждый слот.
    std::vector<TSolution> candidates;
    candidates.reserve(PoolSize * 2);
    for (i64 i = 0; i < PoolSize; i++) {
      candidates.push_back(pool[i]);
    }
    for (i64 i = 0; i < PoolSize; i++) {
      TSolution imp =
          Improve(objects, pool[i], rng, randomPickProb, removeFraction);
      // принимаем только допустимое (DropRedundant вызван внутри, если ок)
      i64 unc = 0;
      for (i64 c : imp.coverCount) {
        if (c == 0) unc++;
      }
      if (unc == 0) candidates.push_back(std::move(imp));
    }

    // 2) Из всех кандидатов оставляем 5 лучших, не слишком похожих.
    std::sort(candidates.begin(), candidates.end(),
              [](const TSolution& a, const TSolution& b) {
                return a.totalCost < b.totalCost;
              });

    std::vector<TSolution> newPool;
    newPool.reserve(PoolSize);
    std::vector<int> taken(candidates.size(), 0);
    for (size_t i = 0;
         i < candidates.size() && static_cast<i64>(newPool.size()) < PoolSize;
         i++) {
      bool ok = true;
      for (const auto& r : newPool) {
        if (Similarity(candidates[i], r) >= simThreshold) {
          ok = false;
          break;
        }
      }
      if (ok) {
        newPool.push_back(candidates[i]);
        taken[i] = 1;
      }
    }
    // Если не хватает — добиваем самыми дешёвыми из оставшихся, игнорируя
    // порог.
    for (size_t i = 0;
         i < candidates.size() && static_cast<i64>(newPool.size()) < PoolSize;
         i++) {
      if (!taken[i]) newPool.push_back(candidates[i]);
    }
    pool = std::move(newPool);
  }

  // Берём лучшее по стоимости в пуле.
  i64 bestIdx = 0;
  for (i64 i = 1; i < PoolSize; i++) {
    if (pool[i].totalCost < pool[bestIdx].totalCost) bestIdx = i;
  }

  std::vector<i64> ans;
  ans.reserve(N);
  for (i64 i = 0; i < N; i++) {
    if (pool[bestIdx].selected[i]) ans.push_back(objects[i].Id);
  }
  return ans;
}
