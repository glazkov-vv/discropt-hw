#include "sol.h"

#include <algorithm>
#include <chrono>
#include <random>
#include <vector>

namespace {

// Жадно достраиваем покрытие. На каждом шаге с вероятностью randomPickProb
// берём случайное множество, покрывающее хотя бы один непокрытый элемент,
// иначе — лучшее по отношению (новых элементов)/стоимость.
void GreedyExtend(const std::vector<TObject>& objects,
                  std::vector<bool>& selected, std::vector<i64>& coverCount,
                  i64& uncovered, i64& totalCost, std::mt19937_64& rng,
                  double randomPickProb) {
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
        if (selected[i]) continue;
        for (i64 s : objects[i].ss) {
          if (coverCount[s] == 0) {
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
        if (selected[i]) continue;
        i64 newCov = 0;
        for (i64 s : objects[i].ss) {
          if (coverCount[s] == 0) newCov++;
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

    selected[chosen] = true;
    totalCost += objects[chosen].c;
    for (i64 s : objects[chosen].ss) {
      if (coverCount[s] == 0) uncovered--;
      coverCount[s]++;
    }
  }
}

// Удаляем избыточные множества: те, у которых каждый элемент уже покрыт
// другим выбранным множеством. Пробуем выбрасывать сначала самые дорогие.
void DropRedundant(const std::vector<TObject>& objects,
                   std::vector<bool>& selected, std::vector<i64>& coverCount,
                   i64& totalCost) {
  const i64 N = static_cast<i64>(objects.size());
  std::vector<i64> order;
  order.reserve(N);
  for (i64 i = 0; i < N; i++) {
    if (selected[i]) order.push_back(i);
  }
  std::sort(order.begin(), order.end(),
            [&](i64 a, i64 b) { return objects[a].c > objects[b].c; });
  for (i64 idx : order) {
    bool redundant = true;
    for (i64 s : objects[idx].ss) {
      if (coverCount[s] <= 1) {
        redundant = false;
        break;
      }
    }
    if (redundant) {
      selected[idx] = false;
      totalCost -= objects[idx].c;
      for (i64 s : objects[idx].ss) coverCount[s]--;
    }
  }
}

}  // namespace

std::vector<i64> Solve(i64 m, std::vector<TObject> objects) {
  const i64 N = static_cast<i64>(objects.size());
  if (N == 0 || m == 0) return {};

  std::vector<bool> selected(N, false);
  std::vector<i64> coverCount(m, 0);
  i64 totalCost = 0;
  i64 uncovered = m;

  std::mt19937_64 rng(123456789);

  // 1) Исходный простой жадный — как было раньше.
  GreedyExtend(objects, selected, coverCount, uncovered, totalCost, rng, 0.0);
  if (uncovered > 0) return {};
  DropRedundant(objects, selected, coverCount, totalCost);

  std::vector<bool> bestSelected = selected;
  i64 bestCost = totalCost;

  const auto startTime = std::chrono::steady_clock::now();
  const auto timeLimit = std::chrono::milliseconds(60000);

  const double randomPickProb = 0.13;

  std::vector<i64> sel;
  while (true) {
    auto now = std::chrono::steady_clock::now();
    if (now - startTime > timeLimit) break;

    // Стартуем с лучшего найденного покрытия.
    selected = bestSelected;
    std::fill(coverCount.begin(), coverCount.end(), 0);
    totalCost = 0;
    uncovered = m;
    for (i64 i = 0; i < N; i++) {
      if (selected[i]) {
        totalCost += objects[i].c;
        for (i64 s : objects[i].ss) {
          if (coverCount[s] == 0) uncovered--;
          coverCount[s]++;
        }
      }
    }

    // Убираем случайную половину выбранных множеств.
    sel.clear();
    for (i64 i = 0; i < N; i++) {
      if (selected[i]) sel.push_back(i);
    }
    if (sel.empty()) break;
    std::shuffle(sel.begin(), sel.end(), rng);
    const i64 removeCnt = sel.size() / 4;
    for (i64 t = 0; t < removeCnt; t++) {
      i64 idx = sel[t];
      selected[idx] = false;
      totalCost -= objects[idx].c;
      for (i64 s : objects[idx].ss) {
        coverCount[s]--;
        if (coverCount[s] == 0) uncovered++;
      }
    }

    // Достраиваем жадно с малой долей случайных шагов.
    GreedyExtend(objects, selected, coverCount, uncovered, totalCost, rng,
                 randomPickProb);
    if (uncovered != 0) continue;
    DropRedundant(objects, selected, coverCount, totalCost);

    if (totalCost < bestCost) {
      bestCost = totalCost;
      bestSelected = selected;
    }
  }

  std::vector<i64> ans;
  ans.reserve(N);
  for (i64 i = 0; i < N; i++) {
    if (bestSelected[i]) ans.push_back(objects[i].Id);
  }
  return ans;
}
