#include "sol.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <numeric>
#include <random>
#include <ranges>
#include <sstream>
#include <unordered_set>
#include <vector>

namespace {

double GetObjectScore(const std::unordered_set<i64>& remaining,
                      const TObject& object) {
  i64 coversNew = 0;
  for (const auto& s : object.ss) {
    if (remaining.contains(s)) {
      coversNew++;
    }
  }
  return static_cast<double>(coversNew) / object.c;
}

// Жадная стартовая инициализация. Возвращает индексы в objects.
std::vector<i64> GreedyInitial(i64 m, const std::vector<TObject>& objects) {
  std::unordered_set<i64> remaining;
  for (i64 i = 0; i < m; i++) {
    remaining.insert(i);
  }

  std::vector<bool> used(objects.size(), false);
  std::vector<i64> picked;
  while (!remaining.empty()) {
    i64 bestIdx = -1;
    double bestScore = -1.0;
    for (i64 i = 0; i < static_cast<i64>(objects.size()); i++) {
      if (used[i]) continue;
      double s = GetObjectScore(remaining, objects[i]);
      if (s > bestScore) {
        bestScore = s;
        bestIdx = i;
      }
    }
    if (bestIdx == -1) return {};
    used[bestIdx] = true;
    for (i64 s : objects[bestIdx].ss) {
      remaining.erase(s);
    }
    picked.push_back(bestIdx);
  }
  return picked;
}

}  // namespace

std::vector<i64> Solve(i64 m, std::vector<TObject> objects) {
  const i64 N = static_cast<i64>(objects.size());
  if (N == 0 || m == 0) return {};

  auto initialIdx = GreedyInitial(m, objects);
  if (initialIdx.empty()) return {};

  // Состояние решения.
  std::vector<bool> selected(N, false);
  std::vector<i64> coverCount(m, 0);
  i64 totalCost = 0;
  i64 uncovered = m;

  for (i64 idx : initialIdx) {
    selected[idx] = true;
    totalCost += objects[idx].c;
    for (i64 s : objects[idx].ss) {
      if (coverCount[s] == 0) uncovered--;
      coverCount[s]++;
    }
  }

  // Списки для быстрой случайной выборки.
  std::vector<i64> selList, unselList;
  std::vector<i64> posInList(N, -1);
  for (i64 i = 0; i < N; i++) {
    if (selected[i]) {
      posInList[i] = selList.size();
      selList.push_back(i);
    } else {
      posInList[i] = unselList.size();
      unselList.push_back(i);
    }
  }

  auto moveSelToUnsel = [&](i64 idx) {
    i64 pos = posInList[idx];
    i64 last = selList.back();
    selList[pos] = last;
    posInList[last] = pos;
    selList.pop_back();
    posInList[idx] = unselList.size();
    unselList.push_back(idx);
  };
  auto moveUnselToSel = [&](i64 idx) {
    i64 pos = posInList[idx];
    i64 last = unselList.back();
    unselList[pos] = last;
    posInList[last] = pos;
    unselList.pop_back();
    posInList[idx] = selList.size();
    selList.push_back(idx);
  };

  // Лучшее найденное допустимое решение.
  std::vector<bool> bestSelected = selected;
  i64 bestCost = totalCost;

  // Большой штраф за непокрытый элемент: гарантированно дороже любого выбора.
  const i64 penalty = std::max<i64>(totalCost, 1) * 4 + 1000;

  std::mt19937_64 rng(123456789);
  std::uniform_real_distribution<double> uni01(0.0, 1.0);

  const auto startTime = std::chrono::steady_clock::now();
  const auto timeLimit = std::chrono::milliseconds(40000);

  // Расписание температуры: экспоненциальное по доле прошедшего времени.
  const double T0 = std::max<double>(1.0, totalCost * 0.02);
  const double Tmin = 1e-2;

  const i64 tournament = 4;  // размер турнира для предпочтительной выборки

  i64 iter = 0;
  while (true) {
    if ((iter++ & 1023) == 0) {
      auto now = std::chrono::steady_clock::now();
      if (now - startTime > timeLimit) break;
    }

    auto elapsed = std::chrono::steady_clock::now() - startTime;
    double frac = std::chrono::duration<double>(elapsed).count() /
                  std::chrono::duration<double>(timeLimit).count();
    if (frac > 1.0) frac = 1.0;
    double T = T0 * std::pow(Tmin / T0, frac);

    // Решаем: добавлять или удалять.
    bool doAdd;
    if (uncovered > 0) {
      // Если есть непокрытые — чаще добавляем.
      doAdd = uni01(rng) < 0.75;
    } else if (selList.empty()) {
      doAdd = true;
    } else {
      doAdd = uni01(rng) < 0.35;
    }
    if (doAdd && unselList.empty()) doAdd = false;
    if (!doAdd && selList.empty()) doAdd = true;
    if ((doAdd && unselList.empty()) || (!doAdd && selList.empty())) continue;

    if (doAdd) {
      // Турнир: выбираем кандидата с максимальной долей новых элементов.
      i64 bestCand = -1;
      double bestFrac = -1.0;
      i64 bestNewCov = 0;
      i64 kk = std::min<i64>(tournament, static_cast<i64>(unselList.size()));
      std::uniform_int_distribution<i64> di(0, unselList.size() - 1);
      for (i64 t = 0; t < kk; t++) {
        i64 idx = unselList[di(rng)];
        const auto& o = objects[idx];
        i64 newCov = 0;
        for (i64 s : o.ss) {
          if (coverCount[s] == 0) newCov++;
        }
        double f =
            o.ss.empty() ? 0.0 : static_cast<double>(newCov) / o.ss.size();
        if (f > bestFrac) {
          bestFrac = f;
          bestCand = idx;
          bestNewCov = newCov;
        }
      }
      if (bestCand == -1) continue;

      const auto& o = objects[bestCand];
      double delta =
          static_cast<double>(o.c) - static_cast<double>(penalty) * bestNewCov;
      bool accept = (delta <= 0.0) || (uni01(rng) < std::exp(-delta / T));
      if (!accept) continue;

      selected[bestCand] = true;
      totalCost += o.c;
      for (i64 s : o.ss) {
        if (coverCount[s] == 0) uncovered--;
        coverCount[s]++;
      }
      moveUnselToSel(bestCand);
    } else {
      // Турнир: выбираем с минимальной долей "исчезающих совсем" элементов.
      i64 bestCand = -1;
      double bestFrac = 2.0;
      i64 bestLose = 0;
      i64 kk = std::min<i64>(tournament, static_cast<i64>(selList.size()));
      std::uniform_int_distribution<i64> di(0, selList.size() - 1);
      for (i64 t = 0; t < kk; t++) {
        i64 idx = selList[di(rng)];
        const auto& o = objects[idx];
        i64 lose = 0;
        for (i64 s : o.ss) {
          if (coverCount[s] == 1) lose++;
        }
        double f = o.ss.empty() ? 0.0 : static_cast<double>(lose) / o.ss.size();
        if (f < bestFrac) {
          bestFrac = f;
          bestCand = idx;
          bestLose = lose;
        }
      }
      if (bestCand == -1) continue;

      const auto& o = objects[bestCand];
      double delta =
          -static_cast<double>(o.c) + static_cast<double>(penalty) * bestLose;
      bool accept = (delta <= 0.0) || (uni01(rng) < std::exp(-delta / T));
      if (!accept) continue;

      selected[bestCand] = false;
      totalCost -= o.c;
      for (i64 s : o.ss) {
        coverCount[s]--;
        if (coverCount[s] == 0) uncovered++;
      }
      moveSelToUnsel(bestCand);
    }

    if (uncovered == 0 && totalCost < bestCost) {
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
