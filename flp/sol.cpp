#include "sol.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>
#include <numeric>
#include <random>
#include <vector>

static constexpr double TimeLimit = 55;
static constexpr i64 PopulationSize = 20;
static constexpr i64 NearestFacilities = 30;
static constexpr i64 NearestCustomers = 100;
static constexpr i64 NumRounds = 5;
static constexpr double StartTempFactor = 1.0;
static constexpr double FinalTempRatio = 1e-3;
static constexpr double Eps = 1e-9;

namespace {

struct TProblem {
  const std::vector<TFacility>& Facilities;
  const std::vector<TCustomer>& Customers;
  // Для каждого покупателя — ближайшие магазины (по возрастанию расстояния)
  std::vector<std::vector<i64>> NearFacilities;
  // Для каждого магазина — ближайшие покупатели (по возрастанию расстояния)
  std::vector<std::vector<i64>> NearCustomers;
};

// Решение с поддержкой инкрементальных изменений: загрузки, списки
// покупателей по магазинам и значение целевой функции.
struct TState {
  std::vector<i64> Assign;
  std::vector<double> Load;
  std::vector<std::vector<i64>> Members;
  std::vector<i64> PosInMembers;
  double Total = 0;
};

TState EmptyState(i64 n, i64 m) {
  TState s;
  s.Assign.assign(m, -1);
  s.Load.assign(n, 0);
  s.Members.assign(n, {});
  s.PosInMembers.assign(m, -1);
  return s;
}

void Attach(const TProblem& p, TState& s, i64 c, i64 f) {
  if (s.Members[f].empty()) {
    s.Total += p.Facilities[f].Cost;
  }
  s.PosInMembers[c] = s.Members[f].size();
  s.Members[f].push_back(c);
  s.Load[f] += p.Customers[c].Demand;
  s.Total += Dist(p.Facilities[f], p.Customers[c]);
  s.Assign[c] = f;
}

void Detach(const TProblem& p, TState& s, i64 c) {
  i64 f = s.Assign[c];
  auto& members = s.Members[f];
  i64 last = members.back();
  members[s.PosInMembers[c]] = last;
  s.PosInMembers[last] = s.PosInMembers[c];
  members.pop_back();
  s.PosInMembers[c] = -1;
  s.Load[f] -= p.Customers[c].Demand;
  s.Total -= Dist(p.Facilities[f], p.Customers[c]);
  if (members.empty()) {
    s.Total -= p.Facilities[f].Cost;
  }
  s.Assign[c] = -1;
}

void MoveCustomer(const TProblem& p, TState& s, i64 c, i64 f) {
  Detach(p, s, c);
  Attach(p, s, c, f);
}

bool Fits(const TProblem& p, const TState& s, i64 c, i64 f) {
  return s.Load[f] + p.Customers[c].Demand <= p.Facilities[f].Capacity + Eps;
}

TState BuildState(const TProblem& p, const std::vector<i64>& assign) {
  TState s = EmptyState(p.Facilities.size(), p.Customers.size());
  for (i64 c = 0; c < std::ssize(assign); c++) {
    Attach(p, s, c, assign[c]);
  }
  return s;
}

// Случайно "предоткрываем" часть магазинов (их стоимость открытия не
// учитывается при выборе), после чего жадно назначаем покупателей (от больших
// требований к меньшим) туда, где прирост целевой функции прямо сейчас
// минимален.
bool RandomGreedy(const TProblem& p, TState& s, double openProb,
                  std::mt19937& rng) {
  i64 n = p.Facilities.size();
  i64 m = p.Customers.size();
  s = EmptyState(n, m);

  std::uniform_real_distribution<double> uni01(0.0, 1.0);
  std::vector<bool> preOpened(n);
  for (i64 f = 0; f < n; f++) {
    preOpened[f] = uni01(rng) < openProb;
  }

  std::vector<i64> order(m);
  std::iota(order.begin(), order.end(), 0);
  std::shuffle(order.begin(), order.end(), rng);
  std::stable_sort(order.begin(), order.end(), [&](i64 a, i64 b) {
    return p.Customers[a].Demand > p.Customers[b].Demand;
  });

  for (i64 c : order) {
    i64 bestF = -1;
    double bestCost = std::numeric_limits<double>::max();
    for (i64 f = 0; f < n; f++) {
      if (!Fits(p, s, c, f)) {
        continue;
      }
      double cost = Dist(p.Facilities[f], p.Customers[c]);
      if (s.Members[f].empty() && !preOpened[f]) {
        cost += p.Facilities[f].Cost;
      }
      if (cost < bestCost) {
        bestCost = cost;
        bestF = f;
      }
    }
    if (bestF == -1) {
      return false;
    }
    Attach(p, s, c, bestF);
  }
  return true;
}

bool Accept(double delta, double temp, std::mt19937& rng) {
  if (delta <= 0) {
    return true;
  }
  std::uniform_real_distribution<double> uni01(0.0, 1.0);
  return uni01(rng) < std::exp(-delta / temp);
}

double RelocateDelta(const TProblem& p, const TState& s, i64 c, i64 to) {
  i64 from = s.Assign[c];
  double delta = Dist(p.Facilities[to], p.Customers[c]) -
                 Dist(p.Facilities[from], p.Customers[c]);
  if (s.Members[from].size() == 1) {
    delta -= p.Facilities[from].Cost;
  }
  if (s.Members[to].empty()) {
    delta += p.Facilities[to].Cost;
  }
  return delta;
}

// Перемещаем одного покупателя в один из ближайших к нему магазинов.
void TryRelocate(const TProblem& p, TState& s, double temp,
                 std::mt19937& rng) {
  std::uniform_int_distribution<i64> cDist(0, std::ssize(p.Customers) - 1);
  i64 c = cDist(rng);
  const auto& near = p.NearFacilities[c];
  std::uniform_int_distribution<i64> fDist(0, std::ssize(near) - 1);
  i64 to = near[fDist(rng)];
  if (to == s.Assign[c] || !Fits(p, s, c, to)) {
    return;
  }
  if (Accept(RelocateDelta(p, s, c, to), temp, rng)) {
    MoveCustomer(p, s, c, to);
  }
}

// Меняем местами двух покупателей из разных (близких) магазинов.
void TrySwap(const TProblem& p, TState& s, double temp, std::mt19937& rng) {
  std::uniform_int_distribution<i64> cDist(0, std::ssize(p.Customers) - 1);
  i64 c1 = cDist(rng);
  i64 f1 = s.Assign[c1];
  const auto& near = p.NearFacilities[c1];
  std::uniform_int_distribution<i64> fDist(0, std::ssize(near) - 1);
  i64 f2 = near[fDist(rng)];
  if (f2 == f1 || s.Members[f2].empty()) {
    return;
  }
  std::uniform_int_distribution<i64> mDist(0, std::ssize(s.Members[f2]) - 1);
  i64 c2 = s.Members[f2][mDist(rng)];

  double d1 = p.Customers[c1].Demand;
  double d2 = p.Customers[c2].Demand;
  if (s.Load[f1] - d1 + d2 > p.Facilities[f1].Capacity + Eps ||
      s.Load[f2] - d2 + d1 > p.Facilities[f2].Capacity + Eps) {
    return;
  }
  double delta = Dist(p.Facilities[f2], p.Customers[c1]) +
                 Dist(p.Facilities[f1], p.Customers[c2]) -
                 Dist(p.Facilities[f1], p.Customers[c1]) -
                 Dist(p.Facilities[f2], p.Customers[c2]);
  if (Accept(delta, temp, rng)) {
    Detach(p, s, c1);
    Detach(p, s, c2);
    Attach(p, s, c1, f2);
    Attach(p, s, c2, f1);
  }
}

void Revert(const TProblem& p, TState& s,
            const std::vector<std::pair<i64, i64>>& moves) {
  for (auto it = moves.rbegin(); it != moves.rend(); ++it) {
    MoveCustomer(p, s, it->first, it->second);
  }
}

// Закрываем магазин, раздавая его покупателей в ближайшие открытые.
void TryClose(const TProblem& p, TState& s, double temp, std::mt19937& rng) {
  std::uniform_int_distribution<i64> fDist(0, std::ssize(p.Facilities) - 1);
  i64 f = -1;
  for (int attempt = 0; attempt < 50; attempt++) {
    i64 cand = fDist(rng);
    if (!s.Members[cand].empty()) {
      f = cand;
      break;
    }
  }
  if (f == -1) {
    return;
  }

  double before = s.Total;
  std::vector<i64> members = s.Members[f];
  std::sort(members.begin(), members.end(), [&](i64 a, i64 b) {
    return p.Customers[a].Demand > p.Customers[b].Demand;
  });

  std::vector<std::pair<i64, i64>> moves;
  for (i64 c : members) {
    i64 target = -1;
    for (i64 g : p.NearFacilities[c]) {
      if (g != f && !s.Members[g].empty() && Fits(p, s, c, g)) {
        target = g;
        break;
      }
    }
    if (target == -1) {
      Revert(p, s, moves);
      return;
    }
    moves.push_back({c, f});
    MoveCustomer(p, s, c, target);
  }
  if (!Accept(s.Total - before, temp, rng)) {
    Revert(p, s, moves);
  }
}

// Открываем магазин и переводим в него покупателей, которым он ближе текущего.
void TryOpen(const TProblem& p, TState& s, double temp, std::mt19937& rng) {
  std::uniform_int_distribution<i64> fDist(0, std::ssize(p.Facilities) - 1);
  i64 f = -1;
  for (int attempt = 0; attempt < 50; attempt++) {
    i64 cand = fDist(rng);
    if (s.Members[cand].empty()) {
      f = cand;
      break;
    }
  }
  if (f == -1) {
    return;
  }

  double before = s.Total;
  std::vector<std::pair<i64, i64>> moves;
  for (i64 c : p.NearCustomers[f]) {
    i64 from = s.Assign[c];
    if (Dist(p.Facilities[f], p.Customers[c]) <
            Dist(p.Facilities[from], p.Customers[c]) &&
        Fits(p, s, c, f)) {
      moves.push_back({c, from});
      MoveCustomer(p, s, c, f);
    }
  }
  if (moves.empty()) {
    return;
  }
  if (!Accept(s.Total - before, temp, rng)) {
    Revert(p, s, moves);
  }
}

// Начальная температура — средний положительный прирост от случайного
// перемещения покупателя.
double EstimateTemperature(const TProblem& p, const TState& s,
                           std::mt19937& rng) {
  std::uniform_int_distribution<i64> cDist(0, std::ssize(p.Customers) - 1);
  double sum = 0;
  i64 cnt = 0;
  for (int i = 0; i < 1000; i++) {
    i64 c = cDist(rng);
    const auto& near = p.NearFacilities[c];
    std::uniform_int_distribution<i64> fDist(0, std::ssize(near) - 1);
    i64 to = near[fDist(rng)];
    if (to == s.Assign[c]) {
      continue;
    }
    double delta = RelocateDelta(p, s, c, to);
    if (delta > 0) {
      sum += delta;
      cnt++;
    }
  }
  return cnt == 0 ? 1.0 : sum / cnt;
}

}  // namespace

std::vector<i64> Solve(const std::vector<TFacility>& facilities,
                       const std::vector<TCustomer>& customers) {
  i64 n = facilities.size();
  i64 m = customers.size();
  std::mt19937 rng(42);

  auto start = std::chrono::steady_clock::now();
  auto elapsed = [&]() {
    return std::chrono::duration<double>(std::chrono::steady_clock::now() -
                                         start)
        .count();
  };

  TProblem p{.Facilities = facilities, .Customers = customers};
  p.NearFacilities.assign(m, {});
  for (i64 c = 0; c < m; c++) {
    std::vector<i64> ids(n);
    std::iota(ids.begin(), ids.end(), 0);
    i64 k = std::min(n, NearestFacilities);
    std::partial_sort(ids.begin(), ids.begin() + k, ids.end(),
                      [&](i64 a, i64 b) {
                        return Dist(facilities[a], customers[c]) <
                               Dist(facilities[b], customers[c]);
                      });
    ids.resize(k);
    p.NearFacilities[c] = std::move(ids);
  }
  p.NearCustomers.assign(n, {});
  for (i64 f = 0; f < n; f++) {
    std::vector<i64> ids(m);
    std::iota(ids.begin(), ids.end(), 0);
    i64 k = std::min(m, NearestCustomers);
    std::partial_sort(ids.begin(), ids.begin() + k, ids.end(),
                      [&](i64 a, i64 b) {
                        return Dist(facilities[f], customers[a]) <
                               Dist(facilities[f], customers[b]);
                      });
    ids.resize(k);
    p.NearCustomers[f] = std::move(ids);
  }

  // Группа стартовых решений: случайно открытые магазины + жадное достраивание.
  std::uniform_real_distribution<double> uni01(0.0, 1.0);
  std::vector<i64> bestAssign;
  double bestTotal = std::numeric_limits<double>::max();
  for (i64 i = 0; (i < PopulationSize || bestAssign.empty()) &&
                  elapsed() < TimeLimit;
       i++) {
    TState s;
    double openProb = (i == 0) ? 0.0 : uni01(rng) * 0.5;
    if (RandomGreedy(p, s, openProb, rng) && s.Total < bestTotal) {
      bestTotal = s.Total;
      bestAssign = s.Assign;
    }
  }
  if (bestAssign.empty()) {
    return std::vector<i64>(m, 0);
  }

  // Отжиг в несколько раундов, каждый раунд стартует с лучшего найденного.
  for (i64 round = 0; round < NumRounds; round++) {
    double roundStart = elapsed();
    double roundDuration = (TimeLimit - roundStart) / (NumRounds - round);
    if (roundDuration <= 0) {
      break;
    }

    TState s = BuildState(p, bestAssign);
    double t0 = EstimateTemperature(p, s, rng) * StartTempFactor;
    double t1 = t0 * FinalTempRatio;
    double temp = t0;

    for (i64 iter = 0;; iter++) {
      if ((iter & 255) == 0) {
        double frac = (elapsed() - roundStart) / roundDuration;
        if (frac >= 1) {
          break;
        }
        temp = t0 * std::pow(t1 / t0, frac);
      }

      double r = uni01(rng);
      if (r < 0.5) {
        TryRelocate(p, s, temp, rng);
      } else if (r < 0.85) {
        TrySwap(p, s, temp, rng);
      } else if (r < 0.95) {
        TryClose(p, s, temp, rng);
      } else {
        TryOpen(p, s, temp, rng);
      }

      if (s.Total < bestTotal - Eps) {
        bestTotal = s.Total;
        bestAssign = s.Assign;
      }
    }
  }

  return bestAssign;
}
