#include "sol.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>
#include <numeric>
#include <random>
#include <vector>

static constexpr double TimeLimit = 55;
// Доля времени, отдаваемая второму этапу (решению TSP на каждой машине)
static constexpr double RouteShare = 0.6;
static constexpr i64 PopulationSize = 20;
static constexpr i64 NumRounds = 5;
static constexpr i64 NearestPoints = 20;
static constexpr double StartTempFactor = 1.0;
static constexpr double FinalTempRatio = 1e-3;
static constexpr double RelocateProb = 0.6;
static constexpr double TwoOptProb = 0.5;
static constexpr i64 OrOptMaxLength = 3;
static constexpr double TwoPi = 6.283185307179586;
static constexpr double Eps = 1e-9;

namespace {

using TClock = std::chrono::steady_clock;

double Seconds(TClock::time_point since) {
  return std::chrono::duration<double>(TClock::now() - since).count();
}

struct TProblem {
  const std::vector<TPoint>& Points;
  i64 Vehicles;
  double Capacity;
  // Матрица расстояний между точками
  std::vector<double> Dists;
  // Для каждого клиента — ближайшие другие клиенты (по возрастанию расстояния)
  std::vector<std::vector<i64>> NearPoints;
};

double PointDist(const TProblem& p, i64 a, i64 b) {
  return p.Dists[a * std::ssize(p.Points) + b];
}

bool Accept(double delta, double temp, std::mt19937& rng) {
  if (delta <= 0) {
    return true;
  }
  std::uniform_real_distribution<double> uni01(0.0, 1.0);
  return uni01(rng) < std::exp(-delta / temp);
}

// Распределение клиентов по машинам. Стоимость машины — жадная оценка длины ее
// маршрута, она поддерживается инкрементально.
struct TState {
  std::vector<i64> Vehicle;
  std::vector<std::vector<i64>> Members;
  std::vector<i64> PosInMembers;
  std::vector<double> Load;
  std::vector<double> Cost;
  double Total = 0;
};

TState EmptyState(i64 n, i64 v) {
  TState s;
  s.Vehicle.assign(n, -1);
  s.Members.assign(v, {});
  s.PosInMembers.assign(n, -1);
  s.Load.assign(v, 0);
  s.Cost.assign(v, 0);
  return s;
}

void Attach(const TProblem& p, TState& s, i64 c, i64 v) {
  s.PosInMembers[c] = s.Members[v].size();
  s.Members[v].push_back(c);
  s.Load[v] += p.Points[c].Demand;
  s.Vehicle[c] = v;
}

void Detach(const TProblem& p, TState& s, i64 c) {
  i64 v = s.Vehicle[c];
  auto& members = s.Members[v];
  i64 last = members.back();
  members[s.PosInMembers[c]] = last;
  s.PosInMembers[last] = s.PosInMembers[c];
  members.pop_back();
  s.PosInMembers[c] = -1;
  s.Load[v] -= p.Points[c].Demand;
  s.Vehicle[c] = -1;
}

bool Fits(const TProblem& p, const TState& s, i64 c, i64 v) {
  return s.Load[v] + p.Points[c].Demand <= p.Capacity + Eps;
}

// Жадная оценка длины маршрута: выезжаем со склада, каждый раз едем к
// ближайшему из еще не посещенных клиентов, затем возвращаемся на склад.
double GreedyCost(const TProblem& p, const std::vector<i64>& members) {
  i64 k = members.size();
  if (k == 0) {
    return 0;
  }
  static thread_local std::vector<i64> order;
  order.assign(members.begin(), members.end());
  double total = 0;
  i64 cur = 0;
  for (i64 i = 0; i < k; i++) {
    i64 best = i;
    double bestDist = PointDist(p, cur, order[i]);
    for (i64 j = i + 1; j < k; j++) {
      double d = PointDist(p, cur, order[j]);
      if (d < bestDist) {
        bestDist = d;
        best = j;
      }
    }
    std::swap(order[i], order[best]);
    total += bestDist;
    cur = order[i];
  }
  return total + PointDist(p, cur, 0);
}

// Оценка маршрута машины v, если из нее убрать клиента del и добавить клиента
// add (-1 — ничего не убирать / не добавлять)
double CostWith(const TProblem& p, const TState& s, i64 v, i64 add, i64 del) {
  static thread_local std::vector<i64> members;
  members.clear();
  for (i64 c : s.Members[v]) {
    if (c != del) {
      members.push_back(c);
    }
  }
  if (add != -1) {
    members.push_back(add);
  }
  return GreedyCost(p, members);
}

void Recost(const TProblem& p, TState& s, i64 v) {
  s.Total -= s.Cost[v];
  s.Cost[v] = GreedyCost(p, s.Members[v]);
  s.Total += s.Cost[v];
}

void RecostAll(const TProblem& p, TState& s) {
  s.Total = 0;
  for (i64 v = 0; v < p.Vehicles; v++) {
    s.Cost[v] = GreedyCost(p, s.Members[v]);
    s.Total += s.Cost[v];
  }
}

TState BuildState(const TProblem& p, const std::vector<i64>& vehicle) {
  TState s = EmptyState(p.Points.size(), p.Vehicles);
  for (i64 c = 1; c < std::ssize(vehicle); c++) {
    Attach(p, s, c, vehicle[c]);
  }
  RecostAll(p, s);
  return s;
}

// Развозка по секторам: клиенты сортируются по углу вокруг склада (начиная с
// направления startAngle) и последовательно набираются в машины.
bool Sweep(const TProblem& p, TState& s, double startAngle) {
  i64 n = p.Points.size();
  s = EmptyState(n, p.Vehicles);

  std::vector<double> angle(n, 0);
  for (i64 c = 1; c < n; c++) {
    double a = std::atan2(p.Points[c].Y - p.Points[0].Y,
                          p.Points[c].X - p.Points[0].X) -
               startAngle;
    angle[c] = a - TwoPi * std::floor(a / TwoPi);
  }
  std::vector<i64> order(n - 1);
  std::iota(order.begin(), order.end(), 1);
  std::sort(order.begin(), order.end(),
            [&](i64 a, i64 b) { return angle[a] < angle[b]; });

  i64 v = 0;
  for (i64 c : order) {
    while (v < p.Vehicles && !Fits(p, s, c, v)) {
      v++;
    }
    if (v == p.Vehicles) {
      return false;
    }
    Attach(p, s, c, v);
  }
  RecostAll(p, s);
  return true;
}

// Резервное начальное решение: клиенты в порядке убывания требований
// раскладываются в первую машину, куда влезают. Если такой нет (суммарного
// объема машин не хватает), кладем в наименее загруженную — решение будет
// недопустимым, но это видно по проверке.
void FirstFitDecreasing(const TProblem& p, TState& s) {
  i64 n = p.Points.size();
  s = EmptyState(n, p.Vehicles);

  std::vector<i64> order(n - 1);
  std::iota(order.begin(), order.end(), 1);
  std::sort(order.begin(), order.end(), [&](i64 a, i64 b) {
    return p.Points[a].Demand > p.Points[b].Demand;
  });

  for (i64 c : order) {
    i64 target = -1;
    for (i64 v = 0; v < p.Vehicles; v++) {
      if (Fits(p, s, c, v)) {
        target = v;
        break;
      }
    }
    if (target == -1) {
      target = std::min_element(s.Load.begin(), s.Load.end()) - s.Load.begin();
    }
    Attach(p, s, c, target);
  }
  RecostAll(p, s);
}

// Перемещаем клиента в машину одного из ближайших к нему клиентов.
void TryRelocate(const TProblem& p, TState& s, double temp, std::mt19937& rng) {
  std::uniform_int_distribution<i64> cDist(1, std::ssize(p.Points) - 1);
  i64 c = cDist(rng);
  const auto& near = p.NearPoints[c];
  std::uniform_int_distribution<i64> nDist(0, std::ssize(near) - 1);
  i64 from = s.Vehicle[c];
  i64 to = s.Vehicle[near[nDist(rng)]];
  if (to == from || !Fits(p, s, c, to)) {
    return;
  }
  double delta = CostWith(p, s, from, -1, c) - s.Cost[from] +
                 CostWith(p, s, to, c, -1) - s.Cost[to];
  if (!Accept(delta, temp, rng)) {
    return;
  }
  Detach(p, s, c);
  Attach(p, s, c, to);
  Recost(p, s, from);
  Recost(p, s, to);
}

// Меняем местами двух близких клиентов из разных машин.
void TrySwap(const TProblem& p, TState& s, double temp, std::mt19937& rng) {
  std::uniform_int_distribution<i64> cDist(1, std::ssize(p.Points) - 1);
  i64 c1 = cDist(rng);
  const auto& near = p.NearPoints[c1];
  std::uniform_int_distribution<i64> nDist(0, std::ssize(near) - 1);
  i64 c2 = near[nDist(rng)];
  i64 v1 = s.Vehicle[c1], v2 = s.Vehicle[c2];
  if (v1 == v2) {
    return;
  }
  double d1 = p.Points[c1].Demand, d2 = p.Points[c2].Demand;
  if (s.Load[v1] - d1 + d2 > p.Capacity + Eps ||
      s.Load[v2] - d2 + d1 > p.Capacity + Eps) {
    return;
  }
  double delta = CostWith(p, s, v1, c2, c1) - s.Cost[v1] +
                 CostWith(p, s, v2, c1, c2) - s.Cost[v2];
  if (!Accept(delta, temp, rng)) {
    return;
  }
  Detach(p, s, c1);
  Detach(p, s, c2);
  Attach(p, s, c1, v2);
  Attach(p, s, c2, v1);
  Recost(p, s, v1);
  Recost(p, s, v2);
}

// Начальная температура — средний положительный прирост от перемещения
// случайного клиента в машину его случайного соседа.
double EstimateTemperature(const TProblem& p, const TState& s,
                           std::mt19937& rng) {
  std::uniform_int_distribution<i64> cDist(1, std::ssize(p.Points) - 1);
  double sum = 0;
  i64 cnt = 0;
  for (i64 i = 0; i < 1000; i++) {
    i64 c = cDist(rng);
    const auto& near = p.NearPoints[c];
    std::uniform_int_distribution<i64> nDist(0, std::ssize(near) - 1);
    i64 from = s.Vehicle[c];
    i64 to = s.Vehicle[near[nDist(rng)]];
    if (to == from) {
      continue;
    }
    double delta = CostWith(p, s, from, -1, c) - s.Cost[from] +
                   CostWith(p, s, to, c, -1) - s.Cost[to];
    if (delta > 0) {
      sum += delta;
      cnt++;
    }
  }
  return cnt == 0 ? 1.0 : sum / cnt;
}

// Обход клиентов в порядке ближайшего соседа, tour[0] — склад.
std::vector<i64> GreedyTour(const TProblem& p,
                            const std::vector<i64>& members) {
  std::vector<i64> tour;
  tour.reserve(members.size() + 1);
  tour.push_back(0);
  std::vector<i64> rest = members;
  while (!rest.empty()) {
    i64 best = 0;
    for (i64 i = 1; i < std::ssize(rest); i++) {
      if (PointDist(p, tour.back(), rest[i]) <
          PointDist(p, tour.back(), rest[best])) {
        best = i;
      }
    }
    tour.push_back(rest[best]);
    rest[best] = rest.back();
    rest.pop_back();
  }
  return tour;
}

double TourLength(const TProblem& p, const std::vector<i64>& tour) {
  i64 size = tour.size();
  double total = 0;
  for (i64 i = 0; i < size; i++) {
    total += PointDist(p, tour[i], tour[(i + 1) % size]);
  }
  return total;
}

// 2-opt: разворачиваем участок маршрута [i, j]. Возвращает изменение длины.
double TryTwoOpt(const TProblem& p, std::vector<i64>& tour, double temp,
                 std::mt19937& rng) {
  i64 size = tour.size();
  std::uniform_int_distribution<i64> posDist(1, size - 1);
  i64 i = posDist(rng), j = posDist(rng);
  if (i > j) {
    std::swap(i, j);
  }
  if (i == j || (i == 1 && j == size - 1)) {
    return 0;
  }
  i64 prev = tour[i - 1], next = tour[(j + 1) % size];
  double delta = PointDist(p, prev, tour[j]) + PointDist(p, tour[i], next) -
                 PointDist(p, prev, tour[i]) - PointDist(p, tour[j], next);
  if (!Accept(delta, temp, rng)) {
    return 0;
  }
  std::reverse(tour.begin() + i, tour.begin() + j + 1);
  return delta;
}

// Or-opt: переносим отрезок длины len (возможно развернув его) после позиции j.
// Возвращает изменение длины.
double TryOrOpt(const TProblem& p, std::vector<i64>& tour, double temp,
                std::mt19937& rng) {
  i64 size = tour.size();
  std::uniform_int_distribution<i64> lenDist(1, OrOptMaxLength);
  i64 len = std::min(lenDist(rng), size - 2);
  std::uniform_int_distribution<i64> startDist(1, size - len);
  i64 i = startDist(rng);
  std::uniform_int_distribution<i64> posDist(0, size - 1);
  i64 j = posDist(rng);
  if (j >= i - 1 && j <= i + len - 1) {
    return 0;
  }

  i64 head = tour[i], tail = tour[i + len - 1];
  i64 prev = tour[i - 1], next = tour[(i + len) % size];
  i64 left = tour[j], right = tour[(j + 1) % size];
  std::uniform_real_distribution<double> uni01(0.0, 1.0);
  bool reversed = uni01(rng) < 0.5;
  double delta = PointDist(p, prev, next) - PointDist(p, prev, head) -
                 PointDist(p, tail, next) - PointDist(p, left, right);
  if (reversed) {
    delta += PointDist(p, left, tail) + PointDist(p, head, right);
  } else {
    delta += PointDist(p, left, head) + PointDist(p, tail, right);
  }
  if (!Accept(delta, temp, rng)) {
    return 0;
  }

  i64 at;
  if (j > i) {
    std::rotate(tour.begin() + i, tour.begin() + i + len, tour.begin() + j + 1);
    at = j + 1 - len;
  } else {
    std::rotate(tour.begin() + j + 1, tour.begin() + i, tour.begin() + i + len);
    at = j + 1;
  }
  if (reversed) {
    std::reverse(tour.begin() + at, tour.begin() + at + len);
  }
  return delta;
}




std::vector<i64> AnnealRoute(const TProblem& p, const std::vector<i64>& members,
                             TClock::time_point start, double from, double to,
                             std::mt19937& rng) {
  std::vector<i64> tour = GreedyTour(p, members);
  if (std::ssize(tour) <= 3) {
    return {tour.begin() + 1, tour.end()};
  }

  double cur = TourLength(p, tour);
  double t0 = cur / std::ssize(tour) * StartTempFactor;
  double t1 = t0 * FinalTempRatio;
  double temp = t0;

  std::vector<i64> best = tour;
  double bestLength = cur;
  std::uniform_real_distribution<double> uni01(0.0, 1.0);
  for (i64 iter = 0;; iter++) {
    if ((iter & 255) == 0) {
      double frac = (Seconds(start) - from) / (to - from);
      if (frac >= 1) {
        break;
      }
      temp = t0 * std::pow(t1 / t0, frac);
    }
    if (uni01(rng) < TwoOptProb) {
      cur += TryTwoOpt(p, tour, temp, rng);
    } else {
      cur += TryOrOpt(p, tour, temp, rng);
    }
    if (cur < bestLength - Eps) {
      bestLength = cur;
      best = tour;
    }
  }
  return {best.begin() + 1, best.end()};
}

}  // namespace

std::vector<std::vector<i64>> Solve(const std::vector<TPoint>& points,
                                    i64 vehicles, double capacity) {
  i64 n = points.size();
  std::mt19937 rng(42);

  auto start = TClock::now();
  auto elapsed = [&]() { return Seconds(start); };

  TProblem p{.Points = points, .Vehicles = vehicles, .Capacity = capacity};
  p.Dists.assign(n * n, 0);
  for (i64 i = 0; i < n; i++) {
    for (i64 j = 0; j < n; j++) {
      p.Dists[i * n + j] = Dist(points[i], points[j]);
    }
  }
  p.NearPoints.assign(n, {});
  for (i64 c = 1; c < n; c++) {
    std::vector<i64> ids(n - 1);
    std::iota(ids.begin(), ids.end(), 1);
    ids.erase(ids.begin() + (c - 1));
    i64 k = std::min<i64>(n - 2, NearestPoints);
    std::partial_sort(
        ids.begin(), ids.begin() + k, ids.end(),
        [&](i64 a, i64 b) { return PointDist(p, c, a) < PointDist(p, c, b); });
    ids.resize(k);
    p.NearPoints[c] = std::move(ids);
  }

  std::vector<i64> bestVehicle;
  double bestTotal = std::numeric_limits<double>::max();
  for (i64 i = 0; i < PopulationSize; i++) {
    TState s;
    if (Sweep(p, s, TwoPi * i / PopulationSize) && s.Total < bestTotal) {
      bestTotal = s.Total;
      bestVehicle = s.Vehicle;
    }
  }
  if (bestVehicle.empty()) {
    TState s;
    FirstFitDecreasing(p, s);
    bestTotal = s.Total;
    bestVehicle = s.Vehicle;
  }



  double annealLimit = TimeLimit * (1 - RouteShare);
  std::uniform_real_distribution<double> uni01(0.0, 1.0);
  for (i64 round = 0; round < NumRounds; round++) {
    double roundStart = elapsed();
    double roundDuration = (annealLimit - roundStart) / (NumRounds - round);
    if (roundDuration <= 0) {
      break;
    }

    TState s = BuildState(p, bestVehicle);
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
      if (uni01(rng) < RelocateProb) {
        TryRelocate(p, s, temp, rng);
      } else {
        TrySwap(p, s, temp, rng);
      }
      if (s.Total < bestTotal - Eps) {
        bestTotal = s.Total;
        bestVehicle = s.Vehicle;
      }
    }
  }


  TState best = BuildState(p, bestVehicle);
  i64 rest = 0;
  for (i64 v = 0; v < vehicles; v++) {
    if (!best.Members[v].empty()) {
      rest++;
    }
  }

  std::vector<std::vector<i64>> routes(vehicles);
  for (i64 v = 0; v < vehicles; v++) {
    if (best.Members[v].empty()) {
      continue;
    }
    double from = elapsed();
    double to = from + (TimeLimit - from) / rest;
    routes[v] = AnnealRoute(p, best.Members[v], start, from, to, rng);
    rest--;
  }
  return routes;
}
