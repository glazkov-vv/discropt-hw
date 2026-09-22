#include "sol.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <functional>
#include <limits>
#include <numeric>
#include <random>
#include <vector>

static constexpr double TimeLimit = 55;
static constexpr i64 PopulationSize = 20;
static constexpr i64 NearestFacilities = 30;
static constexpr i64 NearestCustomers = 100;
static constexpr i64 NearestFacilityFacilities = 10;
static constexpr i64 NumRounds = 5;
static constexpr double StartTempFactor = 1.0;
static constexpr double FinalTempRatio = 1e-3;
static constexpr double KickProb = 7e-9;
static constexpr double KickFraction = 0.1;static constexpr i64 StagnationMode = 0;
static constexpr double StagnationShare = 0.03;
static constexpr double ReheatFactor = 30;
static constexpr double ReheatDecayShare = 0.1;
static constexpr double LocalSearchShare = 0.1;
static constexpr i64 ChainNearest = 10;
static constexpr i64 MaxChainDepth = 6;
static constexpr double LocalSearchOps = 2e8;
static constexpr i64 SwapNearest = 5;
static constexpr i64 SubFacilities = 4;
static constexpr i64 MaxSubCustomers = 12;
static constexpr double SubproblemLeaves = 5e4;
static constexpr i64 ExchangeMaxLength = 10;
static constexpr double ImproveEps = 1e-3;
static constexpr double Eps = 1e-9;

namespace {

using TClock = std::chrono::steady_clock;

double Seconds(TClock::time_point since) {
  return std::chrono::duration<double>(TClock::now() - since).count();
}

double CustomerDist(const TCustomer& a, const TCustomer& b) {
  double dx = a.X - b.X, dy = a.Y - b.Y;
  return std::sqrt(dx * dx + dy * dy);
}

double FacilityDist(const TFacility& a, const TFacility& b) {
  double dx = a.X - b.X, dy = a.Y - b.Y;
  return std::sqrt(dx * dx + dy * dy);
}

struct TProblem {
  const std::vector<TFacility>& Facilities;
  const std::vector<TCustomer>& Customers;
  
  std::vector<std::vector<i64>> NearFacilities;
  
  std::vector<std::vector<i64>> NearCustomers;
  
  std::vector<std::vector<i64>> NearFacilitiesOfFacility;
};
struct TState {
  std::vector<i64> Assign;
  std::vector<double> Load;
  std::vector<std::vector<i64>> Members;
  std::vector<i64> PosInMembers;
  
  std::vector<bool> Banned;
  double Total = 0;
};

TState EmptyState(i64 n, i64 m) {
  TState s;
  s.Assign.assign(m, -1);
  s.Load.assign(n, 0);
  s.Members.assign(n, {});
  s.PosInMembers.assign(m, -1);
  s.Banned.assign(n, false);
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
i64 RandomFacility(const TProblem& p, const TState& s, bool open,
                   std::mt19937& rng) {
  std::uniform_int_distribution<i64> fDist(0, std::ssize(p.Facilities) - 1);
  for (int attempt = 0; attempt < 50; attempt++) {
    i64 f = fDist(rng);
    if (open ? !s.Members[f].empty() : (s.Members[f].empty() && !s.Banned[f])) {
      return f;
    }
  }
  return -1;
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

void TryRelocate(const TProblem& p, TState& s, double temp,
                 std::mt19937& rng) {
  std::uniform_int_distribution<i64> cDist(0, std::ssize(p.Customers) - 1);
  i64 c = cDist(rng);
  const auto& near = p.NearFacilities[c];
  std::uniform_int_distribution<i64> fDist(0, std::ssize(near) - 1);
  i64 to = near[fDist(rng)];
  if (to == s.Assign[c] || s.Banned[to] || !Fits(p, s, c, to)) {
    return;
  }
  if (Accept(RelocateDelta(p, s, c, to), temp, rng)) {
    MoveCustomer(p, s, c, to);
  }
}void TrySwap(const TProblem& p, TState& s, double temp, std::mt19937& rng) {
  std::uniform_int_distribution<i64> cDist(0, std::ssize(p.Customers) - 1);
  i64 c1 = cDist(rng);
  i64 f1 = s.Assign[c1];
  const auto& near = p.NearFacilities[c1];
  std::uniform_int_distribution<i64> fDist(0, std::ssize(near) - 1);
  i64 f2 = near[fDist(rng)];
  if (f2 == f1 || s.Banned[f2] || s.Members[f2].empty()) {
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
bool CloseFacility(const TProblem& p, TState& s, i64 f,
                   std::vector<std::pair<i64, i64>>& moves) {
  std::vector<i64> members = s.Members[f];
  std::sort(members.begin(), members.end(), [&](i64 a, i64 b) {
    return p.Customers[a].Demand > p.Customers[b].Demand;
  });

  for (i64 c : members) {
    i64 target = -1;
    for (i64 g : p.NearFacilities[c]) {
      if (g != f && !s.Banned[g] && !s.Members[g].empty() && Fits(p, s, c, g)) {
        target = g;
        break;
      }
    }
    if (target == -1) {
      Revert(p, s, moves);
      moves.clear();
      return false;
    }
    moves.push_back({c, f});
    MoveCustomer(p, s, c, target);
  }
  return true;
}
void OpenFacility(const TProblem& p, TState& s, i64 f,
                  std::vector<std::pair<i64, i64>>& moves) {
  for (i64 c : p.NearCustomers[f]) {
    i64 from = s.Assign[c];
    if (Dist(p.Facilities[f], p.Customers[c]) <
            Dist(p.Facilities[from], p.Customers[c]) &&
        Fits(p, s, c, f)) {
      moves.push_back({c, from});
      MoveCustomer(p, s, c, f);
    }
  }
}

bool SwapFacilities(const TProblem& p, TState& s, i64 f, i64 g,
                    std::vector<std::pair<i64, i64>>& moves) {
  OpenFacility(p, s, g, moves);
  if (moves.empty()) {
    return false;
  }
  return CloseFacility(p, s, f, moves);
}

void TryClose(const TProblem& p, TState& s, double temp, std::mt19937& rng) {
  i64 f = RandomFacility(p, s, true, rng);
  if (f == -1) {
    return;
  }
  double before = s.Total;
  std::vector<std::pair<i64, i64>> moves;
  if (CloseFacility(p, s, f, moves) &&
      !Accept(s.Total - before, temp, rng)) {
    Revert(p, s, moves);
  }
}

void TryOpen(const TProblem& p, TState& s, double temp, std::mt19937& rng) {
  i64 f = RandomFacility(p, s, false, rng);
  if (f == -1) {
    return;
  }
  double before = s.Total;
  std::vector<std::pair<i64, i64>> moves;
  OpenFacility(p, s, f, moves);
  if (!moves.empty() && !Accept(s.Total - before, temp, rng)) {
    Revert(p, s, moves);
  }
}

void TrySwapFacilities(const TProblem& p, TState& s, double temp,
                       std::mt19937& rng) {
  i64 f = RandomFacility(p, s, true, rng);
  if (f == -1) {
    return;
  }
  const auto& near = p.NearFacilitiesOfFacility[f];
  i64 k = std::min<i64>(std::ssize(near), SwapNearest);
  if (k == 0) {
    return;
  }
  std::uniform_int_distribution<i64> gDist(0, k - 1);
  i64 g = near[gDist(rng)];
  if (!s.Members[g].empty() || s.Banned[g]) {
    return;
  }
  double before = s.Total;
  std::vector<std::pair<i64, i64>> moves;
  if (SwapFacilities(p, s, f, g, moves) &&
      !Accept(s.Total - before, temp, rng)) {
    Revert(p, s, moves);
  }
}
void Kick(const TProblem& p, TState& s, std::mt19937& rng) {
  i64 n = p.Facilities.size();
  std::vector<i64> open;
  for (i64 f = 0; f < n; f++) {
    if (!s.Members[f].empty()) {
      open.push_back(f);
    }
  }
  i64 cnt = std::max<i64>(1, std::ssize(open) * KickFraction);
  if (cnt >= std::ssize(open)) {
    return;
  }
  std::shuffle(open.begin(), open.end(), rng);
  for (i64 i = 0; i < cnt; i++) {
    s.Banned[open[i]] = true;
  }

  for (i64 i = 0; i < cnt; i++) {
    std::vector<i64> members = s.Members[open[i]];
    std::sort(members.begin(), members.end(), [&](i64 a, i64 b) {
      return p.Customers[a].Demand > p.Customers[b].Demand;
    });
    for (i64 c : members) {
      i64 target = -1;
      for (i64 g : p.NearFacilities[c]) {
        if (!s.Banned[g] && !s.Members[g].empty() && Fits(p, s, c, g)) {
          target = g;
          break;
        }
      }
      if (target == -1) {
        double bestCost = std::numeric_limits<double>::max();
        for (i64 g = 0; g < n; g++) {
          if (s.Banned[g] || !Fits(p, s, c, g)) {
            continue;
          }
          double cost = Dist(p.Facilities[g], p.Customers[c]);
          if (s.Members[g].empty()) {
            cost += p.Facilities[g].Cost;
          }
          if (cost < bestCost) {
            bestCost = cost;
            target = g;
          }
        }
      }
      if (target != -1) {
        MoveCustomer(p, s, c, target);
      }
    }
  }
}
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

bool ExtendChain(const TProblem& p, TState& s, i64 c, i64 depth,
                 i64 maxDepth, double startTotal) {
  i64 from = s.Assign[c];
  const auto& near = p.NearFacilities[c];
  i64 k = std::min<i64>(std::ssize(near), ChainNearest);
  for (i64 i = 0; i < k; i++) {
    i64 f = near[i];
    if (f == from) {
      continue;
    }
    MoveCustomer(p, s, c, f);
    if (s.Load[f] <= p.Facilities[f].Capacity + Eps) {
      if (s.Total < startTotal - ImproveEps) {
        return true;
      }
    } else if (depth < maxDepth) {
      std::vector<i64> members = s.Members[f];
      for (i64 next : members) {
        if (next == c || s.Load[f] - p.Customers[next].Demand >
                             p.Facilities[f].Capacity + Eps) {
          continue;
        }
        if (ExtendChain(p, s, next, depth + 1, maxDepth, startTotal)) {
          return true;
        }
      }
    }
    MoveCustomer(p, s, c, from);
  }
  return false;
}
i64 ChainDepth(const TProblem& p, const TState& s) {
  i64 m = p.Customers.size();
  i64 open = 0;
  for (const auto& members : s.Members) {
    if (!members.empty()) {
      open++;
    }
  }
  double branching = ChainNearest * static_cast<double>(m) / std::max<i64>(open, 1);
  double ops = static_cast<double>(m) * ChainNearest;
  i64 depth = 1;
  while (depth < MaxChainDepth && ops * branching <= LocalSearchOps) {
    ops *= branching;
    depth++;
  }
  return depth;
}

bool SolveSubproblem(const TProblem& p, TState& s, i64 seed) {
  std::vector<i64> fs = {s.Assign[seed]};
  for (i64 g : p.NearFacilities[seed]) {
    if (std::ssize(fs) >= SubFacilities) {
      break;
    }
    if (g != fs[0] && !s.Members[g].empty()) {
      fs.push_back(g);
    }
  }
  if (std::ssize(fs) < 2) {
    return false;
  }

  i64 k = 1;
  while (k < MaxSubCustomers &&
         std::pow(static_cast<double>(fs.size()), k + 1) <= SubproblemLeaves) {
    k++;
  }

  std::vector<i64> cands;
  for (i64 f : fs) {
    cands.insert(cands.end(), s.Members[f].begin(), s.Members[f].end());
  }
  k = std::min<i64>(k, std::ssize(cands));
  std::partial_sort(cands.begin(), cands.begin() + k, cands.end(),
                    [&](i64 a, i64 b) {
                      return CustomerDist(p.Customers[a], p.Customers[seed]) <
                             CustomerDist(p.Customers[b], p.Customers[seed]);
                    });
  cands.resize(k);

  std::vector<i64> orig(k), cur(k), best(k);
  for (i64 i = 0; i < k; i++) {
    orig[i] = s.Assign[cands[i]];
  }
  double bestTotal = s.Total - ImproveEps;
  for (i64 c : cands) {
    Detach(p, s, c);
  }

  
  
  std::vector<double> restBound(k + 1, 0);
  for (i64 i = k - 1; i >= 0; i--) {
    double minDist = std::numeric_limits<double>::max();
    for (i64 f : fs) {
      minDist = std::min(minDist, Dist(p.Facilities[f], p.Customers[cands[i]]));
    }
    restBound[i] = restBound[i + 1] + minDist;
  }

  bool found = false;
  std::function<void(i64)> dfs = [&](i64 i) {
    if (s.Total + restBound[i] >= bestTotal) {
      return;
    }
    if (i == k) {
      bestTotal = s.Total;
      best = cur;
      found = true;
      return;
    }
    for (i64 f : fs) {
      if (!Fits(p, s, cands[i], f)) {
        continue;
      }
      Attach(p, s, cands[i], f);
      cur[i] = f;
      dfs(i + 1);
      Detach(p, s, cands[i]);
    }
  };
  dfs(0);

  const auto& assign = found ? best : orig;
  for (i64 i = 0; i < k; i++) {
    Attach(p, s, cands[i], assign[i]);
  }
  return found;
}

bool ChainPass(const TProblem& p, TState& s, i64 maxDepth,
               TClock::time_point start, double deadline) {
  bool improved = false;
  for (i64 c = 0; c < std::ssize(p.Customers); c++) {
    if ((c & 15) == 0 && Seconds(start) >= deadline) {
      break;
    }
    if (ExtendChain(p, s, c, 1, maxDepth, s.Total)) {
      improved = true;
    }
  }
  return improved;
}

bool ClosePass(const TProblem& p, TState& s) {
  bool improved = false;
  for (i64 f = 0; f < std::ssize(p.Facilities); f++) {
    if (s.Members[f].empty()) {
      continue;
    }
    double before = s.Total;
    std::vector<std::pair<i64, i64>> moves;
    if (!CloseFacility(p, s, f, moves)) {
      continue;
    }
    if (s.Total < before - ImproveEps) {
      improved = true;
    } else {
      Revert(p, s, moves);
    }
  }
  return improved;
}

bool OpenPass(const TProblem& p, TState& s) {
  bool improved = false;
  for (i64 f = 0; f < std::ssize(p.Facilities); f++) {
    if (!s.Members[f].empty()) {
      continue;
    }
    double before = s.Total;
    std::vector<std::pair<i64, i64>> moves;
    OpenFacility(p, s, f, moves);
    if (moves.empty()) {
      continue;
    }
    if (s.Total < before - ImproveEps) {
      improved = true;
    } else {
      Revert(p, s, moves);
    }
  }
  return improved;
}

bool SwapFacilitiesPass(const TProblem& p, TState& s) {
  bool improved = false;
  for (i64 f = 0; f < std::ssize(p.Facilities); f++) {
    const auto& near = p.NearFacilitiesOfFacility[f];
    i64 k = std::min<i64>(std::ssize(near), SwapNearest);
    for (i64 i = 0; i < k && !s.Members[f].empty(); i++) {
      i64 g = near[i];
      if (!s.Members[g].empty()) {
        continue;
      }
      double before = s.Total;
      std::vector<std::pair<i64, i64>> moves;
      if (!SwapFacilities(p, s, f, g, moves)) {
        continue;
      }
      if (s.Total < before - ImproveEps) {
        improved = true;
        break;
      }
      Revert(p, s, moves);
    }
  }
  return improved;
}

struct TExchange {
  std::vector<i64> Customers;
  
  i64 EndFacility = -1;
};bool FindMultiExchange(const TProblem& p, const TState& s,
                       TExchange& exchange) {
  i64 m = p.Customers.size();
  const double inf = std::numeric_limits<double>::max();
  std::vector<std::vector<double>> cost(ExchangeMaxLength + 1,
                                        std::vector<double>(m, inf));
  std::vector<std::vector<i64>> pred(ExchangeMaxLength + 1,
                                     std::vector<i64>(m, -1));
  std::vector<std::vector<i64>> origin(ExchangeMaxLength + 1,
                                       std::vector<i64>(m, -1));
  std::vector<i64> active(m);
  std::iota(active.begin(), active.end(), 0);
  for (i64 v = 0; v < m; v++) {
    cost[1][v] = 0;
    origin[1][v] = v;
  }

  
  auto onPath = [&](i64 level, i64 v, i64 g) {
    for (i64 t = level; t >= 1; t--) {
      if (s.Assign[v] == g) {
        return true;
      }
      v = pred[t][v];
    }
    return false;
  };

  double bestGain = -ImproveEps;
  i64 bestLevel = -1, bestNode = -1, bestEnd = -1;

  for (i64 k = 1; k <= ExchangeMaxLength && !active.empty(); k++) {
    std::vector<i64> next;
    for (i64 u : active) {
      double c = cost[k][u];
      i64 gu = s.Assign[u];
      double du = p.Customers[u].Demand;
      double base = Dist(p.Facilities[gu], p.Customers[u]);
      i64 first = origin[k][u];
      i64 gFirst = s.Assign[first];

      
      if (k >= 2 && s.Load[gFirst] - p.Customers[first].Demand + du <=
                        p.Facilities[gFirst].Capacity + Eps) {
        double total = c + Dist(p.Facilities[gFirst], p.Customers[u]) - base;
        if (total < bestGain) {
          bestGain = total;
          bestLevel = k;
          bestNode = u;
          bestEnd = -1;
        }
      }

      const auto& near = p.NearFacilities[u];
      i64 kk = std::min<i64>(std::ssize(near), ChainNearest);

      
      for (i64 i = 0; i < kk; i++) {
        i64 f = near[i];
        if (f == gu || s.Banned[f] || !Fits(p, s, u, f) || onPath(k, u, f)) {
          continue;
        }
        double total = c + Dist(p.Facilities[f], p.Customers[u]) - base;
        if (s.Members[f].empty()) {
          total += p.Facilities[f].Cost;
        }
        if (s.Members[gFirst].size() == 1) {
          total -= p.Facilities[gFirst].Cost;
        }
        if (total < bestGain) {
          bestGain = total;
          bestLevel = k;
          bestNode = u;
          bestEnd = f;
        }
      }

      if (k == ExchangeMaxLength) {
        continue;
      }

      
      for (i64 i = 0; i < kk; i++) {
        i64 g = near[i];
        if (g == gu || s.Members[g].empty()) {
          continue;
        }
        double into = c + Dist(p.Facilities[g], p.Customers[u]) - base;
        if (into >= 0 || onPath(k, u, g)) {
          continue;
        }
        for (i64 v : s.Members[g]) {
          if (s.Load[g] - p.Customers[v].Demand + du >
              p.Facilities[g].Capacity + Eps) {
            continue;
          }
          if (into < cost[k + 1][v]) {
            if (cost[k + 1][v] == inf) {
              next.push_back(v);
            }
            cost[k + 1][v] = into;
            pred[k + 1][v] = u;
            origin[k + 1][v] = first;
          }
        }
      }
    }
    active = std::move(next);
  }

  if (bestNode == -1) {
    return false;
  }
  exchange.Customers.assign(bestLevel, -1);
  i64 v = bestNode;
  for (i64 t = bestLevel; t >= 1; t--) {
    exchange.Customers[t - 1] = v;
    v = pred[t][v];
  }
  exchange.EndFacility = bestEnd;
  return true;
}void ApplyExchange(const TProblem& p, TState& s, const TExchange& e,
                   std::vector<i64>& origShops) {
  i64 k = e.Customers.size();
  origShops.assign(k, -1);
  for (i64 i = 0; i < k; i++) {
    origShops[i] = s.Assign[e.Customers[i]];
  }
  for (i64 c : e.Customers) {
    Detach(p, s, c);
  }
  for (i64 i = 0; i + 1 < k; i++) {
    Attach(p, s, e.Customers[i], origShops[i + 1]);
  }
  Attach(p, s, e.Customers[k - 1],
         e.EndFacility == -1 ? origShops[0] : e.EndFacility);
}

void RevertExchange(const TProblem& p, TState& s, const TExchange& e,
                    const std::vector<i64>& origShops) {
  for (i64 c : e.Customers) {
    Detach(p, s, c);
  }
  for (i64 i = 0; i < std::ssize(e.Customers); i++) {
    Attach(p, s, e.Customers[i], origShops[i]);
  }
}bool MultiExchangePass(const TProblem& p, TState& s, TClock::time_point start,
                       double deadline) {
  bool exchanged = false;
  TExchange e;
  while (Seconds(start) < deadline && FindMultiExchange(p, s, e)) {
    double before = s.Total;
    std::vector<i64> origShops;
    ApplyExchange(p, s, e, origShops);
    if (s.Total >= before - ImproveEps) {
      RevertExchange(p, s, e, origShops);
      break;
    }
    exchanged = true;
  }
  return exchanged;
}void LocalSearch(const TProblem& p, TState& s, TClock::time_point start,
                 double deadline, bool useSubproblems, std::mt19937& rng) {
  while (Seconds(start) < deadline) {
    bool improved = ChainPass(p, s, ChainDepth(p, s), start, deadline);
    improved = MultiExchangePass(p, s, start, deadline) || improved;
    improved = ClosePass(p, s) || improved;
    improved = OpenPass(p, s) || improved;
    improved = SwapFacilitiesPass(p, s) || improved;
    if (!improved) {
      break;
    }
  }
  if (!useSubproblems) {
    return;
  }
  std::uniform_int_distribution<i64> cDist(0, std::ssize(p.Customers) - 1);
  for (i64 iter = 0;; iter++) {
    if ((iter & 15) == 0 && Seconds(start) >= deadline) {
      break;
    }
    SolveSubproblem(p, s, cDist(rng));
  }
}

}  

std::vector<i64> Solve(const std::vector<TFacility>& facilities,
                       const std::vector<TCustomer>& customers) {
  i64 n = facilities.size();
  i64 m = customers.size();
  std::mt19937 rng(42);

  auto start = TClock::now();
  auto elapsed = [&]() { return Seconds(start); };

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
  p.NearFacilitiesOfFacility.assign(n, {});
  for (i64 f = 0; f < n; f++) {
    std::vector<i64> ids(n);
    std::iota(ids.begin(), ids.end(), 0);
    ids.erase(ids.begin() + f);
    i64 k = std::min<i64>(n - 1, NearestFacilityFacilities);
    std::partial_sort(ids.begin(), ids.begin() + k, ids.end(),
                      [&](i64 a, i64 b) {
                        return FacilityDist(facilities[f], facilities[a]) <
                               FacilityDist(facilities[f], facilities[b]);
                      });
    ids.resize(k);
    p.NearFacilitiesOfFacility[f] = std::move(ids);
  }

  
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

  
  
  for (i64 round = 0; round < NumRounds; round++) {
    double roundStart = elapsed();
    double roundDuration = (TimeLimit - roundStart) / (NumRounds - round);
    if (roundDuration <= 0) {
      break;
    }
    double annealDuration = roundDuration * (1 - LocalSearchShare);
    double roundEnd = roundStart + roundDuration;

    TState s = BuildState(p, bestAssign);
    double t0 = EstimateTemperature(p, s, rng) * StartTempFactor;
    double t1 = t0 * FinalTempRatio;
    double temp = t0;

    
    
    
    
    double segmentBest = s.Total;
    bool progressed = false;
    double lastProgress = roundStart;
    double lastCheck = roundStart;
    double boost = 1;

    for (i64 iter = 0;; iter++) {
      if ((iter & 255) == 0) {
        double now = elapsed();
        double frac = (now - roundStart) / annealDuration;
        if (frac >= 1) {
          break;
        }
        if (progressed) {
          lastProgress = now;
          progressed = false;
        }
        boost = std::max(1.0, boost * std::pow(ReheatFactor, -(now - lastCheck) /
                                                                 (ReheatDecayShare * annealDuration)));
        lastCheck = now;
        if (StagnationMode != 0 && boost <= 1.0 &&
            now - lastProgress > StagnationShare * annealDuration) {
          if (StagnationMode == 2) {
            std::fill(s.Banned.begin(), s.Banned.end(), false);
            Kick(p, s, rng);
          } else if (StagnationMode == 3) {
            s = BuildState(p, bestAssign);
            Kick(p, s, rng);
          }
          boost = ReheatFactor;
          segmentBest = std::numeric_limits<double>::max();
          lastProgress = now;
        }
        temp = std::min(t0, t0 * std::pow(t1 / t0, frac) * boost);
      }

      double r = uni01(rng);
      if (r < KickProb) {
        Kick(p, s, rng);
      } else if (r < 0.5) {
        TryRelocate(p, s, temp, rng);
      } else if (r < 0.85) {
        TrySwap(p, s, temp, rng);
      } else if (r < 0.93) {
        TryClose(p, s, temp, rng);
      } else if (r < 0.97) {
        TryOpen(p, s, temp, rng);
      } else {
        TrySwapFacilities(p, s, temp, rng);
      }

      if (s.Total < segmentBest - Eps) {
        segmentBest = s.Total;
        progressed = true;
      }
      if (s.Total < bestTotal - Eps) {
        bestTotal = s.Total;
        bestAssign = s.Assign;
      }
    }

    
    
    TState final = BuildState(p, s.Assign);
    LocalSearch(p, final, start, roundEnd, false, rng);
    if (final.Total < bestTotal - Eps) {
      bestTotal = final.Total;
      bestAssign = final.Assign;
    }
    TState best = BuildState(p, bestAssign);
    LocalSearch(p, best, start, roundEnd, true, rng);
    if (best.Total < bestTotal - Eps) {
      bestTotal = best.Total;
      bestAssign = best.Assign;
    }
  }

  return bestAssign;
}
