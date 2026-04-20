#include "sol.h"

#include <algorithm>
#include <chrono>
#include <numeric>
#include <random>
#include <vector>

static constexpr double TimeLimit = 55;
static constexpr i64 PopulationSize = 50;

using TOrdering = std::vector<i64>;

static std::vector<i64> GreedyColor(const std::vector<std::vector<i64>>& adj, const TOrdering& order)
{
    i64 n = adj.size();
    std::vector<i64> color(n, -1);
    for (i64 v : order) {
        i64 deg = std::ssize(adj[v]);
        std::vector<bool> forbidden(deg + 2, false);
        for (i64 u : adj[v]) {
            if (color[u] >= 0 && color[u] <= deg) {
                forbidden[color[u]] = true;
            }
        }
        for (i64 c = 0; ; c++) {
            if (!forbidden[c]) {
                color[v] = c;
                break;
            }
        }
    }
    return color;
}

static i64 NumColors(const std::vector<i64>& coloring)
{
    if (coloring.empty()) {
        return 0;
    }
    return *std::max_element(coloring.begin(), coloring.end()) + 1;
}

static TOrdering SumPosCrossover(const TOrdering& p1, const TOrdering& p2)
{
    i64 n = p1.size();
    std::vector<i64> pos1(n), pos2(n);
    for (i64 i = 0; i < n; i++) {
        pos1[p1[i]] = i;
    }
    for (i64 i = 0; i < n; i++) {
        pos2[p2[i]] = i;
    }

    std::vector<i64> elems(n);
    std::iota(elems.begin(), elems.end(), 0);
    std::stable_sort(elems.begin(), elems.end(), [&] (i64 a, i64 b) {
        return pos1[a] + pos2[a] < pos1[b] + pos2[b];
    });
    return elems;
}

static void Mutate(TOrdering& order, std::mt19937& rng)
{
    i64 n = order.size();
    std::uniform_int_distribution<i64> dist(0, n - 1);
    for (int k = 0; k < 3; k++) {
        std::swap(order[dist(rng)], order[dist(rng)]);
    }
}

static i64 CountInversions(std::vector<i64> arr)
{
    i64 n = arr.size();
    if (n <= 1) {
        return 0;
    }
    i64 mid = n / 2;
    std::vector<i64> left(arr.begin(), arr.begin() + mid);
    std::vector<i64> right(arr.begin() + mid, arr.end());
    i64 inv = CountInversions(left) + CountInversions(right);
    i64 i = 0, j = 0, k = 0;
    while (i < std::ssize(left) && j < std::ssize(right)) {
        if (left[i] <= right[j]) {
            arr[k++] = left[i++];
        } else {
            arr[k++] = right[j++];
            inv += std::ssize(left) - i;
        }
    }
    while (i < std::ssize(left)) {
        arr[k++] = left[i++];
    }
    while (j < std::ssize(right)) {
        arr[k++] = right[j++];
    }
    return inv;
}

static i64 Dissimilarity(const TOrdering& a, const TOrdering& b)
{

    i64 n = a.size();
    std::vector<i64> rank_a(n);
    for (i64 i = 0; i < n; i++) {
        rank_a[a[i]] = i;
    }
    std::vector<i64> b_renum(n);
    for (i64 i = 0; i < n; i++) {
        b_renum[i] = rank_a[b[i]];
    }
    return CountInversions(b_renum);
}

static std::vector<std::pair<i64, i64>> PairPopulation(
    const std::vector<TOrdering>& pop, std::mt19937& rng)
{
    const i64 N = pop.size();
    std::vector<bool> alreadyPaired(N, false);
    std::vector<std::pair<i64, i64>> pairs;

    std::vector<i64> order(N);
    std::iota(order.begin(), order.end(), 0);
    std::shuffle(order.begin(), order.end(), rng);

    for (auto idx : order) {
        if (alreadyPaired[idx]) {
            continue;
        }
        alreadyPaired[idx] = true;
        i64 bestPartner = -1;
        i64 bestDiss = -1;
        for (i64 j = 0; j < N; j++) {
            if (alreadyPaired[j]) {
                continue;
            }
            i64 d = Dissimilarity(pop[idx], pop[j]);
            if (d > bestDiss) {
                bestDiss = d;
                bestPartner = j;
            }
        }
        if (bestPartner != -1) {
            alreadyPaired[bestPartner] = true;
            pairs.push_back({idx, bestPartner});
        }
    }
    return pairs;
}

std::vector<i64> Solve(i64 n, const std::vector<TEdge>& edges)
{
    std::vector<std::vector<i64>> adj(n);
    for (const auto& e : edges) {
        adj[e.u].push_back(e.v);
        adj[e.v].push_back(e.u);
    }

    std::mt19937 rng(42);

    std::vector<TOrdering> pop(PopulationSize, TOrdering(n));
    for (auto& ind : pop) {
        std::iota(ind.begin(), ind.end(), 0);
        std::shuffle(ind.begin(), ind.end(), rng);
    }

    std::vector<i64> fitnesses(PopulationSize);
    for (i64 i = 0; i < PopulationSize; i++) {
        fitnesses[i] = NumColors(GreedyColor(adj, pop[i]));
    }

    i64 bestFitness = *std::min_element(fitnesses.begin(), fitnesses.end());
    TOrdering bestOrder = pop[std::min_element(fitnesses.begin(), fitnesses.end()) - fitnesses.begin()];

    auto start = std::chrono::steady_clock::now();
    auto elapsed = [&]() {
        return std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    };

    while (elapsed() < TimeLimit) {
        std::vector<i64> order;
        std::iota(order.begin(), order.end(), 0);
        std::shuffle(order.begin(), order.end(), rng);
        for (auto [i, j] : PairPopulation(pop, rng)) {
            TOrdering child = SumPosCrossover(pop[i], pop[j]);
            Mutate(child, rng);
            i64 fitness = NumColors(GreedyColor(adj, child));
            pop.push_back(std::move(child));
            fitnesses.push_back(fitness);
        }

        i64 newSize = std::ssize(pop);
        std::vector<i64> sortedIdx(newSize);
        std::iota(sortedIdx.begin(), sortedIdx.end(), 0);
        std::sort(sortedIdx.begin(), sortedIdx.end(), [&] (auto a, auto b) {
            return fitnesses[a] < fitnesses[b];
        });

        i64 immune = newSize / 6;
        i64 killBottom = newSize / 6;
        i64 killRandom = newSize - PopulationSize - killBottom;

        std::vector<bool> alive(newSize, true);
        for (i64 i = newSize - killBottom; i < newSize; i++) {
            alive[sortedIdx[i]] = false;
        }

        std::vector<i64> candidates;
        candidates.reserve(newSize - immune - killBottom);
        for (i64 i = immune; i < newSize - killBottom; i++) {
            candidates.push_back(sortedIdx[i]);
        }
        std::shuffle(candidates.begin(), candidates.end(), rng);
        killRandom = std::min(killRandom, static_cast<i64>(std::size(candidates)));
        for (i64 i = 0; i < killRandom; i++) {
            alive[candidates[i]] = false;
        }

        std::vector<TOrdering> newPop;
        std::vector<i64> newFits;
        newPop.reserve(PopulationSize);
        newFits.reserve(PopulationSize);
        for (i64 i = 0; i < newSize; i++) {
            if (alive[i]) {
                newPop.push_back(std::move(pop[i]));
                newFits.push_back(fitnesses[i]);
            }
        }
        pop = std::move(newPop);
        fitnesses = std::move(newFits);

        for (i64 i = 0; i < std::ssize(fitnesses); i++) {
            if (fitnesses[i] < bestFitness) {
                bestFitness = fitnesses[i];
                bestOrder = pop[i];
            }
        }
    }

    return GreedyColor(adj, bestOrder);
}
