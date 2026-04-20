#include "check.h"
#include "sol.h"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>

const std::string TestsFilename = "tests.txt";

struct TTestCase
{
    i64 N;
    std::vector<TEdge> Edges;
};

TTestCase ReadFromFile(const std::string& filename)
{
    std::ifstream in(filename);
    i64 n, e;
    {
        std::string line;
        std::getline(in, line);
        std::istringstream iss(line);
        iss >> n >> e;
    }
    std::vector<TEdge> edges;
    std::string line;
    while (std::getline(in, line)) {
        std::istringstream iss(line);
        TEdge edge;
        iss >> edge.u >> edge.v;
        edges.push_back(edge);
    }
    return TTestCase{.N = n, .Edges = edges};
}

int main()
{
    std::ifstream in(TestsFilename);
    std::string line;

    while (std::getline(in, line)) {
        std::istringstream iss(line);
        std::string name;
        i64 easyLimit, hardLimit;
        iss >> name >> easyLimit >> hardLimit;

        auto testCase = ReadFromFile(name);
        auto result = Solve(testCase.N, testCase.Edges);

        try {
            Check(testCase.N, testCase.Edges, result);
        } catch (const char* ex) {
            std::cout << "Wrong result on " << name << ": " << ex << std::endl;
            return 0;
        }

        i64 numColors = *std::max_element(result.begin(), result.end()) + 1;
        std::cout << "On test case " << name << " got " << numColors << " colors. ";
        if (numColors <= hardLimit) {
            std::cout << "Hard limit passed. Score 7" << std::endl;
        } else if (numColors <= easyLimit) {
            std::cout << "Easy limit passed. Score 5" << std::endl;
        } else {
            std::cout << "No limit passed. Score 0" << std::endl;
        }
    }
}
