#include <fstream>
#include <iostream>
#include <sstream>

#include "aco.h"
#include "sa.h"
#include "check.h"

const std::string TestsFilename = "tests.txt";

struct TTestCase {
  std::vector<TPoint> Points;
};

TTestCase ReadFromFile(const std::string& filename) {
  std::ifstream in(filename);
  i64 n;

  {
    std::string line;
    std::getline(in, line);
    std::istringstream iss(line);
    iss >> n;
  }

  std::vector<TPoint> points;
  std::string line;
  i64 lineNumber = 0;
  while (std::getline(in, line)) {
    std::istringstream iss(line);
    TPoint p;
    iss >> p.X >> p.Y;
    p.Id = lineNumber;
    points.push_back(std::move(p));
    lineNumber++;
  }

  return TTestCase{.Points = points};
}

static double RouteLength(const std::vector<TPoint>& points,
                          const std::vector<i64>& route) {
  i64 n = route.size();
  double total = 0;
  for (i64 i = 0; i < n; i++) {
    total += Dist(points[route[i]], points[route[(i + 1) % n]]);
  }
  return total;
}

int main() {
  std::ifstream in(TestsFilename);
  std::string line;

  while (std::getline(in, line)) {
    std::istringstream iss(line);
    std::string name;
    double easyLimit;
    double hardLimit;

    iss >> name >> easyLimit >> hardLimit;
    auto testCase = ReadFromFile(name);
    auto result = SolveSA(testCase.Points);

    try {
      Check(testCase.Points, result);
    } catch (const char* ex) {
      std::cout << "Wrong result: " << ex << std::endl;
      return 0;
    }

    double length = RouteLength(testCase.Points, result);

    std::cout << "On test case " << name << " got length " << length << ". ";
    if (length <= hardLimit) {
      std::cout << "Hard limit passed. Score 5" << std::endl;
    } else if (length <= easyLimit) {
      std::cout << "Easy limit passed. Score 3" << std::endl;
    } else {
      std::cout << "No limit passed. Score 0" << std::endl;
    }
  }
}
