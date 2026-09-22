#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>

#include "check.h"
#include "sol.h"

const std::string TestsFilename = "tests.txt";

struct TTestCase {
  std::vector<TPoint> Points;
  i64 Vehicles;
  double Capacity;
};

TTestCase ReadFromFile(const std::string& filename) {
  std::ifstream in(filename);
  i64 n, v;
  double capacity;
  {
    std::string line;
    std::getline(in, line);
    std::istringstream iss(line);
    iss >> n >> v >> capacity;
  }

  std::vector<TPoint> points;
  std::string line;
  i64 lineNumber = 0;
  while (std::getline(in, line) && lineNumber < n) {
    if (line.find_first_not_of(" \t\r") == std::string::npos) {
      continue;
    }
    std::istringstream iss(line);
    TPoint p;
    iss >> p.Demand >> p.X >> p.Y;
    p.Id = lineNumber;
    points.push_back(p);
    lineNumber++;
  }

  return TTestCase{.Points = points, .Vehicles = v, .Capacity = capacity};
}

static double TotalLength(const std::vector<TPoint>& points,
                          const std::vector<std::vector<i64>>& routes) {
  double total = 0;
  for (const auto& route : routes) {
    if (route.empty()) {
      continue;
    }
    total += Dist(points[0], points[route.front()]);
    for (i64 i = 0; i + 1 < std::ssize(route); i++) {
      total += Dist(points[route[i]], points[route[i + 1]]);
    }
    total += Dist(points[route.back()], points[0]);
  }
  return total;
}

int main() {
  std::ifstream in(TestsFilename);
  std::string line;

  while (std::getline(in, line)) {
    std::istringstream iss(line);
    std::string name;
    double easyLimit, hardLimit;
    iss >> name >> easyLimit >> hardLimit;

    auto testCase = ReadFromFile(name);
    auto result = Solve(testCase.Points, testCase.Vehicles, testCase.Capacity);

    try {
      Check(testCase.Points, testCase.Vehicles, testCase.Capacity, result);
    } catch (const char* ex) {
      std::cout << "Wrong result on " << name << ": " << ex << std::endl;
      return 0;
    }

    double length = TotalLength(testCase.Points, result);
    std::cout << std::fixed << std::setprecision(2);
    std::cout << "On test case " << name << " got length " << length << ". ";
    if (length <= hardLimit) {
      std::cout << "Hard limit passed. Score 7" << std::endl;
    } else if (length <= easyLimit) {
      std::cout << "Easy limit passed. Score 5" << std::endl;
    } else {
      std::cout << "No limit passed. Score 0" << std::endl;
    }
  }
}
