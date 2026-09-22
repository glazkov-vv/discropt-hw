#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>

#include "check.h"
#include "sol.h"

const std::string TestsFilename = "tests.txt";

struct TTestCase {
  std::vector<TFacility> Facilities;
  std::vector<TCustomer> Customers;
};

TTestCase ReadFromFile(const std::string& filename) {
  std::ifstream in(filename);
  i64 n, m;
  {
    std::string line;
    std::getline(in, line);
    std::istringstream iss(line);
    iss >> n >> m;
  }

  std::vector<TFacility> facilities;
  std::vector<TCustomer> customers;
  std::string line;
  i64 lineNumber = 0;
  while (std::getline(in, line)) {
    if (line.find_first_not_of(" \t\r") == std::string::npos) {
      continue;
    }
    std::istringstream iss(line);
    if (lineNumber < n) {
      TFacility f;
      iss >> f.Cost >> f.Capacity >> f.X >> f.Y;
      f.Id = lineNumber;
      facilities.push_back(f);
    } else if (lineNumber < n + m) {
      TCustomer c;
      iss >> c.Demand >> c.X >> c.Y;
      c.Id = lineNumber - n;
      customers.push_back(c);
    }
    lineNumber++;
  }
  return TTestCase{.Facilities = facilities, .Customers = customers};
}

static double TotalCost(const TTestCase& testCase,
                        const std::vector<i64>& assignment) {
  std::vector<bool> open(testCase.Facilities.size(), false);
  double total = 0;
  for (i64 c = 0; c < std::ssize(assignment); c++) {
    total += Dist(testCase.Facilities[assignment[c]], testCase.Customers[c]);
    open[assignment[c]] = true;
  }
  for (i64 f = 0; f < std::ssize(open); f++) {
    if (open[f]) {
      total += testCase.Facilities[f].Cost;
    }
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
    auto result = Solve(testCase.Facilities, testCase.Customers);

    try {
      Check(testCase.Facilities, testCase.Customers, result);
    } catch (const char* ex) {
      std::cout << "Wrong result on " << name << ": " << ex << std::endl;
      return 0;
    }

    double cost = TotalCost(testCase, result);
    std::cout << std::fixed << std::setprecision(2);
    std::cout << "On test case " << name << " got cost " << cost << ". ";
    if (cost <= hardLimit) {
      std::cout << "Hard limit passed. Score 7" << std::endl;
    } else if (cost <= easyLimit) {
      std::cout << "Easy limit passed. Score 5" << std::endl;
    } else {
      std::cout << "No limit passed. Score 0" << std::endl;
    }
  }
}
