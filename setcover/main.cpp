#include "check.h"
#include "sol.h"

#include <iostream>
#include <fstream>
#include <sstream>

const std::string TestsFilename = "tests.txt";

struct TTestCase {
    i64 N;
    std::vector<TObject> Objects;
};

TTestCase ReadFromFile(std::string filename) {
    std::ifstream in(filename);
    i64 n, m;

    {
        std::string line;
        std::getline(in, line);
        std::istringstream iss(line);
        iss >> n >> m;
    }

    std::vector<TObject> objects;
    std::string line;
    i64 lineNumber = 0;
    while (std::getline(in, line)) {
        std::istringstream iss(line);
        TObject obj;
        iss >> obj.c;
        i64 x;
        while (iss >> x) {
            obj.ss.push_back(x);
        }
        obj.Id = lineNumber;
        objects.push_back(std::move(obj));
        lineNumber++;
    }
    
    return TTestCase{
        .N = n,
        .Objects = objects};
}

int main() {
    std::ifstream in(TestsFilename);
    std::string line;

    while (std::getline(in, line)) {
        std::istringstream iss(line);
        std::string name;
        i64 easyLimit;
        i64 hardLimit;
        
        iss >> name >> easyLimit >> hardLimit;
        // std::cout << "Got test " << name << easyLimit << hardLimit << std::endl;
        auto testCase = ReadFromFile(name);
        auto result = Solve(testCase.N, testCase.Objects);
        try {
            Check(testCase.N, testCase.Objects, result);
        } catch (const char* ex) {
            std::cout << "Wrong result " << ex << std::endl;
            return 0;
        }

        i64 score = 0;
        for (i64 id : result) {
            score += testCase.Objects[id].c;
        }

        std::cout << "On test case " << name << " got score " << score << ". ";
        if (score <= hardLimit) {
            std::cout << "Hard limit passed. Score 7" << std::endl;
        } else if (score <= easyLimit) {
            std::cout << "Easy limit passed. Score 5" << std::endl;
        } else {
            std::cout << "No limit passed. Score 0" << std::endl;
        }
    }
}