g++ --std=c++23 -O3 sol.cpp main.cpp aco.cpp sa.cpp check.cpp \
    && ./a.out > report.txt \
    && cat report.txt \
    && rm a.out