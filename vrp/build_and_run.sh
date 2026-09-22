g++ --std=c++23 -O2 sol.cpp main.cpp check.cpp \
    && ./a.out > report.txt \
    && cat report.txt \
    && rm a.out
