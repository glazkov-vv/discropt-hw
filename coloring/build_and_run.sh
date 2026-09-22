g++ --std=c++23 sol_tabucol.cpp main.cpp check.cpp \
    && ./a.out \
    && cat report.txt \
    && rm a.out