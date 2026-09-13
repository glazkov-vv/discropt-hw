#include "check.h"

void Check(const std::vector<TFacility>& facilities, const std::vector<TCustomer>& customers, const std::vector<i64>& assignment)
{
    i64 n = facilities.size();
    i64 m = customers.size();
    if ((i64)assignment.size() != m) throw "Assignment size mismatch";
    std::vector<double> load(n, 0);
    for (i64 c = 0; c < m; c++) {
        i64 f = assignment[c];
        if (f < 0 || f >= n) {
            throw "Invalid facility id";
        }
        load[f] += customers[c].Demand;
    }
    for (i64 f = 0; f < n; f++) {
        if (load[f] > facilities[f].Capacity + 1e-6) {
            throw "Facility capacity exceeded";
        }
    }
}
