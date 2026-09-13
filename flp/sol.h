#pragma once

#include "common.h"

#include <vector>

// Returns assignment: result[c] = facility assigned to customer c (0-indexed)
std::vector<i64> Solve(const std::vector<TFacility>& facilities, const std::vector<TCustomer>& customers);
