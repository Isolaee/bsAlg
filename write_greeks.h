/**
 * @file write_greeks.h
 * @brief CSV output utilities for Greeks validation
 */

#ifndef WRITE_GREEKS_H
#define WRITE_GREEKS_H

#include <string>

// Write CSV file comparing FD and complex-step methods across different step sizes
void write_scenario_csv(const std::string& filename,
                        double S, double K, double r, double q, double sigma, double T);

#endif // WRITE_GREEKS_H
