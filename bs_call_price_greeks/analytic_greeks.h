/**
 * @file analytic_greeks.h
 * @brief Analytic formulas for Black-Scholes Greeks
 *
 * Implements closed-form solutions for delta and gamma
 * of Black-Scholes European call options.
 */

#ifndef ANALYTIC_GREEKS_H
#define ANALYTIC_GREEKS_H

#include <cmath>
#include <algorithm>
#include <string>

// Black-Scholes call delta: Δ_call = e^{-qT} Φ(d1)
double bs_delta_call(double S, double K, double r, double q, double sigma, double T);

// Black-Scholes call gamma: Γ = e^{-qT} φ(d1) / (S σ sqrt(T))
double bs_gamma_call(double S, double K, double r, double q, double sigma, double T);

// Generic interface: returns delta or gamma based on type parameter
double bs_analytic_call(std::string type, double S, double K, double r, double q, double sigma, double T);

#endif // ANALYTIC_GREEKS_H
