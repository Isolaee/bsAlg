/**
 * @file classical_forward_differences.h
 * @brief Classical forward difference methods for numerical derivatives
 *
 * Implements forward difference approximations for delta and gamma
 * of Black-Scholes call options.
 */

#ifndef CLASSICAL_FORWARD_DIFFERENCES_H
#define CLASSICAL_FORWARD_DIFFERENCES_H

#include <cmath>
#include <algorithm>

// Generic forward difference for first derivative
double classical_forward_difference(double (*f)(double), double x, double h);

// Forward difference approximation for delta: Δ_fwd(S; h) = [C(S+h) - C(S)] / h
double delta_fwd(double S, double K, double r, double q, double sigma, double T, double h);

// Forward difference approximation for gamma: Γ_fwd(S; h) = [C(S+2h) - 2C(S+h) + C(S)] / h²
double gamma_fwd(double S, double K, double r, double q, double sigma, double T, double h);

#endif // CLASSICAL_FORWARD_DIFFERENCES_H
