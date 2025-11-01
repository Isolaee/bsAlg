/**
 * @file complex_step.h
 * @brief Complex-step differentiation methods for numerical derivatives
 *
 * Implements complex-step methods for computing derivatives with
 * truncation error O(h²) for first derivative and O(h²) or O(h⁴) 
 * for second derivative.
 */

#ifndef COMPLEX_STEP_H
#define COMPLEX_STEP_H

#include <complex>
#include <cmath>

// Complex-step first derivative: f'(x) ≈ Im[f(x + ih)] / h
// Truncation error: O(h²)
double delta_complex_step(double S, double K, double r, double q, double sigma, double T, double h);

// Complex-step second derivative: f''(x) ≈ -2(Re[f(x + ih)] - f(x)) / h²
// Truncation error: O(h²)
double gamma_complex_step(double S, double K, double r, double q, double sigma, double T, double h);

// 45° imaginary-step alternative for second derivative
// f''(x) ≈ Im[f(x + hω) + f(x - hω)] / h² where ω = e^(iπ/4) = (1+i)/√2
// Truncation error: O(h⁴)
double gamma_complex_step_45deg(double S, double K, double r, double q, double sigma, double T, double h);

#endif // COMPLEX_STEP_H
