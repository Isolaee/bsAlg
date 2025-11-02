#include "classical_forward_differences.h"
#include "../bs_call_price/bs_call_price.h"

double classical_forward_difference(double (*f)(double), double x, double h) {
    /**
     * Computes the first derivative of f at x using the classical forward difference method.
     *
     * @param f Function pointer to the function whose derivative is to be computed.
     * @param x Point at which the derivative is evaluated.
     * @param h Step size.
     * @return Approximation of f'(x).
     */
    return (f(x + h) - f(x)) / h;
}

// Forward difference approximation for delta: Δ_fwd(S; h) = [C(S+h) - C(S)] / h
double delta_fwd(double S, double K, double r, double q, double sigma, double T, double h) {
    /**
     * Computes delta using classical forward difference.
     * Δ_fwd(S; h) = [C(S+h) - C(S)] / h
     * where C(·) = bs_price_call(·, K, r, q, σ, T)
     *
     * @param S     Spot price
     * @param K     Strike price
     * @param r     Risk-free rate
     * @param q     Dividend yield
     * @param sigma Volatility
     * @param T     Time to maturity
     * @param h     Step size
     * @return      Forward difference approximation of delta
     */
    double C_S = bs_price_call(S, K, r, q, sigma, T);
    double C_S_plus_h = bs_price_call(S + h, K, r, q, sigma, T);
    return (C_S_plus_h - C_S) / h;
}

// Forward difference approximation for gamma: Γ_fwd(S; h) = [C(S+2h) - 2C(S+h) + C(S)] / h²
double gamma_fwd(double S, double K, double r, double q, double sigma, double T, double h) {
    /**
     * Computes gamma using classical forward difference.
     * Γ_fwd(S; h) = [C(S+2h) - 2C(S+h) + C(S)] / h²
     * where C(·) = bs_price_call(·, K, r, q, σ, T)
     *
     * @param S     Spot price
     * @param K     Strike price
     * @param r     Risk-free rate
     * @param q     Dividend yield
     * @param sigma Volatility
     * @param T     Time to maturity
     * @param h     Step size
     * @return      Forward difference approximation of gamma
     */
    const double C_S = bs_price_call(S, K, r, q, sigma, T);
    const double C_S_plus_h = bs_price_call(S + h, K, r, q, sigma, T);
    const double C_S_plus_2h = bs_price_call(S + 2.0 * h, K, r, q, sigma, T);
    return (C_S_plus_2h - 2.0 * C_S_plus_h + C_S) / (h * h);
}