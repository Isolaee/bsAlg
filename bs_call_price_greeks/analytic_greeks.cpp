#include "analytic_greeks.h"
#include "../bs_call_price/bs_call_price.h"
#include <cmath>
#include <string>
#include <limits>

using namespace std;

// Black-Scholes call delta: Δ_call = e^{-qT} Φ(d1)
double bs_delta_call(double S, double K, double r, double q, double sigma, double T) {
    /**
     * Calculates the Black-Scholes delta for a European call option.
     * @param S     Spot price
     * @param K     Strike price
     * @param r     Continuously compounded risk-free interest rate
     * @param q     Continuous dividend yield
     * @param sigma Annualized volatility
     * @param T     Time to maturity
     * @return      Delta of the call option
     */

    const double sigmaT = sigma * std::sqrt(std::max(T, 0.0));
    const double DFq = std::exp(-q * T); // e^{-qT}

    // Handle zero-vol / zero-time as in bs_price_call: delta -> e^{-qT} * 1_{F>K}
    const double F = S * std::exp((r - q) * T);

    if (sigmaT == 0.0) {
        return DFq * (F > K ? 1.0 : 0.0);
    }

    double ln_F_over_K;
    if (K > 0.0) {
        const double x = (F - K) / K;
        ln_F_over_K = (std::abs(x) <= 1e-12) ? std::log1p(x) : std::log(F / K);
    } else {
        ln_F_over_K = std::log(F / K);
    }

    const double d1 = (ln_F_over_K + 0.5 * sigma * sigma * T) / sigmaT;
    return DFq * Phi_real(d1);
}

// Black-Scholes call gamma: Γ = e^{-qT} φ(d1) / (S σ sqrt(T))
double bs_gamma_call(double S, double K, double r, double q, double sigma, double T) {
    /**
     * Calculates the Black-Scholes gamma for a European call option.
     * @param S     Spot price
     * @param K     Strike price
     * @param r     Continuously compounded risk-free interest rate
     * @param q     Continuous dividend yield
     * @param sigma Annualized volatility
     * @param T     Time to maturity
     * @return      Gamma of the call option
     */
    const double sigmaT = sigma * std::sqrt(std::max(T, 0.0));

    // If zero volatility or zero time to maturity, classical gamma is zero
    if (sigmaT == 0.0) return 0.0;

    const double F = S * std::exp((r - q) * T);

    double ln_F_over_K;
    if (K > 0.0) {
        const double x = (F - K) / K;
        ln_F_over_K = (std::abs(x) <= 1e-12) ? std::log1p(x) : std::log(F / K);
    } else {
        ln_F_over_K = std::log(F / K);
    }

    const double d1 = (ln_F_over_K + 0.5 * sigma * sigma * T) / sigmaT;

    // Compute phi(d1) via log form to avoid underflow: log φ = -0.5 d1^2 - 0.5 log(2π)
    static constexpr double NEG_HALF_LOG_2PI = -0.91893853320467274178; // -0.5*log(2π)
    const double log_phi = -0.5 * d1 * d1 + NEG_HALF_LOG_2PI;
    const double phi_d1 = std::exp(log_phi);

    return std::exp(-q * T) * phi_d1 / (S * sigmaT);
}

double bs_analytic_call(string type, double S, double K, double r, double q, double sigma, double T) {
    /**
     * Returns an analytic Greek for a European call option.
     * Accepted `type` values: "delta" or "gamma".
     * @param type  Which Greek to return: "delta" | "gamma"
     * @param S     Spot price
     * @param K     Strike price
     * @param r     Continuously compounded risk-free interest rate
     * @param q     Continuous dividend yield
     * @param sigma Annualized volatility
     * @param T     Time to maturity
     * @return      The requested Greek value, or NaN if `type` is invalid
     */

    if (type == "delta") {
        return bs_delta_call(S, K, r, q, sigma, T);
    }

    if (type == "gamma") {
        return bs_gamma_call(S, K, r, q, sigma, T);
    }

    // Invalid type: return NaN to signal error to caller
    return std::numeric_limits<double>::quiet_NaN();
}
