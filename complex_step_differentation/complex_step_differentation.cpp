#include "complex_step.h"
#include "../bs_call_price/bs_call_price.h"
#include <complex>
#include <cmath>

// Complex-valued Black-Scholes call price helper
// Uses first-order complex approximation for normal CDF:
//   Φ(x + iy) ≈ Φ(x) + iy φ(x)
static std::complex<double> bs_price_call_complex(
    const std::complex<double>& S,
    const std::complex<double>& K,
    const std::complex<double>& r,
    const std::complex<double>& q,
    const std::complex<double>& sigma,
    const std::complex<double>& T) {
    
    using std::exp;
    using std::log;
    using std::sqrt;
    
    // Discount factor and forward price
    std::complex<double> DF = exp(-r * T);
    std::complex<double> F = S * exp((r - q) * T);
    
    // Total volatility: σ√T
    std::complex<double> sigmaT = sigma * sqrt(T);
    
    // log(F/K)
    std::complex<double> ln_F_over_K = log(F / K);
    
    // d1 and d2
    std::complex<double> d1 = (ln_F_over_K + std::complex<double>(0.5, 0.0) * sigma * sigma * T) / sigmaT;
    std::complex<double> d2 = d1 - sigmaT;
    
    // Complex normal CDF approximation: Φ(x + iy) ≈ Φ(x) + iy φ(x)
    auto Phi_complex = [](const std::complex<double>& z) -> std::complex<double> {
        double x_real = std::real(z);
        double y_imag = std::imag(z);
        return std::complex<double>(Phi_real(x_real), y_imag * phi(x_real));
    };
    
    return DF * (F * Phi_complex(d1) - K * Phi_complex(d2));
}

// First derivative using complex-step method
// f'(x) ≈ Im[f(x + ih)] / h
// Truncation error: O(h²)
double delta_complex_step(double S, double K, double r, double q, double sigma, double T, double h) {
    /**
     * Computes delta using complex-step differentiation.
     * Formula: Δ ≈ Im[C(S + ih)] / h
     * where C(·) = bs_price_call(·, K, r, q, σ, T)
     *
     * @param S     Spot price
     * @param K     Strike price
     * @param r     Risk-free rate
     * @param q     Dividend yield
     * @param sigma Volatility
     * @param T     Time to maturity
     * @param h     Imaginary step size
     * @return      Complex-step approximation of delta
     */
    
    // Create complex spot price: S + ih
    std::complex<double> S_complex(S, h);
    std::complex<double> K_complex(K, 0.0);
    std::complex<double> r_complex(r, 0.0);
    std::complex<double> q_complex(q, 0.0);
    std::complex<double> sigma_complex(sigma, 0.0);
    std::complex<double> T_complex(T, 0.0);
    
    // Evaluate price at S + ih
    std::complex<double> price_complex = bs_price_call_complex(
        S_complex, K_complex, r_complex, q_complex, sigma_complex, T_complex);
    
    // Extract imaginary part and divide by h
    return std::imag(price_complex) / h;
}

// Second derivative using complex-step method
// f''(x) ≈ -2(Re[f(x + ih)] - f(x)) / h²
// Truncation error: O(h²)
double gamma_complex_step(double S, double K, double r, double q, double sigma, double T, double h) {
    /**
     * Computes gamma using complex-step differentiation.
     * Formula: Γ ≈ -2(Re[C(S + ih)] - C(S)) / h²
     * where C(·) = bs_price_call(·, K, r, q, σ, T)
     *
     * @param S     Spot price
     * @param K     Strike price
     * @param r     Risk-free rate
     * @param q     Dividend yield
     * @param sigma Volatility
     * @param T     Time to maturity
     * @param h     Imaginary step size
     * @return      Complex-step approximation of gamma
     */
    
    // Evaluate f(x)
    double f_x = bs_price_call(S, K, r, q, sigma, T);
    
    // Create complex spot price: S + ih
    std::complex<double> S_complex(S, h);
    std::complex<double> K_complex(K, 0.0);
    std::complex<double> r_complex(r, 0.0);
    std::complex<double> q_complex(q, 0.0);
    std::complex<double> sigma_complex(sigma, 0.0);
    std::complex<double> T_complex(T, 0.0);
    
    // Evaluate f(x + ih)
    std::complex<double> f_x_plus_ih = bs_price_call_complex(
        S_complex, K_complex, r_complex, q_complex, sigma_complex, T_complex);
    
    // Extract real part and apply formula
    double real_part = std::real(f_x_plus_ih);
    return -2.0 * (real_part - f_x) / (h * h);
}

// 45° imaginary-step alternative for second derivative
// f''(x) ≈ Im[f(x + hω) + f(x - hω)] / h²
// where ω = e^(iπ/4) = (1 + i)/√2
// Truncation error: O(h⁴)
double gamma_complex_step_45deg(double S, double K, double r, double q, double sigma, double T, double h) {
    /**
     * Computes gamma using 45° complex-step differentiation.
     * Formula: Γ ≈ Im[C(S + hω) + C(S - hω)] / h²
     * where ω = e^(iπ/4) = (1 + i)/√2
     *
     * @param S     Spot price
     * @param K     Strike price
     * @param r     Risk-free rate
     * @param q     Dividend yield
     * @param sigma Volatility
     * @param T     Time to maturity
     * @param h     Step size
     * @return      45° complex-step approximation of gamma (O(h⁴) accuracy)
     */
    
    // ω = e^(iπ/4) = (1 + i)/√2
    static const double inv_sqrt2 = 0.70710678118654752440; // 1/√2
    const std::complex<double> omega(inv_sqrt2, inv_sqrt2);
    
    std::complex<double> K_complex(K, 0.0);
    std::complex<double> r_complex(r, 0.0);
    std::complex<double> q_complex(q, 0.0);
    std::complex<double> sigma_complex(sigma, 0.0);
    std::complex<double> T_complex(T, 0.0);
    
    // Compute shift: h*ω
    std::complex<double> shift = std::complex<double>(h, 0.0) * omega;
    
    // Evaluate at S + hω
    std::complex<double> S_plus(S, 0.0);
    std::complex<double> f_plus = bs_price_call_complex(
        S_plus + shift, K_complex, r_complex, q_complex, sigma_complex, T_complex);
    
    // Evaluate at S - hω
    std::complex<double> f_minus = bs_price_call_complex(
        S_plus - shift, K_complex, r_complex, q_complex, sigma_complex, T_complex);
    
    // Extract imaginary part of sum and divide by h²
    double imag_sum = std::imag(f_plus + f_minus);
    return imag_sum / (h * h);
}
