#include "../bs_call_price_greeks/analytic_greeks.h"
#include "../classical_forward_differences/classical_forward_differences.h"
#include "../complex_step_differentation/complex_step.h"
#include <iostream>
#include <cmath>
#include <cassert>
#include <iomanip>

// Test counter
int tests_passed = 0;
int tests_failed = 0;

// Helper function to compare doubles
bool approx_equal(double a, double b, double tolerance = 1e-10) {
    return std::abs(a - b) < tolerance;
}

void test_delta_bounds() {
    std::cout << "Testing delta bounds... ";
    
    // Delta should be between 0 and exp(-qT) for calls
    double S = 100.0, K = 100.0, r = 0.05, q = 0.02, sigma = 0.2, T = 1.0;
    double delta = bs_delta_call(S, K, r, q, sigma, T);
    double max_delta = std::exp(-q * T);
    
    assert(delta >= 0.0 && "Delta should be non-negative");
    assert(delta <= max_delta && "Delta should be <= exp(-qT)");
    
    std::cout << "✓ PASSED\n";
    tests_passed++;
}

void test_gamma_positive() {
    std::cout << "Testing gamma is positive... ";
    
    double S = 100.0, K = 100.0, r = 0.05, q = 0.02, sigma = 0.2, T = 1.0;
    double gamma = bs_gamma_call(S, K, r, q, sigma, T);
    
    assert(gamma > 0.0 && "Gamma should be positive for ATM options");
    
    std::cout << "✓ PASSED\n";
    tests_passed++;
}

void test_delta_known_value() {
    std::cout << "Testing delta against known value... ";
    
    // ATM with r=q=0
    double S = 100.0, K = 100.0, r = 0.0, q = 0.0, sigma = 0.2, T = 1.0;
    double delta = bs_delta_call(S, K, r, q, sigma, T);
    double expected = 0.5398278372770;  // Known analytical value
    
    assert(approx_equal(delta, expected, 1e-10) && "Delta should match known value");
    
    std::cout << "✓ PASSED (delta = " << std::setprecision(12) << delta << ")\n";
    tests_passed++;
}

void test_gamma_known_value() {
    std::cout << "Testing gamma against known value... ";
    
    // ATM with r=q=0
    double S = 100.0, K = 100.0, r = 0.0, q = 0.0, sigma = 0.2, T = 1.0;
    double gamma = bs_gamma_call(S, K, r, q, sigma, T);
    double expected = 0.01984762737385;  // Known analytical value
    
    assert(approx_equal(gamma, expected, 1e-10) && "Gamma should match known value");
    
    std::cout << "✓ PASSED (gamma = " << std::setprecision(12) << gamma << ")\n";
    tests_passed++;
}

void test_delta_deep_itm() {
    std::cout << "Testing delta for deep in-the-money... ";
    
    // Deep ITM: S >> K
    double S = 150.0, K = 100.0, r = 0.05, q = 0.02, sigma = 0.2, T = 1.0;
    double delta = bs_delta_call(S, K, r, q, sigma, T);
    
    // Delta should be close to exp(-qT) for deep ITM
    assert(delta > 0.95 && "Deep ITM delta should be high");
    assert(delta <= std::exp(-q * T) && "Delta bounded by exp(-qT)");
    
    std::cout << "✓ PASSED (delta = " << std::setprecision(6) << delta << ")\n";
    tests_passed++;
}

void test_delta_deep_otm() {
    std::cout << "Testing delta for deep out-of-the-money... ";
    
    // Deep OTM: S << K
    double S = 50.0, K = 100.0, r = 0.05, q = 0.02, sigma = 0.2, T = 1.0;
    double delta = bs_delta_call(S, K, r, q, sigma, T);
    
    assert(delta < 0.01 && "Deep OTM delta should be near zero");
    
    std::cout << "✓ PASSED (delta = " << std::setprecision(6) << delta << ")\n";
    tests_passed++;
}

void test_gamma_maximum_atm() {
    std::cout << "Testing gamma is maximum at-the-money... ";
    
    double K = 100.0, r = 0.05, q = 0.02, sigma = 0.2, T = 1.0;
    
    double gamma_atm = bs_gamma_call(100.0, K, r, q, sigma, T);
    double gamma_far_itm = bs_gamma_call(150.0, K, r, q, sigma, T);
    double gamma_far_otm = bs_gamma_call(50.0, K, r, q, sigma, T);
    
    assert(gamma_atm > gamma_far_itm && "ATM gamma should be > far ITM gamma");
    assert(gamma_atm > gamma_far_otm && "ATM gamma should be > far OTM gamma");
    
    std::cout << "✓ PASSED\n";
    tests_passed++;
}

void test_forward_difference_accuracy() {
    std::cout << "Testing forward difference vs analytic... ";
    
    double S = 100.0, K = 100.0, r = 0.05, q = 0.02, sigma = 0.2, T = 1.0;
    double h = 0.01;
    
    double delta_analytic = bs_delta_call(S, K, r, q, sigma, T);
    double delta_fd = delta_fwd(S, K, r, q, sigma, T, h);
    
    double error = std::abs(delta_fd - delta_analytic);
    assert(error < 0.001 && "FD delta error should be < 0.001 for h=0.01");
    
    std::cout << "✓ PASSED (error = " << std::scientific << error << ")\n";
    tests_passed++;
}

void test_complex_step_accuracy() {
    std::cout << "Testing complex-step vs analytic... ";
    
    double S = 100.0, K = 100.0, r = 0.05, q = 0.02, sigma = 0.2, T = 1.0;
    double h = 1e-8;
    
    double delta_analytic = bs_delta_call(S, K, r, q, sigma, T);
    double delta_cs = delta_complex_step(S, K, r, q, sigma, T, h);
    
    double error = std::abs(delta_cs - delta_analytic);
    assert(error < 1e-10 && "Complex-step delta should be nearly exact");
    
    std::cout << "✓ PASSED (error = " << std::scientific << error << ")\n";
    tests_passed++;
}

void test_complex_step_gamma_45deg() {
    std::cout << "Testing 45° complex-step gamma... ";
    
    double S = 100.0, K = 100.0, r = 0.05, q = 0.02, sigma = 0.2, T = 1.0;
    double h = 1e-6;
    
    double gamma_analytic = bs_gamma_call(S, K, r, q, sigma, T);
    double gamma_cs_45 = gamma_complex_step_45deg(S, K, r, q, sigma, T, h);
    
    double error = std::abs(gamma_cs_45 - gamma_analytic);
    assert(error < 1e-6 && "45° complex-step gamma should be accurate");
    
    std::cout << "✓ PASSED (error = " << std::scientific << error << ")\n";
    tests_passed++;
}

void test_zero_volatility() {
    std::cout << "Testing zero volatility edge case... ";
    
    double S = 110.0, K = 100.0, r = 0.05, q = 0.02, sigma = 0.0, T = 1.0;
    
    // For zero vol, ITM option: delta should be exp(-qT), gamma should be 0
    double delta = bs_delta_call(S, K, r, q, sigma, T);
    double gamma = bs_gamma_call(S, K, r, q, sigma, T);
    
    assert(approx_equal(delta, std::exp(-q * T), 1e-10) && "Zero vol ITM: delta = exp(-qT)");
    assert(approx_equal(gamma, 0.0, 1e-10) && "Zero vol: gamma = 0");
    
    std::cout << "✓ PASSED\n";
    tests_passed++;
}

void test_put_call_parity_delta() {
    std::cout << "Testing put-call parity for delta... ";
    
    // Delta_put = Delta_call - exp(-qT)
    double S = 100.0, K = 100.0, r = 0.05, q = 0.02, sigma = 0.2, T = 1.0;
    double delta_call = bs_delta_call(S, K, r, q, sigma, T);
    
    // For put: delta_put = delta_call - exp(-qT)
    double delta_put_expected = delta_call - std::exp(-q * T);
    
    // Delta put should be negative for ATM
    assert(delta_put_expected < 0.0 && "Put delta should be negative");
    assert(delta_put_expected > -1.0 && "Put delta should be > -1");
    
    std::cout << "✓ PASSED\n";
    tests_passed++;
}

void test_convergence_fd_to_cs() {
    std::cout << "Testing FD converges to complex-step... ";
    
    double S = 100.0, K = 100.0, r = 0.05, q = 0.02, sigma = 0.2, T = 1.0;
    
    // As h decreases, FD should approach complex-step accuracy
    double h_large = 0.1;
    double h_small = 0.001;
    
    double delta_cs = delta_complex_step(S, K, r, q, sigma, T, 1e-8);
    double error_large = std::abs(delta_fwd(S, K, r, q, sigma, T, h_large) - delta_cs);
    double error_small = std::abs(delta_fwd(S, K, r, q, sigma, T, h_small) - delta_cs);
    
    assert(error_small < error_large && "Smaller h should give smaller error");
    
    std::cout << "✓ PASSED\n";
    tests_passed++;
}

void test_bs_analytic_call_dispatcher() {
    std::cout << "Testing bs_analytic_call dispatcher... ";
    
    double S = 100.0, K = 100.0, r = 0.05, q = 0.02, sigma = 0.2, T = 1.0;
    
    double delta = bs_analytic_call("delta", S, K, r, q, sigma, T);
    double gamma = bs_analytic_call("gamma", S, K, r, q, sigma, T);
    double invalid = bs_analytic_call("invalid", S, K, r, q, sigma, T);
    
    assert(!std::isnan(delta) && "Delta should be valid");
    assert(!std::isnan(gamma) && "Gamma should be valid");
    assert(std::isnan(invalid) && "Invalid type should return NaN");
    
    std::cout << "✓ PASSED\n";
    tests_passed++;
}

int main() {
    std::cout << "\n=== Running Black-Scholes Greeks Unit Tests ===\n\n";
    
    // Analytic Greeks tests
    std::cout << "--- Analytic Greeks Tests ---\n";
    test_delta_bounds();
    test_gamma_positive();
    test_delta_known_value();
    test_gamma_known_value();
    test_delta_deep_itm();
    test_delta_deep_otm();
    test_gamma_maximum_atm();
    test_zero_volatility();
    test_put_call_parity_delta();
    
    // Numerical methods tests
    std::cout << "\n--- Numerical Methods Tests ---\n";
    test_forward_difference_accuracy();
    test_complex_step_accuracy();
    test_complex_step_gamma_45deg();
    test_convergence_fd_to_cs();
    
    // Dispatcher tests
    std::cout << "\n--- Dispatcher Tests ---\n";
    test_bs_analytic_call_dispatcher();
    
    // Summary
    std::cout << "\n=== Test Summary ===\n";
    std::cout << "Tests passed: " << tests_passed << "\n";
    std::cout << "Tests failed: " << tests_failed << "\n";
    
    if (tests_failed == 0) {
        std::cout << "\n✓ All tests passed!\n\n";
        return 0;
    } else {
        std::cout << "\n✗ Some tests failed!\n\n";
        return 1;
    }
}
