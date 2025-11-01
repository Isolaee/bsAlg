#include "write_greeks.h"
#include <iostream>

int main() {
    std::cout << "=== Black-Scholes Greeks Validation ===\n\n";
    
    // Scenario 1: ATM reference
    // S = 100, K = 100, r = q = 0, σ = 0.20, T = 1
    {
        const double S = 100.0;
        const double K = 100.0;
        const double r = 0.0;
        const double q = 0.0;
        const double sigma = 0.20;
        const double T = 1.0;
        
        std::cout << "Scenario 1 (ATM reference):\n";
        std::cout << "  S = " << S << ", K = " << K << ", r = q = " << r 
                  << ", σ = " << sigma << ", T = " << T << "\n";
        write_scenario_csv("output/bs_fd_vs_complex_scenario1.csv", S, K, r, q, sigma, T);
    }
    
    // Scenario 2: Near-expiry, low-vol, ATM
    // S = K = 100, r = q = 0, σ = 0.01, T = 1/365
    {
        const double S = 100.0;
        const double K = 100.0;
        const double r = 0.0;
        const double q = 0.0;
        const double sigma = 0.01;
        const double T = 1.0 / 365.0;  // 1 day to expiry
        
        std::cout << "\nScenario 2 (Near-expiry, low-vol, ATM):\n";
        std::cout << "  S = K = " << S << ", r = q = " << r 
                  << ", σ = " << sigma << ", T = " << T << " (1/365)\n";
        write_scenario_csv("output/bs_fd_vs_complex_scenario2.csv", S, K, r, q, sigma, T);
    }
    
    std::cout << "\nCSV files generated successfully.\n";
    std::cout << "Each file contains data sweeping h_rel over [10^-16, 10^-1] with 24 logarithmically-spaced points.\n";
    
    return 0;
}
