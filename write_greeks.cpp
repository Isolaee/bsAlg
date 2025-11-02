#include "write_greeks.h"
#include "complex_step_differentation/complex_step_differentation.h"
#include "bs_call_price_greeks/analytic_greeks.h"
#include "classical_forward_differences/classical_forward_differences.h"
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>

void write_scenario_csv(const std::string& filename,
                        double S, double K, double r, double q, double sigma, double T) {
    /**
     * Write CSV file comparing FD and complex-step methods across different step sizes.
     * Sweeps h_rel over logarithmic grid [10^-16, 10^-1].
     */
    
    std::ofstream csv(filename);
    if (!csv.is_open()) {
        std::cerr << "Error: Could not open " << filename << " for writing.\n";
        return;
    }
    
    // Write header
    csv << "h_rel,h,";
    csv << "Delta_analytic,Delta_fd,Delta_cs,err_D_fd,err_D_cs,";
    csv << "Gamma_analytic,Gamma_fd,Gamma_cs_real,Gamma_cs_45,";
    csv << "err_G_fd,err_G_cs_real,err_G_cs_45\n";
    
    // Compute analytic Greeks (reference values)
    double delta_analytic = bs_delta_call(S, K, r, q, sigma, T);
    double gamma_analytic = bs_gamma_call(S, K, r, q, sigma, T);
    
    // Set precision for output
    csv << std::scientific << std::setprecision(12);
    
    // Logarithmic grid of relative step sizes h_rel ∈ [10^-16, 10^-4]
    // 24 points as suggested in the validation scenarios
    const int num_points = 24;
    const double log_min = -16.0;  // 10^-16
    const double log_max = -4.0;   // 10^-4
    
    for (int i = 0; i < num_points; ++i) {
        // Logarithmic spacing
        double log_h_rel = log_min + i * (log_max - log_min) / (num_points - 1);
        double h_rel = std::pow(10.0, log_h_rel);
        double h = h_rel * S;  // Absolute step size: h = h_rel * S
        
        // Forward differences
        double delta_fd = delta_fwd(S, K, r, q, sigma, T, h);
        double gamma_fd = gamma_fwd(S, K, r, q, sigma, T, h);
        
        // Complex-step methods
        double delta_cs = delta_complex_step(S, K, r, q, sigma, T, h);
        double gamma_cs_real = gamma_complex_step(S, K, r, q, sigma, T, h);
        double gamma_cs_45 = gamma_complex_step_45deg(S, K, r, q, sigma, T, h);
        
        // Compute absolute errors
        double err_D_fd = std::abs(delta_fd - delta_analytic);
        double err_D_cs = std::abs(delta_cs - delta_analytic);
        double err_G_fd = std::abs(gamma_fd - gamma_analytic);
        double err_G_cs_real = std::abs(gamma_cs_real - gamma_analytic);
        double err_G_cs_45 = std::abs(gamma_cs_45 - gamma_analytic);
        
        // Write row
        csv << h_rel << "," << h << ",";
        csv << delta_analytic << "," << delta_fd << "," << delta_cs << ",";
        csv << err_D_fd << "," << err_D_cs << ",";
        csv << gamma_analytic << "," << gamma_fd << "," << gamma_cs_real << "," << gamma_cs_45 << ",";
        csv << err_G_fd << "," << err_G_cs_real << "," << err_G_cs_45 << "\n";
    }
    
    csv.close();
    std::cout << "Written: " << filename << " (" << num_points << " points)\n";
}
