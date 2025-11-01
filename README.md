# Black-Scholes Greeks Calculator

[![C++ Unit Tests](https://github.com/Isolaee/bsAlg/actions/workflows/test.yml/badge.svg)](https://github.com/Isolaee/bsAlg/actions/workflows/test.yml)

A comprehensive C++ implementation of Black-Scholes option Greeks with multiple numerical differentiation methods.

## Features

- **Analytic Greeks**: Closed-form solutions for Delta and Gamma
- **Classical Forward Differences**: Standard finite difference approximations
- **Complex-Step Differentiation**: High-precision numerical derivatives with O(h²) and O(h⁴) accuracy

## Project Structure

```
bsAlg/
├── bs_call_price/              # Black-Scholes pricing functions
├── bs_call_price_greeks/       # Analytic Greek formulas
├── classical_forward_differences/  # Finite difference methods
├── complex_step_differentation/    # Complex-step methods
├── tests/                      # Unit tests
├── output/                     # Generated CSV validation results
└── test_greeks.cpp            # Main validation program
```

## Building

### Compile Tests
```bash
g++ -std=c++11 -o tests/test_greeks_simple \
    tests/test_greeks_simple.cpp \
    bs_call_price_greeks/analytic_greeks.cpp \
    classical_forward_differences/classical_forward_differences.cpp \
    complex_step_differentation/complex_step_differentation.cpp \
    -I.
```

### Compile Validation Program
```bash
g++ -std=c++11 -o test_greeks \
    test_greeks.cpp \
    write_greeks.cpp \
    bs_call_price_greeks/analytic_greeks.cpp \
    classical_forward_differences/classical_forward_differences.cpp \
    complex_step_differentation/complex_step_differentation.cpp \
    -I.
```

## Running

### Run Unit Tests
```bash
./tests/test_greeks_simple
```

Output:
```
=== Running Black-Scholes Greeks Unit Tests ===

--- Analytic Greeks Tests ---
Testing delta bounds... ✓ PASSED
Testing gamma is positive... ✓ PASSED
...
✓ All tests passed!
```

### Generate Validation CSVs
```bash
mkdir -p output
./test_greeks
```

This generates comparative analysis CSV files in `output/`:
- `bs_fd_vs_complex_scenario1.csv` - ATM reference case
- `bs_fd_vs_complex_scenario2.csv` - Near-expiry, low-volatility case

## Methods Comparison

### Delta (Δ)

| Method | Formula | Truncation Error |
|--------|---------|------------------|
| Analytic | Δ = e^(-qT) Φ(d₁) | Exact |
| Forward Difference | [C(S+h) - C(S)] / h | O(h) |
| Complex-Step | Im[C(S+ih)] / h | O(h²) |

### Gamma (Γ)

| Method | Formula | Truncation Error |
|--------|---------|------------------|
| Analytic | Γ = e^(-qT) φ(d₁) / (Sσ√T) | Exact |
| Forward Difference | [C(S+2h) - 2C(S+h) + C(S)] / h² | O(h²) |
| Complex-Step (Real) | -2(Re[C(S+ih)] - C(S)) / h² | O(h²) |
| Complex-Step (45°) | Im[C(S+hω) + C(S-hω)] / h² | O(h⁴) |

where ω = e^(iπ/4) = (1+i)/√2

## Test Coverage

The test suite includes 14 tests:

**Analytic Greeks** (9 tests):
- Delta bounds, known values, edge cases
- Gamma positivity, known values, ATM maximum
- Zero volatility handling
- Put-call parity

**Numerical Methods** (4 tests):
- Forward difference accuracy
- Complex-step machine precision
- 45° complex-step high-order accuracy
- Convergence analysis

**Dispatcher** (1 test):
- String-based function dispatch

## Validation Scenarios

### Scenario 1: ATM Reference
- S = K = 100, r = q = 0, σ = 0.20, T = 1

### Scenario 2: Near-Expiry, Low-Vol
- S = K = 100, r = q = 0, σ = 0.01, T = 1/365

Both scenarios sweep step sizes from h_rel ∈ [10^-16, 10^-1] with 24 logarithmically-spaced points.

## CI/CD

Automated testing runs on every push request via GitHub Actions. The workflow:
1. Compiles the test suite
2. Runs all unit tests
3. Generates validation CSVs
4. Uploads results as artifacts

## Visualization

### Plotting Error Analysis with gnuplot

The project includes gnuplot scripts to visualize the error behavior of different numerical differentiation methods. Since C++ does not have native plotting capabilities, we use **gnuplot** - a powerful command-line plotting utility commonly used in scientific computing.

#### Installation
```bash
sudo apt install gnuplot-qt  # Ubuntu/Debian
```

#### Generate Plots
```bash
gnuplot plotting/plot_greeks_png.gnuplot  # Generates PNG files
gnuplot plotting/plot_greeks.gnuplot      # Interactive Qt windows
```

#### Output
The scripts generate 2×2 grid plots for each scenario showing:
- Delta: Finite Difference vs Complex Step errors
- Gamma: Finite Difference vs Complex Step (Real & 45°) errors

Each plot includes a **red reference line** showing the analytic Greek value for comparison.

All plots use log-log scale (h_rel vs absolute error) to visualize the behavior across step sizes from 10^-16 to 10^-1. The plots clearly demonstrate:
- **FD methods**: Classic U-shaped error curve (truncation error at large h, roundoff at small h)
- **CS methods**: Machine-precision accuracy maintained across most of the range (line breaks at zero error are replaced with 10^-17 floor for visualization)
- **CS 45° method**: Superior accuracy for second derivatives (Gamma)

Generated files: `output/scenario1_errors.png` and `output/scenario2_errors.png`

## Requirements

- C++11 or later

## References

- Black, F., & Scholes, M. (1973). The Pricing of Options and Corporate Liabilities
- Martins, J. R. R. A., Sturdza, P., & Alonso, J. J. (2003). The complex-step derivative approximation
