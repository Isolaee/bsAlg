# Black-Scholes Greeks Calculator

[![C++ Unit Tests](https://github.com/Isolaee/bsAlg/actions/workflows/test.yml/badge.svg)](https://github.com/Isolaee/bsAlg/actions/workflows/test.yml)

A comprehensive C++ implementation of Black-Scholes option Greeks with multiple numerical differentiation methods. Purpose of this solution is to compare accuracy of different delta and gamma calculation mehtods.

## Features

- **Analytic Greeks**: Closed-form solutions for Delta and Gamma
- **Classical Forward Differences**: Standard finite difference approximations
- **Complex-Step Differentiation**: High-precision numerical derivatives with O(h²) and O(h⁴) accuracy

## Project Structure

```
bsAlg/
├── bs_call_price/                  # Black-Scholes pricing functions
├── bs_call_price_greeks/           # Analytic Greek formulas
├── classical_forward_differences/  # Finite difference methods
├── complex_step_differentation/    # Complex-step methods
├── tests/                          # Unit tests
├── output/                         # Generated CSV validation results
├── plotting/                       # Gnuplot scripts for plotting
├── test_greeks.cpp                 # Main validation program
└── write_greeks.cpp                # Write CSV program
```

## Building

### Requirements

- C++11 or later

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

### Compile Tests
```bash
g++ -std=c++11 -o tests/test_greeks_simple \
    tests/test_greeks_simple.cpp \
    bs_call_price_greeks/analytic_greeks.cpp \
    classical_forward_differences/classical_forward_differences.cpp \
    complex_step_differentation/complex_step_differentation.cpp \
    -I.
```

## Running

### Run program
```bash
./test_greeks
```

Output:
```
=== Black-Scholes Greeks Validation ===
.
.
.
```

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
The scripts generate grid plots for each scenario showing:
- Delta: Finite Difference vs Complex Step errors
- Gamma: Finite Difference vs Complex Step (Real & 45°) errors

Each plot includes a **red reference line** showing the analytic Greek value for comparison.

All plots use log-log scale (h_rel vs absolute error) to visualize the behavior across step sizes from 10^-16 to 10^-1. The plots clearly demonstrate:
- **FD methods**: Classic U-shaped error curve (truncation error at large h, roundoff at small h)
- **CS methods**: Machine-precision accuracy maintained across most of the range (line breaks at zero error are replaced with 10^-17 floor for visualization)
- **CS 45° method**: Superior accuracy for second derivatives (Gamma)

Generated files: `output/scenario1_errors.png` and `output/scenario2_errors.png`

## Conclusion

This validation study compares numerical differentiation methods for computing Black-Scholes Greeks under both standard and stress conditions.

### Accuracy Analysis

#### Delta Computation

Complex-step differentiation demonstrates superior accuracy compared to finite difference methods in both validation scenarios, maintaining machine-precision errors (≈10^-16) across step sizes from 10^-16 to 10^-4. The method begins to degrade only when the relative step size exceeds 10^-4, at which point truncation error dominates.

In contrast, finite difference methods exhibit optimal accuracy at h_rel ≈ 10^-8 with errors around 10^-9, representing a 1,000-fold degradation compared to complex-step approaches.

#### Gamma Computation

For second derivatives, the 45° complex-step method proves optimal, consistently outperforming both standard complex-step (real part) and finite difference implementations. This method maintains machine precision across the entire practical step-size range, demonstrating its theoretical O(h⁴) truncation error advantage.

Finite difference methods for Gamma show significantly larger errors (≈10^-9 to 10^-6) and require careful step-size selection near h_rel ≈ 7×10^-6 to balance roundoff and truncation errors.

### Step-Size Sensitivity

#### Finite Difference Methods

Both Delta and Gamma implementations using finite differences exhibit the characteristic U-shaped error curve on log-log scale:
- **Left regime (h < 10^-8)**: Roundoff error dominates, scaling as O(ε/h) for first derivatives and O(ε/h²) for second derivatives
- **Right regime (h > 10^-6)**: Truncation error dominates, scaling as O(h) for Delta and O(h²) for Gamma
- **Optimal region**: Narrow band near h_rel ≈ 10^-8 for Delta and 7×10^-6 for Gamma

#### Complex-Step Methods

Complex-step approaches demonstrate remarkable robustness, maintaining near-zero error across 8-10 orders of magnitude in step size (10^-16 to 10^-6). This flat error profile eliminates the need for scenario-specific step-size tuning and provides consistent accuracy regardless of market conditions.

### Numerical Stability

The stress scenario (near-expiry, low-volatility with T = 1/365, σ = 0.01) presents significant challenges for finite difference methods:
- Gamma values increase 733-fold compared to the reference case
- Relative errors amplify proportionally
- Optimal step sizes remain constant, but absolute errors degrade

Complex-step methods maintain their accuracy advantage even under these extreme conditions, with the 45° method continuing to deliver machine-precision results.

## Recommendations

### Practical Implementation Guidance

**Primary recommendation**: Use complex-step differentiation with h_rel ∈ [10^-10, 10^-6] for both Delta and Gamma computations. The optimal practical choice is h_rel ≈ 10^-8, which balances numerical precision with computational stability.

**Method-specific recommendations**:
- **Delta**: Complex-step standard method (imaginary part)
- **Gamma**: Complex-step 45° method for highest accuracy

**Fallback strategy**: If complex arithmetic is unavailable, finite difference methods require careful step-size selection:
- Delta: h_rel ≈ 10^-8
- Gamma: h_rel ≈ 7×10^-6

These values must be validated for each specific scenario, particularly in stress conditions with near-expiry or extreme volatility regimes.

**Production considerations**: For risk management systems requiring high accuracy, complex-step methods should be preferred despite their modest computational overhead (approximately 2× compared to finite differences). The elimination of step-size tuning and robust performance across market regimes justify this cost in production environments.

## Use of AI-tools
During the completion of this assignment, I utilized GitHub Copilot as an assistive tool for coding, information retrieval, and README formatting. All solutions and decisions presented are my own and remain my sole responsibility.

## References

- Black, F., & Scholes, M. (1973). The Pricing of Options and Corporate Liabilities
- Martins, J. R. R. A., Sturdza, P., & Alonso, J. J. (2003). The complex-step derivative approximation
