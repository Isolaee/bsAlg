#!/usr/bin/gnuplot

# Set CSV as the data file format
set datafile separator ","

# Generate PNG output files
set terminal pngcairo size 1400,1000 font "Arial,10"

# Plot Scenario 1
set output "output/scenario1_errors.png"
set multiplot layout 2,2 title "Scenario 1: ATM Reference (S=K=100, r=q=0, σ=0.20, T=1)" font "Arial,12"

# Delta FD Error
set logscale xy
set xlabel "h_{rel}"
set ylabel "Absolute Error"
set title "Delta: Finite Difference Error"
set grid
# Get analytic delta value from first data row
stats "output/bs_fd_vs_complex_scenario1.csv" every ::1::1 using 3 nooutput
delta_analytic = STATS_max
plot "output/bs_fd_vs_complex_scenario1.csv" every ::1 using 1:6 with linespoints pt 7 ps 0.5 lw 2 title "FD Error", \
     delta_analytic with lines dt 2 lw 1.5 lc rgb "red" title sprintf("Analytic Δ = %.4f", delta_analytic)

# Delta CS Error
set logscale xy
set xlabel "h_{rel}"
set ylabel "Absolute Error"
set title "Delta: Complex Step Error"
set grid
plot "output/bs_fd_vs_complex_scenario1.csv" every ::1 using 1:($7 == 0 ? 1e-17 : $7) with linespoints pt 7 ps 0.5 lw 2 title "CS Error", \
     delta_analytic with lines dt 2 lw 1.5 lc rgb "red" title sprintf("Analytic Δ = %.4f", delta_analytic)

# Gamma FD Error
set logscale xy
set xlabel "h_{rel}"
set ylabel "Absolute Error"
set title "Gamma: Finite Difference Error"
set grid
# Get analytic gamma value from first data row
stats "output/bs_fd_vs_complex_scenario1.csv" every ::1::1 using 8 nooutput
gamma_analytic = STATS_max
plot "output/bs_fd_vs_complex_scenario1.csv" every ::1 using 1:12 with linespoints pt 7 ps 0.5 lw 2 title "FD Error", \
     gamma_analytic with lines dt 2 lw 1.5 lc rgb "red" title sprintf("Analytic Γ = %.4f", gamma_analytic)

# Gamma CS Errors
set logscale xy
set xlabel "h_{rel}"
set ylabel "Absolute Error"
set title "Gamma: Complex Step Errors"
set grid
plot "output/bs_fd_vs_complex_scenario1.csv" every ::1 using 1:($13 == 0 ? 1e-17 : $13) with linespoints pt 7 ps 0.5 lw 2 title "CS Real Error", \
     "output/bs_fd_vs_complex_scenario1.csv" every ::1 using 1:($14 == 0 ? 1e-17 : $14) with linespoints pt 5 ps 0.5 lw 2 title "CS 45° Error", \
     gamma_analytic with lines dt 2 lw 1.5 lc rgb "red" title sprintf("Analytic Γ = %.4f", gamma_analytic)

unset multiplot

# Plot Scenario 2
set output "output/scenario2_errors.png"
set multiplot layout 2,2 title "Scenario 2: Near-Expiry, Low-Vol, ATM (S=K=100, r=q=0, σ=0.01, T=1/365)" font "Arial,12"

# Delta FD Error
set logscale xy
set xlabel "h_{rel}"
set ylabel "Absolute Error"
set title "Delta: Finite Difference Error"
set grid
# Get analytic delta value from first data row
stats "output/bs_fd_vs_complex_scenario2.csv" every ::1::1 using 3 nooutput
delta_analytic2 = STATS_max
plot "output/bs_fd_vs_complex_scenario2.csv" every ::1 using 1:6 with linespoints pt 7 ps 0.5 lw 2 title "FD Error", \
     delta_analytic2 with lines dt 2 lw 1.5 lc rgb "red" title sprintf("Analytic Δ = %.4f", delta_analytic2)

# Delta CS Error
set logscale xy
set xlabel "h_{rel}"
set ylabel "Absolute Error"
set title "Delta: Complex Step Error"
set grid
plot "output/bs_fd_vs_complex_scenario2.csv" every ::1 using 1:($7 == 0 ? 1e-17 : $7) with linespoints pt 7 ps 0.5 lw 2 title "CS Error", \
     delta_analytic2 with lines dt 2 lw 1.5 lc rgb "red" title sprintf("Analytic Δ = %.4f", delta_analytic2)

# Gamma FD Error
set logscale xy
set xlabel "h_{rel}"
set ylabel "Absolute Error"
set title "Gamma: Finite Difference Error"
set grid
# Get analytic gamma value from first data row
stats "output/bs_fd_vs_complex_scenario2.csv" every ::1::1 using 8 nooutput
gamma_analytic2 = STATS_max
plot "output/bs_fd_vs_complex_scenario2.csv" every ::1 using 1:12 with linespoints pt 7 ps 0.5 lw 2 title "FD Error", \
     gamma_analytic2 with lines dt 2 lw 1.5 lc rgb "red" title sprintf("Analytic Γ = %.4f", gamma_analytic2)

# Gamma CS Errors
set logscale xy
set xlabel "h_{rel}"
set ylabel "Absolute Error"
set title "Gamma: Complex Step Errors"
set grid
plot "output/bs_fd_vs_complex_scenario2.csv" every ::1 using 1:($13 == 0 ? 1e-17 : $13) with linespoints pt 7 ps 0.5 lw 2 title "CS Real Error", \
     "output/bs_fd_vs_complex_scenario2.csv" every ::1 using 1:($14 == 0 ? 1e-17 : $14) with linespoints pt 5 ps 0.5 lw 2 title "CS 45° Error", \
     gamma_analytic2 with lines dt 2 lw 1.5 lc rgb "red" title sprintf("Analytic Γ = %.4f", gamma_analytic2)

unset multiplot

print "Plots saved to output/scenario1_errors.png and output/scenario2_errors.png"
