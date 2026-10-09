# Suite 17-Layered Foliated Bulk Gravity

Analytic and numerical validation suite for the 17-layered foliated bulk, 8(+) vs 8(-) double octet balance (602.50 GeV radion), 3.33 TeV parity doublet in LHC, PeV neutrino inelastic absorption in IceCube, and 0.19 mHz GW stochastic background in LISA.

Author: Albert Cruanas Pardo (ORCID: 0009-0007-0846-2321)
DOI: 10.5281/zenodo.23135100
- `map_3d_foliation_dynamics.py`: 3D topocentric kinematics and velocity hierarchy across microscopic, galactic, and terrestrial scales.
- `plot_3d_foliation_ellipsoid.py`: Metric sensitivity ellipsoid reconstruction, 3-axis diurnal delay envelopes, and invariant trace conservation (Tr(h) = h0).
- `foliation_3d_trajectory.csv`: Numerical time-series data for the 24-hour sidereal 3D kinematic trajectory and strain tensor components.
- `sim_foliation_frequency_comb_psd.py`: Welch PSD spectral estimation resolving the 17-line discrete harmonic comb (Delta f = 29.41 mHz, SNR = +22.36 dB).
- `sim_chiral_polarization_leakage_cavity.py`: Cavity-finesse resonant amplification of vacuum birefringence (P -> S leakage, +59.97 dB optical contrast).
- `sim_dispersion_and_poisson_tensor.py`: Phononic acoustic dispersion chirp (Delta tau = 41.14 ms) across the 17-line comb and metric Poisson tensor reconstruction (nu_exp = 0.4201, Delta nu_radion ~ 1/sqrt(17)).
- `sim_nonlinear_breather_and_mixing.py`: Non-linear Schrödinger equation (NLSE) propagation via Split-Step Fourier method, resolving higher-order breather soliton dynamics (N_sol = 18.87) and polarization mode dispersion delay (delta_tau_PMD = 600.93 ns).
- `sim_scalar_radion_field_dynamics.py`: 1D Klein-Gordon scalar radion field dynamics, fifth-force spatial gradients across the optical table, and fundamental constant modulation (|Delta alpha / alpha| = 4.99e-17, 50x above optical atomic clock limits).
- `sim_radion_thermodynamics_and_poynting.py`: Scalar radion thermodynamics, Poynting energy flux (S_phi = 3.41e-5 W/m^2), net mirror radiative force (F_phi = 12.1 nN), and Standard Quantum Limit (SQL) margin (+13.66 dB).
- `sim_foliation_four_pillars_dynamics.py`: Complete four-pillar dynamics synthesis: Israel junction stress balance (p_eff = -2.70e14 GeV^4), phononic transfer matrix with Bragg gap at [33.19, 36.20] Hz, geodesic Jacobi scalar force dominance (ratio = 3.33e15), and 500 m distributed network triangulation (tau_max = 87.31 s).
