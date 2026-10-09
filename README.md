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

---

# Sección II: Validación Experimental Acústica y Cotas Métricas (LIGO O3)

Protocolo experimental de bombeo-sonda (*Pump-Probe*) y análisis de límites métricos en ultra-alto vacío para evaluar la transferencia de energía-momento y el desacoplo de halos transitorios.

### Módulos de Análisis Experimental y Metrología

- `test_ligo.py`: Descarga y cálculo de la Densidad Espectral de Amplitud (ASD) del canal de deformación métrica $h(t)$ de LIGO Livingston (L1, Run O3) mediante GWOSC.
- `metricas_ligo_valores.py`: Extracción cuantitativa celda a celda de la sensibilidad métrica en las frecuencias críticas ($13{,}15\text{ Hz}$, $17{,}20\text{ Hz}$, $26{,}30\text{ Hz}$, $210{,}40\text{ Hz}$ y $210{,}42\text{ Hz}$).
- `generar_audios_15min.py`: Síntesis de formas de onda bifrecuencia con corte temporal abrupto ($t = t_{\text{corte}}$) para ensayos de relajación en cavidad macroscópica ($53\text{ m}^3$).
- `ligo_espectro_test.png`: Curva experimental de deformación espectral en LIGO L1 contrastando el muro sísmico frente a la zona de máxima sensibilidad.

### Cotas Cuantitativas Obtenidas

| Canal / Frecuencia | Sensibilidad / Nivel | Régimen Físico |
| :--- | :--- | :--- |
| **Bomba desacoplada ($17{,}20\text{ Hz}$)** | Decaimiento $-8{,}41\text{ dB}$ en $50\text{ s}$ | Relajación acústica lineal del volumen |
| **Límite relativo de mezcla ($\eta$)** | $\le -57{,}0\text{ dBc}$ | Cota superior de modulación cruzada en cavidad |
| **LIGO L1 ($13{,}15\text{ Hz}$)** | $\text{ASD} = 3{,}86 \times 10^{-21}\text{ Hz}^{-1/2}$ | Fondo dominado por aislamiento sísmico |
| **LIGO L1 ($210{,}41\text{ Hz}$)** | $\text{ASD} = 4{,}67 \times 10^{-24}\text{ Hz}^{-1/2}$ | Suelo térmico de espejos en vacío puro |
| **Cota coherente integrada ($T = 1024\text{ s}$)** | $h_0 \le 1{,}46 \times 10^{-25}$ | Cota superior estricta a oscilaciones métricas |

### Modelo de Desacoplo de Halo Transitorio

En medios de alta densidad bariónica ($\rho_{\text{in}} \gg \rho_{\text{ext}}$), el grado de libertad escalar permanece apantallado y canaliza su disipación a través del peine armónico elástico ($16 \times 13{,}15\text{ Hz} = 210{,}40\text{ Hz}$). 

En la discontinuidad de corte ($t = t_{\text{corte}}$), la ruptura del equilibrio adiabático en la frontera interfacial desacopla una componente transitoria hacia el medio difuso circundante:

$$h_{ij}(t, \mathbf{x}) = \Theta(t - t_{\text{corte}}) \, \frac{\mathcal{A}_{\text{trans}}}{\vert{}\mathbf{x}\vert{}} \, \exp\left( -\frac{t - t_{\text{corte}}}{\tau_{\text{halo}}} \right) \cos\big( \omega_{\text{halo}}(t - t_{\text{corte}}) + \delta \big)$$

donde $\tau_{\text{halo}} \gg \tau_{\text{acústico}}$ define la persistencia del halo en ausencia de pérdidas viscosas ordinarias.
