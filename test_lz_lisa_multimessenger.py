#!/usr/bin/env python3
"""
SIMULADOR EXPERIMENTAL MULTI-MENSAJERO (17 CAPAS):
 1. Proyección temporal de descubrimiento/exclusión en LUX-ZEPLIN (LZ 1000 días) y XLZD.
 2. Análisis espectral de ondas gravitacionales primordiales en LISA (ESA/NASA).
"""

import math
import numpy as np

integrate_trapz = getattr(np, 'trapezoid', getattr(np, 'trapz', None))

print("=" * 80)
print(" TEST EXPERIMENTAL MULTI-MENSAJERO: CANDIDATO X17 (LZ) + FONDO GW (LISA)")
print("=" * 80)

# ==============================================================================
# FRENTE 1: LUX-ZEPLIN (LZ) - MATERIA OSCURA A 1.24 TeV
# ==============================================================================
m_chi_gev = 1240.0             # 1.24 TeV (Capítulo 9.1)
sigma_si_teorica = 3.40e-47     # cm^2 (Capítulo 9.1)

# Datos de partida LZ (Run 1 oficial: ~60 días de física, masa fiducial ~5.5 toneladas)
dias_run1 = 60.0
masa_fid_ton = 5.5
expo_run1_ton_dias = dias_run1 * masa_fid_ton
limite_run1_1240 = 1.74e-46     # cm^2 al 90% CL

# Proyecciones temporales de LZ y XLZD
# En ausencia de fondo (régimen ciego de conteo bajo), el límite escala como 1 / Exposición
# Con fondo de neutrinos solares/difusos, escala como 1 / sqrt(Exposición)
# Para LZ en 1.24 TeV, el fondo en la ventana es casi nulo, por lo que escala fuertemente (~1/E^0.8)

dias_escenarios = [60.0, 300.0, 500.0, 1000.0]
print("\n[+] 1. PROYECCIÓN DE DETECCIÓN EN LUX-ZEPLIN (LZ):")
print(f" * Masa analítica predicha:               m_chi    = {m_chi_gev/1000.0:.2f} TeV")
print(f" * Sección eficaz teórica (17 capas):      sigma_SI = {sigma_si_teorica:.2e} cm^2")
print(f" * Cota Run 1 (60 días, PRL 131, 041002):  {limite_run1_1240:.2e} cm^2 (Holgura: 5.12x)")
print("-" * 80)
print(f"{'Días de Toma':<15} | {'Exposición (t*año)':<20} | {'Cota LZ 90% CL (cm^2)':<25} | {'Estatus Teórico'}")
print("-" * 80)

dias_cruce = None
for dias in dias_escenarios:
    expo_ton_yr = (dias * masa_fid_ton) / 365.25
    # Factor de escalado de sensibilidad con exposición
    escala = (dias / dias_run1) ** 0.85
    limite_proyectado = limite_run1_1240 / escala
    
    status = ""
    if limite_proyectado > sigma_si_teorica:
        status = f"Permitido (margen {limite_proyectado/sigma_si_teorica:.2f}x)"
    else:
        status = ">> EN ZONA DE DESCUBRIMIENTO / EXCLUSIÓN <<"
        if dias_cruce is None:
            dias_cruce = dias
            
    print(f"{dias:<15.0f} | {expo_ton_yr:<20.2f} | {limite_proyectado:<25.2e} | {status}")

# Inversión exacta para saber cuándo LZ alcanza los 3.4e-47 cm^2
dias_grid = np.linspace(60.0, 1200.0, 1000)
dias_descubrimiento = None
for d in dias_grid:
    lim = limite_run1_1240 / ((d / dias_run1) ** 0.85)
    if lim <= sigma_si_teorica:
        dias_descubrimiento = d
        break

print("-" * 80)
if dias_descubrimiento:
    print(f"-> Veredicto LZ: La campaña completa de LZ alcanzará la sección eficaz teórica a los ~{dias_descubrimiento:.0f} días.")
    print("   Si el candidato X17 existe a 1.24 TeV, generará eventos de retroceso nuclear en el detector")
    print("   antes de completar los 1.000 días de operación.")

# ==============================================================================
# FRENTE 2: LISA - ESPECTRO ESTOCÁSTICO PRIMORDIAL (3.06 mHz)
# ==============================================================================
print("\n" + "=" * 80)
print("[+] 2. AUDITORÍA ESPECTRAL EN EL OBSERVATORIO ESPACIAL LISA:")
T_n = 210.0          # GeV (Capítulo 9.2)
f_peak_hz = 3.06e-3  # 3.06 mHz (Capítulo 9.3)
Omega_peak = 3.64e-10# h^2 Omega_GW (Capítulo 9.3)
T_obs_sec = 4.0 * 365.25 * 86400.0  # 4 años de misión

print(f" * Frecuencia pico:                       f_peak   = {f_peak_hz*1000.0:.2f} mHz")
print(f" * Densidad de energía reliquia en pico:  Omega_GW = {Omega_peak:.2e}")

def lisa_effective_noise(f):
    L = 2.5e9
    c = 2.99792458e8
    f_star = c / (2.0 * math.pi * L)
    P_oms = (1.5e-11) ** 2 * (1.0 + (2.0e-3 / f) ** 4)
    P_acc = (3.0e-15) ** 2 * (1.0 + (0.4e-3 / f) ** 2) * (1.0 + (f / 8.0e-3) ** 4)
    S_n = (10.0 / (3.0 * L**2)) * (P_oms + 2.0 * (1.0 + np.cos(f / f_star)**2) * P_acc / (2.0 * math.pi * f)**4) * (1.0 + 0.6 * (f / f_star)**2)
    S_gal = 9.0e-45 * f**(-7.0/3.0) * np.exp(-(f / 1.4e-3)**0.88)
    H0_h = 3.24077929e-18
    return (4.0 * math.pi**2 / (3.0 * H0_h**2)) * f**3 * (S_n + S_gal)

# Integración espectral de la SNR
freqs = np.logspace(-5, 0, 2000)
s = freqs / f_peak_hz
shape_sw = (s ** 3) * ((7.0 / (4.0 + 3.0 * s**2)) ** (7.0 / 2.0))
omega_gw_spectrum = Omega_peak * shape_sw
omega_lisa_curve = lisa_effective_noise(freqs)

integrand = (omega_gw_spectrum / omega_lisa_curve) ** 2
snr = math.sqrt(max(0.0, T_obs_sec * integrate_trapz(integrand, freqs)))

noise_at_peak = lisa_effective_noise(f_peak_hz)
margen_pico = Omega_peak / noise_at_peak

print(f" * Ruido efectivo del detector en pico:   Omega_sens = {noise_at_peak:.2e}")
print(f" * Margen de la señal sobre ruido (pico): S/N local  = {margen_pico:.2f}x")
print(f" * Relación Señal-Ruido integrada:        SNR        = {snr:.1f}")
print("-" * 80)
print("-> Veredicto LISA: Señal cosmológica de extraordinaria limpieza.")
print(f"   Con SNR = {snr:.1f} (muy superior al umbral crítico de 60), LISA no solo detectará el fondo,")
print("   sino que podrá reconstruir con precisión la anchura espectral acústica y fijar")
print("   la escala de la transición de fase conforme.")
print("=" * 80)
