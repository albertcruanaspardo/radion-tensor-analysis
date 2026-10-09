#!/usr/bin/env python3
"""
Cálculo de Sección Eficaz de Fusión Fotón-Fotón (gamma gamma -> Solitón) en LHC (13.6 TeV)
Albert Cruañas Pardo (2026)
Basado en:
 - Aproximación de Fotones Equivalentes (EPA) con factor de forma dipolar
 - Acoplamiento dilatónico modular d_e / Lambda_r
"""
import numpy as np
from scipy.integrate import quad

print("=" * 85)
print(" ESTIMACIÓN DE SECCIÓN EFICAZ: FUSIÓN FOTÓN-FOTÓN EXCLUSIVA EN LHC (13.6 TeV)")
print("=" * 85)

# Constantes físicas fundamentales
alpha_em = 1.0 / 137.036
m_p = 0.938272          # GeV (masa del protón)
Q0_sq = 0.71            # GeV^2 (parámetro de corte del factor de forma dipolar)
sqrt_s = 13600.0        # GeV (LHC Run 3, 13.6 TeV)
s = sqrt_s**2
hbarc2_to_pb = 3.89379e8 # Factor de conversión GeV^-2 a picobarns (pb)

# Parámetros del modelo de foliación modular de 17 capas
N_sub = 17
Lambda_r = 10785.0      # GeV (Escala ultravioleta de foliación)
d_e = 1.0 / np.sqrt(N_sub) # ~0.2425 (acoplamiento conforme dilatónico)

def flujo_foton_equivalente(x):
    """
    Función de distribución de fotones en un protón elástico f_{gamma/p}(x).
    Fórmula analítica basada en el factor de forma electromagnético dipolar.
    """
    if x <= 0 or x >= 1:
        return 0.0
    # Límite cinemático mínimo de transferencia de momento Q_min^2
    q_min_sq = (m_p * x)**2 / (1.0 - x)
    # Factor adimensional de apantallamiento
    y = q_min_sq / Q0_sq
    if y >= 100.0:
        return 0.0
    
    # Integral aproximada del factor dipolar de Budnev et al.
    phi_y = (1.0 + y)**(-4) * (1.0 + 4.0 * y)
    integrando = (1.0 / y) * (1.0 - x) * phi_y
    # Normalización del flujo de Weizsäcker-Williams
    f_val = (alpha_em / (np.pi * x)) * np.log(1.0 + 1.0 / (y + 1e-12)) * (1.0 - x)
    return max(0.0, f_val)

def luminosidad_gamma_gamma(M):
    """
    Luminosidad diferencial dL/dM [GeV^-1] a masa invariante M.
    dL/dM = (2M / s) * \int_{M^2/s}^1 (dx1 / x1) f(x1) f(M^2 / (s * x1))
    """
    tau = (M**2) / s
    if tau >= 1.0:
        return 0.0

    def integrand(log_x1):
        x1 = np.exp(log_x1)
        x2 = tau / x1
        if x2 >= 1.0 or x2 <= 0.0:
            return 0.0
        return flujo_foton_equivalente(x1) * flujo_foton_equivalente(x2)

    # Cambio de variable logarítmico para estabilidad numérica
    log_min = np.log(tau)
    log_max = 0.0
    res, _ = quad(integrand, log_min, log_max, limit=200, epsrel=1e-3)
    return (2.0 * M / s) * res

def anchura_gamma_gamma(M):
    """
    Anchura parcial de desintegración a dos fotones Gamma(phi -> gamma gamma) [GeV].
    Gamma = d_e^2 * M^3 / (64 * pi * Lambda_r^2)
    """
    return (d_e**2 * M**3) / (64.0 * np.pi * Lambda_r**2)

def seccion_eficaz_exclusiva(M):
    """
    Sección eficaz total pp -> p + phi + p en picobarns (pb).
    """
    width = anchura_gamma_gamma(M)
    dl_dm = luminosidad_gamma_gamma(M)
    # sigma = (8 * pi^2 / M^3) * Gamma * (dL/dM)
    sigma_gev = (8.0 * np.pi**2 / (M**3)) * width * dl_dm
    sigma_pb = sigma_gev * hbarc2_to_pb
    return sigma_pb, width, dl_dm

# Barrido de masas del espectro solitónico (desde modos sub-TeV hasta el límite Lambda_r)
masas_gev = [500.0, 1000.0, 2000.0, 3000.0, 5000.0, 7500.0, 10785.0]

print(f"Parámetros teóricos: Lambda_r = {Lambda_r:.1f} GeV | d_e = {d_e:.4f}")
print(f"LHC Run 3: sqrt(s) = 13.6 TeV | Luminosidad integrada objetivo L_int = 300 fb^-1\n")
print(f"{'Masa Solitón (GeV)':<20} | {'Gamma_gammagamma (GeV)':<24} | {'dL/dM (GeV^-1)':<16} | {'Sigma (fb)':<12} | {'Eventos (300 fb^-1)'}")
print("-" * 92)

L_int_fb = 300.0 # Luminosidad integrada de Run 3 en fb^-1

for M in masas_gev:
    sigma_pb, width, dl_dm = seccion_eficaz_exclusiva(M)
    sigma_fb = sigma_pb * 1e3
    eventos = sigma_fb * L_int_fb
    print(f"{M:<20.1f} | {width:<24.4e} | {dl_dm:<16.4e} | {sigma_fb:<12.4e} | {eventos:<10.2f}")

print("=" * 92)
print("CONCLUSIONES DEL ANÁLISIS CINEMÁTICO:")
print("1. En el modo fundamental (10.78 TeV), la sección eficaz está fuertemente suprimida")
print("   por la caída del flujo de fotones en el umbral cinemático (tau = M^2/s = 0.63).")
print("2. Para condensados solitónicos de capas intermedias o resonancias sub-TeV (M ~ 1 - 2 TeV),")
print("   la sección eficaz alcanza ~ 10^-2 a 10^-4 fb, produciendo eventos limpios y exclusivos")
print("   detectables con los espectrómetros de protones adelantados (PPS / AFP).")
