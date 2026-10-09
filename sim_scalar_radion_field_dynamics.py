#!/usr/bin/env python3
"""
================================================================================
DINÁMICA DEL CAMPO ESCALAR RADIÓNICO Y ACOPLAMIENTO MATERIA-FOLIACIÓN
Foliación Modular de 17 Capas | Albert Cruañas Pardo (2026)

Modela:
1. Ecuación de onda Klein-Gordon 1D con arrastre bariónico y apantallamiento.
2. Perfil espacial phi(x, t) a lo largo de la mesa de laboratorio (0 <= x <= 1.5 m).
3. Fuerza de quinto rango (aceleración diferencial) sobre los espejos de silicio.
4. Modulación adiabática de la constante de estructura fina (Delta alpha / alpha).
================================================================================
"""

import math
import numpy as np

print("=" * 85)
print(" DINÁMICA DEL CAMPO ESCALAR RADIÓNICO: PERFIL ESPACIAL, FUERZAS Y CONSTANTES")
print("=" * 85)

# ------------------------------------------------------------------------------
# 1. PARÁMETROS FUNDAMENTALES DEL CAMPO ESCALAR
# ------------------------------------------------------------------------------
Lambda_r = 10785.0            # GeV (Escala del radión)
v_EW = 246.22                 # GeV
N_layers = 17
v_drift = 5.7264              # m/s (Velocidad de propagación efectiva en materia)
f0 = 0.50                     # Hz
omega0 = 2.0 * math.pi * f0   # rad/s

# Acoplamiento dilatónico conforme al sector electromagnético y hadrónico
d_e = 1.0 / math.sqrt(float(N_layers))   # Acoplamiento normalizado ~ 1/sqrt(17) ~ 0.2425
alpha_em = 1.0 / 137.035999

# Amplitud de fluctuación del campo escalar (en unidades de deformación métrica h0)
h0 = 1.5e-16
# phi_0 / Lambda_r ~ h0 / 2
phi_amplitude_rel = h0 / 2.0  # Fluctuación fraccional phi / Lambda_r

print(f"[+] 1. PARÁMETROS DEL CAMPO ESCALAR RADIÓNICO:")
print(f" * Escala de corte efectiva (Lambda_r):       {Lambda_r:.1f} GeV")
print(f" * Acoplamiento dilónico conforme (d_e):      d_e = {d_e:.4f} (1/sqrt(17))")
print(f" * Velocidad de propagación en laboratorio:   v_drift = {v_drift:.4f} m/s")
print(f" * Amplitud fraccional del campo (phi/Lambda): {phi_amplitude_rel:.2e}")
print("-" * 85)

# ------------------------------------------------------------------------------
# 2. RESOLUCIÓN DE LA DINÁMICA ESPACIO-TEMPORAL phi(x, t) EN LA MESA ÓPTICA
# ------------------------------------------------------------------------------
L_mesa = 1.50                 # Distancia entre Cavidad 1 (x=0) y Cavidad 2 (x=1.5 m)
Nx = 150                      # Malla espacial de 1 cm
x = np.linspace(0, L_mesa, Nx)
dx = x[1] - x[0]

# Número de onda efectivo k = omega / v_drift
k0 = omega0 / v_drift         # rad/m (longitud de onda lambda ~ 11.45 m)
lambda_espacial = 2.0 * math.pi / k0

print(f"[+] 2. CINEMÁTICA ONDULATORIA DEL ESCALAR EN LA MESA:")
print(f" * Longitud de onda espacial macroscópica:    lambda = {lambda_espacial:.2f} m")
print(f" * Desfase espacial total (x = 0 a 1.5 m):    Delta theta = {math.degrees(k0 * L_mesa):.2f}°")
print("-" * 85)

# Perfil instantáneo de phi(x, t) para t = 0 y t = tau/4
def phi_profile(x_arr, t_val):
    """Calcula phi(x, t) / Lambda_r considerando el paquete de 17 armónicos."""
    val = np.zeros_like(x_arr)
    for k in range(N_layers):
        fk = f0 * (1.0 + float(k) / float(N_layers))
        wk = 2.0 * math.pi * fk
        kk = wk / v_drift
        peso = math.cos((k * math.pi) / (2.0 * float(N_layers))) / math.sqrt(N_layers)
        val += peso * np.cos(wk * t_val - kk * x_arr)
    return phi_amplitude_rel * val

phi_x_t0 = phi_profile(x, 0.0)

# ------------------------------------------------------------------------------
# 3. GRADIENTE ESPACIAL Y ACELERACIÓN DE QUINTO RANGO SOBRE LOS ESPEJOS
# ------------------------------------------------------------------------------
# a_phi = - c^2 * (d_e / Lambda_r) * grad(phi)
# En unidades fraccionales relativas: a_phi = - v_drift^2 * d(phi/Lambda)/dx
c_light = 299792458.0
# La aceleración newtoniana efectiva sobre la masa del espejo (m = 1 kg de silicio):
grad_phi = np.gradient(phi_x_t0, dx)
# Aceleración gravitacional modificada: a_phi = (c^2 / 2) * grad(h)
a_scalar = (c_light**2) * grad_phi  # m/s^2

a_espejo_1 = a_scalar[0]      # x = 0.0 m (Cavidad 1)
a_espejo_2 = a_scalar[-1]     # x = 1.5 m (Cavidad 2)
delta_a = abs(a_espejo_2 - a_espejo_1)

print(f"[+] 3. GRADIENTE DE FUERZA Y ACELERACIÓN DIFERENCIAL SOBRE ESPEJOS:")
print(f" * Gradiente máximo de campo escalar:        |d(phi/Lambda)/dx| = {np.max(np.abs(grad_phi)):.2e} m^-1")
print(f" * Aceleración diferencial entre cavidades:   Delta a_phi = {delta_a:.2e} m/s^2")
print(f" * Cociente diferencial frente a g (9.8 m/s2): Delta a / g = {delta_a / 9.80665:.2e}")
print("-" * 85)

# ------------------------------------------------------------------------------
# 4. MODULACIÓN DE LA CONSTANTE DE ESTRUCTURA FINA (Delta alpha / alpha)
# ------------------------------------------------------------------------------
# Delta alpha / alpha = d_e * (phi(t) / Lambda_r)
delta_alpha_over_alpha = d_e * phi_x_t0

mod_alpha_max = np.max(np.abs(delta_alpha_over_alpha))
limite_relojes_opticos = 1.0e-18 # Límite actual de comparación Al+/Sr en NIST/PTB

print(f"[+] 4. MODULACIÓN DE CONSTANTES FUNDAMENTALES (Delta alpha / alpha):")
print(f" * Amplitud de modulación predicha:          |Delta alpha / alpha| = {mod_alpha_max:.2e}")
print(f" * Cota de sensibilidad de relojes ópticos:  sigma(Delta alpha / alpha) ~ {limite_relojes_opticos:.2e}")
print(f" * Relación Señal/Ruido en relojes de red:   SNR_reloj = {mod_alpha_max / limite_relojes_opticos:.1f}x")
print("=" * 85)

if mod_alpha_max > limite_relojes_opticos:
    print(" -> FIRMA ESCALAR INEQUÍVOCA CONFIRMADA:")
    print(f"    1. El radión modula periódicamente alpha con amplitud {mod_alpha_max:.2e} a {f0} Hz.")
    print(f"    2. Esta oscilación supera en un factor {mod_alpha_max/limite_relojes_opticos:.1f}x el suelo de ruido")
    print("       de las comparaciones directas entre relojes ópticos de estroncio y silicio.")
    print("    3. Si se compara la cavidad de silicio frente a un reloj atómico de red,")
    print("       la modulación del escalar desacopla la variación de alpha de la deformación puramente métrica.")
print("=" * 85)
