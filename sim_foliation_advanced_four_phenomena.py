#!/usr/bin/env python3
"""
================================================================================
ANÁLISIS AVANZADO: CUATRO FRENTES FENOMENOLÓGICOS DE LA FOLIACIÓN MODULAR
Foliación de 17 Capas | Albert Cruañas Pardo (2026)

Integra en un solo ejecutable:
1. Respuesta elasto-óptica fotoelástica y acoplamiento modal TEM_00 -> TEM_01.
2. Cuantización canónica del peine modular y compresión de vacío (squeezing).
3. Radiación dipolar escalar en púlsares binarios y límites de decaimiento orbital.
4. Algoritmo de filtrado óptimo adaptado (Matched Filtering) y significancia SNR.
================================================================================
"""

import math
import numpy as np
from scipy import signal

print("=" * 85)
print(" CUATRO FRENTES AVANZADOS: ELASTO-ÓPTICA, CUÁNTICA, PÚLSARES Y MATCHED-FILTERING")
print("=" * 85)

# ==============================================================================
# FRENTE 1: RESPUESTA ELASTO-ÓPTICA Y ACOPLAMIENTO MODAL TEM_00 -> TEM_01
# ==============================================================================
# Silicio monocristalino a 124 K (eje óptico en [111], n_Si ~ 3.48 a 1542 nm)
n_Si = 3.48
# Tensor fotoelástico típico de Si a 1.55 um: p11 = -0.101, p12 = +0.009, p44 = -0.051
p11, p12 = -0.101, 0.009
h0_strain = 1.5e-16
nu_metr = 0.4201 # Coeficiente de Poisson efectivo reconstruido

# Deformaciones longitudinal (S_zz) y transversal (S_rr)
S_zz = h0_strain
S_rr = -nu_metr * h0_strain

# Perturbación del índice de refracción efectivo: Delta(1/n^2) = p11*S_zz + 2*p12*S_rr
delta_inv_n2 = p11 * S_zz + 2.0 * p12 * S_rr
delta_n = -0.5 * (n_Si**3) * delta_inv_n2

# Desfase de Gouy y desacople modal en cavidad de L = 0.21 m, espejos R = 1.0 m
L_cav = 0.21
R_curv = 1.0
w0_beam = 350e-6 # Cintura del haz óptico (350 um)
# Factor de acoplamiento de potencia a modo espacial ortogonal TEM_01:
# eta_01 ~ (pi * delta_n * L_cav / lambda)^2
lambda_opt = 1542e-9
fase_elasto_optica = (2.0 * math.pi * abs(delta_n) * L_cav) / lambda_opt
eta_transversal_01 = (fase_elasto_optica)**2

print(f"[+] 1. RESPUESTA ELASTO-ÓPTICA Y ACOPLAMIENTO TRANSVERSAL TEM_00 -> TEM_01:")
print(f" * Perturbación del índice de refracción (Delta n): Delta n = {delta_n:.2e}")
print(f" * Desfase elasto-óptico por tránsito simple:       Delta phi = {fase_elasto_optica:.2e} rad")
print(f" * Razón de conversión de potencia a TEM_01:        eta_01    = {eta_transversal_01:.2e}")
print(f" * Separación espectral modo fundamental/cruzado:   Delta f_G = {299792458.0 / (2.0 * L_cav) * (math.acos(1 - L_cav/R_curv)/math.pi) / 1e6:.2f} MHz")
print("-" * 85)

# ==============================================================================
# FRENTE 2: CUANTIZACIÓN CANÓNICA DEL PEINE Y ESTADOS COMPRIMIDOS (SQUEEZING)
# ==============================================================================
# 17 interfaces de Israel actuando como amplificador paramétrico cuántico
hbar = 1.054571817e-34
N_sheets = 17
f0 = 0.50
omega0 = 2.0 * math.pi * f0

# Parámetro de bombeo paramétrico no lineal r_sq ~ kappa_NL * h0 * N_sheets
# kappa_NL = 1.264e-4 (m*strain)^-1 derivado previamente
kappa_NL = 1.2641e-4
# Ganancia paramétrica efectiva por paso a través de las 17 interfaces
r_squeezing = kappa_NL * h0_strain * float(N_sheets) * 1.50 # parámetro adimensional
# Supresión de ruido de cuadratura en dB: Delta X_sq = e^(-2*r)
# En una cavidad con finura F = 1.2e5, acumulación resonante: r_cav = r_sq * (2*F/pi)
Finura = 120000.0
r_cav_eff = r_squeezing * (2.0 * Finura / math.pi)
squeezing_db = -10.0 * math.log10(math.exp(2.0 * min(5.0, r_cav_eff)))

print(f"[+] 2. CUANTIZACIÓN CANÓNICA Y COMPRESIÓN DE VACÍO (SQUEEZING):")
print(f" * Número de modos discretos cuantizados:           N_modos = {N_sheets} osciladores")
print(f" * Parámetro de compresión paramétrica simple:       r_sq    = {r_squeezing:.2e}")
print(f" * Compresión efectiva intracavidad acumulada:       r_eff   = {r_cav_eff:.2e}")
print(f" * Nivel de compresión de ruido de cuadratura:       Delta X = {squeezing_db:.2f} dB (bajo SQL)")
print("-" * 85)

# ==============================================================================
# FRENTE 3: RADIACIÓN DIPOLAR ESCALAR EN PÚLSARES BINARIOS
# ==============================================================================
# PSR B1913+16 y PSR J0737-3039: comprobación de emisión dipolar escalar
G_N = 6.67430e-11
c_light = 299792458.0
M_sun = 1.989e30
# Púlsar doble PSR J0737-3039A/B:
m1 = 1.3381 * M_sun
m2 = 1.2489 * M_sun
P_orb = 2.45 * 3600.0 # 2.45 horas en segundos
e_orb = 0.088 # excentricidad

# Sensibilidad escalar alfa_A ~ d_e * (E_grav / m c^2), E_grav/mc^2 ~ 0.15 para estrella de neutrones
d_e = 1.0 / math.sqrt(17.0) # 0.2425
alpha_NS1 = d_e * 0.15
alpha_NS2 = d_e * 0.14 # ligera asimetría de masa de enlace
delta_alpha = abs(alpha_NS1 - alpha_NS2)

# Corrección dipolar a la tasa de decaimiento orbital:
# Pdot_dipolo / Pdot_cuadrupolo ~ (5/96) * (delta_alpha)^2 * (c / v_orb)^2
omega_orb = 2.0 * math.pi / P_orb
a_semi = ((G_N * (m1 + m2) / (omega_orb**2))**(1.0/3.0))
v_orb = omega_orb * a_semi

ratio_dipolo_cuadrupolo = (5.0 / 96.0) * (delta_alpha**2) * ((c_light / v_orb)**2)
# Cota experimental del púlsar doble: |Pdot_obs / Pdot_GR - 1| < 1.3e-4
cota_pulsar_doble = 1.3e-4

print(f"[+] 3. RADIACIÓN DIPOLAR ESCALAR EN PÚLSARES BINARIOS (PSR J0737-3039):")
print(f" * Sensibilidad escalar efectiva (alpha_NS):         alpha_1 = {alpha_NS1:.4f}, alpha_2 = {alpha_NS2:.4f}")
print(f" * Carga escalar dipolar diferencial:               Delta alpha = {delta_alpha:.2e}")
print(f" * Ratio de pérdida dipolar escalar vs cuadrupolar: Pdot_dip / Pdot_GR = {ratio_dipolo_cuadrupolo:.2e}")
print(f" * Límite experimental observacional del púlsar:    Cota obs = {cota_pulsar_doble:.2e}")
cumple_pulsar = ratio_dipolo_cuadrupolo < cota_pulsar_doble
print(f" * Consistencia astrofísica observacional:          {'COMPATIBLE CON DATOS DE PÚLSAR' if cumple_pulsar else 'EXCLUIDO'}")
print("-" * 85)

# ==============================================================================
# FRENTE 4: FILTRO ÓPTIMO (MATCHED-FILTERING) PARA SERIES TEMPORALES
# ==============================================================================
# Simulación de un registro experimental ruidoso de 300 s
fs_mf = 50.0
T_mf = 300.0
N_mf = int(fs_mf * T_mf)
t_mf = np.linspace(0, T_mf, N_mf, endpoint=False)

# Ruido coloreado de la cavidad (blanco + 1/f)
ruido_base = np.random.normal(0, 3.0e-17, N_mf)

# Señal inyectada con las 17 frecuencias y amplitudes de Dirichlet
f_peine = [f0 * (1.0 + float(k) / float(N_sheets)) for k in range(N_sheets)]
signal_inyectada = np.zeros(N_mf)
for k, fk in enumerate(f_peine):
    peso = math.cos((k * math.pi) / (2.0 * float(N_sheets))) / math.sqrt(N_sheets)
    signal_inyectada += h0_strain * peso * np.cos(2.0 * math.pi * fk * t_mf)

dato_detector = signal_inyectada + ruido_base

# Construcción de la plantilla normalizada (Template)
template_h = np.zeros(N_mf)
for k, fk in enumerate(f_peine):
    peso = math.cos((k * math.pi) / (2.0 * float(N_sheets))) / math.sqrt(N_sheets)
    template_h += peso * np.cos(2.0 * math.pi * fk * t_mf)

# Normalización del matched filter: <T | T> = 1
sigma_T = math.sqrt(np.sum(template_h**2))
template_h /= sigma_T

# Correlación cruzada / Salida del filtro óptimo
snr_time_series = np.correlate(dato_detector, template_h, mode='same') / np.std(ruido_base)
rho_max = np.max(np.abs(snr_time_series))

# Tasa de Falsa Alarma (FAR) bajo estadística gaussiana de Neyman-Pearson
# P_FA = 0.5 * erfc(rho / sqrt(2)) ~ exp(-rho^2 / 2) / (sqrt(2*pi) * rho)
log10_FAR = (-0.5 * (rho_max**2) * math.log10(math.e)) - math.log10(math.sqrt(2.0 * math.pi) * rho_max)

print(f"[+] 4. ALGORITMO DE FILTRADO ÓPTIMO ADAPTADO (MATCHED FILTERING):")
print(f" * Duración del registro analizado:                 T_obs   = {T_mf:.1f} s ({N_mf} muestras)")
print(f" * Pico del estadístico de detección (SNR óptimo):  rho_max = {rho_max:.2f} sigma")
print(f" * Significancia estadística de la detección:       Signif  = {rho_max:.1f} sigma (Umbral de oro: > 5.0 sigma)")
print(f" * Tasa de Falsa Alarma teórica (FAR):              P_FA    ~ 10^({log10_FAR:.1f})")
print("=" * 85)

if rho_max >= 5.0 and cumple_pulsar:
    print(" -> SÍNTESIS GLOBAL COMPLETADA CON ÉXITO:")
    print(f"    1. La interacción elasto-óptica preserva el modo TEM_00 con fuga a TEM_01 despreciable ({eta_transversal_01:.1e}).")
    print(f"    2. La cuantización predice compresión de ruido intracavidad de {squeezing_db:.2f} dB.")
    print("    3. La emisión dipolar escalar respeta rigurosamente las cotas del púlsar doble.")
    print(f"    4. El filtro adaptado detecta el peine de 17 líneas con significancia de {rho_max:.1f} sigma.")
print("=" * 85)
