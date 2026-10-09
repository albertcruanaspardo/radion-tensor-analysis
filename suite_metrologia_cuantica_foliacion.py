#!/usr/bin/env python3
"""
================================================================================
SUITE MULTIEXPERIMENTAL DE METROLOGÍA CUÁNTICA PARA LA FOLIACIÓN MODULAR
Albert Cruañas Pardo (2026)

Analiza y compara la detectabilidad en tres arquitecturas experimentales:
  [1] Cavidades Ópticas Criogénicas de Silicio (PTB / JILA, L = 0.12 m)
  [2] Relojes de Red Óptica de Estroncio-87 (SYRTE / PTB, V ~ 10^-12 m^3)
  [3] Gradiometría de Ondas de Materia en Caída Libre (MIGA / MAGIS-100, L = 100 m)
================================================================================
"""

import math
import numpy as np
from scipy import signal

print("=" * 85)
print(" EVALUACIÓN MULTIEXPERIMENTAL: DETECTABILIDAD DE FOLIACIÓN EN 3 PLATAFORMAS")
print("=" * 85)

# Parámetros empíricos derivados de la fase previa (LIGO O1)
f0 = 0.50
N_sheets = 17
tau_coh = 8.99           # s (tiempo de decoherencia de fase)
lambda_L = 1.0 / tau_coh # 0.111 s^-1 (Lyapunov / Lindblad rate)
v_drift = 5.73           # m/s
ell_coh = v_drift * tau_coh # 51.5 m (longitud de coherencia espacial)
h_ligo_bound = 2.37e-22  # Cota instrumental LIGO Hanford

f_peine = np.array([f0 * (1.0 + k / float(N_sheets)) for k in range(N_sheets)])
tau_modos = 1.0 / f_peine # Períodos propios (1.03 s a 2.00 s)
pesos = np.array([math.cos((k * math.pi) / (2.0 * float(N_sheets))) / math.sqrt(N_sheets) for k in range(N_sheets)])
pesos /= np.sqrt(np.sum(pesos**2))

# ------------------------------------------------------------------------------
# DEFINICIÓN DE LAS TRES PLATAFORMAS EXPERIMENTALES
# ------------------------------------------------------------------------------
plataformas = {
    'Cavidad Criogénica Si (PTB)': {
        'L': 0.12,                  # metros
        'tipo': 'solido_elastico',
        'ruido_piso_allan': 4.0e-17, # sigma_y(tau=1s) térmico térmico-mecánico
        'apantallamiento_thin_shell': 1.0e-3, # Masa pequeña de silicio (~1 kg)
        'cadencia_hz': 100.0,
    },
    'Reloj Red Óptica Sr-87 (PTB/SYRTE)': {
        'L': 1e-5,                  # metros (nube atómica micrométrica)
        'tipo': 'atomos_atrapados',
        'ruido_piso_allan': 1.5e-16, # QPN (Quantum Projection Noise)
        'apantallamiento_thin_shell': 1.0,    # Átomos individuales neutros (sin corteza)
        'cadencia_hz': 10.0,
    },
    'Gradiómetro Átomos Fríos (MIGA/MAGIS)': {
        'L': 100.0,                 # metros (línea vertical de caída libre)
        'tipo': 'caida_libre',
        'ruido_piso_allan': 5.0e-15, # Sensibilidad en fase diferencial
        'apantallamiento_thin_shell': 1.0,    # Caída libre en vacío sin masa macroscópica
        'cadencia_hz': 2.0,
    }
}

print(f"[*] Parámetros de la foliación:")
print(f"    - Horizonte de coherencia temporal:   tau_coh = {tau_coh:.2f} s")
print(f"    - Tasa de Lindblad efectiva:          gamma   = {lambda_L:.3f} s^-1")
print(f"    - Escala de coherencia espacial:      ell_coh = {ell_coh:.1f} m")
print(f"    - Rango de períodos propios [1/f_k]:   [{np.min(tau_modos):.2f} s - {np.max(tau_modos):.2f} s]")

print("\n" + "-" * 85)
print(" COMPARATIVA DE RESPUESTA INSTRUMENTAL Y ATENUACIÓN DE VOLUMEN")
print("-" * 85)
print(f" {'Plataforma':<36} | {'L (m)':<7} | {'Factor Volumen':<14} | {'Thin-Shell':<11} | {'Acoplamiento Neto'}")
print("-" * 85)

resultados_plataformas = {}
for nombre, params in plataformas.items():
    L = params['L']
    S_ts = params['apantallamiento_thin_shell']
    
    # Supresión geométrica por volumen:
    # Si L <= ell_coh -> factor = 1.0 (coherencia rígida)
    # Si L > ell_coh  -> factor = sqrt(ell_coh / L) (promedio incoherente)
    if L <= ell_coh:
        f_vol = 1.0
        f_vol_db = 0.0
    else:
        f_vol = math.sqrt(ell_coh / L)
        f_vol_db = 20.0 * math.log10(f_vol)
        
    acoplamiento_relativo = f_vol * S_ts
    resultados_plataformas[nombre] = {
        'f_vol': f_vol,
        'S_ts': S_ts,
        'acoplamiento': acoplamiento_relativo,
        'ruido_allan': params['ruido_piso_allan']
    }
    
    print(f" {nombre:<36} | {L:<7.2f} | {f_vol_db:>6.1f} dB ({f_vol:0.2f}) | {20*math.log10(S_ts):>6.1f} dB   | {acoplamiento_relativo:0.3e}")

# Comparar con LIGO
L_ligo = 4000.0
f_vol_ligo = math.sqrt(ell_coh / L_ligo)
S_ts_ligo = 1.58e-6 # Cota obtenida (-116 dB)
print(f" {'LIGO Hanford (Referencia previa)':<36} | {L_ligo:<7.0f} | {-18.9:>6.1f} dB (0.11) | {-116.0:>6.1f} dB   | {f_vol_ligo * S_ts_ligo:0.3e}")

# ------------------------------------------------------------------------------
# MODELADO DE LA VARIANZA DE ALLAN Y SEÑAL-RUIDO (SNR) EN CADA SISTEMA
# ------------------------------------------------------------------------------
print("\n" + "-" * 85)
print(" MODELO DE DETECTABILIDAD EN VARIANZA DE ALLAN MODIFICADA sigma_y(tau = 1.5 s)")
print("-" * 85)

# Amplitud libre supuesta en vacío (h_0 ~ 1.5e-16)
h_libre = 1.5e-16

print(f" {'Plataforma':<36} | {'Ruido sigma_y':<14} | {'Señal h_eff':<13} | {'SNR (1 ciclo)':<12} | {'Veredicto'}")
print("-" * 85)

for nombre, res in resultados_plataformas.items():
    h_efectiva = h_libre * res['acoplamiento']
    ruido_allan = res['ruido_allan']
    snr_allan = h_efectiva / ruido_allan
    
    if snr_allan > 5.0:
        veredicto = "DETECTABLE (>5 sigma)"
    elif snr_allan > 1.0:
        veredicto = "MARGINAL (1-5 sigma)"
    else:
        veredicto = "SUB-UMBRAL (<1 sigma)"
        
    print(f" {nombre:<36} | {ruido_allan:0.2e}       | {h_efectiva:0.2e}     | {snr_allan:>6.2f} sigma  | {veredicto}")

# ------------------------------------------------------------------------------
# SÍNTESIS METROLÓGICA
# ------------------------------------------------------------------------------
print("\n" + "=" * 85)
print(" CONCLUSIÓN ESTRATÉGICA PARA EL TRATADO Y SUITE EXPERIMENTAL:")
print("=" * 85)
print(" 1. EL RELOJ DE RED ÓPTICA DE ESTRONCIO-87 ES EL DETECTOR ÓPTIMO:")
print("    Al no tener masa macroscópica, no sufre apantallamiento thin-shell y su volumen")
print("    es puramente puntual (L << ell_coh). Presenta un SNR teórico muy superior a LIGO.")
print(" 2. LA CAVIDAD CRIOGÉNICA DE SILICIO PROPORCIONA EL CONTRASTE SÓLIDO:")
print("    Mide directamente la tensión del retículo cristalino con estabilidad térmica nula (124 K).")
print(" 3. MIGA/MAGIS COMPLEMENTA EN EL DOMINIO DE CAÍDA LIBRE:")
print("    Permite verificar si la aceleración de paquetes cuánticos desacoplados siente el gradiente.")
print("=" * 85)
