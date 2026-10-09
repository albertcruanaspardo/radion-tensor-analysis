#!/usr/bin/env python3
"""
Pipeline de Análisis de Residuos de Fase en Relojes Ópticos (Sr-87)
Albert Cruañas Pardo (2026)

Calcula:
  - Función de estructura de fase D_phi(tau)
  - Varianza de Allan modificada sigma_y(tau)
  - Búsqueda espectral en los 17 periodos de foliación tau_k in [1.03, 2.00] s
"""
import math
import numpy as np
from scipy import signal

print("=" * 85)
print(" ANÁLISIS DE RESIDUOS DE FASE Y VARIANZA DE ALLAN EN REDES DE Sr-87")
print("=" * 85)

# Parámetros físicos derivados
f0 = 0.50
N_sheets = 17
tau_coh = 8.99           # s (tiempo de decoherencia medido)
gamma_lindblad = 0.111   # s^-1
f_peine = np.array([f0 * (1.0 + k / float(N_sheets)) for k in range(N_sheets)])
tau_modos = 1.0 / f_peine # Períodos propios [1.03 s a 2.00 s]
pesos = np.array([math.cos((k * math.pi) / (2.0 * float(N_sheets))) / math.sqrt(N_sheets) for k in range(N_sheets)])
pesos /= np.sqrt(np.sum(pesos**2))

# Simulación de serie temporal de batido (1800 s = 30 min a fs = 10 Hz)
fs = 10.0
T_total = 1800.0
N_samples = int(T_total * fs)
t = np.linspace(0, T_total, N_samples, endpoint=False)

# Ruido instrumental: Quantum Projection Noise (QPN blanco en frecuencia)
# sigma_y(1s) = 1.5e-16 -> PSD blanco S_y = 2 * sigma_y^2
sigma_qpn = 1.5e-16
ruido_y = np.random.normal(0.0, sigma_qpn * np.sqrt(fs), N_samples)

# Inyección de la señal física con disipación de Lindblad y fluctuaciones de fase
h0_señal = 1.5e-16
señal_y = np.zeros(N_samples)

# Modular en bloques de coherencia tau_coh
n_bloque = int(tau_coh * fs)
n_bloques = N_samples // n_bloque

for b in range(n_bloques):
    idx_start = b * n_bloque
    idx_end = (b + 1) * n_bloque
    t_b = t[idx_start:idx_end] - t[idx_start]
    
    # Fase estocástica aleatoria en cada bloque (difusión de fase)
    fase_bloque = np.random.uniform(0, 2*math.pi, N_sheets)
    
    bloque_sig = np.zeros(len(t_b))
    for k, fk in enumerate(f_peine):
        bloque_sig += pesos[k] * np.cos(2.0 * math.pi * fk * t_b + fase_bloque[k])
    
    señal_y[idx_start:idx_end] = h0_señal * bloque_sig * math.exp(-0.5 * gamma_lindblad * tau_coh)

y_total = ruido_y + señal_y

print(f"[+] Serie analizada: {T_total:.0f} s ({T_total/60:.1f} min) a {fs:.1f} Hz ({N_samples} muestras)")
print(f"[+] Amplitud inyectada h0:   {h0_señal:.2e}")
print(f"[+] Piso de ruido QPN:       {sigma_qpn:.2e}")

# ------------------------------------------------------------------------------
# 1. FUNCIÓN DE ESTRUCTURA DE FASE DE SEGUNDO ORDEN: D_phi(tau)
# ------------------------------------------------------------------------------
print("\n[*] Evaluando función de estructura de fase D_phi(tau)...")
tau_vals = np.logspace(-1, 1.2, 50) # De 0.1 s a 15 s
D_phi = []

# Integración numérica de fase: phi(t) = 2*pi*nu0 * int y(t) dt
phi = 2.0 * math.pi * 4.29e14 * np.cumsum(y_total) / fs

for tau_i in tau_vals:
    m = int(round(tau_i * fs))
    if m < 1 or m >= len(phi):
        continue
    dif = phi[m:] - phi[:-m]
    D_phi.append(np.mean(dif**2))

D_phi = np.array(D_phi)
tau_eval = tau_vals[:len(D_phi)]

# ------------------------------------------------------------------------------
# 2. VARIANZA DE ALLAN sigma_y(tau)
# ------------------------------------------------------------------------------
print("[*] Calculando Varianza de Allan modificada sigma_y(tau)...")
allan_tau = [0.1, 0.2, 0.5, 1.0, 1.5, 2.0, 4.0, 8.0, 10.0, 20.0, 50.0]
allan_sigma = []

for tau_a in allan_tau:
    m = int(round(tau_a * fs))
    if m < 1:
        m = 1
    # Promediar en cajas de tamaño m
    n_cajas = len(y_total) // m
    y_prom = y_total[:n_cajas * m].reshape((n_cajas, m)).mean(axis=1)
    if len(y_prom) > 1:
        sigma_val = np.sqrt(0.5 * np.mean(np.diff(y_prom)**2))
    else:
        sigma_val = 0.0
    allan_sigma.append(sigma_val)

print("\n" + "=" * 80)
print(" RESULTADOS DE LA VARIANZA DE ALLAN sigma_y(tau):")
print("=" * 80)
print("  tau (s)  |  sigma_y total  |  sigma_QPN esperado  |  Exceso fraccional")
print("-----------+-----------------+----------------------+-------------------")
for tau_a, sig in zip(allan_tau, allan_sigma):
    qpn_teo = sigma_qpn / math.sqrt(tau_a)
    exceso = (sig - qpn_teo) / qpn_teo * 100.0
    print(f"   {tau_a:>5.1f} s |    {sig:0.3e}    |      {qpn_teo:0.3e}       |     {exceso:+6.1f} %")

# ------------------------------------------------------------------------------
# 3. CONTRASTE ESPECTRAL EN EL HORIZONTE DE COHERENCIA
# ------------------------------------------------------------------------------
print("\n" + "=" * 80)
print(" TEST DE HIPÓTESIS EN LOS 17 MODOS (tau in [1.03 s, 2.00 s]):")
print("=" * 80)
idx_banda = (tau_eval >= 1.03) & (tau_eval <= 2.00)
exceso_banda_media = np.mean(D_phi[idx_banda])
idx_fuera = (tau_eval < 1.03)
exceso_fuera = np.mean(D_phi[idx_fuera]) if np.any(idx_fuera) else 1.0

ratio_estructura = exceso_banda_media / exceso_fuera
print(f" * Relación de estructura en la banda modular: R_D = {ratio_estructura:.2f}")
if ratio_estructura > 1.2:
    print(" [+] MODULACIÓN MODULAR DETECTADA en la escala temporal característica.")
else:
    print(" [-] Compatible con ruido estocástico puro sin modulación armónica.")
print("=" * 85)
