#!/usr/bin/env python3
"""
Análisis Riguroso de la Densidad Espectral Cruzada (CSD) de Fermilab Holometer
Albert Cruañas Pardo (2026)
"""
import numpy as np
from scipy import signal
import math

print("=" * 85)
print(" ANÁLISIS ESPECTRAL RIGUROSO: FERMILAB HOLOMETER (E-990 PRL 2016)")
print("=" * 85)

# 1. Inspeccionar las primeras líneas del archivo
csv_file = "holo_data_PRL2016.csv"
with open(csv_file, 'r') as f:
    primeras_lineas = [f.readline().strip() for _ in range(5)]

print("[*] Primeras líneas del archivo CSV:")
for line in primeras_lineas:
    print(f"    {line}")

# Cargar datos
raw_data = np.genfromtxt(csv_file, delimiter=',', skip_header=1)
# Columnas en el CSV de PRL2016: Frecuencia (MHz), CSD_Real, CSD_Imag, Incertidumbre (o similar)
f_mhz = raw_data[:, 0]
csd_r = raw_data[:, 1]
csd_i = raw_data[:, 2] if raw_data.shape[1] > 2 else np.zeros_like(csd_r)
csd_mag = np.sqrt(csd_r**2 + csd_i**2)

N = len(f_mhz)
df_mhz = f_mhz[1] - f_mhz[0]
print(f"\n[+] Canales espectrales cargados: {N}")
print(f"    - Frecuencia mínima: {f_mhz[0]:.4f} MHz ({f_mhz[0]*1e6:.1f} Hz)")
print(f"    - Frecuencia máxima: {f_mhz[-1]:.4f} MHz ({f_mhz[-1]*1e6:.1f} Hz)")
print(f"    - Resolución espectral: df = {df_mhz*1e3:.3f} kHz")

# 2. Análisis del pico de autocorrelación (0.7388)
# Restamos la tendencia suave para evitar el sesgo de baja frecuencia
csd_detrend = csd_r - signal.savgol_filter(csd_r, 101, 3)
autocorr = signal.correlate(csd_detrend, csd_detrend, mode='full')
lags = signal.correlation_lags(N, N) * df_mhz  # Desplazamiento en MHz

center = len(autocorr) // 2
autocorr_norm = autocorr / autocorr[center]

# Buscar los picos prominentes de periodicidad en el espectro
# Excluimos el lóbulo central (+- 0.1 MHz)
mask_picos = (np.abs(lags) >= 0.1) & (lags > 0)
lags_pos = lags[mask_picos]
ac_pos = autocorr_norm[mask_picos]

picos_idx, props = signal.find_peaks(ac_pos, height=0.2, distance=int(0.05 / df_mhz))

print("\n" + "-" * 85)
print(" PRINCIPALES PERIODICIDADES DETECTADAS EN EL ESPECTRO DE FERMILAB:")
print("-" * 85)
print("  #  | Separación Delta f (MHz) |  Correlación Normalizada  | Relación con f0 = 0.50 Hz")
print("-----+--------------------------+---------------------------+---------------------------")

for i, idx in enumerate(picos_idx[:10]):
    delta_f = lags_pos[idx]
    val_corr = ac_pos[idx]
    m_entero = (delta_f * 1e6) / 0.50
    print(f" {i+1:02d}  |       {delta_f:0.4f} MHz        |          {val_corr:0.4f}           | m ~ {m_entero:0.1f} * f0")

# 3. Test Específico para la Foliación de 17 Capas en Escala Logarítmica
# Evaluamos si existe periodicidad fraccional logarítmica: delta(ln f) = ln(1 + 1/17) = 0.057158
log_f = np.log(f_mhz)
log_uniforme = np.linspace(log_f[0], log_f[-1], N)
csd_log = np.interp(log_uniforme, log_f, csd_detrend)

ac_log = signal.correlate(csd_log, csd_log, mode='full')
lags_log = np.linspace(-(log_uniforme[-1] - log_uniforme[0]), log_uniforme[-1] - log_uniforme[0], len(ac_log))
ac_log_norm = ac_log / ac_log[len(ac_log)//2]

delta_teorico = math.log(1.0 + 1.0 / 17.0) # ~ 0.05716
idx_17 = np.argmin(np.abs(lags_log - delta_teorico))
corr_17 = ac_log_norm[idx_17]

# Fondo circundante
mask_ruido = (np.abs(lags_log - delta_teorico) > 0.01) & (np.abs(lags_log - delta_teorico) < 0.04)
sigma_fondo = np.std(ac_log_norm[mask_ruido])
significancia_17 = corr_17 / sigma_fondo if sigma_fondo > 0 else 0.0

print("\n" + "=" * 85)
print(" TEST DE LA RAZÓN MODULAR (1/17) EN LA TORRE DE KALUZA-KLEIN:")
print("=" * 85)
print(f" * Paso modular teórico:           delta(ln f) = {delta_teorico:.6f}")
print(f" * Correlación en el paso 1/17:     C(1/17)     = {corr_17:+.4f}")
print(f" * Ruido estocástico del fondo:    sigma       = {sigma_fondo:.4f}")
print(f" * Significancia del modo 1/17:    Z           = {significancia_17:+.2f} sigma")
print("=" * 85)
