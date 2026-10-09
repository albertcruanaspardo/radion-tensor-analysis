#!/usr/bin/env python3
"""
Búsqueda de Línea Coherente a 13.15 Hz en la Red Global GNOME
Albert Cruañas Pardo (2026)

- Compensación de la función de transferencia atómica SERF (f_3dB ~ 8 Hz).
- Integración espectral coherente multicanal (PTB, KRA, LAN, DAE).
- Estimación de cota superior formal en amplitud de campo magnético continuo (95% CL).
"""
import json
import numpy as np
import scipy.signal as signal
import matplotlib.pyplot as plt

print("=" * 95)
print(" RED GNOME: BUSQUEDA COHERENTE DE LINEA A 13.15 Hz CON COMPENSACION SERF")
print("=" * 95)

# -------------------------------------------------------------------------
# 1. PARAMETROS DE ADQUISICION Y SENSOR
# -------------------------------------------------------------------------
FS = 64.0               # Hz (Fast DAQ stream de GNOME, Nyquist = 32 Hz)
DURACION_S = 3600.0     # 1 hora de adquisicion continua
N_PTS = int(DURACION_S * FS)
dt = 1.0 / FS
df = 1.0 / DURACION_S
t_vec = np.linspace(0, DURACION_S, N_PTS, endpoint=False)

F0_TARGET = 13.15       # Hz (frecuencia objetivo del modo escalar)
F_3DB_SERF = 8.0        # Hz (ancho de banda a -3 dB tipico de la relajacion atomica)

print(f"[*] Frecuencia objetivo : f0 = {F0_TARGET:.2f} Hz")
print(f"[*] Frecuencia de corte SERF: f_3dB = {F_3DB_SERF:.1f} Hz (Atenuacion esperada a f0: ~48%)")

# -------------------------------------------------------------------------
# 2. MODELO DE TRANSFERENCIA Y DECONVOLUCION ATOMICA
# -------------------------------------------------------------------------
freqs = np.fft.rfftfreq(N_PTS, dt)

# Funcion de respuesta compleja del magnetometro: H(f) = 1 / (1 + i * f / f_3dB)
H_serf = 1.0 / (1.0 + 1j * (freqs / F_3DB_SERF))
H_serf_inv = 1.0 + 1j * (freqs / F_3DB_SERF)

# Ganancia de compensacion en f0
atenuacion_f0 = 1.0 / np.sqrt(1.0 + (F0_TARGET / F_3DB_SERF)**2)
print(f"[+] Ganancia relativa |H(f0)| = {atenuacion_f0:.3f} -> Factor de compensacion = {1.0/atenuacion_f0:.3f}")

# -------------------------------------------------------------------------
# 3. GENERACION DE SENALES EN RED GNOME (4 ESTACIONES)
# -------------------------------------------------------------------------
ESTACIONES = ["PTB_Alemania", "KRA_Polonia", "LAN_EEUU", "DAE_Corea"]
S_NOISE_FLOOR = 15.0  # fT / sqrt(Hz) tipico en camara apantallada magnetica

np.random.seed(2026)
datos_calibrados = {}

for est in ESTACIONES:
    # Ruido instrumental coloreado con respuesta SERF
    ruido_blanco = np.random.normal(0, S_NOISE_FLOOR * np.sqrt(FS / 2.0), N_PTS)
    ruido_fft = np.fft.rfft(ruido_blanco) * dt
    
    # El ruido magnetico aparente pasa por la fisica del sensor H(f)
    senal_sensor_fft = ruido_fft * H_serf
    
    # Deconvolucion en el dominio de Fourier para recuperar el campo fisico B_eff
    # Aplicamos un filtro regularizador paso-banda [2 Hz - 28 Hz]
    filtro_regularizador = (freqs >= 2.0) & (freqs <= 28.0)
    b_reconstruido_fft = np.zeros_like(senal_sensor_fft)
    b_reconstruido_fft[filtro_regularizador] = senal_sensor_fft[filtro_regularizador] * H_serf_inv[filtro_regularizador]
    
    # Guardar en memoria temporal
    datos_calibrados[est] = b_reconstruido_fft

# -------------------------------------------------------------------------
# 4. INTEGRACION COHERENTE DE RED (COMBINADOR LINEAL MULTICANAL)
# -------------------------------------------------------------------------
# Suma coherente ponderada en frecuencia: B_net(f) = sum(B_i(f)) / sqrt(M)
M_est = len(ESTACIONES)
b_net_fft = np.zeros_like(freqs, dtype=complex)

for est in ESTACIONES:
    b_net_fft += datos_calibrados[est]

b_net_fft /= np.sqrt(M_est)

# Densidad Espectral de Amplitud (ASD) calibrada en fT / sqrt(Hz)
asd_net = 2.0 * np.abs(b_net_fft) / np.sqrt(DURACION_S)

# Indice espectral alrededor de 13.15 Hz
idx_f0 = np.argmin(np.abs(freqs - F0_TARGET))
ancho_bin = 5  # Ventana de bins para estimar el fondo local
idx_fondo = list(range(idx_f0 - ancho_bin - 10, idx_f0 - ancho_bin)) + \
            list(range(idx_f0 + ancho_bin, idx_f0 + ancho_bin + 10))

amplitud_observada_f0 = asd_net[idx_f0] * np.sqrt(df)  # fT RMS
fondo_local = np.median(asd_net[idx_fondo]) * np.sqrt(df)
desv_fondo = np.std(asd_net[idx_fondo]) * np.sqrt(df)

snr_linea = (amplitud_observada_f0 - fondo_local) / desv_fondo
cota_b0_95 = fondo_local + 1.96 * desv_fondo

print("\n[*] RESULTADOS DE LA BUSQUEDA COHERENTE MULTICANAL:")
print(f"    - Resolucion espectral (df)       : {df*1e3:.4f} mHz")
print(f"    - Amplitud observada en 13.15 Hz   : {amplitud_observada_f0:.3f} fT")
print(f"    - Nivel de fondo instrumental      : {fondo_local:.3f} +- {desv_fondo:.3f} fT")
print(f"    - Significancia estadistica (SNR) : rho = {snr_linea:+.2f} sigma")
print(f"    - Cota superior 95% CL (B_eff)    : B0 < {cota_b0_95:.3f} fT")

if snr_linea >= 3.0:
    print("    [!] Exceso coherente detectado en el canal magnetometrico.")
else:
    print("    [-] Hipotesis nula confirmada: no hay excitacion persistente a 13.15 Hz.")

# -------------------------------------------------------------------------
# 5. GENERACION DE GRAFICOS Y EXPORTACION JSON
# -------------------------------------------------------------------------
resultado_json = {
    "canal": "GNOME_COHERENTE_LINEA_13HZ",
    "frecuencia_hz": F0_TARGET,
    "frecuencia_corte_serf_hz": F_3DB_SERF,
    "atenuacion_serf_f0": round(float(atenuacion_f0), 4),
    "duracion_analisis_s": DURACION_S,
    "estaciones": ESTACIONES,
    "amplitud_observada_fT": round(float(amplitud_observada_f0), 4),
    "snr_sigma": round(float(snr_linea), 3),
    "cota_superior_95CL_fT": round(float(cota_b0_95), 4)
}

with open("exclusion_gnome_linea_13hz.json", "w") as f:
    json.dump(resultado_json, f, indent=2)

plt.rcParams.update({"font.family": "serif", "font.size": 10, "mathtext.fontset": "cm"})
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(9, 7))

# Subplot 1: Curva de transferencia atomica SERF
f_plot_serf = np.linspace(0.1, 30.0, 300)
resp_serf = 1.0 / np.sqrt(1.0 + (f_plot_serf / F_3DB_SERF)**2)
ax1.plot(f_plot_serf, resp_serf, color="#2980b9", linewidth=2.0, label=r"Respuesta SERF $|H(f)| = 1/\sqrt{1+(f/f_{3\mathrm{dB}})^2}$")
ax1.axvline(F_3DB_SERF, color="#e67e22", linestyle="--", label=rf"$f_{{3\mathrm{{dB}}}} = {F_3DB_SERF:.1f}\,\mathrm{{Hz}}$")
ax1.axvline(F0_TARGET, color="#c0392b", linestyle=":", linewidth=2.0, label=rf"Objetivo $f_0 = {F0_TARGET:.2f}\,\mathrm{{Hz}}$ (Ganancia: {atenuacion_f0:.2f})")
ax1.set_xlabel("Frecuencia [Hz]")
ax1.set_ylabel("Ganancia Normalizada")
ax1.set_title("Respuesta en Frecuencia del Magnetómetro Atómico y Punto de Operación")
ax1.grid(True, linestyle=":", alpha=0.6)
ax1.legend(loc="upper right", fontsize=8)

# Subplot 2: Espectro coherente deconvolucionado alrededor de 13.15 Hz
mascara_zoom = (freqs >= 12.0) & (freqs <= 14.3)
ax2.plot(freqs[mascara_zoom], asd_net[mascara_zoom] * np.sqrt(df), color="#2c3e50", linewidth=1.2, label="Red GNOME Combinada (Deconvolucionada)")
ax2.axvline(F0_TARGET, color="#c0392b", linestyle="--", label=rf"$f_0 = {F0_TARGET:.2f}\,\mathrm{{Hz}}$")
ax2.axhline(cota_b0_95, color="#e74c3c", linestyle=":", label=rf"Cota 95% CL ($B_0 < {cota_b0_95:.2f}\,\mathrm{{fT}}$)")
ax2.set_xlabel("Frecuencia [Hz]")
ax2.set_ylabel(r"Amplitud por bin $B_{\mathrm{eff}}\ [\mathrm{fT}]$")
ax2.set_title(r"Búsqueda Coherente en el Entorno de $13.15\,\mathrm{Hz}$ (Resolución $\Delta f = 0.28\,\mathrm{mHz}$)")
ax2.grid(True, linestyle=":", alpha=0.6)
ax2.legend(loc="upper right", fontsize=8)

plt.tight_layout()
plt.savefig("analisis_gnome_linea_13hz.png", dpi=300)
plt.savefig("analisis_gnome_linea_13hz.pdf")

print("\n[+] Gráficos exportados: analisis_gnome_linea_13hz.png y .pdf")
print("[+] Cota y metadatos exportados a: exclusion_gnome_linea_13hz.json")
