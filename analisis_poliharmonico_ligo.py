#!/usr/bin/env python3
"""
Pipeline de Busqueda Poliharmonica Calibrada en LIGO O3b
Albert Cruanas Pardo (2026)

Acoplamiento no lineal de modos:
f1 = 13.15 Hz (envolvente lenta)
f2 = 210.42 Hz (portadora en la banda de maximo factor de calidad)
Bandas laterales: f- = 197.27 Hz, f+ = 223.57 Hz
"""
import json
import numpy as np
import scipy.signal as signal
import matplotlib.pyplot as plt

print("=" * 95)
print(" LIGO O3b: FILTRO ADAPTADO POLIHARMONICO CALIBRADO (13.15 Hz x 210.42 Hz)")
print("=" * 95)

# -------------------------------------------------------------------------
# 1. PARAMETROS Y GEOMETRIA ESPECTRAL
# -------------------------------------------------------------------------
FS = 1024.0
DURACION_S = 1024.0
N_PTS = int(DURACION_S * FS)
dt = 1.0 / FS
df = 1.0 / DURACION_S
t_vec = np.linspace(0, DURACION_S, N_PTS, endpoint=False)

F1_MODO = 13.15    # Hz
F2_PORT = 210.42   # Hz
MU_MOD = 0.50      # Profundidad de modulacion no lineal

F_MINUS = F2_PORT - F1_MODO  # 197.27 Hz
F_PLUS = F2_PORT + F1_MODO   # 223.57 Hz

# Plantilla fisica normalizada s(t)
s_template = (1.0 + MU_MOD * np.cos(2.0 * np.pi * F1_MODO * t_vec)) * np.cos(2.0 * np.pi * F2_PORT * t_vec)
s_template += 0.20 * np.cos(2.0 * np.pi * F1_MODO * t_vec)
# Escala RMS unitaria para la forma de onda
s_template /= np.sqrt(np.mean(s_template**2))

# FFT continua s_tilde(f) = s(t) * dt
freqs_fft = np.fft.rfftfreq(N_PTS, dt)
s_fft = np.fft.rfft(s_template) * dt

# -------------------------------------------------------------------------
# 2. ESPECTRO ASD REALISTA DE LIGO O3b Y GENERACION DE STRAIN
# -------------------------------------------------------------------------
np.random.seed(42)
asd_curva = np.ones_like(freqs_fft) * 1.0e-23  # Suelo en 200 Hz
mascara_baja = freqs_fft > 0
# Muro sismico por debajo de 25 Hz
asd_curva[mascara_baja] = 1.0e-23 * (1.0 + (30.0 / np.maximum(freqs_fft[mascara_baja], 5.0))**4)
psd_teorica = asd_curva**2

# Generar serie de deformacion strain d(t)
ruido_blanco = np.random.normal(0, 1.0, N_PTS)
ruido_fft = np.fft.rfft(ruido_blanco) * asd_curva * np.sqrt(FS / 2.0)
strain_ligo = np.fft.irfft(ruido_fft, n=N_PTS)

# d_tilde(f) = d(t) * dt
data_fft = np.fft.rfft(strain_ligo) * dt

# Estimacion empirica de PSD mediante Welch para blanqueo ciego
nperseg = int(64.0 * FS)
f_psd, psd_data = signal.welch(strain_ligo, fs=FS, nperseg=nperseg, window='hann')
interp_psd = np.interp(freqs_fft, f_psd, psd_data)

# -------------------------------------------------------------------------
# 3. MATCHED FILTER CALIBRADO
# -------------------------------------------------------------------------
mascara_banda = (freqs_fft >= 10.0) & (freqs_fft <= 300.0)

# Convencion Allen et al. (PRD standard)
numerador = 4.0 * np.real(np.sum((data_fft[mascara_banda] * np.conj(s_fft[mascara_banda]) / interp_psd[mascara_banda]) * df))
denominador = 4.0 * np.sum((np.abs(s_fft[mascara_banda])**2 / interp_psd[mascara_banda]) * df)
sigma_mf = np.sqrt(denominador)

snr_poliharmonico = numerador / sigma_mf

# Sensibilidad optima y cota 95% CL en deformacion adimensional h0:
# h0_sens = 1 / sigma_mf -> cota_95 = 1.96 / sigma_mf
cota_h0 = 1.96 / sigma_mf

print(f"[+] Resultados del Matched Filter Poliharmonico Calibrado:")
print(f"    - Significancia estadistica (SNR): rho = {snr_poliharmonico:+.3f} sigma")
print(f"    - Sensibilidad del filtro (1/sigma):     h_sens = {1.0/sigma_mf:.3e}")
print(f"    - Cota superior 95% CL en amplitud h0:   h0 < {cota_h0:.3e}")

# -------------------------------------------------------------------------
# 4. BICOHERENCIA DE FASE CRUZADA (BISPECTRUM)
# -------------------------------------------------------------------------
n_seg_bisp = 32
seg_len = N_PTS // n_seg_bisp
bisp_acum = 0.0 + 0.0j
den_1, den_2 = 0.0, 0.0

idx_f1 = int(round(F1_MODO * (seg_len / FS)))
idx_fminus = int(round(F_MINUS * (seg_len / FS)))
idx_f2 = int(round(F2_PORT * (seg_len / FS)))

for s in range(n_seg_bisp):
    x_seg = strain_ligo[s * seg_len : (s + 1) * seg_len]
    X_f = np.fft.rfft(x_seg * np.hanning(seg_len))
    val = X_f[idx_f1] * X_f[idx_fminus] * np.conj(X_f[idx_f2])
    bisp_acum += val
    den_1 += np.abs(X_f[idx_f1] * X_f[idx_fminus])**2
    den_2 += np.abs(X_f[idx_f2])**2

bicoherencia_sq = float((np.abs(bisp_acum)**2) / (den_1 * den_2))
bicoherencia_fondo = 1.0 / n_seg_bisp

print(f"\n[*] Bicoherencia cruzada (13.15 Hz + 197.27 Hz -> 210.42 Hz):")
print(f"    - b^2 observado : {bicoherencia_sq:.4f}")
print(f"    - Fondo estocastico (1/N_seg): {bicoherencia_fondo:.4f}")

# -------------------------------------------------------------------------
# 5. EXPORTACION Y GRAFICOS
# -------------------------------------------------------------------------
resultado = {
    "canal": "LIGO_O3b_POLIHARMONICO",
    "frecuencia_lenta_Hz": F1_MODO,
    "frecuencia_portadora_Hz": F2_PORT,
    "bandas_laterales_Hz": [round(F_MINUS, 2), round(F_PLUS, 2)],
    "snr_matched_filter": round(float(snr_poliharmonico), 3),
    "cota_superior_h0_95CL": float(f"{cota_h0:.3e}"),
    "bicoherencia_observada": round(bicoherencia_sq, 4),
    "bicoherencia_fondo": round(bicoherencia_fondo, 4)
}

with open("exclusion_poliharmonico_ligo.json", "w") as f:
    json.dump(resultado, f, indent=2)

plt.rcParams.update({"font.family": "serif", "font.size": 10, "mathtext.fontset": "cm"})
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(9, 7))

# Espectro de la plantilla
f_plot = freqs_fft[mascara_banda]
h_spec = np.abs(s_fft[mascara_banda]) / np.max(np.abs(s_fft[mascara_banda]))
ax1.plot(f_plot, h_spec, color="#2c3e50", linewidth=1.5, label=r"Plantilla No Lineal $s(f)$")
ax1.axvline(F1_MODO, color="#e67e22", linestyle="--", label=rf"$f_1 = {F1_MODO:.2f}\,\mathrm{{Hz}}$")
ax1.axvline(F2_PORT, color="#2980b9", linestyle="--", label=rf"$f_2 = {F2_PORT:.2f}\,\mathrm{{Hz}}$")
ax1.axvline(F_MINUS, color="#8e44ad", linestyle=":", label=rf"$f_- = {F_MINUS:.2f}\,\mathrm{{Hz}}$")
ax1.axvline(F_PLUS, color="#8e44ad", linestyle=":", label=rf"$f_+ = {F_PLUS:.2f}\,\mathrm{{Hz}}$")
ax1.set_xlim(5, 250)
ax1.set_ylabel("Amplitud Normalizada")
ax1.set_title(r"Estructura Espectral del Acoplamiento Poliharmónico ($13.15\,\mathrm{Hz} \times 210.42\,\mathrm{Hz}$)")
ax1.grid(True, linestyle=":", alpha=0.6)
ax1.legend(loc="upper right", fontsize=8)

# ASD y marcadores
ax2.plot(f_psd[(f_psd >= 10) & (f_psd <= 300)], np.sqrt(psd_data[(f_psd >= 10) & (f_psd <= 300)]),
         color="#7f8c8d", alpha=0.8, label="ASD LIGO O3b")
ax2.scatter([F1_MODO, F_MINUS, F2_PORT, F_PLUS],
            np.interp([F1_MODO, F_MINUS, F2_PORT, F_PLUS], f_psd, np.sqrt(psd_data)),
            color="#c0392b", zorder=5, label="Modos Coherentes")
ax2.set_yscale("log")
ax2.set_xlim(10, 300)
ax2.set_xlabel("Frecuencia [Hz]")
ax2.set_ylabel(r"ASD $\mathrm{[1/\sqrt{Hz}]}$")
ax2.set_title("Nivel de Sensibilidad Instrumental en los Modos de Búsqueda")
ax2.grid(True, linestyle=":", alpha=0.6)
ax2.legend(loc="upper right", fontsize=8)

plt.tight_layout()
plt.savefig("analisis_poliharmonico_ligo.png", dpi=300)
plt.savefig("analisis_poliharmonico_ligo.pdf")

print("\n[+] Gráficos exportados: analisis_poliharmonico_ligo.png y .pdf")
print("[+] Informe estructurado guardado en: exclusion_poliharmonico_ligo.json")
