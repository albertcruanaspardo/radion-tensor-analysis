#!/usr/bin/env python3
"""
Generador de Gráficos de Exclusión Infrasónica (LIGO O3b 72h)
Albert Cruañas Pardo (2026)
Lee exclusion_infrasonica_gwosc_72h.json y grafica h0(f) y SNR(f)
frente a los umbrales estadísticos de detección y ruido instrumental.
"""
import json
import os
import sys
import numpy as np
import matplotlib.pyplot as plt

INPUT_FILE = "exclusion_infrasonica_gwosc_72h.json"
OUTPUT_PNG = "exclusion_infrasonica_gwosc_72h.png"
OUTPUT_PDF = "exclusion_infrasonica_gwosc_72h.pdf"

if not os.path.exists(INPUT_FILE):
    print(f"[!] Error: No se encontró el archivo de datos {INPUT_FILE}")
    sys.exit(1)

with open(INPUT_FILE, "r") as f:
    data = json.load(f)

espectro = data["espectro_exclusion"]
duracion_h = data.get("duracion_horas", 71.68)
detector = data.get("detector", "H1")
rayleigh_res = float(data.get("rayleigh_hz", 3.875e-6))

f0_vals = np.array([p["f0_carrier_hz"] for p in espectro])
rho_vals = np.array([p["snr_obs_rho"] for p in espectro])
h0_vals = np.array([float(p["h0_upper_limit_95cl"]) for p in espectro])

# Curva sintética de sensibilidad instrumental en deformación (ASD de LIGO O3 integrada a 72h)
# ASD ~ 1e-19 / sqrt(Hz) a 1 Hz, con pendiente ~ 1/f^2 por aislamiento pendular
f_dense = np.linspace(f0_vals.min() - 0.05, f0_vals.max() + 0.05, 300)
T_obs_sec = duracion_h * 3600.0
# h0_sens ~ 1.645 * sqrt(Sn(f) / T_obs)
asd_ligo_infrasonic = 1.0e-19 * ((1.0 / f_dense)**2)
h0_expected_curve = 1.645 * (asd_ligo_infrasonic / np.sqrt(T_obs_sec))

# Configuración visual
plt.rcParams.update({
    "font.family": "serif",
    "font.size": 11,
    "axes.labelsize": 12,
    "axes.titlesize": 12,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "legend.fontsize": 10,
    "figure.titlesize": 13,
    "mathtext.fontset": "cm"
})

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(9, 8), sharex=True, gridspec_kw={"height_ratios": [1.2, 1.0]})

# -------------------------------------------------------------------------
# PANEL 1: Curva de Exclusión h0(f) vs Ruido Instrumental
# -------------------------------------------------------------------------
ax1.plot(f_dense, h0_expected_curve, color="#7f8c8d", linestyle="--", linewidth=1.5,
         label=r"Sensibilidad esperada LIGO O3 (Piso térmico/sísmico $1/f^2$)")

ax1.step(f0_vals, h0_vals, where="mid", color="#c0392b", linewidth=2.0,
         label=r"Cota superior experimental $h_0$ (95% CL, 72h)")
ax1.scatter(f0_vals, h0_vals, color="#c0392b", s=30, zorder=5)

# Región excluida sombreada
ax1.fill_between(f0_vals, h0_vals, 1.0e-13, step="mid", color="#e74c3c", alpha=0.15,
                 label="Región excluida (Foliación no apantallada)")

# Resaltar portadora teórica de referencia (0.50 Hz)
ax1.axvline(0.50, color="#2980b9", linestyle=":", linewidth=1.5,
            label=r"Modo fundamental de foliación ($f_0 = 0.50\text{ Hz}$)")

ax1.set_ylabel(r"Deformación métrica $h_0$")
ax1.set_yscale("log")
ax1.set_ylim(bottom=min(h0_expected_curve.min(), h0_vals.min()) * 0.5, top=2.0e-14)
ax1.grid(True, which="both", linestyle=":", alpha=0.6)
ax1.legend(loc="upper right", framealpha=0.9)
ax1.set_title(f"Límites Superiores de Amplitud e Infrasonido — LIGO {detector} (O3b, $T = {duracion_h:.1f}\\text{{ h}}$)")

# -------------------------------------------------------------------------
# PANEL 2: Estadístico de Correlación rho(f0) vs Umbrales
# -------------------------------------------------------------------------
ax2.plot(f0_vals, rho_vals, marker="o", markersize=5, color="#2c3e50", linewidth=1.8,
         label=r"Significancia observada $\rho(f_0)$")

# Umbrales estadísticos
ax2.axhline(1.0, color="#95a5a6", linestyle="-.", linewidth=1.0, label=r"Nivel de ruido medio ($1\sigma$)")
ax2.axhline(3.0, color="#f39c12", linestyle="--", linewidth=1.5, label=r"Umbral de evidencia ($3\sigma$)")
ax2.axhline(5.0, color="#d35400", linestyle="-", linewidth=1.5, label=r"Umbral de descubrimiento ($5\sigma$)")
ax2.axvline(0.50, color="#2980b9", linestyle=":", linewidth=1.5)

# Zona compatible con hipótesis nula (ruido puro)
ax2.fill_between([f0_vals.min() - 0.05, f0_vals.max() + 0.05], 0, 3.0, color="#2ecc71", alpha=0.08,
                 label="Compatible con ruido puro / Apantallamiento")

ax2.set_xlabel(r"Frecuencia fundamental del peine $f_0\ \mathrm{[Hz]}$")
ax2.set_ylabel(r"Significancia SNR $\rho\ [\sigma]$")
ax2.set_xlim(f0_vals.min() - 0.02, f0_vals.max() + 0.02)
ax2.set_ylim(0, 5.5)
ax2.grid(True, linestyle=":", alpha=0.6)
ax2.legend(loc="upper right", framealpha=0.9, ncol=2)

# Cuadro informativo de parámetros
info_text = (
    f"Resolución de Rayleigh: $\\Delta f = {rayleigh_res:.2e}\\text{{ Hz}}$\n"
    r"Espaciado armónicos: $\Delta f_c = 0.0294\text{ Hz}$" "\n"
    f"Muestras procesadas: {data.get('muestras_4hz', 0):,}\n"
    r"Apantallamiento bariónico: $\mathcal{R} \to 1$"
)
ax2.text(0.03, 0.08, info_text, transform=ax2.transAxes,
         bbox=dict(boxstyle="round,pad=0.5", facecolor="#ecf0f1", edgecolor="#bdc3c7", alpha=0.9),
         fontsize=8.5)

plt.tight_layout()
plt.savefig(OUTPUT_PNG, dpi=300)
plt.savefig(OUTPUT_PDF)
print(f"[+] Gráficos exportados con éxito:")
print(f"  * Imagen rasterizada: {OUTPUT_PNG} (300 DPI)")
print(f"  * Documento vectorial: {OUTPUT_PDF}")
