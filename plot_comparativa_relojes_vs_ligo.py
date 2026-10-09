#!/usr/bin/env python3
"""
Comparativa de Sensibilidad: Relojes Ópticos (Sin Apantallamiento) vs LIGO (Apantallado)
Albert Cruañas Pardo (2026)
"""
import json
import numpy as np
import matplotlib.pyplot as plt

# Cargar ambos JSON
with open("exclusion_infrasonica_gwosc_72h.json", "r") as f:
    data_ligo = json.load(f)

with open("exclusion_relojes_opticos_bipm.json", "r") as f:
    data_clocks = json.load(f)

f_ligo = np.array([p["f0_carrier_hz"] for p in data_ligo["espectro_exclusion"]])
h0_ligo = np.array([float(p["h0_upper_limit_95cl"]) for p in data_ligo["espectro_exclusion"]])
rho_ligo = np.array([p["snr_obs_rho"] for p in data_ligo["espectro_exclusion"]])

f_clk = np.array([p["f0_carrier_hz"] for p in data_clocks["espectro_relojes"]])
y_clk = np.array([float(p["frac_limit_95cl"]) for p in data_clocks["espectro_relojes"]])
rho_clk = np.array([p["snr_rho"] for p in data_clocks["espectro_relojes"]])

plt.rcParams.update({
    "font.family": "serif",
    "font.size": 11,
    "mathtext.fontset": "cm"
})

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(9, 8), sharex=True, gridspec_kw={"height_ratios": [1.2, 1.0]})

# PANEL 1: Cotas Experimentales
ax1.step(f_ligo, h0_ligo, where="mid", color="#c0392b", linewidth=2.0,
         label=r"LIGO O3b H1 ($h_0$, Apantallamiento bariónico $\mathcal{R} \to 1$)")
ax1.step(f_clk, y_clk, where="mid", color="#2980b9", linewidth=2.0,
         label=r"Red de Relojes Ópticos ($\Delta\nu/\nu_0$, Sin apantallamiento $\mathcal{R} \to 0$)")

ax1.axvline(0.50, color="#8e44ad", linestyle=":", linewidth=1.5, label=r"Modo fundamental ($0.50\text{ Hz}$)")

ax1.set_yscale("log")
ax1.set_ylabel(r"Amplitud fraccionaria / Deformación")
ax1.set_ylim(1.0e-20, 1.0e-13)
ax1.grid(True, which="both", linestyle=":", alpha=0.6)
ax1.legend(loc="upper right", framealpha=0.9)
ax1.set_title("Sensibilidad Cruzada: Supresión Macroscópica vs Sondas Cuánticas (72h)")

# PANEL 2: Estadístico de Señal-Ruido (SNR)
ax2.plot(f_ligo, rho_ligo, marker="s", markersize=5, color="#c0392b", linestyle="--", label=r"SNR LIGO H1")
ax2.plot(f_clk, rho_clk, marker="o", markersize=5, color="#2980b9", linestyle="-", label=r"SNR Red de Relojes")

ax2.axhline(3.0, color="#f39c12", linestyle="--", linewidth=1.2, label=r"Umbral $3\sigma$")
ax2.axvline(0.50, color="#8e44ad", linestyle=":", linewidth=1.5)
ax2.fill_between([0.38, 1.22], 0, 3.0, color="#2ecc71", alpha=0.08)

ax2.set_xlabel(r"Frecuencia fundamental $f_0\ \mathrm{[Hz]}$")
ax2.set_ylabel(r"Significancia $\rho\ [\sigma]$")
ax2.set_xlim(0.38, 1.22)
ax2.set_ylim(0, 4.0)
ax2.grid(True, linestyle=":", alpha=0.6)
ax2.legend(loc="upper right", framealpha=0.9, ncol=3)

plt.tight_layout()
plt.savefig("comparativa_relojes_vs_ligo_72h.png", dpi=300)
plt.savefig("comparativa_relojes_vs_ligo_72h.pdf")
print("[+] Gráfico comparativo generado: comparativa_relojes_vs_ligo_72h.png / .pdf")
