#!/usr/bin/env python3
"""
Modelado del Factor Thin-Shell de la Foliacion (N=17) en la Tierra
Albert Cruañas Pardo (2026)
"""
import numpy as np
import matplotlib.pyplot as plt

# Constantes fisicas
R_EARTH = 6.371e6        # m
M_EARTH = 5.972e24       # kg
G_NEWTON = 6.67430e-11   # m^3 / (kg s^2)
C_LIGHT = 2.99792e8      # m/s
M_PL_RED = 2.435e18      # GeV (M_Pl reducida)
GEV_A_JOULE = 1.60218e-10

# Potencial gravitatorio superficial de la Tierra
phi_N = (G_NEWTON * M_EARTH) / (R_EARTH * C_LIGHT**2)  # ~ 6.95e-10

print("=" * 80)
print(" MODELADO DEL FACTOR THIN-SHELL: TIERRA EN EL HALO GALACTICO")
print("=" * 80)
print(f"[*] Potencial gravitatorio superficial terrestre: Phi_N = {phi_N:.3e}")

# Parametros del modelo de 17 hojas
xi_17 = 1.0e-5
beta_eff = 1.0           # Acoplamiento conformal estandar
phi_out_adim = 1.0e-15   # Phi_out / M_Pl en el halo exterior

# Calculo analitico de thin-shell
theta_screen = phi_out_adim / (6.0 * beta_eff * phi_N)
delta_R = theta_screen * R_EARTH

print(f"[*] Escala de vacio en el halo (Phi_out / M_Pl): {phi_out_adim:.2e}")
print(f"[*] Factor de acoplamiento conformal beta    : {beta_eff:.2f}")
print(f"\n[+] RESULTADOS DEL APANTALLAMIENTO:")
print(f"    - Factor Thin-Shell (Delta R / R_oplus)  : theta = {theta_screen:.3e}")
print(f"    - Grosor de la capa activa superficial    : Delta R = {delta_R*100.0:.2f} cm")
print(f"    - Factor de supresion de fuerza quinta   : F_phi / F_N = {6.0 * (beta_eff**2) * theta_screen:.3e}")

# Perfil radial simulado
r_grid = np.linspace(0.8 * R_EARTH, 1.4 * R_EARTH, 1000)
r_norm = r_grid / R_EARTH

# Perfil analitico del campo normalizado (Phi - Phi_in) / (Phi_out - Phi_in)
phi_perfil = np.zeros_like(r_grid)

# En el interior (r < R_earth - Delta R): satura a 0
# En la cascara (R_earth - Delta R <= r <= R_earth): crece parabolicamente
# En el exterior (r > R_earth): cae asintoticamente como 1 - (1 - theta) * (R/r)
idx_int = r_grid < (R_EARTH - delta_R)
idx_shell = (r_grid >= (R_EARTH - delta_R)) & (r_grid <= R_EARTH)
idx_ext = r_grid > R_EARTH

phi_perfil[idx_int] = 0.0
# Transicion en la capa
r_s = r_grid[idx_shell]
phi_perfil[idx_shell] = ((r_s - (R_EARTH - delta_R)) / delta_R)**2 * (delta_R / R_EARTH)
# Exterior
phi_perfil[idx_ext] = 1.0 - (1.0 - theta_screen) * (R_EARTH / r_grid[idx_ext])

# Graficar
plt.rcParams.update({"font.family": "serif", "font.size": 10, "mathtext.fontset": "cm"})
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.5))

ax1.plot(r_norm, phi_perfil, color="#2c3e50", linewidth=2.0, label=r"Perfil de Campo $\Phi(r) / \Phi_{\mathrm{out}}$")
ax1.axvline(1.0, color="#e74c3c", linestyle="--", label="Superficie Terrestre ($R_\\oplus$)")
ax1.axvspan(1.0 - theta_screen, 1.0, color="#f39c12", alpha=0.3, label=rf"Thin Shell ($\Delta R \approx {delta_R*100:.1f}\text{ cm}$)")
ax1.set_xlim(0.95, 1.15)
ax1.set_xlabel(r"Radio Normalizado $r / R_\oplus$")
ax1.set_ylabel(r"Campo Normalizado $\Phi(r) / \Phi_{\mathrm{out}}$")
ax1.set_title("Apantallamiento del Campo en la Corteza Terrestre")
ax1.grid(True, linestyle=":", alpha=0.6)
ax1.legend(loc="upper left", fontsize=8)

# Comparativa de la fuerza quinta frente a la gravedad newtoniana
r_ext_norm = np.linspace(1.0, 3.0, 500)
fuerza_supresion = 2.0 * (beta_eff**2) * (3.0 * theta_screen) * (1.0 / r_ext_norm**2)
ax2.plot(r_ext_norm, fuerza_supresion, color="#2980b9", linewidth=2.0, label=r"Fuerza Quinta $F_\Phi / F_{\mathrm{Newton}}$")
ax2.axhline(0, color="gray", linestyle=":")
ax2.set_yscale("log")
ax2.set_xlabel(r"Distancia Orbital $r / R_\oplus$")
ax2.set_ylabel(r"$F_\Phi / F_{\mathrm{Newton}}$")
ax2.set_title("Supresión Dinámica de la Interacción Fuera de la Tierra")
ax2.grid(True, linestyle=":", alpha=0.6)
ax2.legend(loc="upper right", fontsize=8)

plt.tight_layout()
plt.savefig("perfil_apantallamiento_tierra.png", dpi=300)
plt.savefig("perfil_apantallamiento_tierra.pdf")
print("\n[+] Gráficos generados: perfil_apantallamiento_tierra.png y .pdf")
