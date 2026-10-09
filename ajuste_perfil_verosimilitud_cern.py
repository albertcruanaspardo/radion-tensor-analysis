#!/usr/bin/env python3
"""
Ajuste de Perfil de Verosimilitud (Profile Likelihood Scan)
y Test de Razón de Verosimilitudes (Wilks / Cowan et al. 2011)
Albert Cruañas Pardo (2026)

Determina la significancia local (z_local) y global (Look-Elsewhere Effect)
del exceso observado en UPC Pb-Pb en torno a 30.5 GeV.
"""
import json
import numpy as np
import scipy.optimize as opt
import scipy.stats as stats
import matplotlib.pyplot as plt

print("=" * 95)
print(" TEST DE HIPÓTESIS Y AJUSTE DE PERFIL DE VEROSIMILITUD (UPC Pb-Pb 5.02 TeV)")
print("=" * 95)

with open("exclusion_cern_opendata_upc_pbpb.json", "r") as f:
    data_cern = json.load(f)

m_inv = np.array([p["masa_inv_gev"] for p in data_cern["espectro_masa"]])
n_obs = np.array([p["eventos_observados"] for p in data_cern["espectro_masa"]], dtype=np.float64)

# -------------------------------------------------------------------------
# 1. MODELO H0: HIPÓTESIS NULA (SOLO FONDO CONTINUO EXPONENCIAL / POWER-LAW)
# -------------------------------------------------------------------------
def bkg_model(m, a, b):
    # Continuo QED decreciente suave
    return a * np.exp(-b * (m - 5.0))

def nll_null(params):
    a, b = params
    if a <= 0 or b <= 0:
        return 1e12
    mu = bkg_model(m_inv, a, b)
    # Poisson NLL: sum(mu - n*log(mu))
    return np.sum(mu - n_obs * np.log(np.maximum(mu, 1e-9)))

res_null = opt.minimize(nll_null, [15000.0, 0.16], bounds=[(100.0, 1e6), (0.01, 2.0)])
nll_h0 = res_null.fun
a_h0, b_h0 = res_null.x
bkg_fit_h0 = bkg_model(m_inv, a_h0, b_h0)

print(f"[*] Hipótesis H0 convergida: a = {a_h0:.2f}, b = {b_h0:.4f}")
print(f"    -Log(L) H0 = {nll_h0:.3f}")

# -------------------------------------------------------------------------
# 2. MODELO H1: FONDO + SEÑAL RESONANTE GAUSSIANA (RESOLUCIÓN DETECTOR)
# -------------------------------------------------------------------------
def sig_model(m, n_sig, m_res, sigma_res):
    return (n_sig / (sigma_res * np.sqrt(2 * np.pi))) * np.exp(-0.5 * ((m - m_res) / sigma_res)**2)

def nll_alt(params):
    a, b, n_sig, m_res = params
    sigma_res = 0.85  # Resolución típica calorimétrica en CMS/ALICE a 30 GeV (~3%)
    if a <= 0 or b <= 0 or n_sig < 0:
        return 1e12
    mu = bkg_model(m_inv, a, b) + sig_model(m_inv, n_sig, m_res, sigma_res)
    return np.sum(mu - n_obs * np.log(np.maximum(mu, 1e-9)))

res_alt = opt.minimize(nll_alt, [a_h0, b_h0, 50.0, 30.5], 
                       bounds=[(100.0, 1e6), (0.01, 2.0), (0.0, 5000.0), (28.0, 33.0)])
nll_h1 = res_alt.fun
a_h1, b_h1, n_sig_fit, m_res_fit = res_alt.x
sigma_res = 0.85
tot_fit_h1 = bkg_model(m_inv, a_h1, b_h1) + sig_model(m_inv, n_sig_fit, m_res_fit, sigma_res)
bkg_fit_h1 = bkg_model(m_inv, a_h1, b_h1)

# Estadístico de Wilks q0 = 2 * (NLL(H0) - NLL(H1))
q0 = 2.0 * max(0.0, nll_h0 - nll_h1)
z_local = np.sqrt(q0)

# Corrección de Look-Elsewhere Effect (Gross & Vitells 2010)
# N_trials ~ (M_max - M_min) / (2 * sigma_res) = (35 - 5) / 1.7 ~ 18 resoluciones independientes
n_trials = (35.0 - 5.0) / (2.0 * sigma_res)
p_local = stats.norm.sf(z_local)
p_global = min(1.0, p_local * n_trials)
z_global = stats.norm.isf(p_global) if p_global > 0 else 0.0

print(f"\n[*] Hipótesis H1 convergida:")
print(f"    Masa central resonante: M_0 = {m_res_fit:.3f} GeV")
print(f"    Eventos de señal:        N_sig = {n_sig_fit:.1f} +/- {np.sqrt(max(1.0, n_sig_fit)):.1f}")
print(f"    -Log(L) H1:             {nll_h1:.3f}")
print(f"    Estadístico Delta q_0:  {q0:.3f}")
print(f"    Significancia Local:    {z_local:.2f} sigma  (p-value = {p_local:.3e})")
print(f"    Significancia Global:   {z_global:.2f} sigma  (p-value = {p_global:.3e}, trials = {n_trials:.1f})")

# -------------------------------------------------------------------------
# 3. ESCANEO DE PERFIL EN MASA [26 a 34 GeV]
# -------------------------------------------------------------------------
m_scan = np.linspace(26.0, 34.0, 81)
q0_scan = []

for m_test in m_scan:
    def nll_fixed_mass(params):
        a, b, n_sig = params
        if a <= 0 or b <= 0 or n_sig < 0:
            return 1e12
        mu = bkg_model(m_inv, a, b) + sig_model(m_inv, n_sig, m_test, sigma_res)
        return np.sum(mu - n_obs * np.log(np.maximum(mu, 1e-9)))
    
    r_fix = opt.minimize(nll_fixed_mass, [a_h0, b_h0, 30.0], bounds=[(100.0, 1e6), (0.01, 2.0), (0.0, 2000.0)])
    q_val = 2.0 * max(0.0, nll_h0 - r_fix.fun)
    q0_scan.append(q_val)

q0_scan = np.array(q0_scan)
z_scan = np.sqrt(q0_scan)

# -------------------------------------------------------------------------
# 4. GENERACIÓN DE GRÁFICO TÉCNICO DE PERFIL DE VEROSIMILITUD
# -------------------------------------------------------------------------
plt.rcParams.update({"font.family": "serif", "font.size": 11, "mathtext.fontset": "cm"})
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(9, 8), sharex=True, gridspec_kw={"height_ratios": [1.2, 1.0]})

# Ajuste local en el rango del exceso
mask_zoom = (m_inv >= 26.0) & (m_inv <= 34.0)
ax1.errorbar(m_inv[mask_zoom], n_obs[mask_zoom], yerr=np.sqrt(n_obs[mask_zoom]), fmt='o',
             color='#2c3e50', capsize=3, label="Datos Exclusivos UPC Pb-Pb")
m_fine = np.linspace(26.0, 34.0, 200)
ax1.plot(m_fine, bkg_model(m_fine, a_h1, b_h1), color='#e74c3c', linestyle='--', label="Fondo continuo H0 (QED)")
ax1.plot(m_fine, bkg_model(m_fine, a_h1, b_h1) + sig_model(m_fine, n_sig_fit, m_res_fit, sigma_res),
         color='#27ae60', linewidth=2.0, label=f"Ajuste H1 (Fondo + Señal, M={m_res_fit:.1f} GeV)")
ax1.set_ylabel(r"Eventos / $1.0\ \mathrm{GeV}$")
ax1.grid(True, linestyle=":", alpha=0.6)
ax1.legend(loc="upper right")
ax1.set_title(r"Ajuste de Perfil de Verosimilitud en Exceso Bifotónico ($M \approx 30.5\ \mathrm{GeV}$)")

# Panel de significancia local
ax2.plot(m_scan, z_scan, color='#2980b9', linewidth=2.2, label=r"Significancia local $z_{\mathrm{local}}(M)$")
ax2.axhline(3.0, color='#f39c12', linestyle='--', label=r"Umbral de evidencia $3\sigma$")
ax2.axhline(z_global, color='#8e44ad', linestyle=':', label=f"Significancia global ({z_global:.2f}$\sigma$)")
ax2.fill_between(m_scan, 0, 3.0, color='#2ecc71', alpha=0.08)
ax2.set_xlabel(r"Masa invariante evaluada $M\ \mathrm{[GeV]}$")
ax2.set_ylabel(r"Significancia $\sqrt{q_0}\ [\sigma]$")
ax2.set_xlim(26.0, 34.0)
ax2.set_ylim(0, max(max(z_scan), 3.5) * 1.15)
ax2.grid(True, linestyle=":", alpha=0.6)
ax2.legend(loc="upper left")

plt.tight_layout()
plt.savefig("perfil_verosimilitud_upc_30gev.png", dpi=300)
plt.savefig("perfil_verosimilitud_upc_30gev.pdf")
print("[+] Gráfico de verosimilitud exportado: perfil_verosimilitud_upc_30gev.png / .pdf")

# Guardar informe
scan_summary = {
    "masa_optimizada_gev": round(float(m_res_fit), 3),
    "eventos_senal_ajustados": round(float(n_sig_fit), 2),
    "delta_q0": round(float(q0), 3),
    "significancia_local_sigma": round(float(z_local), 3),
    "significancia_global_sigma": round(float(z_global), 3),
    "p_value_local": f"{p_local:.4e}",
    "p_value_global": f"{p_global:.4e}"
}
with open("resultado_perfil_verosimilitud_30gev.json", "w") as f:
    json.dump(scan_summary, f, indent=2)
print("[+] Informe de verosimilitud guardado en: resultado_perfil_verosimilitud_30gev.json")
