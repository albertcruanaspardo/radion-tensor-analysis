import numpy as np
from scipy import integrate

N_OCTETO = 8
LAMBDA_R = 3000.0   # GeV
L_BULK = 1.0
M_OBS = 602.50      # GeV

# Nodos de las 8 brechas del semiespacio conforme (z en [0, 0.5])
z_nodes = np.linspace(0.0, 0.5, N_OCTETO + 1)
dz_n = np.diff(z_nodes)

# Deformación exponencial y escalas de corte locales
k_warp = 2.5
A_n = k_warp * z_nodes[:-1]
lambda_local = LAMBDA_R * np.exp(-A_n)

# Tensión de Israel y pesos de degeneración Dynkin para SO(8) (suma = 16)
omega_n = 285.4 * np.exp(-0.5 * A_n)
pesos_dynkin = np.array([1, 1, 2, 2, 2, 2, 3, 3])

delta_m2_evanescente = []
delta_m2_escalar = []

for n in range(N_OCTETO):
    dz = dz_n[n]
    om = omega_n[n]
    cut = lambda_local[n]
    w = pesos_dynkin[n]

    # 1. Fase Evanescente (Negativa: túnel / Coleman-Weinberg)
    def integrand_ev(kappa):
        r = om / (2.0 * kappa + om)
        arg = np.maximum(1.0 - (r**2) * np.exp(-2.0 * kappa * dz), 1e-15)
        return (kappa**3) * np.log(arg)

    val_ev, _ = integrate.quad(integrand_ev, 0.5, cut)
    sigma_ev = - (w / (16.0 * np.pi**2)) * abs(val_ev) * (om / LAMBDA_R)
    delta_m2_evanescente.append(sigma_ev)

    # 2. Fase Escalar Conforme (Positiva: confinamiento geométrico de Israel / Ricci)
    sigma_sc = + (w * (om**2) / (LAMBDA_R**2)) * (cut**2) / (16.0 * np.pi**2) * 0.184
    delta_m2_escalar.append(sigma_sc)

delta_m2_evanescente = np.array(delta_m2_evanescente)
delta_m2_escalar = np.array(delta_m2_escalar)
balance_por_rango = delta_m2_escalar + delta_m2_evanescente

print("=" * 82)
print("SISTEMA DE DOBLE OCTETO: 8 FASES ESCALARES (+) vs 8 FASES EVANESCENTES (-)")
print("=" * 82)
print(f"{'Rango n':<8} | {'Fase Escalar (+)':<18} | {'Fase Evanescente (-)':<20} | {'Balance Neto (GeV^2)'}")
print("-" * 82)
for n in range(N_OCTETO):
    print(f"n = {n+1:<4} | +{delta_m2_escalar[n]:14.2f} GeV^2 | {delta_m2_evanescente[n]:16.2f} GeV^2 | {balance_por_rango[n]:+14.2f}")
print("-" * 82)
tot_sc = np.sum(delta_m2_escalar)
tot_ev = np.sum(delta_m2_evanescente)
tot_neto = np.sum(balance_por_rango)

print(f"{'TOTAL':<8} | +{tot_sc:14.2f} GeV^2 | {tot_ev:16.2f} GeV^2 | {tot_neto:+14.2f}")
print("=" * 82)

# Masa desnuda fundamental para cerrar exactamente el polo observable
m_bare = np.sqrt(M_OBS**2 - tot_neto)
print(f"[*] Polo Observable fijado:           {M_OBS:.2f} GeV")
print(f"[*] Desplazamiento neto de autoenergía: {tot_neto:+.2f} GeV^2")
print(f"[*] Masa desnuda fundamental deducida:  {m_bare:.2f} GeV")
print("=" * 82)
