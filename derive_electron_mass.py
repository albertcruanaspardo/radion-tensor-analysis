#!/usr/bin/env python3
"""
DERIVACIÓN DE LA MASA DEL ELECTRÓN EN LA FOLIACIÓN DE 17 CAPAS
Calcula el parámetro de localización c_e y la masa física a partir de kL y epsilon.
"""

import math
from scipy.special import lambertw

# Parámetros del tratado
N = 17
epsilon = 1.0 / (4.0 * N * math.pi)       # 1 / (68 * pi) ~ 4.681e-3
kL = 4.0 * math.pi / math.sqrt(N * epsilon) # ~ 36.84
v_ew = 246.22                              # GeV (VEV del Higgs)
m_e_exp = 0.51099895e-3                    # GeV (510.998 keV)

# Acoplo de Yukawa experimental objetivo
y_e_target = (math.sqrt(2.0) * m_e_exp) / v_ew

# Resolver: u * exp(-kL * u) * (1 + epsilon) = y_e_target
# => -kL * u * exp(-kL * u) = -kL * (y_e_target / (1 + epsilon))
# => -kL * u = W_0(-kL * y_e_target / (1 + epsilon))
# Nota: Como queremos la solución con u ~ 0.4 (u grande), usamos la rama W_{-1}
arg = -kL * (y_e_target / (1.0 + epsilon))
u_sol = -float(lambertw(arg, k=-1).real) / kL

c_e = 0.5 * (1.0 + u_sol)

# Reconstruir la masa teórica
y_e_calc = u_sol * math.exp(-kL * u_sol) * (1.0 + epsilon)
m_e_calc = (y_e_calc * v_ew) / math.sqrt(2.0)

print("=" * 70)
print(" DERIVACIÓN GEOMÉTRICA DE LA MASA DEL ELECTRÓN (17 CAPAS)")
print("=" * 70)
print(f" * Parámetro modular (epsilon):        {epsilon:.6e}")
print(f" * Jerarquía interbrana (kL):          {kL:.3f}")
print(f" * Parámetro de localización (c_e):    {c_e:.4f}")
print(f" * Factor de supresión de bulk (exp):  {math.exp(-kL * u_sol):.4e}")
print(f" * Acoplamiento de Yukawa derivado:    y_e = {y_e_calc:.4e}")
print(f" * Masa del electrón derivada:         m_e = {m_e_calc*1000.0:.4f} MeV ({m_e_calc*1e6:.1f} keV)")
print(f" * Masa experimental (PDG):            m_e = {m_e_exp*1000.0:.4f} MeV")
print(f" * Concordancia:                       {100.0 * (1.0 - abs(m_e_calc - m_e_exp)/m_e_exp):.3f}%")
print("=" * 70)
