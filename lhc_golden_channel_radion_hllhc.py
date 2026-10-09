#!/usr/bin/env python3
"""
CANAL DEFINITIVO DE BÚSQUEDA DEL RADIÓN: pp -> phi -> ZZ -> 4l (HL-LHC 3000 fb^-1)
Modelado en cuatro leptones (e+e-e+e-, mu+mu-mu+mu-, e+e-mu+mu-)
y canal complementario hh -> bb gammagamma.
"""
import math
import numpy as np

print("=" * 80)
print(" SIMULADOR DEL CANAL DORADO: RADIÓN A 602.5 GeV EN EL HL-LHC (3000 fb^-1)")
print("================================================================================")

m_phi = 602.50       # GeV
lumi_hl = 3000.0     # fb^-1 (14 TeV)

lambda_scenarios = [3.8, 5.0, 7.5, 10.0, 12.0]

br_phi_zz = 0.245
br_phi_hh = 0.178
br_zz_4l = 0.0045        # BR(Z -> ll)^2 para l = e, mu
br_hh_bbgaga = 0.0026

acc_eff_4l = 0.42
acc_eff_bbgaga = 0.18

# Fondo continuo en [592.5, 612.5] GeV
sigma_bkg_4l_fiducial = 0.028 * 20.0 * br_zz_4l
n_bkg_4l = sigma_bkg_4l_fiducial * lumi_hl * acc_eff_4l

sigma_bkg_bbgaga_fiducial = 0.085
n_bkg_bbgaga = sigma_bkg_bbgaga_fiducial * lumi_hl * acc_eff_bbgaga

print(f"Masa del radión:              {m_phi:.2f} GeV")
print(f"Luminosidad HL-LHC (14 TeV):   {lumi_hl:.0f} fb^-1")
print(f"Fondo esperado ZZ -> 4l:      B_4l     = {n_bkg_4l:.2f} eventos")
print(f"Fondo esperado hh -> bbgaga:  B_bbgaga = {n_bkg_bbgaga:.2f} eventos")

print("\n" + "=" * 80)
print(f"{'Escala Lambda':<15} | {'sigma(gg->phi)':<15} | {'S (4l)':<10} | {'Z_local (4l)':<14} | {'S (bbgaga)':<12} | {'Z_comb':<10}")
print("=" * 80)

sigma_ref_38 = 62.0  # fb a Lambda = 3.8 TeV

for l_val in lambda_scenarios:
    ratio_sq = (3.8 / l_val) ** 2
    sigma_prod = sigma_ref_38 * ratio_sq

    sigma_fid_4l = sigma_prod * br_phi_zz * br_zz_4l
    n_sig_4l = sigma_fid_4l * lumi_hl * acc_eff_4l
    z_4l = math.sqrt(2 * ((n_sig_4l + n_bkg_4l) * math.log(1 + n_sig_4l / n_bkg_4l) - n_sig_4l))

    sigma_fid_bbgaga = sigma_prod * br_phi_hh * br_hh_bbgaga
    n_sig_bbgaga = sigma_fid_bbgaga * lumi_hl * acc_eff_bbgaga
    z_bbgaga = math.sqrt(2 * ((n_sig_bbgaga + n_bkg_bbgaga) * math.log(1 + n_sig_bbgaga / n_bkg_bbgaga) - n_sig_bbgaga))

    z_comb = math.sqrt(z_4l**2 + z_bbgaga**2)
    flag = "[DESCUBRIMIENTO 5 sigma]" if z_comb >= 5.0 else ("[EVIDENCIA 3 sigma]" if z_comb >= 3.0 else "[Cota Exclusión]")

    print(f" {l_val:4.1f} TeV        | {sigma_prod:8.2f} fb      | {n_sig_4l:6.1f}     | {z_4l:6.2f} sigma    | {n_sig_bbgaga:6.1f}       | {z_comb:5.2f} sigma {flag}")

print("=" * 80)

# Búsqueda numérica exacta del umbral de 5 sigma
lambda_grid = np.linspace(3.0, 15.0, 1000)
lambda_5sigma = None
for l_test in lambda_grid:
    s_p = sigma_ref_38 * ((3.8 / l_test)**2)
    s4 = s_p * br_phi_zz * br_zz_4l * lumi_hl * acc_eff_4l
    sh = s_p * br_phi_hh * br_hh_bbgaga * lumi_hl * acc_eff_bbgaga
    z1 = math.sqrt(2 * ((s4 + n_bkg_4l) * math.log(1 + s4 / n_bkg_4l) - s4))
    z2 = math.sqrt(2 * ((sh + n_bkg_bbgaga) * math.log(1 + sh / n_bkg_bbgaga) - sh))
    zc = math.sqrt(z1**2 + z2**2)
    if zc <= 5.0:
        lambda_5sigma = l_test
        break

print("\n--- CONCLUSIONES FENOMENOLÓGICAS ---")
if lambda_5sigma:
    print(f"1. Umbral de descubrimiento formal (5 sigma) en HL-LHC (3000 fb^-1):")
    print(f"   Requiere una escala conforme: Lambda_phi <= {lambda_5sigma:.2f} TeV.")
    print(f"2. Para escalas con mayor confinamiento (Lambda_phi > {lambda_5sigma:.2f} TeV):")
    print(f"   La señal no alcanzará 5 sigma y se traducirá en cotas de exclusión al 95% CL.")
print("=" * 80)
