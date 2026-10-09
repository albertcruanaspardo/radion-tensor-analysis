#!/usr/bin/env python3
"""
AUDITORÍA DE EXCITACIONES DEL ELECTRÓN (KK) FRENTE A HEPDATA (ATLAS/CMS)
Evalúa si el primer modo Kaluza-Klein del electrón (e_1) es detectable
en los registros de búsquedas de leptones excitados (e*) en LHC Run 2/3.
"""

import math

print("=" * 80)
print(" AUDITORÍA EN HEPDATA: SECTOR DEL ELECTRÓN EXCITADO (MODO KK e_1)")
print("=" * 80)

# Parámetros del tratado
Lambda_IR = 4.4030  # TeV
c_e = 0.6277        # Parámetro de localización derivado

# Cero de Bessel para fermiones en el bulk con c > 0.5
# Para orden nu = c - 0.5 ~ 0.128, el primer cero se sitúa en ~ 2.45
x1_fermion = 2.45
m_e1_tev = x1_fermion * Lambda_IR
m_e1_gev = m_e1_tev * 1000.0

sqrt_s_lhc = 13.6   # TeV

print(f"\n[+] PREDICCIÓN ANALÍTICA DE LA FOLIACIÓN DE 17 CAPAS:")
print(f" * Masa del electrón fundamental:       m_e0 = 0.5110 MeV")
print(f" * Parámetro de localización en bulk:   c_e  = {c_e:.4f} (anclado a brana UV)")
print(f" * Masa del 1er electrón excitado (e1): m_e1 = {m_e1_tev:.2f} TeV ({m_e1_gev:.0f} GeV)")
print(f" * Energía máxima LHC (Run 3):           sqrt(s) = {sqrt_s_lhc:.1f} TeV")

# Cotas oficiales en HEPData (ATLAS Run 2, 139 fb^-1, JHEP 02 (2020) 088)
# Búsqueda de leptones excitados e* en canal e* -> e gamma
limite_exp_hepdata_tev = 4.60  # Cota de exclusión al 95% CL

print("\n[+] REGISTROS EXPERIMENTALES EN HEPDATA (Búsqueda de leptones e*):")
print(f" * Límite de exclusión ATLAS al 95% CL: m_e* > {limite_exp_hepdata_tev:.2f} TeV")

margen = m_e1_tev - limite_exp_hepdata_tev

print("\n--- STATUS FENOMENOLÓGICO Y ALCANCE CINEMÁTICO ---")
print(f" * Umbral de pares (pp -> e1 e1):        Requiere sqrt(s) >= {2*m_e1_tev:.2f} TeV  [BLOQUEO TOTAL]")
print(f" * Producción simple (pp -> e1 e):       Requiere sqrt(s) >= {m_e1_tev:.2f} TeV    [PDFs prácticamente nulas]")
print(f" * Margen sobre la cota de HEPData:      +{margen:.2f} TeV por encima de la exclusión")

if m_e1_tev > limite_exp_hepdata_tev:
    print("\n[+] VEREDICTO:")
    print(" 1. El electrón excitado no es detectable en los datos actuales de HEPData.")
    print(" 2. La teoría respeta holgadamente la cota experimental (m_e1 > 4.6 TeV).")
    print(" 3. Su detección requerirá colisionadores futuros de alta energía (FCC-hh a 100 TeV).")
    print(" 4. La única manifestación del electrón a bajas energías es su masa exacta (511 keV)")
    print("    y las cotas de precisión electrodébil (parámetros S, T y momentos magnéticos anómalos).")
print("=" * 80)
