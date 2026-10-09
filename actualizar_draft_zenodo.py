#!/usr/bin/env python3
"""
Actualización automática de metadatos en Zenodo (Draft ID: 23260503)
Albert Cruañas Pardo (2026)
"""
import requests
import json

DEPOSITION_ID = "23260503"
ACCESS_TOKEN = "ReFZCQKpjmnxSPCv1oFZwlDUjEKHgix0rbH9rOE4ApfZmtX4lb2CmPrErTKj"
BASE_URL = f"https://zenodo.org/api/deposit/depositions/{DEPOSITION_ID}"

headers = {"Content-Type": "application/json"}
params = {"access_token": ACCESS_TOKEN}

print(f"[*] Conectando con Zenodo para el depósito {DEPOSITION_ID}...")

# 1. Obtener estado actual
r_get = requests.get(BASE_URL, params=params)
if r_get.status_code != 200:
    print(f"[-] Error al acceder al borrador: {r_get.status_code} {r_get.text}")
    exit(1)

dep_data = r_get.json()
print(f"[+] Borrador localizado. Título actual: '{dep_data.get('metadata', {}).get('title', 'Sin título')}'")

# 2. Metadatos limpios y validados para el esquema de Zenodo
metadata_update = {
    "upload_type": "publication",
    "publication_type": "preprint",
    "publication_date": "2026-10-09",
    "title": "Non-Linear Dynamics, Causal Hysteresis, and Experimental Falsification of the 17-Layer Modular Foliation: From Israel Junction Conditions to Lagrangian Closure and Quantum Resonator Metrology",
    "creators": [
        {
            "name": "Cruañas Pardo, Albert",
            "affiliation": "Investigador Independiente, Amer, Girona, Cataluña, España",
            "orcid": "0009-0007-0846-2321"
        }
    ],
    "description": """<p>This treatise formalizes the quantum, tensorial, microscopic, and non-linear dynamics of the 17-layer modular foliation theory coupled to a collective radion scalar field &phi;(x<sup>&mu;</sup>). The stabilization of the compactification radius in an AdS<sub>5</sub> bulk is resolved via generalized Israel junction conditions, establishing an average surface tension &langle;&sigma;<sub>k</sub>&rangle; &approx; 2.70 &times; 10<sup>14</sup> GeV<sup>4</sup>.</p>

<p>We introduce an effective 4D Lagrangian density based on inter-sheet relative deformation &chi;<sub>n</sub>(x<sup>&mu;</sup>), a bistable double-well potential with barrier collapse governed by the accumulated causal exposure functional &Epsilon;<sub>n</sub>, and an intrinsic relativity parameter &real; = T<sup>&mu;</sup><sub>&mu;</sub> / &radic;(T<sub>&alpha;&beta;</sub> T<sup>&alpha;&beta;</sup>) that analytically discriminates between closed hysteresis (macroscopic baryonic density) and open hysteresis (radiative flux mediated by the derivative operator).</p>

<p>Particles are defined as localized phase transitions (condensed solitons of the foliated bulk), with half-lives rigorously derived via the Callan-Coleman Euclidean bounce action S<sub>E</sub>(A) &prop; &sigma;<sub>fol</sub><sup>4</sup> / [&Delta;V<sub>0</sub> - &kappa; A<sup>2</sup>]<sup>3</sup>. Finally, the four dynamic experimental pillars are analytically and numerically verified: intracavity chiral birefringence (+59.97 dB), anomalous metric Poisson ratio (&nu; &approx; 0.4201), distributed baseline sidereal triangulation (&plusmn;87.31 s across 500 m), and multichannel matched filtering reaching &rho;<sub>max</sub> = 312.36&sigma;.</p>""",
    "access_right": "open",
    "license": "cc-by-4.0",
    "keywords": [
        "Modular foliation",
        "Radion field",
        "Israel junction conditions",
        "Non-linear Lagrangian",
        "Bistable potential",
        "Causal hysteresis",
        "Fabry-Perot resonator",
        "Chiral birefringence",
        "Optical lattice clock",
        "Thin-shell screening",
        "Matched filtering"
    ],
    "version": "2.0.0",
    "language": "eng",
    "related_identifiers": [
        {
            "identifier": "10.5281/zenodo.22827063",
            "relation": "continues",
            "scheme": "doi"
        }
    ],
    "references": [
        "Cruañas Pardo, A. (2026). Topodynamic Field Theory and SO(8) Triality in Multi-Layer Foliations. Zenodo Preprints. DOI: 10.5281/zenodo.22827063.",
        "Cruañas Pardo, A. (2026). Scale Spectral Theory and Modular Operator Formalism for Boundary Manifolds. Submitted to Journal of Mathematical Physics, Ref. JMP26-AR-02867.",
        "Israel, W. (1966). Singular hypersurfaces and thin shells in general relativity. Nuovo Cimento B, 44, 434-449.",
        "Goldberger, W. D., & Wise, M. B. (1999). Modulus stabilization with bulk fields. Physical Review Letters, 83(24), 4922-4925.",
        "Connes, A. (1994). Noncommutative Geometry. Academic Press, San Diego.",
        "Kramer, M. et al. (2021). Strong-Field Gravity Tests with the Double Pulsar. Physical Review X, 11, 041050.",
        "Matei, D. G. et al. (2017). 1.5 μm Lasers with Sub-10 mHz Linewidth. Physical Review Letters, 118, 263202.",
        "Safronova, M. S. et al. (2018). Search for new physics with atoms and molecules. Reviews of Modern Physics, 90, 025008."
    ]
}

# 3. Enviar actualización
print("[*] Enviando nuevos metadatos a Zenodo...")
payload = {"metadata": metadata_update}
r_put = requests.put(BASE_URL, params=params, data=json.dumps(payload), headers=headers)

if r_put.status_code == 200:
    print("\n" + "=" * 70)
    print(" [+] ÉXITO: Los metadatos de https://zenodo.org/uploads/23260503")
    print("     han sido actualizados correctamente.")
    print("=" * 70)
    print("Recarga la página en tu navegador para ver los campos actualizados.")
else:
    print(f"[-] Error al guardar metadatos: {r_put.status_code}\n{r_put.text}")
