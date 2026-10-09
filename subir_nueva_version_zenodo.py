#!/usr/bin/env python3
"""
Creación de nueva versión y subida de archivos/metadatos a Zenodo
Albert Cruañas Pardo (2026)
"""
import os
import tarfile
import requests
import json

BASE_RECORD_ID = "23247617"  # ID del registro publicado previo
ACCESS_TOKEN = "ReFZCQKpjmnxSPCv1oFZwlDUjEKHgix0rbH9rOE4ApfZmtX4lb2CmPrErTKj"
ARCHIVE_NAME = "treatise_lagrangian_validation_scripts.tar.gz"

ARCHIVOS = [
    "calibrar_y_probar_tratado.py",
    "probar_cuatro_pilares_tratado.py"
]

params = {"access_token": ACCESS_TOKEN}
headers = {"Content-Type": "application/json"}

# 1. Empaquetar scripts si no existen
if not os.path.exists(ARCHIVE_NAME):
    print(f"[*] Empaquetando scripts en {ARCHIVE_NAME}...")
    with tarfile.open(ARCHIVE_NAME, "w:gz") as tar:
        for archivo in ARCHIVOS:
            if os.path.exists(archivo):
                tar.add(archivo)
                print(f"  [+] Añadido: {archivo}")
            else:
                print(f"  [-] Archivo no encontrado: {archivo}")

# 2. Solicitar nueva versión del registro a la API
print(f"[*] Solicitando creación de nueva versión a partir de {BASE_RECORD_ID}...")
url_newversion = f"https://zenodo.org/api/deposit/depositions/{BASE_RECORD_ID}/actions/newversion"
r_nv = requests.post(url_newversion, params=params)

if r_nv.status_code != 201:
    print(f"[-] Error al crear nueva versión: {r_nv.status_code}\n{r_nv.text}")
    print("[*] Comprobando si el borrador ya existe...")
    # Si ya existía, usamos el borrador actual directamente
    draft_url = f"https://zenodo.org/api/deposit/depositions/23260503"
    r_check = requests.get(draft_url, params=params)
    if r_check.status_code == 200:
        draft_data = r_check.json()
        new_dep_id = "23260503"
        bucket_url = draft_data.get("links", {}).get("bucket")
    else:
        exit(1)
else:
    # Obtener el ID del nuevo borrador generado
    draft_url = r_nv.json().get("links", {}).get("latest_draft")
    r_draft = requests.get(draft_url, params=params)
    draft_data = r_draft.json()
    new_dep_id = draft_data.get("id")
    bucket_url = draft_data.get("links", {}).get("bucket")

print(f"[+] Nuevo borrador activo: ID {new_dep_id}")
print(f"[+] Bucket desbloqueado: {bucket_url}")

# 3. Subir el paquete tar.gz al bucket editable
print(f"[*] Subiendo {ARCHIVE_NAME} al bucket...")
upload_url = f"{bucket_url}/{ARCHIVE_NAME}"

with open(ARCHIVE_NAME, "rb") as fp:
    r_upload = requests.put(upload_url, data=fp, params=params)

if r_upload.status_code in [200, 201]:
    print(f"[+] Archivo {ARCHIVE_NAME} subido correctamente.")
else:
    print(f"[-] Error en subida de archivo: {r_upload.status_code}\n{r_upload.text}")
    exit(1)

# 4. Actualizar metadatos del nuevo borrador
print("[*] Aplicando metadatos actualizados...")
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

url_update = f"https://zenodo.org/api/deposit/depositions/{new_dep_id}"
r_put = requests.put(url_update, params=params, data=json.dumps({"metadata": metadata_update}), headers=headers)

if r_put.status_code == 200:
    print("\n" + "=" * 75)
    print(" [+] PROCESO COMPLETADO:")
    print(f"     Borrador: https://zenodo.org/uploads/{new_dep_id}")
    print(f"     Archivos cargados: {ARCHIVE_NAME}")
    print("=" * 75)
    print("Accede al enlace en tu navegador, adjunta el PDF compilado si lo deseas y pulsa 'Publish'.")
else:
    print(f"[-] Error al guardar metadatos: {r_put.status_code}\n{r_put.text}")
