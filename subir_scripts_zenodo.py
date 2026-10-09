#!/usr/bin/env python3
"""
Empaquetado y subida directa de scripts a Zenodo (Draft ID: 23260503)
Albert Cruañas Pardo (2026)
"""
import os
import tarfile
import requests

DEPOSITION_ID = "23260503"
ACCESS_TOKEN = "ReFZCQKpjmnxSPCv1oFZwlDUjEKHgix0rbH9rOE4ApfZmtX4lb2CmPrErTKj"
ARCHIVE_NAME = "treatise_lagrangian_validation_scripts.tar.gz"

ARCHIVOS = [
    "calibrar_y_probar_tratado.py",
    "probar_cuatro_pilares_tratado.py"
]

# 1. Crear el paquete tar.gz con los scripts
print(f"[*] Empaquetando scripts en {ARCHIVE_NAME}...")
with tarfile.open(ARCHIVE_NAME, "w:gz") as tar:
    for archivo in ARCHIVOS:
        if os.path.exists(archivo):
            tar.add(archivo)
            print(f"  [+] Añadido: {archivo}")
        else:
            print(f"  [-] Archivo no encontrado: {archivo}")

if not os.path.exists(ARCHIVE_NAME):
    print("[-] Error: no se pudo crear el archivo comprimido.")
    exit(1)

# 2. Obtener la URL del bucket del depósito
print(f"[*] Obteniendo enlace del bucket para el depósito {DEPOSITION_ID}...")
url_dep = f"https://zenodo.org/api/deposit/depositions/{DEPOSITION_ID}"
params = {"access_token": ACCESS_TOKEN}

r = requests.get(url_dep, params=params)
if r.status_code != 200:
    print(f"[-] Error al consultar depósito: {r.status_code} {r.text}")
    exit(1)

bucket_url = r.json().get("links", {}).get("bucket")
if not bucket_url:
    print("[-] No se encontró la URL del bucket en la respuesta de Zenodo.")
    exit(1)

print(f"[+] Bucket localizado: {bucket_url}")

# 3. Subir el archivo al bucket
print(f"[*] Subiendo {ARCHIVE_NAME} a Zenodo...")
upload_url = f"{bucket_url}/{ARCHIVE_NAME}"

with open(ARCHIVE_NAME, "rb") as fp:
    r_upload = requests.put(upload_url, data=fp, params=params)

if r_upload.status_code in [200, 201]:
    print("\n" + "=" * 70)
    print(" [+] ÉXITO: Los scripts han sido subidos a Zenodo.")
    print(f"     Archivo: {ARCHIVE_NAME}")
    print(f"     Verifica en: https://zenodo.org/uploads/{DEPOSITION_ID}")
    print("=" * 70)
else:
    print(f"[-] Error en la subida: {r_upload.status_code}\n{r_upload.text}")
