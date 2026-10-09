#!/usr/bin/env python3
"""
Empaquetado y subida directa de scripts de validación a Zenodo
Draft ID: 23261127 (Versión 2.0.0)
Albert Cruañas Pardo (2026)
"""
import os
import tarfile
import requests

DEPOSITION_ID = "23261127"
ACCESS_TOKEN = "ReFZCQKpjmnxSPCv1oFZwlDUjEKHgix0rbH9rOE4ApfZmtX4lb2CmPrErTKj"
ARCHIVE_NAME = "treatise_v2_validation_suite.tar.gz"

ARCHIVOS = [
    "calibrar_y_probar_tratado.py",
    "probar_cuatro_pilares_tratado.py",
    "analisis_multiexperimental_conclusiones.py",
    "falsacion_datos_abiertos.py"
]

params = {"access_token": ACCESS_TOKEN}

# 1. Empaquetar los 4 scripts
print(f"[*] Creando paquete comprimido: {ARCHIVE_NAME}...")
with tarfile.open(ARCHIVE_NAME, "w:gz") as tar:
    for archivo in ARCHIVOS:
        if os.path.exists(archivo):
            tar.add(archivo)
            print(f"  [+] Añadido al paquete: {archivo}")
        else:
            print(f"  [-] Aviso: No se encontró el archivo: {archivo}")

if not os.path.exists(ARCHIVE_NAME):
    print("[-] Error crítico: No se pudo generar el archivo .tar.gz")
    exit(1)

# 2. Obtener URL del bucket del borrador
print(f"[*] Conectando con Zenodo para obtener el bucket del borrador {DEPOSITION_ID}...")
url_dep = f"https://zenodo.org/api/deposit/depositions/{DEPOSITION_ID}"

r = requests.get(url_dep, params=params)
if r.status_code != 200:
    print(f"[-] Error al consultar el borrador: {r.status_code}\n{r.text}")
    exit(1)

bucket_url = r.json().get("links", {}).get("bucket")
if not bucket_url:
    print("[-] No se localizó la URL del bucket en la respuesta de la API.")
    exit(1)

print(f"[+] Bucket localizado: {bucket_url}")

# 3. Subir el archivo al bucket del borrador
print(f"[*] Subiendo {ARCHIVE_NAME} al borrador de Zenodo...")
upload_url = f"{bucket_url}/{ARCHIVE_NAME}"

with open(ARCHIVE_NAME, "rb") as fp:
    r_upload = requests.put(upload_url, data=fp, params=params)

if r_upload.status_code in [200, 201]:
    print("\n" + "=" * 75)
    print(" [+] ÉXITO: El paquete de scripts ha sido subido a Zenodo.")
    print(f"     Archivo:  {ARCHIVE_NAME}")
    print(f"     Borrador: https://zenodo.org/uploads/{DEPOSITION_ID}")
    print("=" * 75)
    print("Al recargar la página verás el archivo en la lista de 'Files'.")
    print("Solo te quedará adjuntar el PDF compilado desde Overleaf y pulsar 'Publish'.")
else:
    print(f"[-] Error durante la subida: {r_upload.status_code}\n{r_upload.text}")
