#!/usr/bin/env python3
import requests
import json

DEPOSITION_ID = "23247617"
ACCESS_TOKEN = "ReFZCQKpjmnxSPCv1oFZwlDUjEKHgix0rbH9rOE4ApfZmtX4lb2CmPrErTKj"

BASE_URL = f"https://zenodo.org/api/deposit/depositions/{DEPOSITION_ID}"
headers = {"Content-Type": "application/json"}
params = {"access_token": ACCESS_TOKEN}

print(f"[*] Conectando con el registro {DEPOSITION_ID} en Zenodo...")

# Paso 1: Abrir registro para edición
r_edit = requests.post(f"{BASE_URL}/actions/edit", params=params)
if r_edit.status_code in [201, 200]:
    print("[+] Registro abierto para edición (estado borrador activo).")
elif r_edit.status_code == 403:
    print("[*] El registro ya se encuentra en estado borrador/edición.")
else:
    print(f"[-] Aviso al abrir edición: {r_edit.status_code} {r_edit.text}")

# Paso 2: Obtener metadatos actuales
r_get = requests.get(BASE_URL, params=params)
if r_get.status_code != 200:
    print(f"[-] Error al obtener metadatos: {r_get.status_code} {r_get.text}")
    exit(1)

record_data = r_get.json()
metadata = record_data.get("metadata", {})

# Descripción técnica con la reformulación del lagrangiano
descripcion_adicional = """
<hr/>
<h3>Actualización Teórica: Reformulación No Lineal del Lagrangiano Efectivo y Cierre Dinámico de Histéresis</h3>

<p>Se formaliza la densidad lagrangiana efectiva 4D acoplada a la foliación modular de 17 capas de Israel, incorporando de manera autoconsistente:</p>

<ol>
  <li><strong>Campo de deformación inter-hoja adimensional:</strong> &chi;<sub>n</sub>(x<sup>&mu;</sup>) &equiv; (&xi;<sub>n+1</sub> - &xi;<sub>n</sub>) / d<sub>0</sub> para n &isin; {0, &hellip;, 15}, definiendo los modos transversales relativos de compresión del bulk pentadimensional.</li>
  <li><strong>Potencial biestable con colapso de barrera dependiente de la exposición:</strong> V(&chi;<sub>n</sub>, &Epsilon;<sub>n</sub>) = &frac12; m<sub>0</sub><sup>2</sup> (1 - &Epsilon;<sub>n</sub>/&Epsilon;<sub>c</sub>) &chi;<sub>n</sub><sup>2</sup> - &frac13; &mu; &chi;<sub>n</sub><sup>3</sup> + &frac14; &lambda; &chi;<sub>n</sub><sup>4</sup>, donde el funcional de exposición derivada acumulada &Epsilon;<sub>n</sub> suprime la barrera ante bombeo resonante armónico continuo.</li>
  <li><strong>Parámetro de relatividad intrínseca (&real;):</strong> &real; = T<sup>&mu;</sup><sub>&mu;</sub> / &radic;(T<sub>&alpha;&beta;</sub> T<sup>&alpha;&beta;</sup>), discriminando analíticamente entre:
    <ul>
      <li><em>Histéresis cerrada (&real; &rarr; 1):</em> Dominio de masa bariónica/corteza densa, donde la divergencia total confina la interacción a términos de contorno superficial con apantallamiento &le; -112 dB.</li>
      <li><em>Histéresis abierta (&real; &rarr; 0):</em> Dominio radiativo y exposición a la función derivada (&part;<sub>&mu;</sub> &chi;<sub>n</sub>)<sup>2</sup>, activando la penetración volumétrica a través de las 17 interfaces de salto.</li>
    </ul>
  </li>
  <li><strong>Partículas como transiciones de fase y semivida del solitón:</strong> Formalización de las partículas como condensados locales de foliación con semivida gobernada por la acción euclídea del bounce S<sub>E</sub>(A) &prop; &sigma;<sub>fol</sub><sup>4</sup> / [&Delta;V<sub>0</sub> - &kappa; A<sup>2</sup>]<sup>3</sup>, fundamentando la estabilidad asintótica de baja capa frente a la inestabilidad de alta superposición.</li>
  <li><strong>Validación experimental sub-0.01%:</strong> Calibración de los cuatro pilares (birrefringencia quiral intracavidad +59.97 dB, tensor de Poisson anómalo &nu; = 0.4201, triangulación en red de 500 m con retardo m&aacute;ximo de &plusmn;87.31 s, y filtrado adaptado multicanal con &rho; = 312.36&sigma;).</li>
</ol>
"""

# Comprobar si ya contiene la actualización para evitar duplicados
desc_actual = metadata.get("description", "")
if "Reformulación No Lineal del Lagrangiano Efectivo" not in desc_actual:
    metadata["description"] = desc_actual + descripcion_adicional
else:
    print("[*] La descripción ya contiene la sección técnica. Procediendo a refrescar...")

# Paso 3: Enviar metadatos actualizados (PUT)
r_put = requests.put(BASE_URL, params=params, data=json.dumps({"metadata": metadata}), headers=headers)
if r_put.status_code == 200:
    print("[+] Metadatos actualizados correctamente en el borrador.")
else:
    print(f"[-] Error al guardar metadatos: {r_put.status_code} {r_put.text}")
    exit(1)

# Paso 4: Publicar los cambios para que sean visibles en records/23247617
print("[*] Publicando cambios en Zenodo...")
r_pub = requests.post(f"{BASE_URL}/actions/publish", params=params)
if r_pub.status_code in [202, 200]:
    print(f"[+] ¡Éxito! Registro publicado y visible en: https://zenodo.org/records/{DEPOSITION_ID}")
else:
    print(f"[-] Error en publicación: {r_pub.status_code} {r_pub.text}")
