#!/usr/bin/env bash
set -euo pipefail

ZENODO_TOKEN="1y5FXoiybVkOdndx07C80fqWShrK8U8Tttn2mZzu9UBdAfPSsQGcqUY4idT1"
BASE_RECORD="23135100"
FILE_TAR="17layer_double_octet_analysis_v2.1.0.tar.gz"
FILE_SHA="17layer_double_octet_analysis_v2.1.0.tar.gz.sha256"

echo "================================================================================"
echo " SUBIENDO RELEASE v2.1.0 A ZENODO (Record Base: ${BASE_RECORD})"
echo "================================================================================"

# 1. Crear nueva versión en Zenodo
echo "[1/4] Creando borrador de nueva versión..."
DRAFT_RES=$(curl -s -f -X POST \
  "https://zenodo.org/api/deposit/depositions/${BASE_RECORD}/actions/newversion?access_token=${ZENODO_TOKEN}")

NEW_ID=$(echo "${DRAFT_RES}" | python3 -c "import sys, json; print(json.load(sys.stdin)['links']['latest_draft'].split('/')[-1])")
echo "      -> ID del nuevo borrador: ${NEW_ID}"

# 2. Obtener bucket de subida
BUCKET_RES=$(curl -s -f "https://zenodo.org/api/deposit/depositions/${NEW_ID}?access_token=${ZENODO_TOKEN}")
BUCKET_URL=$(echo "${BUCKET_RES}" | python3 -c "import sys, json; print(json.load(sys.stdin)['links']['bucket'])")
echo "      -> Bucket URL asignado: ${BUCKET_URL}"

# 3. Subir el archivo .tar.gz y el .sha256
echo "[2/4] Subiendo ${FILE_TAR} al bucket..."
curl -f -X PUT "${BUCKET_URL}/${FILE_TAR}?access_token=${ZENODO_TOKEN}" \
  -H "Content-Type: application/octet-stream" \
  --data-binary "@${FILE_TAR}" \
  --progress-bar

echo "[3/4] Subiendo suma criptográfica SHA-256..."
curl -f -X PUT "${BUCKET_URL}/${FILE_SHA}?access_token=${ZENODO_TOKEN}" \
  -H "Content-Type: text/plain" \
  --data-binary "@${FILE_SHA}"

# 4. Actualizar metadatos del release
echo "[4/4] Actualizando metadatos descriptivos..."
curl -s -f -X PUT "https://zenodo.org/api/deposit/depositions/${NEW_ID}?access_token=${ZENODO_TOKEN}" \
  -H "Content-Type: application/json" \
  -d '{
    "metadata": {
      "title": "Analytic & Numerical Validation Suite for 17-Layered Foliated Bulk Gravity: Double-Octet Balance and Multidomain Phenomenology",
      "version": "2.1.0",
      "upload_type": "software",
      "publication_type": "software",
      "description": "<p>Numerical validation suite providing full reproducibility for the 17-layer foliated 5D bulk gravity model. Validates the 8(+) vs 8(-) quantum double-octet phase balance (602.50 GeV radion), 3.33 TeV parity-doubled graviton resonance in LHC, PeV neutrino inelastic absorption in Earth core (IceCube), and the 0.19 mHz stochastic gravitational wave background (LISA). Global statistical significance: <strong>Z = 4.87 sigma</strong> (p = 5.50e-7).</p>",
      "creators": [
        {
          "name": "Cruañas Pardo, Albert",
          "orcid": "0009-0007-0846-2321",
          "affiliation": "Independent Theoretical Physics Researcher"
        }
      ],
      "access_right": "open",
      "license": "mit",
      "keywords": [
        "Kaluza-Klein Graviton",
        "17-Layered Foliated Bulk",
        "Double Octet Balance",
        "Coleman-Weinberg",
        "LHC Phenomenology",
        "IceCube Neutrino Absorption",
        "LISA Gravitational Waves"
      ]
    }
  }' > /dev/null

echo "================================================================================"
echo " ARCHIVOS Y METADATOS SUBIDOS EXITOSAMENTE EN MODO BORRADOR"
echo " Puedes revisarlo en tu panel de Zenodo:"
echo " https://zenodo.org/deposit/${NEW_ID}"
echo "================================================================================"
