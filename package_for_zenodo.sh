#!/usr/bin/env bash
set -e

# Detectar si estamos dentro o fuera de radion_tensor_analysis
if [ -d "modules" ] && [ -f "run_full_analysis.py" ]; then
    TARGET_DIR="."
    ARCHIVE_DIR="radion_tensor_analysis"
    cd ..
elif [ -d "radion_tensor_analysis" ]; then
    TARGET_DIR="radion_tensor_analysis"
    ARCHIVE_DIR="radion_tensor_analysis"
else
    echo "[-] Error: No se localizó la carpeta de la suite."
    exit 1
fi

VERSION="1.0.0"
DATE="2026-10-08"

echo "[*] Generando metadatos para Zenodo en: ${TARGET_DIR}"

# 1. codemeta.json (Estándar interoperable de metadatos de software científico)
cat << 'META_EOF' > "${TARGET_DIR}/codemeta.json"
{
  "@context": "https://doi.org/10.5063/schema/codemeta-2.0",
  "@type": "SoftwareSourceCode",
  "name": "radion_tensor_analysis",
  "version": "1.0.0",
  "description": "Reproducible profile likelihood and multitensor decomposition suite (Ranks 0 to 4) for the 602.50 GeV modular scalar radion in 17-layer warped bulk foliations.",
  "datePublished": "2026-10-08",
  "author": [
    {
      "@type": "Person",
      "givenName": "Albert",
      "familyName": "Cruañas Pardo",
      "@id": "https://orcid.org/0009-0007-0846-2321",
      "affiliation": {
        "@type": "Organization",
        "name": "Independent Researcher in Theoretical Physics"
      }
    }
  ],
  "license": "https://spdx.org/licenses/MIT",
  "programmingLanguage": "Python",
  "runtimePlatform": "Python 3.8+",
  "keywords": [
    "High Energy Physics",
    "Radion",
    "Warped Extra Dimensions",
    "Profile Likelihood",
    "Multitensor Decomposition",
    "LHC Open Data",
    "SMEFT"
  ]
}
META_EOF

# 2. .zenodo.json (Autocompletado directo de campos en la interfaz de Zenodo)
cat << 'ZEN_EOF' > "${TARGET_DIR}/.zenodo.json"
{
  "title": "radion_tensor_analysis: Multitensor Statistical Validation Suite for the 602.50 GeV Modular Scalar Radion",
  "upload_type": "software",
  "description": "Open-source computational implementation of the multitensor profile likelihood framework for validating the 602.50 GeV scalar radion theoretical model across Ranks 0 to 4 using LHC Run 2/3 and IceCube data.",
  "creators": [
    {
      "name": "Cruañas Pardo, Albert",
      "affiliation": "Independent Researcher in Theoretical Physics",
      "orcid": "0009-0007-0846-2321"
    }
  ],
  "access_right": "open",
  "license": "MIT",
  "version": "1.0.0",
  "keywords": [
    "Theoretical Physics",
    "High Energy Physics",
    "Radion",
    "Warped Extra Dimensions",
    "Standard Model Effective Field Theory",
    "LHC",
    "IceCube"
  ]
}
ZEN_EOF

# 3. Empaquetar tarball excluyendo caches
TAR_NAME="radion_tensor_analysis_v${VERSION}.tar.gz"
echo "[*] Comprimiendo el archivo de distribución: ${TAR_NAME}"

tar --exclude='__pycache__' \
    --exclude='*.pyc' \
    --exclude='.git' \
    -czvf "${TAR_NAME}" "${ARCHIVE_DIR}"

echo ""
echo "======================================================================"
echo "[+] Paquete para Zenodo generado con éxito:"
echo "    -> Archivo: $(pwd)/${TAR_NAME}"
echo "    -> Tamaño:  $(du -h "${TAR_NAME}" | cut -f1)"
echo "======================================================================"
echo "Instrucciones de subida a Zenodo:"
echo "1. Ve a https://zenodo.org/deposit/new"
echo "2. Arrastra '${TAR_NAME}' a la sección 'Files'."
echo "3. Zenodo leerá '.zenodo.json' y rellenará el formulario automáticamente."
echo "4. Haz clic en 'Publish' para obtener el DOI oficial."
echo "======================================================================"
