#!/usr/bin/env bash
set -euo pipefail

# Download and extract a prebuilt CReM fragment database from Zenodo.
# Usage:
#   bash scripts/download_crem_db.sh
#   bash scripts/download_crem_db.sh chembl33_sa2_f5
#   bash scripts/download_crem_db.sh chembl33_sa25_f5 /custom/output/dir

DB_NAME="${1:-chembl33_sa25_f5}"
OUT_DIR="${2:-}"

case "${DB_NAME}" in
  chembl22_sa2)
    FILE_NAME="chembl22_sa2.db.gz"
    ;;
  chembl22_sa2_hac12)
    FILE_NAME="chembl22_sa2_hac12.db.gz"
    ;;
  chembl22_sa25_hac12)
    FILE_NAME="chembl22_sa25_hac12.db.gz"
    ;;
  chembl33_sa2_f5)
    FILE_NAME="chembl33_sa2_f5.db.gz"
    ;;
  chembl33_sa25_f5)
    FILE_NAME="chembl33_sa25_f5.db.gz"
    ;;
  chembl33_f5)
    FILE_NAME="chembl33_f5.db.gz"
    ;;
  enamine2025_sa2_f5)
    FILE_NAME="enamine2025_sa2_f5.db.gz"
    ;;
  *)
    echo "Unsupported DB name: ${DB_NAME}" >&2
    echo "Supported values:" >&2
    echo "  chembl22_sa2" >&2
    echo "  chembl22_sa2_hac12" >&2
    echo "  chembl22_sa25_hac12" >&2
    echo "  chembl33_sa2_f5" >&2
    echo "  chembl33_sa25_f5" >&2
    echo "  chembl33_f5" >&2
    echo "  enamine2025_sa2_f5" >&2
    exit 1
    ;;
esac

if [[ -z "${OUT_DIR}" ]]; then
  SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
  REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
  OUT_DIR="${REPO_ROOT}/data/crem_db"
fi

URL="https://zenodo.org/records/16909329/files/${FILE_NAME}?download=1"
GZ_PATH="${OUT_DIR}/${FILE_NAME}"
DB_PATH="${GZ_PATH%.gz}"

mkdir -p "${OUT_DIR}"

echo "Target DB     : ${DB_NAME}"
echo "Download URL  : ${URL}"
echo "Output folder : ${OUT_DIR}"
echo

if command -v wget >/dev/null 2>&1; then
  wget -c -O "${GZ_PATH}" "${URL}"
elif command -v curl >/dev/null 2>&1; then
  curl -L -C - -o "${GZ_PATH}" "${URL}"
else
  echo "Error: neither wget nor curl is installed." >&2
  exit 1
fi

if [[ ! -s "${GZ_PATH}" ]]; then
  echo "Error: download failed or empty file: ${GZ_PATH}" >&2
  exit 1
fi

echo "Extracting ${GZ_PATH} ..."
gzip -dk "${GZ_PATH}"

if [[ ! -s "${DB_PATH}" ]]; then
  echo "Error: extraction failed: ${DB_PATH} not found." >&2
  exit 1
fi

echo
echo "Done."
echo "DB file: ${DB_PATH}"
echo
ls -lh "${GZ_PATH}" "${DB_PATH}"
