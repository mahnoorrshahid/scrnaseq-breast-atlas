#!/usr/bin/env bash
# scripts/00_download_data.sh
#
# Downloads the GSE176078 supplementary files (Wu et al. 2021 breast cancer
# scRNA-seq atlas) into data/raw/.
#
# GEO supplementary file listings can change their exact filenames over time,
# so this script lists what it finds and lets you confirm before extracting —
# check that against data/README.md's "expected contents" section.

set -euo pipefail

DEST="$(dirname "$0")/../data/raw"
mkdir -p "$DEST"

GEO_ACC="GSE176078"
GEO_FTP="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE176nnn/${GEO_ACC}/suppl/"

echo "Fetching supplementary file listing for ${GEO_ACC}..."
echo "  ${GEO_FTP}"
echo

# List what's available (requires network access to NCBI FTP)
curl -sSL "${GEO_FTP}" | grep -oE 'href="[^"]+"' | sed 's/href="//;s/"$//' || {
  echo "Could not reach NCBI FTP directly from this machine."
  echo "Open ${GEO_FTP} in a browser instead, or use the GEO accession"
  echo "viewer: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=${GEO_ACC}"
  exit 1
}

echo
echo "Review the filenames above, then download the relevant archive, e.g.:"
echo "  curl -sSL -o \"${DEST}/GSE176078_supp.tar.gz\" \"${GEO_FTP}<exact_filename>\""
echo "  tar -xzf \"${DEST}/GSE176078_supp.tar.gz\" -C \"${DEST}\""
echo
echo "Then confirm the extracted files match data/README.md before running 01_qc_filtering.R."
