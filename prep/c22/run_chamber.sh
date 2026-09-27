#!/usr/bin/env bash
set -euo pipefail

usage() {
  echo "Usage: $0 <toppar_dir> <system_prefix>"
  echo
  echo "Example:"
  echo "  $0 ../../toppar capped_gkg_solv_10A_octahedral_50mM"
}

if [[ $# -ne 2 ]]; then
  usage
  exit 1
fi

# Path to c22 param parent directory
TOPPAR="$1"

# Prefix for the system to convert, no file extension
SYS="$2"

if [[ ! -d "$TOPPAR" ]]; then
  echo "Error: TOPPAR directory not found: $TOPPAR"
  exit 1
fi

if [[ ! -f "${SYS}.psf" ]]; then
  echo "Error: Missing PSF file: ${SYS}.psf"
  exit 1
fi

if [[ ! -f "${SYS}.crd" ]]; then
  echo "Error: Missing CRD file: ${SYS}.crd"
  exit 1
fi

# Find matching solvation log.
# Uses case-insensitive matching so these both work:
#   capped_gkg_solv_10A_octahedral_50mM
#   capped_gkg_solv_10a_octahedral_50mm
mapfile -d '' LOG_MATCHES < <(
  find . -maxdepth 1 -type f -iname "solvate_${SYS}.out" -print0
)

if [[ "${#LOG_MATCHES[@]}" -eq 0 ]]; then
  echo "Error: Could not find matching solvation log: solvate_${SYS}.out"
  echo "Available solvation logs:"
  ls -1 solvate_*.out 2>/dev/null || true
  exit 1
fi

if [[ "${#LOG_MATCHES[@]}" -gt 1 ]]; then
  echo "Error: Multiple matching solvation logs found for SYS=${SYS}:"
  printf '  %s\n' "${LOG_MATCHES[@]}"
  exit 1
fi

LOGFILE="${LOG_MATCHES[0]}"

read -r A B C < <(
  awk '
    toupper($1) == "CELL" && toupper($2) == "DIMENSIONS:" {
      a=$3; b=$4; c=$5
    }
    END {
      if (a != "") print a, b, c
    }
  ' "$LOGFILE"
)

read -r ALPHA BETA GAMMA < <(
  awk '
    toupper($1) == "CELL" && toupper($2) == "ANGLES:" {
      a=$3; b=$4; c=$5
    }
    END {
      if (a != "") print a, b, c
    }
  ' "$LOGFILE"
)

is_number() {
  [[ "$1" =~ ^[0-9]+([.][0-9]+)?([eE][-+]?[0-9]+)?$ ]]
}

if [[ -z "${A:-}" || -z "${B:-}" || -z "${C:-}" || -z "${ALPHA:-}" || -z "${BETA:-}" || -z "${GAMMA:-}" ]]; then
  echo "Error: Failed to parse box dimensions/angles from: $LOGFILE"
  exit 1
fi

for value in "$A" "$B" "$C" "$ALPHA" "$BETA" "$GAMMA"; do
  if ! is_number "$value"; then
    echo "Error: Parsed non-numeric box value from $LOGFILE: $value"
    exit 1
  fi
done

BOX="${A},${B},${C},${ALPHA},${BETA},${GAMMA}"

echo "Using solvation log : $LOGFILE"
echo "Using Amber box     : $BOX"

cat > "chamber_${SYS}.in" <<EOF
chamber \\
  -top ${TOPPAR}/prot/top_all22_prot.rtf \\
  -param ${TOPPAR}/prot/par_all22_prot.prm \\
  -str ${TOPPAR}/toppar_water_ions.str \\
  -str ${TOPPAR}/toppar_phmd_c22_foramber.str \\
  -psf ${SYS}.psf \\
  -crd ${SYS}.crd \\
  -box ${BOX}

outparm ${SYS}.prmtop ${SYS}.inpcrd
quit
EOF

parmed -i "chamber_${SYS}.in" > "chamber_${SYS}.log"
