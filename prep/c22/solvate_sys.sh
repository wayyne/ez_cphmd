#!/usr/bin/env bash
set -euo pipefail

# =====================================================================
# PATH DEFINITIONS & INITIALIZATION
# =====================================================================
# Directory hosting this wrapper and its CHARMM input script.
SCRIPT_DIR="$(
    cd -- "$(dirname -- "${BASH_SOURCE[0]}")"
    pwd -P
)"

# Initialize default parameters
INNAME=""
CUSHION=""
BOX_TYPE=""
BOXID=""
XTL=""
ALPHA=""
BETA=""
GAMMA=""

# User-facing salt concentration is now mM.
# Example: --conc 50 means 50 mM NaCl.
CONC_MM="150"     # Default physiological ionic concentration (mM)

# Internal CHARMM-facing concentration in mol/L.
# This is computed after validation and passed to CHARMM as CONC.
CONC_M=""

TITR_FLAG=""

# Optional exact-ion override.
# AUTO  = compute SOD/CLA from concentration and neutralization.
# EXACT = use exactly --sod/--cla values, no concentration/neutralization math.
ION_MODE="AUTO"
SOD_EXACT=""
CLA_EXACT=""
EXACT_IONS_REQUESTED="0"

# By default, remove intermediate files after a successful CHARMM run.
# Use --keep-temp to preserve them for debugging.
KEEP_TEMP="0"

usage() {
  echo "Usage: $0 --in <name> --box <cubic|octahedral|rhombic> --cushion <val> [--conc <mM>] [--titr <-1|0|1>] [--sod <N>] [--cla <N>] [--keep-temp]"
  echo
  echo "Salt concentration:"
  echo "  --conc is in mM."
  echo "  Example: --conc 50 means 50 mM NaCl."
  echo
  echo "Ion modes:"
  echo "  AUTO  mode: omit --sod/--cla; SOD/CLA are computed from --conc and neutralization."
  echo "  EXACT mode: provide --sod and/or --cla; exact values are used and --conc is ignored for bulk ions."
}

# =====================================================================
# 1. ARGUMENT PARSING (COMMAND LINE FLAGS)
# =====================================================================
while [[ $# -gt 0 ]]; do
    case "$1" in
        --in)
            INNAME="${2:-}"; shift 2 ;;
        --cushion)
            CUSHION="${2:-}"; shift 2 ;;
        --box)
            BOX_TYPE="$(echo "${2:-}" | tr '[:upper:]' '[:lower:]')"; shift 2 ;;
        --conc|--conc-mm|--salt-mm)
            CONC_MM="${2:-}"; shift 2 ;;
        --sod)
            SOD_EXACT="${2:-}"; EXACT_IONS_REQUESTED="1"; shift 2 ;;
        --cla)
            CLA_EXACT="${2:-}"; EXACT_IONS_REQUESTED="1"; shift 2 ;;
        --titr)
            TITR_FLAG="${2:-}"; shift 2 ;;
        --keep-temp)
            KEEP_TEMP="1"; shift ;;
        -h|--help)
            usage
            exit 0 ;;
        *)
            echo "Unknown option: $1"
            usage
            exit 1 ;;
    esac
done

# =====================================================================
# 2. INTERACTIVE PROMPTS (TRIGGERED ONLY IF FLAGS ARE MISSING)
# =====================================================================
if [[ -z "$INNAME" ]]; then
    read -r -p "Enter the base name of your input files (without .psf/.crd): " INNAME
fi

if [[ -z "$BOX_TYPE" ]]; then
    echo "-----------------------------------------"
    echo "Select Unit Cell Geometry:"
    echo " 1) Cubic"
    echo " 2) Truncated Octahedral"
    echo " 3) Rhombic Dodecahedron"
    echo "-----------------------------------------"
    read -r -p "Enter choice [1, 2, or 3]: " BOX_TYPE
    BOX_TYPE="$(echo "$BOX_TYPE" | tr '[:upper:]' '[:lower:]')"
fi

if [[ -z "$CUSHION" ]]; then
    read -r -p "Enter the water cushion size in Angstroms (e.g., 15): " CUSHION
fi

if [[ -z "$TITR_FLAG" ]]; then
    echo "-----------------------------------------"
    echo "Enable Titratable Water Co-ions (pHMD)?"
    echo " 0) No (Standard MD / No Co-ions)"
    echo " -1) Yes (Acidic Titration Co-ions: TIPU)"
    echo " 1) Yes (Basic Titration Co-ions: TIPP)"
    echo "-----------------------------------------"
    read -r -p "Enter choice [0, -1, or 1]: " TITR_FLAG
fi

# =====================================================================
# 3. VALIDATION
# =====================================================================
if [[ -z "$INNAME" ]]; then
    echo "Error: --in cannot be empty"
    exit 1
fi

if [[ ! -f "${INNAME}.psf" ]]; then
    echo "Error: Missing input PSF: ${INNAME}.psf"
    exit 1
fi

if [[ ! -f "${INNAME}.crd" ]]; then
    echo "Error: Missing input CRD: ${INNAME}.crd"
    exit 1
fi

if [[ ! -f "${charmm_scripts}/solvate_sys.inp" ]]; then
    echo "Error: Missing CHARMM input script: ${charmm_scripts}/solvate_sys.inp"
    exit 1
fi

if ! [[ "$CUSHION" =~ ^[0-9]+([.][0-9]+)?$ ]]; then
    echo "Error: --cushion must be a non-negative number"
    exit 1
fi

if ! [[ "$CONC_MM" =~ ^[0-9]+([.][0-9]+)?$ ]]; then
    echo "Error: --conc must be a non-negative number in mM"
    exit 1
fi

# Soft warning for the common old-style mistake:
# Old interface: --conc 0.05 meant 50 mM.
# New interface: --conc 50 means 50 mM.
if awk -v c="$CONC_MM" 'BEGIN { exit !(c > 0.0 && c < 1.0) }'; then
    echo "WARNING: --conc now expects mM, not mol/L."
    echo "         You entered --conc ${CONC_MM}, which means ${CONC_MM} mM."
    echo "         For 50 mM, use --conc 50."
fi

# Convert user-facing mM to CHARMM-facing mol/L.
CONC_M="$(awk -v mm="$CONC_MM" 'BEGIN { printf "%.10g", mm / 1000.0 }')"

if [[ -n "$SOD_EXACT" ]] && ! [[ "$SOD_EXACT" =~ ^[0-9]+$ ]]; then
    echo "Error: --sod must be a non-negative integer"
    exit 1
fi

if [[ -n "$CLA_EXACT" ]] && ! [[ "$CLA_EXACT" =~ ^[0-9]+$ ]]; then
    echo "Error: --cla must be a non-negative integer"
    exit 1
fi

if [[ "$EXACT_IONS_REQUESTED" == "1" ]]; then
    ION_MODE="EXACT"
    SOD_EXACT="${SOD_EXACT:-0}"
    CLA_EXACT="${CLA_EXACT:-0}"
else
    ION_MODE="AUTO"
    SOD_EXACT="0"
    CLA_EXACT="0"
fi

case "$TITR_FLAG" in
    -1|0|1) ;;
    *)
        echo "Error: --titr must be one of: -1, 0, 1"
        exit 1 ;;
esac

# =====================================================================
# 4. CANONICAL UNIT CELL MAPPING
# =====================================================================
# BOXID is the canonical flag intended for the updated solvate_sys.inp:
#   1 = Cubic
#   2 = Truncated Octahedral
#   3 = Rhombic Dodecahedron
#
# XTLtype and angles are still passed for compatibility/debug visibility.

case "$BOX_TYPE" in
    cubic|cube|cub|1)
        BOX_TYPE="cubic"
        BOXID="1"
        XTL="CUBIC"
        ALPHA="90.0"
        BETA="90.0"
        GAMMA="90.0"
        ;;

    octahedral|octa|oct|truncated_octahedral|truncated-octahedral|2)
        BOX_TYPE="octahedral"
        BOXID="2"
        XTL="OCTAHEDRAL"
        ALPHA="109.471220634"
        BETA="109.471220634"
        GAMMA="109.471220634"
        ;;

    rhombic|rhdo|dodecahedron|rhombic_dodecahedron|rhombic-dodecahedron|3)
        BOX_TYPE="rhombic"
        BOXID="3"
        XTL="RHDO"
        ALPHA="60.0"
        BETA="90.0"
        GAMMA="60.0"
        ;;

    *)
        echo "Error: Unknown box type '$BOX_TYPE'"
        usage
        exit 1
        ;;
esac

# =====================================================================
# 4b. OUTPUT NAMING
# =====================================================================
# Naming pattern:
#   AUTO mode : <in>_solv_<cushion>A_<box>_<conc>mM
#   EXACT mode: <in>_solv_<cushion>A_<box>_<NxSOD>_<NxCLA>
#
# Examples:
#   capped_gdg_solv_15A_cubic_50mM
#   capped_gdg_solv_15A_octahedral_1xSOD
#   capped_gdg_solv_10A_rhombic_2xSOD_1xCLA

sanitize_num_tag() {
    local value="$1"

    # Strip trailing .0, .00, etc.
    value="$(echo "$value" | sed -E 's/\.0+$//')"

    # Replace decimal point with p for filesystem-safe names.
    value="${value//./p}"

    echo "$value"
}

format_mM_tag() {
    local value="$1"

    awk -v mm="$value" '
        BEGIN {
            if (mm == int(mm)) {
                printf "%dmM", mm
            } else {
                s = sprintf("%.3f", mm)
                sub(/0+$/, "", s)
                sub(/\.$/, "", s)
                gsub(/\./, "p", s)
                printf "%smM", s
            }
        }'
}

CUSHION_TAG="$(sanitize_num_tag "$CUSHION")"

if [[ "$ION_MODE" == "EXACT" ]]; then
    ION_TAG=""

    if [[ "$SOD_EXACT" -gt 0 ]]; then
        ION_TAG="${SOD_EXACT}xsod"
    fi

    if [[ "$CLA_EXACT" -gt 0 ]]; then
        if [[ -n "$ION_TAG" ]]; then
            ION_TAG="${ION_TAG}_${CLA_EXACT}xcla"
        else
            ION_TAG="${CLA_EXACT}xcla"
        fi
    fi

    if [[ -z "$ION_TAG" ]]; then
        ION_TAG="0ions"
    fi
else
    ION_TAG="$(format_mM_tag "$CONC_MM")"
fi

OUTNAME="${INNAME}_solv_${CUSHION_TAG}A_${BOX_TYPE}_${ION_TAG}"
LOGFILE="solvate_${OUTNAME}.out"

# =====================================================================
# 5. RUNTIME LOG DISPLAY & CHARMM EXECUTION
# =====================================================================
echo "========================================="
echo "Solvation System Variables Confirmed:"
echo "Input System     : $INNAME"
echo "Output System    : $OUTNAME"
echo "Box Type         : $BOX_TYPE"
echo "BOXID            : $BOXID"
echo "CHARMM XTL       : $XTL"
echo "Angles           : $ALPHA $BETA $GAMMA"
echo "Cushion Size     : $CUSHION Å"
echo "Salt Conc        : $CONC_MM mM"
echo "Salt Conc CHARMM : $CONC_M M"
echo "Titration Flag   : $TITR_FLAG"
echo "Ion Mode         : $ION_MODE"
echo "Exact SOD        : $SOD_EXACT"
echo "Exact CLA        : $CLA_EXACT"
echo "Keep Temp Files  : $KEEP_TEMP"
echo "Log File         : $LOGFILE"
echo "========================================="

# Do the damn thing.
# The wrapper accepts --conc in mM, but solvate_sys.inp expects CONC in mol/L.
if [[ -z "${C22_TOPPAR:-}" ]]; then
    echo "Error: C22_TOPPAR is not set. Check prep/config.sh."
    exit 1
fi

charmm TOPPAR="$C22_TOPPAR" INNAME="$INNAME" OUTNAME="$OUTNAME" CUSHION="$CUSHION" \
       BOXID="$BOXID" XTLtype="$XTL" A_ANG="$ALPHA" B_ANG="$BETA" G_ANG="$GAMMA" \
       CONC="$CONC_M" CONCMM="$CONC_MM" TITRFLAG="$TITR_FLAG" \
       IONMODE="$ION_MODE" NSOD="$SOD_EXACT" NCLA="$CLA_EXACT" \
       < "${SCRIPT_DIR}/solvate_sys.inp" > "$LOGFILE"

if [[ "$KEEP_TEMP" == "0" ]]; then
    rm -f prot_temp.crd prot_temp.psf water_tmp.crd crystal_image.str
else
    echo "Temporary files preserved:"
    echo "  prot_temp.crd"
    echo "  prot_temp.psf"
    echo "  water_tmp.crd"
    echo "  crystal_image.str"
fi

echo "Process complete. Output saved to ${OUTNAME}.psf/.crd/.pdb"
echo "CHARMM log saved to $LOGFILE"
