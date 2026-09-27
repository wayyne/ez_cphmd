#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# Amber24 CpHMD dry peptide generator.
#
# ff19sb : ff19SB + CpHMD params
# ff14sb : ff14SB + ff14sb_pme_dayhoff_dev.parm + CpHMD params
#
# Capped   : ACE / NHE
# Uncapped : zwitterionic Amber termini
#
# CpHMD residue mapping:
#   D -> AS2
#   E -> GL2
#   H -> HIP
#   C -> CYS
#   K -> LYS
###############################################################################

AMBER_SH="${AMBER_ENV:-}"

RAW_SEQ=""
FF_CHOICE="ff19sb"
CAP_CHOICE=""
CAPPED_FLAG="0"
UNCAPPED_FLAG="0"

usage()
{
    echo "Usage: $0 [--pep <SEQ>] [--ff <ff19sb|ff14sb>] [--capped | --uncapped]"
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --pep)
            if [[ -n "${2:-}" && "${2:-}" != --* ]]; then
                RAW_SEQ="$2"
                shift 2
            else
                echo "Error: --pep requires a sequence argument."
                exit 1
            fi
            ;;
        --ff|--forcefield)
            if [[ -n "${2:-}" && "${2:-}" != --* ]]; then
                FF_CHOICE="$(printf "%s" "$2" | tr '[:upper:]' '[:lower:]')"
                shift 2
            else
                echo "Error: --ff requires ff19sb or ff14sb."
                exit 1
            fi
            ;;
        --capped)
            CAPPED_FLAG="1"
            CAP_CHOICE="1"
            shift
            ;;
        --uncapped)
            UNCAPPED_FLAG="1"
            CAP_CHOICE="2"
            shift
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "Error: unknown argument '$1'"
            usage
            exit 1
            ;;
    esac
done

case "$FF_CHOICE" in
    ff19sb|ff19)
        FF_CHOICE="ff19sb"
        FF_DIR="ff19"
        PROTEIN_LEAPRC="leaprc.protein.ff19SB"
        EXTRA_FF_PARAMS=""
        ;;
    ff14sb|ff14)
        FF_CHOICE="ff14sb"
        FF_DIR="ff14"
        PROTEIN_LEAPRC="leaprc.protein.ff14SB"
        EXTRA_FF_PARAMS="ff14sb_pme_dayhoff_dev.parm"
        ;;
    *)
        echo "Error: --ff must be ff19sb or ff14sb."
        exit 1
        ;;
esac

if [[ "$FF_CHOICE" == "ff14sb" ]]; then
    PHMD_DIR="${FF14_PHMD_DIR:-}"
else
    PHMD_DIR="${FF19_PHMD_DIR:-}"
fi
PHMD_LIB="${PHMD_DIR}/phmd.lib"
PHMD_FRCMOD="${PHMD_DIR}/frcmod.phmd_barrier"

if [[ -n "$EXTRA_FF_PARAMS" ]]; then
    EXTRA_FF_PARAMS_PATH="${PHMD_DIR}/${EXTRA_FF_PARAMS}"
else
    EXTRA_FF_PARAMS_PATH=""
fi

if [[ "$CAPPED_FLAG" -eq 1 && "$UNCAPPED_FLAG" -eq 1 ]]; then
    echo "Error: cannot pass both --capped and --uncapped."
    exit 1
fi

if [[ -z "$RAW_SEQ" ]]; then
    read -r -p "Enter the single-letter peptide sequence: " RAW_SEQ
fi

if [[ -z "$RAW_SEQ" ]]; then
    echo "Error: sequence cannot be empty."
    exit 1
fi

SEQ_IN="$(printf "%s" "$RAW_SEQ" | tr '[:lower:]' '[:upper:]')"
PEPLEN="${#SEQ_IN}"

if [[ -z "$CAP_CHOICE" ]]; then
    echo "-----------------------------------------"
    echo "Select Capping Style:"
    echo " 1) Capped   (Neutral: ACE / NHE)"
    echo " 2) Uncapped (Zwitterionic termini)"
    echo "-----------------------------------------"
    read -r -p "Enter choice [1 or 2]: " CAP_CHOICE
fi

if [[ "$CAP_CHOICE" == "1" ]]; then
    CAP_LABEL="capped"
    NCAP="ACE"
    CCAP="NHE"
elif [[ "$CAP_CHOICE" == "2" ]]; then
    CAP_LABEL="zwitter"
    NCAP=""
    CCAP=""
else
    echo "Error: invalid capping choice."
    exit 1
fi

map_seq()
{
    python3 - "$SEQ_IN" <<'PY'
import sys

seq = sys.argv[1]

aa = {
    "A": "ALA",
    "R": "ARG",
    "N": "ASN",
    "D": "AS2",
    "C": "CYS",
    "Q": "GLN",
    "E": "GL2",
    "G": "GLY",
    "H": "HIP",
    "I": "ILE",
    "L": "LEU",
    "K": "LYS",
    "M": "MET",
    "F": "PHE",
    "P": "PRO",
    "S": "SER",
    "T": "THR",
    "W": "TRP",
    "Y": "TYR",
    "V": "VAL",
}

out = []

for ch in seq:
    if ch not in aa:
        sys.exit("Error: unknown amino acid character '%s'" % ch)
    out.append(aa[ch])

print(" ".join(out))
PY
}

PEPRSQ="$(map_seq)"
OUTNAME="${CAP_LABEL}_${SEQ_IN}"
LEAPIN="gen_${OUTNAME}.leap"
LOGFILE="gen_${OUTNAME}.out"

if [[ "$CAP_CHOICE" == "1" ]]; then
    LEAP_SEQ="${NCAP} ${PEPRSQ} ${CCAP}"
else
    LEAP_SEQ="${PEPRSQ}"
fi

if [[ -z "$PHMD_DIR" ]]; then
    echo "Error: CpHMD parameter directory is not set for $FF_CHOICE. Check prep/config.sh."
    exit 1
fi

if [[ ! -f "$AMBER_SH" ]]; then
    echo "Error: missing Amber environment: $AMBER_SH"
    exit 1
fi

if [[ ! -f "$PHMD_LIB" ]]; then
    echo "Error: missing CpHMD library: $PHMD_LIB"
    exit 1
fi

if [[ ! -f "$PHMD_FRCMOD" ]]; then
    echo "Error: missing CpHMD frcmod: $PHMD_FRCMOD"
    exit 1
fi

if [[ -n "$EXTRA_FF_PARAMS_PATH" && ! -f "$EXTRA_FF_PARAMS_PATH" ]]; then
    echo "Error: missing force-field parameter file: $EXTRA_FF_PARAMS_PATH"
    exit 1
fi

# shellcheck disable=SC1090
set +u
source "$AMBER_SH"
set -u

echo "========================================="
echo "Processing Sequence : $SEQ_IN"
echo "Force Field         : $FF_CHOICE"
echo "Protein Leaprc      : $PROTEIN_LEAPRC"
echo "CpHMD Parm Dir      : $PHMD_DIR"
echo "Extra FF Params     : ${EXTRA_FF_PARAMS_PATH:-none}"
echo "Calculated Length   : $PEPLEN"
echo "Amber Formatted     : $PEPRSQ"
echo "N-Terminus Cap      : ${NCAP:-standard}"
echo "C-Terminus Cap      : ${CCAP:-standard}"
echo "Output Filename     : $OUTNAME"
echo "========================================="

{
    echo "source ${PROTEIN_LEAPRC}"
    echo "loadoff ${PHMD_LIB}"
    echo "loadamberparams ${PHMD_FRCMOD}"

    if [[ -n "$EXTRA_FF_PARAMS_PATH" ]]; then
        echo "loadamberparams ${EXTRA_FF_PARAMS_PATH}"
    fi

    echo
    echo "pep = sequence { ${LEAP_SEQ} }"
    echo
    echo "check pep"
    echo "charge pep"
    echo
    echo "saveamberparm pep ${OUTNAME}.parm7 ${OUTNAME}.rst7"
    echo "savepdb pep ${OUTNAME}.pdb"
    echo
    echo "quit"
} > "$LEAPIN"

tleap -f "$LEAPIN" > "$LOGFILE"

if [[ ! -f "${OUTNAME}.parm7" ]]; then
    echo "Error: failed to write ${OUTNAME}.parm7"
    exit 1
fi

if [[ ! -f "${OUTNAME}.rst7" ]]; then
    echo "Error: failed to write ${OUTNAME}.rst7"
    exit 1
fi

if [[ ! -f "${OUTNAME}.pdb" ]]; then
    echo "Error: failed to write ${OUTNAME}.pdb"
    exit 1
fi

echo "Execution finished. Check $LOGFILE for logs."
