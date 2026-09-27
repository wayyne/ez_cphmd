#!/bin/bash

SCRIPT_DIR="$(
    cd -- "$(dirname -- "${BASH_SOURCE[0]}")"
    pwd -P
)"

CHARMM_INPUT="${SCRIPT_DIR}/genpep.inp"

RAW_SEQ=""
CAP_CHOICE=""
CAPPED_FLAG=0
UNCAPPED_FLAG=0

while [[ $# -gt 0 ]]; do
    case "$1" in
        --pep)
            if [ -n "$2" ] && [[ "$2" != --* ]]; then
                RAW_SEQ="$2"
                shift 2
            else
                echo "Error: --pep requires a sequence argument."
                exit 1
            fi
            ;;
        --capped)
            CAPPED_FLAG=1
            CAP_CHOICE="1"
            shift
            ;;
        --uncapped)
            UNCAPPED_FLAG=1
            CAP_CHOICE="2"
            shift
            ;;
        *)
            echo "Error: Unknown argument '$1'"
            echo "Usage: $0 [--pep <SEQ>] [--capped | --uncapped]"
            exit 1
            ;;
    esac
done

if [ $CAPPED_FLAG -eq 1 ] && [ $UNCAPPED_FLAG -eq 1 ]; then
    echo "Error: Cannot pass both --capped and --uncapped simultaneously."
    exit 1
fi

if [ -z "$RAW_SEQ" ]; then
    read -p "Enter the single-letter peptide sequence (e.g., GGDGG): " RAW_SEQ
fi
if [ -z "$RAW_SEQ" ]; then
    echo "Error: Sequence cannot be empty."
    exit 1
fi

SEQ_IN=$(echo "$RAW_SEQ" | tr '[:lower:]' '[:upper:]')
PEPLEN=${#SEQ_IN}

if [ -z "$CAP_CHOICE" ]; then
    echo "-----------------------------------------"
    echo "Select Capping Style:"
    echo " 1) Capped (Neutral: ACE / CT2)"
    echo " 2) Uncapped (Zwitterionic: NTER / CTER)"
    echo "-----------------------------------------"
    read -p "Enter choice [1 or 2]: " CAP_CHOICE
fi

if [ "$CAP_CHOICE" == "1" ]; then
    NCAP="ACE"
    CCAP="CT2"
    CAP_LABEL="capped"
elif [ "$CAP_CHOICE" == "2" ]; then
    NCAP="NTER"
    CCAP="CTER"
    CAP_LABEL="zwitter"
else
    echo "Error: Invalid selection."
    exit 1
fi

declare -A AA_MAP=(
    [A]="ALA" [R]="ARG" [N]="ASN" [D]="ASP" [C]="CYS"
    [Q]="GLN" [E]="GLU" [G]="GLY" [H]="HSP" [I]="ILE"
    [L]="LEU" [K]="LYS" [M]="MET" [F]="PHE" [P]="PRO"
    [S]="SER" [T]="THR" [W]="TRP" [Y]="TYR" [V]="VAL"
)

PEPRSQ=""
for (( i=0; i<$PEPLEN; i++ )); do
    CHAR="${SEQ_IN:$i:1}"
    THREE_LETTER="${AA_MAP[$CHAR]}"
    if [ -z "$THREE_LETTER" ]; then
        echo "Error: Unknown amino acid character '$CHAR' in sequence."
        exit 1
    fi
    if [ -z "$PEPRSQ" ]; then
        PEPRSQ="$THREE_LETTER"
    else
        PEPRSQ="$PEPRSQ $THREE_LETTER"
    fi
done

OUTNAME="$(printf '%s_%s' "$CAP_LABEL" "$SEQ_IN" | tr '[:upper:]' '[:lower:]')"

echo "========================================="
echo "Processing Sequence : $SEQ_IN"
echo "Calculated Length   : $PEPLEN"
echo "CHARMM Formatted    : $PEPRSQ"
echo "N-Terminus Patch    : $NCAP"
echo "C-Terminus Patch    : $CCAP"
echo "Output Filename     : $OUTNAME"
echo "========================================="

if [[ -z "${C22_TOPPAR:-}" ]]; then
    echo "Error: C22_TOPPAR is not set. Check prep/config.sh."
    exit 1
fi

LOGFILE="gen_${OUTNAME}.out"

charmm TOPPAR="$C22_TOPPAR" PEPLEN="$PEPLEN" PEPRSQ="$PEPRSQ" NCAP="$NCAP" CCAP="$CCAP" OUTNAME="$OUTNAME" < "$CHARMM_INPUT" > "$LOGFILE" 2>&1
CHARMM_STATUS=$?

if [[ "$CHARMM_STATUS" -ne 0 ]]; then
    echo "Error: CHARMM exited with status ${CHARMM_STATUS}." >&2
    echo "See log: $LOGFILE" >&2
    exit "$CHARMM_STATUS"
fi

for output_file in "${OUTNAME}.psf" "${OUTNAME}.crd" "${OUTNAME}.pdb"; do
    if [[ ! -s "$output_file" ]]; then
        echo "Error: expected CHARMM output was not created: $output_file" >&2
        echo "See log: $LOGFILE" >&2
        exit 1
    fi
done

echo "CHARMM peptide preparation completed successfully."
echo "PSF         : ${OUTNAME}.psf"
echo "Coordinates : ${OUTNAME}.crd"
echo "PDB         : ${OUTNAME}.pdb"
echo "CHARMM log  : $LOGFILE"
