#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# GENERATE A CHARMM PSF/CRD FROM A CLEANED SINGLE-CHAIN PDB
###############################################################################

SCRIPT_DIR="$(
    cd -- "$(dirname -- "${BASH_SOURCE[0]}")"
    pwd -P
)"

CHARMM_INPUT="${SCRIPT_DIR}/gen_pdb_system.inp"

INFILE=""
OUTNAME=""
DISU_FILE=""

CAPPED_FLAG="0"
UNCAPPED_FLAG="0"

usage() {
    cat <<EOF
Usage:
  $0 \\
    --in <clean.pdb> \\
    --out <output_prefix> \\
    --disu <disulfides.dat> \\
    [--capped | --uncapped]

Required arguments:
  --in <file>          Cleaned single-chain PDB from prep_pdb.sh
  --out <prefix>       Output prefix without a file extension
  --disu <file>        Disulfide-pair file from prep_pdb.sh

Terminal options:
  --capped             Apply ACE/CT2 terminal patches
  --uncapped           Apply NTER/CTER terminal patches

Exactly one terminal option must be selected.

Example:
  $0 \\
    --in 1ubq_a_clean.pdb \\
    --out uncapped_1ubq_a \\
    --disu 1ubq_a_disulfides.dat \\
    --uncapped
EOF
}

die() {
    echo "Error: $*" >&2
    exit 1
}

require_argument() {
    local option="$1"
    local value="${2:-}"

    [[ -n "$value" ]] || die "$option requires an argument"
}

absolute_existing_file() {
    local input_path="$1"
    local input_dir
    local input_base

    input_dir="$(dirname -- "$input_path")"
    input_base="$(basename -- "$input_path")"

    (
        cd -- "$input_dir"
        printf '%s/%s\n' "$(pwd -P)" "$input_base"
    )
}

absolute_output_path() {
    local output_path="$1"
    local output_dir
    local output_base

    output_dir="$(dirname -- "$output_path")"
    output_base="$(basename -- "$output_path")"

    (
        cd -- "$output_dir"
        printf '%s/%s\n' "$(pwd -P)" "$output_base"
    )
}

###############################################################################
# ARGUMENT PARSING
###############################################################################

while [[ $# -gt 0 ]]; do
    case "$1" in
        --in)
            require_argument "$1" "${2:-}"
            INFILE="$2"
            shift 2
            ;;

        --out)
            require_argument "$1" "${2:-}"
            OUTNAME="$2"
            shift 2
            ;;

        --disu)
            require_argument "$1" "${2:-}"
            DISU_FILE="$2"
            shift 2
            ;;

        --capped)
            CAPPED_FLAG="1"
            shift
            ;;

        --uncapped)
            UNCAPPED_FLAG="1"
            shift
            ;;

        -h|--help)
            usage
            exit 0
            ;;

        *)
            die "unknown option: $1"
            ;;
    esac
done

###############################################################################
# VALIDATION
###############################################################################

[[ -n "$INFILE" ]] || die "--in is required"
[[ -n "$OUTNAME" ]] || die "--out is required"
[[ -n "$DISU_FILE" ]] || die "--disu is required"

[[ -f "$INFILE" ]] || die "cleaned PDB not found: $INFILE"
[[ -s "$INFILE" ]] || die "cleaned PDB is empty: $INFILE"

[[ -f "$DISU_FILE" ]] || die "disulfide file not found: $DISU_FILE"
[[ -f "$CHARMM_INPUT" ]] || die "CHARMM input script not found: $CHARMM_INPUT"
[[ -n "${C22_TOPPAR:-}" ]] || die "C22_TOPPAR is not set. Check prep/config.sh"

if [[ "$CAPPED_FLAG" == "1" && "$UNCAPPED_FLAG" == "1" ]]; then
    die "cannot use --capped and --uncapped together"
fi

if [[ "$CAPPED_FLAG" == "0" && "$UNCAPPED_FLAG" == "0" ]]; then
    die "one of --capped or --uncapped is required"
fi

if ! command -v charmm >/dev/null 2>&1; then
    die "CHARMM executable was not found in PATH"
fi

if [[ "$INFILE" =~ [[:space:]] ]]; then
    die "input PDB path cannot contain whitespace"
fi

if [[ "$OUTNAME" =~ [[:space:]] ]]; then
    die "output prefix cannot contain whitespace"
fi

if [[ "$DISU_FILE" =~ [[:space:]] ]]; then
    die "disulfide-file path cannot contain whitespace"
fi

OUTDIR="$(dirname -- "$OUTNAME")"
OUTBASE="$(basename -- "$OUTNAME")"

[[ -d "$OUTDIR" ]] || die "output directory does not exist: $OUTDIR"
[[ -n "$OUTBASE" ]] || die "invalid output prefix: $OUTNAME"

if [[ "$OUTBASE" == "." || "$OUTBASE" == ".." ]]; then
    die "invalid output prefix: $OUTNAME"
fi

if [[ "$CAPPED_FLAG" == "1" ]]; then
    NCAP="ACE"
    CCAP="CT2"
    TERMINI_LABEL="capped"
else
    NCAP="NTER"
    CCAP="CTER"
    TERMINI_LABEL="uncapped"
fi

###############################################################################
# TEMPORARY WORKSPACE
###############################################################################

TMPDIR="$(mktemp -d "${TMPDIR:-/tmp}/gen_pdb_system.XXXXXX")"

cleanup() {
    rm -rf -- "$TMPDIR"
}

trap cleanup EXIT

RESIDUE_TABLE="${TMPDIR}/residues.dat"

###############################################################################
# VERIFY CLEANED-PDB CONTENT
###############################################################################

ATOM_COUNT="$(
    awk '
        substr($0, 1, 6) == "ATOM  " {
            count++
        }
        END {
            print count + 0
        }
    ' "$INFILE"
)"

[[ "$ATOM_COUNT" -gt 0 ]] || die "input PDB contains no ATOM records"

HETATM_COUNT="$(
    awk '
        substr($0, 1, 6) == "HETATM" {
            count++
        }
        END {
            print count + 0
        }
    ' "$INFILE"
)"

if [[ "$HETATM_COUNT" -ne 0 ]]; then
    die "input PDB contains HETATM records; run prep_pdb.sh first"
fi

MODEL_COUNT="$(
    awk '
        substr($0, 1, 6) == "MODEL " {
            count++
        }
        END {
            print count + 0
        }
    ' "$INFILE"
)"

if [[ "$MODEL_COUNT" -ne 0 ]]; then
    die "input PDB contains MODEL records; run prep_pdb.sh first"
fi

CHAIN_COUNT="$(
    awk '
        substr($0, 1, 6) == "ATOM  " {
            chain = substr($0, 22, 1)
            seen[chain] = 1
        }
        END {
            for (chain in seen) {
                count++
            }

            print count + 0
        }
    ' "$INFILE"
)"

if [[ "$CHAIN_COUNT" -ne 1 ]]; then
    die "input PDB must contain exactly one chain; found ${CHAIN_COUNT}"
fi

###############################################################################
# EXTRACT THE RESIDUE SEQUENCE
###############################################################################

LC_ALL=C awk '
function trim(value) {
    sub(/^[[:space:]]+/, "", value)
    sub(/[[:space:]]+$/, "", value)
    return value
}

substr($0, 1, 6) == "ATOM  " {
    resid = trim(substr($0, 23, 4))
    resname = toupper(trim(substr($0, 18, 3)))

    if (resid !~ /^[0-9]+$/) {
        print "Non-integer residue number: " resid > "/dev/stderr"
        exit 1
    }

    key = resid + 0

    if (!(key in seen)) {
        seen[key] = resname
        order[++count] = key
    } else if (seen[key] != resname) {
        print \
            "Conflicting residue names at residue " key ": " \
            seen[key] " and " resname \
            > "/dev/stderr"
        exit 1
    }
}

END {
    if (count == 0) {
        print "No residues were found" > "/dev/stderr"
        exit 1
    }

    for (i = 1; i <= count; i++) {
        key = order[i]
        print key, seen[key]
    }
}
' "$INFILE" > "$RESIDUE_TABLE"

declare -A ALLOWED_RESIDUES=(
    [ALA]=1
    [ARG]=1
    [ASN]=1
    [ASP]=1
    [CYS]=1
    [GLN]=1
    [GLU]=1
    [GLY]=1
    [HSP]=1
    [ILE]=1
    [LEU]=1
    [LYS]=1
    [MET]=1
    [PHE]=1
    [PRO]=1
    [SER]=1
    [THR]=1
    [TRP]=1
    [TYR]=1
    [VAL]=1
)

declare -A RESNAME_BY_RESID=()
declare -a SEQUENCE=()

EXPECTED_RESID="1"

while read -r RESID RESNAME; do
    [[ -n "$RESID" ]] || continue

    if [[ "$RESID" -ne "$EXPECTED_RESID" ]]; then
        die \
            "residues must be numbered consecutively from 1; expected " \
            "${EXPECTED_RESID}, found ${RESID}"
    fi

    if [[ -z "${ALLOWED_RESIDUES[$RESNAME]+x}" ]]; then
        die "unsupported or unnormalized residue ${RESNAME} at residue ${RESID}"
    fi

    RESNAME_BY_RESID["$RESID"]="$RESNAME"
    SEQUENCE+=("$RESNAME")

    EXPECTED_RESID=$((EXPECTED_RESID + 1))
done < "$RESIDUE_TABLE"

NRES="${#SEQUENCE[@]}"

[[ "$NRES" -gt 0 ]] || die "failed to extract the protein sequence"

###############################################################################
# CREATE THE CHARMM SEQUENCE STREAM
###############################################################################

SEQUENCE_STREAM="${OUTDIR}/sequence_${OUTBASE}.str"

{
    echo "read sequence card"
    echo "* Sequence generated from ${INFILE}"
    echo "*"
    echo "$NRES"

    for ((i = 0; i < NRES; i++)); do
        printf '%-4s' "${SEQUENCE[$i]}"

        if (( (i + 1) % 12 == 0 || i + 1 == NRES )); then
            printf '\n'
        else
            printf ' '
        fi
    done

    echo "RETURN"
} > "$SEQUENCE_STREAM"

###############################################################################
# VALIDATE AND CREATE THE CHARMM DISULFIDE STREAM
###############################################################################

DISULFIDE_STREAM="${OUTDIR}/disulfides_${OUTBASE}.str"

declare -A DISULFIDE_SEEN=()
DISULFIDE_COUNT="0"
LINE_NUMBER="0"

{
    echo "! DISU patches generated by gen_pdb_system.sh"
    echo "! Segment ID is supplied by gen_pdb_system.inp"
} > "$DISULFIDE_STREAM"

while IFS= read -r RAW_LINE || [[ -n "$RAW_LINE" ]]; do
    LINE_NUMBER=$((LINE_NUMBER + 1))

    LINE="${RAW_LINE%%#*}"

    LINE="${LINE#"${LINE%%[![:space:]]*}"}"
    LINE="${LINE%"${LINE##*[![:space:]]}"}"

    [[ -n "$LINE" ]] || continue

    read -r RESID1 RESID2 EXTRA <<< "$LINE"

    if [[ -n "${EXTRA:-}" ]]; then
        die \
            "invalid disulfide line ${LINE_NUMBER}; expected two residue numbers"
    fi

    if ! [[ "$RESID1" =~ ^[0-9]+$ && "$RESID2" =~ ^[0-9]+$ ]]; then
        die \
            "invalid disulfide line ${LINE_NUMBER}; residue numbers must be integers"
    fi

    if [[ "$RESID1" -lt 1 || "$RESID1" -gt "$NRES" ]]; then
        die \
            "disulfide residue ${RESID1} is outside the valid range 1-${NRES}"
    fi

    if [[ "$RESID2" -lt 1 || "$RESID2" -gt "$NRES" ]]; then
        die \
            "disulfide residue ${RESID2} is outside the valid range 1-${NRES}"
    fi

    if [[ "$RESID1" -eq "$RESID2" ]]; then
        die "a disulfide cannot connect residue ${RESID1} to itself"
    fi

    if [[ "${RESNAME_BY_RESID[$RESID1]}" != "CYS" ]]; then
        die \
            "disulfide residue ${RESID1} is " \
            "${RESNAME_BY_RESID[$RESID1]}, not CYS"
    fi

    if [[ "${RESNAME_BY_RESID[$RESID2]}" != "CYS" ]]; then
        die \
            "disulfide residue ${RESID2} is " \
            "${RESNAME_BY_RESID[$RESID2]}, not CYS"
    fi

    if [[ "$RESID1" -lt "$RESID2" ]]; then
        PAIR_KEY="${RESID1}:${RESID2}"
        FIRST_RESID="$RESID1"
        SECOND_RESID="$RESID2"
    else
        PAIR_KEY="${RESID2}:${RESID1}"
        FIRST_RESID="$RESID2"
        SECOND_RESID="$RESID1"
    fi

    if [[ -n "${DISULFIDE_SEEN[$PAIR_KEY]+x}" ]]; then
        die "duplicate disulfide pair: ${FIRST_RESID}-${SECOND_RESID}"
    fi

    DISULFIDE_SEEN["$PAIR_KEY"]="1"
    DISULFIDE_COUNT=$((DISULFIDE_COUNT + 1))

    printf \
        'patch DISU @SEGID %d @SEGID %d\n' \
        "$FIRST_RESID" \
        "$SECOND_RESID" \
        >> "$DISULFIDE_STREAM"

done < "$DISU_FILE"

if [[ "$DISULFIDE_COUNT" -eq 0 ]]; then
    echo "! No disulfide patches requested" >> "$DISULFIDE_STREAM"
fi

echo "RETURN" >> "$DISULFIDE_STREAM"

###############################################################################
# RESOLVE ABSOLUTE PATHS FOR CHARMM
###############################################################################

INFILE_ABS="$(absolute_existing_file "$INFILE")"
OUTNAME_ABS="$(absolute_output_path "$OUTNAME")"
SEQUENCE_STREAM_ABS="$(absolute_existing_file "$SEQUENCE_STREAM")"
DISULFIDE_STREAM_ABS="$(absolute_existing_file "$DISULFIDE_STREAM")"

LOGFILE="${OUTDIR}/gen_${OUTBASE}.out"

###############################################################################
# DISPLAY SETTINGS
###############################################################################

echo "============================================================"
echo "Generating CHARMM protein system"
echo "Input PDB        : $INFILE"
echo "Output prefix    : $OUTNAME"
echo "Residue count    : $NRES"
echo "Terminal mode    : $TERMINI_LABEL"
echo "N-terminal patch : $NCAP"
echo "C-terminal patch : $CCAP"
echo "Disulfide count  : $DISULFIDE_COUNT"
echo "Sequence stream  : $SEQUENCE_STREAM"
echo "Disulfide stream : $DISULFIDE_STREAM"
echo "CHARMM log       : $LOGFILE"
echo "============================================================"

###############################################################################
# RUN CHARMM
###############################################################################

set +e

charmm \
    TOPPAR="$C22_TOPPAR" \
    INPDB="$INFILE_ABS" \
    OUTNAME="$OUTNAME_ABS" \
    NCAP="$NCAP" \
    CCAP="$CCAP" \
    SEQFILE="$SEQUENCE_STREAM_ABS" \
    DISUFILE="$DISULFIDE_STREAM_ABS" \
    NRES="$NRES" \
    < "$CHARMM_INPUT" \
    > "$LOGFILE" 2>&1

CHARMM_STATUS=$?

set -e

if [[ "$CHARMM_STATUS" -ne 0 ]]; then
    echo "Error: CHARMM exited with status ${CHARMM_STATUS}." >&2
    echo "See log: $LOGFILE" >&2
    exit "$CHARMM_STATUS"
fi

###############################################################################
# OUTPUT VALIDATION
###############################################################################

if grep -Eqi \
    'ABNORMAL TERMINATION|terminated due to the detection of a fatal error|Fortran runtime error' \
    "$LOGFILE"; then
    die "CHARMM reported a fatal error; see $LOGFILE"
fi

EXPECTED_OUTPUTS=(
    "${OUTNAME}.psf"
    "${OUTNAME}.crd"
    "${OUTNAME}.pdb"
)

for output_file in "${EXPECTED_OUTPUTS[@]}"; do
    if [[ ! -s "$output_file" ]]; then
        die "expected CHARMM output was not created: $output_file"
    fi
done

echo
echo "CHARMM protein preparation completed successfully."
echo "PSF             : ${OUTNAME}.psf"
echo "Coordinates     : ${OUTNAME}.crd"
echo "PDB             : ${OUTNAME}.pdb"
echo "CHARMM log      : $LOGFILE"
echo "Sequence stream : $SEQUENCE_STREAM"
echo "DISU stream     : $DISULFIDE_STREAM"
