#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# BUILD A DRY AMBER CpHMD SYSTEM FROM A CLEANED SINGLE-CHAIN PDB
#
# Capped:
#   ACE + protein core + NHE
#   Explicit ACE-C/N and C/NHE-N bonds
#
# Uncapped:
#   Protein core only
#
# Both modes use loadPdbUsingSeq for strict residue-template mapping.
# A restrained vacuum minimization repairs generated hydrogens, missing atoms,
# and cap coordinates while restraining protein heavy atoms.
###############################################################################

AMBER_SH="${AMBER_ENV:-}"

INFILE=""
OUTNAME=""
DISU_FILE=""
FF_CHOICE="ff19sb"
CAPPED_FLAG="0"
UNCAPPED_FLAG="0"
KEEP_TEMP="0"

usage() {
    cat <<EOF
Usage:
  $0 \
    --in <clean.pdb> \
    --out <output_prefix> \
    --disu <disulfides.dat> \
    --ff <ff19sb|ff14sb> \
    [--capped | --uncapped] \
    [--keep-temp]

Required arguments:
  --in <file>          Cleaned PDB from prep_pdb_amber.sh
  --out <prefix>       Output prefix without extension
  --disu <file>        Disulfide-pair file from prep_pdb_amber.sh
  --ff <name>          ff19sb or ff14sb

Terminal options:
  --capped             Build ACE/protein/NHE
  --uncapped           Build the protein core without caps

Optional:
  --keep-temp          Preserve pre-minimization files
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
        --ff|--forcefield)
            require_argument "$1" "${2:-}"
            FF_CHOICE="$(printf "%s" "$2" | tr '[:upper:]' '[:lower:]')"
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
        --keep-temp)
            KEEP_TEMP="1"
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

[[ -n "$INFILE" ]] || die "--in is required"
[[ -n "$OUTNAME" ]] || die "--out is required"
[[ -n "$DISU_FILE" ]] || die "--disu is required"

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
        die "--ff must be ff19sb or ff14sb"
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

[[ -f "$INFILE" ]] || die "cleaned PDB not found: $INFILE"
[[ -s "$INFILE" ]] || die "cleaned PDB is empty: $INFILE"
[[ -f "$DISU_FILE" ]] || die "disulfide file not found: $DISU_FILE"

if [[ "$CAPPED_FLAG" == "1" && "$UNCAPPED_FLAG" == "1" ]]; then
    die "cannot use --capped and --uncapped together"
fi

if [[ "$CAPPED_FLAG" == "0" && "$UNCAPPED_FLAG" == "0" ]]; then
    die "one of --capped or --uncapped is required"
fi

for path_value in "$INFILE" "$OUTNAME" "$DISU_FILE"; do
    if [[ "$path_value" =~ [[:space:]] ]]; then
        die "paths and output prefixes cannot contain whitespace: $path_value"
    fi
done

OUTDIR="$(dirname -- "$OUTNAME")"
OUTBASE="$(basename -- "$OUTNAME")"

[[ -d "$OUTDIR" ]] || die "output directory does not exist: $OUTDIR"
[[ -n "$OUTBASE" && "$OUTBASE" != "." && "$OUTBASE" != ".." ]] ||
    die "invalid output prefix: $OUTNAME"

[[ -n "$PHMD_DIR" ]] || die "CpHMD parameter directory is not set for $FF_CHOICE. Check prep/config.sh"
[[ -f "$AMBER_SH" ]] || die "missing Amber environment: $AMBER_SH"
[[ -f "$PHMD_LIB" ]] || die "missing CpHMD library: $PHMD_LIB"
[[ -f "$PHMD_FRCMOD" ]] || die "missing CpHMD frcmod: $PHMD_FRCMOD"

if [[ -n "$EXTRA_FF_PARAMS_PATH" && ! -f "$EXTRA_FF_PARAMS_PATH" ]]; then
    die "missing force-field parameter file: $EXTRA_FF_PARAMS_PATH"
fi

# shellcheck disable=SC1090
set +u
source "$AMBER_SH"
set -u

for required_command in tleap sander cpptraj python3; do
    command -v "$required_command" >/dev/null 2>&1 ||
        die "required command not found in PATH: $required_command"
done

TMPDIR="$(mktemp -d "${TMPDIR:-/tmp}/gen_pdb_amber.XXXXXX")"

cleanup() {
    if [[ "$KEEP_TEMP" == "1" ]]; then
        echo "Temporary files preserved in: $TMPDIR"
    else
        rm -rf -- "$TMPDIR"
    fi
}
trap cleanup EXIT

RESIDUE_TABLE="${TMPDIR}/residues.dat"
DISU_TABLE="${TMPDIR}/disulfides.dat"

python3 - "$INFILE" "$RESIDUE_TABLE" <<'PY'
from __future__ import annotations

import sys
from collections import OrderedDict
from pathlib import Path

infile = Path(sys.argv[1])
outfile = Path(sys.argv[2])

allowed = {
    "ALA",
    "ARG",
    "ASN",
    "AS2",
    "CYS",
    "CYX",
    "GLN",
    "GL2",
    "GLY",
    "HIP",
    "ILE",
    "LEU",
    "LYS",
    "MET",
    "PHE",
    "PRO",
    "SER",
    "THR",
    "TRP",
    "TYR",
    "VAL",
}

residues: OrderedDict[int, str] = OrderedDict()
chains: set[str] = set()
atom_count = 0
hetatm_count = 0
model_count = 0

for line in infile.read_text(encoding="utf-8", errors="replace").splitlines():
    record = line[:6]

    if record == "MODEL ":
        model_count += 1
        continue

    if record == "HETATM":
        hetatm_count += 1
        continue

    if record != "ATOM  ":
        continue

    atom_count += 1

    chain = line[21:22]
    chains.add(chain)

    resid_text = line[22:26].strip()
    resname = line[17:20].strip().upper()

    try:
        resid = int(resid_text)
    except ValueError:
        raise SystemExit(f"Non-integer residue number: {resid_text}")

    if resid not in residues:
        residues[resid] = resname
    elif residues[resid] != resname:
        raise SystemExit(
            f"Conflicting residue names at residue {resid}: "
            f"{residues[resid]} and {resname}"
        )

if atom_count == 0:
    raise SystemExit("Input PDB contains no ATOM records")

if hetatm_count != 0:
    raise SystemExit(
        "Input PDB contains HETATM records; run prep_pdb_amber.sh first"
    )

if model_count != 0:
    raise SystemExit(
        "Input PDB contains MODEL records; run prep_pdb_amber.sh first"
    )

if len(chains) != 1:
    raise SystemExit(
        f"Input PDB must contain exactly one chain; found {len(chains)}"
    )

for expected, (resid, resname) in enumerate(residues.items(), start=1):
    if resid != expected:
        raise SystemExit(
            f"Residues must be numbered consecutively from 1; "
            f"expected {expected}, found {resid}"
        )

    if resname not in allowed:
        raise SystemExit(
            f"Unsupported or unnormalized residue {resname} at residue {resid}"
        )

with outfile.open("w", encoding="utf-8") as handle:
    for resid, resname in residues.items():
        handle.write(f"{resid} {resname}\n")
PY

NRES="$(wc -l < "$RESIDUE_TABLE" | tr -d '[:space:]')"
[[ "$NRES" -gt 0 ]] || die "failed to extract the protein sequence"

declare -A RESNAME_BY_RESID=()
declare -a SEQUENCE=()

while read -r RESID RESNAME; do
    [[ -n "$RESID" ]] || continue
    RESNAME_BY_RESID["$RESID"]="$RESNAME"
    SEQUENCE+=("$RESNAME")
done < "$RESIDUE_TABLE"

declare -A DISULFIDE_SEEN=()
declare -A DISULFIDE_RESIDUE_USED=()
DISULFIDE_COUNT="0"
LINE_NUMBER="0"

: > "$DISU_TABLE"

while IFS= read -r RAW_LINE || [[ -n "$RAW_LINE" ]]; do
    LINE_NUMBER=$((LINE_NUMBER + 1))
    LINE="${RAW_LINE%%#*}"
    LINE="${LINE#"${LINE%%[![:space:]]*}"}"
    LINE="${LINE%"${LINE##*[![:space:]]}"}"

    [[ -n "$LINE" ]] || continue

    read -r RESID1 RESID2 EXTRA <<< "$LINE"

    [[ -z "${EXTRA:-}" ]] ||
        die "invalid disulfide line ${LINE_NUMBER}; expected two residue numbers"

    [[ "$RESID1" =~ ^[0-9]+$ && "$RESID2" =~ ^[0-9]+$ ]] ||
        die "invalid disulfide line ${LINE_NUMBER}; residue numbers must be integers"

    [[ "$RESID1" -ge 1 && "$RESID1" -le "$NRES" ]] ||
        die "disulfide residue ${RESID1} is outside 1-${NRES}"

    [[ "$RESID2" -ge 1 && "$RESID2" -le "$NRES" ]] ||
        die "disulfide residue ${RESID2} is outside 1-${NRES}"

    [[ "$RESID1" -ne "$RESID2" ]] ||
        die "a disulfide cannot connect residue ${RESID1} to itself"

    [[ "${RESNAME_BY_RESID[$RESID1]}" == "CYX" ]] ||
        die "disulfide residue ${RESID1} is ${RESNAME_BY_RESID[$RESID1]}, not CYX"

    [[ "${RESNAME_BY_RESID[$RESID2]}" == "CYX" ]] ||
        die "disulfide residue ${RESID2} is ${RESNAME_BY_RESID[$RESID2]}, not CYX"

    if [[ "$RESID1" -lt "$RESID2" ]]; then
        FIRST_RESID="$RESID1"
        SECOND_RESID="$RESID2"
    else
        FIRST_RESID="$RESID2"
        SECOND_RESID="$RESID1"
    fi

    PAIR_KEY="${FIRST_RESID}:${SECOND_RESID}"

    [[ -z "${DISULFIDE_SEEN[$PAIR_KEY]+x}" ]] ||
        die "duplicate disulfide pair: ${FIRST_RESID}-${SECOND_RESID}"

    [[ -z "${DISULFIDE_RESIDUE_USED[$FIRST_RESID]+x}" ]] ||
        die "CYX residue ${FIRST_RESID} occurs in more than one disulfide"

    [[ -z "${DISULFIDE_RESIDUE_USED[$SECOND_RESID]+x}" ]] ||
        die "CYX residue ${SECOND_RESID} occurs in more than one disulfide"

    DISULFIDE_SEEN["$PAIR_KEY"]="1"
    DISULFIDE_RESIDUE_USED["$FIRST_RESID"]="1"
    DISULFIDE_RESIDUE_USED["$SECOND_RESID"]="1"

    printf '%s %s\n' "$FIRST_RESID" "$SECOND_RESID" >> "$DISU_TABLE"
    DISULFIDE_COUNT=$((DISULFIDE_COUNT + 1))
done < "$DISU_FILE"

for ((RESID = 1; RESID <= NRES; RESID++)); do
    if [[ "${RESNAME_BY_RESID[$RESID]}" == "CYX" ]] &&
       [[ -z "${DISULFIDE_RESIDUE_USED[$RESID]+x}" ]]; then
        die "CYX residue ${RESID} is not assigned to a disulfide pair"
    fi
done

INFILE_ABS="$(absolute_existing_file "$INFILE")"
OUTNAME_ABS="$(absolute_output_path "$OUTNAME")"

LEAPIN="${OUTDIR}/gen_${OUTBASE}.leap"
LEAPLOG="${OUTDIR}/gen_${OUTBASE}.out"
MININ="${OUTDIR}/min_${OUTBASE}.in"
MINOUT="${OUTDIR}/min_${OUTBASE}.out"
CPPTRAJLOG="${OUTDIR}/capcheck_${OUTBASE}.out"
CAPDATA="${OUTDIR}/cap_bonds_${OUTBASE}.dat"

PREMIN_RST="${TMPDIR}/${OUTBASE}_premin.rst7"
PREMIN_PDB="${TMPDIR}/${OUTBASE}_premin.pdb"

SEQUENCE_TEXT="${SEQUENCE[*]}"

{
    echo "source ${PROTEIN_LEAPRC}"
    echo "loadoff ${PHMD_LIB}"
    echo "loadamberparams ${PHMD_FRCMOD}"

    if [[ -n "$EXTRA_FF_PARAMS_PATH" ]]; then
        echo "loadamberparams ${EXTRA_FF_PARAMS_PATH}"
    fi
    echo
    echo "coreseq = { ${SEQUENCE_TEXT} }"
    echo "core = loadPdbUsingSeq ${INFILE_ABS} coreseq"
    echo

    while read -r FIRST_RESID SECOND_RESID; do
        [[ -n "$FIRST_RESID" ]] || continue
        echo "bond core.${FIRST_RESID}.SG core.${SECOND_RESID}.SG"
    done < "$DISU_TABLE"

    if [[ "$CAPPED_FLAG" == "1" ]]; then
        LAST_CORE_SYSTEM_RESID=$((NRES + 1))
        NHE_RESID=$((NRES + 2))

        echo
        echo "ncap = sequence { ACE }"
        echo "ccap = sequence { NHE }"
        echo "system = combine { ncap core ccap }"
        echo "bond system.1.C system.2.N"
        echo "bond system.${LAST_CORE_SYSTEM_RESID}.C system.${NHE_RESID}.N"
    else
        echo
        echo "system = core"
    fi

    echo
    echo "check system"
    echo "charge system"
    echo
    echo "saveamberparm system ${OUTNAME_ABS}.parm7 ${PREMIN_RST}"
    echo "savepdb system ${PREMIN_PDB}"
    echo
    echo "quit"
} > "$LEAPIN"

if [[ "$CAPPED_FLAG" == "1" ]]; then
    TERMINI_LABEL="capped"
    TOTAL_RESIDUES=$((NRES + 2))
    RESTRAINT_MASK="!:1,${TOTAL_RESIDUES} & !@H="
else
    TERMINI_LABEL="uncapped"
    TOTAL_RESIDUES="$NRES"
    RESTRAINT_MASK="!@H="
fi

cat > "$MININ" <<EOF
Restrained dry-system minimization
&cntrl
  imin=1,
  maxcyc=5000,
  ncyc=1000,
  ntb=0,
  igb=6,
  cut=999.0,
  ntr=1,
  restraint_wt=50.0,
  restraintmask='${RESTRAINT_MASK}',
  ntpr=100,
  drms=0.0001,
/
EOF

echo "============================================================"
echo "Generating Amber protein system"
echo "Force field       : $FF_CHOICE"
echo "Protein Leaprc    : $PROTEIN_LEAPRC"
echo "CpHMD Parm Dir    : $PHMD_DIR"
echo "Extra FF Params   : ${EXTRA_FF_PARAMS_PATH:-none}"
echo "Input PDB          : $INFILE"
echo "Output prefix      : $OUTNAME"
echo "Protein residues   : $NRES"
echo "Total residues     : $TOTAL_RESIDUES"
echo "Terminal mode      : $TERMINI_LABEL"
echo "Disulfide count    : $DISULFIDE_COUNT"
echo "Min restraint mask : $RESTRAINT_MASK"
echo "LEaP input         : $LEAPIN"
echo "LEaP log           : $LEAPLOG"
echo "Minimization input : $MININ"
echo "Minimization log   : $MINOUT"
echo "============================================================"

set +e
tleap -f "$LEAPIN" > "$LEAPLOG" 2>&1
TLEAP_STATUS=$?
set -e

if [[ "$TLEAP_STATUS" -ne 0 ]]; then
    echo "Error: tleap exited with status ${TLEAP_STATUS}." >&2
    echo "See log: $LEAPLOG" >&2
    exit "$TLEAP_STATUS"
fi

if grep -Eqi \
    'FATAL:|does not have a type|Could not find bond parameter|Could not find angle parameter|Could not find torsion parameter|Exiting LEaP: Errors = [1-9]' \
    "$LEAPLOG"; then
    die "tleap reported a fatal topology error; see $LEAPLOG"
fi

[[ -s "${OUTNAME}.parm7" ]] ||
    die "tleap did not create ${OUTNAME}.parm7"

[[ -s "$PREMIN_RST" ]] ||
    die "tleap did not create the pre-minimization restart"

set +e
sander -O \
    -i "$MININ" \
    -o "$MINOUT" \
    -p "${OUTNAME}.parm7" \
    -c "$PREMIN_RST" \
    -ref "$PREMIN_RST" \
    -r "${OUTNAME}.rst7"
SANDER_STATUS=$?
set -e

if [[ "$SANDER_STATUS" -ne 0 ]]; then
    echo "Error: sander exited with status ${SANDER_STATUS}." >&2
    echo "See log: $MINOUT" >&2
    exit "$SANDER_STATUS"
fi

if grep -Eqi \
    'FATAL|ERROR|BOMB|NaN|vlimit exceeded|Could not open' \
    "$MINOUT"; then
    die "sander reported a minimization failure; see $MINOUT"
fi

[[ -s "${OUTNAME}.rst7" ]] ||
    die "sander did not create ${OUTNAME}.rst7"

PDBLOG="${OUTDIR}/pdb_${OUTBASE}.out"

rm -f -- "${OUTNAME}.pdb"

set +e
cpptraj -p "${OUTNAME}.parm7" > "$PDBLOG" 2>&1 <<EOF
trajin ${OUTNAME}.rst7
trajout ${OUTNAME}.pdb pdb
run
quit
EOF
CPPTRAJ_PDB_STATUS=$?
set -e

if [[ "$CPPTRAJ_PDB_STATUS" -ne 0 ]]; then
    echo "Error: cpptraj PDB conversion exited with status ${CPPTRAJ_PDB_STATUS}." >&2
    echo "See log: $PDBLOG" >&2
    exit "$CPPTRAJ_PDB_STATUS"
fi

[[ -s "${OUTNAME}.pdb" ]] ||
    die "cpptraj did not create ${OUTNAME}.pdb; see $PDBLOG"

if [[ "$CAPPED_FLAG" == "1" ]]; then
    LAST_CORE_SYSTEM_RESID=$((NRES + 1))
    NHE_RESID=$((NRES + 2))

    cpptraj > "$CPPTRAJLOG" 2>&1 <<EOF
parm ${OUTNAME}.parm7
trajin ${OUTNAME}.rst7
distance ncap_bond :1@C :2@N out ${CAPDATA}
distance ccap_bond :${LAST_CORE_SYSTEM_RESID}@C :${NHE_RESID}@N out ${CAPDATA}
run
quit
EOF

    [[ -s "$CAPDATA" ]] ||
        die "cpptraj did not create cap-bond validation data"

    read -r _FRAME NCAP_DIST CCAP_DIST < <(
        awk '!/^#/ && NF >= 3 {frame=$1; ncap=$2; ccap=$3}
             END {print frame, ncap, ccap}' "$CAPDATA"
    )

    [[ -n "${NCAP_DIST:-}" && -n "${CCAP_DIST:-}" ]] ||
        die "failed to parse minimized cap-bond distances"

    python3 - "$NCAP_DIST" "$CCAP_DIST" <<'PY'
import sys

ncap = float(sys.argv[1])
ccap = float(sys.argv[2])

for label, value in (("ACE-C/N", ncap), ("C/NHE-N", ccap)):
    if not 1.20 <= value <= 1.50:
        raise SystemExit(
            f"Error: minimized {label} bond is outside 1.20-1.50 A: {value:.4f}"
        )
PY

    echo "Minimized ACE C-N bond : ${NCAP_DIST} Å"
    echo "Minimized C-NHE bond   : ${CCAP_DIST} Å"
fi

if [[ "$KEEP_TEMP" == "0" ]]; then
    rm -f -- "$PREMIN_RST" "$PREMIN_PDB"
fi

echo
echo "Amber ${FF_CHOICE} protein preparation completed successfully."
echo "Topology         : ${OUTNAME}.parm7"
echo "Coordinates      : ${OUTNAME}.rst7"
echo "PDB              : ${OUTNAME}.pdb"
echo "LEaP log         : $LEAPLOG"
echo "Minimization log : $MINOUT"
echo "PDB conversion   : $PDBLOG"

if [[ "$CAPPED_FLAG" == "1" ]]; then
    echo "Cap validation   : $CAPDATA"
fi
