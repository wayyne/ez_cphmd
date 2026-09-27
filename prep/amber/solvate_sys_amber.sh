#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# Amber24 CpHMD solvation and ion builder.
#
# ff19sb : ff19SB + OPC water
# ff14sb : ff14SB + TIP3P water + ff14sb_pme_dayhoff_dev.parm
#
# EXACT mode:
#   Passing --sod and/or --cla uses exactly those bulk ions.
#   No neutralization is performed.
#
# AUTO mode:
#   Omit --sod/--cla.
#   Use --conc in mM.
#   Salt pairs are computed from water count:
#
#       npairs = int(nwater * conc_molar / 55.5)
#
#   Then charge compensation is added.
#
# Box handling:
#   --box-edge <A> fits the Amber LEaP cushion until the final cell
#   edge is close to the requested value.
###############################################################################

AMBER_SH="${AMBER_ENV:-}"

INNAME=""
FF_CHOICE="ff19sb"
CUSHION=""
BOX_TYPE=""
BOX_EDGE=""
CONC_MM="150"
TITR_FLAG=""

ION_MODE="AUTO"
SOD_EXACT=""
CLA_EXACT=""
EXACT_IONS_REQUESTED="0"

KEEP_TEMP="0"

FIT_TOL="0.05"
FIT_MAX="10"

usage()
{
    echo "Usage: $0 --in <name> --ff <ff19sb|ff14sb> --box <type> --cushion <A> [opts]"
    echo
    echo "Box types:"
    echo "  cubic"
    echo "  octahedral"
    echo
    echo "Options:"
    echo "  --ff <ff19sb|ff14sb>"
    echo "  --box-edge <A>"
    echo "  --conc <mM>"
    echo "  --sod <N>"
    echo "  --cla <N>"
    echo "  --titr <-1|0|1>"
    echo "  --keep-temp"
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --in)
            INNAME="${2:-}"
            shift 2
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
        --cushion)
            CUSHION="${2:-}"
            shift 2
            ;;
        --box)
            BOX_TYPE="$(printf "%s" "${2:-}" | tr '[:upper:]' '[:lower:]')"
            shift 2
            ;;
        --box-edge|--edge)
            BOX_EDGE="${2:-}"
            shift 2
            ;;
        --conc|--conc-mm|--salt-mm)
            CONC_MM="${2:-}"
            shift 2
            ;;
        --sod)
            SOD_EXACT="${2:-}"
            EXACT_IONS_REQUESTED="1"
            shift 2
            ;;
        --cla)
            CLA_EXACT="${2:-}"
            EXACT_IONS_REQUESTED="1"
            shift 2
            ;;
        --titr)
            TITR_FLAG="${2:-}"
            shift 2
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
            echo "Unknown option: $1"
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
        WATER_LEAPRC="leaprc.water.opc"
        WATER_BOX="OPCBOX"
        EXTRA_FF_PARAMS=""
        ;;
    ff14sb|ff14)
        FF_CHOICE="ff14sb"
        FF_DIR="ff14"
        PROTEIN_LEAPRC="leaprc.protein.ff14SB"
        WATER_LEAPRC="leaprc.water.tip3p"
        WATER_BOX="TIP3PBOX"
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

if [[ -z "$INNAME" ]]; then
    read -r -p "Enter input base name: " INNAME
fi

if [[ -z "$BOX_TYPE" ]]; then
    echo "-----------------------------------------"
    echo "Select Unit Cell Geometry:"
    echo " 1) Cubic"
    echo " 2) Truncated Octahedral"
    echo "-----------------------------------------"
    read -r -p "Enter choice [1 or 2]: " BOX_TYPE
    BOX_TYPE="$(printf "%s" "$BOX_TYPE" | tr '[:upper:]' '[:lower:]')"
fi

if [[ -z "$CUSHION" ]]; then
    read -r -p "Enter water cushion size in Angstroms: " CUSHION
fi

if [[ -z "$TITR_FLAG" ]]; then
    echo "-----------------------------------------"
    echo "Enable titratable co-ions?"
    echo " 0) No"
    echo " -1) Acidic"
    echo " 1) Basic"
    echo "-----------------------------------------"
    read -r -p "Enter choice [0, -1, or 1]: " TITR_FLAG
fi

if [[ -z "$INNAME" ]]; then
    echo "Error: --in cannot be empty."
    exit 1
fi

if [[ ! -f "${INNAME}.pdb" ]]; then
    echo "Error: missing input PDB: ${INNAME}.pdb"
    exit 1
fi

if [[ ! -f "${INNAME}.parm7" ]]; then
    echo "Error: missing input parm7: ${INNAME}.parm7"
    exit 1
fi

if [[ ! -f "${INNAME}.rst7" ]]; then
    echo "Error: missing input rst7: ${INNAME}.rst7"
    exit 1
fi

if ! [[ "$CUSHION" =~ ^[0-9]+([.][0-9]+)?$ ]]; then
    echo "Error: --cushion must be a non-negative number."
    exit 1
fi

if [[ -n "$BOX_EDGE" ]]; then
    if ! [[ "$BOX_EDGE" =~ ^[0-9]+([.][0-9]+)?$ ]]; then
        echo "Error: --box-edge must be a positive number."
        exit 1
    fi
fi

if ! [[ "$CONC_MM" =~ ^[0-9]+([.][0-9]+)?$ ]]; then
    echo "Error: --conc must be a non-negative mM value."
    exit 1
fi

if [[ -n "$SOD_EXACT" ]] && ! [[ "$SOD_EXACT" =~ ^[0-9]+$ ]]; then
    echo "Error: --sod must be a non-negative integer."
    exit 1
fi

if [[ -n "$CLA_EXACT" ]] && ! [[ "$CLA_EXACT" =~ ^[0-9]+$ ]]; then
    echo "Error: --cla must be a non-negative integer."
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
    -1|0|1)
        ;;
    *)
        echo "Error: --titr must be one of -1, 0, or 1."
        exit 1
        ;;
esac

case "$BOX_TYPE" in
    cubic|cube|cub|1)
        BOX_TYPE="cubic"
        ;;
    octahedral|octa|oct|truncated_octahedral|2)
        BOX_TYPE="octahedral"
        ;;
    rhombic|rhdo|dodecahedron|rhombic_dodecahedron|3)
        echo "Error: rhombic dodecahedron is not supported here."
        echo "       LEaP supports cubic and octahedral simply."
        exit 1
        ;;
    *)
        echo "Error: unknown box type '$BOX_TYPE'"
        usage
        exit 1
        ;;
esac

sanitize_num_tag()
{
    local value="$1"

    value="$(printf "%s" "$value" | sed -E 's/\.0+$//')"
    value="${value//./p}"

    printf "%s" "$value"
}

format_mM_tag()
{
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

if [[ -n "$BOX_EDGE" ]]; then
    EDGE_TAG="$(sanitize_num_tag "$BOX_EDGE")"
    OUTNAME="${INNAME}_solv_${EDGE_TAG}Aedge_${BOX_TYPE}_${ION_TAG}"
fi

LOGFILE="solvate_${OUTNAME}.out"

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

TMPDIR="amber_solv_tmp_${OUTNAME}"

if [[ -d "$TMPDIR" ]]; then
    rm -rf "$TMPDIR"
fi

mkdir "$TMPDIR"

cleanup()
{
    if [[ "$KEEP_TEMP" == "0" ]]; then
        rm -rf "$TMPDIR"
    else
        echo "Temporary files preserved in $TMPDIR"
    fi
}

trap cleanup EXIT

write_common_leap()
{
    echo "source ${PROTEIN_LEAPRC}"
    echo "source ${WATER_LEAPRC}"
    echo "loadoff ${PHMD_LIB}"
    echo "loadamberparams ${PHMD_FRCMOD}"

    if [[ -n "$EXTRA_FF_PARAMS_PATH" ]]; then
        echo "loadamberparams ${EXTRA_FF_PARAMS_PATH}"
    fi

    echo
    echo "pep = loadpdb ${INNAME}.pdb"
}

write_solvate_cmd()
{
    local buf="$1"

    if [[ "$BOX_TYPE" == "octahedral" ]]; then
        echo "solvateOct pep ${WATER_BOX} ${buf}"
    else
        echo "solvateBox pep ${WATER_BOX} ${buf}"
    fi
}

get_prmtop_info()
{
    python3 - "${INNAME}.parm7" "$TITR_FLAG" <<'PY'
import sys

parm = sys.argv[1]
titr = int(sys.argv[2])

labels = []
charges = []
reading_charge = False
reading_label = False

with open(parm, "r", encoding="utf-8", errors="replace") as fh:
    for line in fh:
        if line.startswith("%FLAG CHARGE"):
            reading_charge = True
            reading_label = False
            continue

        if line.startswith("%FLAG RESIDUE_LABEL"):
            reading_label = True
            reading_charge = False
            continue

        if line.startswith("%FLAG"):
            reading_charge = False
            reading_label = False
            continue

        if reading_charge and line.startswith("%FORMAT"):
            continue

        if reading_label and line.startswith("%FORMAT"):
            continue

        if reading_charge:
            charges.extend(float(x) for x in line.split())

        if reading_label:
            row = line.rstrip("\n")
            for i in range(0, len(row), 4):
                lab = row[i:i + 4].strip()
                if lab:
                    labels.append(lab)

charge = sum(charges) / 18.2223
netchg = int(charge + 0.1) if charge >= 0.0 else int(charge - 0.1)

nacid = sum(1 for x in labels if x in ("AS2", "GL2", "ASP", "GLU"))
nbase = sum(1 for x in labels if x in ("LYS", "ARG"))
nhis = sum(1 for x in labels if x == "HIP")

nposco = 0
nnegco = 0

if titr < 0:
    nposco = nhis
    nnegco = nacid
elif titr > 0:
    nposco = nbase + nhis
    nnegco = 0

print(netchg, nposco, nnegco)
PY
}

get_box_edge()
{
    python3 - "$1" <<'PY'
import sys

rst = sys.argv[1]

with open(rst, "r", encoding="utf-8", errors="replace") as fh:
    lines = [x.strip() for x in fh if x.strip()]

vals = [float(x) for x in lines[-1].split()]

if len(vals) < 6:
    sys.exit("Error: restart file lacks box information.")

print("%.8f" % vals[0])
PY
}

count_waters()
{
    python3 - "$1" <<'PY'
import sys

parm = sys.argv[1]
labels = []
reading = False

with open(parm, "r", encoding="utf-8", errors="replace") as fh:
    for line in fh:
        if line.startswith("%FLAG RESIDUE_LABEL"):
            reading = True
            continue

        if reading and line.startswith("%FORMAT"):
            continue

        if reading and line.startswith("%FLAG"):
            break

        if reading:
            row = line.rstrip("\n")
            for i in range(0, len(row), 4):
                lab = row[i:i + 4].strip()
                if lab:
                    labels.append(lab)

print(labels.count("WAT"))
PY
}

run_pre_solvate()
{
    local buf="$1"
    local tag="$2"
    local leap="${TMPDIR}/${tag}.leap"
    local parm="${TMPDIR}/${tag}.parm7"
    local rst="${TMPDIR}/${tag}.rst7"
    local log="${TMPDIR}/${tag}.log"

    {
        write_common_leap
        echo
        write_solvate_cmd "$buf"
        echo
        echo "saveamberparm pep ${parm} ${rst}"
        echo "quit"
    } > "$leap"

    tleap -f "$leap" > "$log"

    if [[ ! -f "$rst" ]]; then
        echo "Error: failed pre-solvation with buffer $buf"
        exit 1
    fi
}

fit_cushion_to_edge()
{
    python3 - "$CUSHION" "$BOX_EDGE" <<'PY'
import sys

start = float(sys.argv[1])
target = float(sys.argv[2])

if start <= 0.0:
    start = target / 3.0

print("%.6f" % start)
PY
}

update_cushion_guess()
{
    python3 - "$1" "$2" "$BOX_EDGE" <<'PY'
import sys

buf = float(sys.argv[1])
edge = float(sys.argv[2])
target = float(sys.argv[3])

newbuf = buf + (target - edge)

if newbuf < 0.1:
    newbuf = 0.1

print("%.6f" % newbuf)
PY
}

edge_abs_error()
{
    python3 - "$1" "$BOX_EDGE" <<'PY'
import sys

edge = float(sys.argv[1])
target = float(sys.argv[2])

print("%.8f" % abs(edge - target))
PY
}

read -r NETCHG NPOSCO NNEGCO < <(get_prmtop_info)

FINAL_CUSHION="$CUSHION"
FIT_EDGE=""

if [[ -n "$BOX_EDGE" ]]; then
    guess="$(fit_cushion_to_edge)"
    best_buf="$guess"
    best_edge="0.0"
    best_err="999999.0"

    for iter in $(seq 1 "$FIT_MAX"); do
        tag="fit_${iter}"
        run_pre_solvate "$guess" "$tag"

        edge="$(get_box_edge "${TMPDIR}/${tag}.rst7")"
        err="$(edge_abs_error "$edge")"

        better="$(awk -v e="$err" -v b="$best_err" \
            'BEGIN { if (e < b) print 1; else print 0 }')"

        if [[ "$better" -eq 1 ]]; then
            best_buf="$guess"
            best_edge="$edge"
            best_err="$err"
        fi

        ok="$(awk -v e="$err" -v t="$FIT_TOL" \
            'BEGIN { if (e <= t) print 1; else print 0 }')"

        if [[ "$ok" -eq 1 ]]; then
            break
        fi

        guess="$(update_cushion_guess "$guess" "$edge")"
    done

    FINAL_CUSHION="$best_buf"
    FIT_EDGE="$best_edge"
fi

PREPARM="${TMPDIR}/preion.parm7"
PRERST="${TMPDIR}/preion.rst7"

run_pre_solvate "$FINAL_CUSHION" "preion"

NWAT="$(count_waters "$PREPARM")"
FINAL_EDGE="$(get_box_edge "$PRERST")"

calc_ion_counts()
{
    python3 - "$ION_MODE" "$CONC_MM" "$NWAT" \
        "$NETCHG" "$NPOSCO" "$NNEGCO" \
        "$SOD_EXACT" "$CLA_EXACT" <<'PY'
import sys

mode = sys.argv[1]
conc_mm = float(sys.argv[2])
nwater = int(sys.argv[3])
netchg = int(sys.argv[4])
nposco = int(sys.argv[5])
nnegco = int(sys.argv[6])
sod_exact = int(sys.argv[7])
cla_exact = int(sys.argv[8])

if mode == "EXACT":
    npos = sod_exact
    nneg = cla_exact
    npairs = 0
else:
    conc_m = conc_mm / 1000.0
    npairs = int(nwater * conc_m / 55.5)

    npos = npairs
    nneg = npairs

    coichg = nposco - nnegco
    fixchg = netchg + coichg

    if fixchg < 0:
        npos += -fixchg
    elif fixchg > 0:
        nneg += fixchg

    if npos < 0:
        npos = 0
    if nneg < 0:
        nneg = 0

print(npos + nposco, nneg + nnegco, npairs)
PY
}

read -r NPOS_TOTAL NNEG_TOTAL NPAIRS < <(calc_ion_counts)

LEAPIN="${TMPDIR}/final.leap"

{
    write_common_leap
    echo
    write_solvate_cmd "$FINAL_CUSHION"
    echo

    if [[ "$NPOS_TOTAL" -gt 0 ]]; then
        echo "addionsrand pep Na+ ${NPOS_TOTAL}"
    fi

    if [[ "$NNEG_TOTAL" -gt 0 ]]; then
        echo "addionsrand pep Cl- ${NNEG_TOTAL}"
    fi

    echo
    echo "check pep"
    echo "charge pep"
    echo
    echo "saveamberparm pep ${OUTNAME}.parm7 ${OUTNAME}.rst7"
    echo "savepdb pep ${OUTNAME}.pdb"
    echo
    echo "quit"
} > "$LEAPIN"

{
    echo "========================================="
    echo "Solvation System Variables Confirmed:"
    echo "Input System       : $INNAME"
    echo "Output System      : $OUTNAME"
    echo "Force Field        : $FF_CHOICE"
    echo "Protein Leaprc     : $PROTEIN_LEAPRC"
    echo "Water Leaprc       : $WATER_LEAPRC"
    echo "Water Box          : $WATER_BOX"
    echo "CpHMD Parm Dir     : $PHMD_DIR"
    echo "Extra FF Params    : ${EXTRA_FF_PARAMS_PATH:-none}"
    echo "Box Type           : $BOX_TYPE"
    echo "Input Cushion      : $CUSHION Angstrom"
    echo "Final Cushion      : $FINAL_CUSHION Angstrom"
    echo "Requested Edge     : ${BOX_EDGE:-none}"
    echo "Pre-ion Edge       : $FINAL_EDGE"
    echo "Best Fit Edge      : ${FIT_EDGE:-not-fit}"
    echo "Salt Conc          : $CONC_MM mM"
    echo "Water Count        : $NWAT"
    echo "Salt Pairs         : $NPAIRS"
    echo "Titration Flag     : $TITR_FLAG"
    echo "Ion Mode           : $ION_MODE"
    echo "Exact SOD          : $SOD_EXACT"
    echo "Exact CLA          : $CLA_EXACT"
    echo "Solute Charge      : $NETCHG"
    echo "pHMD + co-ions     : $NPOSCO"
    echo "pHMD - co-ions     : $NNEGCO"
    echo "Final Na+ count    : $NPOS_TOTAL"
    echo "Final Cl- count    : $NNEG_TOTAL"
    echo "Keep Temp Files    : $KEEP_TEMP"
    echo "Log File           : $LOGFILE"
    echo "========================================="
} | tee "$LOGFILE"

tleap -f "$LEAPIN" >> "$LOGFILE"

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

echo "Process complete."
mv "${OUTNAME}.parm7" "${OUTNAME}.prmtop"
mv "${OUTNAME}.rst7" "${OUTNAME}.inpcrd"
echo "Output saved to ${OUTNAME}.prmtop/.inpcrd/.pdb"
echo "Amber log saved to $LOGFILE"
