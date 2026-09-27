#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# STRUCTURE INPUT
###############################################################################

# Force field: c22, ff14sb, or ff19sb
ff="c22"

# Set exactly one:
pdbid="2lzt"          # Fetch this PDB ID from RCSB
pdbfile=""            # Existing PDB, e.g. ./1ubq_original.pdb

chain="A"             # One protein chain
cushion="15"          # Water cushion in Angstrom
boxtype="octahedral"  # C22: cubic/rhombic/octahedral; Amber: cubic/octahedral
boxedge=""            # Optional Amber-only requested final box edge in Angstrom
capped="1"            # 1 = capped; 0 = uncapped

sod="1"               # Exact Na+ count when conc is blank
cla="0"               # Exact Cl- count when conc is blank
conc="150"            # Salt concentration in mM; blank uses sod/cla

titr="0"              # -1 acidic co-ions; 0 none; 1 basic co-ions

###############################################################################

SCRIPT_DIR="$(
    cd -- "$(dirname -- "${BASH_SOURCE[0]}")"
    pwd -P
)"
PREP_DIR="$(
    cd -- "${SCRIPT_DIR}/.."
    pwd -P
)"
CONFIG="${PREP_DIR}/config.sh"

if [[ ! -f "$CONFIG" ]]; then
    echo "Error: missing configuration file: $CONFIG" >&2
    echo "Copy config.example.sh to config.sh and set the local paths." >&2
    exit 1
fi

# shellcheck disable=SC1090
source "$CONFIG"
export AMBER_ENV C22_TOPPAR FF14_PHMD_DIR FF19_PHMD_DIR

die() {
    echo "Error: $*" >&2
    exit 1
}

sanitize_num_tag() {
    local value="$1"
    value="$(printf "%s" "$value" | sed -E 's/\.0+$//')"
    value="${value//./p}"
    printf "%s" "$value"
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

make_ion_tag() {
    if [[ -n "$conc" ]]; then
        format_mM_tag "$conc"
        return
    fi

    local tag=""

    if [[ "$sod" -gt 0 ]]; then
        tag="${sod}xsod"
    fi

    if [[ "$cla" -gt 0 ]]; then
        if [[ -n "$tag" ]]; then
            tag="${tag}_${cla}xcla"
        else
            tag="${cla}xcla"
        fi
    fi

    if [[ -z "$tag" ]]; then
        tag="0ions"
    fi

    printf "%s" "$tag"
}

load_amber_env() {
    [[ -f "${AMBER_ENV:-}" ]] || die "Amber environment not found: ${AMBER_ENV:-unset}"
    set +u
    # shellcheck disable=SC1090
    source "$AMBER_ENV"
    set -u
}

ff="$(printf "%s" "$ff" | tr '[:upper:]' '[:lower:]')"
case "$ff" in
    c22|charmm22) ff="c22" ;;
    ff14sb|ff14) ff="ff14sb" ;;
    ff19sb|ff19) ff="ff19sb" ;;
    *) die "ff must be c22, ff14sb, or ff19sb" ;;
esac

if [[ -n "$pdbid" && -n "$pdbfile" ]]; then
    die "set either pdbid or pdbfile, not both"
fi

if [[ -z "$pdbid" && -z "$pdbfile" ]]; then
    die "either pdbid or pdbfile must be set"
fi

[[ "$chain" =~ ^[[:alnum:]]$ ]] || die "chain must contain exactly one alphanumeric character"
[[ "$capped" == "0" || "$capped" == "1" ]] || die "capped must be 0 or 1"
[[ "$titr" == "-1" || "$titr" == "0" || "$titr" == "1" ]] || die "titr must be -1, 0, or 1"

chain_select="$(printf "%s" "$chain" | tr '[:lower:]' '[:upper:]')"
chain_tag="$(printf "%s" "$chain" | tr '[:upper:]' '[:lower:]')"

if [[ -n "$pdbid" ]]; then
    pdbid="$(printf "%s" "$pdbid" | tr '[:upper:]' '[:lower:]')"
    [[ "$pdbid" =~ ^[[:alnum:]]{4}$ ]] || die "pdbid must contain exactly four alphanumeric characters"

    pdb_original="${pdbid}_original.pdb"
    bash "${SCRIPT_DIR}/get_pdb.sh" --pdbid "$pdbid" --out "$pdb_original"
else
    [[ -f "$pdbfile" ]] || die "PDB file not found: $pdbfile"
    [[ -s "$pdbfile" ]] || die "PDB file is empty: $pdbfile"

    pdb_original="$pdbfile"
    pdb_basename="$(basename -- "$pdbfile")"

    if [[ "$pdb_basename" =~ ^([[:alnum:]]{4})([_\.-]|$) ]]; then
        pdbid="$(printf "%s" "${BASH_REMATCH[1]}" | tr '[:upper:]' '[:lower:]')"
    else
        die "local PDB filename must begin with its four-character PDB ID: $pdb_basename"
    fi
fi

if [[ "$ff" == "c22" ]]; then
    [[ -n "${C22_TOPPAR:-}" ]] || die "C22_TOPPAR is not set in config.sh"
    [[ -z "$boxedge" ]] || die "boxedge is only supported by the Amber backend"

    case "$(printf "%s" "$boxtype" | tr '[:upper:]' '[:lower:]')" in
        cubic|cube|cub|1) boxtype="cubic" ;;
        octahedral|octa|oct|truncated_octahedral|2) boxtype="octahedral" ;;
        rhombic|rhdo|dodecahedron|rhombic_dodecahedron|3) boxtype="rhombic" ;;
        *) die "C22 boxtype must be cubic, octahedral, or rhombic" ;;
    esac

    pdb_clean="${pdbid}_${chain_tag}_clean.pdb"
    disulfides="${pdbid}_${chain_tag}_disulfides.dat"
    prep_report="${pdbid}_${chain_tag}_prep_report.txt"

    bash "${SCRIPT_DIR}/c22/prep_pdb.sh" \
        --in "$pdb_original" \
        --pdbid "$pdbid" \
        --chain "$chain_select" \
        --out "$pdb_clean" \
        --disu-out "$disulfides" \
        --report "$prep_report"

    if [[ "$capped" == "1" ]]; then
        termini="capped"
        bash "${SCRIPT_DIR}/c22/gen_pdb_system.sh" \
            --in "$pdb_clean" \
            --out "${termini}_${pdbid}_${chain_tag}" \
            --disu "$disulfides" \
            --capped
    else
        termini="uncapped"
        bash "${SCRIPT_DIR}/c22/gen_pdb_system.sh" \
            --in "$pdb_clean" \
            --out "${termini}_${pdbid}_${chain_tag}" \
            --disu "$disulfides" \
            --uncapped
    fi

    prepared="${termini}_${pdbid}_${chain_tag}"

    solvate_args=(
        --in "$prepared"
        --box "$boxtype"
        --cushion "$cushion"
        --titr "$titr"
    )

    if [[ -n "$conc" ]]; then
        solvate_args+=(--conc "$conc")
    else
        solvate_args+=(--sod "$sod" --cla "$cla")
    fi

    bash "${PREP_DIR}/c22/solvate_sys.sh" "${solvate_args[@]}"

    cushion_tag="$(sanitize_num_tag "$cushion")"
    ion_tag="$(make_ion_tag)"
    system="${prepared}_solv_${cushion_tag}A_${boxtype}_${ion_tag}"

    load_amber_env
    bash "${PREP_DIR}/c22/run_chamber.sh" "$C22_TOPPAR" "$system"
else
    case "$(printf "%s" "$boxtype" | tr '[:upper:]' '[:lower:]')" in
        cubic|cube|cub|1) boxtype="cubic" ;;
        octahedral|octa|oct|truncated_octahedral|truncated-octahedral|2) boxtype="octahedral" ;;
        *) die "Amber boxtype must be cubic or octahedral" ;;
    esac

    pdb_clean="${pdbid}_${chain_tag}_amber_clean.pdb"
    disulfides="${pdbid}_${chain_tag}_amber_disulfides.dat"
    prep_report="${pdbid}_${chain_tag}_amber_prep_report.txt"

    bash "${SCRIPT_DIR}/amber/prep_pdb_amber.sh" \
        --in "$pdb_original" \
        --pdbid "$pdbid" \
        --chain "$chain_select" \
        --out "$pdb_clean" \
        --disu-out "$disulfides" \
        --report "$prep_report"

    if [[ "$capped" == "1" ]]; then
        termini="capped"
        bash "${SCRIPT_DIR}/amber/gen_pdb_system_amber.sh" \
            --in "$pdb_clean" \
            --out "${termini}_${pdbid}_${chain_tag}" \
            --disu "$disulfides" \
            --ff "$ff" \
            --capped
    else
        termini="uncapped"
        bash "${SCRIPT_DIR}/amber/gen_pdb_system_amber.sh" \
            --in "$pdb_clean" \
            --out "${termini}_${pdbid}_${chain_tag}" \
            --disu "$disulfides" \
            --ff "$ff" \
            --uncapped
    fi

    prepared="${termini}_${pdbid}_${chain_tag}"

    solvate_args=(
        --in "$prepared"
        --ff "$ff"
        --box "$boxtype"
        --cushion "$cushion"
        --titr "$titr"
    )

    if [[ -n "$boxedge" ]]; then
        solvate_args+=(--box-edge "$boxedge")
    fi

    if [[ -n "$conc" ]]; then
        solvate_args+=(--conc "$conc")
    else
        solvate_args+=(--sod "$sod" --cla "$cla")
    fi

    bash "${PREP_DIR}/amber/solvate_sys_amber.sh" "${solvate_args[@]}"

    ion_tag="$(make_ion_tag)"
    if [[ -n "$boxedge" ]]; then
        edge_tag="$(sanitize_num_tag "$boxedge")"
        system="${prepared}_solv_${edge_tag}Aedge_${boxtype}_${ion_tag}"
    else
        cushion_tag="$(sanitize_num_tag "$cushion")"
        system="${prepared}_solv_${cushion_tag}A_${boxtype}_${ion_tag}"
    fi
fi

prmtop="${system}.prmtop"
inpcrd="${system}.inpcrd"
final_pdb="${system}.pdb"

[[ -s "$prmtop" ]] || die "final topology not found or empty: $prmtop"
[[ -s "$inpcrd" ]] || die "final coordinates not found or empty: $inpcrd"
[[ -s "$final_pdb" ]] || die "final PDB not found or empty: $final_pdb"

awk -v sysname="$system" -v ffname="$ff" '
/^%FLAG RESIDUE_LABEL/ {
    reading = 1
    next
}

/^%FLAG/ && reading {
    exit
}

reading && /^%FORMAT/ {
    next
}

reading {
    for (i = 1; i <= length($0); i += 4) {
        residue = toupper(substr($0, i, 4))
        gsub(/[[:space:]]/, "", residue)

        if (residue == "WAT") {
            waters++
        } else if (residue == "NA+" || residue == "NA" || residue == "SOD") {
            sodium++
        } else if (residue == "CL-" || residue == "CL" || residue == "CLA") {
            chloride++
        }
    }
}

END {
    printf "\n"
    printf "============================================================\n"
    printf "Final System Composition\n"
    printf "Force field    : %s\n", ffname
    printf "System         : %s\n", sysname
    printf "Waters         : %d\n", waters + 0
    printf "Sodium ions    : %d\n", sodium + 0
    printf "Chloride ions  : %d\n", chloride + 0
    printf "Total ions     : %d\n", sodium + chloride
    printf "============================================================\n"
}
' "$prmtop"

echo
echo "Final topology    : $prmtop"
echo "Final coordinates : $inpcrd"
echo "Final PDB         : $final_pdb"
echo "Preparation report: $prep_report"
