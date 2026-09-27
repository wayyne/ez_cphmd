#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# SEQUENCE INPUT
###############################################################################

# Force field: c22, ff14sb, or ff19sb
ff="c22"

pep="gdg"
cushion="15"
boxtype="octahedral"  # C22: cubic, rhombic, octahedral; Amber: cubic, octahedral
boxedge=""            # Optional Amber-only requested final box edge in Angstrom
capped="1"            # 1 = capped; 0 = uncapped/zwitterionic

sod="1"               # Exact Na+ count when conc is blank
cla="0"               # Exact Cl- count when conc is blank
conc=""                # Salt concentration in mM; blank uses sod/cla

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
    c22|charmm22)
        ff="c22"
        ;;
    ff14sb|ff14)
        ff="ff14sb"
        ;;
    ff19sb|ff19)
        ff="ff19sb"
        ;;
    *)
        die "ff must be c22, ff14sb, or ff19sb"
        ;;
esac

[[ "$capped" == "0" || "$capped" == "1" ]] || die "capped must be 0 or 1"
[[ "$titr" == "-1" || "$titr" == "0" || "$titr" == "1" ]] || die "titr must be -1, 0, or 1"

seq_up="$(printf "%s" "$pep" | tr '[:lower:]' '[:upper:]')"

if [[ "$ff" == "c22" ]]; then
    [[ -n "${C22_TOPPAR:-}" ]] || die "C22_TOPPAR is not set in config.sh"
    [[ -z "$boxedge" ]] || die "boxedge is only supported by the Amber backend"

    case "$(printf "%s" "$boxtype" | tr '[:upper:]' '[:lower:]')" in
        cubic|cube|cub|1) boxtype="cubic" ;;
        octahedral|octa|oct|truncated_octahedral|2) boxtype="octahedral" ;;
        rhombic|rhdo|dodecahedron|rhombic_dodecahedron|3) boxtype="rhombic" ;;
        *) die "C22 boxtype must be cubic, octahedral, or rhombic" ;;
    esac

    if [[ "$capped" == "1" ]]; then
        bash "${SCRIPT_DIR}/c22/genpep.sh" --pep "$pep" --capped
        prepared="capped_${seq_up}"
    else
        bash "${SCRIPT_DIR}/c22/genpep.sh" --pep "$pep" --uncapped
        prepared="zwitter_${seq_up}"
    fi

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
        octahedral|octa|oct|truncated_octahedral|2) boxtype="octahedral" ;;
        *) die "Amber boxtype must be cubic or octahedral" ;;
    esac

    if [[ "$capped" == "1" ]]; then
        bash "${SCRIPT_DIR}/amber/genpep_amber.sh" --pep "$pep" --ff "$ff" --capped
        prepared="capped_${seq_up}"
    else
        bash "${SCRIPT_DIR}/amber/genpep_amber.sh" --pep "$pep" --ff "$ff" --uncapped
        prepared="zwitter_${seq_up}"
    fi

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

echo
echo "Preparation complete."
echo "Force field       : $ff"
echo "Final topology    : $prmtop"
echo "Final coordinates : $inpcrd"
echo "Final PDB         : $final_pdb"
