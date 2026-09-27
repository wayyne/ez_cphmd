#!/usr/bin/env bash
set -euo pipefail

# analyze_arex.sh
#
# Residue-centric analysis of one or more Amber/AREX trajectories using CPPTRAJ.
# The script processes each trajectory independently so pH-slot identity is
# preserved (e.g., arex_ph*.nc), then optionally calls plot_arex_analysis.py.
#
# Analyses:
#   - backbone RMSD (and optional coordinate fitting)
#   - radius of gyration
#   - per-residue SASA
#   - per-residue side-chain SASA
#   - first/second-shell water counts around a residue/site
#   - residue/site <-> water hydrogen-bond counts
#   - residue/site <-> rest-of-protein hydrogen-bond counts
#
# Default local hydrogen-bond definition:
#   donor-acceptor distance <= 3.0 A and donor-H-acceptor angle >= 135 deg.
#
# The scalar SASA/H-bond/watershell quantities do not require global alignment,
# but fitting is enabled by default so the same processed coordinates can also
# be used for structural analyses and optional aligned trajectory output.

usage() {
    cat <<'EOF'
Usage:
  analyze_arex.sh -p system.parm7 -t 'arex_ph*.nc' [options]
  analyze_arex.sh -p system.parm7 [options] traj1.nc traj2.nc ...

Required:
  -p, --topology FILE          Amber parm7/prmtop topology.
  -t, --traj PATTERN           Trajectory file or shell pattern. May be repeated.
                               Quote globs so the script expands them itself.

Residue selection:
  -r, --residues SPEC          Comma-separated biological/sequence residue
                               numbers or "all". Default: all residues selected
                               by --protein-mask.
  --residue-offset N           Add N to each reported/sequence residue number
                               to obtain the CPPTRAJ/topology residue index.
                               Example: N-terminal cap => --residue-offset 1.
                               Default: 0.
  --protein-mask MASK          CPPTRAJ mask defining the protein/solute.
                               Default: ^1 (first molecule).
  --site-mode MODE             "residue" or "sidechain".
                               Controls the mask used for water/H-bond analyses.
                               Default: sidechain.
  --site-file FILE             Optional tab-separated file:
                                   sequence_residue<TAB>cpptraj_mask_template
                               A custom mask overrides --site-mode. The template
                               may contain {resid} (sequence number) and
                               {topres} (topology/CPPTRAJ residue index).

Trajectory processing:
  --ref FILE                   Optional reference structure. If omitted, frame 1
                               of each trajectory is the fitting reference.
                               The reference should have atom ordering compatible
                               with the supplied topology and fit mask.
  --fit-mask MASK              Mask used for RMS fitting.
                               Default: (^1)&@N,CA,C
  --no-align                   Calculate RMSD but do not rotate/translate frames.
  --write-aligned              Write the processed trajectory in NetCDF format.
  --start N                    First trajectory frame. Default: 1
  --stop N|last               Last trajectory frame. Default: last
  --stride N                   Frame stride. Default: 1
  --dt-ps FLOAT                Time spacing (ps) between STORED INPUT frames.
                               Used only for plotting. Default: 1.0

Analysis:
  --metrics LIST               Comma-separated subset of:
                               sasa,watershell,hbond
                               Default: sasa,watershell,hbond
  --water-mask MASK            Solvent mask for watershell / water donors.
                               Default: :WAT
  --water-acceptor-mask MASK   Water acceptor mask for H-bond analysis.
                               Default: :WAT@O
  --first-shell FLOAT          First water-shell cutoff (A). Default: 3.5
  --second-shell FLOAT         Second water-shell cutoff (A). Default: 5.0
  --hbond-dist FLOAT           H-bond heavy-atom distance cutoff (A). Default: 3.0
  --hbond-angle FLOAT          H-bond angle cutoff (deg). Default: 135

Output:
  -o, --outdir DIR             Output directory. Default: arex_analysis
  --no-plot                    Do not run plot_arex_analysis.py after CPPTRAJ.
  -h, --help                   Show this message.

HEWL titratable-residue example:
  ./analyze_arex.sh \
      -p hewl.parm7 \
      -t 'arex_ph*.nc' \
      -r 7,15,18,35,48,52,66,87,101,119 \
      --residue-offset 0 \
      --site-file hewl_titratable_sites.tsv \
      --dt-ps 10 \
      -o hewl_arex_analysis

All-residue SASA/water-shell example:
  ./analyze_arex.sh \
      -p hewl.parm7 \
      -t 'arex_ph*.nc' \
      -r all \
      --site-mode residue \
      --metrics sasa,watershell \
      -o hewl_all_residues
EOF
}

die() {
    echo "ERROR: $*" >&2
    exit 1
}

warn() {
    echo "WARNING: $*" >&2
}

have_metric() {
    local needle="$1"
    case ",${METRICS}," in
        *",${needle},"*) return 0 ;;
        *) return 1 ;;
    esac
}

extract_ph() {
    local name="$1"
    local value=""

    # Accept both decimal-point and Amber/AREX underscore pH notation.
    # Examples:
    #   arex.ph0_5  -> 0.5
    #   arex.ph5_0  -> 5.0
    #   arex_ph7.0  -> 7.0
    #   pH_7.0      -> 7.0
    #   PH-7.0      -> 7.0
    if [[ "$name" =~ [pP][hH][_=-]?([0-9]+([._][0-9]+)?) ]]; then
        value="${BASH_REMATCH[1]}"
        value="${value/_/.}"
        printf '%s' "$value"
    else
        printf ''
    fi
}

# Defaults
TOP=""
OUTDIR="arex_analysis"
RES_SPEC="all"
RESIDUE_OFFSET=0
PROTEIN_MASK="^1"
FIT_MASK="(^1)&@N,CA,C"
SITE_MODE="sidechain"
SITE_FILE=""
REF=""
ALIGN=1
WRITE_ALIGNED=0
START=1
STOP="last"
STRIDE=1
DT_PS=1.0
METRICS="sasa,watershell,hbond"
WATER_MASK=":WAT"
WATER_ACCEPTOR_MASK=":WAT@O"
FIRST_SHELL=3.5
SECOND_SHELL=5.0
HBOND_DIST=3.0
HBOND_ANGLE=135
DO_PLOT=1

declare -a TRAJ_SPECS=()

while [[ $# -gt 0 ]]; do
    case "$1" in
        -p|--topology)
            [[ $# -ge 2 ]] || die "$1 requires an argument"
            TOP="$2"; shift 2 ;;
        -t|--traj)
            [[ $# -ge 2 ]] || die "$1 requires an argument"
            TRAJ_SPECS+=("$2"); shift 2 ;;
        -r|--residues)
            [[ $# -ge 2 ]] || die "$1 requires an argument"
            RES_SPEC="$2"; shift 2 ;;
        --residue-offset)
            [[ $# -ge 2 ]] || die "$1 requires an argument"
            RESIDUE_OFFSET="$2"; shift 2 ;;
        --protein-mask)
            PROTEIN_MASK="$2"; shift 2 ;;
        --fit-mask)
            FIT_MASK="$2"; shift 2 ;;
        --site-mode)
            SITE_MODE="$2"; shift 2 ;;
        --site-file)
            SITE_FILE="$2"; shift 2 ;;
        --ref)
            REF="$2"; shift 2 ;;
        --no-align)
            ALIGN=0; shift ;;
        --write-aligned)
            WRITE_ALIGNED=1; shift ;;
        --start)
            START="$2"; shift 2 ;;
        --stop)
            STOP="$2"; shift 2 ;;
        --stride)
            STRIDE="$2"; shift 2 ;;
        --dt-ps)
            DT_PS="$2"; shift 2 ;;
        --metrics)
            METRICS="$2"; shift 2 ;;
        --water-mask)
            WATER_MASK="$2"; shift 2 ;;
        --water-acceptor-mask)
            WATER_ACCEPTOR_MASK="$2"; shift 2 ;;
        --first-shell)
            FIRST_SHELL="$2"; shift 2 ;;
        --second-shell)
            SECOND_SHELL="$2"; shift 2 ;;
        --hbond-dist)
            HBOND_DIST="$2"; shift 2 ;;
        --hbond-angle)
            HBOND_ANGLE="$2"; shift 2 ;;
        -o|--outdir)
            OUTDIR="$2"; shift 2 ;;
        --no-plot)
            DO_PLOT=0; shift ;;
        -h|--help)
            usage; exit 0 ;;
        --)
            shift
            while [[ $# -gt 0 ]]; do
                TRAJ_SPECS+=("$1")
                shift
            done
            ;;
        -*)
            die "Unknown option: $1" ;;
        *)
            # Friendly handling of unflagged trajectory files/globs.
            TRAJ_SPECS+=("$1"); shift ;;
    esac
done

[[ -n "$TOP" ]] || { usage >&2; die "Topology is required."; }
[[ -f "$TOP" ]] || die "Topology not found: $TOP"
command -v cpptraj >/dev/null 2>&1 || die "cpptraj not found in PATH."
[[ "$SITE_MODE" == "residue" || "$SITE_MODE" == "sidechain" ]] ||
    die "--site-mode must be 'residue' or 'sidechain'."
[[ "$RESIDUE_OFFSET" =~ ^-?[0-9]+$ ]] ||
    die "--residue-offset must be an integer."
[[ "$STRIDE" =~ ^[1-9][0-9]*$ ]] || die "--stride must be a positive integer."
[[ "$START" =~ ^[1-9][0-9]*$ ]] || die "--start must be a positive integer."
[[ "$STOP" == "last" || "$STOP" =~ ^[1-9][0-9]*$ ]] ||
    die "--stop must be a positive integer or 'last'."
if [[ -n "$SITE_FILE" && ! -f "$SITE_FILE" ]]; then
    die "Site file not found: $SITE_FILE"
fi
if [[ -n "$REF" && ! -f "$REF" ]]; then
    die "Reference not found: $REF"
fi

# Expand trajectory specifications and sort naturally.
declare -a TRAJS=()
for spec in "${TRAJ_SPECS[@]}"; do
    if [[ -f "$spec" ]]; then
        TRAJS+=("$spec")
    else
        mapfile -t matches < <(compgen -G "$spec" || true)
        if [[ ${#matches[@]} -eq 0 ]]; then
            warn "No trajectory matched: $spec"
        else
            TRAJS+=("${matches[@]}")
        fi
    fi
done

[[ ${#TRAJS[@]} -gt 0 ]] || die "No trajectory files were found."
mapfile -t TRAJS < <(printf '%s\n' "${TRAJS[@]}" | awk '!seen[$0]++' | sort -V)

mkdir -p "$OUTDIR"

# Build protein residue map once using CPPTRAJ. This is also provenance for
# how "all" was interpreted.
PROTEIN_RES_MAP="$OUTDIR/protein_residues.tsv"
{
    printf "top_resid\tresname\n"
    cpptraj -p "$TOP" --resmask "$PROTEIN_MASK" 2>/dev/null |
        awk '
            $1 ~ /^[0-9]+$/ && $2 ~ /^[[:alnum:]_+-]+$/ {
                print $1 "\t" $2
            }
        '
} > "$PROTEIN_RES_MAP"

# Deduplicate numeric residue rows if CPPTRAJ printed any repeated records.
{
    head -n 1 "$PROTEIN_RES_MAP"
    tail -n +2 "$PROTEIN_RES_MAP" | awk -F '\t' '!seen[$1]++'
} > "${PROTEIN_RES_MAP}.tmp"
mv "${PROTEIN_RES_MAP}.tmp" "$PROTEIN_RES_MAP"

declare -a RESIDUES=()

top_resid_for() {
    local seq_r="$1"
    printf '%d' "$((seq_r + RESIDUE_OFFSET))"
}

if [[ "$RES_SPEC" == "all" ]]; then
    # protein_residues.tsv stores CPPTRAJ/topology residue indices. Convert
    # them back to reported/sequence numbering by subtracting the offset.
    # Non-positive sequence numbers (e.g., an N-terminal cap at topology
    # residue 1 with offset +1) are omitted.
    while IFS=$'\t' read -r top_r rn; do
        [[ "$top_r" =~ ^[1-9][0-9]*$ ]] || continue
        seq_r=$((top_r - RESIDUE_OFFSET))
        if (( seq_r > 0 )); then
            RESIDUES+=("$seq_r")
        fi
    done < <(tail -n +2 "$PROTEIN_RES_MAP")

    [[ ${#RESIDUES[@]} -gt 0 ]] ||
        die "Could not determine residues for mask '$PROTEIN_MASK'. Try an explicit -r list."
else
    IFS=',' read -r -a RESIDUES <<< "$RES_SPEC"
    for r in "${RESIDUES[@]}"; do
        [[ "$r" =~ ^[1-9][0-9]*$ ]] ||
            die "Invalid residue number in --residues: '$r'"
        top_r="$(top_resid_for "$r")"
        (( top_r > 0 )) ||
            die "Residue $r with offset $RESIDUE_OFFSET maps to invalid topology residue $top_r."
    done
fi

get_resname() {
    local seq_r="$1"
    local top_r
    local name
    top_r="$(top_resid_for "$seq_r")"

    name="$(awk -F '\t' -v r="$top_r" '$1 == r {print $2; exit}' "$PROTEIN_RES_MAP")"
    if [[ -z "$name" ]]; then
        # Fallback for explicitly selected residues outside the protein mask.
        name="$(cpptraj -p "$TOP" --resmask ":$top_r" 2>/dev/null |
            awk -v r="$top_r" '$1 == r {print $2; exit}')"
    fi
    [[ -n "$name" ]] || name="RES"
    printf '%s' "$name"
}

get_custom_site_mask() {
    local seq_r="$1"
    [[ -n "$SITE_FILE" ]] || return 0
    awk -F '\t' -v r="$seq_r" '
        $0 !~ /^[[:space:]]*#/ && $1 == r {
            print $2
            exit
        }
    ' "$SITE_FILE"
}

site_mask_for() {
    local seq_r="$1"
    local top_r
    local custom
    top_r="$(top_resid_for "$seq_r")"
    custom="$(get_custom_site_mask "$seq_r")"

    if [[ -n "$custom" ]]; then
        custom="${custom//\{resid\}/$seq_r}"
        custom="${custom//\{topres\}/$top_r}"
        printf '%s' "$custom"
    elif [[ "$SITE_MODE" == "residue" ]]; then
        printf ':%s' "$top_r"
    else
        # Full side chain, including side-chain donor hydrogens but excluding
        # standard backbone atoms/hydrogens.
        printf '(:%s)&!@N,CA,C,O,OXT,H,H1,H2,H3,HA,HA2,HA3' "$top_r"
    fi
}

# Record residue/site definitions.
RES_FILE="$OUTDIR/residues.tsv"
printf "resid\ttop_resid\tresname\tsite_mask\n" > "$RES_FILE"
for r in "${RESIDUES[@]}"; do
    top_r="$(top_resid_for "$r")"
    rn="$(get_resname "$r")"
    sm="$(site_mask_for "$r")"
    printf "%s\t%s\t%s\t%s\n" "$r" "$top_r" "$rn" "$sm" >> "$RES_FILE"
done

METADATA="$OUTDIR/metadata.tsv"
printf "label\tph\ttrajectory\tstart\tstop\tstride\tdt_ps\tframe_dt_ps\n" > "$METADATA"

# Store exact command-line/provenance.
{
    echo "date_utc=$(date -u +%FT%TZ)"
    echo "cpptraj=$(command -v cpptraj)"
    cpptraj -V 2>/dev/null | head -n 1 | sed 's/^/cpptraj_version=/' || true
    echo "topology=$TOP"
    echo "residue_offset=$RESIDUE_OFFSET"
    echo "protein_mask=$PROTEIN_MASK"
    echo "fit_mask=$FIT_MASK"
    echo "site_mode=$SITE_MODE"
    echo "site_file=$SITE_FILE"
    echo "reference=$REF"
    echo "align=$ALIGN"
    echo "start=$START"
    echo "stop=$STOP"
    echo "stride=$STRIDE"
    echo "dt_ps=$DT_PS"
    echo "metrics=$METRICS"
    echo "water_mask=$WATER_MASK"
    echo "water_acceptor_mask=$WATER_ACCEPTOR_MASK"
    echo "first_shell_A=$FIRST_SHELL"
    echo "second_shell_A=$SECOND_SHELL"
    echo "hbond_dist_A=$HBOND_DIST"
    echo "hbond_angle_deg=$HBOND_ANGLE"
} > "$OUTDIR/settings.txt"

FRAME_DT_PS="$(awk -v dt="$DT_PS" -v s="$STRIDE" 'BEGIN{printf "%.10g", dt*s}')"

for traj in "${TRAJS[@]}"; do
    [[ -f "$traj" ]] || die "Trajectory disappeared: $traj"

    base="$(basename "$traj")"
    label="${base%.*}"
    ph="$(extract_ph "$label")"
    run_dir="$OUTDIR/$label"
    mkdir -p "$run_dir"

    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "$label" "$ph" "$traj" "$START" "$STOP" "$STRIDE" "$DT_PS" "$FRAME_DT_PS" \
        >> "$METADATA"

    cpin="$run_dir/analysis.in"

    {
        echo "# Auto-generated by analyze_arex.sh"
        echo "parm \"$TOP\""
        echo "trajin \"$traj\" $START $STOP $STRIDE"

        if [[ -n "$REF" ]]; then
            echo "reference \"$REF\" [REF]"
        fi

        # Reconstruct molecules across periodic boundaries before geometric
        # analyses. Default autoimage anchor is the first molecule, matching
        # the default protein mask.
        echo "autoimage anchor \"$PROTEIN_MASK\""

        if [[ -n "$REF" ]]; then
            if [[ "$ALIGN" -eq 1 ]]; then
                echo "rmsd Fit \"$FIT_MASK\" ref [REF] out \"$run_dir/rmsd.dat\" mass"
            else
                echo "rmsd Fit \"$FIT_MASK\" ref [REF] out \"$run_dir/rmsd.dat\" mass nomod"
            fi
        else
            if [[ "$ALIGN" -eq 1 ]]; then
                echo "rmsd Fit \"$FIT_MASK\" first out \"$run_dir/rmsd.dat\" mass"
            else
                echo "rmsd Fit \"$FIT_MASK\" first out \"$run_dir/rmsd.dat\" mass nomod"
            fi
        fi

        echo "radgyr Rg \"$PROTEIN_MASK\" out \"$run_dir/rg.dat\" mass nomax"

        # Collect dataset names so CPPTRAJ 'create' writes compact numeric
        # tables with one column per selected residue.
        declare -a SASA_RES_DS=()
        declare -a SASA_SC_DS=()
        declare -a HBW_DS=()
        declare -a HBPD_DS=()
        declare -a HBPA_DS=()

        for r in "${RESIDUES[@]}"; do
            top_r="$(top_resid_for "$r")"
            rn="$(get_resname "$r")"
            site_mask="$(site_mask_for "$r")"
            rest_mask="(${PROTEIN_MASK})&!:${top_r}"

            if have_metric "sasa"; then
                ds="SASA_RES_R${r}"
                echo "surf $ds \":${top_r}\""
                SASA_RES_DS+=("$ds")

                # Gly has no conventional side chain. Skip side-chain SASA
                # rather than silently substituting whole-residue SASA.
                if [[ "$rn" != "GLY" ]]; then
                    sc_mask="(:${top_r})&!@N,CA,C,O,OXT,H,H1,H2,H3,HA,HA2,HA3"
                    ds="SASA_SC_R${r}"
                    echo "surf $ds \"$sc_mask\""
                    SASA_SC_DS+=("$ds")
                fi
            fi

            if have_metric "watershell"; then
                echo "watershell \"$site_mask\" out \"$run_dir/watershell_R${r}.dat\" lower $FIRST_SHELL upper $SECOND_SHELL \"$WATER_MASK\""
            fi

            if have_metric "hbond"; then
                # Solute/site <-> water H-bonds. CPPTRAJ finds donor hydrogens
                # attached to donor heavy atoms selected by the site mask.
                echo "hbond HBW_R${r} \"$site_mask\" dist $HBOND_DIST angle $HBOND_ANGLE solventdonor \"$WATER_MASK\" solventacceptor \"$WATER_ACCEPTOR_MASK\""
                HBW_DS+=("HBW_R${r}[UV]")

                # Protein H-bonds involving the site. Donor and acceptor
                # directions are split so the site is counted against the rest
                # of the protein without including within-residue H-bonds.
                echo "hbond HBPD_R${r} \"$PROTEIN_MASK\" dist $HBOND_DIST angle $HBOND_ANGLE donormask \"$site_mask\" acceptormask \"$rest_mask\""
                echo "hbond HBPA_R${r} \"$PROTEIN_MASK\" dist $HBOND_DIST angle $HBOND_ANGLE donormask \"$rest_mask\" acceptormask \"$site_mask\""
                HBPD_DS+=("HBPD_R${r}[UU]")
                HBPA_DS+=("HBPA_R${r}[UU]")
            fi
        done

        if have_metric "sasa"; then
            if [[ ${#SASA_RES_DS[@]} -gt 0 ]]; then
                printf 'create "%s" ' "$run_dir/sasa_residue.dat"
                printf '%s ' "${SASA_RES_DS[@]}"
                echo
            fi
            if [[ ${#SASA_SC_DS[@]} -gt 0 ]]; then
                printf 'create "%s" ' "$run_dir/sasa_sidechain.dat"
                printf '%s ' "${SASA_SC_DS[@]}"
                echo
            fi
        fi

        if have_metric "hbond"; then
            if [[ ${#HBW_DS[@]} -gt 0 ]]; then
                printf 'create "%s" ' "$run_dir/hbond_water.dat"
                printf '%s ' "${HBW_DS[@]}"
                echo
            fi
            if [[ ${#HBPD_DS[@]} -gt 0 ]]; then
                printf 'create "%s" ' "$run_dir/hbond_protein_donor.dat"
                printf '%s ' "${HBPD_DS[@]}"
                echo
            fi
            if [[ ${#HBPA_DS[@]} -gt 0 ]]; then
                printf 'create "%s" ' "$run_dir/hbond_protein_acceptor.dat"
                printf '%s ' "${HBPA_DS[@]}"
                echo
            fi
        fi

        if [[ "$WRITE_ALIGNED" -eq 1 ]]; then
            echo "trajout \"$run_dir/processed.nc\" netcdf"
        fi

        echo "run"
    } > "$cpin"

    echo "Running CPPTRAJ: $traj -> $run_dir"
    cpptraj -i "$cpin" > "$run_dir/cpptraj.log" 2>&1 || {
        tail -n 80 "$run_dir/cpptraj.log" >&2 || true
        die "CPPTRAJ failed for $traj. See $run_dir/cpptraj.log"
    }
done

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PLOT_SCRIPT="$SCRIPT_DIR/plot_arex_analysis.py"

if [[ "$DO_PLOT" -eq 1 ]]; then
    if [[ -f "$PLOT_SCRIPT" ]]; then
        echo "Building summary tables and figures..."
        python3 "$PLOT_SCRIPT" --analysis-dir "$OUTDIR"
    else
        warn "plot_arex_analysis.py not found next to this script; skipping plots."
    fi
fi

echo
echo "Done."
echo "Analysis directory: $OUTDIR"
echo "Residue definitions: $RES_FILE"
echo "Run metadata: $METADATA"