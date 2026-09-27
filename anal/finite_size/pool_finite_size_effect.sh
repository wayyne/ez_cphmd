#!/usr/bin/env bash

set -euo pipefail
shopt -s nullglob

if [[ $# -lt 2 || $# -gt 3 ]]; then
    echo "Usage: $0 <prmtop> '<trajectory_glob>' [skip_frames_per_file]" >&2
    echo "Example: $0 ../../pre/system.prmtop 'arex.ph*.nc' 100" >&2
    exit 1
fi

prmtop=$1
pattern=$2
skip_frames=${3:-0}

if [[ ! -f "$prmtop" ]]; then
    echo "Error: topology not found: $prmtop" >&2
    exit 1
fi

if ! [[ "$skip_frames" =~ ^[0-9]+$ ]]; then
    echo "Error: skip_frames_per_file must be a nonnegative integer." >&2
    exit 1
fi

# Expand the user-supplied pattern inside the script.
# The pattern should be quoted on the command line.
trajectories=( $pattern )

if [[ ${#trajectories[@]} -eq 0 ]]; then
    echo "Error: no trajectories matched: $pattern" >&2
    exit 1
fi

# Sort naturally when filenames contain numbers.
mapfile -t trajectories < <(printf '%s\n' "${trajectories[@]}" | sort -V)

workdir=$(mktemp -d "${TMPDIR:-/tmp}/finite_size_pool.XXXXXX")
trap 'rm -rf "$workdir"' EXIT

watinfo="$workdir/watinfo.txt"
pooled_volumes="$workdir/pooled_volumes.dat"
per_file_summary="finite_size_per_file.tsv"

# Count waters once from the shared topology.
cpptraj -p "$prmtop" <<'EOF' > "$watinfo"
resinfo :WAT
EOF

NWAT=$(awk '
    $1 ~ /^[0-9]+$/ && $2 == "WAT" {
        n++
    }
    END {
        print n+0
    }
' "$watinfo")

if [[ "$NWAT" -eq 0 ]]; then
    echo "Error: no WAT residues were detected in $prmtop" >&2
    exit 1
fi

: > "$pooled_volumes"
printf "trajectory\tframes_used\tmean_volume_A3\tvolume_sd_A3\tmean_density_A-3\n" \
    > "$per_file_summary"

echo "Topology:                 $prmtop"
echo "Trajectory pattern:       $pattern"
echo "Matched trajectories:     ${#trajectories[@]}"
echo "Frames skipped per file:  $skip_frames"
echo "Waters:                   $NWAT"
echo

for i in "${!trajectories[@]}"; do
    traj=${trajectories[$i]}

    if [[ ! -f "$traj" ]]; then
        echo "Error: trajectory not found after expansion: $traj" >&2
        exit 1
    fi

    volume_file="$workdir/volume_${i}.dat"

    cpptraj -p "$prmtop" > /dev/null <<EOF
trajin "$traj"
volume out "$volume_file"
run
EOF

    # Extract retained volumes and append them to the pooled set.
    awk -v skip="$skip_frames" '
        !/^#/ {
            seen++
            if (seen <= skip) {
                next
            }
            print $2
        }
    ' "$volume_file" >> "$pooled_volumes"

    # Report a per-trajectory summary for quality control.
    awk -v skip="$skip_frames" -v nwat="$NWAT" -v name="$traj" '
        !/^#/ {
            seen++
            if (seen <= skip) {
                next
            }

            n++
            volume=$2
            sum_volume+=volume
            sumsq_volume+=volume*volume
            sum_density+=nwat/volume
        }
        END {
            if (n == 0) {
                printf "Error: no frames remain for %s after skipping %d frames\n",
                       name, skip > "/dev/stderr"
                exit 1
            }

            mean_volume=sum_volume/n
            volume_sd=(n > 1) \
                ? sqrt((sumsq_volume-sum_volume*sum_volume/n)/(n-1)) \
                : 0
            mean_density=sum_density/n

            printf "%s\t%d\t%.6f\t%.6f\t%.10f\n",
                   name, n, mean_volume, volume_sd, mean_density
        }
    ' "$volume_file" >> "$per_file_summary"
done

echo "Per-trajectory results:"
column -t -s $'\t' "$per_file_summary" 2>/dev/null || cat "$per_file_summary"
echo

awk -v nwat="$NWAT" -v nfiles="${#trajectories[@]}" '
    {
        frames++
        volume=$1
        rho=nwat/volume

        sum_volume+=volume
        sumsq_volume+=volume*volume
        sum_rho+=rho
        sumsq_rho+=rho*rho
    }
    END {
        if (frames == 0) {
            print "Error: pooled volume set is empty" > "/dev/stderr"
            exit 1
        }

        mean_volume=sum_volume/frames
        volume_sd=(frames > 1) \
            ? sqrt((sumsq_volume-sum_volume*sum_volume/frames)/(frames-1)) \
            : 0
        mean_rho=sum_rho/frames
        rho_sd=(frames > 1) \
            ? sqrt((sumsq_rho-sum_rho*sum_rho/frames)/(frames-1)) \
            : 0

        printf "Pooled result\n"
        printf "-------------\n"
        printf "Trajectories:               %d\n", nfiles
        printf "Waters:                     %d\n", nwat
        printf "Retained frames:            %d\n", frames
        printf "Mean box volume:            %.3f A^3\n", mean_volume
        printf "Volume SD:                  %.3f A^3\n", volume_sd
        printf "Mean <Nwater/V>:            %.8f A^-3\n", mean_rho
        printf "Mean <Nwater/V>:            %.5f nm^-3\n", mean_rho*1000
        printf "Framewise density SD:       %.8f A^-3\n", rho_sd
        printf "Density from Nwater/<V>:    %.8f A^-3\n", nwat/mean_volume
    }
' "$pooled_volumes"

echo
echo "Per-trajectory table saved to: $per_file_summary"
