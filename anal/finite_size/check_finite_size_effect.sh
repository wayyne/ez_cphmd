#!/usr/bin/env bash

set -euo pipefail

if [[ $# -ne 2 ]]; then
    echo "Usage: $0 <prmtop> <trajectory.nc>" >&2
    exit 1
fi

prmtop=$1
traj=$2

if [[ ! -f "$prmtop" ]]; then
    echo "Error: topology not found: $prmtop" >&2
    exit 1
fi

if [[ ! -f "$traj" ]]; then
    echo "Error: trajectory not found: $traj" >&2
    exit 1
fi

workdir=$(mktemp -d "${TMPDIR:-/tmp}/finite_size_check.XXXXXX")
trap 'rm -rf "$workdir"' EXIT
volume_file="$workdir/volume.dat"
watinfo="$workdir/watinfo.txt"

# Get volume.
cpptraj -p "$prmtop" <<EOF
trajin "$traj"
volume out "$volume_file"
run
EOF

awk '
    !/^#/ {
        n++
        x=$2
        sum+=x
        sumsq+=x*x
    }
    END {
        if (n == 0) {
            print "Error: no volume data were produced" > "/dev/stderr"
            exit 1
        }
        mean=sum/n
        sd=(n > 1) ? sqrt((sumsq-sum*sum/n)/(n-1)) : 0
        printf "Frames: %d\nMean volume: %.3f A^3\nVolume SD: %.3f A^3\n",
               n, mean, sd
    }
' "$volume_file"

# Get number of waters.
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

# Compute water number density.
awk -v n="$NWAT" '
    !/^#/ {
        frames++
        volume=$2
        rho=n/volume

        sum_volume+=volume
        sum_rho+=rho
        sumsq_rho+=rho*rho
    }
    END {
        if (frames == 0) {
            print "Error: no volume data were produced" > "/dev/stderr"
            exit 1
        }
        mean_volume=sum_volume/frames
        mean_rho=sum_rho/frames
        sd_rho=(frames > 1) \
            ? sqrt((sumsq_rho-sum_rho*sum_rho/frames)/(frames-1)) \
            : 0

        printf "Waters:                    %d\n", n
        printf "Frames:                    %d\n", frames
        printf "Mean volume:               %.3f A^3\n", mean_volume
        printf "Mean <Nwater/V>:           %.8f A^-3\n", mean_rho
        printf "Mean <Nwater/V>:           %.5f nm^-3\n", mean_rho*1000
        printf "Framewise density SD:      %.8f A^-3\n", sd_rho
        printf "Density from Nwater/<V>:   %.8f A^-3\n", n/mean_volume
    }
' "$volume_file"
