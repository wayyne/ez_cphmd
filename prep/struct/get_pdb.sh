#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# DOWNLOAD A LEGACY PDB-FORMAT STRUCTURE FROM RCSB
###############################################################################

PDBID=""
OUTFILE=""
REFRESH="0"

usage() {
    cat <<EOF
Usage:
  $0 --pdbid <PDBID> --out <output.pdb> [--refresh]

Required arguments:
  --pdbid <PDBID>    Four-character PDB identifier
  --out <file>       Output filename

Optional arguments:
  --refresh          Redownload and replace an existing output file
  -h, --help         Show this help message
EOF
}

die() {
    echo "Error: $*" >&2
    exit 1
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --pdbid)
            [[ $# -ge 2 ]] || die "--pdbid requires an argument"
            PDBID="$2"
            shift 2
            ;;
        --out)
            [[ $# -ge 2 ]] || die "--out requires an argument"
            OUTFILE="$2"
            shift 2
            ;;
        --refresh)
            REFRESH="1"
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

[[ -n "$PDBID" ]] || die "--pdbid is required"
[[ -n "$OUTFILE" ]] || die "--out is required"

PDBID="$(printf "%s" "$PDBID" | tr '[:upper:]' '[:lower:]')"

if ! [[ "$PDBID" =~ ^[[:alnum:]]{4}$ ]]; then
    die "PDB ID must contain exactly four alphanumeric characters"
fi

OUTDIR="$(dirname -- "$OUTFILE")"
[[ -d "$OUTDIR" ]] || die "output directory does not exist: $OUTDIR"

validate_pdb_file() {
    local pdb_file="$1"

    [[ -s "$pdb_file" ]] || return 1

    if grep -Eqi \
        '^[[:space:]]*<(html|!doctype html)|404[[:space:]]+not[[:space:]]+found' \
        "$pdb_file"; then
        return 1
    fi

    grep -q '^ATOM  ' "$pdb_file"
}

if [[ -e "$OUTFILE" && "$REFRESH" == "0" ]]; then
    if validate_pdb_file "$OUTFILE"; then
        echo "Using existing PDB file: $OUTFILE"
        exit 0
    fi

    die "existing output is not a valid coordinate-bearing PDB file: $OUTFILE"
fi

DOWNLOAD_URL="https://files.rcsb.org/download/${PDBID}.pdb"
TMPFILE="$(mktemp "${OUTDIR}/.${PDBID}.download.XXXXXX")"

cleanup() {
    rm -f -- "$TMPFILE"
}
trap cleanup EXIT

echo "Downloading PDB entry ${PDBID}..."
echo "Source: $DOWNLOAD_URL"

if command -v curl >/dev/null 2>&1; then
    curl \
        --fail \
        --location \
        --silent \
        --show-error \
        --retry 3 \
        --retry-delay 2 \
        --connect-timeout 30 \
        --output "$TMPFILE" \
        "$DOWNLOAD_URL"
elif command -v wget >/dev/null 2>&1; then
    wget \
        --quiet \
        --tries=3 \
        --timeout=30 \
        --output-document="$TMPFILE" \
        "$DOWNLOAD_URL"
else
    die "neither curl nor wget is available"
fi

if ! validate_pdb_file "$TMPFILE"; then
    die "downloaded file is not a valid coordinate-bearing legacy PDB file for ${PDBID}"
fi

mv -f -- "$TMPFILE" "$OUTFILE"
trap - EXIT

echo "Downloaded PDB entry : $PDBID"
echo "Saved original file  : $OUTFILE"
