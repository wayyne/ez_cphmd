#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# EXTRACT, NORMALIZE, AND VALIDATE ONE PROTEIN CHAIN FOR AMBER CpHMD
#
# Output residue mapping:
#   ASP -> AS2
#   GLU -> GL2
#   HIS/HID/HIE/HIP/HSD/HSE/HSP -> HIP
#   CYS participating in retained SSBOND records -> CYX
#   Other standard amino acids are unchanged.
###############################################################################

INFILE=""
PDBID=""
CHAIN=""
OUTFILE=""
DISU_OUT=""
REPORT=""

usage() {
    cat <<EOF
Usage:
  $0 \
    --in <original.pdb> \
    --pdbid <PDBID> \
    --chain <CHAIN> \
    --out <clean.pdb> \
    --disu-out <disulfides.dat> \
    --report <prep_report.txt>

Required arguments:
  --in <file>          Original PDB file
  --pdbid <PDBID>      Four-character PDB identifier
  --chain <CHAIN>      Single chain identifier to retain
  --out <file>         Cleaned Amber-compatible PDB
  --disu-out <file>    Cleaned disulfide residue pairs
  --report <file>      Structure-preparation report

Supported scope:
  - One selected protein chain
  - Zero or one MODEL record
  - Standard amino-acid residues
  - No insertion codes
  - Consecutive original residue numbering
  - Complete N, CA, C, and O backbone atoms
  - Blank alternate locations preferred over alternate location A
  - Existing hydrogens and OXT atoms removed
  - Same-chain SSBOND records retained and renumbered
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

while [[ $# -gt 0 ]]; do
    case "$1" in
        --in)
            require_argument "$1" "${2:-}"
            INFILE="$2"
            shift 2
            ;;
        --pdbid)
            require_argument "$1" "${2:-}"
            PDBID="$2"
            shift 2
            ;;
        --chain)
            require_argument "$1" "${2:-}"
            CHAIN="$2"
            shift 2
            ;;
        --out)
            require_argument "$1" "${2:-}"
            OUTFILE="$2"
            shift 2
            ;;
        --disu-out)
            require_argument "$1" "${2:-}"
            DISU_OUT="$2"
            shift 2
            ;;
        --report)
            require_argument "$1" "${2:-}"
            REPORT="$2"
            shift 2
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
[[ -n "$PDBID" ]] || die "--pdbid is required"
[[ -n "$CHAIN" ]] || die "--chain is required"
[[ -n "$OUTFILE" ]] || die "--out is required"
[[ -n "$DISU_OUT" ]] || die "--disu-out is required"
[[ -n "$REPORT" ]] || die "--report is required"

[[ -f "$INFILE" ]] || die "input PDB file not found: $INFILE"
[[ -s "$INFILE" ]] || die "input PDB file is empty: $INFILE"

PDBID="$(printf "%s" "$PDBID" | tr '[:upper:]' '[:lower:]')"
CHAIN="$(printf "%s" "$CHAIN" | tr '[:lower:]' '[:upper:]')"

[[ "$PDBID" =~ ^[[:alnum:]]{4}$ ]] ||
    die "PDB ID must contain exactly four alphanumeric characters"

[[ "$CHAIN" =~ ^[[:alnum:]]$ ]] ||
    die "chain must contain exactly one alphanumeric character"

for output_path in "$OUTFILE" "$DISU_OUT" "$REPORT"; do
    output_dir="$(dirname -- "$output_path")"
    [[ -d "$output_dir" ]] ||
        die "output directory does not exist: $output_dir"
done

[[ "$INFILE" != "$OUTFILE" ]] ||
    die "input and cleaned output filenames must be different"

command -v python3 >/dev/null 2>&1 ||
    die "python3 is required"

TMPDIR="$(mktemp -d "${TMPDIR:-/tmp}/prep_pdb_amber.XXXXXX")"

cleanup() {
    rm -rf -- "$TMPDIR"
}
trap cleanup EXIT

TMP_PDB="${TMPDIR}/clean.pdb"
TMP_DISU="${TMPDIR}/disulfides.dat"
TMP_REPORT="${TMPDIR}/prep_report.txt"
TMP_STATUS="${TMPDIR}/status"

python3 - \
    "$INFILE" \
    "$PDBID" \
    "$CHAIN" \
    "$TMP_PDB" \
    "$TMP_DISU" \
    "$TMP_REPORT" \
    "$TMP_STATUS" <<'PY'
from __future__ import annotations

import math
import sys
from collections import OrderedDict
from pathlib import Path

(
    infile_s,
    pdbid,
    target_chain,
    outfile_s,
    disu_out_s,
    report_s,
    status_s,
) = sys.argv[1:]

infile = Path(infile_s)
outfile = Path(outfile_s)
disu_out = Path(disu_out_s)
report = Path(report_s)
status = Path(status_s)

lines = infile.read_text(encoding="utf-8", errors="replace").splitlines()

three_to_one = {
    "ALA": "A",
    "ARG": "R",
    "ASN": "N",
    "ASP": "D",
    "CYS": "C",
    "GLN": "Q",
    "GLU": "E",
    "GLY": "G",
    "HIS": "H",
    "HID": "H",
    "HIE": "H",
    "HIP": "H",
    "HSD": "H",
    "HSE": "H",
    "HSP": "H",
    "ILE": "I",
    "LEU": "L",
    "LYS": "K",
    "MET": "M",
    "PHE": "F",
    "PRO": "P",
    "SER": "S",
    "THR": "T",
    "TRP": "W",
    "TYR": "Y",
    "VAL": "V",
}

histidine_names = {"HIS", "HID", "HIE", "HIP", "HSD", "HSE", "HSP"}

errors: list[str] = []
error_seen: set[str] = set()


def add_error(message: str) -> None:
    if message not in error_seen:
        error_seen.add(message)
        errors.append(message)


def field(line: str, start: int, end: int) -> str:
    return line[start:end] if len(line) >= start else ""


def parse_int(text: str) -> int | None:
    text = text.strip()
    try:
        return int(text)
    except ValueError:
        return None


def parse_float(text: str) -> float | None:
    text = text.strip()
    try:
        value = float(text)
    except ValueError:
        return None
    return value if math.isfinite(value) else None


def is_hydrogen(line: str, atom_name: str) -> bool:
    element = field(line, 76, 78).strip().upper()
    if element in {"H", "D"}:
        return True

    compact = "".join(ch for ch in atom_name.upper() if not ch.isdigit())
    return compact.startswith("H") or compact.startswith("D")


model_count = sum(1 for line in lines if line.startswith("MODEL "))
if model_count > 1:
    add_error(f"Multiple structural models are not supported; found {model_count}")

ssbond_records: list[tuple[str, int | None, str, str, int | None, str]] = []

for line in lines:
    if not line.startswith("SSBOND"):
        continue

    chain1 = field(line, 15, 16).strip()
    seq1 = parse_int(field(line, 17, 21))
    icode1 = field(line, 21, 22) or " "

    chain2 = field(line, 29, 30).strip()
    seq2 = parse_int(field(line, 31, 35))
    icode2 = field(line, 35, 36) or " "

    ssbond_records.append((chain1, seq1, icode1, chain2, seq2, icode2))

selected_atom_records = 0
selected_hetatm_records = 0
removed_hydrogens = 0
removed_oxt = 0
selected_blank_altloc = 0
selected_altloc_a = 0
normalized_histidines = 0

# Ordered by first appearance in the selected chain.
residues: OrderedDict[tuple[int, str], dict] = OrderedDict()
completed_residues: set[tuple[int, str]] = set()
last_residue_key: tuple[int, str] | None = None

for line in lines:
    record = field(line, 0, 6)

    if record == "HETATM":
        if field(line, 21, 22) == target_chain:
            selected_hetatm_records += 1
        continue

    if record != "ATOM  ":
        continue

    if field(line, 21, 22) != target_chain:
        continue

    selected_atom_records += 1

    raw_resseq = field(line, 22, 26).strip()
    icode = field(line, 26, 27) or " "
    resseq = parse_int(raw_resseq)
    resname = field(line, 17, 20).strip().upper()
    atom_name = field(line, 12, 16).strip().upper()
    altloc = field(line, 16, 17) or " "

    if resseq is None:
        add_error(
            f"Selected chain contains a non-integer residue number: {raw_resseq}"
        )
        continue

    key = (resseq, icode)

    if icode != " ":
        add_error(f"Insertion code found at original residue {resseq}{icode}")

    if resname not in three_to_one:
        add_error(f"Unsupported residue {resname} at original residue {resseq}{icode}")
        continue

    if key not in residues:
        residues[key] = {
            "original_resname": resname,
            "atoms": OrderedDict(),
        }

        if resname in histidine_names:
            normalized_histidines += 1
    elif residues[key]["original_resname"] != resname:
        add_error(
            f"Residue {resseq} has conflicting residue names: "
            f"{residues[key]['original_resname']} and {resname}"
        )

    if (
        last_residue_key is not None
        and last_residue_key != key
        and key in completed_residues
    ):
        add_error(f"Residue records are noncontiguous for original residue {resseq}")

    if last_residue_key is not None and last_residue_key != key:
        completed_residues.add(last_residue_key)

    last_residue_key = key

    if is_hydrogen(line, atom_name):
        removed_hydrogens += 1
        continue

    if atom_name == "OXT":
        removed_oxt += 1
        continue

    atom_table: OrderedDict[str, dict] = residues[key]["atoms"]

    if atom_name not in atom_table:
        atom_table[atom_name] = {
            "chosen_priority": 0,
            "line": None,
            "altloc": None,
        }

    atom = atom_table[atom_name]

    if altloc == " ":
        priority = 2
    elif altloc == "A":
        priority = 1
    else:
        # Ignore B/C/etc. if a blank or A coordinate exists. Validation below
        # fails if no supported coordinate was selected.
        continue

    if atom["chosen_priority"] == priority:
        label = "blank" if altloc == " " else altloc
        add_error(
            f"Duplicate coordinate records for residue {resseq}, atom "
            f"{atom_name}, alternate location {label}"
        )
        continue

    if priority > atom["chosen_priority"]:
        atom["chosen_priority"] = priority
        atom["line"] = line
        atom["altloc"] = altloc

if selected_atom_records == 0:
    add_error(f"No ATOM records were found for requested chain {target_chain}")

if not residues:
    add_error(f"No supported amino-acid residues were found for chain {target_chain}")

residue_keys = list(residues)
clean_number: dict[tuple[int, str], int] = {
    key: index for index, key in enumerate(residue_keys, start=1)
}

for index, key in enumerate(residue_keys):
    resseq, _ = key

    if index > 0:
        previous_resseq = residue_keys[index - 1][0]
        if resseq != previous_resseq + 1:
            add_error(
                "Nonconsecutive residue numbering between original residues "
                f"{previous_resseq} and {resseq}"
            )

    atoms = residues[key]["atoms"]

    for atom_name, atom in atoms.items():
        if atom["line"] is None:
            add_error(
                f"No blank or A alternate-location coordinate exists for "
                f"original residue {resseq}, atom {atom_name}"
            )
            continue

        line = atom["line"]
        x = parse_float(field(line, 30, 38))
        y = parse_float(field(line, 38, 46))
        z = parse_float(field(line, 46, 54))

        if x is None or y is None or z is None:
            add_error(
                f"Invalid coordinate found for original residue {resseq}, "
                f"atom {atom_name}"
            )

        if atom["altloc"] == "A":
            selected_altloc_a += 1
        else:
            selected_blank_altloc += 1

    for required_atom in ("N", "CA", "C", "O"):
        if required_atom not in atoms or atoms[required_atom]["line"] is None:
            add_error(
                f"Missing backbone atom {required_atom} at original residue "
                f"{resseq} {residues[key]['original_resname']}"
            )

disulfides: list[tuple[int, int, int, int]] = []
disulfide_seen: set[tuple[int, int]] = set()
disulfide_residues: set[int] = set()

for chain1, seq1, icode1, chain2, seq2, icode2 in ssbond_records:
    involves_selected = chain1 == target_chain or chain2 == target_chain
    if not involves_selected:
        continue

    if chain1 != target_chain or chain2 != target_chain:
        add_error(
            f"SSBOND involving selected chain {target_chain} connects to another chain"
        )
        continue

    if icode1 != " " or icode2 != " ":
        add_error("SSBOND uses an insertion-coded cysteine, which is unsupported")
        continue

    if seq1 is None or seq2 is None:
        add_error("SSBOND contains a non-integer residue number")
        continue

    key1 = (seq1, " ")
    key2 = (seq2, " ")

    if key1 not in clean_number or key2 not in clean_number:
        add_error(
            f"SSBOND references a residue missing from selected chain "
            f"{target_chain}: {seq1} or {seq2}"
        )
        continue

    if residues[key1]["original_resname"] != "CYS" or residues[key2][
        "original_resname"
    ] != "CYS":
        add_error(f"SSBOND does not connect two CYS residues: {seq1} and {seq2}")
        continue

    clean1 = clean_number[key1]
    clean2 = clean_number[key2]
    first, second = sorted((clean1, clean2))
    pair = (first, second)

    if pair in disulfide_seen:
        continue

    if first in disulfide_residues or second in disulfide_residues:
        add_error(
            f"A cysteine occurs in more than one disulfide pair: {first}-{second}"
        )
        continue

    disulfide_seen.add(pair)
    disulfide_residues.update(pair)
    disulfides.append((first, second, seq1, seq2))


def output_resname(clean_resid: int, original_resname: str) -> str:
    if original_resname == "ASP":
        return "AS2"
    if original_resname == "GLU":
        return "GL2"
    if original_resname in histidine_names:
        return "HIP"
    if original_resname == "CYS" and clean_resid in disulfide_residues:
        return "CYX"
    return original_resname


sequence = "".join(
    three_to_one[residues[key]["original_resname"]] for key in residue_keys
)

report_lines = [
    "============================================================",
    "Amber ff19SB PDB structure preparation report",
    "============================================================",
    f"PDB ID                  : {pdbid}",
    f"Input file              : {infile}",
    f"Selected chain          : {target_chain}",
    f"MODEL records           : {model_count}",
    f"Selected ATOM records   : {selected_atom_records}",
    f"Selected HETATM records : {selected_hetatm_records}",
    f"Clean residue count     : {len(residue_keys)}",
    f"Removed hydrogens       : {removed_hydrogens}",
    f"Removed OXT atoms       : {removed_oxt}",
    f"Normalized histidines   : {normalized_histidines}",
    f"Blank-altloc atoms      : {selected_blank_altloc}",
    f"Altloc-A atoms          : {selected_altloc_a}",
    f"Disulfide pairs         : {len(disulfides)}",
    f"Sequence                : {sequence}",
    "",
    "Residue mapping:",
    "  Clean  Original  Input  Amber",
]

for key in residue_keys:
    clean_resid = clean_number[key]
    original_resid = key[0]
    original_name = residues[key]["original_resname"]
    amber_name = output_resname(clean_resid, original_name)

    report_lines.append(
        f"  {clean_resid:<6d} {original_resid:<9d} "
        f"{original_name:<6s} {amber_name:<6s}"
    )

report_lines.extend(["", "Disulfide mapping:"])

if not disulfides:
    report_lines.append("  None")
else:
    for clean1, clean2, orig1, orig2 in disulfides:
        report_lines.append(
            f"  Original {orig1}-{orig2} -> Clean {clean1}-{clean2}"
        )

if errors:
    report_lines.extend(["", "Status: FAILED", "", "Errors:"])
    report_lines.extend(f"  - {message}" for message in errors)
    report.write_text("\n".join(report_lines) + "\n", encoding="utf-8")
    status.write_text("FAILED\n", encoding="utf-8")
    sys.exit(0)

report_lines.extend(["", "Status: PASSED"])
report.write_text("\n".join(report_lines) + "\n", encoding="utf-8")

disu_lines = [
    "# Amber CYX disulfide residue pairs",
    "# Columns: clean_residue_1 clean_residue_2",
]
disu_lines.extend(f"{first} {second}" for first, second, _, _ in disulfides)
disu_out.write_text("\n".join(disu_lines) + "\n", encoding="utf-8")

output_lines = [
    "REMARK   Generated by prep_pdb_amber.sh",
    f"REMARK   Source PDB ID {pdbid.upper()}",
    f"REMARK   Original chain {target_chain}",
    "REMARK   Residues renumbered consecutively from 1",
    "REMARK   ASP->AS2, GLU->GL2, histidines->HIP",
    "REMARK   Disulfide cysteines renamed CYX",
]

serial = 0

for key in residue_keys:
    clean_resid = clean_number[key]
    original_name = residues[key]["original_resname"]
    amber_name = output_resname(clean_resid, original_name)

    for atom_name, atom in residues[key]["atoms"].items():
        line = atom["line"]
        if line is None:
            continue

        x = float(field(line, 30, 38).strip())
        y = float(field(line, 38, 46).strip())
        z = float(field(line, 46, 54).strip())

        occupancy = parse_float(field(line, 54, 60))
        bfactor = parse_float(field(line, 60, 66))

        if occupancy is None:
            occupancy = 1.0
        if bfactor is None:
            bfactor = 0.0

        element = field(line, 76, 78).strip().upper()
        if not element:
            compact = "".join(ch for ch in atom_name if not ch.isdigit())
            element = compact[:1].upper()

        atom_field = field(line, 12, 16)
        if len(atom_field) != 4:
            atom_field = f"{atom_name:>4s}"

        serial += 1

        output_lines.append(
            f"ATOM  {serial:5d} {atom_field} {amber_name:>3s} "
            f"{target_chain}{clean_resid:4d}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}"
            f"{occupancy:6.2f}{bfactor:6.2f}"
            f"          {element:>2s}"
        )

last_key = residue_keys[-1]
last_resname = output_resname(
    clean_number[last_key], residues[last_key]["original_resname"]
)
output_lines.append(
    f"TER   {serial + 1:5d}      {last_resname:>3s} "
    f"{target_chain}{len(residue_keys):4d}"
)
output_lines.append("END")

outfile.write_text("\n".join(output_lines) + "\n", encoding="utf-8")
status.write_text("PASSED\n", encoding="utf-8")
PY

mv -f -- "$TMP_REPORT" "$REPORT"

STATUS="$(cat "$TMP_STATUS" 2>/dev/null || printf "FAILED")"

if [[ "$STATUS" != "PASSED" ]]; then
    rm -f -- "$OUTFILE" "$DISU_OUT"
    echo "Structure preparation failed." >&2
    echo "See report: $REPORT" >&2
    exit 1
fi

mv -f -- "$TMP_PDB" "$OUTFILE"
mv -f -- "$TMP_DISU" "$DISU_OUT"

echo "Amber PDB preparation completed successfully."
echo "Input PDB       : $INFILE"
echo "Selected chain  : $CHAIN"
echo "Cleaned PDB     : $OUTFILE"
echo "Disulfide file  : $DISU_OUT"
echo "Preparation log : $REPORT"
