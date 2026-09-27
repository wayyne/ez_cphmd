#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# EXTRACT, NORMALIZE, AND VALIDATE ONE PROTEIN CHAIN FROM A PDB FILE
###############################################################################

INFILE=""
PDBID=""
CHAIN=""
OUTFILE=""
DISU_OUT=""
REPORT=""
SEGID="PROA"

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
  --in <file>          Original downloaded PDB file
  --pdbid <PDBID>      Four-character PDB identifier
  --chain <CHAIN>      Single chain identifier to retain
  --out <file>         Cleaned output PDB
  --disu-out <file>    Output file containing mapped disulfide pairs
  --report <file>      Structure-preparation report

Current supported scope:
  - One selected protein chain
  - Zero or one MODEL record
  - Standard amino-acid residues
  - No insertion codes
  - Consecutive residue numbering
  - Complete backbone atoms: N, CA, C, O
  - Blank alternate locations preferred over alternate location A
  - Existing hydrogens and OXT atoms removed
  - HIS, HSD, and HSE normalized to HSP
  - Same-chain SSBOND records retained and renumbered
  - CHARMM segment ID written as PROA
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

###############################################################################
# INPUT VALIDATION
###############################################################################

[[ -n "$INFILE" ]] || die "--in is required"
[[ -n "$PDBID" ]] || die "--pdbid is required"
[[ -n "$CHAIN" ]] || die "--chain is required"
[[ -n "$OUTFILE" ]] || die "--out is required"
[[ -n "$DISU_OUT" ]] || die "--disu-out is required"
[[ -n "$REPORT" ]] || die "--report is required"

[[ -f "$INFILE" ]] || die "input PDB file not found: $INFILE"
[[ -s "$INFILE" ]] || die "input PDB file is empty: $INFILE"

PDBID="$(echo "$PDBID" | tr '[:upper:]' '[:lower:]')"
CHAIN="$(echo "$CHAIN" | tr '[:lower:]' '[:upper:]')"

if ! [[ "$PDBID" =~ ^[[:alnum:]]{4}$ ]]; then
    die "PDB ID must contain exactly four alphanumeric characters"
fi

if ! [[ "$CHAIN" =~ ^[[:alnum:]]$ ]]; then
    die "chain must contain exactly one alphanumeric character"
fi

for output_path in "$OUTFILE" "$DISU_OUT" "$REPORT"; do
    output_dir="$(dirname -- "$output_path")"

    if [[ ! -d "$output_dir" ]]; then
        die "output directory does not exist: $output_dir"
    fi
done

if [[ "$INFILE" == "$OUTFILE" ]]; then
    die "input and cleaned output filenames must be different"
fi

if ! command -v awk >/dev/null 2>&1; then
    die "awk is required"
fi

###############################################################################
# MODEL VALIDATION
###############################################################################

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

if [[ "$MODEL_COUNT" -gt 1 ]]; then
    die "multiple structural models are not supported; found ${MODEL_COUNT}"
fi

###############################################################################
# TEMPORARY OUTPUTS
###############################################################################

TMPDIR="$(mktemp -d "${TMPDIR:-/tmp}/prep_pdb.XXXXXX")"

cleanup() {
    rm -rf -- "$TMPDIR"
}

trap cleanup EXIT

TMP_PDB="${TMPDIR}/clean.pdb"
TMP_DISU="${TMPDIR}/disulfides.dat"
TMP_REPORT="${TMPDIR}/prep_report.txt"
TMP_STATUS="${TMPDIR}/status"

###############################################################################
# STRUCTURE EXTRACTION AND VALIDATION
###############################################################################

LC_ALL=C awk \
    -v target_chain="$CHAIN" \
    -v pdbid="$PDBID" \
    -v input_file="$INFILE" \
    -v model_count="$MODEL_COUNT" \
    -v output_file="$TMP_PDB" \
    -v disu_file="$TMP_DISU" \
    -v report_file="$TMP_REPORT" \
    -v status_file="$TMP_STATUS" \
    -v output_segid="$SEGID" '
BEGIN {
    standard["ALA"] = "A"
    standard["ARG"] = "R"
    standard["ASN"] = "N"
    standard["ASP"] = "D"
    standard["CYS"] = "C"
    standard["GLN"] = "Q"
    standard["GLU"] = "E"
    standard["GLY"] = "G"
    standard["HIS"] = "H"
    standard["HSD"] = "H"
    standard["HSE"] = "H"
    standard["HSP"] = "H"
    standard["ILE"] = "I"
    standard["LEU"] = "L"
    standard["LYS"] = "K"
    standard["MET"] = "M"
    standard["PHE"] = "F"
    standard["PRO"] = "P"
    standard["SER"] = "S"
    standard["THR"] = "T"
    standard["TRP"] = "W"
    standard["TYR"] = "Y"
    standard["VAL"] = "V"
}

function trim(value) {
    sub(/^[[:space:]]+/, "", value)
    sub(/[[:space:]]+$/, "", value)
    return value
}

function upper(value) {
    return toupper(value)
}

function add_error(message) {
    if (!(message in error_seen)) {
        error_seen[message] = 1
        errors[++error_count] = message
    }
}

function is_integer(value) {
    return value ~ /^-?[0-9]+$/
}

function is_number(value) {
    return value ~ /^[-+]?[0-9]*[.]?[0-9]+([eE][-+]?[0-9]+)?$/
}

function is_hydrogen(line, atom_name, element, test_name) {
    element = upper(trim(substr(line, 77, 2)))

    if (element == "H" || element == "D") {
        return 1
    }

    test_name = upper(atom_name)
    gsub(/[[:space:]]/, "", test_name)
    gsub(/[0-9]/, "", test_name)

    return substr(test_name, 1, 1) == "H" ||
           substr(test_name, 1, 1) == "D"
}

function normalize_residue(resname) {
    if (resname == "HIS" ||
        resname == "HSD" ||
        resname == "HSE" ||
        resname == "HSP") {
        return "HSP"
    }

    return resname
}

function residue_key(resseq, icode) {
    return resseq "|" icode
}

function atom_key(reskey, atom_name) {
    return reskey SUBSEP atom_name
}

function canonical_pair(first, second) {
    if ((first + 0) < (second + 0)) {
        return first ":" second
    }

    return second ":" first
}

substr($0, 1, 6) == "SSBOND" {
    ssbond_count++

    ss_chain1[ssbond_count] = trim(substr($0, 16, 1))
    ss_seq1[ssbond_count] = trim(substr($0, 18, 4))
    ss_icode1[ssbond_count] = substr($0, 22, 1)

    ss_chain2[ssbond_count] = trim(substr($0, 30, 1))
    ss_seq2[ssbond_count] = trim(substr($0, 32, 4))
    ss_icode2[ssbond_count] = substr($0, 36, 1)

    next
}

substr($0, 1, 6) == "HETATM" {
    if (substr($0, 22, 1) == target_chain) {
        selected_hetatm_records++
    }

    next
}

substr($0, 1, 6) != "ATOM  " {
    next
}

{
    chain = substr($0, 22, 1)

    if (chain != target_chain) {
        next
    }

    selected_atom_records++

    raw_resseq = trim(substr($0, 23, 4))
    icode = substr($0, 27, 1)
    resname = upper(trim(substr($0, 18, 3)))
    atom_name = upper(trim(substr($0, 13, 4)))
    altloc = substr($0, 17, 1)

    if (!is_integer(raw_resseq)) {
        add_error("Selected chain contains a non-integer residue number: " raw_resseq)
        next
    }

    resseq = raw_resseq + 0
    reskey = residue_key(resseq, icode)

    if (icode != " ") {
        insertion_codes[reskey] = 1
    }

    if (!(resname in standard)) {
        unsupported_residues[reskey] = resname
        next
    }

    normalized_resname = normalize_residue(resname)

    if (!(reskey in residue_seen)) {
        residue_seen[reskey] = 1
        residue_order[++residue_count] = reskey
        original_resseq[reskey] = resseq
        original_icode[reskey] = icode
        original_resname[reskey] = resname
        normalized_residue[reskey] = normalized_resname

        if (resname != normalized_resname) {
            histidine_residue_count++
        }
    } else {
        if (original_resname[reskey] != resname) {
            add_error("Residue " resseq " has conflicting residue names: " original_resname[reskey] " and " resname)
        }

        if (last_reskey != "" &&
            last_reskey != reskey &&
            (reskey in completed_residue)) {
            add_error("Residue records are noncontiguous for original residue " resseq)
        }
    }

    if (last_reskey != "" && last_reskey != reskey) {
        completed_residue[last_reskey] = 1
    }

    last_reskey = reskey

    if (is_hydrogen($0, atom_name)) {
        removed_hydrogen_count++
        next
    }

    if (atom_name == "OXT") {
        removed_oxt_count++
        next
    }

    akey = atom_key(reskey, atom_name)

    if (!(akey in atom_order_seen)) {
        atom_order_seen[akey] = 1
        atom_order[++atom_order_count] = akey
        atom_residue[akey] = reskey
        atom_names[akey] = atom_name
    }

    if (altloc == " ") {
        priority = 2
    } else if (altloc == "A") {
        priority = 1
    } else {
        unsupported_altloc[akey] = altloc
        next
    }

    if ((akey in chosen_priority) &&
        chosen_priority[akey] == priority) {
        add_error("Duplicate coordinate records for residue " resseq ", atom " atom_name ", alternate location " (altloc == " " ? "blank" : altloc))
        next
    }

    if (!(akey in chosen_priority) ||
        priority > chosen_priority[akey]) {
        chosen_priority[akey] = priority
        chosen_line[akey] = $0
        chosen_altloc[akey] = altloc
    }
}

END {
    if (selected_atom_records == 0) {
        add_error("No ATOM records were found for requested chain " target_chain)
    }

    if (residue_count == 0) {
        add_error("No supported amino-acid residues were found for chain " target_chain)
    }

    for (reskey in insertion_codes) {
        add_error("Insertion code found at original residue " original_resseq[reskey] original_icode[reskey])
    }

    for (reskey in unsupported_residues) {
        add_error("Unsupported residue " unsupported_residues[reskey] " at original residue " original_resseq[reskey] original_icode[reskey])
    }

    for (i = 1; i <= residue_count; i++) {
        reskey = residue_order[i]
        clean_resseq[reskey] = i
        sequence = sequence standard[original_resname[reskey]]

        if (i > 1) {
            previous_key = residue_order[i - 1]
            previous_number = original_resseq[previous_key]
            current_number = original_resseq[reskey]

            if (current_number != previous_number + 1) {
                add_error("Nonconsecutive residue numbering between original " "residues " previous_number " and " current_number)
            }
        }
    }

    for (i = 1; i <= atom_order_count; i++) {
        akey = atom_order[i]

        if (!(akey in chosen_line)) {
            reskey = atom_residue[akey]
            atom_name = atom_names[akey]

            add_error("No blank or A alternate-location coordinate exists for " "original residue " original_resseq[reskey] ", atom " atom_name)
            continue
        }

        line = chosen_line[akey]

        x_text = trim(substr(line, 31, 8))
        y_text = trim(substr(line, 39, 8))
        z_text = trim(substr(line, 47, 8))

        if (!is_number(x_text) ||
            !is_number(y_text) ||
            !is_number(z_text)) {
            reskey = atom_residue[akey]

            add_error("Invalid coordinate found for original residue " original_resseq[reskey] ", atom " atom_names[akey])
        }

        if (chosen_altloc[akey] == "A") {
            selected_altloc_a_count++
        } else {
            selected_blank_altloc_count++
        }
    }

    required_backbone[1] = "N"
    required_backbone[2] = "CA"
    required_backbone[3] = "C"
    required_backbone[4] = "O"

    for (i = 1; i <= residue_count; i++) {
        reskey = residue_order[i]

        for (j = 1; j <= 4; j++) {
            atom_name = required_backbone[j]
            akey = atom_key(reskey, atom_name)

            if (!(akey in chosen_line)) {
                add_error("Missing backbone atom " atom_name " at original residue " original_resseq[reskey] " " original_resname[reskey])
            }
        }
    }

    for (i = 1; i <= ssbond_count; i++) {
        chain1 = ss_chain1[i]
        chain2 = ss_chain2[i]

        involves_selected = \
            (chain1 == target_chain || chain2 == target_chain)

        if (!involves_selected) {
            continue
        }

        if (chain1 != target_chain ||
            chain2 != target_chain) {
            add_error("SSBOND involving selected chain " target_chain " connects to another chain")
            continue
        }

        if (ss_icode1[i] != " " ||
            ss_icode2[i] != " ") {
            add_error("SSBOND uses an insertion-coded cysteine, " "which is unsupported")
            continue
        }

        if (!is_integer(ss_seq1[i]) ||
            !is_integer(ss_seq2[i])) {
            add_error("SSBOND contains a non-integer residue number")
            continue
        }

        key1 = residue_key(ss_seq1[i] + 0, " ")
        key2 = residue_key(ss_seq2[i] + 0, " ")

        if (!(key1 in clean_resseq) ||
            !(key2 in clean_resseq)) {
            add_error("SSBOND references a residue missing from selected chain " target_chain ": " ss_seq1[i] " or " ss_seq2[i])
            continue
        }

        if (normalized_residue[key1] != "CYS" ||
            normalized_residue[key2] != "CYS") {
            add_error("SSBOND does not connect two CYS residues: " ss_seq1[i] " and " ss_seq2[i])
            continue
        }

        clean1 = clean_resseq[key1]
        clean2 = clean_resseq[key2]
        pair_key = canonical_pair(clean1, clean2)

        if (!(pair_key in disulfide_seen)) {
            disulfide_seen[pair_key] = 1

            disulfide_first[++disulfide_count] = \
                (clean1 < clean2 ? clean1 : clean2)

            disulfide_second[disulfide_count] = \
                (clean1 < clean2 ? clean2 : clean1)

            disulfide_orig_first[disulfide_count] = \
                original_resseq[key1]

            disulfide_orig_second[disulfide_count] = \
                original_resseq[key2]
        }
    }

    print "============================================================" \
        > report_file

    print "PDB structure preparation report" \
        > report_file

    print "============================================================" \
        > report_file

    print "PDB ID                  : " pdbid \
        > report_file

    print "Input file              : " input_file \
        > report_file

    print "Selected chain          : " target_chain \
        > report_file

    print "CHARMM segment ID       : " output_segid \
        > report_file

    print "MODEL records           : " model_count \
        > report_file

    print "Selected ATOM records   : " selected_atom_records \
        > report_file

    print "Selected HETATM records : " selected_hetatm_records \
        > report_file

    print "Clean residue count     : " residue_count \
        > report_file

    print "Removed hydrogens       : " removed_hydrogen_count \
        > report_file

    print "Removed OXT atoms       : " removed_oxt_count \
        > report_file

    print "Normalized histidines   : " histidine_residue_count \
        > report_file

    print "Blank-altloc atoms      : " selected_blank_altloc_count \
        > report_file

    print "Altloc-A atoms          : " selected_altloc_a_count \
        > report_file

    print "Disulfide pairs         : " disulfide_count \
        > report_file

    print "Sequence                : " sequence \
        > report_file

    print "" > report_file

    print "Residue mapping:" > report_file
    print "  Clean  Original  Input  Output" > report_file

    for (i = 1; i <= residue_count; i++) {
        reskey = residue_order[i]

        printf "  %-6d %-9d %-6s %-6s\n", \
            clean_resseq[reskey], \
            original_resseq[reskey], \
            original_resname[reskey], \
            normalized_residue[reskey] \
            > report_file
    }

    print "" > report_file
    print "Disulfide mapping:" > report_file

    if (disulfide_count == 0) {
        print "  None" > report_file
    } else {
        for (i = 1; i <= disulfide_count; i++) {
            printf "  Original %d-%d -> Clean %d-%d\n", \
                disulfide_orig_first[i], \
                disulfide_orig_second[i], \
                disulfide_first[i], \
                disulfide_second[i] \
                > report_file
        }
    }

    if (error_count > 0) {
        print "" > report_file
        print "Status: FAILED" > report_file
        print "" > report_file
        print "Errors:" > report_file

        for (i = 1; i <= error_count; i++) {
            print "  - " errors[i] > report_file
        }

        print "FAILED" > status_file

        close(report_file)
        close(status_file)

        exit
    }

    print "" > report_file
    print "Status: PASSED" > report_file

    close(report_file)

    print "# CHARMM DISU residue pairs" > disu_file
    print "# Columns: clean_residue_1 clean_residue_2" > disu_file

    for (i = 1; i <= disulfide_count; i++) {
        print disulfide_first[i], disulfide_second[i] \
            > disu_file
    }

    close(disu_file)

    print "REMARK   Generated by prep_pdb.sh" \
        > output_file

    print "REMARK   Source PDB ID " toupper(pdbid) \
        > output_file

    print "REMARK   Original chain " target_chain \
        > output_file

    print "REMARK   CHARMM segment ID " output_segid \
        > output_file

    print "REMARK   Residues renumbered consecutively from 1" \
        > output_file

    print "REMARK   Histidines normalized to HSP" \
        > output_file

    serial = 0

    for (i = 1; i <= atom_order_count; i++) {
        akey = atom_order[i]

        if (!(akey in chosen_line)) {
            continue
        }

        line = chosen_line[akey]
        reskey = atom_residue[akey]

        atom_field = substr(line, 13, 4)
        resname = normalized_residue[reskey]
        new_resseq = clean_resseq[reskey]

        x = trim(substr(line, 31, 8)) + 0.0
        y = trim(substr(line, 39, 8)) + 0.0
        z = trim(substr(line, 47, 8)) + 0.0

        occupancy_text = trim(substr(line, 55, 6))
        bfactor_text = trim(substr(line, 61, 6))

        occupancy = \
            (is_number(occupancy_text) ? \
                occupancy_text + 0.0 : 1.0)

        bfactor = \
            (is_number(bfactor_text) ? \
                bfactor_text + 0.0 : 0.0)

        element = upper(trim(substr(line, 77, 2)))

        if (element == "") {
            element = atom_names[akey]
            gsub(/[0-9]/, "", element)
            element = substr(element, 1, 1)
        }

        serial++

        # PDB columns:
        #   22       = chain ID
        #   23-26    = residue number
        #   73-76    = CHARMM segment ID
        #   77-78    = element
        printf \
            "ATOM  %5d %s %3s %1s%4d    " \
            "%8.3f%8.3f%8.3f%6.2f%6.2f" \
            "      %-4s%2s\n", \
            serial, \
            atom_field, \
            resname, \
            target_chain, \
            new_resseq, \
            x, \
            y, \
            z, \
            occupancy, \
            bfactor, \
            output_segid, \
            element \
            > output_file
    }

    printf "TER   %5d      %3s %1s%4d\n", \
        serial + 1, \
        normalized_residue[residue_order[residue_count]], \
        target_chain, \
        residue_count \
        > output_file

    print "END" > output_file

    close(output_file)

    print "PASSED" > status_file

    close(status_file)
}
' "$INFILE"

###############################################################################
# INSTALL REPORT AND CHECK STATUS
###############################################################################

mv -f -- "$TMP_REPORT" "$REPORT"

STATUS="$(cat "$TMP_STATUS" 2>/dev/null || echo "FAILED")"

if [[ "$STATUS" != "PASSED" ]]; then
    rm -f -- "$OUTFILE" "$DISU_OUT"

    echo "Structure preparation failed." >&2
    echo "See report: $REPORT" >&2

    exit 1
fi

mv -f -- "$TMP_PDB" "$OUTFILE"
mv -f -- "$TMP_DISU" "$DISU_OUT"

###############################################################################
# COMPLETION SUMMARY
###############################################################################

echo "PDB preparation completed successfully."
echo "Input PDB       : $INFILE"
echo "Selected chain  : $CHAIN"
echo "CHARMM segid    : $SEGID"
echo "Cleaned PDB     : $OUTFILE"
echo "Disulfide file  : $DISU_OUT"
echo "Preparation log : $REPORT"
