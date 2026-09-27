# RMSF

`calc_rmsf.py` calculates per-residue RMSF from Amber trajectories using CPPTRAJ.

Supported modes are:

```text
backbone   heavy backbone atoms N, CA, C, O
sidechain  heavy side-chain atoms
total      all heavy atoms
```

Frames are fitted with heavy backbone atoms by default regardless of the RMSF measurement mode. This removes global translation/rotation without fitting away the side-chain motion being measured.

## Basic usage

```bash
python3 calc_rmsf.py \
    -p system.prmtop \
    -y prod1.nc prod2.nc \
    --rmsf-mode backbone \
    -o backbone_rmsf
```

Multiple trajectories should be supplied in chronological order when they are consecutive simulation segments.

## Side-chain RMSF

```bash
python3 calc_rmsf.py \
    -p system.prmtop \
    -y prod.nc \
    --selection-mask ':1-129' \
    --rmsf-mode sidechain \
    --highlight-residues 35 52 \
    -o sidechain_rmsf
```

Glycine has no conventional side-chain heavy atoms and may therefore be absent from side-chain RMSF output.

The default fit reference is the trajectory average. `--fit-reference first` or an external `--reference` can be used instead.

Residue labels can be adjusted with `--residue-offset`, or original PDB numbering can be requested with `--pdbres` where supported by CPPTRAJ.

## Output

For an output prefix such as `backbone_rmsf`, the script creates:

```text
backbone_rmsf.png
backbone_rmsf.csv
backbone_rmsf.summary.txt
backbone_rmsf.cpptraj.in
backbone_rmsf.cpptraj.log
backbone_rmsf.cpptraj.dat
```

The text summary includes masks, CPPTRAJ provenance, frame counts, RMSF statistics, and the highest-RMSF residues.

## Requirements

- AmberTools / `cpptraj`
- Python 3
- Matplotlib

Run `python3 calc_rmsf.py --help` for all options.
