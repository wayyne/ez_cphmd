# RMSD

`calc_rmsd.py` calculates Amber trajectory RMSD with separate masks for coordinate fitting and RMSD measurement.

The workflow is:

```text
autoimage (optional)
    -> fit coordinates using --fit-mask
    -> measure RMSD using --rmsd-mask without a second fit
```

This makes it possible, for example, to align on the protein backbone while measuring RMSD for a selected domain or region.

## Basic usage

```bash
python3 calc_rmsd.py \
    -p system.prmtop \
    -y prod1.nc prod2.nc \
    -o backbone_rmsd
```

Multiple trajectories are interpreted in the order supplied and should be listed chronologically when they are consecutive simulation segments.

## Separate fit and measurement masks

```bash
python3 calc_rmsd.py \
    -p system.prmtop \
    -y prod.nc \
    --fit-mask ':1-129@N,CA,C,O' \
    --rmsd-mask ':35-52@N,CA,C,O' \
    -o region_rmsd
```

## Time axis

Supply the saved-frame spacing directly:

```bash
--time-per-frame-ps 10
```

or let the script read `dt` and `ntwx` from an Amber mdin/mdout-style text file:

```bash
--mdin production.in
```

If no timing information is supplied, plots use processed frame number.

## Output

For an output prefix such as `backbone_rmsd`, the script creates:

```text
backbone_rmsd.png
backbone_rmsd.csv
backbone_rmsd.summary.txt
backbone_rmsd.cpptraj.in
backbone_rmsd.cpptraj.log
backbone_rmsd.cpptraj.dat
```

The CSV retains source trajectory and source frame provenance. The script also warns about adjacent RMSD jumps larger than the configurable `--jump-threshold`.

## Requirements

- AmberTools / `cpptraj`
- Python 3
- Matplotlib

Run `python3 calc_rmsd.py --help` for all options.
