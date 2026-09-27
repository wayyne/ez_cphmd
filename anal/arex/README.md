# AREX analysis

Tools for checking Amber asynchronous replica exchange (AREX) mixing and for comparing structural observables across pH-slot trajectories.

## Tools

- `arex_acceptance_rate.sh` — acceptance rate for each adjacent pH pair from `arex.history.tsv`.
- `roundtrip_check.py` — per-replica pH-slot occupancy, movement, and complete round-trip counts.
- `analyze_arex.sh` — residue-centric CPPTRAJ analysis of one or more pH-slot trajectories.
- `plot_arex_analysis.py` — summary tables and figures for `analyze_arex.sh`, with optional protonation-state conditioning.

## Exchange acceptance

```bash
./arex_acceptance_rate.sh arex.history.tsv \
    '4.5 5.0 5.5 6.0 6.5 7.0 7.5 8.0'
```

The history file is interpreted as native AREX `slot_of_replica` data with 1-based slot indices. The acceptance calculation assumes the neighboring exchange pattern alternates as:

```text
odd epoch   (1,2) (3,4) ...
even epoch  (2,3) (4,5) ...
```

The script assumes successive valid rows represent successive exchange epochs. If history rows are missing, inspect the file before interpreting the rates.

## Replica mixing and round trips

```bash
python3 roundtrip_check.py arex.history.tsv \
    '4.5 5.0 5.5 6.0 6.5 7.0 7.5 8.0'
```

For each replica the script reports:

- occupancy fraction at each pH slot
- number of slots visited
- minimum and maximum pH reached
- number of slot-changing moves
- average absolute slot jump when a move occurs
- complete endpoint round trips

A round trip is counted only after a replica reaches one endpoint of the pH ladder, reaches the opposite endpoint, and returns to the original endpoint.

## Residue-centric AREX analysis

`analyze_arex.sh` keeps each pH-slot trajectory separate and can calculate:

- backbone RMSD
- radius of gyration
- per-residue SASA
- side-chain SASA
- first- and second-shell water counts
- site-water hydrogen bonds
- site-protein hydrogen bonds

Example:

```bash
./analyze_arex.sh \
    -p system.parm7 \
    -t 'arex_ph*.nc' \
    -r all \
    --site-mode residue \
    --metrics sasa,watershell,hbond \
    --dt-ps 10 \
    -o arex_analysis
```

For a selected set of biological residue numbers:

```bash
./analyze_arex.sh \
    -p system.parm7 \
    -t 'arex_ph*.nc' \
    -r 7,15,18,35,48,52,66,87,101,119 \
    --residue-offset 0 \
    --dt-ps 10 \
    -o titratable_arex_analysis
```

Use `--residue-offset` when biological residue numbering differs from topology/CPPTRAJ numbering, for example because of an N-terminal cap.

The default local hydrogen-bond definition is a donor-acceptor distance of 3.0 Å or less and a donor-H-acceptor angle of at least 135 degrees. Both cutoffs are configurable.

The driver records the selected residues, masks, trajectory metadata, settings, CPPTRAJ input, and CPPTRAJ logs in the output directory.

## Summary tables and figures

`analyze_arex.sh` calls `plot_arex_analysis.py` automatically unless `--no-plot` is supplied.

It creates:

```text
residue_summary.tsv
global_summary.tsv
figures/
```

The summaries include means, sample standard deviations, and a block-based SEM diagnostic. The block SEM should be treated as a convergence diagnostic rather than proof of statistically independent samples.

The plotting script requires Python 3, NumPy, and Matplotlib.

## Protonation-state-conditioned structural analysis

Coupled titration analysis from `../titration/calc_titration.py` can write a frame-resolved state CSV. That file can be joined to the structural analysis by exact MD step:

```bash
python3 plot_arex_analysis.py \
    --analysis-dir arex_analysis \
    --state-file ../titration/coupled_anal/ff19_E35_D52_frame_states.csv \
    --first-frame-step 0 \
    --frame-step-interval 5000
```

Additional outputs include:

```text
dyad_frame_assignments.tsv
dyad_state_summary.tsv
dyad_state_counts.tsv
```

The join is intentionally strict. A coordinate frame whose MD step lies inside the state-file window must have an exact state assignment; nearest-neighbor state matching is not used.

## Requirements

- AmberTools / `cpptraj`
- Bash
- Python 3
- NumPy and Matplotlib for `plot_arex_analysis.py`
