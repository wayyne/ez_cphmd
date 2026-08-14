# `calc_distance.py`

Frame-resolved minimum-distance analysis for Amber trajectories using CPPTRAJ.

`calc_distance.py` is deliberately **state-agnostic**. It calculates structural observables from coordinates and writes a standardized per-frame table that can later be joined to protonation-state data, coupled-state assignments, or any other frame annotation.

The tool is intended to be one small component of a modular trajectory-analysis workflow:

```text
Amber topology + trajectories
            │
            ▼
     calc_distance.py
            │
            ▼
 frame-wise distance TSV
            │
            ├── use directly
            └── join to states with partition_by_state.py
```

It does **not** read lambda coordinates, assign protonation states, interpret coupled titration models, or calculate state-conditioned statistics.

## Features

- Calculates the minimum atom-to-atom distance between two arbitrary CPPTRAJ masks.
- Supports one or many named distance pairs in a single pass.
- Processes multiple trajectory files or quoted trajectory globs.
- Preserves trajectory identity and frame numbering.
- Optionally maps every coordinate frame to an absolute MD step.
- Extracts pH from common CpHMD/AREX trajectory names when available.
- Writes a long-format, frame-wise TSV suitable for downstream joins.
- Writes a provenance/settings file alongside the result.
- Keeps CPPTRAJ input, log, and intermediate series for inspection.

## Requirements

- Python 3.10 or newer
- NumPy
- Amber/AmberTools with `cpptraj` available on `PATH`, or a path supplied with `--cpptraj`

Install the Python dependency with, for example:

```bash
python3 -m pip install numpy
```

Confirm that CPPTRAJ is available:

```bash
cpptraj -V
```

## Basic usage

```bash
python3 calc_distance.py \
    -p system.parm7 \
    -t 'arex.ph*.nc' \
    --pair DYAD ':36@OE1,OE2' ':53@OD1,OD2' \
    --image-anchor '^1' \
    --first-frame-step 510000 \
    --frame-step-interval 10000 \
    -o dyad_distance_frames.tsv
```

`--pair` takes three arguments:

```text
--pair NAME MASK1 MASK2
```

`NAME` becomes the observable label in the output. `MASK1` and `MASK2` are ordinary CPPTRAJ masks and therefore use **topology/CPPTRAJ numbering**, not necessarily biological sequence numbering.

The example above calculates the minimum distance between any atom selected by `:36@OE1,OE2` and any atom selected by `:53@OD1,OD2`.

## Multiple distances in one run

Repeat `--pair` to calculate several observables from the same trajectories:

```bash
python3 calc_distance.py \
    -p system.parm7 \
    -t 'arex.ph*.nc' \
    --pair E35_D52 ':36@OE1,OE2' ':53@OD1,OD2' \
    --pair D52_N46 ':53@OD1,OD2' ':47@ND2' \
    --pair D52_N59 ':53@OD1,OD2' ':60@ND2' \
    --image-anchor '^1' \
    --first-frame-step 510000 \
    --frame-step-interval 10000 \
    -o distances.tsv
```

Pair names must be unique after sanitization. Characters outside letters, numbers, `_`, `.`, and `-` are converted to `_`.

## Trajectory input

`-t/--trajectory` may be repeated and may contain a quoted glob:

```bash
-t 'arex.ph*.nc'
```

or:

```bash
-t arex.ph5_0.nc \
-t arex.ph5_5.nc \
-t arex.ph6_0.nc
```

Matched files are deduplicated and naturally sorted before processing.

### pH extraction

The tool attempts to infer pH from the trajectory basename. Supported examples include:

```text
arex.ph0_5.nc   -> 0.5
arex.ph5_0.nc   -> 5.0
arex_ph7.0.nc   -> 7.0
pH_7.0.nc       -> 7.0
```

If no pH-like token is present, the output `pH` field is left blank. pH is metadata only; it does not affect the distance calculation.

## Periodic imaging

Use `--image-anchor` when the system requires CPPTRAJ imaging before measuring distances:

```bash
--image-anchor '^1'
```

This generates:

```text
autoimage anchor "^1"
```

before the distance actions.

Choose an anchor appropriate for the topology. In systems where the protein is not molecule 1, do not assume `^1` is correct.

## Frame and MD-step mapping

The output always records both:

- `analysis_row`: 1-based row number in the processed trajectory output
- `trajectory_frame`: the corresponding 1-based frame number in the original stored trajectory

When both of the following are supplied:

```bash
--first-frame-step 510000
--frame-step-interval 10000
```

the tool also records the absolute MD step for each analyzed coordinate frame:

```text
step = first_frame_step + (trajectory_frame - 1) * frame_step_interval
```

This mapping is particularly useful when the structural result will later be joined to CpHMD lambda or state assignments.

`--first-frame-step` and `--frame-step-interval` must be supplied together. The frame-step interval must be positive.

### `--start`, `--stop`, and `--stride`

Trajectory processing follows CPPTRAJ frame selection:

```bash
--start 1
--stop last
--stride 1
```

For example:

```bash
--start 101 --stride 5
```

causes the first analyzed row to correspond to original stored trajectory frame 101, and subsequent rows to frames 106, 111, and so on. MD-step mapping is calculated from the **original trajectory frame number**, not the compact analysis-row number.

## Output

The main output is a tab-separated, long-format table with one row per trajectory frame per named distance observable.

### Output columns

| Column | Description |
|---|---|
| `trajectory` | Trajectory label, currently the input filename stem |
| `pH` | pH parsed from the trajectory label, when available |
| `analysis_row` | 1-based row in the processed trajectory |
| `trajectory_frame` | 1-based frame number in the original stored trajectory |
| `step` | Absolute MD step, if frame-step mapping was supplied |
| `observable` | Sanitized name supplied to `--pair` |
| `metric` | Currently always `min_distance` |
| `selection1` | First CPPTRAJ mask |
| `selection2` | Second CPPTRAJ mask |
| `value` | Minimum atom-to-atom distance for the frame |
| `unit` | Currently `A` for Å |

Example:

```text
trajectory    pH   analysis_row  trajectory_frame  step    observable  metric        selection1       selection2       value   unit
arex.ph6_5    6.5  1             1                 510000  E35_D52     min_distance  :36@OE1,OE2     :53@OD1,OD2     7.214   A
arex.ph6_5    6.5  2             2                 520000  E35_D52     min_distance  :36@OE1,OE2     :53@OD1,OD2     7.083   A
```

## Provenance and working files

If the main output is:

```text
dyad_distance_frames.tsv
```

the tool also writes:

```text
dyad_distance_frames.tsv.settings.txt
```

The settings file records the resolved topology and trajectory paths, frame-selection settings, step mapping, image anchor, CPPTRAJ command, and pair definitions.

CPPTRAJ inputs and intermediate outputs are stored under `--workdir`, which defaults to:

```text
distance_work/
```

with one subdirectory per trajectory:

```text
distance_work/
└── arex.ph6_5/
    ├── distance.in
    ├── cpptraj.log
    ├── E35_D52.dat
    ├── D52_N46.dat
    └── D52_N59.dat
```

These files are useful for debugging and reproducibility.

## Caching and `--force`

Existing intermediate pair files are reused unless `--force` is specified.

```bash
--force
```

forces CPPTRAJ to rerun and overwrite/recalculate the requested series.

**Important:** the current cache is file-name based. If you reuse the same pair name but change its masks, topology, trajectory processing options, or other CPPTRAJ-relevant settings, use `--force` or a fresh `--workdir` to avoid reusing stale intermediate data.

## Integration with state analysis

The intended downstream interface is the frame identity, not row order alone. For example:

```bash
python3 ../partition/partition_by_state.py \
    --observations dyad_distance_frames.tsv \
    --states dyad_frame_assignments.tsv \
    --join-on trajectory,trajectory_frame \
    --group-by pH,observable,metric,selection1,selection2,unit \
    -o dyad_distance_by_state.tsv
```

When both tables contain `step`, `partition_by_state.py` can independently validate that the structural frame and state assignment refer to the same absolute MD step.

This separation is intentional: `calc_distance.py` remains reusable for conventional MD, independent CpHMD titration, coupled CpHMD models, or any future state model.

## Validation and failure behavior

The tool fails rather than silently continuing when it detects conditions that could invalidate the output, including:

- missing topology;
- trajectory specifications that match no files;
- invalid `start`, `stop`, or `stride` values;
- only one of the two MD-step mapping arguments being supplied;
- duplicate pair names after sanitization;
- CPPTRAJ failure;
- missing or unreadable CPPTRAJ output;
- non-finite or negative distance values;
- inconsistent frame counts among distance observables from the same trajectory.

Check the per-trajectory `cpptraj.log` when CPPTRAJ fails.

## Current scope and limitations

`calc_distance.py` is intentionally small. In its current form:

- it calculates **minimum atom-to-atom distances only**;
- masks use topology/CPPTRAJ numbering directly;
- it does not translate biological residue numbers or apply residue offsets;
- it does not assign protonation states;
- it does not perform coupled-state analysis;
- it does not calculate summary statistics;
- it assumes the same `--first-frame-step` and `--frame-step-interval` mapping applies to every trajectory supplied in one invocation;
- it does not attempt to infer whether imaging is required or which image anchor is physically appropriate.

These are deliberate boundaries. Residue-centric selection helpers and additional geometric observables can be added later without changing the frame-wise output contract.

## Reproducibility recommendations

For production analyses:

1. Keep the generated `.settings.txt` file with the results.
2. Retain the CPPTRAJ input and log files in the working directory.
3. Record the Amber/AmberTools version used for CPPTRAJ.
4. Supply absolute MD-step mapping when results will be joined to CpHMD state data.
5. Validate the first and last frame-to-step mappings against the original simulation output.
6. Use `--force` or a fresh work directory whenever pair definitions or trajectory-processing settings change.

## Command-line reference

Run:

```bash
python3 calc_distance.py --help
```

Key options:

| Option | Purpose |
|---|---|
| `-p, --topology FILE` | Amber topology; required |
| `-t, --trajectory SPEC` | Trajectory or quoted glob; required, repeatable |
| `--pair NAME MASK1 MASK2` | Named minimum-distance pair; required, repeatable |
| `--image-anchor MASK` | Optional CPPTRAJ `autoimage` anchor |
| `--start N` | First stored input frame; default `1` |
| `--stop N\|last` | Last stored input frame; default `last` |
| `--stride N` | Input-frame stride; default `1` |
| `--first-frame-step N` | Absolute MD step of stored trajectory frame 1 |
| `--frame-step-interval N` | MD-step spacing between stored frames |
| `--cpptraj CMD` | CPPTRAJ executable/command; default `cpptraj` |
| `--workdir DIR` | CPPTRAJ working directory; default `distance_work` |
| `-o, --output FILE` | Frame-wise TSV; default `distance_frames.tsv` |
| `--force` | Recalculate existing CPPTRAJ series |

