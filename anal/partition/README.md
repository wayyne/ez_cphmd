# `partition_by_state.py`

Generic post-processing for joining frame-wise observables to state assignments and producing state-conditioned summaries.

`partition_by_state.py` is intentionally **observable-agnostic** and **state-model-agnostic**. It does not calculate SASA, distances, hydrogen bonds, hydration, or any other structural property. It also does not interpret CpHMD lambda coordinates or construct coupled protonation states.

Instead, it joins two already-defined data products:

1. a frame-wise observation table; and
2. a frame/state assignment table.

```text
frame-wise observations        frame-wise states
         │                           │
         └────────────┬──────────────┘
                      ▼
            partition_by_state.py
                      │
            ┌─────────┴─────────┐
            ▼                   ▼
      joined audit TSV    state summary TSV
```

This separation allows the same post-processing code to be reused for distances, SASA, fSASA-derived inputs, hydrogen bonds, hydration, or future observables, regardless of whether the underlying CpHMD model uses independent titration coordinates or a coupled titration model.

## Features

- Exact TSV join on configurable frame identifiers.
- Default join on `trajectory,trajectory_frame`.
- Supports more specific keys such as `trajectory,trajectory_frame,resid` for per-residue states.
- Validates absolute MD step when both inputs provide a populated `step` column.
- Requires unique state-table join keys.
- Supports configurable state, value, `clean`, and `matched` columns.
- Computes per-state frame count, mean, sample SD, median, 10th/90th percentiles, minimum, and maximum.
- Computes a block-based SEM while preserving the original trajectory timeline.
- Always summarizes each original trajectory independently.
- Writes both the frame-wise joined table and the reduced state summary.

## Requirements

- Python 3.10 or newer
- NumPy

No Amber or CPPTRAJ installation is required by this tool itself.

Install the Python dependency with, for example:

```bash
python3 -m pip install numpy
```

## Basic usage

For a coupled-state table with one state assignment per trajectory frame:

```bash
python3 partition_by_state.py \
    --observations ../distance/dyad_distance_frames.tsv \
    --states dyad_frame_assignments.tsv \
    --join-on trajectory,trajectory_frame \
    --group-by pH,observable,metric,selection1,selection2,unit \
    --blocks 5 \
    --joined-output dyad_distance_with_states.tsv \
    -o dyad_distance_by_state.tsv
```

The two outputs serve different purposes:

- `dyad_distance_with_states.tsv` is the frame-level audit trail after the join.
- `dyad_distance_by_state.tsv` contains the state-conditioned summary statistics.

## Design: state model and observable are separate

The program does not have `stateless`, `protonation`, or `coupled` analysis modes. Those concepts belong to the upstream state representation.

### Stateless analysis

Do not run this tool. Use the structural observable table directly.

### Independent protonation states

Use a state table with a residue identifier and one local state per residue/frame, for example:

```text
trajectory    trajectory_frame    step      resid    state_label    protonation    clean    matched
arex.ph6_5    1                   510000    35       H              H              1        1
arex.ph6_5    1                   510000    52       deprot         deprot         1        1
```

Join with:

```bash
--join-on trajectory,trajectory_frame,resid
```

The observation table must contain the same `resid` field.

### Coupled protonation states

Use one joint state per coupled group/frame, for example:

```text
trajectory    trajectory_frame    step      state_label    clean    matched
arex.ph6_5    1                   510000    S1_E35         1        1
arex.ph6_5    2                   520000    S1_D52         1        1
```

Join with the default:

```bash
--join-on trajectory,trajectory_frame
```

The partitioner does not need to know what `S1_E35` means. It simply treats `state_label` as the categorical state used to partition the observable.

If a downstream derived quantity needs the **local protonation state of an individual residue**—for example selecting the correct protonation-specific reference for fSASA—that local-state projection should be present in, or derived from, the state data before or during the appropriate derived-analysis step.

## Observation-table contract

At minimum, the observation table must contain:

- every column named in `--join-on`;
- the value column specified by `--value-column` (default: `value`); and
- `analysis_row` or `trajectory_frame` so the original timeline can be reconstructed for blocking.

A canonical frame-wise observable table might look like:

```text
trajectory    pH   analysis_row  trajectory_frame  step    observable  metric        value  unit
arex.ph6_5    6.5  1             1                 510000  E35_D52     min_distance  7.214  A
arex.ph6_5    6.5  2             2                 520000  E35_D52     min_distance  7.083  A
```

`calc_distance.py` produces this style of output directly.

The observation table may contain additional descriptor columns. These can be preserved as summary grouping variables with `--group-by`.

## State-table contract

At minimum, the state table must contain:

- every column named in `--join-on`; and
- the categorical state column specified by `--state-column` (default: `state_label`).

The following columns are optional but receive special handling when present:

| Column | Default role |
|---|---|
| `step` | Independent validation of frame identity after the join |
| `clean` | State-quality filter unless `--include-unclean` is supplied |
| `matched` | Frame/state-match filter unless `--include-unmatched` is supplied |

Column names can be changed with:

```bash
--state-column STATE
--clean-column CLEAN
--matched-column MATCHED
```

Accepted truthy values for `clean` and `matched` are case-insensitive:

```text
1, true, yes, y
```

## Exact joins and frame identity

The default join key is:

```text
trajectory,trajectory_frame
```

The state table must contain **exactly one row per join key** after the selected state filters are applied. Duplicate state keys are treated as an error.

The tool does not perform nearest-neighbor frame matching, time interpolation, or fuzzy matching.

### MD-step validation

If both the observation row and matched state row contain a non-empty `step` value, the two steps must agree exactly as integers.

A disagreement causes the program to fail rather than silently accept a potentially shifted frame/state mapping.

For state-conditioned CpHMD analysis, retaining `step` in both upstream products is strongly recommended even when the primary join key is trajectory/frame.

## Clean and matched states

By default, when the state table contains `clean` and/or `matched` columns, rows that are not truthy are removed from the state lookup before joining.

Use:

```bash
--include-unclean
```

to retain state rows not marked clean, and:

```bash
--include-unmatched
```

to retain rows not marked matched.

### Important current behavior

The current implementation requires **every observation row to find a state row after these filters are applied**. Therefore, if the observation table contains frames labeled `mixed`, `outside_state_window`, or otherwise unclean/unmatched in the state table, the default run will fail with a missing-state error.

For the current implementation, either:

1. provide an observation table restricted to the same clean/matched frame set; or
2. use `--include-unclean --include-unmatched` when you intentionally want those rows present in the join and subsequent state categories.

This behavior is intentionally documented because it affects how mixed and out-of-window CpHMD frames are handled. A future refinement may separate “join all auditable rows” from “include rows in statistical summaries.”

## Grouping behavior

`--group-by` controls which observation descriptors define an independent summary group.

Example:

```bash
--group-by pH,observable,metric,selection1,selection2,unit
```

`trajectory` is always added to the grouping key if not already present. This is deliberate: blocking and state summaries are performed within each original trajectory rather than across concatenated pH slots or independent runs.

If `--group-by` is omitted, the tool derives grouping columns from the observation table by excluding frame-identification columns, the value column, and copied state metadata.

For production analysis, an explicit `--group-by` is recommended because it makes the intended reduction unambiguous.

## State-conditioned statistics

For each trajectory/group/state combination, the tool reports:

| Field | Definition |
|---|---|
| `nframes` | Number of finite observation values assigned to the state |
| `mean` | Arithmetic mean over state frames |
| `sd` | Sample standard deviation (`ddof=1`) |
| `block_sem` | SEM of state-specific means across occupied temporal blocks |
| `nblocks` | Number of temporal blocks that contained at least one frame in the state |
| `median` | Median over state frames |
| `q10` | 10th percentile |
| `q90` | 90th percentile |
| `min` | Minimum |
| `max` | Maximum |

Non-finite observation values are excluded from the numerical summaries.

## Block-SEM definition

The block calculation is designed for protonation states that may be visited intermittently.

For each trajectory/group:

1. rows are sorted by the original trajectory timeline (`analysis_row`, falling back to `trajectory_frame`);
2. the **full trajectory timeline** is divided into `N` contiguous blocks;
3. frames belonging to the target state are selected within each block;
4. a state-specific mean is calculated for each block that contains at least one frame in that state;
5. `block_sem` is the sample standard deviation of those block means divided by the square root of the number of occupied blocks.

This is intentionally different from concatenating all frames in one state and only then dividing the concatenated state trajectory into blocks.

The default is:

```bash
--blocks 5
```

`nblocks` in the output may be smaller than the requested block count when a state is absent from some temporal blocks.

### Statistical interpretation

`block_sem` is a **block-based diagnostic**, not an assertion that the resulting block means are statistically independent. For manuscript-quality uncertainty estimates, block length should be comfortably larger than the relevant correlation time, and block-size/convergence sensitivity should be checked explicitly.

## Outputs

### Joined audit table

Set with:

```bash
--joined-output observations_with_states.tsv
```

This contains every successfully joined observation row plus state-table metadata. If a state column name conflicts with an existing observation column, the copied state version is prefixed with `state_`.

The tool also creates a stable `partition_state` field internally/in the joined output so the requested state column can be used consistently during summary generation.

This table is the preferred place to audit frame-to-state assignments before interpreting summary statistics.

### Summary table

Set with:

```bash
-o state_partition_summary.tsv
```

The table contains the selected grouping descriptors, state label, frame count, and summary statistics described above.

## Examples

### Coupled dyad distance

```bash
python3 partition_by_state.py \
    --observations E35_D52_distance_frames.tsv \
    --states dyad_frame_assignments.tsv \
    --join-on trajectory,trajectory_frame \
    --group-by pH,observable,metric,selection1,selection2,unit \
    --blocks 5 \
    --joined-output E35_D52_distance_with_states.tsv \
    -o E35_D52_distance_by_state.tsv
```

### Alternative state column

```bash
python3 partition_by_state.py \
    --observations observable.tsv \
    --states states.tsv \
    --state-column protonation \
    --join-on trajectory,trajectory_frame,resid \
    --group-by resid,observable,unit \
    -o observable_by_protonation.tsv
```

### Alternative value column

```bash
--value-column sasa_A2
```

This allows a compatible table to be analyzed without renaming its numeric observation field to `value`.

### Include all state-table rows

```bash
python3 partition_by_state.py \
    --observations observable.tsv \
    --states states.tsv \
    --include-unclean \
    --include-unmatched \
    -o observable_all_states.tsv
```

Use this only when the meaning of the additional state labels is understood and they should appear in the statistical partitioning.

## Validation and failure behavior

The tool fails rather than silently guessing when it encounters ambiguous or inconsistent input, including:

- empty observation or state tables;
- missing value or state columns;
- missing join-key columns;
- non-integral frame/step/residue identifiers where exact integers are required;
- duplicate state-table keys after filtering;
- no state rows remaining after filtering;
- observation rows with no corresponding state assignment;
- disagreement between observation and state MD steps;
- missing usable timeline identifier for block statistics.

This behavior is intentional for frame-resolved CpHMD analysis, where a one-frame offset can invalidate a state-conditioned structural interpretation.

## Current scope and limitations

The current implementation:

- expects tab-separated input files;
- performs an exact in-memory join;
- summarizes one numeric value column at a time;
- treats state labels as categorical strings and does not interpret their chemistry;
- does not derive protonation states from lambda coordinates;
- does not convert coupled joint states into local per-residue protonation states;
- always keeps trajectory identity in the summary grouping;
- does not pool replicas or pH slots automatically;
- does not estimate autocorrelation times or choose statistically optimal block sizes;
- currently filters `clean`/`matched` state rows before joining, as described above.

These boundaries are deliberate. The tool is intended to provide one auditable state-partitioning layer shared by many structural observables.

## Reproducibility recommendations

For production state-conditioned analysis:

1. Preserve the original frame-wise observation table.
2. Preserve the original state-assignment table.
3. Preserve the joined audit table generated here.
4. Keep absolute MD `step` in both input tables whenever possible.
5. Use an explicit `--join-on` key appropriate to the state model.
6. Use an explicit `--group-by` rather than relying on automatic grouping for final analyses.
7. Inspect state frame counts before interpreting state-conditioned means.
8. Test block-size sensitivity before using `block_sem` as a reported uncertainty.

## Command-line reference

Run:

```bash
python3 partition_by_state.py --help
```

Key options:

| Option | Purpose |
|---|---|
| `--observations FILE` | Frame-wise observation TSV; required |
| `--states FILE` | Frame/state assignment TSV; required |
| `--join-on COLS` | Comma-separated exact join key; default `trajectory,trajectory_frame` |
| `--value-column COL` | Numeric observation column; default `value` |
| `--state-column COL` | Categorical state column; default `state_label` |
| `--clean-column COL` | Clean-state flag column; default `clean` |
| `--matched-column COL` | Matched-state flag column; default `matched` |
| `--include-unclean` | Do not remove unclean state rows before joining |
| `--include-unmatched` | Do not remove unmatched state rows before joining |
| `--group-by COLS` | Comma-separated summary descriptors |
| `--blocks N` | Number of contiguous temporal blocks; default `5` |
| `--joined-output FILE` | Joined frame-wise audit TSV |
| `-o, --output FILE` | State-conditioned summary TSV |

