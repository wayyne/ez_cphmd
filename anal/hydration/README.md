# calc_hydration.py

General first-shell hydration observable generator for Amber trajectories.

`calc_hydration.py` calculates the number of water molecules in the first hydration shell of one or more user-defined solute sites using CPPTRAJ `watershell`. The output is a long-format TSV designed to integrate directly with `partition_by_state.py`.

The analyzer is deliberately **state-agnostic**: it measures hydration geometry only. Protonation-state or microscopic-state conditioning is performed downstream by joining the output to an authoritative frame-state table.

The default first-shell definition is a **unique-water union count within 3.4 Å of any atom in the supplied site mask**. With a water-oxygen mask such as `:WAT@O`, a water that lies within 3.4 Å of multiple atoms in the same site is counted once.

## Requirements

- Python 3
- AmberTools / CPPTRAJ available on `$PATH`
- Amber topology (`.prmtop`)
- One or more Amber trajectories, typically NetCDF (`.nc`)

## Basic usage

```bash
python3 calc_hydration.py \
    -p system.prmtop \
    -t "arex.ph*.nc" \
    --site E35 ':36@OE1,OE2' \
    --site D52 ':53@OD1,OD2' \
    --water-mask ':WAT@O' \
    --cutoff 3.4 \
    --image-anchor '^1' \
    --first-frame-step 510000 \
    --frame-step-interval 10000 \
    --jobs 4 \
    --workdir hydration_work \
    -o hydration.tsv
```

## HEWL ff19SB example

```bash
ff19_hewl=/home/wayyne/cphmd/amber24/ff19/ti/harris_2022/prot/hewl/val-take2/rex-cont

python3 calc_hydration.py \
  -p ${ff19_hewl}/../../pre/*.prmtop \
  -t "${ff19_hewl}/arex.ph*.nc" \
  --site H15 ':16@ND1,NE2' \
  --site E35 ':36@OE1,OE2' \
  --site D52 ':53@OD1,OD2' \
  --water-mask ':WAT@O' \
  --cutoff 3.4 \
  --image-anchor '^1' \
  --first-frame-step 510000 \
  --frame-step-interval 10000 \
  --jobs 9 \
  --workdir ff19_hydration_work_v1_0 \
  -o ff19_hewl_hydration_v1_0.tsv
```

## Hydration-site semantics

Each `--site` consists of:

```text
LABEL MASK
```

For example:

```bash
--site E35 ':36@OE1,OE2'
```

defines one hydration observable named `E35` around the union of the two carboxylate oxygens.

The code passes the entire site mask to CPPTRAJ `watershell`. Therefore, the first-shell count is a **unique solvent-molecule count around the whole site**, not a sum of separate per-atom hydration numbers.

Example:

```text
E35 site = OE1 OR OE2
```

If one water oxygen lies within 3.4 Å of both OE1 and OE2, that water contributes:

```text
1 water
```

not:

```text
2 waters
```

This union-count behavior is important for carboxylate hydration.

## Historical 3.4 Å definition

The default first-shell cutoff is:

```text
3.4 Å
```

The script was designed to reproduce the first-shell water-count convention used in recent Shen-lab structural analyses.

Examples used in the HEWL analysis are:

```text
H15 : water O within 3.4 Å of ND1 OR NE2
E35 : water O within 3.4 Å of OE1 OR OE2
D52 : water O within 3.4 Å of OD1 OR OD2
```

The cutoff can be changed explicitly with:

```bash
--cutoff VALUE
```

but changing it creates a different hydration observable and should be documented accordingly.

## Water selection

The default solvent atom mask is:

```text
:WAT@O
```

specified by:

```bash
--water-mask ':WAT@O'
```

This selects water oxygen atoms.

The water mask must select at least one atom in the supplied topology. The script validates both site masks and the water mask before launching trajectory analysis.

## CPPTRAJ `watershell`

For each trajectory and each requested site, the script generates a CPPTRAJ command of the form:

```text
watershell SITE_MASK out FILE lower 3.4 upper 5.0 WATER_MASK
```

The lower-shell count is retained as the structural observable.

CPPTRAJ also calculates an upper-shell dataset because `watershell` requires an upper cutoff, but this program does **not** emit that value.

By default:

```text
first-shell cutoff = 3.4 Å
upper cutoff       = 5.0 Å
```

The upper cutoff can be changed with:

```bash
--upper-cutoff VALUE
```

provided:

```text
upper_cutoff > first_shell_cutoff
```

The upper-shell values remain intermediate-only and are not written to the final TSV.

## Imaging

CPPTRAJ `watershell` uses periodic imaging unless `noimage` is specified. This program does not disable imaging.

An optional global CPPTRAJ autoimage step may also be requested:

```bash
--image-anchor '^1'
```

which adds:

```text
autoimage anchor ^1
```

before the hydration calculations.

For the capped ff14SB/ff19SB HEWL systems, `^1` is the protein molecule used as the imaging anchor.

## Output format

The output is a tab-separated long-format table with one row per trajectory frame and hydration site.

Columns are:

```text
trajectory
pH
analysis_row
trajectory_frame
step
observable
site_mask
water_mask
cutoff_A
metric
value
unit
```

### trajectory

The trajectory identifier is written using the trajectory filename stem.

For example:

```text
arex.ph6_5.nc
```

becomes:

```text
arex.ph6_5
```

This allows exact joining to state-assignment tables that use the same trajectory identifier.

### trajectory_frame

`trajectory_frame` is one-based:

```text
first trajectory frame  -> 1
second trajectory frame -> 2
...
```

### analysis_row

`analysis_row` is a **frame-level identifier**, not an output-row identifier.

For every site generated from trajectory frame 1:

```text
analysis_row = 1
```

For every site generated from trajectory frame 2:

```text
analysis_row = 2
```

and so forth.

This matches the canonical observable-table interface used by `partition_by_state.py`.

### step

The MD step is reconstructed as:

```text
step = first_frame_step
     + (trajectory_frame - 1) * frame_step_interval
```

For the ff19SB HEWL trajectories:

```text
first_frame_step     = 510000
frame_step_interval  = 10000
```

## Metric

The program currently emits one metric:

```text
first_shell_water_count
```

with unit:

```text
water
```

Each value is the number of unique selected solvent molecules within the first-shell cutoff of any atom in the site mask.

No normalization is applied.

The program does not calculate:

```text
fractional hydration
per-oxygen summed hydration
RDF-derived coordination numbers
residence times
water identities
hydration free energies
```

Those would be distinct observables and should not be inferred from `first_shell_water_count`.

## State conditioning

`calc_hydration.py` does not read CpHMD lambda coordinates and does not assign protonation states.

This is intentional.

Hydration is first calculated as a pure structural observable:

```text
trajectory
    |
    v
calc_hydration.py
    |
    v
state-agnostic hydration TSV
```

The resulting table is then joined to an authoritative state-assignment table:

```text
hydration TSV
    +
state assignments
    |
    v
partition_by_state.py
    |
    v
state-conditioned hydration statistics
```

## HEWL E35/D52 dyad states

For the dyad analysis:

```text
S0       E35- / D52-
S1_E35   E35H / D52-
S1_D52   E35- / D52H
S2       E35H / D52H
```

Therefore:

```text
E35 hydration in S1_E35
```

means hydration of protonated E35 in the `E35H/D52-` microstate.

Likewise:

```text
D52 hydration in S1_D52
```

means hydration of protonated D52 in the `E35-/D52H` microstate.

## H15 interpretation

If H15 hydration is partitioned with the E35/D52 dyad state table, the result means:

```text
H15 hydration conditional on E35/D52 dyad state
```

It does **not** mean H15 hydration conditioned on H15 protonation or tautomer.

For a historical H15 analysis conditioned on HIP/HID/HIE, an H15-specific authoritative state assignment should be used downstream.

## Partitioning example

For ff19SB HEWL:

```bash
partition=/home/wayyne/cphmd/ez_cphmd/anal/partition/partition_by_state.py

ff19_hewl=/home/wayyne/cphmd/amber24/ff19/ti/harris_2022/prot/hewl/val-take2/rex-cont
ff19_states=${ff19_hewl}/hewl_ff19_analysis/dyad_frame_assignments.tsv

python3 ${partition} \
    --observations ff19_hewl_hydration_v1_0.tsv \
    --states "${ff19_states}" \
    --join-on trajectory,trajectory_frame \
    --value-column value \
    --state-column state_label \
    --blocks 5 \
    --joined-output ff19_hewl_hydration_v1_0.by_state.joined.tsv \
    -o ff19_hewl_hydration_v1_0.by_state.tsv
```

## Block statistics

Temporal blocking is not performed inside `calc_hydration.py`.

Blocking is handled by `partition_by_state.py`.

The structural suite forms temporal blocks on the original trajectory timeline before state selection, so the reported block SEM is a diagnostic of temporal variability rather than a claim of fully independent sampling.

## Parallel execution

Trajectory-level parallelism is controlled by:

```bash
--jobs N
```

CPPTRAJ OpenMP threads per worker are controlled separately:

```bash
--cpptraj-threads N
```

Example:

```bash
--jobs 9 --cpptraj-threads 1
```

runs nine trajectory jobs concurrently with one CPPTRAJ thread each.

The script also constrains common BLAS thread pools to one thread per worker to reduce accidental CPU oversubscription.

## pH parsing

By default, the pH is extracted from trajectory names using:

```text
ph(?P<ph>\d+(?:[_\.]\d+)?)
```

Examples:

```text
arex.ph5.nc
arex.ph5_5.nc
arex.ph6.5.nc
```

are interpreted as:

```text
5.0
5.5
6.5
```

A custom expression can be supplied with:

```bash
--ph-regex ...
```

The regular expression must contain a named capture group:

```text
(?P<ph>...)
```

## Work directory

CPPTRAJ inputs, logs, and intermediate `watershell` series are retained under:

```bash
--workdir hydration_work
```

For each trajectory/site combination, an intermediate file similar to:

```text
arex.ph6_5.E35.watershell.dat
```

is retained.

CPPTRAJ input and log files are also preserved for audit and debugging.

## Settings file

Every successful run writes:

```text
OUTPUT.tsv.settings.txt
```

The settings file records:

```text
program
version
cpptraj executable
topology
number of trajectories
number of sites
first-frame MD step
frame-step interval
image anchor
water mask
first-shell cutoff
CPPTRAJ upper cutoff
output-row count
metric name
count semantics
historical-definition label
site labels and masks
```

## Expected output accounting

Because there is one emitted metric per site per frame:

```text
output_rows =
    number_of_trajectory_frames
    × number_of_sites
```

For the ff19SB HEWL run:

```text
18 trajectories
× 5150 frames per trajectory
= 92,700 trajectory frames
```

with three sites:

```text
H15
E35
D52
```

there are:

```text
92,700 × 3 = 278,100 output rows
```

plus the TSV header.

The validated ff19SB run reported:

```text
Trajectories analyzed       : 18
Hydration sites             : 3
Trajectory frames           : 92700
Frame-wise output rows      : 278100
Program version             : 1.0.0
First-shell cutoff          : 3.4 A
Water mask                  : :WAT@O
Metric                      : first_shell_water_count
```

## State-partition accounting

For the same ff19SB run, exact state joining produced:

```text
Joined rows              : 278100
Rows used in summaries   : 249600
Excluded unmatched rows  : 8100
Excluded unclean rows    : 20400
```

These values are consistent with the underlying state table:

```text
2700 unmatched frames × 3 sites = 8100 rows
6800 unclean frames   × 3 sites = 20400 rows
```

which provides an end-to-end bookkeeping check.

## Validation checks

Recommended checks before using results are:

1. Confirm that all site masks select at least one atom.
2. Confirm that the water mask selects the intended solvent oxygen atoms.
3. Confirm that the trajectory count and frame count match expectations.
4. Confirm that:

   ```text
   output_rows = total_frames × number_of_sites
   ```

5. Confirm that `trajectory` identifiers exactly match the state-assignment table.
6. Confirm that all site rows from one trajectory frame share the same `analysis_row`.
7. Run `partition_by_state.py` and verify exact joined/unmatched/unclean accounting.
8. Inspect CPPTRAJ logs if any trajectory fails.
9. Keep the same site mask, water mask, and cutoff when comparing force fields.

## Force-field comparisons

For cross-force-field comparisons, the physical site definition should remain chemically equivalent even if topology residue numbering differs.

For capped ff14SB and ff19SB HEWL:

```text
H15 = topology residue 16
E35 = topology residue 36
D52 = topology residue 53
```

The corresponding c22 masks must be mapped to the actual c22 topology before running.

Do not assume the ff14SB/ff19SB residue indices apply directly to c22.

## Version

Current program version:

```text
1.0.0
```

The source defines:

```python
PROGRAM_VERSION = "1.0.0"
```

## Design philosophy

`calc_hydration.py` follows the same separation used by the other structural-analysis tools:

```text
calc_distance.py
calc_sasa.py
calc_hbond.py
calc_dihedral.py
calc_hydration.py
        |
        v
frame-level structural observables
        |
        + authoritative state assignments
        |
        v
partition_by_state.py
        |
        v
state-conditioned structural statistics
```

This keeps observable generation independent of a particular thermodynamic-state definition and allows the same hydration table to be reused for different mechanistic analyses.

