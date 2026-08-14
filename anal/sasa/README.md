# SASA and fSASA Analysis Tools

This directory contains three complementary tools for solvent-accessibility
analysis:

- `calc_sasa.py` — state-agnostic frame-wise raw SASA calculation from protein
  or model-system trajectories.
- `calc_refsasa.py` — protonation-state-specific model-compound reference SASA
  calculation used to define the fSASA denominator.
- `calc_fsasa.py` — table-based normalization of raw protein SASA to
  protonation-appropriate reference SASA.

The tools are intentionally separated. Coordinate analysis, reference
generation, protonation-state lookup, normalization, and state-conditioned
statistics are different operations and should remain independently auditable.

For state-conditioned statistics, the canonical output from `calc_sasa.py` or
`calc_fsasa.py` is passed to the separate generic tool
`partition_by_state.py`.

---

## 1. Conceptual workflow

For a selected chemical group, the analysis is:

```text
protein trajectories
      |
      v
calc_sasa.py
      |
      |  frame-wise raw LCPO SASA
      v
protein_sasa.tsv
      |
      |                         model-compound trajectories + lambda files
      |                                      |
      |                                      v
      |                               calc_refsasa.py
      |                                      |
      |                                      | chemistry/protonation-specific
      |                                      | reference SASA
      |                                      v
      |                                  refsasa.tsv
      |                                      |
      +-------------------+------------------+
                          |
             authoritative state table
             with local lambda columns
                          |
                          v
                    calc_fsasa.py
                          |
                          | frame-wise fSASA
                          v
                     fsasa.tsv
                          |
                          v
                partition_by_state.py
                          |
                          v
             state-conditioned statistics
```

The frame-wise normalization is

\[
fSASA_i =
\frac{SASA_{\mathrm{protein},i}}
     {\langle SASA_{\mathrm{reference}}\rangle_{\mathrm{chemistry,protonation}}}
\]

where the denominator is selected according to the **local protonation state**
of the site in that frame.

The reference key is:

```text
forcefield + chemistry + site_definition + protonation
```

For example:

```text
ff19SB + ASP + carboxylate_oxygens + H
ff19SB + ASP + carboxylate_oxygens + deprot
ff19SB + GLU + carboxylate_oxygens + H
ff19SB + GLU + carboxylate_oxygens + deprot
```

Biological identities such as `D52` and `E35` are deliberately **not** part of
the reference key.

---

# 2. `calc_sasa.py`

## Purpose

`calc_sasa.py` is the state-agnostic coordinate-analysis layer.

It runs CPPTRAJ `surf` for one or more named atom selections and writes a
canonical long-format frame-wise TSV. It does **not** read lambda files,
assign protonation states, compute fSASA, or calculate state-conditioned
statistics.

Typical uses include:

- catalytic carboxylate oxygen SASA,
- arbitrary residue or atom-group SASA,
- protein SASA observables that will later be partitioned by protonation state,
- the raw numerator for fSASA.

## SASA definition

The tool uses CPPTRAJ `surf`, i.e. the LCPO approximation, with a selected atom
group evaluated in the context of a specified solute:

```text
surf <dataset> "<selection>" solutemask "<solutemask>"
```

The default settings used by these tools are:

```text
offset = 1.4 Å
nbrcut = 2.5 Å
```

For the HEWL catalytic carboxylates, the selections are:

```text
E35_COO  :36@OE1,OE2
D52_COO  :53@OD1,OD2
```

for the capped ff14SB/ff19SB HEWL topology used in this project.

## Important: partial LCPO contributions can be negative

When `surf` is applied to a small selected subset in the context of a larger
`solutemask`, the reported value is the selected atoms' **LCPO contribution**
to the solute SASA. It is not guaranteed to be a nonnegative geometrical area
for every individual frame.

Therefore:

- negative finite values must be preserved;
- they must not be clamped to zero;
- they must not be converted to absolute values;
- they must not be discarded simply because they are negative.

The production `calc_sasa.py` must reject non-finite values but must **not**
reject or clamp finite negative values.

This distinction matters especially for buried carboxylate oxygen selections.
Ensemble means can remain positive and chemically useful even when some
instantaneous partial LCPO contributions are negative.

For manuscript language, a precise description is:

> LCPO surface-area contribution of the selected carboxylate oxygen atoms in
> the context of the protein solute.

## Command-line interface

```text
python3 calc_sasa.py \
    -p TOPOLOGY \
    -t TRAJECTORY_OR_GLOB \
    --site NAME MASK \
    --solutemask MASK \
    [--image-anchor MASK] \
    [--offset 1.4] \
    [--nbrcut 2.5] \
    [--start 1] \
    [--stop last] \
    [--stride 1] \
    [--first-frame-step STEP] \
    [--frame-step-interval N] \
    [--jobs N] \
    [--workdir DIR] \
    -o OUTPUT.tsv
```

`-t/--trajectory` and `--site` may be repeated.

`--jobs`, `--fork`, and `-j` are aliases for trajectory-level concurrency.
Each worker launches an independent CPPTRAJ process. There is no frame-level
parallelism.

## Example: HEWL ff19SB carboxylate oxygen SASA

```bash
python3 /home/wayyne/cphmd/ez_cphmd/anal/sasa/calc_sasa.py \
    -p ../../pre/*.prmtop \
    -t 'arex.ph*.nc' \
    --site E35_COO ':36@OE1,OE2' \
    --site D52_COO ':53@OD1,OD2' \
    --solutemask '^1' \
    --image-anchor '^1' \
    --first-frame-step 510000 \
    --frame-step-interval 10000 \
    --jobs 9 \
    --workdir sasa_work_allph \
    -o test_COO_sasa_allph.tsv
```

For the validated 18-pH HEWL ff19SB dataset this produced:

```text
18 trajectories
5150 frames per trajectory
2 observables per frame

18 x 5150 x 2 = 185400 data rows
```

## Canonical output schema

```text
trajectory
pH
analysis_row
trajectory_frame
step
observable
metric
selection
solutemask
value
unit
```

For this tool:

```text
metric = sasa_lcpo
unit   = A2
```

The first five columns are the canonical frame provenance used throughout the
structural-analysis toolkit.

### Frame fields

- `trajectory` — trajectory stem, e.g. `arex.ph6_5`
- `pH` — parsed from the trajectory name when available
- `analysis_row` — 1-based row number in the analyzed trajectory subset
- `trajectory_frame` — 1-based frame number in the original stored trajectory
- `step` — absolute MD step when frame-step mapping was supplied

## Work directory and reproducibility

The work directory stores per-trajectory CPPTRAJ input, log, and intermediate
data files. A settings file is also written next to the combined output.

Parallel job completion order does not affect final TSV ordering.

If any trajectory-level CPPTRAJ job fails, the combined final TSV should not be
treated as a valid completed analysis.

---

# 3. `calc_refsasa.py`

## Purpose

`calc_refsasa.py` generates the **denominator** used for fSASA normalization.

It is run on a model compound or model peptide containing the chemical group
of interest. For this project, ASP and GLU model tripeptides are analyzed
separately.

The script:

1. calculates frame-wise partial LCPO SASA for the reference selection;
2. maps every coordinate frame to an absolute MD step;
3. exact-joins that step to the corresponding lambda record;
4. classifies the local protonation state;
5. excludes intermediate/mixed lambda frames from each clean reference
   ensemble;
6. computes a protonated (`H`) and deprotonated (`deprot`) ensemble mean;
7. updates a persistent chemistry-keyed reference table.

It does **not** compute protein fSASA.

## Reference identity

References are keyed by chemistry rather than biological residue identity:

```text
forcefield
chemistry
site_definition
protonation
```

Examples:

```text
ff19SB  ASP  carboxylate_oxygens  H
ff19SB  ASP  carboxylate_oxygens  deprot
ff19SB  GLU  carboxylate_oxygens  H
ff19SB  GLU  carboxylate_oxygens  deprot
```

This makes the reference table reusable for any ASP or GLU site analyzed with
the same force field and site definition.

## Local protonation classification

The default clean-state thresholds are:

```text
lambda <= 0.2  -> H
lambda >= 0.8  -> deprot
0.2 < lambda < 0.8 -> mixed
```

These can be changed with `--low` and `--high`, but the thresholds used for
reference generation and subsequent `calc_fsasa.py` analysis must remain
consistent.

## Exact frame/lambda matching

Coordinate frames are mapped to MD steps as:

```text
step =
    first_frame_step
    + (trajectory_frame - 1) * frame_step_interval
```

Matching to the lambda file is exact.

Behavior is:

- frame outside the lambda-file step range:
  retained in the audit as `outside_lambda_window`;
- exact lambda record present:
  assigned normally;
- step lies inside the lambda-file range but has no exact record:
  fatal error;
- nearest-neighbor lambda matching is never allowed.

## Block-SEM diagnostic

Reference statistics include a block SEM.

Blocks are formed on the **complete original coordinate timeline for each
trajectory first**. Protonation-state frames are selected only after block
boundaries have been established.

Thus the block definition is independent of the protonation-state residence
pattern.

The block SEM is a diagnostic of temporal variability; the denominator itself
is the full clean-state ensemble mean.

## Negative reference-frame contributions

As in `calc_sasa.py`, finite negative subset LCPO contributions are valid
possible outputs of `surf` and are preserved.

The final **ensemble-mean reference denominator**, however, must be finite and
strictly positive. A zero or negative reference mean cannot be used as an
fSASA denominator and is treated as an error.

## Command-line interface

```text
python3 calc_refsasa.py \
    --forcefield FORCEFIELD \
    --chemistry ASP_OR_GLU \
    --site-definition SITE_DEFINITION \
    --selection CPPTRAJ_MASK \
    -p TOPOLOGY \
    -t TRAJECTORY_OR_GLOB \
    -l LAMBDA_OR_GLOB \
    --solutemask MASK \
    [--image-anchor MASK] \
    --first-frame-step STEP \
    --frame-step-interval N \
    [--start 1] \
    [--stop last] \
    [--stride 1] \
    [--lambda-step-col 0] \
    [--lambda-value-col 1] \
    [--low 0.2] \
    [--high 0.8] \
    [--blocks 10] \
    [--jobs N] \
    [--workdir DIR] \
    -o refsasa.tsv
```

Trajectory and lambda files are paired by matching basename stems.

## Example: ff19SB ASP reference

```bash
python3 /home/wayyne/cphmd/ez_cphmd/anal/sasa/calc_refsasa.py \
    --forcefield ff19SB \
    --chemistry ASP \
    --site-definition carboxylate_oxygens \
    --selection ':3@OD1,OD2' \
    -p ../../pre/*.prmtop \
    -t 'arex.ph*.nc' \
    -l 'arex.ph*.lambda' \
    --solutemask '^1-5' \
    --first-frame-step 505000 \
    --frame-step-interval 5000 \
    --jobs 6 \
    --workdir model_refsasa_work \
    -o /home/wayyne/cphmd/ez_cphmd/proj/fstpme/tripep/refsasa.tsv \
    --force
```

Validated ff19SB ASP references:

```text
H       90.57428367 A2
deprot  91.88665728 A2
```

The validated ASP accounting was:

```text
coordinate frames      12600
exact matched          12000
outside lambda range     600
mixed/intermediate       384
clean H/deprot         11616

H frames                6117
deprot frames           5499
```

## Example: ff19SB GLU reference

Run from the GLU model-system directory and point `-o` to the **same absolute
reference table**:

```bash
python3 /home/wayyne/cphmd/ez_cphmd/anal/sasa/calc_refsasa.py \
    --forcefield ff19SB \
    --chemistry GLU \
    --site-definition carboxylate_oxygens \
    --selection ':3@OE1,OE2' \
    -p ../../pre/*.prmtop \
    -t 'arex.ph*.nc' \
    -l 'arex.ph*.lambda' \
    --solutemask '^1-5' \
    --first-frame-step 505000 \
    --frame-step-interval 5000 \
    --jobs 6 \
    --workdir model_refsasa_work \
    -o /home/wayyne/cphmd/ez_cphmd/proj/fstpme/tripep/refsasa.tsv \
    --force
```

Validated ff19SB GLU references:

```text
H       90.78260805 A2
deprot  92.45699580 A2
```

The validated GLU accounting was:

```text
coordinate frames      12600
exact matched          12000
outside lambda range     600
mixed/intermediate       453
clean H/deprot         11547

H frames                6423
deprot frames           5124
```

## Persistent reference-table behavior

The same `refsasa.tsv` can be updated by multiple model-chemistry runs.

A new calculation replaces rows with the same key:

```text
forcefield + chemistry + site_definition + protonation
```

and leaves unrelated rows intact.

After the validated ff19SB ASP and GLU calculations, the table contains four
rows:

```text
ff19SB  ASP  carboxylate_oxygens  H
ff19SB  ASP  carboxylate_oxygens  deprot
ff19SB  GLU  carboxylate_oxygens  H
ff19SB  GLU  carboxylate_oxygens  deprot
```

## Reference-table schema

The current table contains 26 columns:

```text
forcefield
chemistry
site_definition
protonation
reference_sasa_A2
sd_A2
block_sem_A2
nframes
nblocks
median_A2
q10_A2
q90_A2
min_A2
max_A2
nnegative
negative_fraction
selection
solutemask
method
offset_A
nbrcut_A
low_cutoff
high_cutoff
topology_file
n_trajectories
blocks_per_trajectory
```

The primary denominator field is:

```text
reference_sasa_A2
```

The additional columns provide provenance and diagnostic statistics.

## Audit output

Each chemistry/reference calculation writes:

```text
<workdir>/<forcefield>/<chemistry>/<site_definition>/
    reference_frame_assignments.tsv
    refsasa.settings.txt
```

The frame audit includes the coordinate/frame identifiers, exact lambda
assignment, protonation classification, clean/matched status, and raw
reference SASA contribution.

---

# 4. `calc_fsasa.py`

## Purpose

`calc_fsasa.py` combines:

- raw frame-wise protein SASA from `calc_sasa.py`;
- local protonation coordinates from an authoritative state table;
- chemistry/protonation-specific reference values from `calc_refsasa.py`.

It then writes frame-wise fSASA in the same canonical observable format used by
the rest of the analysis toolkit.

`calc_fsasa.py` does **not** read coordinates and does **not** run CPPTRAJ.
Changing a reference table or state assignment therefore does not require
rerunning the expensive coordinate-level SASA calculation.

## Site mapping

Each protein observable must be mapped to:

```text
OBSERVABLE
CHEMISTRY
SITE_DEFINITION
LOCAL_LAMBDA_COLUMN
```

For the HEWL dyad state table:

```text
E35_COO -> GLU -> carboxylate_oxygens -> primary_lambda
D52_COO -> ASP -> carboxylate_oxygens -> coupled_lambda
```

This mapping is supplied with repeated `--site` arguments:

```text
--site E35_COO GLU carboxylate_oxygens primary_lambda
--site D52_COO ASP carboxylate_oxygens coupled_lambda
```

`primary_label` and `coupled_label` are **site identity labels** such as `E35`
and `D52`; they are not protonation-state columns and must not be passed here.

## Local protonation assignment

For matched state rows, `calc_fsasa.py` classifies the specified local lambda
coordinate using:

```text
lambda <= 0.2  -> H
lambda >= 0.8  -> deprot
otherwise      -> local mixed
```

The reference denominator is then looked up using:

```text
forcefield + chemistry + site_definition + local_protonation
```

## Exact state-table join

The default join key is:

```text
trajectory + trajectory_frame
```

The script then independently validates shared frame metadata such as
`analysis_row`, `trajectory_frame`, `step`, and pH when populated on both
sides.

There is no nearest-frame or nearest-step matching.

## Treatment of matched, mixed, and globally unclean frames

The distinction between **local** and **global/coupled** cleanliness is
important.

### Unmatched/outside state window

The row is retained in the canonical output, but `value` is blank.

Audit status:

```text
outside_state_window
```

### Matched but local lambda is intermediate

The row is retained, but `value` is blank.

Audit status:

```text
local_mixed
```

### Local lambda is clean

fSASA is calculated from the appropriate protonation-specific denominator.

This is true even if the **global coupled state** is unclean because the other
site is intermediate.

For example, if E35 has a clean local lambda but D52 is intermediate, the E35
fSASA value is still chemically well defined and is retained.

Later, `partition_by_state.py` excludes globally unclean coupled-state frames
from coupled-state summaries.

This separation maximizes preservation of valid local information while
keeping the state-conditioned analysis rigorous.

### Consistency check

If the authoritative state table says a frame is globally `clean` while the
mapped local lambda is intermediate, the script stops with an error. Such a
case indicates inconsistent state definitions or lambda thresholds.

## Negative values

The reference denominators are positive ensemble means.

Therefore normalization does not change the sign of a raw value:

```text
negative raw partial LCPO contribution
    -> negative frame-wise fSASA
```

Such values are preserved.

Do not clamp, take absolute values, or remove them at the fSASA stage.

## Command-line interface

```text
python3 calc_fsasa.py \
    --sasa RAW_SASA.tsv \
    --states STATE_ASSIGNMENTS.tsv \
    --references refsasa.tsv \
    --forcefield FORCEFIELD \
    --site OBSERVABLE CHEMISTRY SITE_DEFINITION LAMBDA_COLUMN \
    [--site ...] \
    [--join-on trajectory,trajectory_frame] \
    [--low 0.2] \
    [--high 0.8] \
    [--clean-column clean] \
    [--matched-column matched] \
    [--sasa-metric sasa_lcpo] \
    [--audit-output FILE] \
    -o FSASA.tsv
```

## Example: HEWL ff19SB E35/D52 fSASA

```bash
python3 /home/wayyne/cphmd/ez_cphmd/anal/sasa/calc_fsasa.py \
    --sasa /home/wayyne/cphmd/amber24/ff19/ti/harris_2022/prot/hewl/val-take2/rex-cont/test_COO_sasa_allph.tsv \
    --states /home/wayyne/cphmd/amber24/ff19/ti/harris_2022/prot/hewl/val-take2/rex-cont/hewl_ff19_analysis/dyad_frame_assignments.tsv \
    --references /home/wayyne/cphmd/ez_cphmd/proj/fstpme/tripep/refsasa.tsv \
    --forcefield ff19SB \
    --site E35_COO GLU carboxylate_oxygens primary_lambda \
    --site D52_COO ASP carboxylate_oxygens coupled_lambda \
    --low 0.2 \
    --high 0.8 \
    -o test_COO_fsasa_allph.tsv
```

Validated accounting for the 18-pH ff19SB HEWL calculation:

```text
Input SASA rows                    : 185400
Exact state joins                  : 185400
Outside state window               : 5400
Local mixed rows                   : 6995
fSASA values computed              : 173005
Global-unclean but local-clean rows: 6605
Negative raw SASA rows             : 5735
Negative fSASA rows                : 5389
```

The accounting identities are:

```text
185400 - 5400 - 6995 = 173005
```

and

```text
173005 - 6605 = 166400
```

Thus every globally clean coupled-state observation has a valid fSASA value.

The globally mixed accounting also closes:

```text
6800 globally mixed dyad frames x 2 sites = 13600 site observations

6995 locally mixed
6605 locally clean while the other site is mixed
----
13600
```

## Canonical output schema

The primary fSASA output has the same canonical observable structure:

```text
trajectory
pH
analysis_row
trajectory_frame
step
observable
metric
selection
solutemask
value
unit
```

For this tool:

```text
metric = fsasa
unit   = 1
```

Rows with unavailable fSASA remain in the file with a blank `value`.

Keeping these rows is important because downstream timeline-based block
construction should operate on the original trajectory timeline rather than a
state-filtered or value-filtered compressed timeline.

## Audit output

By default:

```text
<output_stem>.audit.tsv
```

For example:

```text
test_COO_fsasa_allph.audit.tsv
```

The audit records:

```text
trajectory
pH
analysis_row
trajectory_frame
step
observable
raw_metric
selection
solutemask
raw_sasa_A2
chemistry
site_definition
lambda_column
local_lambda
local_protonation
state_clean
state_matched
state_code
state_label
reference_forcefield
reference_sasa_A2
reference_method
fsasa
fsasa_status
```

A settings file is also written:

```text
<output>.settings.txt
```

It records the input files, join key, lambda thresholds, site mappings, and
actual denominator values used.

---

# 5. Complete fSASA workflow

This section shows the three SASA tools together using the validated ff19SB
HEWL analysis as an example.

## Step 1 — Calculate the protein numerator

From the HEWL production trajectory directory:

```bash
python3 /home/wayyne/cphmd/ez_cphmd/anal/sasa/calc_sasa.py \
    -p ../../pre/*.prmtop \
    -t 'arex.ph*.nc' \
    --site E35_COO ':36@OE1,OE2' \
    --site D52_COO ':53@OD1,OD2' \
    --solutemask '^1' \
    --image-anchor '^1' \
    --first-frame-step 510000 \
    --frame-step-interval 10000 \
    --jobs 9 \
    --workdir sasa_work_allph \
    -o test_COO_sasa_allph.tsv
```

Expected validated size:

```text
185400 data rows
= 18 pH trajectories x 5150 frames x 2 sites
```

## Step 2 — Generate the ASP denominator

From the ASP model-peptide trajectory directory:

```bash
python3 /home/wayyne/cphmd/ez_cphmd/anal/sasa/calc_refsasa.py \
    --forcefield ff19SB \
    --chemistry ASP \
    --site-definition carboxylate_oxygens \
    --selection ':3@OD1,OD2' \
    -p ../../pre/*.prmtop \
    -t 'arex.ph*.nc' \
    -l 'arex.ph*.lambda' \
    --solutemask '^1-5' \
    --first-frame-step 505000 \
    --frame-step-interval 5000 \
    --jobs 6 \
    --workdir model_refsasa_work \
    -o /home/wayyne/cphmd/ez_cphmd/proj/fstpme/tripep/refsasa.tsv \
    --force
```

## Step 3 — Generate the GLU denominator

From the GLU model-peptide trajectory directory:

```bash
python3 /home/wayyne/cphmd/ez_cphmd/anal/sasa/calc_refsasa.py \
    --forcefield ff19SB \
    --chemistry GLU \
    --site-definition carboxylate_oxygens \
    --selection ':3@OE1,OE2' \
    -p ../../pre/*.prmtop \
    -t 'arex.ph*.nc' \
    -l 'arex.ph*.lambda' \
    --solutemask '^1-5' \
    --first-frame-step 505000 \
    --frame-step-interval 5000 \
    --jobs 6 \
    --workdir model_refsasa_work \
    -o /home/wayyne/cphmd/ez_cphmd/proj/fstpme/tripep/refsasa.tsv \
    --force
```

After Steps 2 and 3, inspect:

```bash
column -t -s $'\t' /home/wayyne/cphmd/ez_cphmd/proj/fstpme/tripep/refsasa.tsv
```

The validated ff19SB denominators are:

```text
ASP H       90.57428367 A2
ASP deprot  91.88665728 A2
GLU H       90.78260805 A2
GLU deprot  92.45699580 A2
```

## Step 4 — Normalize protein SASA to frame-wise fSASA

```bash
python3 /home/wayyne/cphmd/ez_cphmd/anal/sasa/calc_fsasa.py \
    --sasa /home/wayyne/cphmd/amber24/ff19/ti/harris_2022/prot/hewl/val-take2/rex-cont/test_COO_sasa_allph.tsv \
    --states /home/wayyne/cphmd/amber24/ff19/ti/harris_2022/prot/hewl/val-take2/rex-cont/hewl_ff19_analysis/dyad_frame_assignments.tsv \
    --references /home/wayyne/cphmd/ez_cphmd/proj/fstpme/tripep/refsasa.tsv \
    --forcefield ff19SB \
    --site E35_COO GLU carboxylate_oxygens primary_lambda \
    --site D52_COO ASP carboxylate_oxygens coupled_lambda \
    --low 0.2 \
    --high 0.8 \
    -o test_COO_fsasa_allph.tsv
```

## Step 5 — Partition fSASA by coupled protonation state

`partition_by_state.py` is not part of the SASA calculation itself. It is the
generic downstream statistics layer.

```bash
python3 /home/wayyne/cphmd/ez_cphmd/anal/partition/partition_by_state.py \
    --observations test_COO_fsasa_allph.tsv \
    --states /home/wayyne/cphmd/amber24/ff19/ti/harris_2022/prot/hewl/val-take2/rex-cont/hewl_ff19_analysis/dyad_frame_assignments.tsv \
    --join-on trajectory,trajectory_frame \
    --group-by pH,observable,metric,selection,solutemask,unit \
    --blocks 5 \
    --joined-output test_COO_fsasa_allph_with_states.tsv \
    -o test_COO_fsasa_allph_by_state.tsv
```

Validated accounting:

```text
Joined rows              : 185400
Rows used in summaries   : 166400
Excluded unmatched rows  : 5400
Excluded unclean rows    : 13600
```

This is exactly the same globally clean observation count obtained when the
raw SASA table was partitioned.

## Step 6 — Inspect state-conditioned results

```bash
column -t -s $'\t' test_COO_fsasa_allph_by_state.tsv | less -S
```

The summary table contains:

```text
trajectory
pH
observable
metric
selection
solutemask
unit
state
nframes
mean
sd
block_sem
nblocks
median
q10
q90
min
max
```

---

# 6. Why the tools are separate

The separation is deliberate and should be maintained.

## `calc_sasa.py`

Knows about:

- coordinates,
- topology,
- CPPTRAJ masks,
- LCPO settings,
- trajectory/frame mapping.

Does not know about:

- lambda,
- protonation states,
- coupled states,
- reference normalization.

## `calc_refsasa.py`

Knows about:

- model-system coordinates,
- one local titration coordinate,
- exact lambda-step matching,
- chemistry-level reference identity,
- state-specific reference statistics.

Does not know about:

- HEWL E35 or D52 identities,
- protein coupled-state labels,
- protein fSASA.

## `calc_fsasa.py`

Knows about:

- canonical raw SASA tables,
- authoritative local lambda columns,
- chemistry/site mappings,
- reference-table lookup,
- frame-wise normalization.

Does not know about:

- coordinates,
- CPPTRAJ,
- how the state table itself was generated,
- state-conditioned averaging.

## `partition_by_state.py`

Knows about:

- authoritative state assignments,
- exact frame joins,
- globally clean/matched filtering,
- timeline-first block statistics.

Does not know whether an observable is:

- distance,
- SASA,
- fSASA,
- hydrogen bonding,
- hydration,
- or another structural quantity.

This architecture allows each layer to be tested and reused independently.

---

# 7. Data-contract requirements

## Raw SASA to fSASA

`calc_fsasa.py` expects the canonical raw SASA fields:

```text
trajectory
analysis_row
trajectory_frame
step
observable
metric
selection
solutemask
value
unit
```

`pH` or `ph` may also be present and are validated/canonicalized.

The expected raw metric is:

```text
sasa_lcpo
```

unless changed with `--sasa-metric`.

## State table

The default exact join is:

```text
trajectory + trajectory_frame
```

The current HEWL state table additionally contains fields including:

```text
analysis_row
step
pH
primary_label
primary_lambda
coupled_label
coupled_lambda
state_code
state_label
clean
matched
```

For the current dyad analysis:

```text
primary_label   = E35
primary_lambda  = E35 local lambda
coupled_label   = D52
coupled_lambda  = D52 local lambda
```

Do not infer protonation from the label columns.

## Reference table

At minimum, `calc_fsasa.py` requires:

```text
forcefield
chemistry
site_definition
protonation
reference_sasa_A2
```

The production `calc_refsasa.py` writes substantially more provenance.

---

# 8. Validation checklist

Before accepting a new fSASA analysis, check all of the following.

### Raw SASA

- every expected trajectory completed;
- expected frame count is present for every observable;
- raw output ordering is deterministic;
- finite negative partial-LCPO values are preserved;
- `offset`, `nbrcut`, selections, solutemask, topology, and frame-step mapping
  are recorded.

### Reference SASA

- every coordinate frame is accounted for as matched or outside the lambda
  range;
- every matched frame is classified as `H`, `deprot`, or `mixed`;
- `H + deprot + mixed = matched`;
- the reference mean is finite and > 0 for both clean protonation states;
- model/reference selection matches the chemical definition of the protein
  numerator;
- ASP and GLU rows are written to the intended shared reference table;
- reference force field matches the protein force field.

### fSASA

- every raw SASA observation exact-joins to one authoritative state row;
- no nearest-neighbor frame matching occurs;
- locally mixed/unmatched rows remain in the canonical timeline with blank
  fSASA;
- locally clean rows use the correct chemistry/protonation denominator;
- negative raw values remain negative after division by positive references;
- globally clean coupled-state rows all have valid fSASA values.

### State partition

- joined row count equals the canonical fSASA row count;
- blocks are formed from the complete original timeline before state
  selection;
- globally unmatched and unclean counts agree with the authoritative state
  assignment;
- rare one-frame states retain blank SD/SEM rather than fabricated
  uncertainty.

---

# 9. Current ff19SB HEWL validation summary

The validated production chain for the catalytic E35/D52 carboxylate oxygen
analysis is:

```text
raw protein SASA observations       185400
fSASA exact state joins             185400
outside state-window observations     5400
locally mixed observations            6995
computed frame-wise fSASA           173005
global-unclean/local-clean            6605

globally clean observations         166400
```

The downstream state partition gives:

```text
joined                            185400
used in summaries                166400
excluded unmatched                5400
excluded globally unclean        13600
```

The validated ff19SB reference values are:

```text
ASP carboxylate_oxygens H       90.57428367 A2
ASP carboxylate_oxygens deprot  91.88665728 A2
GLU carboxylate_oxygens H       90.78260805 A2
GLU carboxylate_oxygens deprot  92.45699580 A2
```

---

# 10. Interpretation notes

The quantity produced here is named `fSASA` for continuity with the
solvent-exposure framework, but its precise computational definition should be
reported in Methods.

Because both numerator and denominator use CPPTRAJ `surf`, the numerator is a
partial LCPO contribution of the selected atoms to the SASA of the specified
solute context. Individual frame values therefore need not be bounded between
0 and 1 and can be negative.

Recommended interpretation:

- use ensemble/state-conditioned means and distributions;
- define the selected atoms and reference model explicitly;
- state that numerator and denominator use the same LCPO/site definition;
- avoid interpreting an individual negative value as a negative physical
  surface area;
- do not silently truncate the distribution.

For the current HEWL analysis, state-conditioned mean fSASA values are positive
and well behaved.

---

# 11. Suggested directory layout

A simple layout is:

```text
sasa/
├── calc_sasa.py
├── calc_refsasa.py
├── calc_fsasa.py
└── README.md
```

Reference tables and project-specific run scripts should normally live with
the project data rather than inside the analysis-code directory, for example:

```text
proj/fstpme/
├── tripep/
│   ├── get_refsasa.sh
│   └── refsasa.tsv
└── hewl/
    └── fsasa_analysis.sh
```

This keeps the analysis tools reusable while preserving project-specific
commands and generated data with the project.

---

# 12. Quick reference

Raw protein SASA:

```bash
python3 calc_sasa.py \
    -p topology.prmtop \
    -t 'arex.ph*.nc' \
    --site E35_COO ':36@OE1,OE2' \
    --site D52_COO ':53@OD1,OD2' \
    --solutemask '^1' \
    --first-frame-step 510000 \
    --frame-step-interval 10000 \
    --jobs 9 \
    -o protein_sasa.tsv
```

Reference denominator:

```bash
python3 calc_refsasa.py \
    --forcefield ff19SB \
    --chemistry ASP \
    --site-definition carboxylate_oxygens \
    --selection ':3@OD1,OD2' \
    -p model.prmtop \
    -t 'arex.ph*.nc' \
    -l 'arex.ph*.lambda' \
    --solutemask '^1-5' \
    --first-frame-step 505000 \
    --frame-step-interval 5000 \
    --jobs 6 \
    -o refsasa.tsv
```

Frame-wise fSASA:

```bash
python3 calc_fsasa.py \
    --sasa protein_sasa.tsv \
    --states frame_assignments.tsv \
    --references refsasa.tsv \
    --forcefield ff19SB \
    --site E35_COO GLU carboxylate_oxygens primary_lambda \
    --site D52_COO ASP carboxylate_oxygens coupled_lambda \
    -o protein_fsasa.tsv
```

State-conditioned fSASA:

```bash
python3 ../partition/partition_by_state.py \
    --observations protein_fsasa.tsv \
    --states frame_assignments.tsv \
    --join-on trajectory,trajectory_frame \
    --group-by pH,observable,metric,selection,solutemask,unit \
    --blocks 5 \
    --joined-output protein_fsasa_with_states.tsv \
    -o protein_fsasa_by_state.tsv
```

