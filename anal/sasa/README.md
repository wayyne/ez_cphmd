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

For a selected chemical group, the analysis is deliberately split into three
layers:

```text
protein trajectories
      |
      v
calc_sasa.py
      |
      | frame-wise raw LCPO SASA
      v
protein_sasa.tsv
      |
      |                         model-compound trajectories + lambda files
      |                                      |
      |                                      v
      |                               calc_refsasa.py
      |                                      |
      |                                      | chemistry/state-specific
      |                                      | reference SASA
      |                                      v
      |                                  refsasa.tsv
      |                                      |
      +-------------------+------------------+
                          |
              authoritative state table
              with local lambda coordinates
              (+ tautomer x for His)
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
     {\langle SASA_{\mathrm{reference}}\rangle_{\mathrm{chemistry,state}}}
\]

where the denominator is selected according to the **local chemical state** of
the site in that frame.

The reference key remains:

```text
forcefield + chemistry + site_definition + protonation
```

The field name `protonation` is retained for backward compatibility, but its
allowed values depend on the lambda mode.

For ordinary one-coordinate titratable groups:

```text
ff19SB + ASP + carboxylate_oxygens + H
ff19SB + ASP + carboxylate_oxygens + deprot
ff19SB + GLU + carboxylate_oxygens + H
ff19SB + GLU + carboxylate_oxygens + deprot
```

For histidine:

```text
ff19SB + HIS + <site_definition> + HIP
ff19SB + HIS + <site_definition> + HID
ff19SB + HIS + <site_definition> + HIE
```

Biological identities such as `D52`, `E35`, or `H15` are deliberately **not**
part of the reference key. The same chemistry-level reference can therefore be
reused for another biological site only when the force field and atom/site
definition are identical.

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
- histidine side-chain or ring-atom SASA,
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

For the validated HEWL catalytic carboxylate analysis, the selections are:

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

It is run on a model compound or model peptide containing the chemical group of
interest. It supports two state-classification modes:

- `--lambda-mode single` — legacy one-coordinate `H`/`deprot` references for
  ASP, GLU, and other ordinary single-coordinate titratable groups;
- `--lambda-mode his` — histidine-aware `HIP`/`HID`/`HIE` references using the
  Amber protonation and neutral-tautomer coordinates.

The script:

1. calculates frame-wise partial LCPO SASA for the reference selection;
2. maps every coordinate frame to an absolute MD step;
3. exact-joins that step to the corresponding lambda record;
4. classifies the local chemical state;
5. excludes intermediate/mixed frames from clean reference ensembles;
6. computes one ensemble-mean denominator for every required clean state;
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

Examples for single-coordinate groups:

```text
ff19SB  ASP  carboxylate_oxygens  H
ff19SB  ASP  carboxylate_oxygens  deprot
ff19SB  GLU  carboxylate_oxygens  H
ff19SB  GLU  carboxylate_oxygens  deprot
```

Examples for histidine:

```text
ff19SB  HIS  histidine_sidechain  HIP
ff19SB  HIS  histidine_sidechain  HID
ff19SB  HIS  histidine_sidechain  HIE
```

The exact `site_definition` string is project-defined and must describe the
same atom selection used for both model-reference and protein SASA.

## Single-coordinate mode

`--lambda-mode single` is the default and preserves the original behavior.

The default clean-state thresholds are:

```text
lambda <= 0.2       -> H
lambda >= 0.8       -> deprot
0.2 < lambda < 0.8  -> mixed
```

The lambda step and value columns are selected explicitly with:

```text
--lambda-step-col
--lambda-value-col
```

`--chemistry HIS` is intentionally rejected in this mode so neutral HID and HIE
cannot be accidentally pooled into one denominator.

## Histidine mode

Use:

```text
--lambda-mode his
```

for histidine model-compound trajectories.

In this mode the script reads the Amber lambda-file `ires` and `itauto` headers.
For the selected His `ires`:

```text
itauto = 1  -> protonation coordinate lambda
itauto = 2  -> neutral-tautomer coordinate x
```

The His `ires` can be supplied explicitly:

```text
--lambda-resid N
```

or, if omitted, the script requires that exactly one `ires` in the lambda file
contains one `itauto=1` variable and one `itauto=2` variable.

Physical-state assignment is:

```text
lambda <= low                         -> HIP
lambda >= high and x <= low           -> x-low tautomer
lambda >= high and x >= high          -> x-high tautomer
low < lambda < high                   -> mixed
lambda >= high and low < x < high     -> neutral_tautomer_mixed
```

`HIP` is assigned independently of `x`.

The mapping of the two neutral `x` end states is supplied explicitly:

```text
--x-low-tautomer HID|HIE
--x-high-tautomer HID|HIE
```

For the force-field parameterizations used in the HEWL/H15 work:

```text
c22:
    x low  -> HIE
    x high -> HID

ff14SB / ff19SB:
    x low  -> HID
    x high -> HIE
```

The script requires the two endpoints to map to different tautomers.

Frames with intermediate protonation lambda or neutral-intermediate `x` remain
in the audit table but do not contribute to any `HIP`, `HID`, or `HIE`
denominator.

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

In His mode, reset/restarted lambda step counters are first unwrapped into a
strictly increasing step series before exact matching. The unwrapping logic is
the same monotonic-counter strategy used by the H15 analysis.

## Block-SEM diagnostic

Reference statistics include a block SEM.

Blocks are formed on the **complete original coordinate timeline for each
trajectory first**. State frames are selected only after block boundaries have
been established.

Thus the block definition is independent of state residence.

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
    --chemistry CHEMISTRY \
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
    [--lambda-mode single|his] \
    [--lambda-step-col 0] \
    [--lambda-value-col 1] \
    [--lambda-resid N] \
    [--x-low-tautomer HID|HIE] \
    [--x-high-tautomer HID|HIE] \
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

## Histidine reference command template

The exact His model-peptide residue and molecule masks must be verified from
each topology before running. Once known, the general command is:

```bash
python3 calc_refsasa.py \
    --forcefield ff19SB \
    --chemistry HIS \
    --site-definition histidine_sidechain \
    --selection 'HIS_REFERENCE_SELECTION' \
    -p MODEL_HIS.prmtop \
    -t 'arex.ph*.nc' \
    -l 'arex.ph*.lambda' \
    --solutemask 'MODEL_PEPTIDE_SOLUTEMASK' \
    --image-anchor 'MODEL_PEPTIDE_SOLUTEMASK' \
    --first-frame-step FIRST_STEP \
    --frame-step-interval STEP_INTERVAL \
    --lambda-mode his \
    --x-low-tautomer HID \
    --x-high-tautomer HIE \
    --jobs N \
    --workdir model_refsasa_work \
    -o refsasa.tsv \
    --force
```

For c22, reverse the two `x` endpoint mappings:

```text
--x-low-tautomer HIE
--x-high-tautomer HID
```

Do not copy the placeholder His masks into production commands without topology
verification.

## Persistent reference-table behavior

The same `refsasa.tsv` can be updated by multiple model-chemistry runs.

A new calculation replaces rows with the same key:

```text
forcefield + chemistry + site_definition + protonation
```

and leaves unrelated rows intact.

A shared table may therefore contain both ordinary two-state references and
His three-state references.

## Reference-table schema

The current table contains 30 columns:

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
lambda_mode
lambda_resid
x_low_tautomer
x_high_tautomer
```

The primary denominator field is:

```text
reference_sasa_A2
```

The last four fields make the state-classification provenance explicit. They
are blank where not applicable to legacy single-coordinate calculations.

## Audit output

Each chemistry/reference calculation writes:

```text
<workdir>/<forcefield>/<chemistry>/<site_definition>/
    reference_frame_assignments.tsv
    refsasa.settings.txt
```

The frame audit contains:

```text
forcefield
chemistry
site_definition
trajectory
pH
lambda_file
analysis_row
trajectory_frame
step
lambda
x
protonation
tautomer
clean
matched
sasa_A2
```

In single-coordinate mode, `x` and `tautomer` are blank. In His mode they make
the HIP/HID/HIE assignment independently auditable.

---

# 4. `calc_fsasa.py`

## Purpose

`calc_fsasa.py` combines:

- raw frame-wise protein SASA from `calc_sasa.py`;
- local state coordinates from an authoritative state table;
- chemistry/state-specific reference values from `calc_refsasa.py`.

It then writes frame-wise fSASA in the same canonical observable format used by
the rest of the analysis toolkit.

`calc_fsasa.py` does **not** read coordinates and does **not** run CPPTRAJ.
Changing a reference table or state assignment therefore does not require
rerunning the expensive coordinate-level SASA calculation.

Plain TSV and `.tsv.gz` state/input tables are supported transparently.

## Site mapping

Each protein observable must be mapped to:

```text
OBSERVABLE
CHEMISTRY
SITE_DEFINITION
LOCAL_PROTONATION_LAMBDA_COLUMN
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

For His mode, the fourth field remains the **protonation** lambda column. The
neutral-tautomer coordinate is supplied separately with:

```text
--tautomer-column x
```

## Single-coordinate normalization

`--lambda-mode single` is the default and preserves the original behavior:

```text
lambda <= 0.2       -> H
lambda >= 0.8       -> deprot
otherwise           -> local mixed
```

The denominator is looked up with:

```text
forcefield + chemistry + site_definition + H|deprot
```

A `HIS` site is rejected in single mode so HID/HIE cannot be silently pooled.

## Histidine normalization

Use:

```text
--lambda-mode his
```

for histidine.

The authoritative state table must contain the protonation coordinate named in
the `--site` mapping plus a tautomer coordinate, default:

```text
x
```

The physical-state rules are:

```text
lambda <= low                         -> HIP
lambda >= high and x <= low           -> x-low tautomer
lambda >= high and x >= high          -> x-high tautomer
low < lambda < high                   -> mixed
lambda >= high and low < x < high     -> neutral_tautomer_mixed
```

The endpoint mapping must match the corresponding reference calculation:

```text
c22:
    --x-low-tautomer HIE
    --x-high-tautomer HID

ff14SB / ff19SB:
    --x-low-tautomer HID
    --x-high-tautomer HIE
```

Clean His frames use one of three reference denominators:

```text
HIP
HID
HIE
```

Intermediate protonation-lambda frames and neutral frames with intermediate `x`
remain on the canonical timeline with blank fSASA.

By default, when the state table contains `state` and `tautomer` columns, the
newly classified His state is cross-checked against those labels. A mismatch is
fatal. This validation can be disabled only explicitly with:

```text
--no-validate-his-state-labels
```

## Exact state-table join

The canonical raw-SASA join fields are selected with:

```text
--join-on trajectory,trajectory_frame
```

by default.

State-side column names can be remapped without changing the canonical raw SASA
schema:

```text
--state-trajectory-column
--state-analysis-row-column
--state-trajectory-frame-column
--state-step-column
```

This is specifically useful for the H15 `analyze_his.py` output, whose
equivalent fields are:

```text
canonical raw SASA        H15 state table
------------------        ----------------
trajectory                trajectory
analysis_row               analysis_frame
trajectory_frame           raw_coord_frame
step                       target_md_step
```

If one table stores `arex.ph5_0` and the other stores `arex.ph5_0.nc`, use:

```text
--trajectory-key-mode stem
```

The normalized key removes only recognized trajectory-file extensions; it does
not blindly truncate labels containing a dot.

Shared frame metadata and pH are independently checked when available on both
tables. There is no nearest-frame or nearest-step matching.

## Matched and clean-state inference

If the requested `--matched-column` exists, it is respected.

If it does not exist, a row is considered matched when the mapped local
protonation-lambda field is populated. This allows direct use of state tables
such as the H15 frame table that encode unmatched rows through missing lambda
values rather than a separate `matched` flag.

Likewise, if the requested `--clean-column` exists, it is respected and
cross-checked against the locally classified state.

If it does not exist, local clean/mixed status is inferred directly from lambda
(and, for His neutral frames, `x`).

## Treatment of unavailable states

### Outside/unmatched state window

The row is retained in the canonical output, but `value` is blank.

Audit status:

```text
outside_state_window
```

### Intermediate protonation lambda

The row is retained, but `value` is blank.

Audit status:

```text
local_mixed
```

### His neutral frame with intermediate tautomer coordinate

The row is retained, but `value` is blank because neither the pure-HID nor
pure-HIE denominator applies.

Audit status:

```text
neutral_tautomer_mixed
```

### Clean local state

fSASA is calculated from the matching state-specific denominator.

For single-coordinate mode this is `H` or `deprot`.

For His mode this is `HIP`, `HID`, or `HIE`.

If a global `clean` flag exists and says a locally mixed frame is clean, the
script stops with an error.

A locally clean site may still retain a valid fSASA when a broader coupled-state
definition is globally unclean because another site is intermediate. This
preserves the original dyad behavior.

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
    --states STATE_ASSIGNMENTS.tsv[.gz] \
    --references refsasa.tsv \
    --forcefield FORCEFIELD \
    --site OBSERVABLE CHEMISTRY SITE_DEFINITION LAMBDA_COLUMN \
    [--site ...] \
    [--lambda-mode single|his] \
    [--tautomer-column x] \
    [--x-low-tautomer HID|HIE] \
    [--x-high-tautomer HID|HIE] \
    [--join-on trajectory,trajectory_frame] \
    [--trajectory-key-mode exact|stem] \
    [--state-trajectory-column trajectory] \
    [--state-analysis-row-column analysis_row] \
    [--state-trajectory-frame-column trajectory_frame] \
    [--state-step-column step] \
    [--low 0.2] \
    [--high 0.8] \
    [--clean-column clean] \
    [--matched-column matched] \
    [--state-class-column state] \
    [--tautomer-label-column tautomer] \
    [--no-validate-his-state-labels] \
    [--sasa-metric sasa_lcpo] \
    [--audit-output FILE] \
    -o FSASA.tsv
```

## Example: HEWL ff19SB E35/D52 fSASA

The validated single-coordinate dyad workflow remains unchanged:

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

These validated numbers refer to the catalytic-dyad analysis and are not
expected values for the H15 His analysis.

## H15 His normalization template

Once the H15 raw-SASA observable and reference table have been generated, the
H15 `analyze_his.py` state table can be used directly with column remapping.

The general pattern is:

```bash
python3 calc_fsasa.py \
    --sasa H15_raw_sasa.tsv \
    --states H15_FRAME_TABLE.tsv.gz \
    --references refsasa.tsv \
    --forcefield ff19SB \
    --site H15_SIDECHAIN HIS histidine_sidechain lambda \
    --lambda-mode his \
    --tautomer-column x \
    --x-low-tautomer HID \
    --x-high-tautomer HIE \
    --join-on trajectory,trajectory_frame \
    --trajectory-key-mode stem \
    --state-analysis-row-column analysis_frame \
    --state-trajectory-frame-column raw_coord_frame \
    --state-step-column target_md_step \
    -o H15_fsasa.tsv
```

For c22, reverse the two `x` endpoint mappings.

The `H15_SIDECHAIN` observable name and `histidine_sidechain` site-definition
string are examples; the production run script should use the exact names
chosen for the final H15 SASA definition.

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

The current audit records:

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
lambda_mode
lambda_column
local_lambda
tautomer_column
local_tautomer_coordinate
local_protonation
local_tautomer
reference_state
state_clean
state_matched
state_code
state_label
source_state
source_tautomer
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

It records the input files, join definition, state-column remapping, lambda
mode, threshold/mapping choices, site mappings, state-label validation setting,
and actual denominator values used.

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

## Histidine workflow template

The His workflow uses the same three-layer architecture but three clean
reference states rather than two.

### Step H1 — Calculate raw protein SASA

Run `calc_sasa.py` on the protein trajectories using the final verified H15 atom
selection and complete-protein `solutemask`.

`calc_sasa.py` remains state-agnostic. No His lambda or tautomer information is
used at this stage.

### Step H2 — Generate force-field-specific His references

Run `calc_refsasa.py --lambda-mode his` on the blocked His model-peptide
trajectories.

Generate `HIP`, `HID`, and `HIE` rows separately for each force field using the
same atom/site definition as the protein numerator.

The model-peptide masks and first-frame/step interval must be verified from the
actual His reference simulations before production use.

### Step H3 — Normalize H15 frame-wise SASA

Run `calc_fsasa.py --lambda-mode his` using the H15 frame/state table.

For the current H15 analysis, the state-side frame provenance can be remapped
from:

```text
analysis_frame
raw_coord_frame
target_md_step
```

to the canonical raw-SASA identity fields.

### Step H4 — Condition the resulting fSASA downstream

The fSASA tools stop at frame-wise normalization. H15-specific comparisons such
as:

```text
fSASA(HIE, chi2+)
fSASA(HIE, chi2-)
fSASA(HID, chi2+)
fSASA(HID, chi2-)
delta fSASA = fSASA(chi2+) - fSASA(chi2-)
```

belong in the downstream H15 analysis/plotting layer, not inside
`calc_fsasa.py`.

This separation is important because `HIP`/`HID`/`HIE` normalization and
`chi2` conditioning answer different questions.

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
- His tautomer coordinates,
- coupled states,
- reference normalization.

## `calc_refsasa.py`

Knows about:

- model-system coordinates,
- exact lambda-step matching,
- chemistry-level reference identity,
- state-specific reference statistics,
- legacy one-coordinate `H`/`deprot` classification,
- His `ires`/`itauto` parsing and `HIP`/`HID`/`HIE` classification.

Does not know about:

- biological identities such as HEWL E35, D52, or H15,
- protein rotamers,
- protein coupled-state labels,
- protein fSASA.

## `calc_fsasa.py`

Knows about:

- canonical raw SASA tables,
- authoritative local state coordinates,
- chemistry/site mappings,
- reference-table lookup,
- exact/remapped frame joins,
- single-coordinate and His state classification,
- frame-wise normalization.

Does not know about:

- coordinates,
- CPPTRAJ,
- how the state table itself was generated,
- `chi2` rotamer definitions,
- state-conditioned averaging.

## `partition_by_state.py` and project-specific analysis

Downstream statistics tools know about:

- authoritative state assignments,
- exact frame joins,
- timeline-first block statistics,
- requested grouping/conditioning variables.

They do not need to know whether the normalized observable came from:

- distance,
- SASA,
- fSASA,
- hydrogen bonding,
- hydration,
- or another structural quantity.

For specialized analyses such as H15 rotamer-conditioned fSASA, a
project-specific plotting/analysis layer may be more appropriate than the
generic coupled-state partitioner.

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

but the state-side names of the canonical identity columns can be remapped with
the `--state-*-column` options.

### Dyad-style state table

The validated HEWL catalytic-dyad table contains fields including:

```text
analysis_row
trajectory_frame
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

For that analysis:

```text
primary_label   = E35
primary_lambda  = E35 local lambda
coupled_label   = D52
coupled_lambda  = D52 local lambda
```

Do not infer protonation from the label columns.

### H15 `analyze_his.py` frame table

The H15 table uses different frame-column names and contains the physical His
state directly:

```text
trajectory
analysis_frame
raw_coord_frame
target_md_step
ph
lambda
x
state
tautomer
chi2_deg
...
```

For fSASA joining, the relevant remapping is:

```text
canonical                  H15 state table
---------                  ---------------
trajectory                 trajectory
analysis_row                analysis_frame
trajectory_frame            raw_coord_frame
step                        target_md_step
```

The H15 trajectory field may include a trajectory extension while
`calc_sasa.py` stores the trajectory stem. Use `--trajectory-key-mode stem`
when necessary.

If the H15 table has no explicit `matched` or `clean` columns,
`calc_fsasa.py` infers these from the availability and endpoint classification
of the local state coordinates.

## Reference table

At minimum, `calc_fsasa.py` requires:

```text
forcefield
chemistry
site_definition
protonation
reference_sasa_A2
```

Required state rows are:

```text
single mode: H, deprot
His mode:    HIP, HID, HIE
```

The production `calc_refsasa.py` writes substantially more provenance,
including lambda mode and His endpoint mapping.

---

# 8. Validation checklist

Before accepting a new fSASA analysis, check all of the following.

### Raw SASA

- every expected trajectory completed;
- expected frame count is present for every observable;
- raw output ordering is deterministic;
- finite negative partial-LCPO values are preserved;
- `offset`, `nbrcut`, selections, solutemask, topology, and frame-step mapping
  are recorded;
- the protein and reference calculations use the **same chemical atom/site
  definition**.

### Reference SASA — all modes

- every coordinate frame is accounted for as matched or outside the lambda
  range;
- every step inside the lambda range has an exact lambda record;
- no nearest-neighbor lambda matching occurs;
- every required reference mean is finite and > 0;
- reference force field matches the protein force field;
- state thresholds match those used by `calc_fsasa.py`.

### Reference SASA — single mode

- every matched frame is classified as `H`, `deprot`, or `mixed`;
- `H + deprot + mixed = matched`;
- both clean `H` and `deprot` ensembles contain frames.

### Reference SASA — His mode

- the intended His `ires` is resolved;
- exactly one `itauto=1` and one `itauto=2` variable are selected for that
  `ires`;
- `x` endpoint mapping is explicitly correct for the force field;
- every matched frame is classified as `HIP`, `HID`, `HIE`,
  protonation-`mixed`, or `neutral_tautomer_mixed`;
- all three `HIP`, `HID`, and `HIE` reference ensembles contain frames;
- neutral tautomer-mixed frames are excluded from the pure HID/HIE
  denominators;
- any lambda counter resets are reported and the unwrapped step series remains
  strictly increasing.

### fSASA — all modes

- every raw SASA observation exact-joins to one authoritative state row;
- no nearest-neighbor frame matching occurs;
- shared frame/step/pH metadata agree;
- unavailable-state rows remain in the canonical timeline with blank fSASA;
- locally clean rows use the correct chemistry/state denominator;
- negative raw values remain negative after division by positive references.

### fSASA — His mode

- `lambda` and `x` endpoint mappings match the reference calculation;
- `HIP` uses the HIP denominator independent of `x`;
- pure neutral frames use HID or HIE denominators separately;
- neutral tautomer-mixed frames do not receive an fSASA value;
- source `state`/`tautomer` labels agree with the lambda/x classification when
  validation is enabled;
- H15 frame-column remapping and trajectory stem normalization are recorded in
  the settings file.

### Downstream state/rotamer partition

- blocks are formed from the intended original timeline before conditional
  state selection;
- sampling thresholds are applied to the actual microstate being summarized;
- for H15, HID/HIE and `chi2+`/`chi2-` conditioning is performed downstream,
  not by changing the reference denominator definition;
- rare one-frame states retain blank SD/SEM rather than fabricated uncertainty.

---

# 9. Current ff19SB HEWL carboxylate validation summary

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

## Raw protein SASA

```bash
python3 calc_sasa.py \
    -p topology.prmtop \
    -t 'arex.ph*.nc' \
    --site OBSERVABLE 'CPPTRAJ_SELECTION' \
    --solutemask 'SOLUTE_MASK' \
    --first-frame-step FIRST_STEP \
    --frame-step-interval STEP_INTERVAL \
    --jobs N \
    -o protein_sasa.tsv
```

## Single-coordinate reference denominator

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

## Histidine reference denominator

```bash
python3 calc_refsasa.py \
    --forcefield ff19SB \
    --chemistry HIS \
    --site-definition histidine_sidechain \
    --selection 'HIS_REFERENCE_SELECTION' \
    -p model_his.prmtop \
    -t 'arex.ph*.nc' \
    -l 'arex.ph*.lambda' \
    --solutemask 'MODEL_PEPTIDE_SOLUTEMASK' \
    --first-frame-step FIRST_STEP \
    --frame-step-interval STEP_INTERVAL \
    --lambda-mode his \
    --x-low-tautomer HID \
    --x-high-tautomer HIE \
    --jobs N \
    -o refsasa.tsv
```

For c22, use:

```text
--x-low-tautomer HIE
--x-high-tautomer HID
```

## Single-coordinate frame-wise fSASA

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

## H15-style His frame-wise fSASA

```bash
python3 calc_fsasa.py \
    --sasa H15_raw_sasa.tsv \
    --states H15_frames.tsv.gz \
    --references refsasa.tsv \
    --forcefield ff19SB \
    --site H15_SIDECHAIN HIS histidine_sidechain lambda \
    --lambda-mode his \
    --tautomer-column x \
    --x-low-tautomer HID \
    --x-high-tautomer HIE \
    --trajectory-key-mode stem \
    --state-analysis-row-column analysis_frame \
    --state-trajectory-frame-column raw_coord_frame \
    --state-step-column target_md_step \
    -o H15_fsasa.tsv
```

The His masks, observable name, site-definition string, model-peptide masks, and
frame-step mapping in these templates are placeholders until verified for the
specific production system.

## Generic state-conditioned fSASA

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

For H15 rotamer/tautomer-conditioned analysis, use the project-specific H15
analysis layer so the same `chi2` basin definitions and sampling thresholds are
used as in the rest of the H15 figure analysis.

