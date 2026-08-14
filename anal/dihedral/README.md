# calc_dihedral.py

General protein side-chain dihedral and rotamer observable generator for Amber trajectories.

`calc_dihedral.py` calculates standard amino-acid side-chain χ torsions or user-defined four-atom dihedrals using CPPTRAJ and writes a long-format TSV designed to integrate directly with `partition_by_state.py`.

The analyzer is deliberately **state-agnostic**: it measures geometry only. Protonation-state or microscopic-state conditioning is performed downstream by joining the output to an authoritative frame-state table.

The current implementation uses CPPTRAJ for the actual torsion calculation, wraps angles to `[-180, 180)`, emits rotamer indicator variables, and provides sine/cosine components for statistically correct circular analysis.

## Requirements

- Python 3
- AmberTools / CPPTRAJ available on `$PATH`
- Amber topology (`.prmtop`)
- One or more Amber trajectories, typically NetCDF (`.nc`)

## Basic usage

```bash
python3 calc_dihedral.py \
    -p system.prmtop \
    -t "arex.ph*.nc" \
    --residue E35 ':36' \
    --residue D52 ':53' \
    --torsions all \
    --image-anchor '^1' \
    --first-frame-step 510000 \
    --frame-step-interval 10000 \
    --jobs 4 \
    --workdir dihedral_work \
    -o dihedrals.tsv
```

## HEWL ff19SB example

```bash
ff19_hewl=/home/wayyne/cphmd/amber24/ff19/ti/harris_2022/prot/hewl/val-take2/rex-cont

python3 calc_dihedral.py \
  -p ${ff19_hewl}/../../pre/*.prmtop \
  -t "${ff19_hewl}/arex.ph*.nc" \
  --residue H15 ':16' \
  --residue E35 ':36' \
  --residue D52 ':53' \
  --residue N46 ':47' \
  --residue N59 ':60' \
  --torsions all \
  --image-anchor '^1' \
  --first-frame-step 510000 \
  --frame-step-interval 10000 \
  --jobs 9 \
  --workdir ff19_dihedral_work_v1_0 \
  -o ff19_hewl_dihedrals_v1_0.tsv
```

For these five residues, the program requests 11 torsions:

```text
H15   HIS   chi1 chi2
E35   GLU   chi1 chi2 chi3
D52   ASP   chi1 chi2
N46   ASN   chi1 chi2
N59   ASN   chi1 chi2
```

## Standard χ definitions

Examples:

```text
HIS chi1 = N-CA-CB-CG
HIS chi2 = CA-CB-CG-ND1

ASP chi1 = N-CA-CB-CG
ASP chi2 = CA-CB-CG-OD1

GLU chi1 = N-CA-CB-CG
GLU chi2 = CA-CB-CG-CD
GLU chi3 = CB-CG-CD-OE1

ASN chi1 = N-CA-CB-CG
ASN chi2 = CA-CB-CG-OD1
```

The H15 definitions intentionally reproduce the prior HEWL H15 analysis:

```text
chi1 = N-CA-CB-CG
chi2 = CA-CB-CG-ND1
```

Supported residue classes include standard and several CpHMD-specific Amber residue names such as ASP/ASH/AS2/AS4, GLU/GLH/GL2/GL4, HIS/HID/HIE/HIP, CYS/CYM, LYS/LYN, TYR/TYM, and common non-titratable amino acids.

## Selecting torsions

Analyze every standard χ torsion:

```bash
--torsions all
```

Analyze selected torsions:

```bash
--torsions chi1,chi2
```

The same `--torsions` selection is applied to every `--residue` argument.

## Arbitrary dihedrals

Custom four-atom torsions are supported:

```bash
--dihedral NAME MASK1 MASK2 MASK3 MASK4
```

Example:

```bash
--dihedral custom_tor ':36@CA' ':36@CB' ':36@CG' ':36@CD'
```

To also classify a custom torsion into coarse g−/g+/trans bins:

```bash
--classify-custom
```

## Output format

Columns:

```text
trajectory
pH
analysis_row
trajectory_frame
step
observable
residue_label
residue_mask
residue_number
residue_name
torsion
atom1
atom2
atom3
atom4
symmetry_period_deg
metric
value
unit
```

`trajectory` is written using the trajectory filename stem, e.g. `arex.ph6_5.nc` becomes `arex.ph6_5`.

`trajectory_frame` is one-based.

`analysis_row` is a **frame-level identifier**. All torsions and metrics belonging to the same trajectory frame share the same `analysis_row`.

The MD step is reconstructed as:

```text
step = first_frame_step
     + (trajectory_frame - 1) * frame_step_interval
```

## Metrics

For each standard torsion:

```text
angle_deg
sin_angle
cos_angle
rotamer_gminus
rotamer_gplus
rotamer_trans
```

Rotamer metrics are binary 0/1 values. Their state-conditioned means therefore equal rotamer populations.

## Angle convention

All angles are wrapped to:

```text
[-180°, 180°)
```

## Rotamer convention

```text
g-      [-120°,   0°)
g+      [   0°, 120°)
trans   [ 120°, 180°) U [-180°, -120°)
```

These are general coarse bins and should not replace residue-specific historical definitions where those are needed.

## Circular statistics

Do **not** interpret the arithmetic mean of `angle_deg` as a mean torsion angle.

Use the state-conditioned means of:

```text
sin_angle
cos_angle
```

and compute:

```text
theta_mean = atan2(mean_sin, mean_cos)
R = sqrt(mean_sin^2 + mean_cos^2)
```

`R` near 1 indicates a tightly localized angular distribution; smaller `R` indicates a broader or potentially multimodal distribution.

## Terminal-group symmetry

The registry can mark selected terminal torsions as having twofold atom-label symmetry. Such torsions may emit:

```text
symmetry_angle_deg
```

Examples include standard ASP chi2, GLU chi3, PHE chi2, TYR chi2, and TYM chi2.

### CpHMD carboxylates

For CpHMD carboxylate residue types such as `AS2` and `GL2`, terminal-group symmetry should be treated carefully. A protonated carboxylic acid does not have the same instantaneous oxygen equivalence as a deprotonated carboxylate. For these cases, symmetry-aware analysis is best performed downstream using protonation-state information.

## State conditioning

`calc_dihedral.py` does not read lambda coordinates and does not assign protonation states.

Use the authoritative state table downstream with `partition_by_state.py`.

For the HEWL E35/D52 dyad:

```text
S0       E35- / D52-
S1_E35   E35H / D52-
S1_D52   E35- / D52H
S2       E35H / D52H
```

An H15, N46, or N59 result partitioned with this table means the conformation of that residue conditional on the E35/D52 dyad state. It does not report that residue's own protonation state.

## Parallel execution

Trajectory-level parallelism:

```bash
--jobs N
```

CPPTRAJ threading per worker:

```bash
--cpptraj-threads N
```

Example:

```bash
--jobs 9 --cpptraj-threads 1
```

## Imaging

Optional CPPTRAJ autoimage anchor:

```bash
--image-anchor '^1'
```

## pH parsing

Default regex:

```text
ph(?P<ph>\d+(?:[_\.]\d+)?)
```

Examples:

```text
arex.ph5.nc
arex.ph5_5.nc
arex.ph6.5.nc
```

parse as pH 5, 5.5, and 6.5.

A custom regex can be supplied with `--ph-regex`; it must contain a named capture group `ph`.

## Work directory

Intermediate CPPTRAJ inputs, logs, topology information, and torsion series are retained under:

```bash
--workdir dihedral_work
```

This provides an audit trail and simplifies debugging.

## Settings file

Every successful run also writes:

```text
OUTPUT.tsv.settings.txt
```

This records program/version information, topology, trajectory count, torsion count, frame-step mapping, image anchor, output-row count, rotamer boundaries, circular-mean formula, torsion definitions, and symmetry settings.

## Expected output accounting

For `Ntraj` trajectories, `Nframe` frames per trajectory, `Ntors` torsions, and 6 standard metrics:

```text
Ntraj × Nframe × Ntors × 6
```

For the ff19 HEWL run:

```text
18 × 5150 × 11 × 6 = 6,118,200 rows
```

plus the TSV header.

## Validation

Recommended checks:

1. Confirm output row count matches expected frame × torsion × metric accounting.
2. Confirm trajectory identifiers exactly match the authoritative state table.
3. Confirm all metric rows from one trajectory frame share the same `analysis_row`.
4. Run `partition_by_state.py` and verify exact joined/unmatched/unclean accounting.
5. For H15, compare raw chi1/chi2 values against the previous dedicated H15 analyzer on identical frames.

For HEWL H15:

```text
chi1 = N-CA-CB-CG
chi2 = CA-CB-CG-ND1
```

## Version

The uploaded source currently reports version `1.0.1`.

If the corrected frame-level `analysis_row` and corrected output-row counter are designated as the finalized release, bump the source consistently before archival, e.g. to `1.0.2`.

## Design philosophy

```text
trajectory geometry
        |
        v
calc_dihedral.py
        |
        v
state-agnostic observable TSV
        |
        + authoritative state assignments
        |
        v
partition_by_state.py
        |
        v
state-conditioned structural statistics
```

This keeps geometry generation independent of a particular CpHMD state definition and allows the same observable table to be reused for different mechanistic questions.

