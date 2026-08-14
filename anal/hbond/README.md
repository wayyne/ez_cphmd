# `calc_hbond.py`

Lambda-aware hydrogen-bond analysis for Amber continuous constant-pH molecular dynamics (CpHMD).

`calc_hbond.py` is a general hydrogen-bond observable generator designed for CpHMD trajectories. It uses CPPTRAJ to detect **geometric candidate hydrogen bonds**, then applies the instantaneous CpHMD protonation and tautomer state to decide whether each donor/acceptor assignment is chemically valid in each frame.

The production implementation described here is **version 5.1**. The installed production filename should be:

```text
calc_hbond.py
```

The versioned development file was `calc_hbond_v5_1.py`.

---

## 1. Why this tool is CpHMD-aware

A conventional hydrogen-bond analysis assumes that donor and acceptor identities are fixed by the topology. That assumption is not valid for constant-pH MD.

A titratable atom can change chemical role during the trajectory. For example, a carboxylate oxygen can be an acceptor in the deprotonated state and, depending on the tautomer, a donor in the protonated state. Histidine ring nitrogens likewise change donor/acceptor roles with protonation and tautomer state.

Therefore, `calc_hbond.py` requires both:

```text
coordinate trajectory (.nc)
CpHMD lambda trajectory (.lambda)
```

for normal CpHMD use.

CPPTRAJ is deliberately used only to answer the geometric question:

> Does this donor-H-acceptor triplet satisfy the requested distance and angle criteria in this frame?

Python then answers the chemical question:

> Is the donor/acceptor assignment valid for the instantaneous CpHMD state of the residues involved?

This separation avoids treating hydrogens present in a constant-pH topology as if they necessarily represent the physical protonation state of every frame.

---

## 2. Design overview

The production workflow is:

```text
Amber topology
     +
trajectory (.nc)
     +
matching lambda file
     |
     v
calc_hbond.py
     |
     |-- expand Amber masks using CPPTRAJ
     |-- enumerate candidate donor heavy atoms
     |-- detect geometric H-bond series with CPPTRAJ
     |-- exact-join each coordinate frame to lambda MD step
     |-- apply frame-specific protonation/tautomer chemistry
     |-- retain definite H-bonds
     |-- separately count chemically unresolved candidates
     v
frame-wise canonical TSV
     +
atom-level interaction TSV
```

The candidate-generation engine is reported as:

```text
donor-enumerated-v5.1
```

This is important. Earlier broad-mask implementations could miss valid atom-level H-bond series that were recovered when the same interaction was specified explicitly. The production engine therefore enumerates candidate donor heavy atoms and builds the broad observables from atom-level directed searches.

For example,

```bash
--focus DYAD ':36,53' '^1'
```

is conceptually evaluated as the union of:

```text
every possible donor in :36,53 -> possible acceptors in ^1
every possible donor in ^1     -> possible acceptors in :36,53
```

with physical donor-H-acceptor triplets deduplicated afterward.

---

## 3. Requirements

The script requires:

- Python 3
- Amber topology (`prmtop`/`parm7`)
- Amber trajectory readable by CPPTRAJ
- CPPTRAJ on `PATH`, or supplied explicitly with `--cpptraj`
- matching CpHMD lambda files for CpH-filtered analysis

In the Shen-lab Amber environment, source Amber before running if necessary, for example:

```bash
source /opt/shen/shared/apps/amber/ambertools25_gcc12_avx2_sm75_86/amber.sh
which cpptraj
```

Alternatively:

```bash
--cpptraj /full/path/to/cpptraj
```

The script checks that CPPTRAJ exists before beginning analysis.

---

## 4. Default H-bond definition

The defaults are:

```text
heavy-atom distance cutoff = 3.0 A
angle cutoff               = 135 degrees
```

These are controlled with:

```bash
--distance 3.0
--angle 135
```

The distance criterion is applied by CPPTRAJ to the donor-heavy/acceptor geometry, and the angular criterion is the CPPTRAJ hydrogen-bond angular definition.

The angle may be disabled using CPPTRAJ-compatible behavior with:

```bash
--angle -1
```

if a geometry-only control is desired.

---

## 5. Exact frame/lambda synchronization

CpH filtering requires exact synchronization between the coordinate trajectory and lambda trajectory.

Supply:

```bash
--first-frame-step STEP
--frame-step-interval INTERVAL
```

For analyzed trajectory frame `k`, the corresponding MD step is constructed exactly from these values and the frame-selection stride.

The lambda file is indexed by its explicit MD-step column. The join is **exact**.

The program never uses nearest-neighbor lambda matching.

Behavior is:

```text
exact lambda step exists
    -> use that lambda row

frame lies outside lambda-file step range
    -> retain frame; any H-bond requiring unavailable CpH chemistry is
       counted as outside-lambda uncertainty

step lies inside lambda-file range but exact row is absent
    -> error
```

Both `--first-frame-step` and `--frame-step-interval` are required in CpH-filtered mode.

---

## 6. Lambda-file interpretation

The parser reads the Amber CpHMD lambda headers:

```text
# ititr ...
# ires  ...
# itauto ...
```

and associates lambda coordinates with topology residue IDs through `ires`.

Default clean-state thresholds are:

```text
lambda <= 0.2  -> low state
lambda >= 0.8  -> high state
0.2 < lambda < 0.8 -> mixed
```

These are configurable:

```bash
--lambda-low 0.2
--lambda-high 0.8
```

For the built-in CpH models, low/high are interpreted according to the residue-state model below.

---

## 7. Built-in CpH residue chemistry

The chemistry registry is general and residue-model driven. It is not specific to HEWL, E35/D52, Asp, or Glu.

Supported residue families are:

```text
ASP family : ASP, ASH, AS2, AS4
GLU family : GLU, GLH, GL2, GL4
HIS family : HIS, HID, HIE, HIP
CYS family : CYS, CYM
LYS family : LYS, LYN
TYR family : TYR, TYM
```

Only atoms whose donor/acceptor chemistry actually changes with CpH state are lambda-dependent. Backbone atoms on titratable residues remain ordinary fixed-chemistry protein atoms.

### 7.1 Asp and Glu

Amber/Shen continuous-CpH convention used by the code:

```text
itauto = 3 -> protonation lambda
itauto = 4 -> carboxylate tautomer coordinate
```

For Asp:

```text
O1 = OD1
O2 = OD2
```

For Glu:

```text
O1 = OE1
O2 = OE2
```

Chemistry:

```text
protonation lambda high:
    O1 acceptor
    O2 acceptor
    neither is donor

protonation lambda low, tautomer low:
    O2-H donor
    O1 acceptor

protonation lambda low, tautomer high:
    O1-H donor
    O2 acceptor
```

Thus the carboxylate tautomer convention is:

```text
tautomer low  -> O2-H
tautomer high -> O1-H
```

If the protonation coordinate is mixed, the side-chain donor/acceptor assignment is unresolved. If the residue is protonated but the tautomer coordinate is mixed, the side-chain assignment is likewise unresolved.

### 7.2 Histidine

Convention:

```text
itauto = 1 -> protonation lambda
itauto = 2 -> tautomer coordinate
```

Chemistry:

```text
protonation lambda low:
    HIP-like state
    ND1 donor
    NE2 donor
    neither ring N is an acceptor

protonation lambda high, tautomer low:
    HID-like state
    ND1 donor
    NE2 acceptor

protonation lambda high, tautomer high:
    HIE-like state
    NE2 donor
    ND1 acceptor
```

### 7.3 Cys

```text
protonated:
    SG donor

deprotonated:
    SG acceptor
```

### 7.4 Lys

```text
protonated LysH+:
    NZ donor
    NZ not acceptor

neutral Lys:
    NZ donor
    NZ acceptor
```

### 7.5 Tyr

```text
protonated Tyr-OH:
    OH donor
    OH acceptor under the protein chemistry profile

deprotonated Tyr-O-:
    OH/O- site acceptor only
```

### 7.6 Backbone atoms on titratable residues

A titratable residue does **not** imply that all atoms in that residue depend on lambda.

Examples:

```text
GLU/GL2 backbone N-H -> fixed donor chemistry
GLU/GL2 backbone O   -> fixed acceptor chemistry
ASP/AS2 backbone N-H -> fixed donor chemistry
ASP/AS2 backbone O   -> fixed acceptor chemistry
```

Only the actual titrating side-chain atoms listed in the registry are state-dependent.

---

## 8. Candidate chemistry profiles

Default:

```bash
--chemistry-profile protein
```

This uses standard residue/atom-specific protein donor and acceptor candidates plus the superset of possible CpH roles.

A more permissive fallback is available:

```bash
--chemistry-profile fon
```

which uses an F/O/N-heavy candidate mask before CpH filtering.

The `protein` profile is recommended for protein CpHMD analysis.

Methionine sulfur is excluded by default. To include fixed Met `SD` as an acceptor candidate:

```bash
--include-met-sulfur
```

---

## 9. Analysis modes

Multiple observables may be requested in a single invocation. The options are repeatable.

### 9.1 `--scan`

```bash
--scan NAME MASK
```

Find H-bonds within a selected region.

Example:

```bash
--scan PROTEIN '^1'
```

This is useful for broad H-bond counts or discovery.

### 9.2 `--focus`

```bash
--focus NAME FOCUS_MASK PARTNER_MASK
```

Find all H-bonds for which one endpoint belongs to `FOCUS_MASK` and the other belongs to `PARTNER_MASK`, in either donor/acceptor direction.

Example:

```bash
--focus DYAD ':36,53' '^1'
```

This is the preferred broad mechanism-discovery mode for a residue/site of interest.

### 9.3 `--interaction`

```bash
--interaction NAME MASK1 MASK2
```

Analyze H-bonds between two masks in either direction.

This is useful when two groups are known but donor/acceptor direction should not be imposed beforehand.

### 9.4 `--directed`

```bash
--directed NAME DONOR_MASK ACCEPTOR_MASK
```

Analyze one explicit donor-to-acceptor direction.

Example:

```bash
--directed N46_to_D52 ':47@ND2' ':53@OD1,OD2'
```

This is useful for targeted mechanistic observables or validation of a broader `--focus` calculation.

### 9.5 `--directed-h`

```bash
--directed-h NAME DONOR_MASK DONOR_H_MASK ACCEPTOR_MASK
```

Explicitly specify donor heavy atom(s), donor hydrogen(s), and acceptor(s).

Use this only when an atom-specific hydrogen definition is required. In ordinary directed analysis, allowing CPPTRAJ to identify hydrogens bonded to the donor heavy atom is usually preferable.

### 9.6 `--solvent-site`

```bash
--solvent-site NAME SOLUTE_MASK
```

Analyze solute-site hydrogen bonds involving the configured solvent donor and/or acceptor masks.

Defaults:

```text
solvent donor mask    : :WAT
solvent acceptor mask : :WAT@O
```

Override with:

```bash
--solvent-donor-mask MASK
--solvent-acceptor-mask MASK
```

Use `calc_hbond.py` for hydrogen-bond geometry/occupancy. Water-shell coordination or hydration-number analysis should remain a separate observable rather than being inferred from H-bond counts.

---

## 10. Frame selection

Coordinate frames may be selected with:

```bash
--start 1
--stop last
--stride 1
```

`--start` and numeric `--stop` use trajectory-frame numbering expected by the script/CPPTRAJ workflow.

The canonical output records both:

```text
analysis_row
trajectory_frame
```

so a subsampled calculation remains exactly auditable.

When using a stride, the calculated `step` follows the selected original trajectory frame, not merely the compact output row number.

---

## 11. Imaging

Optional imaging controls:

```bash
--image-anchor MASK
--hbond-image
```

For protein systems, an image anchor such as:

```bash
--image-anchor '^1'
```

may be appropriate to keep the protein molecule whole/consistently imaged before H-bond analysis.

Imaging should be chosen consistently with the trajectory-preparation and other structural analyses used in the project.

---

## 12. Fixed-protonation controls

To run without CpH filtering:

```bash
--no-cph-filter
```

In this mode:

- `--lambda` is not required;
- exact lambda step mapping is not required;
- residue chemistry is evaluated using fixed topology/residue-name chemistry.

This mode is useful for controls and conventional MD trajectories. It should not be used as a substitute for lambda-aware analysis of a genuine continuous-CpHMD trajectory.

---

## 13. Frame-wise output

The main TSV begins with the canonical frame identifiers:

```text
trajectory
pH
analysis_row
trajectory_frame
step
observable
metric
interaction_class
mode
selection1
selection2
value
unit
chemistry_status
lambda_file
```

Each requested observable produces **six rows per analyzed frame**, one for each metric.

### 13.1 `hbond_count`

```text
unit: count
```

Number of definitely chemically valid H-bond triplets present in the frame.

This value is never erased merely because another candidate in the same frame has unresolved CpH chemistry.

### 13.2 `hbond_present`

```text
unit: 1
```

Binary indicator:

```text
1 -> at least one definitely valid H-bond is present
0 -> no definitely valid H-bond is present
```

### 13.3 `ambiguous_hbond_count`

```text
unit: count
```

Number of geometrically present candidates whose donor/acceptor chemistry cannot be resolved because a required protonation or tautomer coordinate lies between the clean-state thresholds.

### 13.4 `outside_lambda_hbond_count`

```text
unit: count
```

Number of geometrically present candidates whose chemical validity depends on lambda data unavailable because the coordinate frame lies outside the matched lambda-file step range.

### 13.5 `hbond_count_min`

```text
unit: count
```

Lower bound on the true H-bond count.

Currently:

```text
hbond_count_min = hbond_count
```

### 13.6 `hbond_count_max`

```text
unit: count
```

Upper bound if every unresolved candidate were chemically valid:

```text
hbond_count_max =
    hbond_count
  + ambiguous_hbond_count
  + outside_lambda_hbond_count
```

---

## 14. `chemistry_status`

Each observable/frame receives one status:

```text
clean
partial_ambiguity
outside_lambda_window
```

Definitions:

```text
clean
    no present H-bond candidate requires unresolved CpH chemistry

partial_ambiguity
    at least one present candidate requires a mixed protonation/tautomer
    coordinate, but no candidate requires unavailable lambda data

outside_lambda_window
    at least one present candidate requires lambda information outside the
    available lambda-file step range
```

`outside_lambda_window` takes precedence over `partial_ambiguity`.

The status belongs to the **frame-observable**, so it is duplicated identically across all six metrics for that frame.

For a completely clean observable/frame:

```text
ambiguous_hbond_count      = 0
outside_lambda_hbond_count = 0
hbond_count_min            = hbond_count_max = hbond_count
```

---

## 15. Atom-level interaction report

By default, the script also writes:

```text
<output_stem>.interactions.tsv
```

or a custom location supplied through:

```bash
--detail-output FILE
```

Fields include:

```text
trajectory
pH
observable
mode
interaction_class
selection1
selection2
acceptor
donor_h
donor
count
count_semantics
fraction
raw_count
raw_fraction
ambiguous_present_frames
outside_lambda_present_frames
raw_avg_distance_A
raw_avg_angle_deg
geometry_average_scope
```

### 15.1 Chemically filtered occupancy

`count` is:

```text
chemically_valid_frames_present
```

That is, the number of frames in which the triplet is both geometrically present **and** chemically valid according to the matched CpH state.

`fraction` is:

```text
count / number_of_analyzed_frames
```

### 15.2 Raw geometric occupancy

`raw_count` and `raw_fraction` report the CPPTRAJ geometric candidate before CpH filtering.

Comparing raw and filtered values is a useful audit of the effect of protonation/tautomer chemistry.

### 15.3 Ambiguous and outside-lambda occurrences

For each atom-level triplet:

```text
ambiguous_present_frames
outside_lambda_present_frames
```

report how often the geometry was present but the CpH chemistry could not be definitively classified.

### 15.4 Geometry averages

The columns:

```text
raw_avg_distance_A
raw_avg_angle_deg
```

are CPPTRAJ averages over **geometric candidate-present frames before CpH filtering**.

They are deliberately labelled `raw_...` because they are not recomputed only over chemically valid frames.

Do not describe these values in a manuscript as CpH-filtered distance/angle averages.

If chemically filtered geometry distributions are needed, they should be generated explicitly from frame-resolved geometry rather than inferred from these raw CPPTRAJ averages.

---

## 16. Output settings and reproducibility

The script writes:

```text
<output>.settings.txt
```

including:

```text
Program version
analysis type
chemistry profile
candidate engine
CpH filtering state
lambda thresholds
distance cutoff
angle cutoff
number of trajectories
number of requested observables
trajectory -> lambda-file mapping
```

Per-trajectory work directories additionally contain job settings used to detect accidental reuse of incompatible intermediate results.

If an existing work directory contains different settings, the program fails rather than silently mixing analyses.

Use either:

```bash
--force
```

or, preferably for materially different analyses, a new `--workdir`.

---

## 17. Parallel execution

Use:

```bash
-j N
```

or equivalently:

```bash
--jobs N
--fork N
```

Parallelism is across trajectories. Output ordering is deterministic and does not depend on completion order.

For one trajectory, increasing `--jobs` provides no trajectory-level speedup.

---

## 18. Multiple trajectories and lambda matching

`-t/--trajectory` and `--lambda` may be repeated and may contain quoted globs.

Example:

```bash
-t 'arex.ph*.nc' \
--lambda 'rex/arex.ph*.lambda'
```

Trajectory/lambda matching is performed by normalized filename identity/pH-aware matching logic in the script. The final mapping is recorded in the settings file.

Always inspect the settings file for a new production system to confirm that each trajectory was paired with the intended lambda file.

---

## 19. Production example: HEWL dyad

Example for one ff19SB HEWL trajectory:

```bash
ff19_hewl=/home/wayyne/cphmd/amber24/ff19/ti/harris_2022/prot/hewl/val-take2/rex-cont

python3 /home/wayyne/cphmd/ez_cphmd/anal/hbond/calc_hbond.py \
    -p ${ff19_hewl}/../../pre/*.prmtop \
    -t "${ff19_hewl}/arex.ph6_5.nc" \
    --lambda "${ff19_hewl}/rex/arex.ph6_5.lambda" \
    --scan PROTEIN '^1' \
    --focus DYAD ':36,53' '^1' \
    --image-anchor '^1' \
    --first-frame-step 510000 \
    --frame-step-interval 10000 \
    --jobs 1 \
    --workdir hbond_work \
    -o hewl_hbonds.tsv
```

Adjust the lambda path to the actual directory layout if needed.

For this topology convention:

```text
biological E35 -> topology residue 36
biological D52 -> topology residue 53
biological N46 -> topology residue 47
biological N59 -> topology residue 60
```

Those offsets are system-specific and are **not** hardcoded in `calc_hbond.py`.

---

## 20. Targeted mechanistic examples

### N46 -> D52

```bash
--directed N46_to_D52 ':47@ND2' ':53@OD1,OD2'
```

### N59 -> D52

```bash
--directed N59_to_D52 ':60@ND2' ':53@OD1,OD2'
```

The broad focus mode should recover these interactions automatically when they satisfy the geometric criteria:

```bash
--focus DYAD ':36,53' '^1'
```

The targeted forms remain useful as diagnostic or manuscript-specific observables because they produce a scalar frame-wise count/presence for precisely defined interactions.

---

## 21. Important counting semantics

The atom-level report represents individual donor-H-acceptor triplets.

Do **not** blindly add atom-level occupancy fractions to obtain a residue-pair occupancy. Multiple triplets can occur in the same frame.

For example, an Asn side-chain donor has two donor hydrogens, and a carboxylate has two oxygens. Several triplets can therefore represent the same residue-pair interaction.

For residue-pair questions, prefer the frame-wise observable generated by a targeted `--directed` or `--interaction` definition and use:

```text
hbond_count
hbond_present
```

according to the scientific question.

Typical interpretations:

```text
hbond_present
    fraction of frames with at least one H-bond between the defined groups

hbond_count
    instantaneous number of valid donor-H-acceptor triplets between the groups
```

The distinction should be stated explicitly in methods and figure legends.

---

## 22. Integration with `partition_by_state.py`

`calc_hbond.py` performs **local chemical validation** using CpH lambda values. It does not replace the coupled-state analysis layer.

The responsibilities are:

```text
calc_hbond.py
    decides whether a geometric H-bond is chemically possible in each frame

partition_by_state.py
    conditions the resulting scalar observable on the desired protonation/
    coupled-state classification and computes state-specific statistics
```

For example, after generating a targeted observable such as `N46_to_D52`, select a scalar metric such as:

```text
hbond_present
```

or:

```text
hbond_count
```

and partition it using the authoritative state-assignment table.

The frame-level H-bond output already contains:

```text
trajectory
trajectory_frame
step
```

which support the exact joining contract used by the structural-analysis toolkit.

`partition_by_state.py` should still validate `step` independently and should form statistical blocks on the full original trajectory timeline before state selection.

### Which H-bond metric should normally be partitioned?

For a defined residue-pair H-bond, `hbond_present` is often the cleanest occupancy observable.

`hbond_count` is useful when simultaneous multiple triplets are mechanistically meaningful.

For production state-conditioned analysis, a conservative default is to restrict summary statistics to:

```text
chemistry_status == clean
```

if the downstream partitioning script does not explicitly understand the uncertainty-bound metrics.

Alternatively, `hbond_count_min` and `hbond_count_max` provide lower/upper sensitivity bounds for frames with unresolved local chemistry.

---

## 23. Validated v5.1 behavior

The production candidate was validated on a 5,150-frame ff19SB HEWL pH 6.5 trajectory with:

```text
--scan PROTEIN '^1'
--focus DYAD ':36,53' '^1'
```

Expected accounting for two observables is:

```text
5150 frames
x 2 observables
x 6 metrics
= 61800 frame-wise rows
```

The validation run produced:

```text
Trajectories analyzed       : 1
Requested observables       : 2
Frame-wise output rows      : 61800
Atom-level report rows      : 878
Program version             : 5.1
Chemistry profile           : protein
Candidate engine            : donor-enumerated-v5.1
CpH filtering               : enabled
Partial-ambiguity frame obs : 1240
Outside-lambda frame obs    : 93
Distance cutoff (A)         : 3
Angle cutoff (deg)          : 135
```

Status accounting closed exactly for every metric.

For the DYAD observable:

```text
clean                  4716
partial_ambiguity       391
outside_lambda_window    43
---------------------------
total                  5150
```

For the whole-protein observable:

```text
clean                  4251
partial_ambiguity       849
outside_lambda_window    50
---------------------------
total                  5150
```

The sums across observables match the run summary:

```text
391 + 849 = 1240 partial-ambiguity frame-observations
43  + 50  =   93 outside-lambda frame-observations
```

### Candidate-generation validation

The broad DYAD focus recovered the targeted N46/N59 -> D52 H-bond candidates that were previously detected only by explicit directed tests.

Representative chemically valid triplet counts at pH 6.5 included:

```text
N59 ND2-HD22 -> D52 OD1 : 1387 frames
N46 ND2-HD21 -> D52 OD1 : 1073 frames
N59 ND2-HD22 -> D52 OD2 : 1013 frames
N46 ND2-HD21 -> D52 OD2 :  848 frames
```

This validates the donor-enumerated broad-focus engine against the corresponding targeted interactions.

### Backbone chemistry validation

Backbone H-bonds on titratable E35/D52 were verified to remain independent of side-chain lambda chemistry. Representative interactions retained identical raw and chemically valid counts, confirming that only the true titrating atoms are CpH-dependent.

---

## 24. Recommended production validation checklist

Before launching a new system or full pH series, verify:

1. The program reports `Program version : 5.1`.
2. The candidate engine reports `donor-enumerated-v5.1`.
3. Every trajectory is paired with the intended lambda file in `<output>.settings.txt`.
4. Frame count matches the expected trajectory selection.
5. For every observable, each of the six metrics has the same number of frame rows.
6. `clean + partial_ambiguity + outside_lambda_window` equals the number of analyzed frames for each observable.
7. A few known protein H-bonds have chemically sensible donor and acceptor atom identities.
8. Known titratable-site interactions show sensible differences between `raw_count` and chemically filtered `count`.
9. Broad `--focus` results recover equivalent candidates found by one or two targeted `--directed` controls.
10. `hbond_count_min <= hbond_count_max` for every frame, with equality in clean frames.

A useful accounting command is:

```bash
awk -F $'\t' '
NR==1 {
    for (i=1;i<=NF;i++) h[$i]=i
    next
}
{
    key=$h["observable"] "\t" $h["metric"] "\t" $h["chemistry_status"]
    n[key]++
}
END {
    for (k in n) print k "\t" n[k]
}' output.tsv | sort
```

---

## 25. Common mistakes

### Treating raw CPPTRAJ candidates as CpH-valid H-bonds

Do not use `raw_count` as the CpHMD occupancy. Use `count` in the interaction report or the frame-wise `hbond_count`/`hbond_present` metrics.

### Summing atom-triplet occupancy fractions

Do not sum fractions from different donor-H/acceptor triplets to obtain residue-pair occupancy. Use a frame-wise targeted observable.

### Ignoring the tautomer coordinate

For Asp/Glu and His, protonation state alone is insufficient to assign all atom-specific donor/acceptor roles in the protonated/neutral tautomeric states.

### Making all atoms in a titratable residue lambda-dependent

Only the titrating atoms change donor/acceptor role. Backbone chemistry remains fixed.

### Using nearest lambda rows

The script intentionally forbids nearest-neighbor synchronization. Coordinate and lambda trajectories must be mapped by exact MD step.

### Reusing an incompatible work directory

Use a new `--workdir` for changed selections/settings or deliberately invoke `--force` after confirming that overwriting is appropriate.

### Forgetting to source Amber

If `cpptraj` is not found, source the Amber environment or use `--cpptraj /full/path/cpptraj`.

---

## 26. Command-line reference

Current production options:

```text
-p, --topology TOPOLOGY
-t, --trajectory TRAJECTORY
--lambda LAMBDA
--no-cph-filter
--lambda-low FLOAT
--lambda-high FLOAT

--interaction NAME MASK1 MASK2
--directed NAME DONOR_MASK ACCEPTOR_MASK
--directed-h NAME DONOR_MASK DONOR_H_MASK ACCEPTOR_MASK
--scan NAME MASK
--focus NAME FOCUS_MASK PARTNER_MASK
--solvent-site NAME SOLUTE_MASK
--solvent-donor-mask MASK
--solvent-acceptor-mask MASK

--chemistry-profile {protein,fon}
--include-met-sulfur
--distance FLOAT
--angle FLOAT
--image-anchor MASK
--hbond-image

--start INT
--stop INT|last
--stride INT
--first-frame-step INT
--frame-step-interval INT

-j, --jobs, --fork INT
--cpptraj PATH
--workdir DIR
--detail-output FILE
-o, --output FILE
--force
```

Use:

```bash
python3 calc_hbond.py --help
```

for the executable's authoritative CLI summary.

---

## 27. Methods-language guidance

A concise methods description for CpHMD use is:

> Hydrogen bonds were identified using a 3.0 A donor-heavy-atom/acceptor distance cutoff and a 135 degree angular cutoff. CPPTRAJ was used to generate frame-resolved geometric donor-H-acceptor candidates. Candidate donor/acceptor assignments involving titratable atoms were subsequently filtered using the exactly matched CpHMD protonation and tautomer coordinates for each trajectory frame. Asp/Glu and His tautomer-specific donor/acceptor identities were resolved from the corresponding Amber continuous-CpH coordinates, whereas backbone and other non-titrating atoms retained fixed protein chemistry. Frames containing intermediate lambda/tautomer values were retained with explicit lower/upper H-bond-count bounds rather than being discarded.

For atom-level geometry averages, state explicitly that the reported CPPTRAJ distance and angle averages are over geometrically present candidates before CpH filtering unless a separate filtered-geometry analysis has been performed.

---

## 28. Scope

`calc_hbond.py` is intended to be a general CpHMD hydrogen-bond observable generator.

It is **not**:

- hardcoded to HEWL;
- hardcoded to the E35/D52 dyad;
- restricted to Asp/Glu;
- a coupled-state classifier;
- a hydration-number calculator;
- a hydrogen-bond network/community-analysis package.

Its responsibility is narrower and explicit:

> generate auditable frame-resolved hydrogen-bond observables whose donor/acceptor chemistry is consistent with the instantaneous CpHMD state.

Coupled-state conditioning belongs downstream in `partition_by_state.py`, and hydration-shell observables should be generated by a dedicated hydration tool.

---

## 29. Production status

The version documented here is:

```text
calc_hbond.py
Program version: 5.1
Candidate engine: donor-enumerated-v5.1
```

