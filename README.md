# ez_cphmd

Tools for streamlined constant-pH molecular dynamics (CpHMD) system preparation and trajectory analysis.

The goal of `ez_cphmd` is to make common CpHMD workflows easier to set up and reproduce while keeping the underlying CHARMM and Amber preparation steps accessible.

## Features

### System preparation

Prepare CpHMD systems from either:

- an amino-acid sequence, or
- an existing protein structure / RCSB PDB ID.

Supported force fields:

- CHARMM22 CpHMD (`c22`)
- Amber `ff14SB`
- Amber `ff19SB`

Preparation includes:

- capped or uncapped termini
- CpHMD residue preparation
- PDB cleaning and disulfide handling
- solvation
- salt concentration or explicit ion counts
- cubic and truncated-octahedral boxes
- rhombic boxes for the CHARMM22 workflow
- final Amber `prmtop` / `inpcrd` generation

### Analysis

Analysis utilities are available for:

| Directory | Purpose |
| --- | --- |
| `anal/arex/` | AREX exchange diagnostics and residue-centric multi-pH analysis |
| `anal/titration/` | CpHMD lambda trajectories, titration curves, pKa, and coupled-site analysis |
| `anal/finite_size/` | Box-volume and water-density diagnostics for finite-size comparisons |
| `anal/rmsd/` | Trajectory RMSD with separate fit and measurement masks |
| `anal/rmsf/` | Per-residue backbone, side-chain, or total heavy-atom RMSF |
| `anal/dihedral/` | Dihedral-angle analysis |
| `anal/distance/` | Interatomic distance analysis |
| `anal/hbond/` | Hydrogen-bond analysis |
| `anal/hydration/` | Hydration analysis |
| `anal/partition/` | Partitioning trajectories/data by state |
| `anal/sasa/` | Solvent-accessible surface area analysis |
| `anal/water_wire/` | Water proximity and water-wire analysis |

Each analysis directory is intended to remain usable on its own, with a local `README.md` for workflow-specific details.

## Repository layout

```text
ez_cphmd/
├── prep/
│   ├── config.example.sh
│   │
│   ├── seq/
│   │   ├── prep_seq.sh
│   │   ├── c22/
│   │   └── amber/
│   │
│   ├── struct/
│   │   ├── prep_struct.sh
│   │   ├── get_pdb.sh
│   │   ├── c22/
│   │   └── amber/
│   │
│   ├── c22/
│   └── amber/
│
└── anal/
    ├── arex/
    ├── titration/
    ├── finite_size/
    ├── rmsd/
    ├── rmsf/
    ├── dihedral/
    ├── distance/
    ├── hbond/
    ├── hydration/
    ├── partition/
    ├── sasa/
    └── water_wire/
```

## Installation

Clone the repository:

```bash
git clone https://github.com/wayyne/ez_cphmd.git
cd ez_cphmd
```

There is no package installation step. The preparation workflows are driven by Bash scripts and use existing CHARMM / AmberTools installations.

## Configuration

Copy the example configuration:

```bash
cp prep/config.example.sh prep/config.sh
```

Then edit `prep/config.sh` for your local installation:

```bash
AMBER_ENV="/path/to/amber.sh"

C22_TOPPAR="/path/to/c22/toppar"

FF14_PHMD_DIR="/path/to/ff14/parmfiles"
FF19_PHMD_DIR="/path/to/ff19/parmfiles"
```

`prep/config.sh` contains machine-specific paths and is not intended to be committed.

Runtime settings such as sequence, force field, box size, salt concentration, and capping are set in the preparation entry scripts rather than in `config.sh`.

## Preparing a system from sequence

Edit the configuration block near the top of:

```text
prep/seq/prep_seq.sh
```

For example:

```bash
ff="ff19sb"

pep="GDG"
capped="1"

cushion="15"
boxtype="octahedral"

sod="1"
cla="0"
conc=""

titr="0"
```

Available force-field selections are:

```text
c22
ff14sb
ff19sb
```

Run the script from the directory where the system files should be created:

```bash
mkdir my_system
cd my_system

bash /path/to/ez_cphmd/prep/seq/prep_seq.sh
```

The general workflow is:

```text
sequence
   ↓
dry CpHMD system
   ↓
solvation
   ↓
ions
   ↓
prmtop / inpcrd / PDB
```

## Preparing a system from a protein structure

Edit:

```text
prep/struct/prep_struct.sh
```

For example:

```bash
ff="ff19sb"

pdbid="2lzt"
pdbfile=""
chain="A"

cushion="15"
boxtype="octahedral"
capped="1"

sod="0"
cla="0"
conc="150"

titr="0"
```

Set exactly one of `pdbid` or `pdbfile`.

To download a structure from the RCSB PDB:

```bash
pdbid="2lzt"
pdbfile=""
```

To use a local structure:

```bash
pdbid=""
pdbfile="/path/to/2lzt_original.pdb"
```

Then run:

```bash
mkdir my_system
cd my_system

bash /path/to/ez_cphmd/prep/struct/prep_struct.sh
```

The structure workflow performs:

```text
PDB input
   ↓
chain extraction and validation
   ↓
CpHMD residue preparation
   ↓
disulfide handling
   ↓
dry system construction
   ↓
restrained initial minimization
   ↓
solvation and ions
   ↓
prmtop / inpcrd / PDB
```

## Force-field details

### CHARMM22

The `c22` workflow uses the CHARMM22 CpHMD topology and parameter set.

The dry CHARMM system is solvated in CHARMM and converted to Amber topology/coordinate format using ParmEd.

Supported box geometries include:

- cubic
- octahedral
- rhombic

### Amber ff14SB

The `ff14sb` workflow uses:

- `leaprc.protein.ff14SB`
- TIP3P water
- CpHMD parameters for ff14SB
- `ff14sb_pme_dayhoff_dev.parm`

### Amber ff19SB

The `ff19sb` workflow uses:

- `leaprc.protein.ff19SB`
- OPC water
- CpHMD parameters for ff19SB

## Ion setup

Two ion modes are supported.

For a target salt concentration:

```bash
conc="150"
```

with concentration specified in mM.

For explicit ion counts:

```bash
conc=""
sod="10"
cla="12"
```

When `conc` is set, the explicit `sod` and `cla` settings are not used as the bulk-salt specification.

## Structure input scope

The current structure-preparation workflow is intentionally conservative and expects:

- one selected protein chain
- standard amino-acid residues
- zero or one structural model
- consecutive residue numbering
- no insertion codes
- complete backbone `N`, `CA`, `C`, and `O` atoms
- blank or `A` alternate locations
- same-chain disulfide bonds

Existing hydrogens and terminal `OXT` atoms are removed and rebuilt during preparation.

If a structure falls outside the supported scope, the preparation scripts generally stop rather than silently guessing how the structure should be interpreted.

## Analysis tools

Each analysis utility is kept in a focused directory. See the local README for full options and workflow-specific assumptions.

### AREX

```text
anal/arex/analyze_arex.sh
anal/arex/plot_arex_analysis.py
anal/arex/arex_acceptance_rate.sh
anal/arex/roundtrip_check.py
```

Exchange diagnostics report neighboring-pH acceptance, replica slot occupancy, movement, and complete round trips. The structural AREX workflow analyzes each pH-slot trajectory independently and summarizes residue SASA, hydration shells, hydrogen bonds, RMSD, and radius of gyration.

See [`anal/arex/README.md`](anal/arex/README.md).

### Titration

```text
anal/titration/calc_titration.py
```

Analyzes Amber CpHMD lambda trajectories in single-site, histidine, or coupled two-site modes. Coupled mode includes joint protonation-state populations, macroscopic stepwise pKa fitting, residue-specific conditional fits, four-state thermodynamic modeling, and optional block/bootstrap/convergence diagnostics.

See [`anal/titration/README.md`](anal/titration/README.md).

### Finite-size diagnostics

```text
anal/finite_size/check_finite_size_effect.sh
anal/finite_size/pool_finite_size_effect.sh
```

Reports box-volume and water-number-density statistics for one trajectory or a pooled trajectory set. These quantities are intended to support finite-size comparisons across systems rather than serve as a finite-size test by themselves.

See [`anal/finite_size/README.md`](anal/finite_size/README.md).

### RMSD

```text
anal/rmsd/calc_rmsd.py
```

Fits trajectories using one CPPTRAJ mask and can measure RMSD with a different mask on the aligned coordinates.

See [`anal/rmsd/README.md`](anal/rmsd/README.md).

### RMSF

```text
anal/rmsf/calc_rmsf.py
```

Calculates per-residue backbone, side-chain, or total heavy-atom RMSF after backbone fitting.

See [`anal/rmsf/README.md`](anal/rmsf/README.md).

### Dihedral

```text
anal/dihedral/calc_dihedral.py
```

See [`anal/dihedral/README.md`](anal/dihedral/README.md).

### Distance

```text
anal/distance/calc_distance.py
```

See [`anal/distance/README.md`](anal/distance/README.md).

### Hydrogen bonds

```text
anal/hbond/calc_hbond.py
```

See [`anal/hbond/README.md`](anal/hbond/README.md).

### Hydration

```text
anal/hydration/calc_hydration.py
```

See [`anal/hydration/README.md`](anal/hydration/README.md).

### State partitioning

```text
anal/partition/partition_by_state.py
```

See [`anal/partition/README.md`](anal/partition/README.md).

### SASA

```text
anal/sasa/calc_sasa.py
anal/sasa/calc_fsasa.py
anal/sasa/calc_refsasa.py
```

See [`anal/sasa/README.md`](anal/sasa/README.md).

### Water wires

```text
anal/water_wire/calc_water_prox.py
anal/water_wire/calc_water_wire.py
```

The directory also contains smoke tests and reference outputs for the supported force-field workflows.

See [`anal/water_wire/README.md`](anal/water_wire/README.md).

## Requirements

Exact requirements depend on the workflow being used.

Preparation may require:

- CHARMM
- AmberTools
- ParmEd
- `tleap`
- `sander`
- `cpptraj`
- Python 3

Some analysis scripts additionally use NumPy and Matplotlib. `calc_titration.py` can use SciPy for numerical optimization but includes fallback fitting paths when SciPy is unavailable. See the individual analysis READMEs for tool-specific requirements.

## Output and reproducibility

Preparation scripts retain logs and intermediate reports so that generated systems can be inspected.

Running each system in its own working directory is recommended.

Before starting production simulations, users should inspect the generated structure, topology, ion composition, CpHMD residue setup, and preparation logs.

## Citation

If these tools are used in published work, please cite the relevant CpHMD methodology, molecular dynamics software, and force field used in the simulations.

Repository-specific citation information can be added here as appropriate.
