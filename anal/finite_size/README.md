# Finite-size diagnostics

These scripts report box-volume and water-number-density quantities useful when checking or comparing simulation box sizes.

They do not, by themselves, establish a finite-size effect. A finite-size comparison requires comparing the resulting observables across systems prepared with different box sizes or other relevant conditions.

## Single trajectory

```bash
./check_finite_size_effect.sh system.prmtop trajectory.nc
```

The script reports:

- number of analyzed frames
- mean box volume and volume SD
- number of `WAT` residues in the topology
- mean framewise water number density, `<Nwater/V>`
- framewise density SD
- `Nwater/<V>`

Temporary CPPTRAJ files are removed automatically.

## Pool multiple trajectories

```bash
./pool_finite_size_effect.sh \
    system.prmtop \
    'arex.ph*.nc' \
    100
```

The optional third argument is the number of initial frames skipped independently from each trajectory.

The script writes:

```text
finite_size_per_file.tsv
```

and prints an overall pooled result to standard output.

The pooled result is **frame-weighted**: every retained frame enters the pooled statistics once. If trajectories contribute different numbers of retained frames, trajectories with more frames therefore contribute more strongly.

## Assumptions

- all pooled trajectories use the same topology
- water residues are named `WAT`
- topology water count is fixed during the analyzed trajectories
- the supplied trajectory glob should be quoted so it is expanded inside the script

## Requirements

- AmberTools / `cpptraj`
- Bash
- `awk`
