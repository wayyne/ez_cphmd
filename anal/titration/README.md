# Titration analysis

`calc_titration.py` analyzes Amber CpHMD lambda trajectories and supports ordinary single-site, histidine, and coupled two-site titration analysis.

Input pH values are read from filenames using Amber/AREX-style names such as `arex.ph5_0.lambda`.

## Single-site titration

```bash
python3 calc_titration.py \
    -t 3:43 \
    -R 35 \
    -o E35_titration.png \
    arex.ph*.lambda
```

The time window is given in nanoseconds as `START:END`.

By default, lambda values at or below `0.2` are treated as protonated and values at or above `0.8` as deprotonated. Frames between the two cutoffs are tracked as mixed rather than assigned to either clean state.

Ordinary mode writes `titcurve.dat` and reports fitted pKa and Hill values. Optional block/convergence tables can be requested with `--write-stats-csv` and `--write-convergence-csv`.

## Histidine

Use `-H/--his-resid` for histidine analysis:

```bash
python3 calc_titration.py \
    -t 3:43 \
    -H 15 \
    --his-map hid-low \
    -o H15_titration.png \
    arex.ph*.lambda
```

The script uses the `itauto` header to distinguish the histidine protonation and tautomer coordinates.

## Coupled two-site titration

Use `-R` for the primary residue and `--coupled` for the second residue:

```bash
python3 calc_titration.py \
    -t 3:43 \
    -R 35 \
    --coupled 52 \
    --primary-label E35 \
    --coupled-label D52 \
    --run-label ff19 \
    --stats both \
    arex.ph*.lambda
```

Coupled mode classifies the four clean joint protonation states plus a mixed state, then performs several complementary analyses:

1. pair-level macroscopic stepwise pKa fitting
2. residue-specific marginal and conditional titration curves
3. independent and interacting four-state thermodynamic fits
4. optional block, bootstrap, convergence, and cutoff-sensitivity diagnostics

Macroscopic `pK1` and `pK2` are properties of the two-site pair. Residue-specific assignments are reported separately by the four-state model.

All coupled-mode outputs are routed into `coupled_anal/` by default. A typical run produces files such as:

```text
ff19_E35_D52_settings.json
ff19_E35_D52_frame_states.csv
ff19_E35_D52_joint_states.csv
ff19_E35_D52_marginal_and_conditional_curves.csv
ff19_E35_D52_macroscopic_fit.csv
ff19_E35_D52_four_state_model_comparison.csv
ff19_E35_D52_summary.csv
ff19_E35_D52_summary.json
ff19_E35_D52_report.txt
```

Depending on the selected statistics and convergence options, block-fit, bootstrap, convergence, cutoff-sensitivity, and plot files are also generated.

The frame-state CSV can be passed to `../arex/plot_arex_analysis.py` for exact MD-step state-conditioned structural analysis.

## Lambda-coordinate selection

Amber CpHMD lambda files may contain more than one coordinate for a titratable residue. The script uses the `ires` and `itauto` headers to identify protonation coordinates rather than silently selecting the first matching column. Ambiguous duplicated residue mappings cause that file to be rejected rather than guessed.

## Fitting and dependencies

The default ordinary/conditional fit is nonlinear Henderson-Hasselbalch fitting with a free Hill coefficient. `--hill fixed1` fixes Hill to 1, and `--fit-method logit` is also available.

- Python 3 is required.
- SciPy is used when available for numerical optimization; built-in fallback searches are provided when it is unavailable.
- Matplotlib is needed for plots.
- Coupled bootstrap analysis can use multiple CPU processes with `--jobs`.

Run the full option list with:

```bash
python3 calc_titration.py --help
```
