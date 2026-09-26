#!/usr/bin/env python3
"""
calc_fsasa.py

Calculate frame-wise fractional SASA (fSASA) by exact-joining raw SASA to a
state table and normalizing each frame by the chemistry/state-specific reference
SASA denominator.

Two state-classification modes are supported.

1. ``--lambda-mode single`` (default; backward compatible)

   For ordinary one-coordinate titratable groups such as Asp/Glu:

       lambda <= --low   -> H
       lambda >= --high  -> deprot
       otherwise         -> mixed/excluded

   Reference rows are therefore keyed by H/deprot exactly as in the original
   workflow.

2. ``--lambda-mode his``

   For histidine, the joined state table must contain both the protonation
   coordinate and neutral-tautomer coordinate. The physical state is assigned
   with the same rules used by the H15 analysis and the His-aware
   ``calc_refsasa.py``:

       lambda <= --low   -> HIP, independent of x
       lambda >= --high and x <= --low  -> --x-low-tautomer
       lambda >= --high and x >= --high -> --x-high-tautomer

   Intermediate lambda values are excluded. Neutral frames with intermediate x
   are also excluded from fSASA normalization because no pure HID/HIE reference
   denominator applies. Reference rows are keyed separately as HIP, HID, HIE.

The raw SASA table remains state-agnostic. All joins are exact on the requested
frame identity. Optional state-column mappings make the tool compatible with
both the earlier dyad frame tables and the H15 ``analyze_his.py`` output.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import math
from pathlib import Path


OUT_FIELDS = [
    "trajectory", "pH", "analysis_row", "trajectory_frame", "step",
    "observable", "metric", "selection", "solutemask", "value", "unit",
]

AUDIT_FIELDS = [
    "trajectory", "pH", "analysis_row", "trajectory_frame", "step",
    "observable", "raw_metric", "selection", "solutemask", "raw_sasa_A2",
    "chemistry", "site_definition", "lambda_mode", "lambda_column",
    "local_lambda", "tautomer_column", "local_tautomer_coordinate",
    "local_protonation", "local_tautomer", "reference_state",
    "state_clean", "state_matched", "state_code", "state_label",
    "source_state", "source_tautomer", "reference_forcefield",
    "reference_sasa_A2", "reference_method", "fsasa", "fsasa_status",
]


def open_text(path: Path, mode: str):
    """Open plain-text or .gz text files transparently."""
    if "b" in mode:
        raise ValueError("open_text only supports text modes")
    if path.suffix.lower() == ".gz":
        return gzip.open(path, mode, newline="")
    return path.open(mode, newline="")


def read_tsv(path: Path):
    if not path.is_file():
        raise FileNotFoundError(path)
    with open_text(path, "rt") as fh:
        r = csv.DictReader(fh, delimiter="\t")
        if r.fieldnames is None:
            raise RuntimeError(f"TSV has no header: {path}")
        return list(r.fieldnames), list(r)


def write_tsv(path: Path, rows, fields):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open_text(path, "wt") as fh:
        w = csv.DictWriter(
            fh, fieldnames=fields, delimiter="\t", extrasaction="ignore"
        )
        w.writeheader()
        for row in rows:
            out = dict(row)
            for key, value in list(out.items()):
                if isinstance(value, float):
                    out[key] = "" if not math.isfinite(value) else f"{value:.10g}"
            w.writerow(out)


def nonempty(value):
    return value is not None and bool(str(value).strip())


def truthy(value):
    return str(value).strip().lower() in {"1", "true", "yes", "y"}


def number(value, context):
    if not nonempty(value):
        raise RuntimeError(f"Missing numeric {context}")
    try:
        x = float(str(value).strip())
    except ValueError as exc:
        raise RuntimeError(f"Non-numeric {context}: {value!r}") from exc
    if not math.isfinite(x):
        raise RuntimeError(f"Non-finite {context}: {value!r}")
    return x


def integer(value, context):
    x = number(value, context)
    i = int(round(x))
    if abs(x - i) > 1.0e-8:
        raise RuntimeError(f"Non-integral {context}: {value!r}")
    return i


def canonical_ph(row):
    vals = []
    for key in ("pH", "ph"):
        if nonempty(row.get(key)):
            vals.append(number(row[key], key))
    if not vals:
        return ""
    if any(
        not math.isclose(vals[0], x, rel_tol=0.0, abs_tol=1.0e-8)
        for x in vals[1:]
    ):
        raise RuntimeError(f"Conflicting pH/ph values in row: {row}")
    return f"{vals[0]:.10g}"


def trajectory_key(value: str, mode: str) -> str:
    text = str(value).strip()
    if mode == "exact":
        return text
    if mode == "stem":
        # Do not use Path.stem blindly: labels such as ``arex.ph5_0`` contain
        # a dot but have no filename extension. Strip only common trajectory
        # extensions from the basename.
        name = Path(text).name
        lower = name.lower()
        for ext in (".nc", ".netcdf", ".mdcrd", ".crd", ".dcd", ".xtc", ".trr"):
            if lower.endswith(ext):
                return name[:-len(ext)]
        return name
    raise ValueError(f"Unknown trajectory key mode: {mode}")


def canonical_value(row, canonical_col, *, source, args, context):
    """Return a canonical join/validation field from observation or state row."""
    if source == "observation":
        actual = canonical_col
    elif source == "state":
        mapping = {
            "trajectory": args.state_trajectory_column,
            "analysis_row": args.state_analysis_row_column,
            "trajectory_frame": args.state_trajectory_frame_column,
            "step": args.state_step_column,
        }
        actual = mapping.get(canonical_col, canonical_col)
    else:
        raise ValueError(source)

    value = row.get(actual)
    if not nonempty(value):
        raise RuntimeError(
            f"Missing {source} column/value {actual!r} for canonical "
            f"{canonical_col!r} in {context}"
        )

    if canonical_col == "trajectory":
        return trajectory_key(str(value), args.trajectory_key_mode)
    if canonical_col in {"analysis_row", "trajectory_frame", "step"}:
        return str(integer(value, f"{context} {actual}"))
    return str(value).strip()


def make_key(row, columns, *, source, args, context):
    return tuple(
        canonical_value(
            row, col, source=source, args=args, context=context
        )
        for col in columns
    )


def optional_int(row, column, context):
    if column is None or not nonempty(row.get(column)):
        return None
    return integer(row[column], f"{context} {column}")


def validate_same_frame(obs, state, args, context):
    comparisons = (
        ("analysis_row", args.state_analysis_row_column),
        ("trajectory_frame", args.state_trajectory_frame_column),
        ("step", args.state_step_column),
    )
    for obs_col, state_col in comparisons:
        if not nonempty(obs.get(obs_col)) or not state_col or not nonempty(state.get(state_col)):
            continue
        a = integer(obs[obs_col], f"observation {obs_col} ({context})")
        b = integer(state[state_col], f"state {state_col} ({context})")
        if a != b:
            raise RuntimeError(
                f"Frame metadata mismatch for {context}: "
                f"{obs_col}={a}, state {state_col}={b}"
            )

    op = canonical_ph(obs)
    sp = canonical_ph(state)
    if op and sp and not math.isclose(float(op), float(sp), rel_tol=0.0, abs_tol=1.0e-8):
        raise RuntimeError(
            f"pH mismatch for {context}: observation={op}, state={sp}"
        )


def ref_key(forcefield, chemistry, site_definition, protonation):
    return (
        forcefield.strip(),
        chemistry.strip().upper(),
        site_definition.strip(),
        protonation.strip(),
    )


def load_refs(path):
    fields, rows = read_tsv(path)
    required = {
        "forcefield", "chemistry", "site_definition",
        "protonation", "reference_sasa_A2",
    }
    missing = sorted(required - set(fields))
    if missing:
        raise RuntimeError(
            f"Reference table {path} missing: {', '.join(missing)}"
        )

    lookup = {}
    for lineno, row in enumerate(rows, 2):
        key = ref_key(
            row["forcefield"], row["chemistry"],
            row["site_definition"], row["protonation"]
        )
        den = number(
            row["reference_sasa_A2"],
            f"reference_sasa_A2 at {path}:{lineno}",
        )
        if den <= 0:
            raise RuntimeError(
                f"Reference denominator must be > 0 at {path}:{lineno}"
            )
        if key in lookup:
            raise RuntimeError(
                f"Duplicate reference key at {path}:{lineno}: {key}"
            )
        lookup[key] = row
    return lookup


def build_sites(site_args):
    out = {}
    for obs, chemistry, site_definition, lambda_column in site_args:
        if obs in out:
            raise ValueError(f"Duplicate --site mapping for {obs!r}")
        out[obs] = {
            "chemistry": chemistry.strip().upper(),
            "site_definition": site_definition.strip(),
            "lambda_column": lambda_column.strip(),
        }
    return out


def classify_lambda(lam, low, high):
    if lam <= low:
        return "H"
    if lam >= high:
        return "deprot"
    return "mixed"


def classify_his(lam, x, low, high, x_low_tautomer, x_high_tautomer):
    """Return (reference_state, physical_tautomer, status_class)."""
    if lam <= low:
        return "HIP", "HIP", "clean"
    if lam < high:
        return "mixed", "mixed", "lambda_mixed"

    if x <= low:
        return x_low_tautomer, x_low_tautomer, "clean"
    if x >= high:
        return x_high_tautomer, x_high_tautomer, "clean"
    return "neutral_tautomer_mixed", "mixed", "tautomer_mixed"


def required_reference_states(args, chemistry):
    if args.lambda_mode == "his":
        if chemistry != "HIS":
            raise ValueError(
                "--lambda-mode his currently requires every configured site to "
                "use chemistry HIS"
            )
        return ("HIP", "HID", "HIE")
    if chemistry == "HIS":
        raise ValueError(
            "HIS sites require --lambda-mode his so HID and HIE are not pooled"
        )
    return ("H", "deprot")


def expected_his_source_labels(reference_state):
    if reference_state == "HIP":
        return "protonated", "HIP"
    if reference_state in {"HID", "HIE"}:
        return "neutral", reference_state
    if reference_state == "mixed":
        return "mixed", "mixed"
    if reference_state == "neutral_tautomer_mixed":
        return "neutral", "mixed"
    return None, None


def validate_his_source_labels(state, reference_state, args, context):
    if not args.validate_his_state_labels:
        return

    expected_state, expected_taut = expected_his_source_labels(reference_state)

    if args.state_class_column in state and nonempty(state.get(args.state_class_column)):
        observed = str(state[args.state_class_column]).strip()
        if expected_state is not None and observed != expected_state:
            raise RuntimeError(
                f"His state-label mismatch for {context}: classified "
                f"{reference_state}, but {args.state_class_column}={observed!r}; "
                f"expected {expected_state!r}"
            )

    if args.tautomer_label_column in state and nonempty(state.get(args.tautomer_label_column)):
        observed = str(state[args.tautomer_label_column]).strip()
        if expected_taut is not None and observed != expected_taut:
            raise RuntimeError(
                f"His tautomer-label mismatch for {context}: classified "
                f"{reference_state}, but {args.tautomer_label_column}={observed!r}; "
                f"expected {expected_taut!r}"
            )


def parse_args():
    ap = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=(
            "Calculate frame-wise fSASA by exact-joining raw SASA to local "
            "state coordinates and chemistry/state-specific references."
        ),
    )
    ap.add_argument("--sasa", type=Path, required=True)
    ap.add_argument("--states", type=Path, required=True)
    ap.add_argument("--references", type=Path, required=True)
    ap.add_argument("--forcefield", required=True)
    ap.add_argument(
        "--site", action="append", nargs=4, required=True,
        metavar=("OBSERVABLE", "CHEMISTRY", "SITE_DEFINITION", "LAMBDA_COLUMN"),
        help=(
            "Map a raw SASA observable to its chemistry/site definition and "
            "state-table protonation coordinate. Example: "
            "--site E35_COO GLU carboxylate_oxygens primary_lambda. "
            "For His, LAMBDA_COLUMN is the protonation lambda column; the "
            "tautomer coordinate is supplied globally with --tautomer-column."
        ),
    )
    ap.add_argument(
        "--lambda-mode", choices=("single", "his"), default="single",
        help=(
            "single: legacy H/deprot normalization; "
            "his: HIP/HID/HIE normalization using lambda plus tautomer x."
        ),
    )
    ap.add_argument(
        "--tautomer-column", default="x",
        help="State-table tautomer-coordinate column used in --lambda-mode his.",
    )
    ap.add_argument(
        "--x-low-tautomer", choices=("HID", "HIE"), default=None,
        help="Physical tautomer at x <= --low in His mode.",
    )
    ap.add_argument(
        "--x-high-tautomer", choices=("HID", "HIE"), default=None,
        help="Physical tautomer at x >= --high in His mode.",
    )

    ap.add_argument(
        "--join-on", default="trajectory,trajectory_frame",
        help=(
            "Canonical raw-SASA identity columns used for the exact join. "
            "State-side column names may be remapped with the --state-*-column options."
        ),
    )
    ap.add_argument(
        "--trajectory-key-mode", choices=("exact", "stem"), default="exact",
        help=(
            "exact requires identical trajectory labels; stem compares Path(...).stem "
            "and is useful when one table stores 'arex.ph5_0' and the other "
            "'arex.ph5_0.nc'."
        ),
    )
    ap.add_argument("--state-trajectory-column", default="trajectory")
    ap.add_argument("--state-analysis-row-column", default="analysis_row")
    ap.add_argument("--state-trajectory-frame-column", default="trajectory_frame")
    ap.add_argument("--state-step-column", default="step")

    ap.add_argument("--low", type=float, default=0.2)
    ap.add_argument("--high", type=float, default=0.8)
    ap.add_argument("--clean-column", default="clean")
    ap.add_argument("--matched-column", default="matched")
    ap.add_argument(
        "--state-class-column", default="state",
        help="Optional source state label used to validate His classification if present.",
    )
    ap.add_argument(
        "--tautomer-label-column", default="tautomer",
        help="Optional source tautomer label used to validate His classification if present.",
    )
    ap.add_argument(
        "--no-validate-his-state-labels", dest="validate_his_state_labels",
        action="store_false",
        help="Do not cross-check His lambda/x classification against state/tautomer labels.",
    )
    ap.set_defaults(validate_his_state_labels=True)

    ap.add_argument("--sasa-metric", default="sasa_lcpo")
    ap.add_argument("--audit-output", type=Path)
    ap.add_argument("-o", "--output", type=Path, default=Path("fsasa_frames.tsv"))
    return ap.parse_args()


def validate_args(args):
    if not (0.0 <= args.low < args.high <= 1.0):
        raise ValueError("Require 0 <= --low < --high <= 1")

    if args.lambda_mode == "his":
        if args.x_low_tautomer is None or args.x_high_tautomer is None:
            raise ValueError(
                "--lambda-mode his requires --x-low-tautomer and --x-high-tautomer"
            )
        if args.x_low_tautomer == args.x_high_tautomer:
            raise ValueError(
                "His x endpoints must map to different tautomers (HID and HIE)"
            )

    for name in (
        "state_trajectory_column", "state_analysis_row_column",
        "state_trajectory_frame_column", "state_step_column",
    ):
        if not str(getattr(args, name)).strip():
            raise ValueError(f"--{name.replace('_', '-')} must not be empty")


def main():
    args = parse_args()
    validate_args(args)

    join_cols = [x.strip() for x in args.join_on.split(",") if x.strip()]
    if not join_cols:
        raise ValueError("--join-on must contain at least one column")

    sites = build_sites(args.site)
    refs = load_refs(args.references)
    sasa_fields, sasa_rows = read_tsv(args.sasa)
    state_fields, state_rows = read_tsv(args.states)

    required_sasa = {
        "trajectory", "analysis_row", "trajectory_frame", "step",
        "observable", "metric", "selection", "solutemask", "value", "unit",
    }
    missing = sorted(required_sasa - set(sasa_fields))
    if missing:
        raise RuntimeError(
            f"SASA table {args.sasa} missing: {', '.join(missing)}"
        )

    state_join_mapping = {
        "trajectory": args.state_trajectory_column,
        "analysis_row": args.state_analysis_row_column,
        "trajectory_frame": args.state_trajectory_frame_column,
        "step": args.state_step_column,
    }

    for col in join_cols:
        if col not in sasa_fields:
            raise RuntimeError(f"Join column {col!r} absent from SASA table")
        state_col = state_join_mapping.get(col, col)
        if state_col not in state_fields:
            raise RuntimeError(
                f"Canonical join column {col!r} maps to state column "
                f"{state_col!r}, but that column is absent from {args.states}"
            )

    for obs, spec in sites.items():
        if spec["lambda_column"] not in state_fields:
            raise RuntimeError(
                f"{obs!r} requests state lambda column "
                f"{spec['lambda_column']!r}, but it is absent"
            )
        if args.lambda_mode == "his" and args.tautomer_column not in state_fields:
            raise RuntimeError(
                f"His mode requires tautomer coordinate column "
                f"{args.tautomer_column!r}, but it is absent"
            )

        for reference_state in required_reference_states(args, spec["chemistry"]):
            key = ref_key(
                args.forcefield, spec["chemistry"],
                spec["site_definition"], reference_state,
            )
            if key not in refs:
                raise RuntimeError(f"No reference row for {obs!r}: {key}")

    state_lookup = {}
    for lineno, row in enumerate(state_rows, 2):
        key = make_key(
            row, join_cols, source="state", args=args,
            context=f"state row {args.states}:{lineno}",
        )
        if key in state_lookup:
            raise RuntimeError(
                f"Duplicate state join key at {args.states}:{lineno}: {key}"
            )
        state_lookup[key] = row

    output = []
    audit = []
    seen = set()
    observables_seen = set()

    counts = {
        "input": 0,
        "joined": 0,
        "outside": 0,
        "local_mixed": 0,
        "tautomer_mixed": 0,
        "computed": 0,
        "global_unclean_local_clean": 0,
        "negative_raw": 0,
        "negative_fsasa": 0,
        "HIP": 0,
        "HID": 0,
        "HIE": 0,
        "H": 0,
        "deprot": 0,
    }

    for lineno, obs in enumerate(sasa_rows, 2):
        counts["input"] += 1
        observable = str(obs.get("observable", "")).strip()
        observables_seen.add(observable)

        if observable not in sites:
            raise RuntimeError(
                f"No --site mapping for observable {observable!r} "
                f"at {args.sasa}:{lineno}"
            )

        if (
            args.sasa_metric
            and str(obs.get("metric", "")).strip() != args.sasa_metric
        ):
            raise RuntimeError(
                f"Unexpected metric at {args.sasa}:{lineno}: "
                f"{obs.get('metric')!r}"
            )

        if str(obs.get("unit", "")).strip() not in {
            "A2", "A^2", "Angstrom^2", "angstrom^2",
        }:
            raise RuntimeError(
                f"Unexpected SASA unit at {args.sasa}:{lineno}: "
                f"{obs.get('unit')!r}"
            )

        raw = number(obs.get("value"), f"raw SASA at {args.sasa}:{lineno}")
        if raw < 0.0:
            counts["negative_raw"] += 1

        key = make_key(
            obs, join_cols, source="observation", args=args,
            context=f"SASA row {args.sasa}:{lineno}",
        )
        state = state_lookup.get(key)
        if state is None:
            raise RuntimeError(
                f"No exact state row for {args.sasa}:{lineno}; key={key}"
            )
        counts["joined"] += 1

        context = f"{args.sasa}:{lineno}, key={key}"
        validate_same_frame(obs, state, args, context)

        ident = key + (observable,)
        if ident in seen:
            raise RuntimeError(f"Duplicate SASA observation: {ident}")
        seen.add(ident)

        spec = sites[observable]
        lam_col = spec["lambda_column"]

        # If an explicit matched column exists, respect it. Otherwise infer
        # matched status from whether the local protonation coordinate exists.
        if args.matched_column in state_fields:
            matched = truthy(state.get(args.matched_column))
        else:
            matched = nonempty(state.get(lam_col))

        explicit_clean = None
        if args.clean_column in state_fields:
            explicit_clean = truthy(state.get(args.clean_column))

        lam = math.nan
        x = math.nan
        local_prot = ""
        local_taut = ""
        reference_state = ""
        den = math.nan
        ref_method = ""
        fsasa = math.nan

        if not matched:
            status = "outside_state_window"
            clean = False if explicit_clean is None else explicit_clean
            counts["outside"] += 1
        else:
            lam = number(state.get(lam_col), f"{lam_col} for {context}")
            if not 0.0 <= lam <= 1.0:
                raise RuntimeError(
                    f"{lam_col} outside [0,1] for {context}: {lam}"
                )

            if args.lambda_mode == "single":
                reference_state = classify_lambda(lam, args.low, args.high)
                local_prot = reference_state
                clean_local = reference_state != "mixed"

                if reference_state == "mixed":
                    status = "local_mixed"
                    counts["local_mixed"] += 1
                else:
                    rkey = ref_key(
                        args.forcefield, spec["chemistry"],
                        spec["site_definition"], reference_state,
                    )
                    rrow = refs[rkey]
                    den = number(
                        rrow["reference_sasa_A2"],
                        f"reference denominator {rkey}",
                    )
                    ref_method = str(rrow.get("method", "")).strip()
                    fsasa = raw / den
                    if not math.isfinite(fsasa):
                        raise RuntimeError(f"Non-finite fSASA for {context}")
                    status = "computed"
                    counts["computed"] += 1
                    counts[reference_state] += 1

            else:
                # HIP and lambda-intermediate frames do not chemically require
                # a tautomer coordinate. Only neutral lambda-end-state frames
                # need x to distinguish HID from HIE.
                if lam >= args.high:
                    x = number(
                        state.get(args.tautomer_column),
                        f"{args.tautomer_column} for {context}",
                    )
                    if not 0.0 <= x <= 1.0:
                        raise RuntimeError(
                            f"{args.tautomer_column} outside [0,1] for {context}: {x}"
                        )
                else:
                    x = math.nan

                reference_state, local_taut, cls = classify_his(
                    lam, x, args.low, args.high,
                    args.x_low_tautomer, args.x_high_tautomer,
                )
                local_prot = "HIP" if reference_state == "HIP" else (
                    "neutral" if reference_state in {"HID", "HIE", "neutral_tautomer_mixed"}
                    else "mixed"
                )
                clean_local = cls == "clean"

                validate_his_source_labels(
                    state, reference_state, args, context
                )

                if cls == "lambda_mixed":
                    status = "local_mixed"
                    counts["local_mixed"] += 1
                elif cls == "tautomer_mixed":
                    status = "neutral_tautomer_mixed"
                    counts["tautomer_mixed"] += 1
                else:
                    rkey = ref_key(
                        args.forcefield, spec["chemistry"],
                        spec["site_definition"], reference_state,
                    )
                    rrow = refs[rkey]
                    den = number(
                        rrow["reference_sasa_A2"],
                        f"reference denominator {rkey}",
                    )
                    ref_method = str(rrow.get("method", "")).strip()
                    fsasa = raw / den
                    if not math.isfinite(fsasa):
                        raise RuntimeError(f"Non-finite fSASA for {context}")
                    status = "computed"
                    counts["computed"] += 1
                    counts[reference_state] += 1

            if explicit_clean is None:
                clean = clean_local
            else:
                clean = explicit_clean
                if clean and not clean_local:
                    raise RuntimeError(
                        f"Global clean flag conflicts with local state for {context}: "
                        f"reference_state={reference_state!r}"
                    )

            if clean_local and not clean:
                counts["global_unclean_local_clean"] += 1

            if math.isfinite(fsasa) and fsasa < 0.0:
                counts["negative_fsasa"] += 1

        ph = canonical_ph(obs)
        if not ph:
            ph = canonical_ph(state)

        output.append({
            "trajectory": str(obs.get("trajectory", "")).strip(),
            "pH": ph,
            "analysis_row": str(obs.get("analysis_row", "")).strip(),
            "trajectory_frame": str(obs.get("trajectory_frame", "")).strip(),
            "step": str(obs.get("step", "")).strip(),
            "observable": observable,
            "metric": "fsasa",
            "selection": str(obs.get("selection", "")).strip(),
            "solutemask": str(obs.get("solutemask", "")).strip(),
            "value": fsasa,
            "unit": "1",
        })

        audit.append({
            "trajectory": str(obs.get("trajectory", "")).strip(),
            "pH": ph,
            "analysis_row": str(obs.get("analysis_row", "")).strip(),
            "trajectory_frame": str(obs.get("trajectory_frame", "")).strip(),
            "step": str(obs.get("step", "")).strip(),
            "observable": observable,
            "raw_metric": str(obs.get("metric", "")).strip(),
            "selection": str(obs.get("selection", "")).strip(),
            "solutemask": str(obs.get("solutemask", "")).strip(),
            "raw_sasa_A2": raw,
            "chemistry": spec["chemistry"],
            "site_definition": spec["site_definition"],
            "lambda_mode": args.lambda_mode,
            "lambda_column": lam_col,
            "local_lambda": lam,
            "tautomer_column": args.tautomer_column if args.lambda_mode == "his" else "",
            "local_tautomer_coordinate": x,
            "local_protonation": local_prot,
            "local_tautomer": local_taut,
            "reference_state": reference_state,
            "state_clean": str(clean).lower(),
            "state_matched": str(matched).lower(),
            "state_code": str(state.get("state_code", "")).strip(),
            "state_label": str(state.get("state_label", "")).strip(),
            "source_state": str(state.get(args.state_class_column, "")).strip(),
            "source_tautomer": str(state.get(args.tautomer_label_column, "")).strip(),
            "reference_forcefield": args.forcefield,
            "reference_sasa_A2": den,
            "reference_method": ref_method,
            "fsasa": fsasa,
            "fsasa_status": status,
        })

    unused = sorted(set(sites) - observables_seen)
    if unused:
        raise RuntimeError(
            "Configured --site mappings absent from SASA table: "
            + ", ".join(unused)
        )

    audit_path = (
        args.audit_output
        if args.audit_output is not None
        else args.output.with_name(args.output.stem + ".audit.tsv")
    )

    write_tsv(args.output, output, OUT_FIELDS)
    write_tsv(audit_path, audit, AUDIT_FIELDS)

    settings = args.output.with_suffix(args.output.suffix + ".settings.txt")
    with settings.open("w") as fh:
        fh.write(f"sasa={args.sasa.resolve()}\n")
        fh.write(f"states={args.states.resolve()}\n")
        fh.write(f"references={args.references.resolve()}\n")
        fh.write(f"forcefield={args.forcefield}\n")
        fh.write(f"lambda_mode={args.lambda_mode}\n")
        fh.write(f"join_on={','.join(join_cols)}\n")
        fh.write(f"trajectory_key_mode={args.trajectory_key_mode}\n")
        fh.write(f"state_trajectory_column={args.state_trajectory_column}\n")
        fh.write(f"state_analysis_row_column={args.state_analysis_row_column}\n")
        fh.write(f"state_trajectory_frame_column={args.state_trajectory_frame_column}\n")
        fh.write(f"state_step_column={args.state_step_column}\n")
        fh.write(f"low={args.low:.10g}\n")
        fh.write(f"high={args.high:.10g}\n")
        fh.write(f"clean_column={args.clean_column}\n")
        fh.write(f"matched_column={args.matched_column}\n")
        fh.write(f"sasa_metric={args.sasa_metric}\n")

        if args.lambda_mode == "his":
            fh.write(f"tautomer_column={args.tautomer_column}\n")
            fh.write(f"x_low_tautomer={args.x_low_tautomer}\n")
            fh.write(f"x_high_tautomer={args.x_high_tautomer}\n")
            fh.write(f"state_class_column={args.state_class_column}\n")
            fh.write(f"tautomer_label_column={args.tautomer_label_column}\n")
            fh.write(
                f"validate_his_state_labels={str(args.validate_his_state_labels).lower()}\n"
            )

        for obs, spec in sites.items():
            fh.write(f"site.{obs}.chemistry={spec['chemistry']}\n")
            fh.write(
                f"site.{obs}.site_definition={spec['site_definition']}\n"
            )
            fh.write(f"site.{obs}.lambda_column={spec['lambda_column']}\n")
            for reference_state in required_reference_states(
                args, spec["chemistry"]
            ):
                key = ref_key(
                    args.forcefield, spec["chemistry"],
                    spec["site_definition"], reference_state,
                )
                fh.write(
                    f"reference.{obs}.{reference_state}.A2="
                    f"{refs[key]['reference_sasa_A2']}\n"
                )

    print(f"Input SASA rows                    : {counts['input']}")
    print(f"Exact state joins                  : {counts['joined']}")
    print(f"Outside state window               : {counts['outside']}")
    print(f"Local lambda-mixed rows            : {counts['local_mixed']}")
    if args.lambda_mode == "his":
        print(f"Neutral tautomer-mixed rows        : {counts['tautomer_mixed']}")
        print(f"Computed HIP frames                : {counts['HIP']}")
        print(f"Computed HID frames                : {counts['HID']}")
        print(f"Computed HIE frames                : {counts['HIE']}")
    else:
        print(f"Computed H frames                  : {counts['H']}")
        print(f"Computed deprot frames             : {counts['deprot']}")
    print(f"fSASA values computed              : {counts['computed']}")
    print(
        "Global-unclean but local-clean rows : "
        f"{counts['global_unclean_local_clean']}"
    )
    print(f"Negative raw SASA rows             : {counts['negative_raw']}")
    print(f"Negative fSASA rows                : {counts['negative_fsasa']}")
    print(f"Wrote {args.output}")
    print(f"Wrote {audit_path}")
    print(f"Wrote {settings}")


if __name__ == "__main__":
    main()

