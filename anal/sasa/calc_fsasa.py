#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path


OUT_FIELDS = [
    "trajectory", "pH", "analysis_row", "trajectory_frame", "step",
    "observable", "metric", "selection", "solutemask", "value", "unit",
]

AUDIT_FIELDS = [
    "trajectory", "pH", "analysis_row", "trajectory_frame", "step",
    "observable", "raw_metric", "selection", "solutemask", "raw_sasa_A2",
    "chemistry", "site_definition", "lambda_column", "local_lambda",
    "local_protonation", "state_clean", "state_matched", "state_code",
    "state_label", "reference_forcefield", "reference_sasa_A2",
    "reference_method", "fsasa", "fsasa_status",
]


def read_tsv(path: Path):
    if not path.is_file():
        raise FileNotFoundError(path)
    with path.open(newline="") as fh:
        r = csv.DictReader(fh, delimiter="\t")
        if r.fieldnames is None:
            raise RuntimeError(f"TSV has no header: {path}")
        return list(r.fieldnames), list(r)


def write_tsv(path: Path, rows, fields):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(
            fh, fieldnames=fields, delimiter="\t", extrasaction="ignore"
        )
        w.writeheader()
        for row in rows:
            out = dict(row)
            for k, v in list(out.items()):
                if isinstance(v, float):
                    out[k] = "" if not math.isfinite(v) else f"{v:.10g}"
            w.writerow(out)


def nonempty(v):
    return v is not None and bool(str(v).strip())


def truthy(v):
    return str(v).strip().lower() in {"1", "true", "yes", "y"}


def number(v, context):
    if not nonempty(v):
        raise RuntimeError(f"Missing numeric {context}")
    try:
        x = float(str(v).strip())
    except ValueError as exc:
        raise RuntimeError(f"Non-numeric {context}: {v!r}") from exc
    if not math.isfinite(x):
        raise RuntimeError(f"Non-finite {context}: {v!r}")
    return x


def integer(v, context):
    x = number(v, context)
    i = int(round(x))
    if abs(x - i) > 1e-8:
        raise RuntimeError(f"Non-integral {context}: {v!r}")
    return i


def canonical_ph(row):
    vals = []
    for key in ("pH", "ph"):
        if nonempty(row.get(key)):
            vals.append(number(row[key], key))
    if not vals:
        return ""
    if any(not math.isclose(vals[0], x, rel_tol=0, abs_tol=1e-8) for x in vals[1:]):
        raise RuntimeError(f"Conflicting pH/ph values in row: {row}")
    return f"{vals[0]:.10g}"


def make_key(row, columns, context):
    vals = []
    for col in columns:
        if not nonempty(row.get(col)):
            raise RuntimeError(f"Missing join key {col!r} in {context}")
        text = str(row[col]).strip()
        if col in {"analysis_row", "trajectory_frame", "step"}:
            text = str(integer(text, f"{context} {col}"))
        vals.append(text)
    return tuple(vals)


def validate_same_frame(obs, state, context):
    for col in ("analysis_row", "trajectory_frame", "step"):
        if nonempty(obs.get(col)) and nonempty(state.get(col)):
            a = integer(obs[col], f"observation {col} ({context})")
            b = integer(state[col], f"state {col} ({context})")
            if a != b:
                raise RuntimeError(
                    f"Frame metadata mismatch for {context}: "
                    f"{col} observation={a}, state={b}"
                )
    op = canonical_ph(obs)
    sp = canonical_ph(state)
    if op and sp and not math.isclose(float(op), float(sp), rel_tol=0, abs_tol=1e-8):
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
            raise RuntimeError(f"Duplicate reference key at {path}:{lineno}: {key}")
        lookup[key] = row
    return lookup


def build_sites(site_args):
    out = {}
    for obs, chemistry, site_definition, lambda_column in site_args:
        if obs in out:
            raise ValueError(f"Duplicate --site mapping for {obs!r}")
        out[obs] = {
            "chemistry": chemistry.upper(),
            "site_definition": site_definition,
            "lambda_column": lambda_column,
        }
    return out


def classify_lambda(lam, low, high):
    if lam <= low:
        return "H"
    if lam >= high:
        return "deprot"
    return "mixed"


def parse_args():
    ap = argparse.ArgumentParser(
        description=(
            "Calculate frame-wise fSASA by exact-joining raw SASA to local "
            "lambda coordinates and chemistry/protonation-specific references."
        )
    )
    ap.add_argument("--sasa", type=Path, required=True)
    ap.add_argument("--states", type=Path, required=True)
    ap.add_argument("--references", type=Path, required=True)
    ap.add_argument("--forcefield", required=True)
    ap.add_argument(
        "--site", action="append", nargs=4, required=True,
        metavar=("OBSERVABLE", "CHEMISTRY", "SITE_DEFINITION", "LAMBDA_COLUMN"),
        help=(
            "Example: --site E35_COO GLU carboxylate_oxygens primary_lambda"
        ),
    )
    ap.add_argument(
        "--join-on", default="trajectory,trajectory_frame",
        help="Default: trajectory,trajectory_frame",
    )
    ap.add_argument("--low", type=float, default=0.2)
    ap.add_argument("--high", type=float, default=0.8)
    ap.add_argument("--clean-column", default="clean")
    ap.add_argument("--matched-column", default="matched")
    ap.add_argument("--sasa-metric", default="sasa_lcpo")
    ap.add_argument("--audit-output", type=Path)
    ap.add_argument("-o", "--output", type=Path, default=Path("fsasa_frames.tsv"))
    return ap.parse_args()


def main():
    args = parse_args()

    if not (0 <= args.low < args.high <= 1):
        raise ValueError("Require 0 <= --low < --high <= 1")

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

    for col in join_cols:
        if col not in sasa_fields:
            raise RuntimeError(f"Join column {col!r} absent from SASA table")
        if col not in state_fields:
            raise RuntimeError(f"Join column {col!r} absent from state table")

    for obs, spec in sites.items():
        if spec["lambda_column"] not in state_fields:
            raise RuntimeError(
                f"{obs!r} requests state lambda column "
                f"{spec['lambda_column']!r}, but it is absent"
            )
        for prot in ("H", "deprot"):
            key = ref_key(
                args.forcefield, spec["chemistry"],
                spec["site_definition"], prot
            )
            if key not in refs:
                raise RuntimeError(f"No reference row for {obs!r}: {key}")

    state_lookup = {}
    for lineno, row in enumerate(state_rows, 2):
        key = make_key(
            row, join_cols, f"state row {args.states}:{lineno}"
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
        "computed": 0,
        "global_unclean_local_clean": 0,
        "negative_raw": 0,
        "negative_fsasa": 0,
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

        if args.sasa_metric and str(obs.get("metric", "")).strip() != args.sasa_metric:
            raise RuntimeError(
                f"Unexpected metric at {args.sasa}:{lineno}: "
                f"{obs.get('metric')!r}"
            )

        if str(obs.get("unit", "")).strip() not in {
            "A2", "A^2", "Angstrom^2", "angstrom^2"
        }:
            raise RuntimeError(
                f"Unexpected SASA unit at {args.sasa}:{lineno}: "
                f"{obs.get('unit')!r}"
            )

        raw = number(obs.get("value"), f"raw SASA at {args.sasa}:{lineno}")
        if raw < 0:
            counts["negative_raw"] += 1

        key = make_key(obs, join_cols, f"SASA row {args.sasa}:{lineno}")
        state = state_lookup.get(key)
        if state is None:
            raise RuntimeError(
                f"No exact state row for {args.sasa}:{lineno}; key={key}"
            )
        counts["joined"] += 1

        context = f"{args.sasa}:{lineno}, key={key}"
        validate_same_frame(obs, state, context)

        ident = key + (observable,)
        if ident in seen:
            raise RuntimeError(f"Duplicate SASA observation: {ident}")
        seen.add(ident)

        spec = sites[observable]
        lam_col = spec["lambda_column"]
        matched = (
            truthy(state.get(args.matched_column))
            if args.matched_column in state_fields else True
        )
        clean = (
            truthy(state.get(args.clean_column))
            if args.clean_column in state_fields else False
        )

        lam = math.nan
        local_prot = ""
        den = math.nan
        ref_method = ""
        fsasa = math.nan

        if not matched:
            status = "outside_state_window"
            counts["outside"] += 1
        else:
            lam = number(state.get(lam_col), f"{lam_col} for {context}")
            if not 0 <= lam <= 1:
                raise RuntimeError(
                    f"{lam_col} outside [0,1] for {context}: {lam}"
                )

            local_prot = classify_lambda(lam, args.low, args.high)

            if local_prot == "mixed":
                status = "local_mixed"
                counts["local_mixed"] += 1
                if clean:
                    raise RuntimeError(
                        f"Global clean flag conflicts with intermediate "
                        f"{lam_col}={lam:.10g} for {context}"
                    )
            else:
                rkey = ref_key(
                    args.forcefield, spec["chemistry"],
                    spec["site_definition"], local_prot
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

                if not clean:
                    counts["global_unclean_local_clean"] += 1

                if fsasa < 0:
                    counts["negative_fsasa"] += 1

        ph = canonical_ph(obs)

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
            "lambda_column": lam_col,
            "local_lambda": lam,
            "local_protonation": local_prot,
            "state_clean": str(clean).lower(),
            "state_matched": str(matched).lower(),
            "state_code": str(state.get("state_code", "")).strip(),
            "state_label": str(state.get("state_label", "")).strip(),
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
        fh.write(f"join_on={','.join(join_cols)}\n")
        fh.write(f"low={args.low:.10g}\n")
        fh.write(f"high={args.high:.10g}\n")
        fh.write(f"clean_column={args.clean_column}\n")
        fh.write(f"matched_column={args.matched_column}\n")
        fh.write(f"sasa_metric={args.sasa_metric}\n")

        for obs, spec in sites.items():
            fh.write(f"site.{obs}.chemistry={spec['chemistry']}\n")
            fh.write(
                f"site.{obs}.site_definition={spec['site_definition']}\n"
            )
            fh.write(f"site.{obs}.lambda_column={spec['lambda_column']}\n")
            for prot in ("H", "deprot"):
                key = ref_key(
                    args.forcefield, spec["chemistry"],
                    spec["site_definition"], prot
                )
                fh.write(
                    f"reference.{obs}.{prot}.A2="
                    f"{refs[key]['reference_sasa_A2']}\n"
                )

    print(f"Input SASA rows                    : {counts['input']}")
    print(f"Exact state joins                  : {counts['joined']}")
    print(f"Outside state window               : {counts['outside']}")
    print(f"Local mixed rows                   : {counts['local_mixed']}")
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
