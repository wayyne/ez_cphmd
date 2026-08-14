#!/usr/bin/env python3
"""
calc_refsasa.py

Generate protonation-state-specific reference SASA denominators for later
fSASA normalization from model-compound/model-peptide trajectories.

The reference identity is chemistry-level, not biological-site-level:

    forcefield + chemistry + site_definition + protonation

Example keys:
    ff19SB + ASP + carboxylate_oxygens + H
    ff19SB + ASP + carboxylate_oxygens + deprot
    ff19SB + GLU + carboxylate_oxygens + H
    ff19SB + GLU + carboxylate_oxygens + deprot

CPPTRAJ `surf` is used for the partial LCPO contribution of --selection to
--solutemask. Negative instantaneous subset contributions are preserved.
The ensemble-mean reference denominator must be finite and > 0.

Lambda state assignment is by exact absolute-MD-step matching only:
    lambda <= --low   -> H
    lambda >= --high  -> deprot
    otherwise         -> mixed

Frames outside the lambda-file range are retained in the audit table as
outside_lambda_window. A missing exact lambda record inside the range is fatal.

For block-SEM diagnostics, contiguous blocks are created on each trajectory's
complete analyzed coordinate timeline before H/deprot frames are selected.
"""

from __future__ import annotations

import argparse
import csv
from concurrent.futures import ThreadPoolExecutor, as_completed
import glob
import math
import re
import subprocess
from pathlib import Path
from typing import Iterable

import numpy as np


def natural_key(value: str | Path):
    return [int(x) if x.isdigit() else x.lower()
            for x in re.split(r"([0-9]+)", str(value))]


def safe_component(value: str) -> str:
    out = re.sub(r"[^A-Za-z0-9_.-]+", "_", value.strip())
    if not out:
        raise ValueError(f"Invalid empty label derived from {value!r}")
    return out


def expand_specs(specs: Iterable[str]) -> list[Path]:
    found: dict[str, Path] = {}
    for spec in specs:
        matches = glob.glob(spec)
        if matches:
            for match in matches:
                p = Path(match).resolve()
                if p.is_file():
                    found[str(p)] = p
        else:
            p = Path(spec).resolve()
            if not p.is_file():
                raise FileNotFoundError(f"No file matched: {spec}")
            found[str(p)] = p
    paths = sorted(found.values(), key=natural_key)
    if not paths:
        raise RuntimeError("No files were found.")
    return paths


def extract_ph(name: str) -> float:
    m = re.search(r"[pP][hH][_=-]?([0-9]+(?:[._][0-9]+)?)", name)
    return float(m.group(1).replace("_", ".")) if m else math.nan


def read_numeric(path: Path) -> np.ndarray:
    rows = []
    with path.open() as fh:
        for line in fh:
            text = line.strip()
            if not text or text.startswith("#"):
                continue
            fields = text.replace(",", " ").split()
            try:
                vals = [float(x) for x in fields]
            except ValueError:
                continue
            rows.append(vals)
    if not rows:
        return np.empty((0, 0), dtype=float)
    ncol = min(len(row) for row in rows)
    return np.asarray([row[:ncol] for row in rows], dtype=float)


def write_tsv(path: Path, rows: list[dict], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t",
                           extrasaction="ignore")
        w.writeheader()
        for row in rows:
            out = dict(row)
            for key, value in out.items():
                if isinstance(value, (float, np.floating)):
                    out[key] = ("" if not math.isfinite(float(value))
                                else f"{float(value):.10g}")
            w.writerow(out)


def qmask(mask: str) -> str:
    return '"' + mask.replace('"', r'\"') + '"'


def read_lambda_file(path: Path, step_col: int,
                     value_col: int) -> dict[int, float]:
    data: dict[int, float] = {}
    with path.open() as fh:
        for lineno, line in enumerate(fh, 1):
            text = line.strip()
            if not text or text.startswith("#"):
                continue
            fields = text.replace(",", " ").split()
            if len(fields) <= max(step_col, value_col):
                continue
            try:
                step_f = float(fields[step_col])
                lam = float(fields[value_col])
            except ValueError:
                continue
            step = int(round(step_f))
            if abs(step_f - step) > 1.0e-6:
                raise RuntimeError(
                    f"{path}:{lineno}: non-integer step {step_f}")
            if not math.isfinite(lam):
                raise RuntimeError(
                    f"{path}:{lineno}: non-finite lambda {lam}")
            if step in data and abs(data[step] - lam) > 1.0e-12:
                raise RuntimeError(
                    f"{path}:{lineno}: conflicting lambda at step {step}: "
                    f"{data[step]} vs {lam}")
            data[step] = lam
    if not data:
        raise RuntimeError(f"No lambda records parsed from {path}")
    return data


def classify_lambda(lam: float, low: float, high: float) -> str:
    if lam <= low:
        return "H"
    if lam >= high:
        return "deprot"
    return "mixed"


def run_cpptraj_sasa(*, topology: Path, trajectory: Path, selection: str,
                     solutemask: str, image_anchor: str, offset: float,
                     nbrcut: float, start: int, stop: str, stride: int,
                     cpptraj: str, case_dir: Path, force: bool) -> np.ndarray:
    run_dir = case_dir / trajectory.stem
    run_dir.mkdir(parents=True, exist_ok=True)
    sasa_file = run_dir / "reference_sasa.dat"
    input_file = run_dir / "refsasa.in"
    log_file = run_dir / "cpptraj.log"

    input_file.write_text("\n".join([
        f'parm "{topology.resolve()}"',
        f'trajin "{trajectory.resolve()}" {start} {stop} {stride}',
        f"autoimage anchor {qmask(image_anchor)}",
        (f"surf REF_SASA {qmask(selection)} out \"{sasa_file.resolve()}\" "
         f"solutemask {qmask(solutemask)} offset {offset:.10g} "
         f"nbrcut {nbrcut:.10g}"),
        "run", "quit", ""
    ]))

    if force or not sasa_file.exists():
        with log_file.open("w") as log_fh:
            proc = subprocess.run([cpptraj, "-i", str(input_file)],
                                  stdout=log_fh,
                                  stderr=subprocess.STDOUT)
        if proc.returncode != 0:
            raise RuntimeError(
                f"CPPTRAJ failed for {trajectory.name}; see {log_file}")

    if not sasa_file.exists():
        raise RuntimeError(f"Expected CPPTRAJ output not found: {sasa_file}")

    arr = read_numeric(sasa_file)
    if arr.ndim != 2 or arr.shape[1] < 2:
        raise RuntimeError(f"Could not read SASA series from {sasa_file}")
    values = np.asarray(arr[:, 1], dtype=float)
    if np.any(~np.isfinite(values)):
        raise RuntimeError(f"Non-finite SASA values in {sasa_file}")

    # Partial LCPO contributions can legitimately be negative. Preserve them.
    return values


def summarize(timeline_by_run: dict[str, list[dict]], protonation: str,
              nblocks: int) -> dict:
    selected: list[float] = []
    block_means: list[float] = []
    nnegative = 0

    for run in sorted(timeline_by_run, key=natural_key):
        rr = sorted(timeline_by_run[run],
                    key=lambda row: row["trajectory_frame"])

        for row in rr:
            if row["state"] == protonation:
                value = row["sasa_A2"]
                selected.append(value)
                if value < 0.0:
                    nnegative += 1

        # Full-timeline blocks first; state selection second.
        blocks = np.array_split(np.arange(len(rr)),
                                min(nblocks, len(rr)))
        for block in blocks:
            vals = [rr[int(i)]["sasa_A2"] for i in block
                    if rr[int(i)]["state"] == protonation]
            if vals:
                block_means.append(float(np.mean(vals)))

    x = np.asarray(selected, dtype=float)
    bm = np.asarray(block_means, dtype=float)
    if x.size == 0:
        return {"nframes": 0, "mean": np.nan, "sd": np.nan,
                "block_sem": np.nan, "nblocks": 0, "median": np.nan,
                "q10": np.nan, "q90": np.nan, "min": np.nan,
                "max": np.nan, "nnegative": 0,
                "negative_fraction": np.nan}

    sem = (float(bm.std(ddof=1) / math.sqrt(len(bm)))
           if len(bm) >= 2 else np.nan)
    return {
        "nframes": int(x.size),
        "mean": float(np.mean(x)),
        "sd": float(np.std(x, ddof=1)) if x.size >= 2 else np.nan,
        "block_sem": sem,
        "nblocks": int(len(bm)),
        "median": float(np.median(x)),
        "q10": float(np.quantile(x, 0.10)),
        "q90": float(np.quantile(x, 0.90)),
        "min": float(np.min(x)),
        "max": float(np.max(x)),
        "nnegative": int(nnegative),
        "negative_fraction": float(nnegative / x.size),
    }


REFERENCE_FIELDS = [
    "forcefield", "chemistry", "site_definition", "protonation",
    "reference_sasa_A2", "sd_A2", "block_sem_A2", "nframes", "nblocks",
    "median_A2", "q10_A2", "q90_A2", "min_A2", "max_A2",
    "nnegative", "negative_fraction", "selection", "solutemask", "method",
    "offset_A", "nbrcut_A", "low_cutoff", "high_cutoff", "topology_file",
    "n_trajectories", "blocks_per_trajectory",
]


def update_reference_table(path: Path, new_rows: list[dict]) -> None:
    rows: list[dict] = []
    if path.exists():
        with path.open(newline="") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            required = {"forcefield", "chemistry", "site_definition",
                        "protonation"}
            if reader.fieldnames is None or not required.issubset(reader.fieldnames):
                raise RuntimeError(
                    f"Existing {path} uses an incompatible schema. "
                    "Write the chemistry-keyed references to a new file.")
            rows = list(reader)

    replacement_keys = {
        (r["forcefield"], r["chemistry"], r["site_definition"],
         r["protonation"]) for r in new_rows
    }
    rows = [r for r in rows if
            (r.get("forcefield", ""), r.get("chemistry", ""),
             r.get("site_definition", ""), r.get("protonation", ""))
            not in replacement_keys]
    rows.extend(new_rows)
    state_order = {"H": 0, "deprot": 1}
    rows.sort(key=lambda r: (
        str(r.get("forcefield", "")), str(r.get("chemistry", "")),
        str(r.get("site_definition", "")),
        state_order.get(str(r.get("protonation", "")), 99)))
    write_tsv(path, rows, REFERENCE_FIELDS)


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(
        description="Generate model-compound reference SASA denominators for fSASA.")
    ap.add_argument("--forcefield", required=True,
                    help="Reference key, e.g. ff19SB, ff14SB, c22.")
    ap.add_argument("--chemistry", required=True,
                    help="Chemical identity, e.g. ASP or GLU (not D52/E35).")
    ap.add_argument("--site-definition", required=True,
                    help="Stable definition, e.g. carboxylate_oxygens.")
    ap.add_argument("--selection", required=True,
                    help="CPPTRAJ reference selection, e.g. ':3@OD1,OD2'.")
    ap.add_argument("-p", "--topology", type=Path, required=True)
    ap.add_argument("-t", "--trajectory", action="append", required=True,
                    help="Trajectory path/glob; may be repeated.")
    ap.add_argument("-l", "--lambda-file", action="append", required=True,
                    help="Lambda path/glob; paired by basename stem.")
    ap.add_argument("--solutemask", "--solute-mask", dest="solutemask",
                    required=True,
                    help="Complete model compound/peptide for `surf`.")
    ap.add_argument("--image-anchor", default=None,
                    help="autoimage anchor; default is --solutemask.")
    ap.add_argument("--offset", type=float, default=1.4,
                    help="CPPTRAJ surf offset in A (default 1.4).")
    ap.add_argument("--nbrcut", type=float, default=2.5,
                    help="CPPTRAJ surf neighbor cutoff in A (default 2.5).")
    ap.add_argument("--first-frame-step", type=int, required=True)
    ap.add_argument("--frame-step-interval", type=int, required=True)
    ap.add_argument("--start", type=int, default=1)
    ap.add_argument("--stop", default="last")
    ap.add_argument("--stride", type=int, default=1)
    ap.add_argument("--lambda-step-col", type=int, default=0)
    ap.add_argument("--lambda-value-col", type=int, default=1)
    ap.add_argument("--low", type=float, default=0.2)
    ap.add_argument("--high", type=float, default=0.8)
    ap.add_argument("--blocks", type=int, default=10)
    ap.add_argument("-j", "--jobs", "--fork", dest="jobs", type=int,
                    default=1, metavar="N")
    ap.add_argument("--cpptraj", default="cpptraj")
    ap.add_argument("--workdir", type=Path,
                    default=Path("model_refsasa_work"))
    ap.add_argument("-o", "--output", type=Path,
                    default=Path("refsasa.tsv"))
    ap.add_argument("--force", action="store_true",
                    help="Re-run CPPTRAJ intermediates.")
    return ap.parse_args()


def validate_args(args: argparse.Namespace) -> None:
    if not args.topology.is_file():
        raise FileNotFoundError(args.topology)
    for name in ("forcefield", "chemistry", "site_definition", "selection",
                 "solutemask"):
        if not str(getattr(args, name)).strip():
            raise ValueError(f"--{name.replace('_', '-')} must not be empty")
    if args.start < 1 or args.stride < 1:
        raise ValueError("--start and --stride must be >= 1")
    if args.frame_step_interval <= 0:
        raise ValueError("--frame-step-interval must be > 0")
    if args.blocks < 1 or args.jobs < 1:
        raise ValueError("--blocks and --jobs must be >= 1")
    if args.lambda_step_col < 0 or args.lambda_value_col < 0:
        raise ValueError("lambda column indices must be >= 0")
    if not (0.0 <= args.low < args.high <= 1.0):
        raise ValueError("Require 0 <= --low < --high <= 1")
    if not math.isfinite(args.offset) or args.offset < 0.0:
        raise ValueError("--offset must be finite and >= 0")
    if not math.isfinite(args.nbrcut) or args.nbrcut < 0.0:
        raise ValueError("--nbrcut must be finite and >= 0")
    if args.stop != "last":
        try:
            stop = int(args.stop)
        except ValueError as exc:
            raise ValueError("--stop must be a positive integer or 'last'") from exc
        if stop < args.start:
            raise ValueError("--stop must be >= --start")


def main() -> None:
    args = parse_args()
    validate_args(args)
    topology = args.topology.resolve()
    trajectories = expand_specs(args.trajectory)
    lambda_files = expand_specs(args.lambda_file)
    image_anchor = args.image_anchor or args.solutemask

    lambda_by_stem: dict[str, Path] = {}
    for path in lambda_files:
        if path.stem in lambda_by_stem:
            raise RuntimeError(f"Duplicate lambda-file stem: {path.stem}")
        lambda_by_stem[path.stem] = path

    pairs: list[tuple[Path, Path]] = []
    for trajectory in trajectories:
        lam = lambda_by_stem.get(trajectory.stem)
        if lam is None:
            raise RuntimeError(
                f"No lambda file with matching stem for {trajectory.name}")
        pairs.append((trajectory, lam))

    used = {str(lam) for _, lam in pairs}
    unused = [p for p in lambda_files if str(p) not in used]
    if unused:
        raise RuntimeError("Lambda files without matching trajectories: " +
                           ", ".join(str(p) for p in unused))

    stems = [traj.stem for traj, _ in pairs]
    if len(stems) != len(set(stems)):
        raise RuntimeError("Trajectory filename stems must be unique.")

    case_dir = (args.workdir / safe_component(args.forcefield) /
                safe_component(args.chemistry) /
                safe_component(args.site_definition))
    case_dir.mkdir(parents=True, exist_ok=True)

    # Parse lambda first so malformed state input fails before expensive SASA.
    lambda_maps = {
        traj.stem: read_lambda_file(lam, args.lambda_step_col,
                                    args.lambda_value_col)
        for traj, lam in pairs
    }

    sasa_by_stem: dict[str, np.ndarray] = {}
    failures: list[tuple[Path, BaseException]] = []
    workers = min(args.jobs, len(pairs))

    def calculate(pair):
        traj, _ = pair
        values = run_cpptraj_sasa(
            topology=topology, trajectory=traj, selection=args.selection,
            solutemask=args.solutemask, image_anchor=image_anchor,
            offset=args.offset, nbrcut=args.nbrcut, start=args.start,
            stop=args.stop, stride=args.stride, cpptraj=args.cpptraj,
            case_dir=case_dir, force=args.force)
        return traj.stem, values

    if workers == 1:
        for pair in pairs:
            traj, _ = pair
            print(f"[{traj.stem}] running")
            try:
                stem, values = calculate(pair)
                sasa_by_stem[stem] = values
            except BaseException as exc:
                failures.append((traj, exc))
                break
            print(f"[{traj.stem}] complete")
    else:
        print(f"Processing {len(pairs)} trajectories with {workers} "
              "concurrent CPPTRAJ jobs.")
        with ThreadPoolExecutor(max_workers=workers) as ex:
            futures = {ex.submit(calculate, pair): pair for pair in pairs}
            for fut in as_completed(futures):
                traj, _ = futures[fut]
                try:
                    stem, values = fut.result()
                    sasa_by_stem[stem] = values
                except BaseException as exc:
                    failures.append((traj, exc))
                    print(f"[{traj.stem}] FAILED")
                else:
                    print(f"[{traj.stem}] complete")

    if failures:
        details = "\n".join(
            f"  {traj}: {type(exc).__name__}: {exc}"
            for traj, exc in failures)
        raise RuntimeError(
            f"{len(failures)} trajectory job(s) failed; reference table "
            f"was not updated:\n{details}")

    lambda_file_by_stem = {traj.stem: lam for traj, lam in pairs}
    audit_rows: list[dict] = []
    timeline_by_run: dict[str, list[dict]] = {}
    total_frames = total_matched = total_outside = total_mixed = total_clean = 0

    for traj, _ in pairs:
        stem = traj.stem
        values = sasa_by_stem[stem]
        lmap = lambda_maps[stem]
        lmin, lmax = min(lmap), max(lmap)
        timeline = []

        for i, value in enumerate(values):
            analysis_row = i + 1
            trajectory_frame = args.start + i * args.stride
            step = (args.first_frame_step +
                    (trajectory_frame - 1) * args.frame_step_interval)
            matched = 0
            clean = 0
            lam_out: float | str = ""

            if step < lmin or step > lmax:
                state = "outside_lambda_window"
                total_outside += 1
            else:
                if step not in lmap:
                    raise RuntimeError(
                        f"{stem}: frame {trajectory_frame} maps to step {step}, "
                        f"inside lambda range {lmin}--{lmax}, but has no exact "
                        "lambda record. No nearest-neighbor matching is allowed.")
                matched = 1
                total_matched += 1
                lam_out = lmap[step]
                state = classify_lambda(lam_out, args.low, args.high)
                if state == "mixed":
                    total_mixed += 1
                else:
                    clean = 1
                    total_clean += 1

            audit_rows.append({
                "forcefield": args.forcefield,
                "chemistry": args.chemistry,
                "site_definition": args.site_definition,
                "trajectory": stem,
                "pH": extract_ph(stem),
                "lambda_file": lambda_file_by_stem[stem].name,
                "analysis_row": analysis_row,
                "trajectory_frame": trajectory_frame,
                "step": step,
                "lambda": lam_out,
                "protonation": state,
                "clean": clean,
                "matched": matched,
                "sasa_A2": float(value),
            })
            timeline.append({"trajectory_frame": trajectory_frame,
                             "state": state,
                             "sasa_A2": float(value)})
            total_frames += 1
        timeline_by_run[stem] = timeline

    audit_path = case_dir / "reference_frame_assignments.tsv"
    write_tsv(audit_path, audit_rows, [
        "forcefield", "chemistry", "site_definition", "trajectory", "pH",
        "lambda_file", "analysis_row", "trajectory_frame", "step", "lambda",
        "protonation", "clean", "matched", "sasa_A2"])

    print("\nFrame join:")
    print(f"  coordinate frames     : {total_frames}")
    print(f"  exact matched frames  : {total_matched}")
    print(f"  outside lambda range  : {total_outside}")
    print(f"  mixed/intermediate    : {total_mixed}")
    print(f"  clean H/deprot frames : {total_clean}")

    new_rows = []
    print("\nReference SASA:")
    for protonation in ("H", "deprot"):
        stat = summarize(timeline_by_run, protonation, args.blocks)
        if stat["nframes"] == 0:
            raise RuntimeError(
                f"No clean {protonation} frames; denominator cannot be generated.")
        if not math.isfinite(stat["mean"]) or stat["mean"] <= 0.0:
            raise RuntimeError(
                f"{protonation} reference mean must be finite and > 0; "
                f"observed {stat['mean']}")

        print(f"  {protonation:6s}: {stat['mean']:.6f} A^2  "
              f"SD={stat['sd']:.6f}  block_SEM={stat['block_sem']:.6f}  "
              f"n={stat['nframes']}  blocks={stat['nblocks']}  "
              f"negative={stat['nnegative']} "
              f"({100.0 * stat['negative_fraction']:.3f}%)")

        new_rows.append({
            "forcefield": args.forcefield,
            "chemistry": args.chemistry,
            "site_definition": args.site_definition,
            "protonation": protonation,
            "reference_sasa_A2": stat["mean"],
            "sd_A2": stat["sd"],
            "block_sem_A2": stat["block_sem"],
            "nframes": stat["nframes"],
            "nblocks": stat["nblocks"],
            "median_A2": stat["median"],
            "q10_A2": stat["q10"],
            "q90_A2": stat["q90"],
            "min_A2": stat["min"],
            "max_A2": stat["max"],
            "nnegative": stat["nnegative"],
            "negative_fraction": stat["negative_fraction"],
            "selection": args.selection,
            "solutemask": args.solutemask,
            "method": "cpptraj_surf_lcpo",
            "offset_A": args.offset,
            "nbrcut_A": args.nbrcut,
            "low_cutoff": args.low,
            "high_cutoff": args.high,
            "topology_file": topology.name,
            "n_trajectories": len(pairs),
            "blocks_per_trajectory": args.blocks,
        })

    update_reference_table(args.output, new_rows)

    settings_path = case_dir / "refsasa.settings.txt"
    with settings_path.open("w") as fh:
        fh.write(f"forcefield={args.forcefield}\n")
        fh.write(f"chemistry={args.chemistry}\n")
        fh.write(f"site_definition={args.site_definition}\n")
        fh.write(f"selection={args.selection}\n")
        fh.write(f"topology={topology}\n")
        fh.write("trajectories=" + ",".join(str(t) for t, _ in pairs) + "\n")
        fh.write("lambda_files=" + ",".join(str(l) for _, l in pairs) + "\n")
        fh.write(f"solutemask={args.solutemask}\n")
        fh.write(f"image_anchor={image_anchor}\n")
        fh.write(f"offset={args.offset:.10g}\n")
        fh.write(f"nbrcut={args.nbrcut:.10g}\n")
        fh.write(f"first_frame_step={args.first_frame_step}\n")
        fh.write(f"frame_step_interval={args.frame_step_interval}\n")
        fh.write(f"start={args.start}\nstop={args.stop}\nstride={args.stride}\n")
        fh.write(f"lambda_step_col={args.lambda_step_col}\n")
        fh.write(f"lambda_value_col={args.lambda_value_col}\n")
        fh.write(f"low={args.low:.10g}\nhigh={args.high:.10g}\n")
        fh.write(f"blocks={args.blocks}\njobs={args.jobs}\ncpptraj={args.cpptraj}\n")

    print(f"\nUpdated: {args.output}")
    print(f"Audit:   {audit_path}")
    print(f"Settings: {settings_path}")


if __name__ == "__main__":
    main()
