#!/usr/bin/env python3
"""
calc_refsasa.py

Generate chemistry/state-specific reference SASA denominators for later fSASA
normalization from model-compound/model-peptide trajectories.

The reference identity is chemistry-level, not biological-site-level:

    forcefield + chemistry + site_definition + protonation

For ordinary one-coordinate titratable groups (ASP/GLU/etc.), the existing
reference states are preserved:

    lambda <= --low   -> H
    lambda >= --high  -> deprot
    otherwise         -> mixed/excluded

For histidine, use ``--lambda-mode his``. The Amber lambda-file ``ires`` and
``itauto`` headers are parsed exactly as in the H15 analysis: for the selected
His ires, ``itauto=1`` is the protonation coordinate lambda and ``itauto=2`` is
the neutral-tautomer coordinate x. The physical reference states are then
assigned as:

    lambda <= --low   -> HIP, independent of x
    lambda >= --high and x <= --low  -> --x-low-tautomer
    lambda >= --high and x >= --high -> --x-high-tautomer

Intermediate lambda values and neutral frames with intermediate x are retained
in the audit table but excluded from reference denominators. Thus His reference
rows are generated separately for HIP, HID, and HIE.

CPPTRAJ ``surf`` is used for the partial LCPO contribution of --selection to
--solutemask. Negative instantaneous subset contributions are preserved. The
ensemble-mean reference denominator for every required state must be finite and
> 0.

Coordinate/SASA frames are joined to lambda records by exact absolute MD step.
Frames outside the lambda-file range are retained in the audit table. A missing
exact lambda record inside the range is fatal. His lambda step counters are
unwrapped using the same monotonic-counter logic as the H15 analysis before the
exact join.

For block-SEM diagnostics, contiguous blocks are created on each trajectory's
complete analyzed coordinate timeline before state selection.
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
from typing import Iterable, Sequence

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


def read_single_lambda_file(path: Path, step_col: int,
                            value_col: int) -> dict[int, float]:
    """Read the legacy one-coordinate lambda format without changing behavior."""
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


def _header_values(line: str, key: str) -> list[int] | None:
    toks = line.strip().split()
    norm = [tok.lstrip("#").lower() for tok in toks]
    if key.lower() not in norm:
        return None
    i = norm.index(key.lower())
    vals: list[int] = []
    for tok in toks[i + 1:]:
        try:
            vals.append(int(tok))
        except ValueError:
            pass
    return vals


def read_lambda_header(path: Path) -> tuple[list[int], list[int]]:
    """Read Amber CpHMD ires/itauto metadata from a lambda file."""
    ires: list[int] | None = None
    itauto: list[int] | None = None
    with path.open(errors="replace") as fh:
        for line in fh:
            text = line.strip()
            if not text:
                continue
            if ires is None:
                vals = _header_values(text, "ires")
                if vals is not None:
                    ires = vals
            if itauto is None:
                vals = _header_values(text, "itauto")
                if vals is not None:
                    itauto = vals
            if ires is not None and itauto is not None:
                break
    if ires is None or itauto is None:
        raise RuntimeError(f"{path}: missing ires and/or itauto lambda header")
    if len(ires) != len(itauto):
        raise RuntimeError(f"{path}: ires and itauto header lengths differ")
    return ires, itauto


def lambda_pair_columns(ires: Sequence[int], itauto: Sequence[int],
                        resid: int) -> tuple[list[int], list[int]]:
    pvars = [
        i for i, (r, a) in enumerate(zip(ires, itauto))
        if r == resid and a == 1
    ]
    xvars = [
        i for i, (r, a) in enumerate(zip(ires, itauto))
        if r == resid and a == 2
    ]
    return pvars, xvars


def resolve_his_lambda_resid(ires: Sequence[int], itauto: Sequence[int],
                             requested: int | None, path: Path) -> int:
    """Resolve one His protonation/tautomer pair from ires/itauto metadata."""
    if requested is not None:
        pvars, xvars = lambda_pair_columns(ires, itauto, requested)
        if len(pvars) != 1 or len(xvars) != 1:
            raise RuntimeError(
                f"{path}: --lambda-resid={requested} does not identify exactly "
                f"one His pair; itauto=1 columns={pvars}, itauto=2 columns={xvars}"
            )
        return int(requested)

    candidates: list[int] = []
    for resid in sorted(set(int(x) for x in ires)):
        pvars, xvars = lambda_pair_columns(ires, itauto, resid)
        if len(pvars) == 1 and len(xvars) == 1:
            candidates.append(resid)

    if len(candidates) != 1:
        raise RuntimeError(
            f"{path}: could not auto-resolve a unique His lambda pair from "
            f"ires/itauto; candidates={candidates}. Supply --lambda-resid N."
        )
    return candidates[0]


def unwrap_steps(raw_steps: np.ndarray) -> tuple[np.ndarray, int]:
    """Make restarted/reset lambda step counters strictly increasing."""
    raw = np.asarray(raw_steps, dtype=np.int64)
    if raw.size == 0:
        return raw.copy(), 0

    out = np.empty_like(raw)
    out[0] = raw[0]
    offset = 0
    resets = 0
    positive_diffs: list[int] = []

    for i in range(1, len(raw)):
        d = int(raw[i] - raw[i - 1])
        if d > 0:
            positive_diffs.append(d)
        else:
            resets += 1
            typical = (
                int(round(float(np.median(positive_diffs))))
                if positive_diffs else 1
            )
            typical = max(1, typical)
            offset = int(out[i - 1] + typical - raw[i])
        out[i] = raw[i] + offset
        if out[i] <= out[i - 1]:
            raise RuntimeError(f"{i}: could not unwrap lambda step counters")

    return out, resets


def read_his_lambda_file(path: Path, lambda_resid: int | None
                         ) -> tuple[dict[int, tuple[float, float]], dict]:
    """
    Read a two-coordinate Amber histidine lambda file.

    The selected ires must have exactly one itauto=1 (protonation lambda) and
    one itauto=2 (tautomer x) variable. Numeric columns are step + lambda vars.
    """
    ires, itauto = read_lambda_header(path)
    resolved_resid = resolve_his_lambda_resid(
        ires, itauto, lambda_resid, path
    )
    pvars, xvars = lambda_pair_columns(ires, itauto, resolved_resid)
    pcol = pvars[0] + 1
    xcol = xvars[0] + 1
    nvar = len(ires)

    rows: list[list[float]] = []
    with path.open(errors="replace") as fh:
        for line in fh:
            text = line.strip()
            if not text or text.startswith("#"):
                continue
            try:
                vals = [float(x) for x in text.replace(",", " ").split()]
            except ValueError:
                continue
            if len(vals) >= nvar + 1:
                rows.append(vals[: nvar + 1])

    if not rows:
        raise RuntimeError(f"{path}: no usable numeric lambda rows")

    arr = np.asarray(rows, dtype=float)
    raw_step_float = arr[:, 0]
    raw_steps = np.rint(raw_step_float).astype(np.int64)
    if np.any(np.abs(raw_step_float - raw_steps) > 1.0e-6):
        bad = int(np.flatnonzero(np.abs(raw_step_float - raw_steps) > 1.0e-6)[0])
        raise RuntimeError(
            f"{path}: non-integer lambda step {raw_step_float[bad]} at numeric row {bad + 1}"
        )

    steps, resets = unwrap_steps(raw_steps)
    lam = np.asarray(arr[:, pcol], dtype=float)
    x = np.asarray(arr[:, xcol], dtype=float)
    if np.any(~np.isfinite(lam)) or np.any(~np.isfinite(x)):
        raise RuntimeError(f"{path}: non-finite His lambda/tautomer coordinate")

    data: dict[int, tuple[float, float]] = {}
    for step, lv, xv in zip(steps, lam, x):
        step_i = int(step)
        pair = (float(lv), float(xv))
        if step_i in data:
            old = data[step_i]
            if (abs(old[0] - pair[0]) > 1.0e-12 or
                    abs(old[1] - pair[1]) > 1.0e-12):
                raise RuntimeError(
                    f"{path}: conflicting His lambda values at step {step_i}: "
                    f"{old} vs {pair}"
                )
        data[step_i] = pair

    return data, {
        "lambda_resid": resolved_resid,
        "protonation_data_column_0based": pcol,
        "tautomer_data_column_0based": xcol,
        "counter_resets": resets,
        "ires": ires,
        "itauto": itauto,
    }


def classify_lambda(lam: float, low: float, high: float) -> str:
    if lam <= low:
        return "H"
    if lam >= high:
        return "deprot"
    return "mixed"


def classify_his(lam: float, x: float, low: float, high: float,
                 x_low_tautomer: str, x_high_tautomer: str
                 ) -> tuple[str, str]:
    """Return (reference_state, tautomer_label) for one His lambda pair."""
    if lam <= low:
        return "HIP", "HIP"
    if lam < high:
        return "mixed", "mixed"

    if x <= low:
        tautomer = x_low_tautomer
        return tautomer, tautomer
    if x >= high:
        tautomer = x_high_tautomer
        return tautomer, tautomer
    return "neutral_tautomer_mixed", "mixed"


def reference_states(args: argparse.Namespace) -> tuple[str, ...]:
    return ("HIP", "HID", "HIE") if args.lambda_mode == "his" else ("H", "deprot")

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
    "n_trajectories", "blocks_per_trajectory", "lambda_mode",
    "lambda_resid", "x_low_tautomer", "x_high_tautomer",
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
    state_order = {"H": 0, "deprot": 1, "HIP": 0, "HID": 1, "HIE": 2}
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
    ap.add_argument("--lambda-step-col", type=int, default=0,
                    help="Step column for --lambda-mode single.")
    ap.add_argument("--lambda-value-col", type=int, default=1,
                    help="Lambda-value column for --lambda-mode single.")
    ap.add_argument(
        "--lambda-mode", choices=("single", "his"), default="single",
        help=(
            "single: legacy one-coordinate H/deprot classification; "
            "his: parse Amber ires/itauto headers and generate HIP/HID/HIE references."
        ),
    )
    ap.add_argument(
        "--lambda-resid", type=int, default=None,
        help=(
            "His CpHMD ires for --lambda-mode his. If omitted, auto-resolve "
            "the unique ires having exactly one itauto=1 and one itauto=2 variable."
        ),
    )
    ap.add_argument(
        "--x-low-tautomer", choices=("HID", "HIE"), default=None,
        help="Physical tautomer at x <= --low for --lambda-mode his.",
    )
    ap.add_argument(
        "--x-high-tautomer", choices=("HID", "HIE"), default=None,
        help="Physical tautomer at x >= --high for --lambda-mode his.",
    )
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
    if args.lambda_resid is not None and args.lambda_resid < 1:
        raise ValueError("--lambda-resid must be >= 1")

    chemistry = args.chemistry.strip().upper()
    if args.lambda_mode == "his":
        if chemistry != "HIS":
            raise ValueError("--lambda-mode his requires --chemistry HIS")
        if args.x_low_tautomer is None or args.x_high_tautomer is None:
            raise ValueError(
                "--lambda-mode his requires --x-low-tautomer and --x-high-tautomer"
            )
        if args.x_low_tautomer == args.x_high_tautomer:
            raise ValueError(
                "His x endpoints must map to different tautomers (HID and HIE)"
            )
    elif chemistry == "HIS":
        raise ValueError(
            "--chemistry HIS requires --lambda-mode his so HID and HIE are not pooled"
        )
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
    lambda_maps: dict[str, dict] = {}
    lambda_meta: dict[str, dict] = {}
    for traj, lam in pairs:
        if args.lambda_mode == "his":
            lmap, meta = read_his_lambda_file(lam, args.lambda_resid)
            lambda_maps[traj.stem] = lmap
            lambda_meta[traj.stem] = meta
        else:
            lambda_maps[traj.stem] = read_single_lambda_file(
                lam, args.lambda_step_col, args.lambda_value_col
            )
            lambda_meta[traj.stem] = {
                "lambda_resid": "",
                "protonation_data_column_0based": args.lambda_value_col,
                "tautomer_data_column_0based": "",
                "counter_resets": 0,
            }

    if args.lambda_mode == "his":
        resolved = {int(meta["lambda_resid"]) for meta in lambda_meta.values()}
        if len(resolved) != 1:
            raise RuntimeError(
                f"His lambda ires differs among paired files: {sorted(resolved)}"
            )
        resolved_his_lambda_resid = next(iter(resolved))
        print("Resolved His lambda mapping:")
        print(f"  ires                    : {resolved_his_lambda_resid}")
        print(f"  x <= {args.low:g}              : {args.x_low_tautomer}")
        print(f"  x >= {args.high:g}              : {args.x_high_tautomer}")
        resets = sum(int(meta["counter_resets"]) for meta in lambda_meta.values())
        print(f"  lambda counter resets   : {resets}")
    else:
        resolved_his_lambda_resid = None

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
    total_frames = 0
    total_matched = 0
    total_outside = 0
    total_lambda_mixed = 0
    total_tautomer_mixed = 0
    total_clean = 0

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
            x_out: float | str = ""
            tautomer = ""

            if step < lmin or step > lmax:
                state = "outside_lambda_window"
                total_outside += 1
            else:
                if step not in lmap:
                    raise RuntimeError(
                        f"{stem}: frame {trajectory_frame} maps to step {step}, "
                        f"inside lambda range {lmin}--{lmax}, but has no exact "
                        "lambda record. No nearest-neighbor matching is allowed."
                    )
                matched = 1
                total_matched += 1

                if args.lambda_mode == "his":
                    lam_out, x_out = lmap[step]
                    state, tautomer = classify_his(
                        float(lam_out), float(x_out), args.low, args.high,
                        args.x_low_tautomer, args.x_high_tautomer,
                    )
                    if state == "mixed":
                        total_lambda_mixed += 1
                    elif state == "neutral_tautomer_mixed":
                        total_tautomer_mixed += 1
                    else:
                        clean = 1
                        total_clean += 1
                else:
                    lam_out = lmap[step]
                    state = classify_lambda(float(lam_out), args.low, args.high)
                    if state == "mixed":
                        total_lambda_mixed += 1
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
                "x": x_out,
                "protonation": state,
                "tautomer": tautomer,
                "clean": clean,
                "matched": matched,
                "sasa_A2": float(value),
            })
            timeline.append({
                "trajectory_frame": trajectory_frame,
                "state": state,
                "sasa_A2": float(value),
            })
            total_frames += 1
        timeline_by_run[stem] = timeline

    audit_path = case_dir / "reference_frame_assignments.tsv"
    write_tsv(audit_path, audit_rows, [
        "forcefield", "chemistry", "site_definition", "trajectory", "pH",
        "lambda_file", "analysis_row", "trajectory_frame", "step", "lambda",
        "x", "protonation", "tautomer", "clean", "matched", "sasa_A2",
    ])

    print("\nFrame join:")
    print(f"  coordinate frames       : {total_frames}")
    print(f"  exact matched frames    : {total_matched}")
    print(f"  outside lambda range    : {total_outside}")
    print(f"  lambda mixed/intermed.  : {total_lambda_mixed}")
    if args.lambda_mode == "his":
        print(f"  neutral x intermediate  : {total_tautomer_mixed}")
        print(f"  clean HIP/HID/HIE frames: {total_clean}")
    else:
        print(f"  clean H/deprot frames   : {total_clean}")

    new_rows = []
    print("\nReference SASA:")
    for protonation in reference_states(args):
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
            "lambda_mode": args.lambda_mode,
            "lambda_resid": (
                resolved_his_lambda_resid
                if resolved_his_lambda_resid is not None else ""
            ),
            "x_low_tautomer": args.x_low_tautomer or "",
            "x_high_tautomer": args.x_high_tautomer or "",
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
        fh.write(f"lambda_mode={args.lambda_mode}\n")
        fh.write(f"lambda_step_col={args.lambda_step_col}\n")
        fh.write(f"lambda_value_col={args.lambda_value_col}\n")
        fh.write(
            "lambda_resid="
            + ("" if resolved_his_lambda_resid is None
               else str(resolved_his_lambda_resid))
            + "\n"
        )
        fh.write(f"x_low_tautomer={args.x_low_tautomer or ''}\n")
        fh.write(f"x_high_tautomer={args.x_high_tautomer or ''}\n")
        fh.write(f"low={args.low:.10g}\nhigh={args.high:.10g}\n")
        fh.write(f"blocks={args.blocks}\njobs={args.jobs}\ncpptraj={args.cpptraj}\n")

    print(f"\nUpdated: {args.output}")
    print(f"Audit:   {audit_path}")
    print(f"Settings: {settings_path}")


if __name__ == "__main__":
    main()
