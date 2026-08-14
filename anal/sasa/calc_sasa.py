#!/usr/bin/env python3
"""
calc_sasa.py

State-agnostic, frame-resolved solvent-accessible surface-area analysis for
Amber trajectories using CPPTRAJ's LCPO `surf` action.

Design goals
------------
1. Know nothing about protonation, lambda coordinates, or coupled states.
2. Accept arbitrary named atom masks; one or many sites may be analyzed.
3. Require an explicit solutemask so the structural SASA definition is
   unambiguous and recorded in provenance.
4. Emit the same canonical frame identity used by the other observable tools.
5. Optionally process independent trajectories concurrently with --jobs/--fork.
6. Calculate raw SASA only. fSASA normalization belongs in a later derived-data
   layer after local protonation/reference information is available.

Example
-------
python3 calc_sasa.py \
    -p hewl.parm7 \
    -t 'arex.ph*.nc' \
    --site E35_COO ':36@OE1,OE2' \
    --site D52_COO ':53@OD1,OD2' \
    --solutemask '^1' \
    --image-anchor '^1' \
    --first-frame-step 510000 \
    --frame-step-interval 10000 \
    --jobs 8 \
    -o carboxylate_sasa_frames.tsv

The --site and --solutemask values are CPPTRAJ/topology masks. Biological-to-
topology residue mapping is intentionally outside this low-level calculator.
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
    return [
        int(x) if x.isdigit() else x.lower()
        for x in re.split(r"([0-9]+)", str(value))
    ]


def expand_trajectories(specs: Iterable[str]) -> list[Path]:
    found: dict[str, Path] = {}

    for spec in specs:
        matches = glob.glob(spec)

        if matches:
            for match in matches:
                path = Path(match).resolve()
                if path.is_file():
                    found[str(path)] = path
            continue

        path = Path(spec).resolve()

        if path.is_file():
            found[str(path)] = path
        else:
            raise FileNotFoundError(
                f"No trajectory matched: {spec}"
            )

    paths = sorted(
        found.values(),
        key=natural_key,
    )

    if not paths:
        raise RuntimeError(
            "No trajectory files were found."
        )

    return paths


def extract_ph(name: str) -> float:
    """
    Accept common AREX naming conventions:
      arex.ph0_5.nc -> 0.5
      arex.ph5_0.nc -> 5.0
      arex_ph7.0.nc -> 7.0
      pH_7.0.nc     -> 7.0
    """
    match = re.search(
        r"[pP][hH][_=-]?([0-9]+(?:[._][0-9]+)?)",
        name,
    )

    if not match:
        return math.nan

    return float(
        match.group(1).replace("_", ".")
    )


def read_numeric(path: Path) -> np.ndarray:
    rows: list[list[float]] = []

    with path.open() as fh:
        for line in fh:
            text = line.strip()

            if not text or text.startswith("#"):
                continue

            fields = text.replace(",", " ").split()

            try:
                values = [float(x) for x in fields]
            except ValueError:
                continue

            rows.append(values)

    if not rows:
        return np.empty(
            (0, 0),
            dtype=float,
        )

    ncolumns = min(
        len(row)
        for row in rows
    )

    return np.asarray(
        [
            row[:ncolumns]
            for row in rows
        ],
        dtype=float,
    )


def qmask(mask: str) -> str:
    """
    Quote a CPPTRAJ mask in generated input.
    """
    return '"' + mask.replace('"', r'\"') + '"'


def safe_name(name: str) -> str:
    clean = re.sub(
        r"[^A-Za-z0-9_.-]+",
        "_",
        name.strip(),
    )

    if not clean:
        raise ValueError(
            f"Invalid empty site name derived from {name!r}"
        )

    return clean


def write_tsv(
    path: Path,
    rows: list[dict],
    fields: list[str],
) -> None:
    path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=fields,
            delimiter="\t",
            extrasaction="ignore",
        )

        writer.writeheader()

        for row in rows:
            out = dict(row)

            for key, value in out.items():
                if isinstance(
                    value,
                    (float, np.floating),
                ):
                    out[key] = (
                        ""
                        if not math.isfinite(float(value))
                        else f"{float(value):.10g}"
                    )

            writer.writerow(out)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Calculate state-agnostic frame-wise LCPO SASA with CPPTRAJ "
            "and write a canonical long-format TSV."
        )
    )

    parser.add_argument(
        "-p",
        "--topology",
        type=Path,
        required=True,
        help="Amber topology file.",
    )

    parser.add_argument(
        "-t",
        "--trajectory",
        action="append",
        required=True,
        help=(
            "Trajectory path or glob; may be repeated. "
            "Quote shell globs."
        ),
    )

    parser.add_argument(
        "--site",
        action="append",
        nargs=2,
        required=True,
        metavar=("NAME", "MASK"),
        help=(
            "Named SASA observable: NAME MASK. May be repeated. "
            "MASK uses CPPTRAJ/topology numbering."
        ),
    )

    parser.add_argument(
        "--solutemask",
        required=True,
        help=(
            "CPPTRAJ mask defining the solute context used for partial SASA, "
            "for example '^1' for the protein molecule."
        ),
    )

    parser.add_argument(
        "--image-anchor",
        default=None,
        help=(
            "Optional CPPTRAJ autoimage anchor mask, "
            "for example '^1'."
        ),
    )

    parser.add_argument(
        "--offset",
        type=float,
        default=1.4,
        help=(
            "CPPTRAJ surf van-der-Waals radius offset in Angstrom "
            "(default: 1.4)."
        ),
    )

    parser.add_argument(
        "--nbrcut",
        type=float,
        default=2.5,
        help=(
            "CPPTRAJ surf neighbor-radius cutoff in Angstrom "
            "(default: 2.5)."
        ),
    )

    parser.add_argument(
        "--start",
        type=int,
        default=1,
        help="First stored trajectory frame to analyze (default: 1).",
    )

    parser.add_argument(
        "--stop",
        default="last",
        help="Last stored trajectory frame to analyze or 'last'.",
    )

    parser.add_argument(
        "--stride",
        type=int,
        default=1,
        help="Stored-frame stride (default: 1).",
    )

    parser.add_argument(
        "--first-frame-step",
        type=int,
        default=None,
        help=(
            "Absolute MD step corresponding to stored trajectory frame 1."
        ),
    )

    parser.add_argument(
        "--frame-step-interval",
        type=int,
        default=None,
        help=(
            "MD-step spacing between consecutive stored trajectory frames."
        ),
    )

    parser.add_argument(
        "--cpptraj",
        default="cpptraj",
        help="CPPTRAJ executable (default: cpptraj).",
    )

    parser.add_argument(
        "-j",
        "--jobs",
        "--fork",
        dest="jobs",
        type=int,
        default=1,
        metavar="N",
        help=(
            "Number of trajectories to process concurrently using independent "
            "CPPTRAJ processes (default: 1). --fork is an alias for --jobs."
        ),
    )

    parser.add_argument(
        "--workdir",
        type=Path,
        default=Path("sasa_work"),
        help=(
            "Directory for CPPTRAJ inputs, logs, and intermediate series "
            "(default: sasa_work)."
        ),
    )

    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=Path("sasa_frames.tsv"),
        help="Combined canonical frame-wise TSV.",
    )

    parser.add_argument(
        "--force",
        action="store_true",
        help="Re-run CPPTRAJ even when expected intermediate files exist.",
    )

    return parser.parse_args()


def validate_args(
    args: argparse.Namespace,
) -> None:
    if not args.topology.is_file():
        raise FileNotFoundError(
            args.topology
        )

    if args.start < 1:
        raise ValueError(
            "--start must be >= 1"
        )

    if args.stride < 1:
        raise ValueError(
            "--stride must be >= 1"
        )

    if args.jobs < 1:
        raise ValueError(
            "--jobs/--fork must be >= 1"
        )

    if not math.isfinite(args.offset) or args.offset < 0.0:
        raise ValueError(
            "--offset must be a finite value >= 0"
        )

    if not math.isfinite(args.nbrcut) or args.nbrcut < 0.0:
        raise ValueError(
            "--nbrcut must be a finite value >= 0"
        )

    if not args.solutemask.strip():
        raise ValueError(
            "--solutemask must not be empty"
        )

    if args.stop != "last":
        try:
            stop = int(args.stop)
        except ValueError as exc:
            raise ValueError(
                "--stop must be a positive integer or 'last'"
            ) from exc

        if stop < 1:
            raise ValueError(
                "--stop must be a positive integer or 'last'"
            )

        if stop < args.start:
            raise ValueError(
                "--stop must be >= --start"
            )

    have_first_step = (
        args.first_frame_step is not None
    )
    have_interval = (
        args.frame_step_interval is not None
    )

    if have_first_step != have_interval:
        raise ValueError(
            "--first-frame-step and --frame-step-interval "
            "must be supplied together"
        )

    if (
        have_interval
        and args.frame_step_interval <= 0
    ):
        raise ValueError(
            "--frame-step-interval must be > 0"
        )

    names = [
        safe_name(site[0])
        for site in args.site
    ]

    if len(names) != len(set(names)):
        raise ValueError(
            "Site NAME values must be unique."
        )


def run_one_trajectory(
    args: argparse.Namespace,
    trajectory: Path,
    sites: list[tuple[str, str]],
) -> dict[str, np.ndarray]:
    label = trajectory.stem

    run_dir = (
        args.workdir
        / label
    )

    run_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    outputs: dict[str, Path] = {}

    lines = [
        f'parm "{args.topology.resolve()}"',
        (
            f'trajin "{trajectory.resolve()}" '
            f'{args.start} {args.stop} {args.stride}'
        ),
    ]

    if args.image_anchor:
        lines.append(
            f"autoimage anchor {qmask(args.image_anchor)}"
        )

    for name, mask in sites:
        output = (
            run_dir
            / f"{safe_name(name)}.dat"
        )

        outputs[name] = output

        # CPPTRAJ `surf` calculates the LCPO SASA contribution of the
        # selected mask in the context of the explicit solutemask.
        lines.append(
            "surf "
            f"{safe_name(name)} "
            f"{qmask(mask)} "
            f'out "{output.resolve()}" '
            f"solutemask {qmask(args.solutemask)} "
            f"offset {args.offset:.10g} "
            f"nbrcut {args.nbrcut:.10g}"
        )

    lines.extend(
        [
            "run",
            "quit",
            "",
        ]
    )

    input_path = (
        run_dir
        / "sasa.in"
    )

    log_path = (
        run_dir
        / "cpptraj.log"
    )

    input_path.write_text(
        "\n".join(lines)
    )

    need_run = (
        args.force
        or any(
            not path.exists()
            for path in outputs.values()
        )
    )

    if need_run:
        with log_path.open("w") as log_fh:
            process = subprocess.run(
                [
                    args.cpptraj,
                    "-i",
                    str(input_path),
                ],
                stdout=log_fh,
                stderr=subprocess.STDOUT,
            )

        if process.returncode != 0:
            raise RuntimeError(
                f"CPPTRAJ failed for {trajectory.name}; "
                f"see {log_path}"
            )

    result: dict[str, np.ndarray] = {}
    counts: dict[str, int] = {}

    for name, path in outputs.items():
        if not path.exists():
            raise RuntimeError(
                f"Expected CPPTRAJ output not found: {path}"
            )

        array = read_numeric(
            path
        )

        if (
            array.ndim != 2
            or array.shape[1] < 2
        ):
            raise RuntimeError(
                f"Could not read SASA series from {path}"
            )

        values = np.asarray(
            array[:, 1],
            dtype=float,
        )

        if np.any(
            ~np.isfinite(values)
        ):
            raise RuntimeError(
                f"Non-finite SASA values found in {path}"
            )

#       if np.any(
#           values < -1.0e-8
#       ):
#           raise RuntimeError(
#               f"Negative SASA values found in {path}"
#           )
#
#       # Avoid retaining tiny negative roundoff if CPPTRAJ ever emits it.
#       values = np.where(
#           values < 0.0,
#           0.0,
#           values,
#       )

        result[name] = values
        counts[name] = len(values)

    if len(set(counts.values())) != 1:
        details = ", ".join(
            f"{name}={count}"
            for name, count in sorted(
                counts.items()
            )
        )

        raise RuntimeError(
            f"Inconsistent frame counts for {trajectory.name}: "
            f"{details}"
        )

    return result


def main() -> None:
    args = parse_args()
    validate_args(args)

    trajectories = expand_trajectories(
        args.trajectory
    )

    args.workdir.mkdir(
        parents=True,
        exist_ok=True,
    )

    # One trajectory stem maps to one work directory and one trajectory label
    # in the final TSV. Reject ambiguous inputs instead of risking collisions.
    labels = [
        trajectory.stem
        for trajectory in trajectories
    ]

    if len(labels) != len(set(labels)):
        duplicates = sorted(
            label
            for label in set(labels)
            if labels.count(label) > 1
        )

        raise RuntimeError(
            "Trajectory filename stems must be unique because each stem is "
            "used as the trajectory label and work-directory name. "
            "Duplicate stem(s): "
            + ", ".join(duplicates)
        )

    sites = [
        (
            safe_name(name),
            mask,
        )
        for name, mask in args.site
    ]

    site_lookup = {
        name: mask
        for name, mask in sites
    }

    # Each worker merely launches/waits for one independent external CPPTRAJ
    # process. The scientific calculation itself is not performed by Python.
    results: dict[
        Path,
        dict[str, np.ndarray],
    ] = {}

    failures: list[
        tuple[Path, BaseException]
    ] = []

    workers = min(
        args.jobs,
        len(trajectories),
    )

    if workers == 1:
        for trajectory in trajectories:
            label = trajectory.stem

            print(
                f"[{label}] running"
            )

            try:
                results[trajectory] = run_one_trajectory(
                    args,
                    trajectory,
                    sites,
                )
            except BaseException as exc:
                failures.append(
                    (
                        trajectory,
                        exc,
                    )
                )
                break

            print(
                f"[{label}] complete"
            )

    else:
        print(
            f"Processing {len(trajectories)} trajectories "
            f"with {workers} concurrent CPPTRAJ jobs."
        )

        with ThreadPoolExecutor(
            max_workers=workers
        ) as executor:
            future_to_trajectory = {
                executor.submit(
                    run_one_trajectory,
                    args,
                    trajectory,
                    sites,
                ): trajectory
                for trajectory in trajectories
            }

            for future in as_completed(
                future_to_trajectory
            ):
                trajectory = (
                    future_to_trajectory[
                        future
                    ]
                )

                label = trajectory.stem

                try:
                    results[trajectory] = (
                        future.result()
                    )
                except BaseException as exc:
                    failures.append(
                        (
                            trajectory,
                            exc,
                        )
                    )
                    print(
                        f"[{label}] FAILED"
                    )
                else:
                    print(
                        f"[{label}] complete"
                    )

    # Do not write a combined table that could be mistaken for a complete
    # dataset if any trajectory failed.
    if failures:
        details = "\n".join(
            (
                f"  {trajectory}: "
                f"{type(exc).__name__}: {exc}"
            )
            for trajectory, exc in failures
        )

        raise RuntimeError(
            f"{len(failures)} trajectory job(s) failed; "
            f"combined output was not written:\n"
            f"{details}"
        )

    rows: list[dict] = []

    # Assemble only after every trajectory succeeds, in the original natural
    # trajectory order. --jobs therefore affects execution only, not output
    # ordering.
    for trajectory in trajectories:
        label = trajectory.stem
        ph = extract_ph(
            label
        )

        series = results[
            trajectory
        ]

        nframes = len(
            next(
                iter(
                    series.values()
                )
            )
        )

        for index in range(
            nframes
        ):
            analysis_row = (
                index + 1
            )

            trajectory_frame = (
                args.start
                + index * args.stride
            )

            if (
                args.first_frame_step
                is None
            ):
                step = ""
            else:
                step = (
                    args.first_frame_step
                    + (
                        trajectory_frame - 1
                    )
                    * args.frame_step_interval
                )

            for name, values in series.items():
                rows.append(
                    {
                        "trajectory": label,
                        "pH": ph,
                        "analysis_row": analysis_row,
                        "trajectory_frame": trajectory_frame,
                        "step": step,
                        "observable": name,
                        "metric": "sasa_lcpo",
                        "selection": site_lookup[
                            name
                        ],
                        "solutemask": args.solutemask,
                        "value": float(
                            values[index]
                        ),
                        "unit": "A2",
                    }
                )

    fields = [
        "trajectory",
        "pH",
        "analysis_row",
        "trajectory_frame",
        "step",
        "observable",
        "metric",
        "selection",
        "solutemask",
        "value",
        "unit",
    ]

    write_tsv(
        args.output,
        rows,
        fields,
    )

    settings = (
        args.output.with_suffix(
            args.output.suffix
            + ".settings.txt"
        )
    )

    with settings.open("w") as fh:
        fh.write(
            f"topology={args.topology.resolve()}\n"
        )

        fh.write(
            "trajectories="
            + ",".join(
                str(path)
                for path in trajectories
            )
            + "\n"
        )

        fh.write(
            f"start={args.start}\n"
        )
        fh.write(
            f"stop={args.stop}\n"
        )
        fh.write(
            f"stride={args.stride}\n"
        )
        fh.write(
            f"first_frame_step={args.first_frame_step}\n"
        )
        fh.write(
            f"frame_step_interval={args.frame_step_interval}\n"
        )
        fh.write(
            f"image_anchor={args.image_anchor}\n"
        )
        fh.write(
            f"solutemask={args.solutemask}\n"
        )
        fh.write(
            f"offset={args.offset:.10g}\n"
        )
        fh.write(
            f"nbrcut={args.nbrcut:.10g}\n"
        )
        fh.write(
            f"cpptraj={args.cpptraj}\n"
        )
        fh.write(
            f"jobs={args.jobs}\n"
        )

        for name, mask in sites:
            fh.write(
                f"site.{name}.mask={mask}\n"
            )

    print(
        f"Wrote {args.output}"
    )
    print(
        f"Wrote {settings}"
    )


if __name__ == "__main__":
    main()
