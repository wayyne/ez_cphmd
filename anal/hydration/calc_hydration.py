#!/usr/bin/env python3
"""
calc_hydration.py
=================
General first-shell hydration observable generator for Amber trajectories.

The primary metric reproduces the Shen-lab water-shell definition used in
Peeples, Liu, and Shen (J. Phys. Chem. B 2024): the number of water molecules
whose selected solvent atom(s) lie within 3.4 A of any atom in the specified
solute-site mask. With a water-oxygen mask, this is a unique-water union count
around the complete site mask (e.g. both carboxylate oxygens).

Design:
* State-agnostic: no CpH lambda/state information is required.
* One cpptraj pass per trajectory for all requested hydration sites.
* Uses cpptraj `watershell`; periodic imaging is enabled by cpptraj by default.
* Long-format frame-level output compatible with partition_by_state.py.
* analysis_row is frame-level (identical to trajectory_frame), not metric-row-level.
"""

from __future__ import annotations

import argparse
import csv
import glob
import os
import re
import shutil
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

PROGRAM_VERSION = "1.0.0"


@dataclass(frozen=True)
class SiteSpec:
    observable: str
    site_mask: str


@dataclass
class TrajectoryResult:
    trajectory: Path
    ph: float
    nframes: int
    data: Dict[str, List[float]]


def expand_paths(items: Sequence[str]) -> List[Path]:
    out: List[Path] = []
    for item in items:
        hits = [Path(p) for p in glob.glob(item)]
        if hits:
            out.extend(hits)
        else:
            p = Path(item)
            if p.exists():
                out.append(p)
    seen = set()
    uniq: List[Path] = []
    for p in out:
        rp = p.resolve()
        if rp not in seen:
            seen.add(rp)
            uniq.append(rp)
    return sorted(uniq, key=lambda p: p.name)


def run_command(cmd: Sequence[str], input_text: Optional[str] = None, env=None) -> str:
    proc = subprocess.run(
        list(map(str, cmd)), input=input_text, text=True,
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT, env=env,
    )
    if proc.returncode != 0:
        raise RuntimeError(
            "Command failed:\n  " + " ".join(map(str, cmd)) + "\n\n" + proc.stdout
        )
    return proc.stdout


def qpath(path: Path) -> str:
    s = str(path.resolve())
    if any(c.isspace() for c in s):
        return '"' + s.replace('"', '\\"') + '"'
    return s


def parse_ph(path: Path, regex: re.Pattern) -> float:
    m = regex.search(path.name) or regex.search(str(path))
    if not m:
        raise RuntimeError(
            f"Could not parse pH from trajectory {path}; regex={regex.pattern!r}"
        )
    return float(m.group("ph").replace("_", "."))


def selected_atom_ids(cpptraj: str, topology: Path, mask: str) -> List[int]:
    text = run_command([cpptraj, "-p", str(topology), "-ms", mask])
    ids: List[int] = []
    for line in text.splitlines():
        m = re.search(r"Selected\s*=\s*(.*)$", line)
        if m:
            ids.extend(int(x) for x in re.findall(r"\d+", m.group(1)))
    if ids:
        return sorted(set(ids))

    # Some cpptraj builds print a detailed table rather than Selected=.
    for line in text.splitlines():
        m = re.match(r"\s*(\d+)\s+[A-Za-z0-9'+*_-]+\s+\d+\s+[A-Za-z0-9'+*_-]+\b", line)
        if m:
            ids.append(int(m.group(1)))
    return sorted(set(ids))


def validate_mask_nonempty(cpptraj: str, topology: Path, mask: str, context: str) -> None:
    ids = selected_atom_ids(cpptraj, topology, mask)
    if not ids:
        raise RuntimeError(f"{context}: mask {mask!r} selected zero atoms")


def cpptraj_env(threads: int) -> Dict[str, str]:
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = str(threads)
    env["OMP_DYNAMIC"] = "FALSE"
    env["OPENBLAS_NUM_THREADS"] = "1"
    env["MKL_NUM_THREADS"] = "1"
    env["NUMEXPR_NUM_THREADS"] = "1"
    return env


def work_tag(path: Path) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", path.stem)


def site_tag(name: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", name)


def build_cpptraj_input(
    args, topology: Path, traj: Path, sites: Sequence[SiteSpec], workdir: Path
) -> Tuple[str, Dict[str, Path]]:
    tag = work_tag(traj)
    outputs: Dict[str, Path] = {}
    lines = [f"parm {qpath(topology)}", f"trajin {qpath(traj)}"]
    if args.image_anchor:
        lines.append(f"autoimage anchor {args.image_anchor}")
    lines.append("")

    for site in sites:
        out = workdir / f"{tag}.{site_tag(site.observable)}.watershell.dat"
        outputs[site.observable] = out
        # watershell itself performs periodic imaging unless 'noimage' is given.
        lines.append(
            f"watershell {site.site_mask} out {qpath(out)} "
            f"lower {args.cutoff:.10g} upper {args.upper_cutoff:.10g} {args.water_mask}"
        )

    lines.extend(["", "run", "quit", ""])
    return "\n".join(lines), outputs


def read_watershell_first_shell(path: Path) -> List[float]:
    """Read the first-shell count from cpptraj watershell output.

    cpptraj watershell output contains a frame/index column followed by the
    lower-cutoff (first-shell) and upper-cutoff datasets. We retain only the
    first-shell count because that is the historical 3.4-A hydration metric.
    """
    vals: List[float] = []
    with path.open("r", errors="replace") as fh:
        for line in fh:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            numeric: List[float] = []
            for tok in s.split():
                try:
                    numeric.append(float(tok))
                except ValueError:
                    pass
            if len(numeric) >= 2:
                vals.append(numeric[1])
    if not vals:
        raise RuntimeError(f"No first-shell watershell values read from {path}")
    return vals


def run_trajectory(
    args, topology: Path, traj: Path, ph_regex: re.Pattern,
    sites: Sequence[SiteSpec], workdir: Path,
) -> TrajectoryResult:
    ph = parse_ph(traj, ph_regex)
    print(f"[{traj.stem}] running", flush=True)
    workdir.mkdir(parents=True, exist_ok=True)
    text, outputs = build_cpptraj_input(args, topology, traj, sites, workdir)
    inp = workdir / f"{work_tag(traj)}.in"
    log = workdir / f"{work_tag(traj)}.log"
    inp.write_text(text)

    proc = subprocess.run(
        [args.cpptraj, "-i", str(inp)],
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True,
        env=cpptraj_env(args.cpptraj_threads),
    )
    log.write_text(proc.stdout)
    explicit_error = any(
        re.match(r"\s*Error:", line, flags=re.IGNORECASE)
        for line in proc.stdout.splitlines()
    )
    if proc.returncode != 0 or explicit_error:
        raise RuntimeError(f"cpptraj failed for {traj}; inspect {log}")

    data: Dict[str, List[float]] = {}
    lengths: List[int] = []
    for site in sites:
        out = outputs[site.observable]
        if not out.exists() or out.stat().st_size == 0:
            raise RuntimeError(f"Missing/empty cpptraj output {out}; inspect {log}")
        series = read_watershell_first_shell(out)
        data[site.observable] = series
        lengths.append(len(series))

    if len(set(lengths)) != 1:
        raise RuntimeError(f"Hydration series lengths differ for {traj}: {lengths}")

    print(f"[{traj.stem}] complete", flush=True)
    return TrajectoryResult(traj, ph, lengths[0], data)


def float_text(x: float) -> str:
    return f"{x:.10g}"


def write_output(
    args, results: Sequence[TrajectoryResult], sites: Sequence[SiteSpec], output: Path
) -> int:
    fields = [
        "trajectory", "pH", "analysis_row", "trajectory_frame", "step",
        "observable", "site_mask", "water_mask", "cutoff_A",
        "metric", "value", "unit",
    ]

    output_rows = 0
    with output.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", lineterminator="\n")
        w.writeheader()
        for result in sorted(results, key=lambda r: (r.ph, r.trajectory.name)):
            for frame0 in range(result.nframes):
                trajectory_frame = frame0 + 1
                analysis_row = trajectory_frame
                step = args.first_frame_step + frame0 * args.frame_step_interval
                for site in sites:
                    value = result.data[site.observable][frame0]
                    w.writerow({
                        "trajectory": result.trajectory.stem,
                        "pH": f"{result.ph:g}",
                        "analysis_row": analysis_row,
                        "trajectory_frame": trajectory_frame,
                        "step": step,
                        "observable": site.observable,
                        "site_mask": site.site_mask,
                        "water_mask": args.water_mask,
                        "cutoff_A": float_text(args.cutoff),
                        "metric": "first_shell_water_count",
                        "value": float_text(value),
                        "unit": "water",
                    })
                    output_rows += 1
    return output_rows


def write_settings(
    args, topology: Path, trajectories: Sequence[Path], sites: Sequence[SiteSpec],
    output: Path, nrows: int,
) -> None:
    p = Path(str(output) + ".settings.txt")
    lines = [
        "program=calc_hydration.py",
        f"version={PROGRAM_VERSION}",
        f"cpptraj={args.cpptraj}",
        f"topology={topology}",
        f"n_trajectories={len(trajectories)}",
        f"n_sites={len(sites)}",
        f"first_frame_step={args.first_frame_step}",
        f"frame_step_interval={args.frame_step_interval}",
        f"image_anchor={args.image_anchor or ''}",
        f"water_mask={args.water_mask}",
        f"first_shell_cutoff_A={args.cutoff:g}",
        f"cpptraj_upper_cutoff_A={args.upper_cutoff:g}",
        f"output_rows={nrows}",
        "metric=first_shell_water_count",
        "count_semantics=unique solvent molecules within cutoff of any atom in site_mask",
        "historical_definition=Peeples_Liu_Shen_JPCB_2024_3.4A_first_shell",
        "",
        "sites:",
    ]
    for s in sites:
        lines.append(f"{s.observable}\t{s.site_mask}")
    p.write_text("\n".join(lines) + "\n")


def make_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="General Amber/cpptraj first-shell hydration analyzer.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("-p", "--topology", nargs="+", required=True,
                   help="Topology path/glob. Must resolve to exactly one file.")
    p.add_argument("-t", "--trajectory", nargs="+", required=True,
                   help="Trajectory paths/globs.")
    p.add_argument(
        "--site", nargs=2, action="append", metavar=("LABEL", "MASK"), required=True,
        help="Hydration site mask; repeatable. Use a union mask for unique-water counts, "
             "e.g. --site E35 ':36@OE1,OE2'.",
    )
    p.add_argument("--water-mask", default=":WAT@O",
                   help="Solvent atom mask used by watershell; normally water oxygens.")
    p.add_argument("--cutoff", type=float, default=3.4,
                   help="First-shell cutoff in Angstrom; 3.4 A reproduces Shen-lab analyses.")
    p.add_argument("--upper-cutoff", type=float, default=5.0,
                   help="cpptraj watershell upper cutoff; upper-shell values are not emitted.")
    p.add_argument("--image-anchor", default=None,
                   help="Optional cpptraj autoimage anchor mask, e.g. '^1'.")
    p.add_argument("--first-frame-step", type=int, required=True,
                   help="MD step associated with trajectory frame 1.")
    p.add_argument("--frame-step-interval", type=int, required=True,
                   help="MD step increment between consecutive trajectory frames.")
    p.add_argument(
        "--ph-regex", default=r"ph(?P<ph>\d+(?:[_\.]\d+)?)",
        help="Regex containing named group 'ph' for parsing pH from trajectory name.",
    )
    p.add_argument("--cpptraj", default="cpptraj", help="cpptraj executable.")
    p.add_argument("--cpptraj-threads", type=int, default=1,
                   help="OMP threads per cpptraj worker.")
    p.add_argument("--jobs", type=int, default=1,
                   help="Number of trajectories processed concurrently.")
    p.add_argument("--workdir", type=Path, default=Path("hydration_work"),
                   help="Directory for cpptraj inputs/logs/intermediate series.")
    p.add_argument("-o", "--output", type=Path, required=True,
                   help="Long-format TSV output.")
    p.add_argument("--version", action="version", version=f"%(prog)s {PROGRAM_VERSION}")
    return p


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = make_parser().parse_args(argv)
    if args.jobs < 1:
        raise RuntimeError("--jobs must be >=1")
    if args.cpptraj_threads < 1:
        raise RuntimeError("--cpptraj-threads must be >=1")
    if args.frame_step_interval <= 0:
        raise RuntimeError("--frame-step-interval must be >0")
    if args.cutoff <= 0:
        raise RuntimeError("--cutoff must be >0")
    if args.upper_cutoff <= args.cutoff:
        raise RuntimeError("--upper-cutoff must be greater than --cutoff")
    if shutil.which(args.cpptraj) is None:
        raise RuntimeError(f"cpptraj executable not found: {args.cpptraj}")

    topologies = expand_paths(args.topology)
    if len(topologies) != 1:
        raise RuntimeError(
            f"Topology arguments must resolve to exactly one file; found {len(topologies)}: "
            + ", ".join(map(str, topologies))
        )
    topology = topologies[0]

    trajectories = expand_paths(args.trajectory)
    if not trajectories:
        raise RuntimeError("No trajectories matched -t/--trajectory")

    try:
        ph_regex = re.compile(args.ph_regex, flags=re.IGNORECASE)
    except re.error as exc:
        raise RuntimeError(f"Invalid --ph-regex: {exc}") from exc
    if "ph" not in ph_regex.groupindex:
        raise RuntimeError("--ph-regex must contain a named capture group (?P<ph>...)")

    labels = [x[0] for x in args.site]
    if len(set(labels)) != len(labels):
        raise RuntimeError("Duplicate --site labels detected; labels must be unique")
    sites = [SiteSpec(label, mask) for label, mask in args.site]

    for site in sites:
        validate_mask_nonempty(args.cpptraj, topology, site.site_mask, f"site {site.observable}")
    validate_mask_nonempty(args.cpptraj, topology, args.water_mask, "water mask")

    args.workdir.mkdir(parents=True, exist_ok=True)
    args.output.parent.mkdir(parents=True, exist_ok=True)

    results: List[TrajectoryResult] = []
    errors: List[str] = []
    with ThreadPoolExecutor(max_workers=args.jobs) as ex:
        futs = {
            ex.submit(run_trajectory, args, topology, traj, ph_regex, sites, args.workdir): traj
            for traj in trajectories
        }
        for fut in as_completed(futs):
            traj = futs[fut]
            try:
                results.append(fut.result())
            except Exception as exc:
                errors.append(f"{traj}: {exc}")
    if errors:
        raise RuntimeError("One or more trajectories failed:\n  " + "\n  ".join(errors))

    nframes = sum(r.nframes for r in results)
    nrows = write_output(args, results, sites, args.output)
    write_settings(args, topology, trajectories, sites, args.output, nrows)

    print()
    print(f"Trajectories analyzed       : {len(results)}")
    print(f"Hydration sites             : {len(sites)}")
    print(f"Trajectory frames           : {nframes}")
    print(f"Frame-wise output rows      : {nrows}")
    print(f"Program version             : {PROGRAM_VERSION}")
    print(f"First-shell cutoff          : {args.cutoff:g} A")
    print(f"Water mask                  : {args.water_mask}")
    print("Metric                      : first_shell_water_count")
    print(f"Wrote {args.output}")
    print(f"Wrote {args.output}.settings.txt")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except KeyboardInterrupt:
        print("Interrupted", file=sys.stderr)
        raise SystemExit(130)
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1)
