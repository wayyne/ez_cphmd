#!/usr/bin/env python3
"""
calc_water_wire.py
==================

Frame-level distance-defined water-connectivity analysis for Amber trajectories.

Dependencies
------------
- cpptraj
- numpy
- Python standard library

No MDAnalysis.
No pytraj.

Purpose
-------
Detect short solvent-mediated structural connectivity paths between two solute
endpoint sites:

    A -- W -- B
    A -- W1 -- W2 -- B
    A -- W1 -- W2 -- W3 -- B

Edges are distance-defined:

    endpoint atom -- water O <= endpoint_cutoff
    water O       -- water O <= water_cutoff

Default:
    endpoint_cutoff = 3.4 A
    water_cutoff    = 3.5 A
    max_waters      = 3

This is intentionally a STRUCTURAL water-connectivity metric, not yet a
donor/acceptor/angle-qualified hydrogen-bond wire.

Implementation
--------------
CPPTRAJ writes three synchronized reduced XYZ streams:

    site-A.xyz
        only site-A endpoint atoms

    site-B.xyz
        only site-B endpoint atoms

    water.xyz
        oxygen atoms from the N closest water molecules to the combined dyad

This avoids relying on:
    - PDB atom naming
    - CHARMM/Amber naming conventions
    - retained topology ordering
    - residue numbering in coordinate output

Python reads the three streams frame-by-frame and performs the local graph
analysis.

Parallelism
-----------
--jobs parallelizes independent trajectories.

Output
------
One TSV row per analyzed trajectory frame.

Important columns:
    trajectory
    pH
    trajectory_frame
    step

    min_endpoint_distance_A

    n_local_waters
    n_water_site_a
    n_water_site_b

    has_path_1w
    has_path_2w
    has_path_3w
    connected_le3w

    shortest_path_waters

    n_paths_1w
    n_paths_2w
    n_paths_3w

shortest_path_waters semantics
------------------------------
     0 : endpoint sites themselves are within endpoint_cutoff
     1 : A-W-B
     2 : A-W-W-B
     3 : A-W-W-W-B
    -1 : no water-mediated path with <= max_waters

The has_path_* fields are existence indicators and are not mutually exclusive.
shortest_path_waters gives the mutually exclusive shortest-path class.
"""

from __future__ import annotations

import argparse
import csv
import glob
import math
import os
import re
import shutil
import subprocess
import sys

from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np


PROGRAM_VERSION = "1.2.0"


# =============================================================================
# DATA STRUCTURES
# =============================================================================

@dataclass
class FrameResult:
    trajectory_frame: int

    min_endpoint_distance_A: float

    n_local_waters: int
    n_water_site_a: int
    n_water_site_b: int

    has_path_1w: int
    has_path_2w: int
    has_path_3w: int

    connected_le3w: int
    shortest_path_waters: int

    n_paths_1w: int
    n_paths_2w: int
    n_paths_3w: int


@dataclass
class TrajectoryResult:
    trajectory: Path
    ph: float
    frames: List[FrameResult]


# =============================================================================
# GENERIC HELPERS
# =============================================================================

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

    return sorted(
        uniq,
        key=lambda p: p.name,
    )


def qpath(path: Path) -> str:
    s = str(path.resolve())

    if any(c.isspace() for c in s):
        return '"' + s.replace('"', '\\"') + '"'

    return s


def run_command(
    cmd: Sequence[str],
    *,
    input_text: Optional[str] = None,
    env: Optional[Dict[str, str]] = None,
) -> str:

    proc = subprocess.run(
        list(map(str, cmd)),
        input=input_text,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        env=env,
    )

    if proc.returncode != 0:
        raise RuntimeError(
            "Command failed:\n  "
            + " ".join(map(str, cmd))
            + "\n\n"
            + proc.stdout
        )

    return proc.stdout


def cpptraj_env(threads: int) -> Dict[str, str]:
    env = os.environ.copy()

    env["OMP_NUM_THREADS"] = str(threads)
    env["OMP_DYNAMIC"] = "FALSE"

    # Avoid nested BLAS threading when several trajectory workers run.
    env["OPENBLAS_NUM_THREADS"] = "1"
    env["MKL_NUM_THREADS"] = "1"
    env["NUMEXPR_NUM_THREADS"] = "1"

    return env


def parse_ph(
    path: Path,
    regex: re.Pattern,
) -> float:

    m = (
        regex.search(path.name)
        or regex.search(str(path))
    )

    if not m:
        raise RuntimeError(
            f"Could not parse pH from trajectory {path}; "
            f"regex={regex.pattern!r}"
        )

    return float(
        m.group("ph").replace("_", ".")
    )


def work_tag(path: Path) -> str:
    return re.sub(
        r"[^A-Za-z0-9_.-]+",
        "_",
        path.stem,
    )


# =============================================================================
# CPPTRAJ MASK VALIDATION
# =============================================================================

def selected_atom_ids(
    cpptraj: str,
    topology: Path,
    mask: str,
) -> List[int]:

    text = run_command(
        [
            cpptraj,
            "-p",
            str(topology),
            "-ms",
            mask,
        ]
    )

    ids: List[int] = []

    for line in text.splitlines():

        m = re.search(
            r"Selected\s*=\s*(.*)$",
            line,
        )

        if m:
            ids.extend(
                int(x)
                for x in re.findall(
                    r"\d+",
                    m.group(1),
                )
            )

    if ids:
        return sorted(set(ids))

    # Fallback for cpptraj builds that emit detailed selection tables.
    for line in text.splitlines():

        m = re.match(
            r"\s*(\d+)\s+"
            r"[A-Za-z0-9'+*_-]+\s+"
            r"\d+\s+"
            r"[A-Za-z0-9'+*_-]+\b",
            line,
        )

        if m:
            ids.append(
                int(m.group(1))
            )

    return sorted(set(ids))


def validate_mask_nonempty(
    cpptraj: str,
    topology: Path,
    mask: str,
    context: str,
) -> List[int]:

    ids = selected_atom_ids(
        cpptraj,
        topology,
        mask,
    )

    if not ids:
        raise RuntimeError(
            f"{context}: mask {mask!r} selected zero atoms"
        )

    return ids


# =============================================================================
# XYZ READER
# =============================================================================

def read_xyz_frames(
    path: Path,
) -> Iterable[Tuple[List[str], np.ndarray]]:
    """
    Read a standard multi-frame XYZ trajectory.

    Expected format per frame:

        NATOM
        comment line
        LABEL X Y Z
        LABEL X Y Z
        ...

    Labels are returned but are not used for endpoint identity.
    """

    with path.open(
        "r",
        errors="replace",
    ) as fh:

        frame_index = 0

        while True:

            line = fh.readline()

            if not line:
                break

            line = line.strip()

            if not line:
                continue

            frame_index += 1

            try:
                natoms = int(line)

            except ValueError as exc:
                raise RuntimeError(
                    f"{path}: frame {frame_index}: "
                    f"expected XYZ atom count, got {line!r}"
                ) from exc

            comment = fh.readline()

            if not comment:
                raise RuntimeError(
                    f"{path}: frame {frame_index}: "
                    "truncated after XYZ atom count"
                )

            labels: List[str] = []

            coords = np.empty(
                (natoms, 3),
                dtype=np.float64,
            )

            for i in range(natoms):

                row = fh.readline()

                if not row:
                    raise RuntimeError(
                        f"{path}: frame {frame_index}: "
                        "truncated XYZ frame"
                    )

                fields = row.split()

                if len(fields) < 4:
                    raise RuntimeError(
                        f"{path}: frame {frame_index}: "
                        f"malformed XYZ row: {row!r}"
                    )

                labels.append(
                    fields[0]
                )

                try:
                    coords[i, 0] = float(fields[1])
                    coords[i, 1] = float(fields[2])
                    coords[i, 2] = float(fields[3])

                except ValueError as exc:
                    raise RuntimeError(
                        f"{path}: frame {frame_index}: "
                        f"malformed XYZ coordinates: {row!r}"
                    ) from exc

            yield labels, coords


# =============================================================================
# CPPTRAJ INPUT
# =============================================================================

def build_cpptraj_input(
    args,
    topology: Path,
    traj: Path,
    site_a_out: Path,
    site_b_out: Path,
    water_out: Path,
    closest_out: Path,
) -> str:
    """
    Create three synchronized reduced coordinate streams.

    Important:
        Endpoint identity is defined entirely by the CPPTRAJ masks.
        Python does not infer site identity from output atom names.
    """

    combined_site = (
        f"({args.site_a})|({args.site_b})"
    )

    if args.trajin_stop is None:

        trajin = (
            f"trajin {qpath(traj)} "
            f"{args.trajin_start} last "
            f"{args.trajin_stride}"
        )

    else:

        trajin = (
            f"trajin {qpath(traj)} "
            f"{args.trajin_start} "
            f"{args.trajin_stop} "
            f"{args.trajin_stride}"
        )

    lines = [
        f"parm {qpath(topology)}",
        trajin,
        "",
        f"solvent {args.solvent_mask}",
        "",
    ]

    # =========================================================================
    # SITE A STREAM
    # =========================================================================

    if args.image_anchor:
        lines.append(
            f"autoimage anchor {args.image_anchor}"
        )

    lines.extend([
        f"strip !({args.site_a}) nobox",
        f"outtraj {qpath(site_a_out)} xyz",
        "unstrip",
        "",
    ])

    # =========================================================================
    # SITE B STREAM
    # =========================================================================

    if args.image_anchor:
        lines.append(
            f"autoimage anchor {args.image_anchor}"
        )

    lines.extend([
        f"strip !({args.site_b}) nobox",
        f"outtraj {qpath(site_b_out)} xyz",
        "unstrip",
        "",
    ])

    # =========================================================================
    # LOCAL WATER STREAM
    # =========================================================================

    if args.image_anchor:
        lines.append(
            f"autoimage anchor {args.image_anchor}"
        )

    lines.extend([
        (
            f"closest {args.closest_waters} "
            f"{combined_site} "
            f"oxygen "
            f"closestout {qpath(closest_out)}"
        ),
        f"strip !({args.water_oxygen_mask}) nobox",
        f"outtraj {qpath(water_out)} xyz",
        "",
        "run",
        "quit",
        "",
    ])

    return "\n".join(lines)


# =============================================================================
# DISTANCE ROUTINES
# =============================================================================

def pairwise_distances(
    a: np.ndarray,
    b: np.ndarray,
) -> np.ndarray:
    """
    Return Euclidean distance matrix with shape:

        (len(a), len(b))
    """

    if (
        len(a) == 0
        or len(b) == 0
    ):
        return np.empty(
            (len(a), len(b)),
            dtype=np.float64,
        )

    diff = (
        a[:, None, :]
        -
        b[None, :, :]
    )

    return np.sqrt(
        np.einsum(
            "ijk,ijk->ij",
            diff,
            diff,
        )
    )


def minimum_distance(
    a: np.ndarray,
    b: np.ndarray,
) -> float:

    d = pairwise_distances(
        a,
        b,
    )

    if d.size == 0:
        return math.nan

    return float(
        np.min(d)
    )


# =============================================================================
# GRAPH ANALYSIS
# =============================================================================

def analyze_frame(
    site_a_xyz: np.ndarray,
    site_b_xyz: np.ndarray,
    water_xyz: np.ndarray,
    endpoint_cutoff: float,
    water_cutoff: float,
    max_waters: int,
) -> FrameResult:

    nw = len(water_xyz)

    min_ab = minimum_distance(
        site_a_xyz,
        site_b_xyz,
    )

    # =========================================================================
    # ENDPOINT -> WATER EDGES
    # =========================================================================

    d_a_w = pairwise_distances(
        site_a_xyz,
        water_xyz,
    )

    d_b_w = pairwise_distances(
        site_b_xyz,
        water_xyz,
    )

    if nw:

        a_connected = (
            np.min(
                d_a_w,
                axis=0,
            )
            <= endpoint_cutoff
        )

        b_connected = (
            np.min(
                d_b_w,
                axis=0,
            )
            <= endpoint_cutoff
        )

    else:

        a_connected = np.zeros(
            0,
            dtype=bool,
        )

        b_connected = np.zeros(
            0,
            dtype=bool,
        )

    a_idx = np.flatnonzero(
        a_connected
    )

    b_idx = np.flatnonzero(
        b_connected
    )

    # =========================================================================
    # WATER -> WATER EDGES
    # =========================================================================

    if nw:

        d_ww = pairwise_distances(
            water_xyz,
            water_xyz,
        )

        adjacency = (
            d_ww
            <= water_cutoff
        )

        np.fill_diagonal(
            adjacency,
            False,
        )

    else:

        adjacency = np.zeros(
            (0, 0),
            dtype=bool,
        )

    # =========================================================================
    # ONE-WATER PATH
    #
    # A -- Wi -- B
    # =========================================================================

    n_paths_1w = int(
        np.count_nonzero(
            a_connected
            &
            b_connected
        )
    )

    # =========================================================================
    # TWO-WATER PATH
    #
    # A -- Wi -- Wj -- B
    # =========================================================================

    n_paths_2w = 0

    if max_waters >= 2:

        for i in a_idx:

            n_paths_2w += int(
                np.count_nonzero(
                    adjacency[i, :]
                    &
                    b_connected
                )
            )

    # =========================================================================
    # THREE-WATER PATH
    #
    # A -- Wi -- Wj -- Wk -- B
    #
    # Water nodes must be distinct.
    # =========================================================================

    n_paths_3w = 0

    if max_waters >= 3:

        for i in a_idx:

            middle_candidates = np.flatnonzero(
                adjacency[i, :]
            )

            for j in middle_candidates:

                if j == i:
                    continue

                end_candidates = np.flatnonzero(
                    adjacency[j, :]
                    &
                    b_connected
                )

                for k in end_candidates:

                    if (
                        k == i
                        or
                        k == j
                    ):
                        continue

                    n_paths_3w += 1

    has_1 = int(
        n_paths_1w > 0
    )

    has_2 = int(
        n_paths_2w > 0
    )

    has_3 = int(
        n_paths_3w > 0
    )

    # =========================================================================
    # SHORTEST PATH CLASS
    # =========================================================================

    if (
        np.isfinite(min_ab)
        and
        min_ab <= endpoint_cutoff
    ):

        shortest = 0

    elif has_1:

        shortest = 1

    elif (
        max_waters >= 2
        and has_2
    ):

        shortest = 2

    elif (
        max_waters >= 3
        and has_3
    ):

        shortest = 3

    else:

        shortest = -1

    connected = int(
        has_1
        or has_2
        or has_3
    )

    return FrameResult(
        trajectory_frame=-1,

        min_endpoint_distance_A=min_ab,

        n_local_waters=nw,
        n_water_site_a=len(a_idx),
        n_water_site_b=len(b_idx),

        has_path_1w=has_1,
        has_path_2w=has_2,
        has_path_3w=has_3,

        connected_le3w=connected,

        shortest_path_waters=shortest,

        n_paths_1w=n_paths_1w,
        n_paths_2w=n_paths_2w,
        n_paths_3w=n_paths_3w,
    )


# =============================================================================
# TRAJECTORY WORKER
# =============================================================================

def run_trajectory(
    args,
    topology: Path,
    traj: Path,
    ph_regex: re.Pattern,
    expected_a: int,
    expected_b: int,
) -> TrajectoryResult:

    ph = parse_ph(
        traj,
        ph_regex,
    )

    tag = work_tag(
        traj
    )

    traj_workdir = (
        args.workdir
        /
        tag
    )

    traj_workdir.mkdir(
        parents=True,
        exist_ok=True,
    )

    site_a_out = (
        traj_workdir
        /
        f"{tag}.siteA.xyz"
    )

    site_b_out = (
        traj_workdir
        /
        f"{tag}.siteB.xyz"
    )

    water_out = (
        traj_workdir
        /
        f"{tag}.water.xyz"
    )

    closest_out = (
        traj_workdir
        /
        f"{tag}.closest.dat"
    )

    cpptraj_in = (
        traj_workdir
        /
        f"{tag}.in"
    )

    cpptraj_log = (
        traj_workdir
        /
        f"{tag}.log"
    )

    text = build_cpptraj_input(
        args,
        topology,
        traj,
        site_a_out,
        site_b_out,
        water_out,
        closest_out,
    )

    cpptraj_in.write_text(
        text
    )

    print(
        f"[{traj.stem}] cpptraj preprocessing...",
        flush=True,
    )

    proc = subprocess.run(
        [
            args.cpptraj,
            "-i",
            str(cpptraj_in),
        ],
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        env=cpptraj_env(
            args.cpptraj_threads
        ),
    )

    cpptraj_log.write_text(
        proc.stdout
    )

    explicit_error = any(
        re.match(
            r"\s*Error:",
            line,
            flags=re.IGNORECASE,
        )
        for line in proc.stdout.splitlines()
    )

    if (
        proc.returncode != 0
        or explicit_error
    ):
        raise RuntimeError(
            f"cpptraj failed for {traj}; "
            f"inspect {cpptraj_log}"
        )

    # -------------------------------------------------------------------------
    # Validate coordinate outputs.
    # -------------------------------------------------------------------------

    for label, path in [
        ("site A", site_a_out),
        ("site B", site_b_out),
        ("water", water_out),
    ]:

        if (
            not path.exists()
            or path.stat().st_size == 0
        ):

            raise RuntimeError(
                f"{traj}: missing/empty {label} coordinate output "
                f"{path}; inspect {cpptraj_log}"
            )

    # -------------------------------------------------------------------------
    # Read all three streams in lockstep.
    # -------------------------------------------------------------------------

    a_iter = iter(
        read_xyz_frames(
            site_a_out
        )
    )

    b_iter = iter(
        read_xyz_frames(
            site_b_out
        )
    )

    w_iter = iter(
        read_xyz_frames(
            water_out
        )
    )

    frames: List[FrameResult] = []

    selected_frame_index = 0

    while True:

        try:
            _, site_a_xyz = next(a_iter)
            a_done = False

        except StopIteration:
            site_a_xyz = None
            a_done = True

        try:
            _, site_b_xyz = next(b_iter)
            b_done = False

        except StopIteration:
            site_b_xyz = None
            b_done = True

        try:
            _, water_xyz = next(w_iter)
            w_done = False

        except StopIteration:
            water_xyz = None
            w_done = True

        # All streams ended simultaneously.
        if (
            a_done
            and b_done
            and w_done
        ):
            break

        # Any partial termination is an error.
        if (
            a_done
            or b_done
            or w_done
        ):
            raise RuntimeError(
                f"{traj}: synchronized XYZ streams "
                "have different frame counts"
            )

        selected_frame_index += 1

        assert site_a_xyz is not None
        assert site_b_xyz is not None
        assert water_xyz is not None

        # ---------------------------------------------------------------------
        # Strong per-frame atom-count validation.
        # ---------------------------------------------------------------------

        if len(site_a_xyz) != expected_a:

            raise RuntimeError(
                f"{traj}: frame {selected_frame_index}: "
                f"site-A stream contains {len(site_a_xyz)} atoms; "
                f"expected {expected_a}"
            )

        if len(site_b_xyz) != expected_b:

            raise RuntimeError(
                f"{traj}: frame {selected_frame_index}: "
                f"site-B stream contains {len(site_b_xyz)} atoms; "
                f"expected {expected_b}"
            )

        if len(water_xyz) != args.closest_waters:

            raise RuntimeError(
                f"{traj}: frame {selected_frame_index}: "
                f"water stream contains {len(water_xyz)} atoms; "
                f"expected exactly {args.closest_waters}"
            )

        fr = analyze_frame(
            site_a_xyz,
            site_b_xyz,
            water_xyz,
            args.endpoint_cutoff,
            args.water_cutoff,
            args.max_waters,
        )

        # Map analyzed frame back to the ORIGINAL trajectory frame index.
        trajectory_frame = (
            args.trajin_start
            +
            (
                selected_frame_index
                -
                1
            )
            *
            args.trajin_stride
        )

        fr.trajectory_frame = (
            trajectory_frame
        )

        frames.append(
            fr
        )

    if not frames:

        raise RuntimeError(
            f"{traj}: zero synchronized frames analyzed"
        )

    print(
        f"[{traj.stem}] complete: "
        f"{len(frames)} frames",
        flush=True,
    )

    return TrajectoryResult(
        trajectory=traj,
        ph=ph,
        frames=frames,
    )


# =============================================================================
# OUTPUT
# =============================================================================

def write_output(
    args,
    results: Sequence[TrajectoryResult],
    output: Path,
) -> int:

    fields = [
        "trajectory",
        "pH",
        "analysis_row",
        "trajectory_frame",
        "step",

        "site_a_mask",
        "site_b_mask",

        "endpoint_cutoff_A",
        "water_cutoff_A",

        "closest_waters",
        "max_waters",

        "min_endpoint_distance_A",

        "n_local_waters",
        "n_water_site_a",
        "n_water_site_b",

        "has_path_1w",
        "has_path_2w",
        "has_path_3w",

        "connected_le3w",
        "shortest_path_waters",

        "n_paths_1w",
        "n_paths_2w",
        "n_paths_3w",
    ]

    nrows = 0

    with output.open(
        "w",
        newline="",
    ) as fh:

        writer = csv.DictWriter(
            fh,
            fieldnames=fields,
            delimiter="\t",
            lineterminator="\n",
        )

        writer.writeheader()

        for result in sorted(
            results,
            key=lambda r: (
                r.ph,
                r.trajectory.name,
            ),
        ):

            for fr in result.frames:

                tf = (
                    fr.trajectory_frame
                )

                step = (
                    args.first_frame_step
                    +
                    (tf - 1)
                    *
                    args.frame_step_interval
                )

                writer.writerow({
                    "trajectory":
                        result.trajectory.stem,

                    "pH":
                        f"{result.ph:g}",

                    "analysis_row":
                        tf,

                    "trajectory_frame":
                        tf,

                    "step":
                        step,

                    "site_a_mask":
                        args.site_a,

                    "site_b_mask":
                        args.site_b,

                    "endpoint_cutoff_A":
                        f"{args.endpoint_cutoff:.6g}",

                    "water_cutoff_A":
                        f"{args.water_cutoff:.6g}",

                    "closest_waters":
                        args.closest_waters,

                    "max_waters":
                        args.max_waters,

                    "min_endpoint_distance_A":
                        f"{fr.min_endpoint_distance_A:.8g}",

                    "n_local_waters":
                        fr.n_local_waters,

                    "n_water_site_a":
                        fr.n_water_site_a,

                    "n_water_site_b":
                        fr.n_water_site_b,

                    "has_path_1w":
                        fr.has_path_1w,

                    "has_path_2w":
                        fr.has_path_2w,

                    "has_path_3w":
                        fr.has_path_3w,

                    "connected_le3w":
                        fr.connected_le3w,

                    "shortest_path_waters":
                        fr.shortest_path_waters,

                    "n_paths_1w":
                        fr.n_paths_1w,

                    "n_paths_2w":
                        fr.n_paths_2w,

                    "n_paths_3w":
                        fr.n_paths_3w,
                })

                nrows += 1

    return nrows


def write_settings(
    args,
    topology: Path,
    trajectories: Sequence[Path],
    output: Path,
    nrows: int,
    expected_a: int,
    expected_b: int,
) -> None:

    settings = Path(
        str(output)
        +
        ".settings.txt"
    )

    lines = [
        "program=calc_water_wire.py",
        f"version={PROGRAM_VERSION}",

        f"cpptraj={args.cpptraj}",
        f"topology={topology}",

        f"n_trajectories={len(trajectories)}",

        f"site_a={args.site_a}",
        f"site_b={args.site_b}",

        f"site_a_atom_count={expected_a}",
        f"site_b_atom_count={expected_b}",

        f"solvent_mask={args.solvent_mask}",
        f"water_oxygen_mask={args.water_oxygen_mask}",

        f"endpoint_cutoff_A={args.endpoint_cutoff:g}",
        f"water_cutoff_A={args.water_cutoff:g}",

        f"closest_waters={args.closest_waters}",
        f"max_waters={args.max_waters}",

        f"image_anchor={args.image_anchor or ''}",

        f"first_frame_step={args.first_frame_step}",
        f"frame_step_interval={args.frame_step_interval}",

        f"trajin_start={args.trajin_start}",
        f"trajin_stop={args.trajin_stop}",
        f"trajin_stride={args.trajin_stride}",

        f"output_rows={nrows}",

        "coordinate_interchange=three_synchronized_cpptraj_xyz_streams",

        "site_identity=defined_by_cpptraj_mask_not_coordinate_atom_names",

        "metric_definition=distance_defined_water_connectivity",

        (
            "shortest_path_semantics="
            "0=direct_endpoint_contact,"
            "1=one_water,"
            "2=two_waters,"
            "3=three_waters,"
            "-1=no_path_le_max"
        ),

        (
            "IMPORTANT="
            "Structural distance connectivity; "
            "not donor-angle-qualified hydrogen-bond connectivity."
        ),
    ]

    settings.write_text(
        "\n".join(lines)
        +
        "\n"
    )


# =============================================================================
# TERMINAL SUMMARY
# =============================================================================

def print_summary(
    results: Sequence[TrajectoryResult],
) -> None:

    print()
    print("=" * 92)
    print("WATER-WIRE SUMMARY")
    print("=" * 92)

    for result in sorted(
        results,
        key=lambda r: (
            r.ph,
            r.trajectory.name,
        ),
    ):

        frames = result.frames

        n = len(frames)

        p1 = np.mean(
            [
                f.has_path_1w
                for f in frames
            ]
        )

        p2 = np.mean(
            [
                f.has_path_2w
                for f in frames
            ]
        )

        p3 = np.mean(
            [
                f.has_path_3w
                for f in frames
            ]
        )

        pany = np.mean(
            [
                f.connected_le3w
                for f in frames
            ]
        )

        shortest = np.array(
            [
                f.shortest_path_waters
                for f in frames
            ],
            dtype=int,
        )

        d = np.array(
            [
                f.min_endpoint_distance_A
                for f in frames
            ],
            dtype=float,
        )

        nlocal = np.array(
            [
                f.n_local_waters
                for f in frames
            ],
            dtype=float,
        )

        na = np.array(
            [
                f.n_water_site_a
                for f in frames
            ],
            dtype=float,
        )

        nb = np.array(
            [
                f.n_water_site_b
                for f in frames
            ],
            dtype=float,
        )

        print()
        print(
            f"{result.trajectory.name}  "
            f"pH={result.ph:g}  "
            f"frames={n}"
        )

        print(
            "  any path exists: "
            f"1W={p1:.4f}  "
            f"2W={p2:.4f}  "
            f"3W={p3:.4f}  "
            f"<=3W={pany:.4f}"
        )

        print(
            "  shortest path:   "
            f"1W={np.mean(shortest == 1):.4f}  "
            f"2W={np.mean(shortest == 2):.4f}  "
            f"3W={np.mean(shortest == 3):.4f}  "
            f"none={np.mean(shortest == -1):.4f}"
        )

        print(
            "  endpoint O-O:    "
            f"median={np.nanmedian(d):.3f} A  "
            f"q10={np.nanquantile(d, 0.10):.3f}  "
            f"q90={np.nanquantile(d, 0.90):.3f}"
        )

        print(
            "  local waters:    "
            f"median={np.median(nlocal):.1f}  "
            f"range={int(np.min(nlocal))}-{int(np.max(nlocal))}"
        )

        print(
            "  endpoint shells: "
            f"A mean={np.mean(na):.3f}  "
            f"B mean={np.mean(nb):.3f}"
        )


# =============================================================================
# CLI
# =============================================================================

def make_parser() -> argparse.ArgumentParser:

    p = argparse.ArgumentParser(
        description=(
            "Distance-defined water-wire analysis using CPPTRAJ "
            "coordinate reduction and NumPy graph search."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    p.add_argument(
        "-p",
        "--topology",
        nargs="+",
        required=True,
        help=(
            "Topology path/glob; must resolve to exactly one file."
        ),
    )

    p.add_argument(
        "-t",
        "--trajectory",
        nargs="+",
        required=True,
        help="Trajectory paths/globs.",
    )

    p.add_argument(
        "--site-a",
        required=True,
        help=(
            "Endpoint-A atom mask, e.g. ':36@OE1,OE2'."
        ),
    )

    p.add_argument(
        "--site-b",
        required=True,
        help=(
            "Endpoint-B atom mask, e.g. ':53@OD1,OD2'."
        ),
    )

    p.add_argument(
        "--solvent-mask",
        default=":WAT",
        help="CPPTRAJ solvent molecule mask.",
    )

    p.add_argument(
        "--water-oxygen-mask",
        default=":WAT@O",
        help=(
            "Water oxygen mask retained after closest-solvent selection."
        ),
    )

    p.add_argument(
        "--endpoint-cutoff",
        type=float,
        default=3.4,
        help=(
            "Maximum endpoint-atom/water-O edge distance in Angstrom."
        ),
    )

    p.add_argument(
        "--water-cutoff",
        type=float,
        default=3.5,
        help=(
            "Maximum water-O/water-O edge distance in Angstrom."
        ),
    )

    p.add_argument(
        "--max-waters",
        type=int,
        default=3,
        choices=[
            1,
            2,
            3,
        ],
        help=(
            "Maximum number of intervening waters searched."
        ),
    )

    p.add_argument(
        "--closest-waters",
        type=int,
        default=32,
        help=(
            "Number of nearest waters retained around the combined dyad."
        ),
    )

    p.add_argument(
        "--image-anchor",
        default=None,
        help=(
            "CPPTRAJ autoimage anchor, e.g. '^1' or '^17731'."
        ),
    )

    p.add_argument(
        "--first-frame-step",
        type=int,
        required=True,
        help=(
            "MD step corresponding to original trajectory frame 1."
        ),
    )

    p.add_argument(
        "--frame-step-interval",
        type=int,
        required=True,
        help=(
            "MD-step increment between consecutive original trajectory frames."
        ),
    )

    p.add_argument(
        "--trajin-start",
        type=int,
        default=1,
        help="First original trajectory frame read.",
    )

    p.add_argument(
        "--trajin-stop",
        type=int,
        default=None,
        help=(
            "Last original trajectory frame read; default is last."
        ),
    )

    p.add_argument(
        "--trajin-stride",
        type=int,
        default=1,
        help="Original trajectory-frame stride.",
    )

    p.add_argument(
        "--ph-regex",
        default=r"ph(?P<ph>\d+(?:[_\.]\d+)?)",
        help=(
            "Regex containing named group (?P<ph>...)."
        ),
    )

    p.add_argument(
        "--cpptraj",
        default="cpptraj",
        help="CPPTRAJ executable.",
    )

    p.add_argument(
        "--cpptraj-threads",
        type=int,
        default=1,
        help="OMP threads per CPPTRAJ worker.",
    )

    p.add_argument(
        "--jobs",
        type=int,
        default=1,
        help=(
            "Number of independent trajectories processed concurrently."
        ),
    )

    p.add_argument(
        "--workdir",
        type=Path,
        default=Path(
            "water_wire_work"
        ),
        help="Intermediate CPPTRAJ working directory.",
    )

    p.add_argument(
        "-o",
        "--output",
        type=Path,
        required=True,
        help="Frame-level TSV output.",
    )

    p.add_argument(
        "--keep-intermediates",
        action="store_true",
        help=(
            "Keep synchronized reduced XYZ files after analysis."
        ),
    )

    p.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {PROGRAM_VERSION}",
    )

    return p


# =============================================================================
# MAIN
# =============================================================================

def main(
    argv: Optional[Sequence[str]] = None,
) -> int:

    args = (
        make_parser()
        .parse_args(argv)
    )

    # -------------------------------------------------------------------------
    # Validate scalar arguments.
    # -------------------------------------------------------------------------

    if args.jobs < 1:
        raise RuntimeError(
            "--jobs must be >= 1"
        )

    if args.cpptraj_threads < 1:
        raise RuntimeError(
            "--cpptraj-threads must be >= 1"
        )

    if args.closest_waters < 3:
        raise RuntimeError(
            "--closest-waters must be >= 3"
        )

    if args.endpoint_cutoff <= 0:
        raise RuntimeError(
            "--endpoint-cutoff must be > 0"
        )

    if args.water_cutoff <= 0:
        raise RuntimeError(
            "--water-cutoff must be > 0"
        )

    if args.first_frame_step < 0:
        raise RuntimeError(
            "--first-frame-step must be >= 0"
        )

    if args.frame_step_interval <= 0:
        raise RuntimeError(
            "--frame-step-interval must be > 0"
        )

    if args.trajin_start < 1:
        raise RuntimeError(
            "--trajin-start must be >= 1"
        )

    if args.trajin_stride < 1:
        raise RuntimeError(
            "--trajin-stride must be >= 1"
        )

    if (
        args.trajin_stop is not None
        and args.trajin_stop < args.trajin_start
    ):
        raise RuntimeError(
            "--trajin-stop must be >= --trajin-start"
        )

    if shutil.which(
        args.cpptraj
    ) is None:
        raise RuntimeError(
            f"CPPTRAJ executable not found: {args.cpptraj}"
        )

    # -------------------------------------------------------------------------
    # Resolve topology/trajectories.
    # -------------------------------------------------------------------------

    topologies = expand_paths(
        args.topology
    )

    if len(topologies) != 1:
        raise RuntimeError(
            "Topology arguments must resolve to exactly one file; "
            f"found {len(topologies)}:\n  "
            +
            "\n  ".join(
                map(str, topologies)
            )
        )

    topology = topologies[0]

    trajectories = expand_paths(
        args.trajectory
    )

    if not trajectories:
        raise RuntimeError(
            "No trajectories matched -t/--trajectory"
        )

    # -------------------------------------------------------------------------
    # Parse pH regex.
    # -------------------------------------------------------------------------

    try:
        ph_regex = re.compile(
            args.ph_regex,
            flags=re.IGNORECASE,
        )

    except re.error as exc:
        raise RuntimeError(
            f"Invalid --ph-regex: {exc}"
        ) from exc

    if (
        "ph"
        not in ph_regex.groupindex
    ):
        raise RuntimeError(
            "--ph-regex must contain named group (?P<ph>...)"
        )

    # -------------------------------------------------------------------------
    # Validate atom masks directly against topology.
    # -------------------------------------------------------------------------

    a_ids = validate_mask_nonempty(
        args.cpptraj,
        topology,
        args.site_a,
        "site A",
    )

    b_ids = validate_mask_nonempty(
        args.cpptraj,
        topology,
        args.site_b,
        "site B",
    )

    validate_mask_nonempty(
        args.cpptraj,
        topology,
        args.water_oxygen_mask,
        "water oxygen",
    )

    overlap = (
        set(a_ids)
        &
        set(b_ids)
    )

    if overlap:
        raise RuntimeError(
            "Site A and B masks overlap: "
            f"{sorted(overlap)}"
        )

    print(
        f"Site A atoms : {len(a_ids)}"
    )

    print(
        f"Site B atoms : {len(b_ids)}"
    )

    print(
        f"Trajectories : {len(trajectories)}"
    )

    print(
        f"Closest waters retained : {args.closest_waters}"
    )

    print(
        f"Endpoint cutoff : {args.endpoint_cutoff:g} A"
    )

    print(
        f"Water-water cutoff : {args.water_cutoff:g} A"
    )

    print(
        f"Max intervening waters : {args.max_waters}"
    )

    # -------------------------------------------------------------------------
    # Prepare directories.
    # -------------------------------------------------------------------------

    args.workdir.mkdir(
        parents=True,
        exist_ok=True,
    )

    args.output.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    # -------------------------------------------------------------------------
    # Parallel trajectory processing.
    # -------------------------------------------------------------------------

    results: List[
        TrajectoryResult
    ] = []

    errors: List[str] = []

    with ThreadPoolExecutor(
        max_workers=args.jobs
    ) as executor:

        futures = {
            executor.submit(
                run_trajectory,
                args,
                topology,
                traj,
                ph_regex,
                len(a_ids),
                len(b_ids),
            ):
            traj

            for traj in trajectories
        }

        for fut in as_completed(
            futures
        ):

            traj = futures[fut]

            try:
                results.append(
                    fut.result()
                )

            except Exception as exc:
                errors.append(
                    f"{traj}:\n{exc}"
                )

    if errors:
        raise RuntimeError(
            "One or more trajectories failed:\n\n"
            +
            "\n\n".join(errors)
        )

    # -------------------------------------------------------------------------
    # Output.
    # -------------------------------------------------------------------------

    nrows = write_output(
        args,
        results,
        args.output,
    )

    write_settings(
        args,
        topology,
        trajectories,
        args.output,
        nrows,
        len(a_ids),
        len(b_ids),
    )

    # -------------------------------------------------------------------------
    # Optional cleanup.
    # -------------------------------------------------------------------------

    if not args.keep_intermediates:

        for result in results:

            tag = work_tag(
                result.trajectory
            )

            traj_workdir = (
                args.workdir
                /
                tag
            )

            for suffix in [
                "siteA.xyz",
                "siteB.xyz",
                "water.xyz",
            ]:

                p = (
                    traj_workdir
                    /
                    f"{tag}.{suffix}"
                )

                try:
                    p.unlink()

                except FileNotFoundError:
                    pass

    # -------------------------------------------------------------------------
    # Summary.
    # -------------------------------------------------------------------------

    print_summary(
        results
    )

    print()
    print("=" * 92)

    print(
        f"Trajectories analyzed : {len(results)}"
    )

    print(
        f"Frame rows written    : {nrows}"
    )

    print(
        f"Wrote                 : {args.output}"
    )

    print(
        f"Wrote settings        : {args.output}.settings.txt"
    )

    print("=" * 92)

    return 0


if __name__ == "__main__":

    try:
        raise SystemExit(
            main()
        )

    except KeyboardInterrupt:
        print(
            "Interrupted",
            file=sys.stderr,
        )
        raise SystemExit(130)

    except Exception as exc:
        print(
            f"ERROR: {exc}",
            file=sys.stderr,
        )
        raise SystemExit(1)
