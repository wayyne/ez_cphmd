#!/usr/bin/env python3
"""
calc_water_wire.py
==================

Lambda-aware hydrogen-bonded water-wire analysis for Amber continuous CpHMD.

Dependencies
------------
- cpptraj
- numpy
- Python standard library

No MDAnalysis.
No pytraj.

Purpose
-------
Detect genuine simultaneous H-bond paths between two CpHMD endpoint sites:

    A -- W -- B
    A -- W1 -- W2 -- B
    A -- W1 -- W2 -- W3 -- B

A "real" H-bond edge requires:

    donor-heavy ... acceptor distance <= --distance  (default 3.0 A)
    donor-H ... acceptor angle       >= --angle     (default 135 deg)

Every edge in a path must satisfy the criterion in the same trajectory frame.

CpH chemistry
-------------
This script follows the same Asp/Glu convention as calc_hbond_v5_1.py:

    itauto=3 : protonation coordinate
      lambda <= --lambda-low   -> protonated
      lambda >= --lambda-high  -> deprotonated
      between thresholds       -> mixed/ambiguous

    itauto=4 : tautomer coordinate (only needed when protonated)
      lambda <= --lambda-low   -> O2-H
      lambda >= --lambda-high  -> O1-H
      between thresholds       -> mixed/ambiguous

For a deprotonated Asp/Glu, both carboxylate oxygens are acceptors.
For a protonated Asp/Glu, the protonated oxygen is a donor and the other
oxygen is an acceptor.

The endpoint masks must each select exactly the two carboxylate oxygens of one
Asp/Glu site, e.g.

    --site-a ':36@OE1,OE2'
    --site-b ':53@OD1,OD2'

The script expands those masks with CPPTRAJ, determines whether each endpoint
is Asp or Glu, and constructs explicit O1/O2 and H1/H2 masks from topology
residue identity.

Water handling
--------------
CPPTRAJ retains the N closest waters to the combined endpoint site. Separate,
synchronized XYZ streams are written for:

    endpoint A O1
    endpoint A O2
    endpoint A H1
    endpoint A H2
    endpoint B O1
    endpoint B O2
    endpoint B H1
    endpoint B H2
    local water O/H triplets (O,H1,H2) from the same post-closest topology

This avoids dependence on PDB naming and avoids MDAnalysis/pytraj.

For each frame, retained waters are read as topology-order O,H1,H2 triplets.
Every O-H assignment is validated geometrically against --water-oh-max
(default 1.25 A). This works for 3-site waters and 4-site waters such as OPC
because virtual/extra sites are excluded from the retained stream.

Distance-only candidate graph
-----------------------------
For diagnostic comparison, the script also retains the previously validated
distance-only structural-connectivity metric:

    endpoint O -- water O <= --candidate-endpoint-cutoff (default 3.4 A)
    water O    -- water O <= --candidate-water-cutoff    (default 3.5 A)

The H-bond-wire metric is primary. Candidate connectivity is diagnostic.

Frame/lambda synchronization
----------------------------
Analyzed trajectory frame k maps exactly to:

    step = first_frame_step
           + (trajectory_frame - 1) * frame_step_interval

Lambda rows are matched by exact MD step. Missing steps inside the lambda-file
range are fatal. Frames outside the lambda range are retained with
chemistry_status=outside_lambda_window and blank H-bond-wire fields.

Output
------
One TSV row per analyzed trajectory frame.

Primary fields:

    chemistry_status
    hbond_wire_1w
    hbond_wire_2w
    hbond_wire_3w
    hbond_connected_le3w
    hbond_shortest_path_waters

shortest-path semantics:
     1 = A-W-B
     2 = A-W-W-B
     3 = A-W-W-W-B
    -1 = no H-bonded path with <= max_waters
    blank = mixed/outside lambda chemistry

The has-wire fields are not mutually exclusive. shortest path is.

Representative path identity
----------------------------
For each shortest H-bond path, shortest_path_water_slots stores the local
water-slot indices (0-based within the retained closest-water set), e.g.

    4
    2,7
    1,5,9

These are sufficient to validate path geometry within the reduced streams.
They are frame-local slots, not persistent topology residue IDs.
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


PROGRAM_VERSION = "2.1.0"


# =============================================================================
# DATA STRUCTURES
# =============================================================================

@dataclass(frozen=True)
class TopologyAtom:
    index: int
    name: str
    resid: int
    resname: str


@dataclass(frozen=True)
class LambdaFile:
    path: Path
    ires: Tuple[str, ...]
    itauto: Tuple[str, ...]
    step_to_values: Dict[int, Tuple[float, ...]]
    min_step: int
    max_step: int


@dataclass(frozen=True)
class EndpointSpec:
    user_mask: str
    resid: int
    resname: str
    family: str
    o1_name: str
    o2_name: str
    h1_name: str
    h2_name: str
    o1_mask: str
    o2_mask: str
    h1_mask: str
    h2_mask: str


@dataclass
class EndpointRoles:
    status: str
    protonated: Optional[bool]
    donor_o: Tuple[int, ...]
    acceptor_o: Tuple[int, ...]
    donor_h_for_o: Dict[int, int]


@dataclass
class FrameResult:
    trajectory_frame: int
    step: int
    chemistry_status: str

    min_endpoint_distance_A: float
    n_local_waters: int
    n_water_site_a_candidate: int
    n_water_site_b_candidate: int

    candidate_has_path_1w: int
    candidate_has_path_2w: int
    candidate_has_path_3w: int
    candidate_connected_le3w: int
    candidate_shortest_path_waters: int

    hbond_wire_1w: Optional[int]
    hbond_wire_2w: Optional[int]
    hbond_wire_3w: Optional[int]
    hbond_connected_le3w: Optional[int]
    hbond_shortest_path_waters: Optional[int]

    n_hbond_paths_1w: Optional[int]
    n_hbond_paths_2w: Optional[int]
    n_hbond_paths_3w: Optional[int]

    shortest_path_water_slots: str


@dataclass
class TrajectoryResult:
    trajectory: Path
    ph: float
    lambda_file: Path
    frames: List[FrameResult]


# =============================================================================
# GENERIC HELPERS
# =============================================================================

def natural_key(value: str | Path):
    return [
        int(x) if x.isdigit() else x.lower()
        for x in re.split(r"([0-9]+)", str(value))
    ]


def expand_paths(items: Sequence[str], what: str) -> List[Path]:
    found: Dict[str, Path] = {}

    for item in items:
        matches = glob.glob(item)
        if matches:
            for match in matches:
                p = Path(match).resolve()
                if p.is_file():
                    found[str(p)] = p
            continue

        p = Path(item).resolve()
        if p.is_file():
            found[str(p)] = p
        else:
            raise FileNotFoundError(f"No {what} file matched: {item}")

    out = sorted(found.values(), key=natural_key)
    if not out:
        raise RuntimeError(f"No {what} files found")
    return out


def qpath(path: Path) -> str:
    return '"' + str(path.resolve()).replace('"', r'\"') + '"'


def qmask(mask: str) -> str:
    return '"' + mask.replace('"', r'\"') + '"'


def cpptraj_env(threads: int) -> Dict[str, str]:
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = str(threads)
    env["OMP_DYNAMIC"] = "FALSE"
    env["OPENBLAS_NUM_THREADS"] = "1"
    env["MKL_NUM_THREADS"] = "1"
    env["NUMEXPR_NUM_THREADS"] = "1"
    return env


def extract_ph(name: str) -> float:
    m = re.search(r"[pP][hH][_=-]?([0-9]+(?:[._][0-9]+)?)", name)
    return float(m.group(1).replace("_", ".")) if m else math.nan


def safe_component(value: str) -> str:
    text = re.sub(r"[^A-Za-z0-9_.-]+", "_", value.strip()).strip("._")
    return text or "item"


# =============================================================================
# TOPOLOGY MASK EXPANSION
# =============================================================================

def cpptraj_mask_atoms(
    topology: Path,
    mask: str,
    cpptraj: str,
) -> Tuple[TopologyAtom, ...]:
    """
    Expand an Amber mask with CPPTRAJ itself.
    """

    proc = subprocess.run(
        [cpptraj, "-p", str(topology.resolve()), "--mask", mask],
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )

    if proc.returncode != 0:
        raise RuntimeError(
            f"CPPTRAJ could not expand mask {mask!r}:\n"
            f"{proc.stdout}\n{proc.stderr}"
        )

    atoms: List[TopologyAtom] = []

    for text in (proc.stdout, proc.stderr):
        for raw in text.splitlines():
            parts = raw.split()
            if len(parts) < 4:
                continue
            try:
                index = int(parts[0])
                resid = int(parts[2])
            except ValueError:
                continue

            atoms.append(
                TopologyAtom(
                    index=index,
                    name=parts[1].upper(),
                    resid=resid,
                    resname=parts[3].upper(),
                )
            )

        if atoms:
            break

    # Remove duplicates while preserving topology order.
    seen = set()
    uniq = []
    for a in atoms:
        if a.index not in seen:
            seen.add(a.index)
            uniq.append(a)

    return tuple(sorted(uniq, key=lambda a: a.index))


def canonical_family(resname: str) -> Optional[str]:
    name = resname.upper()
    if name in {"ASP", "ASH", "AS2", "AS4"}:
        return "ASP"
    if name in {"GLU", "GLH", "GL2", "GL4"}:
        return "GLU"
    return None


def build_endpoint_spec(
    topology: Path,
    mask: str,
    cpptraj: str,
) -> EndpointSpec:
    atoms = cpptraj_mask_atoms(topology, mask, cpptraj)

    if len(atoms) != 2:
        raise RuntimeError(
            f"Endpoint mask {mask!r} must select exactly two carboxylate "
            f"oxygen atoms; selected {len(atoms)}: {atoms}"
        )

    residues = {(a.resid, a.resname) for a in atoms}
    if len(residues) != 1:
        raise RuntimeError(
            f"Endpoint mask {mask!r} spans multiple residues: {residues}"
        )

    resid, resname = next(iter(residues))
    family = canonical_family(resname)

    if family == "GLU":
        o1, o2 = "OE1", "OE2"
        h1, h2 = "HE1", "HE2"
    elif family == "ASP":
        o1, o2 = "OD1", "OD2"
        h1, h2 = "HD1", "HD2"
    else:
        raise RuntimeError(
            f"Endpoint mask {mask!r} resolves to {resname} {resid}; "
            "this water-wire implementation currently requires Asp/Glu endpoints"
        )

    selected_names = {a.name for a in atoms}
    if selected_names != {o1, o2}:
        raise RuntimeError(
            f"Endpoint mask {mask!r} selected {sorted(selected_names)}, "
            f"expected exactly {o1},{o2}"
        )

    # Explicit residue/name masks. These are independent streams, so no output
    # atom-name inference is needed.
    return EndpointSpec(
        user_mask=mask,
        resid=resid,
        resname=resname,
        family=family,
        o1_name=o1,
        o2_name=o2,
        h1_name=h1,
        h2_name=h2,
        o1_mask=f":{resid}@{o1}",
        o2_mask=f":{resid}@{o2}",
        h1_mask=f":{resid}@{h1}",
        h2_mask=f":{resid}@{h2}",
    )


def validate_single_atom_mask(
    topology: Path,
    mask: str,
    cpptraj: str,
    context: str,
) -> None:
    atoms = cpptraj_mask_atoms(topology, mask, cpptraj)
    if len(atoms) != 1:
        raise RuntimeError(
            f"{context}: mask {mask!r} must select exactly one atom; "
            f"selected {len(atoms)}"
        )


# =============================================================================
# LAMBDA PARSING / EXACT MATCHING
# =============================================================================

def _residue_matches(value: str, resid: int) -> bool:
    try:
        return int(float(value)) == int(resid)
    except (ValueError, TypeError):
        return str(value) == str(resid)


def read_lambda_file(path: Path) -> LambdaFile:
    ires: List[str] = []
    itauto: List[str] = []
    expected: Optional[int] = None
    rows: Dict[int, Tuple[float, ...]] = {}

    with path.open() as fh:
        for lineno, raw in enumerate(fh, 1):
            s = raw.strip()
            if not s:
                continue

            parts = s.split()

            if parts[0] == "#":
                if len(parts) >= 3 and parts[1] == "ititr":
                    expected = len(parts) - 2
                elif len(parts) >= 3 and parts[1] == "ires":
                    ires = parts[2:]
                elif len(parts) >= 3 and parts[1] == "itauto":
                    itauto = parts[2:]
                continue

            try:
                step_f = float(parts[0])
                step = int(round(step_f))
            except ValueError as exc:
                raise RuntimeError(
                    f"{path}:{lineno}: invalid MD step"
                ) from exc

            if abs(step_f - step) > 1e-6:
                raise RuntimeError(
                    f"{path}:{lineno}: non-integral MD step {step_f}"
                )

            try:
                vals = tuple(float(x) for x in parts[1:])
            except ValueError as exc:
                raise RuntimeError(
                    f"{path}:{lineno}: invalid lambda value"
                ) from exc

            if expected is not None and len(vals) != expected:
                raise RuntimeError(
                    f"{path}:{lineno}: expected {expected} lambda columns, "
                    f"found {len(vals)}"
                )

            if step in rows:
                raise RuntimeError(
                    f"{path}:{lineno}: duplicate lambda step {step}"
                )

            rows[step] = vals

    if expected is None:
        raise RuntimeError(f"{path}: missing '# ititr' header")

    if not ires:
        raise RuntimeError(f"{path}: missing '# ires' header")

    # Same compact-header compatibility as calc_hbond_v5_1.py.
    if len(ires) != expected:
        compact = "".join(ires)
        if (
            not compact.isdigit()
            or len(compact) % expected != 0
        ):
            raise RuntimeError(
                f"{path}: cannot interpret ires header "
                f"({len(ires)} tokens for {expected} columns)"
            )

        width = len(compact) // expected
        ires = [
            compact[i:i + width]
            for i in range(0, len(compact), width)
        ]

    if len(ires) != expected:
        raise RuntimeError(f"{path}: ires column count mismatch")

    if itauto and len(itauto) != expected:
        raise RuntimeError(f"{path}: itauto column count mismatch")

    if not itauto:
        itauto = [""] * expected

    if not rows:
        raise RuntimeError(f"{path}: no lambda rows")

    steps = sorted(rows)

    return LambdaFile(
        path=path.resolve(),
        ires=tuple(ires),
        itauto=tuple(itauto),
        step_to_values=rows,
        min_step=steps[0],
        max_step=steps[-1],
    )


def lambda_columns_for_residue(
    lam: LambdaFile,
    resid: int,
) -> List[int]:
    return [
        i
        for i, r in enumerate(lam.ires)
        if _residue_matches(r, resid)
    ]


def itauto_int(
    lam: LambdaFile,
    col: int,
) -> Optional[int]:
    if col >= len(lam.itauto):
        return None

    try:
        return int(float(lam.itauto[col]))
    except (ValueError, TypeError):
        return None


def columns_by_flag(
    lam: LambdaFile,
    resid: int,
) -> Dict[int, List[int]]:
    out: Dict[int, List[int]] = {}

    for col in lambda_columns_for_residue(lam, resid):
        flag = itauto_int(lam, col)
        if flag is not None:
            out.setdefault(flag, []).append(col)

    return out


def threshold(
    x: float,
    low: float,
    high: float,
) -> str:
    if x <= low:
        return "low"
    if x >= high:
        return "high"
    return "mixed"


def lambda_row_exact(
    lam: LambdaFile,
    step: int,
) -> Tuple[str, Optional[Tuple[float, ...]]]:
    if step in lam.step_to_values:
        return "matched", lam.step_to_values[step]

    if step < lam.min_step or step > lam.max_step:
        return "outside_lambda_window", None

    raise RuntimeError(
        f"{lam.path.name}: trajectory step {step} lies inside lambda range "
        f"{lam.min_step}:{lam.max_step} but has no exact lambda row"
    )


def match_lambda_files(
    trajectories: List[Path],
    lambda_files: List[Path],
) -> Dict[Path, Path]:

    by_stem = {p.stem: p for p in lambda_files}
    mapping: Dict[Path, Path] = {}

    for traj in trajectories:
        if traj.stem in by_stem:
            mapping[traj] = by_stem[traj.stem]
            continue

        tph = extract_ph(traj.stem)

        ph_matches = [
            p
            for p in lambda_files
            if (
                math.isfinite(tph)
                and math.isfinite(extract_ph(p.stem))
                and abs(extract_ph(p.stem) - tph) < 1e-8
            )
        ]

        if len(ph_matches) == 1:
            mapping[traj] = ph_matches[0]
            continue

        if len(lambda_files) == 1 and len(trajectories) == 1:
            mapping[traj] = lambda_files[0]
            continue

        raise RuntimeError(
            f"Could not uniquely match trajectory {traj.name} "
            "to a lambda file"
        )

    return mapping


# =============================================================================
# CpH ENDPOINT CHEMISTRY
# =============================================================================

def endpoint_roles_for_frame(
    spec: EndpointSpec,
    lam: LambdaFile,
    values: Tuple[float, ...],
    low: float,
    high: float,
) -> EndpointRoles:
    """
    Asp/Glu chemistry copied semantically from calc_hbond_v5_1.py.

    Oxygen index convention in coordinate arrays:
        0 = O1 (OD1/OE1)
        1 = O2 (OD2/OE2)

    Hydrogen index convention:
        0 = H1 (HD1/HE1)
        1 = H2 (HD2/HE2)
    """

    flags = columns_by_flag(lam, spec.resid)

    pcols = flags.get(3, [])
    tcols = flags.get(4, [])

    if len(pcols) != 1 or len(tcols) != 1:
        raise RuntimeError(
            f"Residue {spec.resid} {spec.resname}: expected one itauto=3 "
            f"protonation and one itauto=4 tautomer coordinate; "
            f"found itauto3={pcols}, itauto4={tcols}"
        )

    pstate = threshold(
        values[pcols[0]],
        low,
        high,
    )

    if pstate == "mixed":
        return EndpointRoles(
            status="mixed",
            protonated=None,
            donor_o=(),
            acceptor_o=(),
            donor_h_for_o={},
        )

    if pstate == "high":
        # Deprotonated: both carboxylate oxygens are acceptors.
        return EndpointRoles(
            status="clean",
            protonated=False,
            donor_o=(),
            acceptor_o=(0, 1),
            donor_h_for_o={},
        )

    # Protonated; tautomer determines O-H identity.
    tstate = threshold(
        values[tcols[0]],
        low,
        high,
    )

    if tstate == "mixed":
        return EndpointRoles(
            status="mixed",
            protonated=True,
            donor_o=(),
            acceptor_o=(),
            donor_h_for_o={},
        )

    if tstate == "low":
        # Shen/Khandogin convention: tautomer low -> O2-H.
        return EndpointRoles(
            status="clean",
            protonated=True,
            donor_o=(1,),
            acceptor_o=(0,),
            donor_h_for_o={1: 1},
        )

    # tautomer high -> O1-H.
    return EndpointRoles(
        status="clean",
        protonated=True,
        donor_o=(0,),
        acceptor_o=(1,),
        donor_h_for_o={0: 0},
    )


# =============================================================================
# XYZ READER
# =============================================================================

def read_xyz_frames(
    path: Path,
) -> Iterable[np.ndarray]:

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
                    f"{path}: frame {frame_index}: expected XYZ atom count, "
                    f"got {line!r}"
                ) from exc

            comment = fh.readline()

            if not comment:
                raise RuntimeError(
                    f"{path}: frame {frame_index}: truncated after atom count"
                )

            xyz = np.empty(
                (natoms, 3),
                dtype=np.float64,
            )

            for i in range(natoms):
                row = fh.readline()

                if not row:
                    raise RuntimeError(
                        f"{path}: frame {frame_index}: truncated XYZ frame"
                    )

                fields = row.split()

                if len(fields) < 4:
                    raise RuntimeError(
                        f"{path}: frame {frame_index}: malformed XYZ row "
                        f"{row!r}"
                    )

                try:
                    xyz[i, 0] = float(fields[1])
                    xyz[i, 1] = float(fields[2])
                    xyz[i, 2] = float(fields[3])
                except ValueError as exc:
                    raise RuntimeError(
                        f"{path}: frame {frame_index}: malformed coordinates"
                    ) from exc

            yield xyz


# =============================================================================
# CPPTRAJ COORDINATE EXTRACTION
# =============================================================================

def endpoint_stream_paths(
    run_dir: Path,
    tag: str,
    label: str,
) -> Dict[str, Path]:

    return {
        "o1": run_dir / f"{tag}.{label}.O1.xyz",
        "o2": run_dir / f"{tag}.{label}.O2.xyz",
        "h1": run_dir / f"{tag}.{label}.H1.xyz",
        "h2": run_dir / f"{tag}.{label}.H2.xyz",
    }


def build_cpptraj_input(
    args,
    topology: Path,
    traj: Path,
    site_a: EndpointSpec,
    site_b: EndpointSpec,
    a_paths: Dict[str, Path],
    b_paths: Dict[str, Path],
    water_out: Path,
    closest_out: Path,
) -> str:

    if args.stop == "last":
        trajin = (
            f"trajin {qpath(traj)} "
            f"{args.start} last {args.stride}"
        )
    else:
        trajin = (
            f"trajin {qpath(traj)} "
            f"{args.start} {args.stop} {args.stride}"
        )

    combined_site = (
        f"({site_a.user_mask})|({site_b.user_mask})"
    )

    lines = [
        f"parm {qpath(topology)}",
        trajin,
        "",
        f"solvent {args.solvent_mask}",
        "",
    ]

    def add_single_atom_stream(mask: str, out: Path):
        if args.image_anchor:
            lines.append(
                f"autoimage anchor {qmask(args.image_anchor)}"
            )

        lines.extend([
            f"strip !({mask}) nobox",
            f"outtraj {qpath(out)} xyz",
            "unstrip",
            "",
        ])

    # Endpoint streams. Each stream contains exactly one atom per frame.
    for key, mask in (
        ("o1", site_a.o1_mask),
        ("o2", site_a.o2_mask),
        ("h1", site_a.h1_mask),
        ("h2", site_a.h2_mask),
    ):
        add_single_atom_stream(mask, a_paths[key])

    for key, mask in (
        ("o1", site_b.o1_mask),
        ("o2", site_b.o2_mask),
        ("h1", site_b.h1_mask),
        ("h2", site_b.h2_mask),
    ):
        add_single_atom_stream(mask, b_paths[key])

    # Local water selection. IMPORTANT: write O + the two physical H atoms
    # together from the SAME post-closest topology. Do not unstrip between
    # separate O/H outputs: CPPTRAJ unstrip can restore the pre-closest solvent
    # set, which would make the H stream contain all box waters.
    if args.image_anchor:
        lines.append(
            f"autoimage anchor {qmask(args.image_anchor)}"
        )

    water_triplet_mask = (
        f"({args.water_oxygen_mask})|({args.water_hydrogen_mask})"
    )

    lines.extend([
        (
            f"closest {args.closest_waters} "
            f"{qmask(combined_site)} oxygen "
            f"closestout {qpath(closest_out)}"
        ),
        "",
        f"strip !({water_triplet_mask}) nobox",
        f"outtraj {qpath(water_out)} xyz",
        "",
        "run",
        "quit",
        "",
    ])

    return "\n".join(lines)


# =============================================================================
# GEOMETRY
# =============================================================================

def pairwise_distances(
    a: np.ndarray,
    b: np.ndarray,
) -> np.ndarray:

    if len(a) == 0 or len(b) == 0:
        return np.empty(
            (len(a), len(b)),
            dtype=np.float64,
        )

    d = a[:, None, :] - b[None, :, :]

    return np.sqrt(
        np.einsum(
            "ijk,ijk->ij",
            d,
            d,
        )
    )


def minimum_distance(
    a: np.ndarray,
    b: np.ndarray,
) -> float:

    d = pairwise_distances(a, b)

    if d.size == 0:
        return math.nan

    return float(np.min(d))


def dha_angle_deg(
    donor: np.ndarray,
    hydrogen: np.ndarray,
    acceptor: np.ndarray,
) -> float:
    """
    D-H...A angle at H.
    """

    v1 = donor - hydrogen
    v2 = acceptor - hydrogen

    n1 = np.linalg.norm(v1)
    n2 = np.linalg.norm(v2)

    if n1 <= 0.0 or n2 <= 0.0:
        return math.nan

    c = float(
        np.dot(v1, v2) / (n1 * n2)
    )

    c = min(1.0, max(-1.0, c))

    return math.degrees(
        math.acos(c)
    )


def hbond_geometry(
    donor: np.ndarray,
    hydrogen: np.ndarray,
    acceptor: np.ndarray,
    distance_cutoff: float,
    angle_cutoff: float,
) -> bool:

    da = float(
        np.linalg.norm(
            donor - acceptor
        )
    )

    if da > distance_cutoff:
        return False

    if angle_cutoff < 0:
        return True

    angle = dha_angle_deg(
        donor,
        hydrogen,
        acceptor,
    )

    return (
        math.isfinite(angle)
        and angle >= angle_cutoff
    )


# =============================================================================
# WATER O/H PAIRING
# =============================================================================

def split_water_triplets(
    water_xyz: np.ndarray,
    nwater: int,
    max_oh: float,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Split one post-closest water coordinate stream into:

        water_o : (N, 3)
        water_h : (N, 2, 3)

    CPPTRAJ preserves retained topology atom order when stripping. With
    --water-oxygen-mask :WAT@O and --water-hydrogen-mask :WAT@H1,H2, each
    retained water contributes O,H1,H2 in topology order (virtual sites, if
    present, are excluded by the strip mask).

    We DO NOT trust the ordering silently: every putative O-H distance is
    validated against --water-oh-max. If the topology does not use O,H1,H2
    ordering, the analysis fails rather than silently constructing false wires.
    """

    expected = 3 * nwater

    if len(water_xyz) != expected:
        raise RuntimeError(
            f"Combined retained-water stream contains {len(water_xyz)} atoms; "
            f"expected exactly {expected} = 3 x {nwater}. "
            "Check --water-oxygen-mask and --water-hydrogen-mask."
        )

    grouped = water_xyz.reshape(
        nwater,
        3,
        3,
    )

    water_o = grouped[:, 0, :].copy()
    water_h = grouped[:, 1:3, :].copy()

    oh = np.linalg.norm(
        water_h - water_o[:, None, :],
        axis=2,
    )

    if not np.all(oh <= max_oh):
        bad = np.argwhere(
            oh > max_oh
        )
        examples = [
            (
                int(i),
                int(h),
                float(oh[i, h]),
            )
            for i, h in bad[:10]
        ]

        raise RuntimeError(
            "Retained-water O/H topology ordering validation failed: "
            f"some putative O-H distances exceed {max_oh:g} A. "
            f"Examples (water_slot,H_slot,distance_A): {examples}. "
            "Do not continue; adjust water atom masks/order for this topology."
        )

    return water_o, water_h


# =============================================================================
# DISTANCE-ONLY CANDIDATE GRAPH
# =============================================================================

def candidate_graph(
    site_a_o: np.ndarray,
    site_b_o: np.ndarray,
    water_o: np.ndarray,
    endpoint_cutoff: float,
    water_cutoff: float,
    max_waters: int,
) -> Dict[str, int]:

    nw = len(water_o)

    da = pairwise_distances(
        site_a_o,
        water_o,
    )

    db = pairwise_distances(
        site_b_o,
        water_o,
    )

    if nw:
        a_conn = (
            np.min(da, axis=0)
            <= endpoint_cutoff
        )
        b_conn = (
            np.min(db, axis=0)
            <= endpoint_cutoff
        )

        dww = pairwise_distances(
            water_o,
            water_o,
        )

        adj = (
            dww <= water_cutoff
        )

        np.fill_diagonal(
            adj,
            False,
        )

    else:
        a_conn = np.zeros(0, dtype=bool)
        b_conn = np.zeros(0, dtype=bool)
        adj = np.zeros((0, 0), dtype=bool)

    a_idx = np.flatnonzero(a_conn)

    n1 = int(
        np.count_nonzero(
            a_conn & b_conn
        )
    )

    n2 = 0
    if max_waters >= 2:
        for i in a_idx:
            n2 += int(
                np.count_nonzero(
                    adj[i] & b_conn
                )
            )

    n3 = 0
    if max_waters >= 3:
        for i in a_idx:
            for j in np.flatnonzero(adj[i]):
                for k in np.flatnonzero(
                    adj[j] & b_conn
                ):
                    if k == i or k == j:
                        continue
                    n3 += 1

    h1 = int(n1 > 0)
    h2 = int(n2 > 0)
    h3 = int(n3 > 0)

    if h1:
        shortest = 1
    elif h2:
        shortest = 2
    elif h3:
        shortest = 3
    else:
        shortest = -1

    return {
        "n_water_site_a": int(np.count_nonzero(a_conn)),
        "n_water_site_b": int(np.count_nonzero(b_conn)),
        "has1": h1,
        "has2": h2,
        "has3": h3,
        "connected": int(h1 or h2 or h3),
        "shortest": shortest,
    }


# =============================================================================
# H-BOND GRAPH
# =============================================================================

def endpoint_water_edges(
    endpoint_o: np.ndarray,
    endpoint_h: np.ndarray,
    roles: EndpointRoles,
    water_o: np.ndarray,
    water_h: np.ndarray,
    distance_cutoff: float,
    angle_cutoff: float,
) -> np.ndarray:
    """
    Return bool[Nwater]. True if at least one chemically allowed H bond exists
    between endpoint and that water in either direction.
    """

    nw = len(water_o)
    connected = np.zeros(nw, dtype=bool)

    # Endpoint donor -> water acceptor.
    for oi in roles.donor_o:
        hi = roles.donor_h_for_o[oi]

        donor = endpoint_o[oi]
        hydrogen = endpoint_h[hi]

        for w in range(nw):
            if hbond_geometry(
                donor,
                hydrogen,
                water_o[w],
                distance_cutoff,
                angle_cutoff,
            ):
                connected[w] = True

    # Water donor -> endpoint acceptor.
    for oi in roles.acceptor_o:
        acceptor = endpoint_o[oi]

        for w in range(nw):
            if connected[w]:
                # Still valid regardless of direction. No need to test more.
                continue

            donor = water_o[w]

            for h in water_h[w]:
                if hbond_geometry(
                    donor,
                    h,
                    acceptor,
                    distance_cutoff,
                    angle_cutoff,
                ):
                    connected[w] = True
                    break

    return connected


def water_water_adjacency(
    water_o: np.ndarray,
    water_h: np.ndarray,
    distance_cutoff: float,
    angle_cutoff: float,
) -> np.ndarray:
    """
    Undirected H-bond adjacency: i-j is True if either water can donate a valid
    H bond to the other in the frame.
    """

    n = len(water_o)

    adj = np.zeros(
        (n, n),
        dtype=bool,
    )

    if n < 2:
        return adj

    d = pairwise_distances(
        water_o,
        water_o,
    )

    # Cheap distance prefilter.
    ii, jj = np.where(
        np.triu(
            d <= distance_cutoff,
            k=1,
        )
    )

    for i, j in zip(ii.tolist(), jj.tolist()):
        valid = False

        # i donor -> j acceptor.
        for h in water_h[i]:
            if hbond_geometry(
                water_o[i],
                h,
                water_o[j],
                distance_cutoff,
                angle_cutoff,
            ):
                valid = True
                break

        # j donor -> i acceptor.
        if not valid:
            for h in water_h[j]:
                if hbond_geometry(
                    water_o[j],
                    h,
                    water_o[i],
                    distance_cutoff,
                    angle_cutoff,
                ):
                    valid = True
                    break

        if valid:
            adj[i, j] = True
            adj[j, i] = True

    return adj


def enumerate_hbond_paths(
    a_conn: np.ndarray,
    b_conn: np.ndarray,
    adj: np.ndarray,
    max_waters: int,
) -> Tuple[int, int, int, str]:

    a_idx = np.flatnonzero(a_conn)

    paths1: List[Tuple[int, ...]] = []
    paths2: List[Tuple[int, ...]] = []
    paths3: List[Tuple[int, ...]] = []

    # 1W: A-W-B
    for i in np.flatnonzero(a_conn & b_conn):
        paths1.append((int(i),))

    # 2W: A-Wi-Wj-B
    if max_waters >= 2:
        for i in a_idx:
            for j in np.flatnonzero(
                adj[i] & b_conn
            ):
                if j == i:
                    continue
                paths2.append(
                    (int(i), int(j))
                )

    # 3W: A-Wi-Wj-Wk-B
    if max_waters >= 3:
        for i in a_idx:
            for j in np.flatnonzero(adj[i]):
                if j == i:
                    continue

                for k in np.flatnonzero(
                    adj[j] & b_conn
                ):
                    if k == i or k == j:
                        continue

                    paths3.append(
                        (int(i), int(j), int(k))
                    )

    # Deduplicate exact oriented paths.
    paths1 = sorted(set(paths1))
    paths2 = sorted(set(paths2))
    paths3 = sorted(set(paths3))

    if paths1:
        representative = ",".join(map(str, paths1[0]))
    elif paths2:
        representative = ",".join(map(str, paths2[0]))
    elif paths3:
        representative = ",".join(map(str, paths3[0]))
    else:
        representative = ""

    return (
        len(paths1),
        len(paths2),
        len(paths3),
        representative,
    )


# =============================================================================
# SYNCHRONIZED STREAM READER
# =============================================================================

def next_or_done(it):
    try:
        return False, next(it)
    except StopIteration:
        return True, None


# =============================================================================
# PER-TRAJECTORY EXECUTION
# =============================================================================

def run_trajectory(
    args,
    topology: Path,
    trajectory: Path,
    lambda_path: Path,
    site_a: EndpointSpec,
    site_b: EndpointSpec,
) -> TrajectoryResult:

    tag = safe_component(
        trajectory.stem
    )

    run_dir = (
        args.workdir
        /
        tag
    )

    run_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    a_paths = endpoint_stream_paths(
        run_dir,
        tag,
        "siteA",
    )

    b_paths = endpoint_stream_paths(
        run_dir,
        tag,
        "siteB",
    )

    water_out = (
        run_dir
        /
        f"{tag}.water_OH.xyz"
    )

    closest_out = (
        run_dir
        /
        f"{tag}.closest.dat"
    )

    cpptraj_in = (
        run_dir
        /
        "water_wire.in"
    )

    cpptraj_log = (
        run_dir
        /
        "water_wire.log"
    )

    text = build_cpptraj_input(
        args,
        topology,
        trajectory,
        site_a,
        site_b,
        a_paths,
        b_paths,
        water_out,
        closest_out,
    )

    cpptraj_in.write_text(
        text
    )

    print(
        f"[{trajectory.stem}] cpptraj coordinate extraction...",
        flush=True,
    )

    proc = subprocess.run(
        [
            args.cpptraj,
            "-i",
            str(cpptraj_in),
        ],
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
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

    if proc.returncode != 0 or explicit_error:
        raise RuntimeError(
            f"CPPTRAJ failed for {trajectory}; inspect {cpptraj_log}"
        )

    expected_files = (
        list(a_paths.values())
        + list(b_paths.values())
        + [water_out]
    )

    for path in expected_files:
        if (
            not path.is_file()
            or path.stat().st_size == 0
        ):
            raise RuntimeError(
                f"Missing/empty coordinate stream {path}; "
                f"inspect {cpptraj_log}"
            )

    lam = read_lambda_file(
        lambda_path
    )

    # Ensure both endpoints are actually present in lambda metadata.
    for spec in (site_a, site_b):
        cols = lambda_columns_for_residue(
            lam,
            spec.resid,
        )

        if not cols:
            raise RuntimeError(
                f"{lambda_path}: endpoint residue {spec.resid} "
                f"({spec.resname}) has no lambda columns"
            )

    streams: Dict[str, Iterable[np.ndarray]] = {}

    for prefix, paths in (
        ("a", a_paths),
        ("b", b_paths),
    ):
        for key, path in paths.items():
            streams[f"{prefix}_{key}"] = iter(
                read_xyz_frames(path)
            )

    streams["water"] = iter(
        read_xyz_frames(
            water_out
        )
    )

    frames: List[FrameResult] = []

    selected_index = 0

    while True:
        read = {}

        done_flags = []

        for key, it in streams.items():
            done, xyz = next_or_done(it)
            done_flags.append(done)
            read[key] = xyz

        if all(done_flags):
            break

        if any(done_flags):
            ended = [
                key
                for key, done in zip(
                    streams.keys(),
                    done_flags,
                )
                if done
            ]

            raise RuntimeError(
                f"{trajectory}: synchronized coordinate streams ended "
                f"at different times; ended={ended}"
            )

        selected_index += 1

        # Every endpoint stream must contain exactly one atom.
        for key in (
            "a_o1", "a_o2", "a_h1", "a_h2",
            "b_o1", "b_o2", "b_h1", "b_h2",
        ):
            if len(read[key]) != 1:
                raise RuntimeError(
                    f"{trajectory}: frame {selected_index}: stream {key} "
                    f"contains {len(read[key])} atoms; expected 1"
                )

        water_o, water_h = split_water_triplets(
            read["water"],
            args.closest_waters,
            args.water_oh_max,
        )

        a_o = np.vstack([
            read["a_o1"][0],
            read["a_o2"][0],
        ])

        a_h = np.vstack([
            read["a_h1"][0],
            read["a_h2"][0],
        ])

        b_o = np.vstack([
            read["b_o1"][0],
            read["b_o2"][0],
        ])

        b_h = np.vstack([
            read["b_h1"][0],
            read["b_h2"][0],
        ])

        trajectory_frame = (
            args.start
            +
            (selected_index - 1)
            *
            args.stride
        )

        step = (
            args.first_frame_step
            +
            (trajectory_frame - 1)
            *
            args.frame_step_interval
        )

        min_ab = minimum_distance(
            a_o,
            b_o,
        )

        cgraph = candidate_graph(
            a_o,
            b_o,
            water_o,
            args.candidate_endpoint_cutoff,
            args.candidate_water_cutoff,
            args.max_waters,
        )

        lambda_status, values = lambda_row_exact(
            lam,
            step,
        )

        if lambda_status == "outside_lambda_window":
            frames.append(
                FrameResult(
                    trajectory_frame=trajectory_frame,
                    step=step,
                    chemistry_status="outside_lambda_window",

                    min_endpoint_distance_A=min_ab,
                    n_local_waters=len(water_o),
                    n_water_site_a_candidate=cgraph["n_water_site_a"],
                    n_water_site_b_candidate=cgraph["n_water_site_b"],

                    candidate_has_path_1w=cgraph["has1"],
                    candidate_has_path_2w=cgraph["has2"],
                    candidate_has_path_3w=cgraph["has3"],
                    candidate_connected_le3w=cgraph["connected"],
                    candidate_shortest_path_waters=cgraph["shortest"],

                    hbond_wire_1w=None,
                    hbond_wire_2w=None,
                    hbond_wire_3w=None,
                    hbond_connected_le3w=None,
                    hbond_shortest_path_waters=None,

                    n_hbond_paths_1w=None,
                    n_hbond_paths_2w=None,
                    n_hbond_paths_3w=None,

                    shortest_path_water_slots="",
                )
            )
            continue

        assert values is not None

        roles_a = endpoint_roles_for_frame(
            site_a,
            lam,
            values,
            args.lambda_low,
            args.lambda_high,
        )

        roles_b = endpoint_roles_for_frame(
            site_b,
            lam,
            values,
            args.lambda_low,
            args.lambda_high,
        )

        if (
            roles_a.status != "clean"
            or roles_b.status != "clean"
        ):
            frames.append(
                FrameResult(
                    trajectory_frame=trajectory_frame,
                    step=step,
                    chemistry_status="mixed_lambda",

                    min_endpoint_distance_A=min_ab,
                    n_local_waters=len(water_o),
                    n_water_site_a_candidate=cgraph["n_water_site_a"],
                    n_water_site_b_candidate=cgraph["n_water_site_b"],

                    candidate_has_path_1w=cgraph["has1"],
                    candidate_has_path_2w=cgraph["has2"],
                    candidate_has_path_3w=cgraph["has3"],
                    candidate_connected_le3w=cgraph["connected"],
                    candidate_shortest_path_waters=cgraph["shortest"],

                    hbond_wire_1w=None,
                    hbond_wire_2w=None,
                    hbond_wire_3w=None,
                    hbond_connected_le3w=None,
                    hbond_shortest_path_waters=None,

                    n_hbond_paths_1w=None,
                    n_hbond_paths_2w=None,
                    n_hbond_paths_3w=None,

                    shortest_path_water_slots="",
                )
            )
            continue

        a_conn = endpoint_water_edges(
            a_o,
            a_h,
            roles_a,
            water_o,
            water_h,
            args.distance,
            args.angle,
        )

        b_conn = endpoint_water_edges(
            b_o,
            b_h,
            roles_b,
            water_o,
            water_h,
            args.distance,
            args.angle,
        )

        ww_adj = water_water_adjacency(
            water_o,
            water_h,
            args.distance,
            args.angle,
        )

        n1, n2, n3, representative = enumerate_hbond_paths(
            a_conn,
            b_conn,
            ww_adj,
            args.max_waters,
        )

        h1 = int(n1 > 0)
        h2 = int(n2 > 0)
        h3 = int(n3 > 0)

        if h1:
            shortest = 1
        elif h2:
            shortest = 2
        elif h3:
            shortest = 3
        else:
            shortest = -1

        frames.append(
            FrameResult(
                trajectory_frame=trajectory_frame,
                step=step,
                chemistry_status="clean",

                min_endpoint_distance_A=min_ab,
                n_local_waters=len(water_o),
                n_water_site_a_candidate=cgraph["n_water_site_a"],
                n_water_site_b_candidate=cgraph["n_water_site_b"],

                candidate_has_path_1w=cgraph["has1"],
                candidate_has_path_2w=cgraph["has2"],
                candidate_has_path_3w=cgraph["has3"],
                candidate_connected_le3w=cgraph["connected"],
                candidate_shortest_path_waters=cgraph["shortest"],

                hbond_wire_1w=h1,
                hbond_wire_2w=h2,
                hbond_wire_3w=h3,
                hbond_connected_le3w=int(h1 or h2 or h3),
                hbond_shortest_path_waters=shortest,

                n_hbond_paths_1w=n1,
                n_hbond_paths_2w=n2,
                n_hbond_paths_3w=n3,

                shortest_path_water_slots=representative,
            )
        )

    if not frames:
        raise RuntimeError(
            f"{trajectory}: zero frames analyzed"
        )

    print(
        f"[{trajectory.stem}] complete: {len(frames)} frames",
        flush=True,
    )

    return TrajectoryResult(
        trajectory=trajectory,
        ph=extract_ph(trajectory.stem),
        lambda_file=lambda_path,
        frames=frames,
    )


# =============================================================================
# OUTPUT
# =============================================================================

OUTPUT_FIELDS = [
    "trajectory",
    "pH",
    "analysis_row",
    "trajectory_frame",
    "step",
    "observable",
    "metric",
    "interaction_class",
    "mode",
    "selection1",
    "selection2",
    "value",
    "unit",
    "chemistry_status",
    "lambda_file",
]


def iter_long_rows(
    args,
    result: TrajectoryResult,
):
    """
    Emit partition_by_state.py-ready long-format rows.

    Strict lambda-aware H-bond-wire metrics are emitted only for chemically
    clean frames. Distance-only candidate metrics are emitted for all frames
    because they do not require an endpoint donor/acceptor assignment.
    """

    for fr in result.frames:
        common = {
            "trajectory": result.trajectory.stem,
            "pH": "" if not math.isfinite(result.ph) else f"{result.ph:g}",
            "analysis_row": fr.trajectory_frame,
            "trajectory_frame": fr.trajectory_frame,
            "step": fr.step,
            "observable": "E35_D52_water_wire",
            "interaction_class": "solute_solvent_network",
            "mode": "water_wire",
            "selection1": "E35",
            "selection2": "D52",
            "chemistry_status": fr.chemistry_status,
            "lambda_file": str(result.lambda_file),
        }

        diagnostics = [
            ("min_endpoint_distance_A", fr.min_endpoint_distance_A, "angstrom"),
            ("candidate_n_water_site_a", fr.n_water_site_a_candidate, "count"),
            ("candidate_n_water_site_b", fr.n_water_site_b_candidate, "count"),
            ("candidate_wire_1W", fr.candidate_has_path_1w, "fraction"),
            ("candidate_wire_2W", fr.candidate_has_path_2w, "fraction"),
            ("candidate_wire_3W", fr.candidate_has_path_3w, "fraction"),
            ("candidate_wire_le3W", fr.candidate_connected_le3w, "fraction"),
            (
                "candidate_shortest_1W",
                int(fr.candidate_shortest_path_waters == 1),
                "fraction",
            ),
            (
                "candidate_shortest_2W",
                int(fr.candidate_shortest_path_waters == 2),
                "fraction",
            ),
            (
                "candidate_shortest_3W",
                int(fr.candidate_shortest_path_waters == 3),
                "fraction",
            ),
            (
                "candidate_shortest_none",
                int(fr.candidate_shortest_path_waters == -1),
                "fraction",
            ),
        ]

        for metric, value, unit in diagnostics:
            row = dict(common)
            row["metric"] = metric
            row["value"] = value
            row["unit"] = unit
            yield row

        # Do not turn mixed/outside chemistry into a numerical "no wire".
        if fr.chemistry_status != "clean":
            continue

        strict = [
            ("wire_1W", fr.hbond_wire_1w, "fraction"),
            ("wire_2W", fr.hbond_wire_2w, "fraction"),
            ("wire_3W", fr.hbond_wire_3w, "fraction"),
            ("wire_le3W", fr.hbond_connected_le3w, "fraction"),
            (
                "shortest_1W",
                int(fr.hbond_shortest_path_waters == 1),
                "fraction",
            ),
            (
                "shortest_2W",
                int(fr.hbond_shortest_path_waters == 2),
                "fraction",
            ),
            (
                "shortest_3W",
                int(fr.hbond_shortest_path_waters == 3),
                "fraction",
            ),
            (
                "shortest_none",
                int(fr.hbond_shortest_path_waters == -1),
                "fraction",
            ),
            ("n_paths_1W", fr.n_hbond_paths_1w, "count"),
            ("n_paths_2W", fr.n_hbond_paths_2w, "count"),
            ("n_paths_3W", fr.n_hbond_paths_3w, "count"),
        ]

        for metric, value, unit in strict:
            row = dict(common)
            row["metric"] = metric
            row["value"] = value
            row["unit"] = unit
            yield row


def write_output(
    args,
    results: Sequence[TrajectoryResult],
) -> int:

    nrows = 0

    with args.output.open("w", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=OUTPUT_FIELDS,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()

        for result in sorted(
            results,
            key=lambda r: (r.ph, r.trajectory.name),
        ):
            for row in iter_long_rows(args, result):
                out = dict(row)
                for key, value in list(out.items()):
                    if isinstance(value, float):
                        out[key] = (
                            ""
                            if not math.isfinite(value)
                            else f"{value:.10g}"
                        )
                    elif value is None:
                        out[key] = ""
                writer.writerow(out)
                nrows += 1

    return nrows


def write_settings(
    args,
    topology: Path,
    site_a: EndpointSpec,
    site_b: EndpointSpec,
    lambda_map: Dict[Path, Path],
    nrows: int,
) -> None:

    path = Path(str(args.output) + ".settings.txt")

    lines = [
        "program=calc_water_wire.py",
        f"version={PROGRAM_VERSION}",
        f"topology={topology}",
        f"site_a={site_a}",
        f"site_b={site_b}",
        f"distance_A={args.distance:g}",
        f"angle_deg={args.angle:g}",
        f"lambda_low={args.lambda_low:g}",
        f"lambda_high={args.lambda_high:g}",
        f"candidate_endpoint_cutoff_A={args.candidate_endpoint_cutoff:g}",
        f"candidate_water_cutoff_A={args.candidate_water_cutoff:g}",
        f"closest_waters={args.closest_waters}",
        f"max_waters={args.max_waters}",
        f"water_oxygen_mask={args.water_oxygen_mask}",
        f"water_hydrogen_mask={args.water_hydrogen_mask}",
        f"water_oh_max_A={args.water_oh_max:g}",
        f"image_anchor={args.image_anchor or ''}",
        f"start={args.start}",
        f"stop={args.stop}",
        f"stride={args.stride}",
        f"first_frame_step={args.first_frame_step}",
        f"frame_step_interval={args.frame_step_interval}",
        f"jobs={args.jobs}",
        f"cpptraj_threads={args.cpptraj_threads}",
        f"output_rows={nrows}",
        "output_format=partition_ready_long_tsv",
        "observable=E35_D52_water_wire",
        "strict_metric_policy=clean_chemistry_only",
        "candidate_metric_policy=all_frames",
        (
            "strict_fraction_metrics="
            "wire_1W,wire_2W,wire_3W,wire_le3W,"
            "shortest_1W,shortest_2W,shortest_3W,shortest_none"
        ),
        "strict_count_metrics=n_paths_1W,n_paths_2W,n_paths_3W",
    ]

    for traj, lam in sorted(
        lambda_map.items(),
        key=lambda x: natural_key(x[0]),
    ):
        lines.append(f"lambda_map[{traj.name}]={lam}")

    path.write_text("\n".join(lines) + "\n")


# =============================================================================
# TERMINAL SUMMARY
# =============================================================================

def print_summary(
    results: Sequence[TrajectoryResult],
) -> None:

    print()
    print("=" * 96)
    print("H-BONDED WATER-WIRE SUMMARY")
    print("=" * 96)

    for result in sorted(
        results,
        key=lambda r: (
            r.ph,
            r.trajectory.name,
        ),
    ):

        clean = [
            f
            for f in result.frames
            if f.chemistry_status == "clean"
        ]

        nmixed = sum(
            f.chemistry_status == "mixed_lambda"
            for f in result.frames
        )

        noutside = sum(
            f.chemistry_status == "outside_lambda_window"
            for f in result.frames
        )

        print()
        print(
            f"{result.trajectory.name}  "
            f"pH={result.ph:g}  "
            f"frames={len(result.frames)}  "
            f"clean={len(clean)}  mixed={nmixed}  outside={noutside}"
        )

        if clean:
            shortest = np.array(
                [
                    f.hbond_shortest_path_waters
                    for f in clean
                ],
                dtype=int,
            )

            h1 = np.mean(
                [
                    f.hbond_wire_1w
                    for f in clean
                ]
            )

            h2 = np.mean(
                [
                    f.hbond_wire_2w
                    for f in clean
                ]
            )

            h3 = np.mean(
                [
                    f.hbond_wire_3w
                    for f in clean
                ]
            )

            hc = np.mean(
                [
                    f.hbond_connected_le3w
                    for f in clean
                ]
            )

            print(
                "  H-bond path exists: "
                f"1W={h1:.4f}  "
                f"2W={h2:.4f}  "
                f"3W={h3:.4f}  "
                f"<=3W={hc:.4f}"
            )

            print(
                "  shortest H-bond:    "
                f"1W={np.mean(shortest == 1):.4f}  "
                f"2W={np.mean(shortest == 2):.4f}  "
                f"3W={np.mean(shortest == 3):.4f}  "
                f"none={np.mean(shortest == -1):.4f}"
            )

        d = np.array(
            [
                f.min_endpoint_distance_A
                for f in result.frames
            ],
            dtype=float,
        )

        print(
            "  endpoint O-O:         "
            f"median={np.nanmedian(d):.3f} A  "
            f"q10={np.nanquantile(d, 0.10):.3f}  "
            f"q90={np.nanquantile(d, 0.90):.3f}"
        )


# =============================================================================
# CLI
# =============================================================================

def parse_args() -> argparse.Namespace:

    ap = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=(
            "Lambda-aware H-bonded water-wire analysis using CPPTRAJ "
            "coordinate extraction and NumPy graph search."
        ),
    )

    ap.add_argument(
        "-p",
        "--topology",
        type=Path,
        required=True,
    )

    ap.add_argument(
        "-t",
        "--trajectory",
        action="append",
        required=True,
        help="Trajectory path or quoted glob; may be repeated",
    )

    ap.add_argument(
        "--lambda",
        dest="lambda_specs",
        action="append",
        required=True,
        help="Matching CpHMD lambda path or quoted glob; may be repeated",
    )

    ap.add_argument(
        "--site-a",
        required=True,
        help="Exactly two Asp/Glu carboxylate O atoms, e.g. ':36@OE1,OE2'",
    )

    ap.add_argument(
        "--site-b",
        required=True,
        help="Exactly two Asp/Glu carboxylate O atoms, e.g. ':53@OD1,OD2'",
    )

    ap.add_argument(
        "--lambda-low",
        type=float,
        default=0.2,
    )

    ap.add_argument(
        "--lambda-high",
        type=float,
        default=0.8,
    )

    ap.add_argument(
        "--distance",
        type=float,
        default=3.0,
        help="H-bond donor-heavy/acceptor cutoff in A",
    )

    ap.add_argument(
        "--angle",
        type=float,
        default=135.0,
        help="D-H...A angle cutoff in degrees; -1 disables angle filtering",
    )

    ap.add_argument(
        "--candidate-endpoint-cutoff",
        type=float,
        default=3.4,
        help="Distance-only diagnostic endpoint-water cutoff in A",
    )

    ap.add_argument(
        "--candidate-water-cutoff",
        type=float,
        default=3.5,
        help="Distance-only diagnostic water-water cutoff in A",
    )

    ap.add_argument(
        "--max-waters",
        type=int,
        default=3,
        choices=[1, 2, 3],
    )

    ap.add_argument(
        "--closest-waters",
        type=int,
        default=32,
    )

    ap.add_argument(
        "--solvent-mask",
        default=":WAT",
    )

    ap.add_argument(
        "--water-oxygen-mask",
        default=":WAT@O",
    )

    ap.add_argument(
        "--water-hydrogen-mask",
        default=":WAT@H1,H2",
        help="Must select exactly two physical H atoms per retained water",
    )

    ap.add_argument(
        "--water-oh-max",
        type=float,
        default=1.25,
        help="Maximum O-H distance used to validate water O/H grouping",
    )

    ap.add_argument(
        "--image-anchor",
        default=None,
    )

    ap.add_argument(
        "--start",
        type=int,
        default=1,
    )

    ap.add_argument(
        "--stop",
        default="last",
    )

    ap.add_argument(
        "--stride",
        type=int,
        default=1,
    )

    ap.add_argument(
        "--first-frame-step",
        type=int,
        required=True,
    )

    ap.add_argument(
        "--frame-step-interval",
        type=int,
        required=True,
    )

    ap.add_argument(
        "-j",
        "--jobs",
        type=int,
        default=1,
    )

    ap.add_argument(
        "--cpptraj-threads",
        type=int,
        default=1,
    )

    ap.add_argument(
        "--cpptraj",
        default="cpptraj",
    )

    ap.add_argument(
        "--workdir",
        type=Path,
        default=Path("water_wire_hbond_work"),
    )

    ap.add_argument(
        "-o",
        "--output",
        type=Path,
        default=Path("water_wire_hbond.tsv"),
    )

    ap.add_argument(
        "--keep-intermediates",
        action="store_true",
    )

    return ap.parse_args()


def parse_stop(value: str) -> str:
    if value == "last":
        return value

    try:
        n = int(value)
    except ValueError as exc:
        raise ValueError(
            "--stop must be a positive integer or 'last'"
        ) from exc

    if n < 1:
        raise ValueError(
            "--stop must be >= 1"
        )

    return str(n)


def validate_args(args: argparse.Namespace) -> None:

    if not args.topology.is_file():
        raise FileNotFoundError(
            args.topology
        )

    if args.start < 1:
        raise ValueError(
            "--start must be >= 1"
        )

    args.stop = parse_stop(
        args.stop
    )

    if (
        args.stop != "last"
        and int(args.stop) < args.start
    ):
        raise ValueError(
            "--stop must be >= --start"
        )

    if args.stride < 1:
        raise ValueError(
            "--stride must be >= 1"
        )

    if args.jobs < 1:
        raise ValueError(
            "--jobs must be >= 1"
        )

    if args.cpptraj_threads < 1:
        raise ValueError(
            "--cpptraj-threads must be >= 1"
        )

    if args.closest_waters < 3:
        raise ValueError(
            "--closest-waters must be >= 3"
        )

    if args.frame_step_interval <= 0:
        raise ValueError(
            "--frame-step-interval must be > 0"
        )

    if not (
        0 <= args.lambda_low
        < args.lambda_high
        <= 1
    ):
        raise ValueError(
            "lambda thresholds must satisfy 0 <= low < high <= 1"
        )

    if not (
        math.isfinite(args.distance)
        and args.distance > 0
    ):
        raise ValueError(
            "--distance must be > 0"
        )

    if not (
        math.isfinite(args.angle)
        and (
            args.angle == -1
            or 0 <= args.angle <= 180
        )
    ):
        raise ValueError(
            "--angle must be -1 or within 0..180"
        )

    if args.water_oh_max <= 0:
        raise ValueError(
            "--water-oh-max must be > 0"
        )

    exe = (
        shutil.which(args.cpptraj)
        if os.sep not in args.cpptraj
        else (
            args.cpptraj
            if Path(args.cpptraj).is_file()
            else None
        )
    )

    if exe is None:
        raise FileNotFoundError(
            f"CPPTRAJ executable not found: {args.cpptraj!r}"
        )

    args.cpptraj = str(exe)


# =============================================================================
# MAIN
# =============================================================================

def main() -> int:

    args = parse_args()

    validate_args(
        args
    )

    topology = (
        args.topology
        .resolve()
    )

    trajectories: List[Path] = []

    for spec in args.trajectory:
        trajectories.extend(
            expand_paths(
                [spec],
                "trajectory",
            )
        )

    # Deduplicate after repeated globs.
    trajectories = sorted(
        {
            p.resolve()
            for p in trajectories
        },
        key=natural_key,
    )

    lambda_files: List[Path] = []

    for spec in args.lambda_specs:
        lambda_files.extend(
            expand_paths(
                [spec],
                "lambda",
            )
        )

    lambda_files = sorted(
        {
            p.resolve()
            for p in lambda_files
        },
        key=natural_key,
    )

    lambda_map = match_lambda_files(
        trajectories,
        lambda_files,
    )

    site_a = build_endpoint_spec(
        topology,
        args.site_a,
        args.cpptraj,
    )

    site_b = build_endpoint_spec(
        topology,
        args.site_b,
        args.cpptraj,
    )

    if site_a.resid == site_b.resid:
        raise RuntimeError(
            "site A and site B resolve to the same residue"
        )

    # Validate generated endpoint atom masks. Constant-pH topologies used by the
    # existing calc_hbond.py are expected to carry both tautomeric donor H atoms.
    for label, spec in (
        ("site A", site_a),
        ("site B", site_b),
    ):
        for atom_label, mask in (
            ("O1", spec.o1_mask),
            ("O2", spec.o2_mask),
            ("H1", spec.h1_mask),
            ("H2", spec.h2_mask),
        ):
            validate_single_atom_mask(
                topology,
                mask,
                args.cpptraj,
                f"{label} {atom_label}",
            )

    # Validate water masks on original topology.
    water_o_atoms = cpptraj_mask_atoms(
        topology,
        args.water_oxygen_mask,
        args.cpptraj,
    )

    water_h_atoms = cpptraj_mask_atoms(
        topology,
        args.water_hydrogen_mask,
        args.cpptraj,
    )

    if not water_o_atoms:
        raise RuntimeError(
            f"--water-oxygen-mask {args.water_oxygen_mask!r} selected zero atoms"
        )

    if not water_h_atoms:
        raise RuntimeError(
            f"--water-hydrogen-mask {args.water_hydrogen_mask!r} selected zero atoms"
        )

    # Global topology check: for ordinary water masks we expect exactly two
    # physical hydrogens per O. This also catches wrong OPC hydrogen naming.
    ratio = len(water_h_atoms) / len(water_o_atoms)

    if abs(ratio - 2.0) > 1e-12:
        raise RuntimeError(
            f"Water masks select {len(water_o_atoms)} O and "
            f"{len(water_h_atoms)} H atoms (H/O={ratio:g}), not exactly 2. "
            "Adjust --water-hydrogen-mask to select the two physical water H atoms."
        )

    print(
        f"Site A : {site_a.resname} {site_a.resid} "
        f"{site_a.o1_name}/{site_a.o2_name}"
    )

    print(
        f"Site B : {site_b.resname} {site_b.resid} "
        f"{site_b.o1_name}/{site_b.o2_name}"
    )

    print(
        f"Trajectories : {len(trajectories)}"
    )

    print(
        f"H-bond geometry : D-A <= {args.distance:g} A, "
        f"D-H...A >= {args.angle:g} deg"
    )

    print(
        f"Closest waters : {args.closest_waters}"
    )

    args.workdir.mkdir(
        parents=True,
        exist_ok=True,
    )

    args.output.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    results: List[TrajectoryResult] = []

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
                lambda_map[traj],
                site_a,
                site_b,
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

    nrows = write_output(
        args,
        results,
    )

    write_settings(
        args,
        topology,
        site_a,
        site_b,
        lambda_map,
        nrows,
    )

    print_summary(
        results
    )

    if not args.keep_intermediates:
        for result in results:
            tag = safe_component(
                result.trajectory.stem
            )
            run_dir = (
                args.workdir
                /
                tag
            )

            for path in run_dir.glob("*.xyz"):
                try:
                    path.unlink()
                except FileNotFoundError:
                    pass

    print()
    print("=" * 96)
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
    print("=" * 96)

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
