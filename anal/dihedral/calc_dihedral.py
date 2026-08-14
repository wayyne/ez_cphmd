#!/usr/bin/env python3
"""
calc_dihedral.py
================
General protein side-chain dihedral / rotamer observable generator for Amber
trajectories. Designed to integrate with the structural-analysis suite and
partition_by_state.py.

Key design choices
------------------
* State-agnostic: no CpH lambda/state information is needed to measure geometry.
* One cpptraj pass per trajectory for all requested torsions.
* Standard protein chi definitions are registry-driven.
* Arbitrary four-mask dihedrals are supported.
* Angles are wrapped to [-180, 180).
* Long-format output is directly usable by partition_by_state.py.
* Rotamer populations are represented by 0/1 indicator metrics, so means from
  partition_by_state.py are true state-conditioned populations.
* sin/cos metrics are emitted so circular means can be reconstructed safely as
      atan2(mean_sin, mean_cos)
  instead of taking an arithmetic mean of angles.
* Symmetry-corrected terminal torsions are emitted where appropriate.

For histidine, the definitions intentionally reproduce the prior HEWL H15 code:
    chi1 = N-CA-CB-CG
    chi2 = CA-CB-CG-ND1

The program uses cpptraj for the actual dihedral calculation.
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


PROGRAM_VERSION = "1.0.2"


# Registry values are atom names in order. symmetry_period=180 means the
# terminal group has a twofold atom-label symmetry for this torsion and a
# symmetry-canonical angle is useful for cross-force-field comparison.
#
# Standard definitions follow conventional protein chi nomenclature and are
# chosen to match cpptraj/Amber atom names used by the user's HEWL systems.
CHI_REGISTRY: Dict[str, Dict[str, Tuple[Tuple[str, str, str, str], Optional[int]]]] = {
    "ARG": {
        "chi1": (("N", "CA", "CB", "CG"), None),
        "chi2": (("CA", "CB", "CG", "CD"), None),
        "chi3": (("CB", "CG", "CD", "NE"), None),
        "chi4": (("CG", "CD", "NE", "CZ"), None),
    },
    "ASN": {
        "chi1": (("N", "CA", "CB", "CG"), None),
        "chi2": (("CA", "CB", "CG", "OD1"), None),
    },
    "ASP": {
        "chi1": (("N", "CA", "CB", "CG"), None),
        "chi2": (("CA", "CB", "CG", "OD1"), 180),
    },
    "ASH": {
        "chi1": (("N", "CA", "CB", "CG"), None),
        "chi2": (("CA", "CB", "CG", "OD1"), None),
    },
    "AS2": {
        "chi1": (("N", "CA", "CB", "CG"), None),
        "chi2": (("CA", "CB", "CG", "OD1"), None),
    },
    "AS4": {
        "chi1": (("N", "CA", "CB", "CG"), None),
        "chi2": (("CA", "CB", "CG", "OD1"), None),
    },
    "CYS": {"chi1": (("N", "CA", "CB", "SG"), None)},
    "CYM": {"chi1": (("N", "CA", "CB", "SG"), None)},
    "CYX": {"chi1": (("N", "CA", "CB", "SG"), None)},
    "GLN": {
        "chi1": (("N", "CA", "CB", "CG"), None),
        "chi2": (("CA", "CB", "CG", "CD"), None),
        "chi3": (("CB", "CG", "CD", "OE1"), None),
    },
    "GLU": {
        "chi1": (("N", "CA", "CB", "CG"), None),
        "chi2": (("CA", "CB", "CG", "CD"), None),
        "chi3": (("CB", "CG", "CD", "OE1"), 180),
    },
    "GLH": {
        "chi1": (("N", "CA", "CB", "CG"), None),
        "chi2": (("CA", "CB", "CG", "CD"), None),
        "chi3": (("CB", "CG", "CD", "OE1"), None),
    },
    "GL2": {
        "chi1": (("N", "CA", "CB", "CG"), None),
        "chi2": (("CA", "CB", "CG", "CD"), None),
        "chi3": (("CB", "CG", "CD", "OE1"), None),
    },
    "GL4": {
        "chi1": (("N", "CA", "CB", "CG"), None),
        "chi2": (("CA", "CB", "CG", "CD"), None),
        "chi3": (("CB", "CG", "CD", "OE1"), None),
    },
    "HIS": {
        "chi1": (("N", "CA", "CB", "CG"), None),
        "chi2": (("CA", "CB", "CG", "ND1"), None),
    },
    "HID": {
        "chi1": (("N", "CA", "CB", "CG"), None),
        "chi2": (("CA", "CB", "CG", "ND1"), None),
    },
    "HIE": {
        "chi1": (("N", "CA", "CB", "CG"), None),
        "chi2": (("CA", "CB", "CG", "ND1"), None),
    },
    "HIP": {
        "chi1": (("N", "CA", "CB", "CG"), None),
        "chi2": (("CA", "CB", "CG", "ND1"), None),
    },
    "HSD": {
        "chi1": (("N", "CA", "CB", "CG"), None),
        "chi2": (("CA", "CB", "CG", "ND1"), None),
    },
    "HSE": {
        "chi1": (("N", "CA", "CB", "CG"), None),
        "chi2": (("CA", "CB", "CG", "ND1"), None),
    },
    "HSP": {
        "chi1": (("N", "CA", "CB", "CG"), None),
        "chi2": (("CA", "CB", "CG", "ND1"), None),
    },
    "ILE": {
        "chi1": (("N", "CA", "CB", "CG1"), None),
        "chi2": (("CA", "CB", "CG1", "CD1"), None),
    },
    "LEU": {
        "chi1": (("N", "CA", "CB", "CG"), None),
        "chi2": (("CA", "CB", "CG", "CD1"), None),
    },
    "LYS": {
        "chi1": (("N", "CA", "CB", "CG"), None),
        "chi2": (("CA", "CB", "CG", "CD"), None),
        "chi3": (("CB", "CG", "CD", "CE"), None),
        "chi4": (("CG", "CD", "CE", "NZ"), None),
    },
    "LYN": {
        "chi1": (("N", "CA", "CB", "CG"), None),
        "chi2": (("CA", "CB", "CG", "CD"), None),
        "chi3": (("CB", "CG", "CD", "CE"), None),
        "chi4": (("CG", "CD", "CE", "NZ"), None),
    },
    "MET": {
        "chi1": (("N", "CA", "CB", "CG"), None),
        "chi2": (("CA", "CB", "CG", "SD"), None),
        "chi3": (("CB", "CG", "SD", "CE"), None),
    },
    "PHE": {
        "chi1": (("N", "CA", "CB", "CG"), None),
        "chi2": (("CA", "CB", "CG", "CD1"), 180),
    },
    "PRO": {
        "chi1": (("N", "CA", "CB", "CG"), None),
        "chi2": (("CA", "CB", "CG", "CD"), None),
    },
    "SER": {"chi1": (("N", "CA", "CB", "OG"), None)},
    "THR": {"chi1": (("N", "CA", "CB", "OG1"), None)},
    "TRP": {
        "chi1": (("N", "CA", "CB", "CG"), None),
        "chi2": (("CA", "CB", "CG", "CD1"), None),
    },
    "TYR": {
        "chi1": (("N", "CA", "CB", "CG"), None),
        "chi2": (("CA", "CB", "CG", "CD1"), 180),
    },
    "TYM": {
        "chi1": (("N", "CA", "CB", "CG"), None),
        "chi2": (("CA", "CB", "CG", "CD1"), 180),
    },
    "VAL": {"chi1": (("N", "CA", "CB", "CG1"), None)},
}

# N/C-terminal Amber residue-name prefixes are normalized before registry lookup.
_BASE_NAMES = set(CHI_REGISTRY)


@dataclass(frozen=True)
class ResidueInfo:
    number: int
    name: str
    first_atom: int = -1
    last_atom: int = -1


@dataclass(frozen=True)
class TorsionSpec:
    observable: str
    residue_label: str
    residue_mask: str
    residue_number: Optional[int]
    residue_name: str
    torsion: str
    atom_names: Tuple[str, str, str, str]
    atom_masks: Tuple[str, str, str, str]
    symmetry_period_deg: Optional[int]
    classify_rotamer: bool


@dataclass
class TrajectoryResult:
    trajectory: Path
    ph: float
    nframes: int
    data: Dict[str, List[float]]


def die(msg: str) -> "None":
    raise RuntimeError(msg)


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
    # stable unique order
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


def normalize_resname(name: str) -> str:
    n = name.strip().upper()
    if n in _BASE_NAMES:
        return n
    if len(n) >= 4 and n[0] in {"N", "C"} and n[1:] in _BASE_NAMES:
        return n[1:]
    return n


def parse_ph(path: Path, regex: re.Pattern) -> float:
    m = regex.search(path.name) or regex.search(str(path))
    if not m:
        raise RuntimeError(
            f"Could not parse pH from trajectory {path}; regex={regex.pattern!r}"
        )
    return float(m.group("ph").replace("_", "."))


def selected_atoms(cpptraj: str, topology: Path, mask: str) -> List[Tuple[int, str, int, str]]:
    """Return selected atom records parsed from cpptraj `-ms` output when available.

    The parser is intentionally permissive across cpptraj versions. If detailed
    records are not printed, a fallback returns only atom indices with blank names.
    """
    text = run_command([cpptraj, "-p", str(topology), "-ms", mask])
    recs: List[Tuple[int, str, int, str]] = []
    # Common cpptraj mask table lines include atom#, atomname, res#, resname.
    for line in text.splitlines():
        m = re.match(r"\s*(\d+)\s+([A-Za-z0-9'+*_-]+)\s+(\d+)\s+([A-Za-z0-9'+*_-]+)\b", line)
        if m:
            recs.append((int(m.group(1)), m.group(2), int(m.group(3)), m.group(4)))
    if recs:
        return recs
    # Fallback used by some builds: "Selected = 1 2 3 ..."
    ids: List[int] = []
    for line in text.splitlines():
        m = re.search(r"Selected\s*=\s*(.*)$", line)
        if m:
            ids.extend(int(x) for x in re.findall(r"\d+", m.group(1)))
    return [(i, "", -1, "") for i in ids]


def inspect_topology_residues(cpptraj: str, topology: Path, workdir: Path) -> List[ResidueInfo]:
    """Read cpptraj resinfo for the topology once, including atom ranges."""
    workdir.mkdir(parents=True, exist_ok=True)
    out = workdir / "topology_resinfo.dat"
    log = workdir / "topology_resinfo.log"
    script = f"resinfo :* out {qpath(out)}\nquit\n"
    text = run_command([cpptraj, "-p", str(topology)], input_text=script)
    log.write_text(text)
    if not out.exists() or out.stat().st_size == 0:
        raise RuntimeError(f"cpptraj did not create topology resinfo file {out}; inspect {log}")
    records: List[ResidueInfo] = []
    with out.open("r", errors="replace") as fh:
        for line in fh:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            toks = s.split()
            if len(toks) < 4:
                continue
            try:
                number = int(toks[0])
                first_atom = int(toks[2])
                last_atom = int(toks[3])
            except ValueError:
                continue
            records.append(ResidueInfo(number, toks[1], first_atom, last_atom))
    if not records:
        raise RuntimeError(f"Could not parse cpptraj resinfo output {out}")
    return records


def topology_residue_info(
    cpptraj: str, topology: Path, residue_mask: str, residues: Sequence[ResidueInfo]
) -> ResidueInfo:
    """Resolve an arbitrary Amber mask to exactly one topology residue."""
    recs = selected_atoms(cpptraj, topology, residue_mask)
    atom_ids = sorted({r[0] for r in recs})
    if not atom_ids:
        raise RuntimeError(f"Residue mask {residue_mask!r} selected zero atoms")
    hits = [
        r for r in residues
        if any(r.first_atom <= atom <= r.last_atom for atom in atom_ids)
    ]
    uniq = {r.number: r for r in hits}
    if len(uniq) != 1:
        raise RuntimeError(
            f"Residue mask {residue_mask!r} must select exactly one residue; "
            f"resolved residues={[(r.number, r.name) for r in uniq.values()]}"
        )
    return next(iter(uniq.values()))


def atom_mask(residue_mask: str, atom_name: str) -> str:
    return f"({residue_mask})&@{atom_name}"


def validate_single_atom(cpptraj: str, topology: Path, mask: str, context: str) -> None:
    recs = selected_atoms(cpptraj, topology, mask)
    if len(recs) != 1:
        raise RuntimeError(
            f"{context}: expected exactly one atom for mask {mask!r}; selected {len(recs)}"
        )


def build_torsions(args, topology: Path) -> List[TorsionSpec]:
    specs: List[TorsionSpec] = []
    requested = [x.strip().lower() for x in args.torsions.split(",") if x.strip()]
    if args.residue and not requested:
        raise RuntimeError("--residue requires at least one name in --torsions")
    residues = inspect_topology_residues(args.cpptraj, topology, args.workdir)

    for label, rmask in args.residue or []:
        rinfo = topology_residue_info(args.cpptraj, topology, rmask, residues)
        base = normalize_resname(rinfo.name)
        if base not in CHI_REGISTRY:
            raise RuntimeError(
                f"{label} ({rmask}) resolved to {rinfo.name} residue {rinfo.number}, "
                "which has no standard chi definition in this program. Use --dihedral "
                "for an explicit four-atom definition."
            )
        available = CHI_REGISTRY[base]
        torsion_names = list(available) if requested == ["all"] else requested
        if "all" in requested and requested != ["all"]:
            raise RuntimeError("Use --torsions all by itself, or list explicit chi names.")
        for torsion in torsion_names:
            if torsion not in available:
                raise RuntimeError(
                    f"{label}: requested {torsion}, but {rinfo.name} supports "
                    f"{','.join(available)}"
                )
            atoms, symmetry = available[torsion]
            masks = tuple(atom_mask(rmask, a) for a in atoms)
            for m, a in zip(masks, atoms):
                validate_single_atom(
                    args.cpptraj, topology, m,
                    f"{label} {torsion} atom {a}"
                )
            specs.append(TorsionSpec(
                observable=label,
                residue_label=label,
                residue_mask=rmask,
                residue_number=rinfo.number,
                residue_name=rinfo.name,
                torsion=torsion,
                atom_names=atoms,
                atom_masks=masks,
                symmetry_period_deg=symmetry,
                classify_rotamer=True,
            ))

    for item in args.dihedral or []:
        name, m1, m2, m3, m4 = item
        masks = (m1, m2, m3, m4)
        for i, m in enumerate(masks, 1):
            validate_single_atom(args.cpptraj, topology, m, f"custom {name} atom{i}")
        specs.append(TorsionSpec(
            observable=name,
            residue_label="",
            residue_mask="",
            residue_number=None,
            residue_name="",
            torsion="custom",
            atom_names=("", "", "", ""),
            atom_masks=masks,
            symmetry_period_deg=None,
            classify_rotamer=bool(args.classify_custom),
        ))

    if not specs:
        raise RuntimeError("No torsions requested. Use --residue and/or --dihedral.")

    keys = [(s.observable, s.torsion) for s in specs]
    if len(set(keys)) != len(keys):
        raise RuntimeError(
            "Duplicate observable/torsion names detected. Give each residue/custom "
            "observable a unique label."
        )
    return specs


def cpptraj_env(threads: int) -> Dict[str, str]:
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = str(threads)
    env["OMP_DYNAMIC"] = "FALSE"
    env["OPENBLAS_NUM_THREADS"] = "1"
    env["MKL_NUM_THREADS"] = "1"
    env["NUMEXPR_NUM_THREADS"] = "1"
    return env


def dataset_name(i: int) -> str:
    return f"D{i:04d}"


def work_tag(path: Path) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", path.stem)


def build_cpptraj_input(
    args, topology: Path, traj: Path, specs: Sequence[TorsionSpec], workdir: Path
) -> Tuple[str, Dict[str, Path]]:
    tag = work_tag(traj)
    outputs: Dict[str, Path] = {}
    lines = [f"parm {qpath(topology)}", f"trajin {qpath(traj)}"]
    if args.image_anchor:
        lines.append(f"autoimage anchor {args.image_anchor}")
    lines.append("")

    for i, spec in enumerate(specs):
        ds = dataset_name(i)
        out = workdir / f"{tag}.{spec.observable}.{spec.torsion}.dat"
        outputs[ds] = out
        m1, m2, m3, m4 = spec.atom_masks
        lines.append(
            f"dihedral {ds} {m1} {m2} {m3} {m4} out {qpath(out)}"
        )
    lines.extend(["", "run", "quit", ""])
    return "\n".join(lines), outputs


def read_cpptraj_series(path: Path) -> List[float]:
    vals: List[float] = []
    with path.open("r", errors="replace") as fh:
        for line in fh:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            toks = s.split()
            # cpptraj dihedral output is normally frame/value; accept last numeric.
            numeric: List[float] = []
            for tok in toks:
                try:
                    numeric.append(float(tok))
                except ValueError:
                    pass
            if len(numeric) >= 2:
                vals.append(numeric[-1])
    if not vals:
        raise RuntimeError(f"No numeric dihedral values read from {path}")
    return vals


def run_trajectory(
    args, topology: Path, traj: Path, ph_regex: re.Pattern,
    specs: Sequence[TorsionSpec], workdir: Path
) -> TrajectoryResult:
    ph = parse_ph(traj, ph_regex)
    print(f"[{traj.stem}] running", flush=True)
    workdir.mkdir(parents=True, exist_ok=True)
    text, outputs = build_cpptraj_input(args, topology, traj, specs, workdir)
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
    lengths = []
    for i, spec in enumerate(specs):
        ds = dataset_name(i)
        out = outputs[ds]
        if not out.exists() or out.stat().st_size == 0:
            raise RuntimeError(f"Missing/empty cpptraj output {out}; inspect {log}")
        series = read_cpptraj_series(out)
        data[ds] = series
        lengths.append(len(series))
    if len(set(lengths)) != 1:
        raise RuntimeError(f"Dihedral series lengths differ for {traj}: {lengths}")

    print(f"[{traj.stem}] complete", flush=True)
    return TrajectoryResult(traj, ph, lengths[0], data)


def wrap180(angle: float) -> float:
    x = (angle + 180.0) % 360.0 - 180.0
    # Avoid a positive 180 endpoint from floating-point noise.
    if x >= 180.0:
        x -= 360.0
    return x


def symmetry_canonical(angle: float, period: int) -> float:
    if period <= 0 or 360 % period != 0:
        raise ValueError(f"Unsupported symmetry period {period}")
    half = period / 2.0
    return (angle + half) % period - half


def rotamer_label(angle: float) -> str:
    """Three-state side-chain rotamer assignment on [-180,180).

    g-    : [-120, 0)
    g+    : [0, 120)
    trans : [120,180) U [-180,-120)

    Boundaries are explicit and deterministic. These labels are intended as a
    coarse conventional rotamer partition; raw angles are always retained.
    """
    a = wrap180(angle)
    if -120.0 <= a < 0.0:
        return "g-"
    if 0.0 <= a < 120.0:
        return "g+"
    return "trans"


def float_text(x: float) -> str:
    # 10 significant decimals is ample relative to cpptraj text precision while
    # keeping files manageable.
    return f"{x:.10g}"


def metric_rows_for_angle(spec: TorsionSpec, angle_raw: float) -> Iterable[Tuple[str, float, str]]:
    angle = wrap180(angle_raw)
    rad = math.radians(angle)
    yield "angle_deg", angle, "degree"
    yield "sin_angle", math.sin(rad), "dimensionless"
    yield "cos_angle", math.cos(rad), "dimensionless"

    if spec.classify_rotamer:
        label = rotamer_label(angle)
        yield "rotamer_gminus", 1.0 if label == "g-" else 0.0, "fraction"
        yield "rotamer_gplus", 1.0 if label == "g+" else 0.0, "fraction"
        yield "rotamer_trans", 1.0 if label == "trans" else 0.0, "fraction"

    if spec.symmetry_period_deg:
        yield (
            "symmetry_angle_deg",
            symmetry_canonical(angle, spec.symmetry_period_deg),
            "degree",
        )


def write_output(
    args, results: Sequence[TrajectoryResult], specs: Sequence[TorsionSpec], output: Path
) -> int:
    fields = [
        "trajectory", "pH", "analysis_row", "trajectory_frame", "step",
        "observable", "residue_label", "residue_mask", "residue_number",
        "residue_name", "torsion", "atom1", "atom2", "atom3", "atom4",
        "symmetry_period_deg", "metric", "value", "unit",
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
                for i, spec in enumerate(specs):
                    angle = result.data[dataset_name(i)][frame0]
                    for metric, value, unit in metric_rows_for_angle(spec, angle):
                        w.writerow({
                            "trajectory": result.trajectory.stem,
                            "pH": f"{result.ph:g}",
                            "analysis_row": analysis_row,
                            "trajectory_frame": trajectory_frame,
                            "step": step,
                            "observable": spec.observable,
                            "residue_label": spec.residue_label,
                            "residue_mask": spec.residue_mask,
                            "residue_number": "" if spec.residue_number is None else spec.residue_number,
                            "residue_name": spec.residue_name,
                            "torsion": spec.torsion,
                            "atom1": spec.atom_masks[0],
                            "atom2": spec.atom_masks[1],
                            "atom3": spec.atom_masks[2],
                            "atom4": spec.atom_masks[3],
                            "symmetry_period_deg": "" if spec.symmetry_period_deg is None else spec.symmetry_period_deg,
                            "metric": metric,
                            "value": float_text(value),
                            "unit": unit,
                        })
                        output_rows += 1
    return output_rows


def write_settings(
    args, topology: Path, trajectories: Sequence[Path], specs: Sequence[TorsionSpec],
    output: Path, nrows: int
) -> None:
    p = Path(str(output) + ".settings.txt")
    lines = [
        f"program=calc_dihedral.py",
        f"version={PROGRAM_VERSION}",
        f"cpptraj={args.cpptraj}",
        f"topology={topology}",
        f"n_trajectories={len(trajectories)}",
        f"n_torsions={len(specs)}",
        f"first_frame_step={args.first_frame_step}",
        f"frame_step_interval={args.frame_step_interval}",
        f"image_anchor={args.image_anchor or ''}",
        f"output_rows={nrows}",
        "angle_wrap=[-180,180)",
        "rotamer_gminus=[-120,0)",
        "rotamer_gplus=[0,120)",
        "rotamer_trans=[120,180)U[-180,-120)",
        "circular_mean_formula=atan2(mean_sin_angle,mean_cos_angle)",
        "",
        "torsions:",
    ]
    for s in specs:
        lines.append(
            "\t".join([
                s.observable, s.torsion, s.residue_name,
                "-".join(s.atom_names) if any(s.atom_names) else "custom",
                "symmetry=" + (str(s.symmetry_period_deg) if s.symmetry_period_deg else "none"),
            ])
        )
    p.write_text("\n".join(lines) + "\n")


def make_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="General Amber/cpptraj protein dihedral and rotamer analyzer.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("-p", "--topology", nargs="+", required=True,
                   help="Topology path/glob. Must resolve to exactly one file.")
    p.add_argument("-t", "--trajectory", nargs="+", required=True,
                   help="Trajectory paths/globs.")
    p.add_argument(
        "--residue", nargs=2, action="append", metavar=("LABEL", "MASK"),
        help="Analyze standard chi torsions for one residue; repeatable."
    )
    p.add_argument(
        "--torsions", default="chi1,chi2",
        help="Comma-separated standard chi names applied to every --residue, or 'all'."
    )
    p.add_argument(
        "--dihedral", nargs=5, action="append",
        metavar=("NAME", "MASK1", "MASK2", "MASK3", "MASK4"),
        help="Arbitrary explicit four-atom dihedral; repeatable."
    )
    p.add_argument(
        "--classify-custom", action="store_true",
        help="Also apply coarse g-/g+/trans indicators to custom dihedrals."
    )
    p.add_argument("--image-anchor", default=None,
                   help="Optional cpptraj autoimage anchor mask, e.g. '^1'.")
    p.add_argument("--first-frame-step", type=int, required=True,
                   help="MD step associated with trajectory frame 1.")
    p.add_argument("--frame-step-interval", type=int, required=True,
                   help="MD step increment between consecutive trajectory frames.")
    p.add_argument(
        "--ph-regex", default=r"ph(?P<ph>\d+(?:[_\.]\d+)?)",
        help="Regex containing named group 'ph' for parsing pH from trajectory name."
    )
    p.add_argument("--cpptraj", default="cpptraj", help="cpptraj executable.")
    p.add_argument("--cpptraj-threads", type=int, default=1,
                   help="OMP threads per cpptraj worker.")
    p.add_argument("--jobs", type=int, default=1,
                   help="Number of trajectories processed concurrently.")
    p.add_argument("--workdir", type=Path, default=Path("dihedral_work"),
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

    specs = build_torsions(args, topology)
    args.workdir.mkdir(parents=True, exist_ok=True)
    args.output.parent.mkdir(parents=True, exist_ok=True)

    results: List[TrajectoryResult] = []
    errors: List[str] = []
    with ThreadPoolExecutor(max_workers=args.jobs) as ex:
        futs = {
            ex.submit(run_trajectory, args, topology, traj, ph_regex, specs, args.workdir): traj
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
    nrows = write_output(args, results, specs, args.output)
    write_settings(args, topology, trajectories, specs, args.output, nrows)

    print()
    print(f"Trajectories analyzed       : {len(results)}")
    print(f"Requested torsions          : {len(specs)}")
    print(f"Trajectory frames           : {nframes}")
    print(f"Frame-wise output rows      : {nrows}")
    print(f"Program version             : {PROGRAM_VERSION}")
    print("Angle convention            : [-180, 180)")
    print("Rotamer indicators          : g-, g+, trans")
    print("Circular components         : sin_angle, cos_angle")
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
