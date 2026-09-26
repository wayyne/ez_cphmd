#!/usr/bin/env python3
"""
calc_hbond_v5_1.py

Lambda-aware hydrogen-bond analysis for Amber continuous constant-pH MD.

Design
------
CPPTRAJ is used only as the geometric candidate detector.  For broad protein
scan/focus/interaction modes, v5 first asks CPPTRAJ to expand each user mask
against the topology, enumerates every candidate donor heavy atom explicitly,
and runs one H-bond action per donor against an explicit acceptor atom list.
This avoids incomplete per-bond series that can occur when many donors and
acceptors are packed into one broad hbond action.  Every detected
donor-H-acceptor candidate has a 0/1 time series.  Python then exact-joins each
analyzed trajectory frame to the corresponding lambda-file MD step and applies
frame-specific protonation/tautomer chemistry before counting H-bonds.

This matters for CpHMD because a titratable atom may be a donor in one state
and an acceptor in another.  A fixed topology alone cannot make that decision.

Built-in CpH residue models
---------------------------
The default registry covers the titratable amino-acid models used in Amber /
Shen-lab continuous CpHMD work:

  ASP-family : ASP, ASH, AS2, AS4
  GLU-family : GLU, GLH, GL2, GL4
  HIS-family : HIS, HID, HIE, HIP
  CYS-family : CYS, CYM
  LYS-family : LYS, LYN
  TYR-family : TYR, TYM

Coordinate conventions used by the built-in registry:

  Asp/Glu: itauto=3 protonation lambda; itauto=4 tautomer coordinate.
           protonation lambda low  -> protonated
           protonation lambda high -> deprotonated
           when protonated: tautomer low -> O2-H; high -> O1-H

  His:     itauto=1 protonation lambda; itauto=2 tautomer coordinate.
           protonation lambda low  -> HIP (both ring nitrogens protonated)
           protonation lambda high -> neutral
           when neutral: tautomer low -> HID; high -> HIE

  Cys/Lys/Tyr: one protonation coordinate; low -> protonated,
               high -> deprotonated.

Lambda values between --lambda-low and --lambda-high are treated as mixed.
For Asp/Glu and neutral His, a mixed tautomer coordinate is also ambiguous.

The implementation is deliberately registry-driven: residue-state chemistry
is centralized in `site_roles_for_frame()` rather than scattered through the
analysis modes.  Additional titratable residue models can be added there
without changing the H-bond machinery.

Exact frame/lambda synchronization
----------------------------------
With --first-frame-step and --frame-step-interval, analyzed frame k maps to
an exact MD step.  Lambda rows are joined by exact step; nearest-neighbor
matching is never performed.  Missing steps inside the lambda-file step range
are fatal.  Frames outside that range are retained.  If a *present* H-bond candidate
requires unavailable lambda data, it is counted in `outside_lambda_hbond_count`
and the observable/frame is marked `chemistry_status=outside_lambda_window`;
definitely valid bonds in that same frame are still retained.

Output
------
Canonical frame-wise TSV fields:

  trajectory, pH, analysis_row, trajectory_frame, step,
  observable, metric, interaction_class, mode, selection1, selection2,
  value, unit, chemistry_status, lambda_file

Each observable emits:

  hbond_count                 number of definitely chemically valid H-bonds
  hbond_present               1 if at least one definitely valid H-bond is present
  ambiguous_hbond_count       present candidates unresolved by mixed lambda/tautomer
  outside_lambda_hbond_count  present candidates requiring lambda outside available range
  hbond_count_min             lower bound; identical to hbond_count
  hbond_count_max             upper bound if every unresolved candidate were chemically valid

Known information is never erased because another candidate is unresolved.
`chemistry_status` is `clean`, `partial_ambiguity`, or
`outside_lambda_window`.  `outside_lambda_window` takes precedence if at least
one present candidate requires unavailable lambda data.  Otherwise any mixed
lambda/tautomer candidate gives `partial_ambiguity`.

The interaction report `<stem>.interactions.tsv` is reconstructed from the
per-bond time series after CpH filtering.  `count` and `fraction` are therefore
chemically filtered.  CPPTRAJ distance/angle averages are retained as
`raw_avg_distance_A` and `raw_avg_angle_deg`; they are averages over geometric
candidate-present frames before lambda filtering and are explicitly labelled
as such.

Examples
--------

  python3 calc_hbond_v5_1.py \
    -p top.prmtop \
    -t 'arex.ph6_5.nc' \
    --lambda 'rex/arex.ph6_5.lambda' \
    --scan PROTEIN '^1' \
    --focus DYAD ':36,53' '^1' \
    --image-anchor '^1' \
    --first-frame-step 510000 \
    --frame-step-interval 10000 \
    --jobs 1 \
    --workdir hbond_work_v5_1 \
    -o hewl_hbonds_v5_1.tsv

For fixed-protonation controls, use --no-cph-filter.  In normal CpHMD use,
--lambda is required.
"""

from __future__ import annotations

import argparse
import csv
import glob
import hashlib
import json
import math
import os
import re
import shutil
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Optional, Sequence

PROGRAM_VERSION = "5.1"
SETTINGS_VERSION = 6

OUTPUT_FIELDS = [
    "trajectory", "pH", "analysis_row", "trajectory_frame", "step",
    "observable", "metric", "interaction_class", "mode",
    "selection1", "selection2", "value", "unit",
    "chemistry_status", "lambda_file",
]

DETAIL_FIELDS = [
    "trajectory", "pH", "observable", "mode", "interaction_class",
    "selection1", "selection2", "acceptor", "donor_h", "donor",
    "count", "count_semantics", "fraction", "raw_count",
    "raw_fraction", "ambiguous_present_frames",
    "outside_lambda_present_frames", "raw_avg_distance_A",
    "raw_avg_angle_deg", "geometry_average_scope",
]


@dataclass(frozen=True)
class UserObservable:
    name: str
    mode: str
    interaction_class: str
    selection1: str
    selection2: str
    donor_h_mask: str = ""


@dataclass(frozen=True)
class ActionSpec:
    action_id: str
    parent_name: str
    role: str
    kind: str                 # UU or UV
    command: str
    count_file: Path
    avg_file: Path
    series_file: Path


@dataclass(frozen=True)
class AtomRef:
    raw: str
    resname: str
    resid: int
    atom: str


@dataclass(frozen=True)
class TopologyAtom:
    index: int
    name: str
    resid: int
    resname: str

    @property
    def atom_ref(self) -> AtomRef:
        raw = f"{self.resname}_{self.resid}@{self.name}"
        return AtomRef(raw, self.resname.upper(), self.resid, self.name.upper())


@dataclass(frozen=True)
class BondSeries:
    key: tuple[str, str, str]
    acceptor: str
    donor_h: str
    donor: str
    values: tuple[int, ...]
    kind: str
    # For UV rows, one side can be V.  The non-V role is represented in
    # acceptor/donor/donor_h as far as CPPTRAJ exposes it.


@dataclass(frozen=True)
class LambdaFile:
    path: Path
    ires: tuple[str, ...]
    itauto: tuple[str, ...]
    step_to_values: dict[int, tuple[float, ...]]
    min_step: int
    max_step: int


# ---------------------------------------------------------------------------
# Generic helpers
# ---------------------------------------------------------------------------

def natural_key(value: str | Path):
    return [int(x) if x.isdigit() else x.lower()
            for x in re.split(r"([0-9]+)", str(value))]


def safe_component(value: str) -> str:
    text = re.sub(r"[^A-Za-z0-9_.-]+", "_", value.strip()).strip("._")
    if not text:
        raise ValueError(f"Invalid empty name derived from {value!r}")
    return text


def expand_specs(specs: Iterable[str], what: str) -> list[Path]:
    found: dict[str, Path] = {}
    for spec in specs:
        matches = glob.glob(spec)
        if matches:
            for match in matches:
                p = Path(match).resolve()
                if p.is_file():
                    found[str(p)] = p
            continue
        p = Path(spec).resolve()
        if p.is_file():
            found[str(p)] = p
        else:
            raise FileNotFoundError(f"No {what} file matched: {spec}")
    out = sorted(found.values(), key=natural_key)
    if not out:
        raise RuntimeError(f"No {what} files were found")
    return out


def extract_ph(name: str) -> float:
    m = re.search(r"[pP][hH][_=-]?([0-9]+(?:[._][0-9]+)?)", name)
    return float(m.group(1).replace("_", ".")) if m else math.nan


def qmask(mask: str) -> str:
    return '"' + mask.replace('"', r'\"') + '"'


def qpath(path: Path) -> str:
    return '"' + str(path).replace('"', r'\"') + '"'


def format_float(x: float) -> str:
    return f"{x:.10g}"


def stable_hash(payload: dict) -> str:
    return hashlib.sha256(json.dumps(payload, sort_keys=True,
                                     separators=(",", ":")).encode()).hexdigest()


def write_tsv(path: Path, rows: list[dict], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t",
                           extrasaction="ignore")
        w.writeheader()
        for row in rows:
            out = dict(row)
            for k, v in list(out.items()):
                if isinstance(v, float):
                    out[k] = "" if not math.isfinite(v) else f"{v:.10g}"
                elif v is None:
                    out[k] = ""
            w.writerow(out)


def parse_stop(stop: str) -> str:
    if stop == "last":
        return stop
    try:
        n = int(stop)
    except ValueError as exc:
        raise ValueError("--stop must be a positive integer or 'last'") from exc
    if n < 1:
        raise ValueError("--stop must be >= 1")
    return str(n)


# ---------------------------------------------------------------------------
# Protein candidate masks
# ---------------------------------------------------------------------------
# These masks deliberately form a SUPERSET of possible chemistry.  Lambda
# filtering later decides whether a titratable atom's role is active.

FIXED_DONOR = (
    "@N"
    "|:ARG@NE,NH1,NH2"
    "|:ASN@ND2"
    "|:GLN@NE2"
    "|:TRP@NE1"
    "|:SER@OG"
    "|:THR@OG1"
)

FIXED_ACCEPTOR = (
    "@O,OXT"
    "|:ASN@OD1"
    "|:GLN@OE1"
    "|:SER@OG"
    "|:THR@OG1"
)

# All roles that can become active in at least one CpH state.
TITRATABLE_DONOR = (
    "|:ASP,ASH,AS2,AS4@OD1,OD2"
    "|:GLU,GLH,GL2,GL4@OE1,OE2"
    "|:HIS,HID,HIE,HIP@ND1,NE2"
    "|:CYS,CYM@SG"
    "|:LYS,LYN@NZ"
    "|:TYR,TYM@OH"
)

TITRATABLE_ACCEPTOR = (
    "|:ASP,ASH,AS2,AS4@OD1,OD2"
    "|:GLU,GLH,GL2,GL4@OE1,OE2"
    "|:HIS,HID,HIE,HIP@ND1,NE2"
    "|:CYS,CYM@SG"
    "|:LYS,LYN@NZ"
    "|:TYR,TYM@OH"
)


def protein_candidate_donor_mask(mask: str, include_met_sulfur: bool) -> str:
    chemistry = FIXED_DONOR + TITRATABLE_DONOR
    # Met sulfur is not proton-titratable and is not a conventional protein
    # H-bond donor; no donor addition is made here.
    return f"({mask})&({chemistry})"


def protein_candidate_acceptor_mask(mask: str, include_met_sulfur: bool) -> str:
    chemistry = FIXED_ACCEPTOR + TITRATABLE_ACCEPTOR
    if include_met_sulfur:
        chemistry += "|:MET@SD"
    return f"({mask})&({chemistry})"


def fon_heavy_mask(mask: str) -> str:
    return f"({mask})&(@/F|@/O|@/N)"


def candidate_donor_mask(mask: str, args: argparse.Namespace) -> str:
    if args.chemistry_profile == "protein":
        return protein_candidate_donor_mask(mask, args.include_met_sulfur)
    return fon_heavy_mask(mask)


def candidate_acceptor_mask(mask: str, args: argparse.Namespace) -> str:
    if args.chemistry_profile == "protein":
        return protein_candidate_acceptor_mask(mask, args.include_met_sulfur)
    return fon_heavy_mask(mask)


def explicit_atom_index_mask(indices: Sequence[int]) -> str:
    vals = sorted({int(i) for i in indices})
    if not vals:
        raise ValueError("Cannot construct an explicit CPPTRAJ mask from zero atoms")
    return "@" + ",".join(str(i) for i in vals)


def cpptraj_mask_atoms(
    topology: Path,
    mask: str,
    args: argparse.Namespace,
    cache: dict[str, tuple[TopologyAtom, ...]],
) -> tuple[TopologyAtom, ...]:
    """Expand an Amber mask to explicit topology atoms with CPPTRAJ itself.

    Using CPPTRAJ for mask expansion means arbitrary Amber masks (^ molecule,
    residue ranges, atom-name masks, boolean expressions, etc.) retain native
    semantics; Python never attempts to reimplement the mask language.
    """
    if mask in cache:
        return cache[mask]
    proc = subprocess.run(
        [args.cpptraj, "-p", str(topology.resolve()), "--mask", mask],
        text=True, capture_output=True,
    )
    if proc.returncode != 0:
        raise RuntimeError(
            f"CPPTRAJ could not expand mask {mask!r}:\n{proc.stdout}\n{proc.stderr}"
        )
    atoms: list[TopologyAtom] = []
    # Typical --mask/atominfo row:
    #   218 N 13 NHE 1 N -0.4630 14.0100 ... N ...
    # We only require the stable first four fields: atom#, atom name,
    # residue#, residue name.
    for raw in proc.stdout.splitlines():
        parts = raw.split()
        if len(parts) < 4:
            continue
        try:
            index = int(parts[0])
            resid = int(parts[2])
        except ValueError:
            continue
        name = parts[1].upper()
        resname = parts[3].upper()
        atoms.append(TopologyAtom(index, name, resid, resname))
    # Some CPPTRAJ builds may send atominfo-style output to stderr.
    if not atoms:
        for raw in proc.stderr.splitlines():
            parts = raw.split()
            if len(parts) < 4:
                continue
            try:
                index = int(parts[0])
                resid = int(parts[2])
            except ValueError:
                continue
            atoms.append(TopologyAtom(index, parts[1].upper(), resid, parts[3].upper()))
    cache[mask] = tuple(atoms)
    return cache[mask]


def enumerated_role_atoms(
    topology: Path,
    user_mask: str,
    role: str,
    args: argparse.Namespace,
    cache: dict[str, tuple[TopologyAtom, ...]],
) -> tuple[TopologyAtom, ...]:
    """Return explicit candidate heavy atoms for one donor/acceptor role."""
    if role == "donor":
        mask = candidate_donor_mask(user_mask, args)
    elif role == "acceptor":
        mask = candidate_acceptor_mask(user_mask, args)
    else:
        raise AssertionError(role)
    return cpptraj_mask_atoms(topology, mask, args, cache)


# ---------------------------------------------------------------------------
# Lambda parsing / exact matching
# ---------------------------------------------------------------------------

def _residue_matches(value: str, resid: str | int) -> bool:
    try:
        return int(float(value)) == int(float(resid))
    except (ValueError, TypeError):
        return str(value) == str(resid)


def read_lambda_file(path: Path) -> LambdaFile:
    ires: list[str] = []
    itauto: list[str] = []
    expected: Optional[int] = None
    rows: dict[int, tuple[float, ...]] = {}

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
                raise RuntimeError(f"{path}:{lineno}: invalid MD step") from exc
            if abs(step_f - step) > 1e-6:
                raise RuntimeError(f"{path}:{lineno}: non-integral MD step {step_f}")
            try:
                vals = tuple(float(x) for x in parts[1:])
            except ValueError as exc:
                raise RuntimeError(f"{path}:{lineno}: invalid lambda value") from exc
            if expected is not None and len(vals) != expected:
                raise RuntimeError(
                    f"{path}:{lineno}: expected {expected} lambda columns, found {len(vals)}"
                )
            if step in rows:
                raise RuntimeError(f"{path}:{lineno}: duplicate lambda step {step}")
            rows[step] = vals

    if expected is None:
        raise RuntimeError(f"{path}: missing '# ititr' header")
    if not ires:
        raise RuntimeError(f"{path}: missing '# ires' header")

    # Preserve compatibility with compact ires headers used by some files.
    if len(ires) != expected:
        compact = "".join(ires)
        if not compact.isdigit() or len(compact) % expected != 0:
            raise RuntimeError(
                f"{path}: cannot interpret ires header ({len(ires)} tokens for {expected} columns)"
            )
        width = len(compact) // expected
        ires = [compact[i:i+width] for i in range(0, len(compact), width)]

    if len(ires) != expected:
        raise RuntimeError(f"{path}: ires column count mismatch")
    if itauto and len(itauto) != expected:
        raise RuntimeError(f"{path}: itauto column count mismatch")
    if not itauto:
        itauto = [""] * expected
    if not rows:
        raise RuntimeError(f"{path}: no lambda data rows")

    steps = sorted(rows)
    return LambdaFile(path.resolve(), tuple(ires), tuple(itauto), rows,
                      steps[0], steps[-1])


def lambda_columns_for_residue(lam: LambdaFile, resid: int) -> list[int]:
    return [i for i, r in enumerate(lam.ires) if _residue_matches(r, resid)]


def itauto_int(lam: LambdaFile, col: int) -> Optional[int]:
    if col >= len(lam.itauto):
        return None
    try:
        return int(float(lam.itauto[col]))
    except (ValueError, TypeError):
        return None


def lambda_row_exact(lam: LambdaFile, step: int) -> tuple[str, Optional[tuple[float, ...]]]:
    if step in lam.step_to_values:
        return "matched", lam.step_to_values[step]
    if step < lam.min_step or step > lam.max_step:
        return "outside_lambda_window", None
    raise RuntimeError(
        f"{lam.path.name}: trajectory step {step} lies inside lambda range "
        f"{lam.min_step}:{lam.max_step} but has no exact lambda row"
    )


def match_lambda_files(trajectories: list[Path], lambda_files: list[Path]) -> dict[Path, Path]:
    if not lambda_files:
        return {}
    by_stem = {p.stem: p for p in lambda_files}
    mapping: dict[Path, Path] = {}
    for traj in trajectories:
        if traj.stem in by_stem:
            mapping[traj] = by_stem[traj.stem]
            continue
        tph = extract_ph(traj.stem)
        ph_matches = [p for p in lambda_files
                      if math.isfinite(tph) and math.isfinite(extract_ph(p.stem))
                      and abs(extract_ph(p.stem) - tph) < 1e-8]
        if len(ph_matches) == 1:
            mapping[traj] = ph_matches[0]
            continue
        if len(lambda_files) == 1 and len(trajectories) == 1:
            mapping[traj] = lambda_files[0]
            continue
        raise RuntimeError(
            f"Could not uniquely match trajectory {traj.name} to a lambda file"
        )
    return mapping


# ---------------------------------------------------------------------------
# Atom-label parsing and CpH chemistry registry
# ---------------------------------------------------------------------------

ATOM_RE = re.compile(
    r"^(?P<resname>[A-Za-z0-9+]+)_(?P<resid>[0-9]+)@(?P<atom>[^\s]+)$"
)


def parse_atom_ref(text: str) -> Optional[AtomRef]:
    if text in {"", "V", "v", "SOLVENT"}:
        return None
    m = ATOM_RE.match(text.strip())
    if not m:
        return None
    return AtomRef(text.strip(), m.group("resname").upper(),
                   int(m.group("resid")), m.group("atom").upper())


def canonical_family(resname: str) -> Optional[str]:
    name = resname.upper()
    if name in {"ASP", "ASH", "AS2", "AS4"}:
        return "ASP"
    if name in {"GLU", "GLH", "GL2", "GL4"}:
        return "GLU"
    if name in {"HIS", "HID", "HIE", "HIP"}:
        return "HIS"
    if name in {"CYS", "CYM"}:
        return "CYS"
    if name in {"LYS", "LYN"}:
        return "LYS"
    if name in {"TYR", "TYM"}:
        return "TYR"
    return None


def _threshold(x: float, low: float, high: float) -> str:
    if x <= low:
        return "low"
    if x >= high:
        return "high"
    return "mixed"


def _columns_by_flag(lam: LambdaFile, resid: int) -> dict[int, list[int]]:
    out: dict[int, list[int]] = {}
    for col in lambda_columns_for_residue(lam, resid):
        flag = itauto_int(lam, col)
        if flag is not None:
            out.setdefault(flag, []).append(col)
    return out


def site_roles_for_frame(
    *,
    atom: AtomRef,
    lam: LambdaFile,
    values: tuple[float, ...],
    low: float,
    high: float,
) -> tuple[str, set[str], set[str]]:
    """Return (status, donor-heavy-atoms, acceptor-atoms) for one CpH site.

    Status is `clean` or `mixed`.  This routine is the central CpH chemistry
    registry.  It is intentionally residue-model driven rather than HEWL- or
    Asp/Glu-specific.
    """
    family = canonical_family(atom.resname)
    cols = lambda_columns_for_residue(lam, atom.resid)
    if not cols:
        # Residue is not titrated in this lambda file.  The caller should use
        # fixed chemistry for it rather than this function.
        return "not_titratable", set(), set()
    if family is None:
        raise RuntimeError(
            f"Lambda file marks topology residue {atom.resid} as titratable, "
            f"but residue name {atom.resname!r} has no built-in CpH chemistry model"
        )

    flags = _columns_by_flag(lam, atom.resid)

    if family in {"ASP", "GLU"}:
        pcols = flags.get(3, [])
        tcols = flags.get(4, [])
        if len(pcols) != 1 or len(tcols) != 1:
            raise RuntimeError(
                f"Residue {atom.resid} {atom.resname}: expected one itauto=3 "
                f"protonation and one itauto=4 tautomer coordinate; found "
                f"itauto3={pcols}, itauto4={tcols}"
            )
        pstate = _threshold(values[pcols[0]], low, high)
        if family == "ASP":
            o1, o2 = "OD1", "OD2"
        else:
            o1, o2 = "OE1", "OE2"
        if pstate == "high":
            return "clean", set(), {o1, o2}
        if pstate == "mixed":
            return "mixed", set(), set()
        tstate = _threshold(values[tcols[0]], low, high)
        if tstate == "mixed":
            return "mixed", set(), set()
        # Shen/Khandogin convention: tautomer low -> O2-H; high -> O1-H.
        if tstate == "low":
            return "clean", {o2}, {o1}
        return "clean", {o1}, {o2}

    if family == "HIS":
        pcols = flags.get(1, [])
        tcols = flags.get(2, [])
        if len(pcols) != 1 or len(tcols) != 1:
            raise RuntimeError(
                f"Residue {atom.resid} {atom.resname}: expected one itauto=1 "
                f"protonation and one itauto=2 tautomer coordinate"
            )
        pstate = _threshold(values[pcols[0]], low, high)
        if pstate == "low":
            return "clean", {"ND1", "NE2"}, set()  # HIP
        if pstate == "mixed":
            return "mixed", set(), set()
        tstate = _threshold(values[tcols[0]], low, high)
        if tstate == "mixed":
            return "mixed", set(), set()
        if tstate == "low":
            return "clean", {"ND1"}, {"NE2"}        # HID
        return "clean", {"NE2"}, {"ND1"}            # HIE

    # One-coordinate models.  If duplicates exist, require a unique
    # protonation-like column (itauto=1 or 3); otherwise a single column.
    if len(cols) == 1:
        pcol = cols[0]
    else:
        preferred = flags.get(3, []) or flags.get(1, [])
        if len(preferred) != 1:
            raise RuntimeError(
                f"Residue {atom.resid} {atom.resname}: cannot identify unique "
                f"protonation coordinate among columns {cols}"
            )
        pcol = preferred[0]
    state = _threshold(values[pcol], low, high)
    if state == "mixed":
        return "mixed", set(), set()

    if family == "CYS":
        # Protonated thiol: donor.  Deprotonated thiolate: acceptor.
        return ("clean", {"SG"}, set()) if state == "low" else ("clean", set(), {"SG"})
    if family == "LYS":
        # LysH+ is donor/non-acceptor; neutral Lys is donor and acceptor.
        return ("clean", {"NZ"}, set()) if state == "low" else ("clean", {"NZ"}, {"NZ"})
    if family == "TYR":
        # Phenol is conventionally allowed as both donor and acceptor in the
        # protein H-bond profile; phenolate is acceptor only.
        return ("clean", {"OH"}, {"OH"}) if state == "low" else ("clean", set(), {"OH"})

    raise AssertionError(family)


def fixed_role_valid(atom: AtomRef, role: str, include_met_sulfur: bool) -> bool:
    """Chemistry for residues not titrated in the matched lambda file."""
    rn, a = atom.resname, atom.atom
    if role == "donor":
        if a == "N": return True
        if rn == "ARG" and a in {"NE", "NH1", "NH2"}: return True
        if rn == "ASN" and a == "ND2": return True
        if rn == "GLN" and a == "NE2": return True
        if rn == "TRP" and a == "NE1": return True
        if rn == "SER" and a == "OG": return True
        if rn == "THR" and a == "OG1": return True
        # Fixed protonation residue aliases.
        if rn in {"ASH", "AS2", "AS4"} and a in {"OD1", "OD2"}: return True
        if rn in {"GLH", "GL2", "GL4"} and a in {"OE1", "OE2"}: return True
        if rn in {"HID", "HIE", "HIP", "HIS"} and a in {"ND1", "NE2"}: return True
        if rn == "CYS" and a == "SG": return True
        if rn in {"LYS", "LYN"} and a == "NZ": return True
        if rn == "TYR" and a == "OH": return True
        return False
    if role == "acceptor":
        if a in {"O", "OXT"}: return True
        if rn in {"ASP", "ASH", "AS2", "AS4"} and a in {"OD1", "OD2"}: return True
        if rn in {"GLU", "GLH", "GL2", "GL4"} and a in {"OE1", "OE2"}: return True
        if rn == "ASN" and a == "OD1": return True
        if rn == "GLN" and a == "OE1": return True
        if rn == "SER" and a == "OG": return True
        if rn == "THR" and a == "OG1": return True
        if rn in {"HID", "HIE", "HIS"} and a in {"ND1", "NE2"}: return True
        if rn in {"CYS", "CYM"} and a == "SG": return True
        if rn in {"LYN"} and a == "NZ": return True
        if rn in {"TYR", "TYM"} and a == "OH": return True
        if include_met_sulfur and rn == "MET" and a == "SD": return True
        return False
    raise AssertionError(role)


def is_cph_role_atom(atom: AtomRef) -> bool:
    """Return True only for atoms whose donor/acceptor role changes with CpH state.

    A residue can be titratable while most of its atoms are not.  In particular,
    backbone N/O atoms on ASP/GLU/HIS/CYS/LYS/TYR must retain ordinary fixed
    protein chemistry and must never become lambda-ambiguous merely because the
    side chain titrates.
    """
    family = canonical_family(atom.resname)
    a = atom.atom
    if family == "ASP":
        return a in {"OD1", "OD2"}
    if family == "GLU":
        return a in {"OE1", "OE2"}
    if family == "HIS":
        return a in {"ND1", "NE2"}
    if family == "CYS":
        return a == "SG"
    if family == "LYS":
        return a == "NZ"
    if family == "TYR":
        return a == "OH"
    return False


def role_valid(
    atom: Optional[AtomRef],
    role: str,
    lam: Optional[LambdaFile],
    values: Optional[tuple[float, ...]],
    args: argparse.Namespace,
) -> str:
    """Return valid / invalid / mixed / outside for one heavy-atom role.

    Lambda filtering is applied only to atoms whose chemical donor/acceptor role
    actually changes with protonation/tautomer state.  All other atoms, including
    backbone atoms belonging to titratable residues, use fixed protein chemistry.
    """
    if atom is None:
        return "valid"  # solvent V endpoint; chemistry handled on solute side
    if args.no_cph_filter or lam is None:
        return "valid" if fixed_role_valid(atom, role, args.include_met_sulfur) else "invalid"

    tit_cols = lambda_columns_for_residue(lam, atom.resid)
    if not tit_cols or not is_cph_role_atom(atom):
        return "valid" if fixed_role_valid(atom, role, args.include_met_sulfur) else "invalid"
    if values is None:
        return "outside"
    status, donors, acceptors = site_roles_for_frame(
        atom=atom, lam=lam, values=values,
        low=args.lambda_low, high=args.lambda_high,
    )
    if status == "mixed":
        return "mixed"
    allowed = donors if role == "donor" else acceptors
    return "valid" if atom.atom in allowed else "invalid"


def donor_heavy_from_hydrogen(h: AtomRef) -> Optional[AtomRef]:
    """Infer titratable donor heavy atom from a CpH donor-H atom label.

    Needed for non-specific solute-solvent CPPTRAJ series, whose legend stores
    the solute donor H rather than the donor heavy atom when solvent is the
    acceptor.
    """
    fam = canonical_family(h.resname)
    a = h.atom
    heavy: Optional[str] = None
    if fam == "ASP":
        if a == "HD1": heavy = "OD1"
        elif a == "HD2": heavy = "OD2"
    elif fam == "GLU":
        if a == "HE1": heavy = "OE1"
        elif a == "HE2": heavy = "OE2"
    elif fam == "HIS":
        if a == "HD1": heavy = "ND1"
        elif a == "HE2": heavy = "NE2"
    elif fam == "CYS" and a in {"HG", "HSG"}:
        heavy = "SG"
    elif fam == "TYR" and a in {"HH", "HO"}:
        heavy = "OH"
    elif fam == "LYS" and a.startswith("HZ"):
        heavy = "NZ"
    if heavy is None:
        return None
    raw = f"{h.resname}_{h.resid}@{heavy}"
    return AtomRef(raw, h.resname, h.resid, heavy)


# ---------------------------------------------------------------------------
# CPPTRAJ output parsing
# ---------------------------------------------------------------------------

def read_hbond_count_series(path: Path, aspect: str) -> list[int]:
    if not path.is_file():
        raise RuntimeError(f"Missing CPPTRAJ count file: {path}")
    header: Optional[list[str]] = None
    data: list[list[str]] = []
    with path.open() as fh:
        for line in fh:
            s = line.strip()
            if not s: continue
            if s.startswith("#"):
                toks = s.lstrip("#").split()
                if toks and toks[0].lower() == "frame": header = toks
                continue
            data.append(s.split())
    if not data:
        raise RuntimeError(f"No count data parsed from {path}")
    wanted = f"[{aspect}]"
    ci = None
    if header:
        for i, tok in enumerate(header):
            if wanted in tok:
                ci = i; break
    if ci is None:
        ci = 1 if aspect == "UU" else 2
    out: list[int] = []
    for row in data:
        if len(row) <= ci:
            raise RuntimeError(f"Malformed count row in {path}")
        x = float(row[ci]); n = int(round(x))
        if abs(x-n) > 1e-6 or n < 0:
            raise RuntimeError(f"Non-integral/negative hbond count {x} in {path}")
        out.append(n)
    return out


def read_hbond_average_report(path: Path) -> dict[tuple[str,str,str], dict]:
    rows: dict[tuple[str,str,str], dict] = {}
    if not path.is_file():
        return rows
    with path.open() as fh:
        for line in fh:
            s = line.strip()
            if not s or s.startswith("#"): continue
            f = s.split()
            if len(f) < 7: continue
            try:
                count = int(round(float(f[-4])))
                frac = float(f[-3]); dist = float(f[-2]); ang = float(f[-1])
            except ValueError:
                continue
            acc, dh, don = f[0], f[1], f[2]
            rows[(acc, dh, don)] = {
                "raw_count": count,
                "raw_fraction": frac,
                "raw_avg_distance_A": dist,
                "raw_avg_angle_deg": ang,
            }
    return rows


def _human_uu_legend_to_triplet(legend: str) -> Optional[tuple[str,str,str]]:
    # CPPTRAJ CreateHBlegend: acceptor-full - donor-full - donorH-atom-name
    parts = legend.strip().strip('"').split("-")
    if len(parts) < 3:
        return None
    acc = parts[0]
    donor = parts[1]
    hatom = "-".join(parts[2:])
    dref = parse_atom_ref(donor)
    if dref is None:
        return None
    dh = f"{dref.resname}_{dref.resid}@{hatom.upper()}"
    if parse_atom_ref(acc) is None:
        return None
    return (acc, dh, donor)


def read_series_matrix(path: Path, nframes: int, kind: str) -> list[tuple[str, tuple[int,...]]]:
    """Read CPPTRAJ generic `.dat` per-H-bond series.

    We intentionally write the series with a `.dat` extension so CPPTRAJ uses
    its plain column format.  The `#Frame` header contains data-set legends.
    """
    if not path.is_file() or path.stat().st_size == 0:
        return []
    header: Optional[list[str]] = None
    numeric: list[list[str]] = []
    with path.open() as fh:
        for raw in fh:
            s = raw.strip()
            if not s: continue
            if s.startswith("#"):
                toks = s.lstrip("#").split()
                if toks and toks[0].lower() == "frame":
                    header = toks
                continue
            # Ignore xmgrace directives defensively; .dat should not have them.
            if s.startswith("@"): continue
            fields = s.split()
            try:
                float(fields[0])
            except (ValueError, IndexError):
                continue
            numeric.append(fields)
    if not numeric:
        return []
    if len(numeric) != nframes:
        raise RuntimeError(
            f"{path}: series has {len(numeric)} frames, expected {nframes}"
        )
    ncol = len(numeric[0])
    if any(len(r) != ncol for r in numeric):
        raise RuntimeError(f"{path}: inconsistent series column counts")
    if ncol == 1:
        return []
    if header is None or len(header) != ncol:
        raise RuntimeError(
            f"{path}: expected a #Frame header with {ncol} columns; found "
            f"{0 if header is None else len(header)}. Keep the uuseries/uvseries "
            "filename extension as .dat."
        )
    out: list[tuple[str, tuple[int,...]]] = []
    for ci in range(1, ncol):
        vals: list[int] = []
        for row in numeric:
            x = float(row[ci])
            n = int(round(x))

            if abs(x - n) > 1e-6 or n < 0:
                raise RuntimeError(
                    f"{path}: non-integral/negative H-bond series value {x}"
                )

            # Solute-solute (UU) series represent a specific physical H-bond
            # and must therefore be binary.  For non-specific solute-solvent
            # (UV) series, CPPTRAJ collapses solvent identity to V; multiple
            # solvent molecules can interact with the same solute endpoint in
            # one frame, so values may be integer counts > 1.
            if kind == "UU" and n not in {0, 1}:
                raise RuntimeError(
                    f"{path}: non-binary UU H-bond series value {x}"
                )

            vals.append(n)

        out.append((header[ci].strip('"'), tuple(vals)))
    return out


def series_to_bonds(
    path: Path,
    nframes: int,
    kind: str,
    avg_rows: dict[tuple[str,str,str], dict],
) -> list[BondSeries]:
    raw = read_series_matrix(path, nframes, kind)
    bonds: list[BondSeries] = []
    avg_keys = set(avg_rows)

    if kind == "UU":
        for legend, vals in raw:
            trip = _human_uu_legend_to_triplet(legend)
            if trip is None:
                # Some CPPTRAJ writers may prefix dataset metadata.  Try to
                # find an avg-row key whose canonical human legend occurs.
                matches = []
                for acc, dh, don in avg_keys:
                    dref = parse_atom_ref(don); href = parse_atom_ref(dh)
                    if dref and href:
                        human = f"{acc}-{don}-{href.atom}"
                        if human in legend:
                            matches.append((acc,dh,don))
                if len(matches) == 1:
                    trip = matches[0]
                else:
                    raise RuntimeError(
                        f"Cannot parse CPPTRAJ UU series legend {legend!r} in {path}"
                    )
            bonds.append(BondSeries(trip, trip[0], trip[1], trip[2], vals, kind))
        return bonds

    # UV legend from CPPTRAJ CreateHBlegend is `<solute atom>-V`.
    # If the atom is H, solute is donor and solvent is acceptor.  Otherwise
    # solute is acceptor and solvent is donor.
    for legend, vals in raw:
        text = legend.strip().strip('"')
        if text.endswith("-V"):
            atom_text = text[:-2]
        elif "-V" in text:
            atom_text = text.split("-V",1)[0].split()[-1]
        else:
            raise RuntimeError(f"Cannot parse CPPTRAJ UV series legend {legend!r}")
        ref = parse_atom_ref(atom_text)
        if ref is None:
            raise RuntimeError(f"Cannot parse solute atom in UV legend {legend!r}")
        if ref.atom.startswith("H"):
            heavy = donor_heavy_from_hydrogen(ref)
            donor = heavy.raw if heavy else ""
            trip = ("V", ref.raw, donor)
        else:
            trip = (ref.raw, "V", "V")
        bonds.append(BondSeries(trip, trip[0], trip[1], trip[2], vals, kind))
    return bonds


# ---------------------------------------------------------------------------
# User observables / CPPTRAJ actions
# ---------------------------------------------------------------------------

def parse_user_observables(args: argparse.Namespace) -> list[UserObservable]:
    items: list[UserObservable] = []
    seen: set[str] = set()
    def add(o: UserObservable):
        if o.name in seen: raise ValueError(f"Duplicate observable name: {o.name}")
        seen.add(o.name); items.append(o)
    for x in args.interaction or []:
        add(UserObservable(x[0], "interaction", "solute_solute", x[1], x[2]))
    for x in args.directed or []:
        add(UserObservable(x[0], "directed", "solute_solute", x[1], x[2]))
    for x in args.directed_h or []:
        add(UserObservable(x[0], "directed_h", "solute_solute", x[1], x[3], x[2]))
    for x in args.scan or []:
        add(UserObservable(x[0], "scan", "solute_solute", x[1], x[1]))
    for x in args.focus or []:
        add(UserObservable(x[0], "focus", "solute_solute", x[1], x[2]))
    for x in args.solvent_site or []:
        desc = ",".join(y for y in (args.solvent_donor_mask,
                                    args.solvent_acceptor_mask) if y)
        add(UserObservable(x[0], "solvent_site", "solute_solvent", x[1], desc))
    if not items:
        raise ValueError("No H-bond observable requested")
    return items


def hbond_geometry_keywords(args: argparse.Namespace) -> str:
    p = [f"dist {format_float(args.distance)}", f"angle {format_float(args.angle)}"]
    if args.hbond_image: p.append("image")
    return " ".join(p)


def build_actions(
    o: UserObservable,
    args: argparse.Namespace,
    run_dir: Path,
    start_index: int,
    topology: Path,
    mask_cache: dict[str, tuple[TopologyAtom, ...]],
) -> tuple[list[ActionSpec], int]:
    """Build CPPTRAJ actions for one observable.

    v5 broad-mode invariant
    -----------------------
    `interaction`, `scan`, and `focus` never put many donor heavy atoms into
    one CPPTRAJ action.  Instead, donor atoms are expanded from the topology
    and each donor gets its own action against an explicit acceptor atom list.
    This makes broad observables the union of atom-level directed searches and
    prevents the N46/N59->D52 class of candidates from disappearing merely
    because they were embedded in a broad donor mask.

    `directed` and `directed_h` remain literal expert overrides.
    """
    actions: list[ActionSpec] = []
    idx = start_index
    geom = hbond_geometry_keywords(args)

    def make(role: str, kind: str, core: str, solvent: bool = False):
        nonlocal idx
        aid = f"HB{idx:04d}"
        idx += 1
        # Include action ID because enumerated modes can generate hundreds of
        # actions with the same conceptual role.
        prefix = run_dir / f"{safe_component(o.name)}.{safe_component(role)}.{aid}"
        countf = Path(str(prefix) + ".count.dat")
        avgf = Path(str(prefix) + ".avg.dat")
        seriesf = Path(str(prefix) + ".series.dat")
        if solvent:
            cmd = (f"hbond {aid} {core} out {qpath(countf.resolve())} "
                   f"solvout {qpath(avgf.resolve())} series "
                   f"uvseries {qpath(seriesf.resolve())} {geom}")
        else:
            cmd = (f"hbond {aid} {core} out {qpath(countf.resolve())} "
                   f"avgout {qpath(avgf.resolve())} series "
                   f"uuseries {qpath(seriesf.resolve())} {geom}")
        actions.append(ActionSpec(aid, o.name, role, kind, cmd,
                                  countf, avgf, seriesf))

    def enumerate_direction(label: str, donor_user_mask: str,
                            acceptor_user_mask: str):
        donors = enumerated_role_atoms(topology, donor_user_mask, "donor",
                                       args, mask_cache)
        acceptors = enumerated_role_atoms(topology, acceptor_user_mask,
                                          "acceptor", args, mask_cache)
        if not donors or not acceptors:
            return
        all_acc = {a.index for a in acceptors}
        for donor in donors:
            # Exclude the donor atom itself if it can also act as an acceptor
            # (Ser/Thr/Tyr/His/CpH sites). CPPTRAJ should not count self bonds,
            # but removing it makes the candidate definition explicit.
            acc_idx = sorted(all_acc - {donor.index})
            if not acc_idx:
                continue
            make(
                f"{label}_D{donor.index}",
                "UU",
                f"donormask {qmask('@'+str(donor.index))} "
                f"acceptormask {qmask(explicit_atom_index_mask(acc_idx))}",
            )

    if o.mode == "interaction":
        enumerate_direction("selection1_donor", o.selection1, o.selection2)
        enumerate_direction("selection2_donor", o.selection2, o.selection1)
    elif o.mode == "directed":
        make("directed", "UU",
             f"donormask {qmask(o.selection1)} "
             f"acceptormask {qmask(o.selection2)}")
    elif o.mode == "directed_h":
        make("directed_h", "UU",
             f"donormask {qmask(o.selection1)} "
             f"donorhmask {qmask(o.donor_h_mask)} "
             f"acceptormask {qmask(o.selection2)}")
    elif o.mode == "scan":
        enumerate_direction("scan", o.selection1, o.selection1)
    elif o.mode == "focus":
        # Union of every focus donor -> partner acceptor and every partner
        # donor -> focus acceptor.  Overlap/focus-focus triplets are merged and
        # deduplicated downstream by physical (acceptor, donor-H, donor) key.
        enumerate_direction("focus_donor", o.selection1, o.selection2)
        enumerate_direction("focus_acceptor", o.selection2, o.selection1)
    elif o.mode == "solvent_site":
        # Solute-solvent detection is already asymmetric in CPPTRAJ and the UV
        # series identifies the solute endpoint. Keep one UV action here; the
        # CpH filter is applied per solute atom afterward.  The v5 donor-level
        # enumeration fix concerns UU broad-mask incompleteness.
        solvent_terms = []
        if args.solvent_donor_mask:
            solvent_terms.append(f"solventdonor {qmask(args.solvent_donor_mask)}")
        if args.solvent_acceptor_mask:
            solvent_terms.append(f"solventacceptor {qmask(args.solvent_acceptor_mask)}")
        make("solute_solvent", "UV",
             f"donormask {qmask(candidate_donor_mask(o.selection1,args))} "
             f"acceptormask {qmask(candidate_acceptor_mask(o.selection1,args))} "
             + " ".join(solvent_terms), solvent=True)
    else:
        raise AssertionError(o.mode)

    if not actions:
        raise RuntimeError(
            f"Observable {o.name!r} ({o.mode}) produced no donor/acceptor "
            "candidate actions after topology mask expansion"
        )
    return actions, idx


# ---------------------------------------------------------------------------
# Chemical filtering of per-bond series
# ---------------------------------------------------------------------------

def bond_chemistry_status(
    bond: BondSeries,
    lam: Optional[LambdaFile],
    lambda_values: Optional[tuple[float,...]],
    lambda_match_status: str,
    args: argparse.Namespace,
) -> str:
    if args.no_cph_filter:
        return "valid"

    if bond.kind == "UU":
        acc = parse_atom_ref(bond.acceptor)
        don = parse_atom_ref(bond.donor)
        if acc is None or don is None:
            raise RuntimeError(f"Cannot parse UU bond atom labels: {bond.key}")
        a = role_valid(acc, "acceptor", lam, lambda_values, args)
        d = role_valid(don, "donor", lam, lambda_values, args)
    else:
        # Solute acceptor <- solvent donor
        if bond.acceptor != "V":
            acc = parse_atom_ref(bond.acceptor)
            a = role_valid(acc, "acceptor", lam, lambda_values, args)
            d = "valid"
        else:
            # Solute donor-H -> solvent acceptor.  Infer donor heavy atom for
            # titratable groups.  For fixed groups, topology-based CPPTRAJ
            # donor detection already ensures a bonded H; no CpH filtering is
            # needed unless its residue is titratable.
            href = parse_atom_ref(bond.donor_h)
            if href is None:
                raise RuntimeError(f"Cannot parse UV donor-H label: {bond.donor_h}")
            heavy = parse_atom_ref(bond.donor) if bond.donor else donor_heavy_from_hydrogen(href)
            if lam is not None and lambda_columns_for_residue(lam, href.resid):
                if heavy is None:
                    raise RuntimeError(
                        f"Cannot infer donor heavy atom for titratable solvent H-bond {href.raw}"
                    )
                d = role_valid(heavy, "donor", lam, lambda_values, args)
            else:
                d = "valid"
            a = "valid"

    statuses = {a,d}
    if "outside" in statuses:
        return "outside_lambda_window"
    if "mixed" in statuses:
        return "mixed_lambda"
    if "invalid" in statuses:
        return "invalid"
    return "valid"


def merge_bonds(bond_lists: list[list[BondSeries]]) -> dict[tuple[str,str,str], BondSeries]:
    merged: dict[tuple[str,str,str], BondSeries] = {}
    for bonds in bond_lists:
        for b in bonds:
            if b.key not in merged:
                merged[b.key] = b
            else:
                old = merged[b.key]
                if old.values != b.values:
                    raise RuntimeError(
                        f"Duplicate physical H-bond {b.key} has inconsistent series"
                    )
    return merged


# ---------------------------------------------------------------------------
# Per-trajectory execution
# ---------------------------------------------------------------------------

def job_payload(topology: Path, trajectory: Path, lambda_path: Optional[Path],
                observables: list[UserObservable], actions: list[ActionSpec],
                args: argparse.Namespace) -> dict:
    return {
        "settings_version": SETTINGS_VERSION,
        "program_version": PROGRAM_VERSION,
        "topology": str(topology.resolve()),
        "trajectory": str(trajectory.resolve()),
        "lambda": str(lambda_path.resolve()) if lambda_path else None,
        "observables": [o.__dict__ for o in observables],
        "actions": [{"id":a.action_id,"parent":a.parent_name,"role":a.role,
                     "kind":a.kind,"command":a.command} for a in actions],
        "distance": args.distance, "angle": args.angle,
        "image_anchor": args.image_anchor, "hbond_image": args.hbond_image,
        "start": args.start, "stop": args.stop, "stride": args.stride,
        "first_frame_step": args.first_frame_step,
        "frame_step_interval": args.frame_step_interval,
        "lambda_low": args.lambda_low, "lambda_high": args.lambda_high,
        "chemistry_profile": args.chemistry_profile,
        "candidate_engine": "donor-enumerated-v5",
        "include_met_sulfur": args.include_met_sulfur,
        "no_cph_filter": args.no_cph_filter,
        "solvent_donor_mask": args.solvent_donor_mask,
        "solvent_acceptor_mask": args.solvent_acceptor_mask,
        "cpptraj": args.cpptraj,
    }


def run_one_trajectory(*, topology: Path, trajectory: Path,
                       lambda_path: Optional[Path],
                       observables: list[UserObservable],
                       args: argparse.Namespace) -> dict:
    stem = trajectory.stem
    run_dir = args.workdir / stem
    run_dir.mkdir(parents=True, exist_ok=True)
    actions: list[ActionSpec] = []
    nxt = 1
    mask_cache: dict[str, tuple[TopologyAtom, ...]] = {}
    for o in observables:
        aa, nxt = build_actions(o, args, run_dir, nxt, topology, mask_cache)
        actions.extend(aa)

    inputf = run_dir/"hbond.in"; logf=run_dir/"hbond.log"; settingsf=run_dir/"job.settings.json"
    payload = job_payload(topology,trajectory,lambda_path,observables,actions,args)
    payload["fingerprint"] = stable_hash(payload)
    expected = [p for a in actions for p in (a.count_file,a.avg_file)]

    reuse = False
    if not args.force and settingsf.is_file():
        prev=json.loads(settingsf.read_text())
        if prev.get("fingerprint") == payload["fingerprint"] and all(p.is_file() for p in expected):
            reuse=True
        elif prev.get("fingerprint") != payload["fingerprint"]:
            raise RuntimeError(f"Existing work settings differ for {stem}; use --force or new --workdir")
    if args.force:
        for a in actions:
            for p in (a.count_file,a.avg_file,a.series_file):
                if p.exists(): p.unlink()

    if not reuse:
        lines=[f"parm {qpath(topology.resolve())}",
               f"trajin {qpath(trajectory.resolve())} {args.start} {args.stop} {args.stride}"]
        if args.image_anchor: lines.append(f"autoimage anchor {qmask(args.image_anchor)}")
        lines.extend(a.command for a in actions); lines += ["run","quit",""]
        inputf.write_text("\n".join(lines))
        settingsf.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")
        with logf.open("w") as fh:
            proc=subprocess.run([args.cpptraj,"-i",str(inputf)],stdout=fh,stderr=subprocess.STDOUT)
        if proc.returncode != 0:
            raise RuntimeError(f"CPPTRAJ failed for {trajectory.name}; see {logf}")
        missing=[p for p in expected if not p.is_file()]
        if missing:
            raise RuntimeError(f"CPPTRAJ missing expected files for {trajectory.name}: {missing}")

    count_by_action: dict[str,list[int]]={}
    avg_by_action: dict[str,dict]={}
    for a in actions:
        count_by_action[a.action_id]=read_hbond_count_series(a.count_file,a.kind)
        avg_by_action[a.action_id]=read_hbond_average_report(a.avg_file)
    lengths={len(v) for v in count_by_action.values()}
    if len(lengths)!=1: raise RuntimeError(f"{stem}: inconsistent action frame counts {sorted(lengths)}")
    nframes=next(iter(lengths))
    if nframes==0: raise RuntimeError(f"{stem}: no analyzed frames")

    # Series files may not be emitted if no bond was ever detected.  If the
    # aggregate count ever exceeds zero, however, a missing/empty series file
    # is a hard error because chemistry-aware reconstruction would be impossible.
    bonds_by_action: dict[str,list[BondSeries]]={}
    for a in actions:
        if a.series_file.is_file() and a.series_file.stat().st_size>0:
            bonds_by_action[a.action_id]=series_to_bonds(a.series_file,nframes,a.kind,avg_by_action[a.action_id])
        else:
            if max(count_by_action[a.action_id])>0:
                raise RuntimeError(
                    f"{stem}/{a.action_id}: geometric H-bonds were detected but per-bond series file "
                    f"{a.series_file} is missing/empty"
                )
            bonds_by_action[a.action_id]=[]

    lam = read_lambda_file(lambda_path) if lambda_path else None

    # Frame metadata follows the original stored trajectory indexing.
    trajectory_frames=[args.start + i*args.stride for i in range(nframes)]
    steps: list[Optional[int]]=[]
    for tf in trajectory_frames:
        if args.first_frame_step is None:
            steps.append(None)
        else:
            steps.append(args.first_frame_step + (tf-1)*args.frame_step_interval)

    if not args.no_cph_filter and lam is not None and any(s is None for s in steps):
        raise RuntimeError("CpH filtering requires --first-frame-step and --frame-step-interval")

    actions_by_parent: dict[str,list[ActionSpec]]={}
    for a in actions: actions_by_parent.setdefault(a.parent_name,[]).append(a)

    combined: dict[str,list[int]]={}
    ambiguous_series: dict[str,list[int]]={}
    outside_series: dict[str,list[int]]={}
    statuses: dict[str,list[str]]={}
    detail_rows: list[dict]=[]
    ph=extract_ph(stem)

    for o in observables:
        parent=actions_by_parent[o.name]
        merged=merge_bonds([bonds_by_action[a.action_id] for a in parent])
        counts=[0]*nframes
        ambiguous_counts=[0]*nframes
        outside_counts=[0]*nframes
        frame_status=["clean"]*nframes

        # per-bond filtered statistics
        bond_valid_counts={k:0 for k in merged}
        bond_ambig_counts={k:0 for k in merged}
        bond_outside_counts={k:0 for k in merged}

        for fi in range(nframes):
            step=steps[fi]
            if args.no_cph_filter:
                lstatus="matched"; vals=None
            elif lam is None:
                raise RuntimeError("Internal error: CpH filtering requested without lambda")
            else:
                assert step is not None
                lstatus, vals = lambda_row_exact(lam,step)

            ambiguous_present=False
            outside_present=False
            for key, b in merged.items():
                multiplicity = b.values[fi]
                if multiplicity == 0:
                    continue

                cstat = bond_chemistry_status(b, lam, vals, lstatus, args)

                if cstat == "valid":
                    counts[fi] += multiplicity
                    bond_valid_counts[key] += multiplicity

                elif cstat == "invalid":
                    pass

                elif cstat == "mixed_lambda":
                    ambiguous_present = True
                    ambiguous_counts[fi] += multiplicity
                    bond_ambig_counts[key] += multiplicity

                elif cstat == "outside_lambda_window":
                    outside_present = True
                    outside_counts[fi] += multiplicity
                    bond_outside_counts[key] += multiplicity

                else:
                    raise AssertionError(cstat)
            if outside_present:
                frame_status[fi]="outside_lambda_window"
            elif ambiguous_present:
                frame_status[fi]="partial_ambiguity"

        combined[o.name]=counts
        ambiguous_series[o.name]=ambiguous_counts
        outside_series[o.name]=outside_counts
        statuses[o.name]=frame_status

        # Merge raw avg rows across the parent actions by physical triplet.
        raw_avg: dict[tuple[str,str,str],dict]={}
        for a in parent:
            for k,row in avg_by_action[a.action_id].items():
                if k not in raw_avg: raw_avg[k]=row
        for key,b in sorted(merged.items(), key=lambda kv:(-bond_valid_counts[kv[0]],kv[0])):
            valid=bond_valid_counts[key]
            raw_count=sum(b.values)
            raw=raw_avg.get(key,{})
            # For UV, key may not match solvout formatting.  Keep geometry blank
            # rather than inventing an association.

            if b.kind == "UU":
                count_semantics = "chemically_valid_frames_present"
            else:
                count_semantics = "chemically_valid_interactions"

            detail_rows.append({
                "trajectory": stem,
                "pH": ph,
                "observable": o.name,
                "mode": o.mode,
                "interaction_class": o.interaction_class,
                "selection1": o.selection1,
                "selection2": o.selection2,
                "acceptor": b.acceptor,
                "donor_h": b.donor_h,
                "donor": b.donor,
                "count": valid,
                "count_semantics": count_semantics,
                "fraction": valid / nframes,
                "raw_count": raw_count,
                "raw_fraction": raw_count / nframes,
                "ambiguous_present_frames": bond_ambig_counts[key],
                "outside_lambda_present_frames": bond_outside_counts[key],
                "raw_avg_distance_A": raw.get("raw_avg_distance_A", math.nan),
                "raw_avg_angle_deg": raw.get("raw_avg_angle_deg", math.nan),
                "geometry_average_scope":
                    "CPPTRAJ geometric-candidate-present frames before CpH filtering",
            })
    return {"trajectory":trajectory,"nframes":nframes,"trajectory_frames":trajectory_frames,
            "steps":steps,"combined_series":combined,
            "ambiguous_series":ambiguous_series,"outside_series":outside_series,
            "chemistry_status":statuses,"detail_rows":detail_rows,
            "lambda_file":lambda_path,"reused":reuse}


# ---------------------------------------------------------------------------
# CLI / main
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    ap=argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
        description="Lambda-aware CpHMD hydrogen-bond analysis using CPPTRAJ geometric series.")
    ap.add_argument("-p","--topology",type=Path,required=True)
    ap.add_argument("-t","--trajectory",action="append",required=True,
                    help="Trajectory path or quoted glob; may be repeated")
    ap.add_argument("--lambda",dest="lambda_specs",action="append",
                    help="Matching CpHMD lambda path or quoted glob; may be repeated")
    ap.add_argument("--no-cph-filter",action="store_true",
                    help="Fixed-protonation control: do not read/apply lambda chemistry")
    ap.add_argument("--lambda-low",type=float,default=0.2)
    ap.add_argument("--lambda-high",type=float,default=0.8)

    ap.add_argument("--interaction",action="append",nargs=3,metavar=("NAME","MASK1","MASK2"))
    ap.add_argument("--directed",action="append",nargs=3,metavar=("NAME","DONOR_MASK","ACCEPTOR_MASK"))
    ap.add_argument("--directed-h",action="append",nargs=4,
                    metavar=("NAME","DONOR_MASK","DONOR_H_MASK","ACCEPTOR_MASK"))
    ap.add_argument("--scan",action="append",nargs=2,metavar=("NAME","MASK"))
    ap.add_argument("--focus",action="append",nargs=3,metavar=("NAME","FOCUS_MASK","PARTNER_MASK"))
    ap.add_argument("--solvent-site",action="append",nargs=2,metavar=("NAME","SOLUTE_MASK"))
    ap.add_argument("--solvent-donor-mask",default=":WAT")
    ap.add_argument("--solvent-acceptor-mask",default=":WAT@O")

    ap.add_argument("--chemistry-profile",choices=("protein","fon"),default="protein",
                    help="Candidate-generation profile. CpH state filtering is applied afterward.")
    ap.add_argument("--include-met-sulfur",action="store_true",
                    help="Include fixed Met SD as an acceptor candidate")
    ap.add_argument("--distance",type=float,default=3.0)
    ap.add_argument("--angle",type=float,default=135.0)
    ap.add_argument("--image-anchor",default=None)
    ap.add_argument("--hbond-image",action="store_true")
    ap.add_argument("--start",type=int,default=1)
    ap.add_argument("--stop",default="last")
    ap.add_argument("--stride",type=int,default=1)
    ap.add_argument("--first-frame-step",type=int,default=None)
    ap.add_argument("--frame-step-interval",type=int,default=None)
    ap.add_argument("-j","--jobs","--fork",dest="jobs",type=int,default=1)
    ap.add_argument("--cpptraj",default="cpptraj")
    ap.add_argument("--workdir",type=Path,default=Path("hbond_work_v5"))
    ap.add_argument("--detail-output",type=Path,default=None)
    ap.add_argument("-o","--output",type=Path,default=Path("hbonds_v5.tsv"))
    ap.add_argument("--force",action="store_true")
    return ap.parse_args()


def validate_args(args: argparse.Namespace) -> None:
    if not args.topology.is_file(): raise FileNotFoundError(args.topology)
    if args.start<1: raise ValueError("--start must be >=1")
    args.stop=parse_stop(args.stop)
    if args.stop!="last" and int(args.stop)<args.start: raise ValueError("--stop must be >= --start")
    if args.stride<1: raise ValueError("--stride must be >=1")
    if args.jobs<1: raise ValueError("--jobs must be >=1")
    if not (math.isfinite(args.distance) and args.distance>0): raise ValueError("--distance must be >0")
    if not (math.isfinite(args.angle) and (args.angle==-1 or 0<=args.angle<=180)):
        raise ValueError("--angle must be -1 or within 0..180")
    if not (0<=args.lambda_low<args.lambda_high<=1):
        raise ValueError("lambda thresholds must satisfy 0 <= low < high <= 1")
    supplied=(args.first_frame_step is not None,args.frame_step_interval is not None)
    if supplied[0]!=supplied[1]:
        raise ValueError("--first-frame-step and --frame-step-interval must be supplied together")
    if args.frame_step_interval is not None and args.frame_step_interval<=0:
        raise ValueError("--frame-step-interval must be >0")
    if not args.no_cph_filter:
        if not args.lambda_specs:
            raise ValueError("CpHMD mode requires --lambda (or use --no-cph-filter for a fixed-state control)")
        if args.first_frame_step is None:
            raise ValueError("CpHMD mode requires exact step mapping via --first-frame-step and --frame-step-interval")
    args.solvent_donor_mask=args.solvent_donor_mask.strip()
    args.solvent_acceptor_mask=args.solvent_acceptor_mask.strip()
    if not (args.solvent_donor_mask or args.solvent_acceptor_mask) and args.solvent_site:
        raise ValueError("--solvent-site requires at least one solvent donor/acceptor mask")
    exe=shutil.which(args.cpptraj) if os.sep not in args.cpptraj else (args.cpptraj if Path(args.cpptraj).is_file() else None)
    if exe is None:
        raise FileNotFoundError(
            f"CPPTRAJ executable not found: {args.cpptraj!r}. Source Amber or pass --cpptraj /full/path/cpptraj"
        )
    args.cpptraj=str(exe)


def settings_text(args: argparse.Namespace, trajectories: list[Path], lambda_map: dict[Path,Path],
                  observables: list[UserObservable]) -> str:
    lines=[
        f"Program version             : {PROGRAM_VERSION}",
        "Analysis                    : lambda-aware CpHMD H-bonds",
        f"Chemistry profile           : {args.chemistry_profile}",
        "Candidate engine            : donor-enumerated-v5.1",
        f"CpH filtering               : {'disabled' if args.no_cph_filter else 'enabled'}",
        f"Lambda cutoffs              : <= {args.lambda_low:g} protonated; >= {args.lambda_high:g} deprotonated",
        f"Distance cutoff (A)         : {args.distance:g}",
        f"Angle cutoff (deg)          : {args.angle:g}",
        f"Trajectories                : {len(trajectories)}",
        f"Requested observables       : {len(observables)}",
    ]
    for t in trajectories:
        lines.append(f"trajectory {t.name} lambda {lambda_map.get(t,'NONE')}")
    return "\n".join(lines)+"\n"


def main() -> int:
    args=parse_args()
    try:
        validate_args(args)
        trajectories=expand_specs(args.trajectory,"trajectory")
        lambda_files=[] if args.no_cph_filter else expand_specs(args.lambda_specs,"lambda")
        lambda_map=match_lambda_files(trajectories,lambda_files)
        observables=parse_user_observables(args)
        args.workdir=args.workdir.resolve(); args.workdir.mkdir(parents=True,exist_ok=True)

        results: dict[Path,dict]={}
        def work(t:Path):
            print(f"[{t.stem}] running",flush=True)
            r=run_one_trajectory(topology=args.topology.resolve(),trajectory=t,
                                 lambda_path=lambda_map.get(t),observables=observables,args=args)
            print(f"[{t.stem}] complete",flush=True)
            return r
        if args.jobs==1:
            for t in trajectories: results[t]=work(t)
        else:
            with ThreadPoolExecutor(max_workers=min(args.jobs,len(trajectories))) as ex:
                fut={ex.submit(work,t):t for t in trajectories}
                for f in as_completed(fut): results[fut[f]]=f.result()

        frame_rows: list[dict]=[]; detail_rows: list[dict]=[]
        for t in trajectories:
            r=results[t]; ph=extract_ph(t.stem); n=r["nframes"]
            for o in observables:
                vals=r["combined_series"][o.name]
                amb=r["ambiguous_series"][o.name]
                outside=r["outside_series"][o.name]
                stats=r["chemistry_status"][o.name]
                for i in range(n):
                    base={"trajectory":t.stem,"pH":ph,"analysis_row":i+1,
                          "trajectory_frame":r["trajectory_frames"][i],
                          "step":r["steps"][i],"observable":o.name,
                          "interaction_class":o.interaction_class,"mode":o.mode,
                          "selection1":o.selection1,"selection2":o.selection2,
                          "chemistry_status":stats[i],
                          "lambda_file":Path(r["lambda_file"]).name if r["lambda_file"] else ""}
                    v=vals[i]
                    unresolved=amb[i]+outside[i]
                    frame_rows.append({**base,"metric":"hbond_count","value":v,"unit":"count"})
                    frame_rows.append({**base,"metric":"hbond_present","value":int(v>0),"unit":"1"})
                    frame_rows.append({**base,"metric":"ambiguous_hbond_count","value":amb[i],"unit":"count"})
                    frame_rows.append({**base,"metric":"outside_lambda_hbond_count","value":outside[i],"unit":"count"})
                    frame_rows.append({**base,"metric":"hbond_count_min","value":v,"unit":"count"})
                    frame_rows.append({**base,"metric":"hbond_count_max","value":v+unresolved,"unit":"count"})
            detail_rows.extend(r["detail_rows"])

        # Deterministic order irrespective of jobs.
        order={t.stem:i for i,t in enumerate(trajectories)}
        obsorder={o.name:i for i,o in enumerate(observables)}
        metricorder={"hbond_count":0,"hbond_present":1,"ambiguous_hbond_count":2,
                     "outside_lambda_hbond_count":3,"hbond_count_min":4,"hbond_count_max":5}
        frame_rows.sort(key=lambda x:(order[x["trajectory"]],obsorder[x["observable"]],
                                      x["analysis_row"],metricorder[x["metric"]]))
        detail_rows.sort(key=lambda x:(order[x["trajectory"]],obsorder[x["observable"]],
                                       -int(x["count"]),x["acceptor"],x["donor"],x["donor_h"]))

        args.output=args.output.resolve()
        if args.detail_output is None:
            args.detail_output=args.output.with_name(args.output.stem+".interactions.tsv")
        else: args.detail_output=args.detail_output.resolve()
        write_tsv(args.output,frame_rows,OUTPUT_FIELDS)
        write_tsv(args.detail_output,detail_rows,DETAIL_FIELDS)
        settings_path=args.output.with_name(args.output.name+".settings.txt")
        settings_path.write_text(settings_text(args,trajectories,lambda_map,observables))

        n_partial=sum(1 for row in frame_rows
                      if row["metric"]=="hbond_count" and row["chemistry_status"]=="partial_ambiguity")
        n_outside=sum(1 for row in frame_rows
                      if row["metric"]=="hbond_count" and row["chemistry_status"]=="outside_lambda_window")
        print()
        print(f"Trajectories analyzed       : {len(trajectories)}")
        print(f"Requested observables       : {len(observables)}")
        print(f"Frame-wise output rows      : {len(frame_rows)}")
        print(f"Atom-level report rows      : {len(detail_rows)}")
        print(f"Program version             : {PROGRAM_VERSION}")
        print(f"Chemistry profile           : {args.chemistry_profile}")
        print("Candidate engine            : donor-enumerated-v5.1")
        print(f"CpH filtering               : {'disabled' if args.no_cph_filter else 'enabled'}")
        print(f"Partial-ambiguity frame obs : {n_partial}")
        print(f"Outside-lambda frame obs    : {n_outside}")
        print(f"Distance cutoff (A)         : {args.distance:g}")
        print(f"Angle cutoff (deg)          : {args.angle:g}")
        print(f"Wrote {args.output.name}")
        print(f"Wrote {args.detail_output.name}")
        print(f"Wrote {settings_path.name}")
        return 0
    except Exception as exc:
        print(f"ERROR: {exc}",file=sys.stderr)
        return 1


if __name__=="__main__":
    sys.exit(main())
