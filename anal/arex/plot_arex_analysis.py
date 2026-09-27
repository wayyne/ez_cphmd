#!/usr/bin/env python3
"""
Summarize and plot residue-centric CPPTRAJ output from analyze_arex.sh.

The script:
  1. reads each pH-slot trajectory directory independently;
  2. computes mean, sample SD, and a block-based SEM estimate;
  3. writes tidy TSV summary tables;
  4. optionally joins exact MD-step dyad protonation-state assignments and writes
     state-conditioned structural summaries;
  5. creates publication-oriented PDF/SVG/PNG figures with Matplotlib.

No statistical independence is implied by the block SEM. Increase --blocks only
if each block remains comfortably longer than the correlation time of the
observable. For final manuscript uncertainty estimates, convergence/block-size
sensitivity should be checked explicitly.
"""

from __future__ import annotations

import argparse
import csv
import math
import re
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--analysis-dir", type=Path, required=True)
    p.add_argument(
        "--blocks",
        type=int,
        default=5,
        help="Number of contiguous blocks used for block-SEM estimation (default: 5).",
    )
    p.add_argument(
        "--max-ph-lines",
        type=int,
        default=15,
        help=(
            "If <= this many residues are analyzed, plot property vs pH with one "
            "line per residue. Above this threshold, plot property vs residue "
            "number with one line per pH/trajectory."
        ),
    )
    p.add_argument(
        "--state-file",
        type=Path,
        default=None,
        help=(
            "Optional frame-resolved coupled-state CSV written by "
            "calc_titration.py. When supplied, state-conditioned "
            "dyad tables are generated in addition to the ordinary summaries."
        ),
    )
    p.add_argument(
        "--first-frame-step",
        type=int,
        default=None,
        help=(
            "MD step corresponding to frame 1 of each input coordinate trajectory. "
            "Required with --state-file."
        ),
    )
    p.add_argument(
        "--frame-step-interval",
        type=int,
        default=None,
        help=(
            "Number of MD steps between consecutive stored input coordinate frames. "
            "Required with --state-file."
        ),
    )
    return p.parse_args()


def read_tsv(path: Path) -> List[dict]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_csv(path: Path) -> List[dict]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def read_numeric(path: Path) -> np.ndarray:
    if not path.exists() or path.stat().st_size == 0:
        return np.empty((0, 0), dtype=float)

    arr = np.genfromtxt(path, comments="#", dtype=float)
    if arr.size == 0:
        return np.empty((0, 0), dtype=float)
    if arr.ndim == 1:
        arr = arr.reshape(1, -1)
    # Drop rows that are completely non-finite, if any.
    arr = arr[np.any(np.isfinite(arr), axis=1)]
    return arr

def determine_processed_frame_count(run_dir: Path) -> int:
    """
    Determine the number of processed coordinate frames from available
    per-frame CPPTRAJ outputs.

    The dyad state join does not require RMSD specifically. Any ordinary
    per-frame structural output can establish the processed frame count.
    When multiple outputs are available, require them to agree.
    """
    candidate_paths = [
        run_dir / "rmsd.dat",
        run_dir / "rg.dat",
        run_dir / "sasa_residue.dat",
        run_dir / "sasa_sidechain.dat",
        run_dir / "hbond_water.dat",
        run_dir / "hbond_protein_donor.dat",
        run_dir / "hbond_protein_acceptor.dat",
    ]

    candidate_paths.extend(sorted(run_dir.glob("watershell_R*.dat")))

    counts: Dict[str, int] = {}

    for path in candidate_paths:
        arr = read_numeric(path)
        if arr.shape[0] > 0:
            counts[path.name] = int(arr.shape[0])

    if not counts:
        raise RuntimeError(
            f"Cannot determine processed frame count for {run_dir.name}: "
            f"no readable per-frame CPPTRAJ outputs were found in {run_dir}"
        )

    unique_counts = sorted(set(counts.values()))

    if len(unique_counts) != 1:
        details = ", ".join(
            f"{name}={count}"
            for name, count in sorted(counts.items())
        )
        raise RuntimeError(
            f"Inconsistent processed frame counts for {run_dir.name}: {details}"
        )

    return unique_counts[0]

def summarize(values: np.ndarray, nblocks: int) -> dict:
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    n = int(values.size)

    if n == 0:
        return {
            "nframes": 0,
            "mean": np.nan,
            "sd": np.nan,
            "block_sem": np.nan,
            "nblocks": 0,
        }

    mean = float(np.mean(values))
    sd = float(np.std(values, ddof=1)) if n > 1 else np.nan

    nb = min(max(1, nblocks), n)
    chunks = [x for x in np.array_split(values, nb) if x.size]
    block_means = np.array([np.mean(x) for x in chunks], dtype=float)

    if block_means.size > 1:
        block_sem = float(np.std(block_means, ddof=1) / np.sqrt(block_means.size))
    else:
        block_sem = np.nan

    return {
        "nframes": n,
        "mean": mean,
        "sd": sd,
        "block_sem": block_sem,
        "nblocks": int(block_means.size),
    }


def summarize_state_masked(
    values: np.ndarray,
    state_labels: np.ndarray,
    target_state: str,
    nblocks: int,
) -> dict:
    """Summarize one observable within one protonation state.

    The mean and SD are calculated over frames assigned to target_state. For the
    block-SEM diagnostic, the original trajectory timeline is divided into
    contiguous blocks first, and the target-state mean is calculated within each
    block. This preserves the temporal spacing of intermittent state visits instead
    of concatenating all target-state frames before blocking.
    """
    values = np.asarray(values, dtype=float)
    state_labels = np.asarray(state_labels, dtype=object)

    if values.ndim != 1:
        values = values.reshape(-1)
    if state_labels.ndim != 1:
        state_labels = state_labels.reshape(-1)
    if values.size != state_labels.size:
        raise RuntimeError(
            f"State/observable length mismatch: {state_labels.size} states for "
            f"{values.size} observable values"
        )

    mask = (state_labels == target_state) & np.isfinite(values)
    selected = values[mask]
    n = int(selected.size)

    if n == 0:
        return {
            "nframes": 0,
            "mean": np.nan,
            "sd": np.nan,
            "block_sem": np.nan,
            "nblocks": 0,
        }

    mean_value = float(np.mean(selected))
    sd = float(np.std(selected, ddof=1)) if n > 1 else np.nan

    nb = min(max(1, nblocks), int(values.size))
    index_blocks = [x for x in np.array_split(np.arange(values.size), nb) if x.size]
    block_means: List[float] = []
    for indices in index_blocks:
        block_mask = mask[indices]
        if not np.any(block_mask):
            continue
        block_values = values[indices][block_mask]
        if block_values.size:
            block_means.append(float(np.mean(block_values)))

    if len(block_means) > 1:
        block_sem = float(
            np.std(np.asarray(block_means, dtype=float), ddof=1)
            / np.sqrt(len(block_means))
        )
    else:
        block_sem = np.nan

    return {
        "nframes": n,
        "mean": mean_value,
        "sd": sd,
        "block_sem": block_sem,
        "nblocks": len(block_means),
    }


def configure_matplotlib() -> None:
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
            "font.size": 8.5,
            "axes.labelsize": 9.5,
            "axes.linewidth": 0.8,
            "xtick.labelsize": 8.0,
            "ytick.labelsize": 8.0,
            "legend.fontsize": 7.5,
            "lines.linewidth": 1.25,
            "lines.markersize": 4.5,
            "xtick.direction": "out",
            "ytick.direction": "out",
            "xtick.major.width": 0.8,
            "ytick.major.width": 0.8,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
        }
    )


def clean_axis(ax: plt.Axes) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(False)


def export_figure(fig: plt.Figure, stem: Path) -> None:
    stem.parent.mkdir(parents=True, exist_ok=True)
    for ext, kwargs in {
        "pdf": {},
        "svg": {},
        "png": {"dpi": 600},
    }.items():
        fig.savefig(
            stem.with_suffix(f".{ext}"),
            bbox_inches="tight",
            pad_inches=0.03,
            **kwargs,
        )
    plt.close(fig)


def trajectory_order(metadata: List[dict]) -> List[str]:
    def key(row: dict):
        ph = row.get("ph", "").strip()
        if ph:
            try:
                return (0, float(ph), row["label"])
            except ValueError:
                pass
        return (1, math.inf, row["label"])

    return [row["label"] for row in sorted(metadata, key=key)]


def _parse_integral_step(value: str, context: str) -> int:
    try:
        numeric = float(value)
    except (TypeError, ValueError) as error:
        raise RuntimeError(f"Invalid MD step in {context}: {value!r}") from error
    if not math.isfinite(numeric):
        raise RuntimeError(f"Non-finite MD step in {context}: {value!r}")
    rounded = int(round(numeric))
    if abs(numeric - rounded) > 1.0e-6:
        raise RuntimeError(f"Non-integral MD step in {context}: {value!r}")
    return rounded


def load_state_file(
    path: Path,
) -> Tuple[
    Dict[Tuple[str, int], dict],
    Dict[str, Tuple[int, int]],
    Dict[str, float],
]:
    """Read the frame-resolved coupled-state CSV.

    Returns an exact (slot_label, MD step) lookup, the inclusive step range
    represented for each slot, and the pH associated with each slot.
    """
    rows = read_csv(path)
    required = {
        "slot_label",
        "pH",
        "step",
        "absolute_time_ns",
        "primary_label",
        "primary_lambda",
        "coupled_label",
        "coupled_lambda",
        "state_code",
        "state_label",
        "clean",
    }
    if not rows:
        raise RuntimeError(f"State file is empty: {path}")
    missing = required.difference(rows[0])
    if missing:
        raise RuntimeError(
            f"State file {path} is missing required columns: "
            + ", ".join(sorted(missing))
        )

    lookup: Dict[Tuple[str, int], dict] = {}
    steps_by_slot: Dict[str, List[int]] = {}
    ph_by_slot: Dict[str, float] = {}

    for row_number, row in enumerate(rows, start=2):
        slot = row["slot_label"].strip()
        if not slot:
            raise RuntimeError(f"Blank slot_label in {path} line {row_number}")
        step = _parse_integral_step(row["step"], f"{path} line {row_number}")
        try:
            ph = float(row["pH"])
            absolute_time_ns = float(row["absolute_time_ns"])
            state_code = int(row["state_code"])
            primary_lambda = float(row["primary_lambda"])
            coupled_lambda = float(row["coupled_lambda"])
        except ValueError as error:
            raise RuntimeError(
                f"Malformed numeric value in {path} line {row_number}"
            ) from error

        clean_text = row["clean"].strip().lower()
        if clean_text in {"1", "true", "yes"}:
            clean = True
        elif clean_text in {"0", "false", "no"}:
            clean = False
        else:
            raise RuntimeError(
                f"Invalid clean flag in {path} line {row_number}: {row['clean']!r}"
            )

        key = (slot, step)
        if key in lookup:
            raise RuntimeError(
                f"Duplicate coupled-state assignment for slot {slot!r}, step {step} "
                f"in {path}"
            )

        normalized = {
            "slot_label": slot,
            "pH": ph,
            "step": step,
            "absolute_time_ns": absolute_time_ns,
            "primary_label": row["primary_label"].strip(),
            "primary_lambda": primary_lambda,
            "coupled_label": row["coupled_label"].strip(),
            "coupled_lambda": coupled_lambda,
            "state_code": state_code,
            "state_label": row["state_label"].strip(),
            "clean": clean,
        }
        lookup[key] = normalized
        steps_by_slot.setdefault(slot, []).append(step)

        if slot in ph_by_slot and abs(ph_by_slot[slot] - ph) > 1.0e-8:
            raise RuntimeError(
                f"State file assigns multiple pH values to slot {slot!r}: "
                f"{ph_by_slot[slot]} and {ph}"
            )
        ph_by_slot[slot] = ph

    ranges = {
        slot: (min(steps), max(steps))
        for slot, steps in steps_by_slot.items()
        if steps
    }
    return lookup, ranges, ph_by_slot


def build_frame_state_assignments(
    analysis_dir: Path,
    metadata: List[dict],
    state_lookup: Dict[Tuple[str, int], dict],
    state_ranges: Dict[str, Tuple[int, int]],
    state_ph: Dict[str, float],
    first_frame_step: int,
    frame_step_interval: int,
) -> Tuple[Dict[str, List[dict]], List[dict], int, int]:
    """Map each processed CPPTRAJ row to an exact MD-step state assignment.

    Frames before or after the state CSV's step range are labeled
    outside_state_window and excluded from state-conditioned summaries. Any frame
    whose MD step lies inside the represented range but lacks an exact state row is
    treated as an error; no nearest-neighbor matching is performed.
    """
    by_trajectory: Dict[str, List[dict]] = {}
    flat_rows: List[dict] = []
    missing_inside: List[Tuple[str, int, int]] = []
    matched_count = 0
    outside_count = 0

    for meta in metadata:
        label = meta["label"]
        if label not in state_ranges:
            raise RuntimeError(
                f"No state-file rows were found for trajectory/slot {label!r}. "
                "The NetCDF basename and state-file slot_label must match exactly."
            )

        ph_text = meta.get("ph", "").strip()
        ph = float(ph_text) if ph_text else np.nan
        if np.isfinite(ph) and label in state_ph:
            if abs(ph - state_ph[label]) > 1.0e-8:
                raise RuntimeError(
                    f"pH mismatch for {label}: metadata.tsv gives {ph:g}, "
                    f"state file gives {state_ph[label]:g}"
                )

        run_dir = analysis_dir / label
        nrows = determine_processed_frame_count(run_dir)
        try:
            start = int(meta.get("start", "1") or 1)
            stride = int(meta.get("stride", "1") or 1)
        except ValueError as error:
            raise RuntimeError(
                f"Invalid start/stride values in metadata.tsv for {label}"
            ) from error
        if start < 1 or stride < 1:
            raise RuntimeError(
                f"Invalid start/stride values for {label}: start={start}, stride={stride}"
            )

        min_step, max_step = state_ranges[label]
        trajectory_rows: List[dict] = []

        for analysis_index in range(nrows):
            trajectory_frame = start + analysis_index * stride
            step = first_frame_step + (trajectory_frame - 1) * frame_step_interval
            key = (label, step)

            if step < min_step or step > max_step:
                row = {
                    "trajectory": label,
                    "ph": ph,
                    "analysis_row": analysis_index + 1,
                    "trajectory_frame": trajectory_frame,
                    "step": step,
                    "absolute_time_ns": np.nan,
                    "primary_label": "",
                    "primary_lambda": np.nan,
                    "coupled_label": "",
                    "coupled_lambda": np.nan,
                    "state_code": "",
                    "state_label": "outside_state_window",
                    "clean": 0,
                    "matched": 0,
                }
                outside_count += 1
            else:
                state = state_lookup.get(key)
                if state is None:
                    missing_inside.append((label, trajectory_frame, step))
                    continue
                row = {
                    "trajectory": label,
                    "ph": ph,
                    "analysis_row": analysis_index + 1,
                    "trajectory_frame": trajectory_frame,
                    "step": step,
                    "absolute_time_ns": state["absolute_time_ns"],
                    "primary_label": state["primary_label"],
                    "primary_lambda": state["primary_lambda"],
                    "coupled_label": state["coupled_label"],
                    "coupled_lambda": state["coupled_lambda"],
                    "state_code": state["state_code"],
                    "state_label": state["state_label"],
                    "clean": int(state["clean"]),
                    "matched": 1,
                }
                matched_count += 1

            trajectory_rows.append(row)
            flat_rows.append(row)

        by_trajectory[label] = trajectory_rows

    if missing_inside:
        examples = "; ".join(
            f"{label} frame {frame} -> step {step}"
            for label, frame, step in missing_inside[:10]
        )
        raise RuntimeError(
            f"Found {len(missing_inside)} coordinate frames whose MD steps fall "
            "inside the state-file window but have no exact state assignment. "
            f"Examples: {examples}. Check --first-frame-step and "
            "--frame-step-interval. Nearest-neighbor state matching is intentionally "
            "not used."
        )

    return by_trajectory, flat_rows, matched_count, outside_count


def clean_state_labels(
    assignments_by_trajectory: Dict[str, List[dict]],
) -> List[str]:
    code_to_label: Dict[int, str] = {}
    for rows in assignments_by_trajectory.values():
        for row in rows:
            if not row.get("matched") or not row.get("clean"):
                continue
            code = int(row["state_code"])
            label = str(row["state_label"])
            if code in code_to_label and code_to_label[code] != label:
                raise RuntimeError(
                    f"State code {code} maps to both {code_to_label[code]!r} and "
                    f"{label!r}"
                )
            code_to_label[code] = label
    return [code_to_label[code] for code in sorted(code_to_label)]


def state_labels_for_rows(
    assignments: List[dict],
    nrows: int,
    trajectory: str,
    metric: str,
) -> np.ndarray:
    if len(assignments) != nrows:
        raise RuntimeError(
            f"Frame-count mismatch for {trajectory} {metric}: "
            f"{nrows} CPPTRAJ rows but {len(assignments)} frame assignments"
        )
    return np.asarray([row["state_label"] for row in assignments], dtype=object)


def build_state_conditioned_summary(
    analysis_dir: Path,
    metadata: List[dict],
    residues: List[dict],
    assignments_by_trajectory: Dict[str, List[dict]],
    blocks: int,
) -> List[dict]:
    rows: List[dict] = []
    residue_ids = [int(r["resid"]) for r in residues]
    resname = {int(r["resid"]): r["resname"] for r in residues}
    states = clean_state_labels(assignments_by_trajectory)

    def summarize_metric(
        trajectory: str,
        ph: float,
        metric: str,
        resid: int,
        values: np.ndarray,
        assignments: List[dict],
    ) -> None:
        labels = state_labels_for_rows(
            assignments, int(values.size), trajectory, metric
        )
        code_by_label = {
            str(row["state_label"]): int(row["state_code"])
            for row in assignments
            if row.get("matched") and row.get("clean")
        }
        for state in states:
            stat = summarize_state_masked(values, labels, state, blocks)
            if stat["nframes"] == 0:
                continue
            rows.append(
                {
                    "trajectory": trajectory,
                    "ph": ph,
                    "state": state,
                    "state_code": code_by_label.get(state, ""),
                    "resid": resid,
                    "resname": resname.get(resid, "RES"),
                    "metric": metric,
                    **stat,
                }
            )

    for meta in metadata:
        label = meta["label"]
        ph_text = meta.get("ph", "").strip()
        ph = float(ph_text) if ph_text else np.nan
        run_dir = analysis_dir / label
        assignments = assignments_by_trajectory[label]

        metric_files = {
            "sasa_residue_A2": run_dir / "sasa_residue.dat",
            "sasa_sidechain_A2": run_dir / "sasa_sidechain.dat",
            "hbond_water_count": run_dir / "hbond_water.dat",
        }

        for metric, path in metric_files.items():
            arr = read_numeric(path)
            if arr.shape[1] < 2:
                continue

            if metric == "sasa_sidechain_A2":
                expected = [r for r in residue_ids if resname.get(r, "") != "GLY"]
            else:
                expected = residue_ids

            ncols = min(len(expected), arr.shape[1] - 1)
            for i in range(ncols):
                summarize_metric(
                    label,
                    ph,
                    metric,
                    expected[i],
                    arr[:, i + 1],
                    assignments,
                )

        donor = read_numeric(run_dir / "hbond_protein_donor.dat")
        acceptor = read_numeric(run_dir / "hbond_protein_acceptor.dat")
        if donor.shape[1] >= 2 and acceptor.shape[1] >= 2:
            ncols = min(len(residue_ids), donor.shape[1] - 1, acceptor.shape[1] - 1)
            nframes = min(donor.shape[0], acceptor.shape[0])
            if donor.shape[0] != acceptor.shape[0]:
                raise RuntimeError(
                    f"Protein H-bond donor/acceptor frame-count mismatch for {label}: "
                    f"{donor.shape[0]} vs {acceptor.shape[0]}"
                )
            for i in range(ncols):
                vals = donor[:nframes, i + 1] + acceptor[:nframes, i + 1]
                summarize_metric(
                    label,
                    ph,
                    "hbond_protein_count",
                    residue_ids[i],
                    vals,
                    assignments,
                )

        for resid in residue_ids:
            ws = read_numeric(run_dir / f"watershell_R{resid}.dat")
            if ws.shape[1] >= 2:
                summarize_metric(
                    label,
                    ph,
                    "water_first_shell_count",
                    resid,
                    ws[:, 1],
                    assignments,
                )
            if ws.shape[1] >= 3:
                summarize_metric(
                    label,
                    ph,
                    "water_second_shell_count",
                    resid,
                    ws[:, 2],
                    assignments,
                )

    return rows


def build_state_count_rows(
    metadata: List[dict],
    assignments_by_trajectory: Dict[str, List[dict]],
) -> List[dict]:
    state_info: Dict[str, Tuple[int, int]] = {}
    for rows in assignments_by_trajectory.values():
        for row in rows:
            if not row.get("matched"):
                continue
            label = str(row["state_label"])
            state_info[label] = (int(row["state_code"]), int(row["clean"]))

    ordered_states = [
        label
        for label, _ in sorted(state_info.items(), key=lambda item: item[1][0])
    ]
    result: List[dict] = []

    for meta in metadata:
        trajectory = meta["label"]
        ph_text = meta.get("ph", "").strip()
        ph = float(ph_text) if ph_text else np.nan
        assignments = assignments_by_trajectory[trajectory]
        matched = [row for row in assignments if row.get("matched")]
        clean_rows = [row for row in matched if row.get("clean")]
        matched_total = len(matched)
        clean_total = len(clean_rows)

        for state in ordered_states:
            code, clean = state_info[state]
            count = sum(1 for row in matched if row["state_label"] == state)
            result.append(
                {
                    "trajectory": trajectory,
                    "ph": ph,
                    "state": state,
                    "state_code": code,
                    "clean": clean,
                    "nframes": count,
                    "matched_frames": matched_total,
                    "clean_frames": clean_total,
                    "fraction_of_matched": (
                        count / matched_total if matched_total else np.nan
                    ),
                    "fraction_of_clean": (
                        count / clean_total if clean and clean_total else np.nan
                    ),
                }
            )

    return result


def build_summary(
    analysis_dir: Path,
    metadata: List[dict],
    residues: List[dict],
    blocks: int,
) -> Tuple[List[dict], List[dict]]:
    rows: List[dict] = []
    global_rows: List[dict] = []

    residue_ids = [int(r["resid"]) for r in residues]
    resname = {int(r["resid"]): r["resname"] for r in residues}

    for meta in metadata:
        label = meta["label"]
        ph_text = meta.get("ph", "").strip()
        ph = float(ph_text) if ph_text else np.nan
        run_dir = analysis_dir / label

        def add_metric(metric: str, resid: int, values: np.ndarray) -> None:
            stat = summarize(values, blocks)
            rows.append(
                {
                    "trajectory": label,
                    "ph": ph,
                    "resid": resid,
                    "resname": resname.get(resid, "RES"),
                    "metric": metric,
                    **stat,
                }
            )

        # Compact CPPTRAJ 'create' files: col 0 is frame/index; subsequent
        # columns follow the residue order used by analyze_arex.sh.
        metric_files = {
            "sasa_residue_A2": run_dir / "sasa_residue.dat",
            "sasa_sidechain_A2": run_dir / "sasa_sidechain.dat",
            "hbond_water_count": run_dir / "hbond_water.dat",
        }

        for metric, path in metric_files.items():
            arr = read_numeric(path)
            if arr.shape[1] < 2:
                continue

            if metric == "sasa_sidechain_A2":
                # Gly side-chain SASA is intentionally omitted by the driver.
                expected = [r for r in residue_ids if resname.get(r, "") != "GLY"]
            else:
                expected = residue_ids

            ncols = min(len(expected), arr.shape[1] - 1)
            for i in range(ncols):
                add_metric(metric, expected[i], arr[:, i + 1])

        # Protein H-bonds are the sum of site-as-donor and site-as-acceptor
        # directions. The site/rest masks are disjoint, so the two directions
        # do not double-count a single H-bond.
        donor = read_numeric(run_dir / "hbond_protein_donor.dat")
        acceptor = read_numeric(run_dir / "hbond_protein_acceptor.dat")
        if donor.shape[1] >= 2 and acceptor.shape[1] >= 2:
            ncols = min(len(residue_ids), donor.shape[1] - 1, acceptor.shape[1] - 1)
            nframes = min(donor.shape[0], acceptor.shape[0])
            for i in range(ncols):
                vals = donor[:nframes, i + 1] + acceptor[:nframes, i + 1]
                add_metric("hbond_protein_count", residue_ids[i], vals)

        # Watershell output is one file per residue. Expected numeric columns:
        # frame/index, first-shell count, second-shell count.
        for resid in residue_ids:
            ws = read_numeric(run_dir / f"watershell_R{resid}.dat")
            if ws.shape[1] >= 2:
                add_metric("water_first_shell_count", resid, ws[:, 1])
            if ws.shape[1] >= 3:
                add_metric("water_second_shell_count", resid, ws[:, 2])

        # Global structural controls.
        frame_dt_ps = float(meta.get("frame_dt_ps", "1.0") or 1.0)
        for metric, filename in [
            ("backbone_rmsd_A", "rmsd.dat"),
            ("radius_of_gyration_A", "rg.dat"),
        ]:
            arr = read_numeric(run_dir / filename)
            if arr.shape[1] >= 2:
                vals = arr[:, 1]
                stat = summarize(vals, blocks)
                global_rows.append(
                    {
                        "trajectory": label,
                        "ph": ph,
                        "metric": metric,
                        "frame_dt_ps": frame_dt_ps,
                        **stat,
                    }
                )

    return rows, global_rows


def write_summary(path: Path, rows: List[dict], fields: List[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for row in rows:
            out = dict(row)
            for key in ("ph", "mean", "sd", "block_sem"):
                if key in out and isinstance(out[key], float):
                    out[key] = "" if not np.isfinite(out[key]) else f"{out[key]:.8g}"
            writer.writerow(out)


def write_rows_tsv(
    path: Path,
    rows: List[dict],
    fields: List[str],
    float_fields: Optional[Iterable[str]] = None,
) -> None:
    float_keys = set(float_fields or [])
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for row in rows:
            out = dict(row)
            for key in float_keys:
                if key in out and isinstance(out[key], (float, np.floating)):
                    value = float(out[key])
                    out[key] = "" if not math.isfinite(value) else f"{value:.10g}"
            writer.writerow({field: out.get(field, "") for field in fields})


METRIC_LABELS = {
    "sasa_residue_A2": r"Residue SASA ($\mathrm{\AA^2}$)",
    "sasa_sidechain_A2": r"Side-chain SASA ($\mathrm{\AA^2}$)",
    "water_first_shell_count": "First-shell waters",
    "water_second_shell_count": "Second-shell waters",
    "hbond_water_count": "Site-water H-bonds",
    "hbond_protein_count": "Site-protein H-bonds",
}


def plot_residue_metric(
    rows: List[dict],
    metric: str,
    trajectory_labels: List[str],
    residues: List[dict],
    figures_dir: Path,
    max_ph_lines: int,
) -> None:
    data = [r for r in rows if r["metric"] == metric and np.isfinite(r["mean"])]
    if not data:
        return

    residue_ids = [int(r["resid"]) for r in residues]
    resname = {int(r["resid"]): r["resname"] for r in residues}

    # Lookup table
    lookup: Dict[Tuple[str, int], dict] = {
        (r["trajectory"], int(r["resid"])): r for r in data
    }

    use_ph_axis = all(
        np.isfinite(r["ph"]) for r in data
    ) and len(residue_ids) <= max_ph_lines

    fig, ax = plt.subplots(figsize=(6.8, 4.1))

    if use_ph_axis:
        # One line per residue; best for targeted HEWL titratable sites.
        for resid in residue_ids:
            pts = [
                lookup[(lab, resid)]
                for lab in trajectory_labels
                if (lab, resid) in lookup
            ]
            if not pts:
                continue
            pts = sorted(pts, key=lambda x: x["ph"])
            x = np.array([p["ph"] for p in pts], dtype=float)
            y = np.array([p["mean"] for p in pts], dtype=float)
            e = np.array([p["block_sem"] for p in pts], dtype=float)
            e[~np.isfinite(e)] = 0.0
            ax.errorbar(
                x,
                y,
                yerr=e,
                marker="o",
                capsize=2.0,
                label=f"{resname.get(resid, 'RES')}{resid}",
            )
        ax.set_xlabel("pH")
        ax.legend(
            frameon=False,
            ncol=2 if len(residue_ids) > 6 else 1,
            title="Residue",
            title_fontsize=7.5,
        )
    else:
        # One line per pH/trajectory; better for all-residue profiles.
        for lab in trajectory_labels:
            pts = [
                lookup[(lab, resid)]
                for resid in residue_ids
                if (lab, resid) in lookup
            ]
            if not pts:
                continue
            x = np.array([p["resid"] for p in pts], dtype=float)
            y = np.array([p["mean"] for p in pts], dtype=float)
            ph = pts[0]["ph"]
            legend = f"pH {ph:g}" if np.isfinite(ph) else lab
            ax.plot(x, y, marker="o", label=legend)
        ax.set_xlabel("Residue number")
        ax.legend(
            frameon=False,
            ncol=2,
            title="Trajectory",
            title_fontsize=7.5,
        )

    ax.set_ylabel(METRIC_LABELS.get(metric, metric))
    clean_axis(ax)
    export_figure(fig, figures_dir / metric)


def plot_global_timeseries(
    analysis_dir: Path,
    metadata: List[dict],
    filename: str,
    ylabel: str,
    stem: str,
    figures_dir: Path,
) -> None:
    fig, ax = plt.subplots(figsize=(6.8, 3.8))
    plotted = 0

    ordered = sorted(
        metadata,
        key=lambda r: (
            0 if r.get("ph", "").strip() else 1,
            float(r["ph"]) if r.get("ph", "").strip() else math.inf,
            r["label"],
        ),
    )

    for meta in ordered:
        arr = read_numeric(analysis_dir / meta["label"] / filename)
        if arr.shape[1] < 2:
            continue
        dt_ns = float(meta.get("frame_dt_ps", "1.0") or 1.0) / 1000.0
        t = np.arange(arr.shape[0], dtype=float) * dt_ns
        ph_text = meta.get("ph", "").strip()
        label = f"pH {float(ph_text):g}" if ph_text else meta["label"]
        ax.plot(t, arr[:, 1], label=label)
        plotted += 1

    if plotted == 0:
        plt.close(fig)
        return

    ax.set_xlabel("Analysis time (ns)")
    ax.set_ylabel(ylabel)
    clean_axis(ax)
    if plotted <= 12:
        ax.legend(frameon=False, ncol=2)
    export_figure(fig, figures_dir / stem)


def main() -> None:
    args = parse_args()
    if args.blocks < 1:
        raise SystemExit("--blocks must be >= 1")

    if args.state_file is not None:
        if args.first_frame_step is None or args.frame_step_interval is None:
            raise SystemExit(
                "--state-file requires both --first-frame-step and "
                "--frame-step-interval"
            )
        if args.first_frame_step < 0:
            raise SystemExit("--first-frame-step must be >= 0")
        if args.frame_step_interval <= 0:
            raise SystemExit("--frame-step-interval must be > 0")
        if not args.state_file.exists():
            raise SystemExit(f"State file not found: {args.state_file}")
    elif args.first_frame_step is not None or args.frame_step_interval is not None:
        raise SystemExit(
            "--first-frame-step/--frame-step-interval are only meaningful with "
            "--state-file"
        )

    analysis_dir = args.analysis_dir
    metadata_path = analysis_dir / "metadata.tsv"
    residues_path = analysis_dir / "residues.tsv"

    if not metadata_path.exists() or not residues_path.exists():
        raise SystemExit(
            "metadata.tsv and residues.tsv are required. Run analyze_arex.sh first."
        )

    configure_matplotlib()

    metadata = read_tsv(metadata_path)
    residues = read_tsv(residues_path)
    labels = trajectory_order(metadata)

    rows, global_rows = build_summary(
        analysis_dir=analysis_dir,
        metadata=metadata,
        residues=residues,
        blocks=args.blocks,
    )

    write_summary(
        analysis_dir / "residue_summary.tsv",
        rows,
        [
            "trajectory",
            "ph",
            "resid",
            "resname",
            "metric",
            "nframes",
            "mean",
            "sd",
            "block_sem",
            "nblocks",
        ],
    )
    write_summary(
        analysis_dir / "global_summary.tsv",
        global_rows,
        [
            "trajectory",
            "ph",
            "metric",
            "frame_dt_ps",
            "nframes",
            "mean",
            "sd",
            "block_sem",
            "nblocks",
        ],
    )

    state_outputs: List[Path] = []
    if args.state_file is not None:
        try:
            state_lookup, state_ranges, state_ph = load_state_file(args.state_file)
            assignments_by_trajectory, assignment_rows, matched_count, outside_count = (
                build_frame_state_assignments(
                    analysis_dir=analysis_dir,
                    metadata=metadata,
                    state_lookup=state_lookup,
                    state_ranges=state_ranges,
                    state_ph=state_ph,
                    first_frame_step=args.first_frame_step,
                    frame_step_interval=args.frame_step_interval,
                )
            )

            assignment_path = analysis_dir / "dyad_frame_assignments.tsv"
            write_rows_tsv(
                assignment_path,
                assignment_rows,
                [
                    "trajectory",
                    "ph",
                    "analysis_row",
                    "trajectory_frame",
                    "step",
                    "absolute_time_ns",
                    "primary_label",
                    "primary_lambda",
                    "coupled_label",
                    "coupled_lambda",
                    "state_code",
                    "state_label",
                    "clean",
                    "matched",
                ],
                float_fields=[
                    "ph",
                    "absolute_time_ns",
                    "primary_lambda",
                    "coupled_lambda",
                ],
            )
            state_outputs.append(assignment_path)

            state_rows = build_state_conditioned_summary(
                analysis_dir=analysis_dir,
                metadata=metadata,
                residues=residues,
                assignments_by_trajectory=assignments_by_trajectory,
                blocks=args.blocks,
            )
            state_summary_path = analysis_dir / "dyad_state_summary.tsv"
            write_rows_tsv(
                state_summary_path,
                state_rows,
                [
                    "trajectory",
                    "ph",
                    "state",
                    "state_code",
                    "resid",
                    "resname",
                    "metric",
                    "nframes",
                    "mean",
                    "sd",
                    "block_sem",
                    "nblocks",
                ],
                float_fields=["ph", "mean", "sd", "block_sem"],
            )
            state_outputs.append(state_summary_path)

            state_count_rows = build_state_count_rows(
                metadata, assignments_by_trajectory
            )
            state_counts_path = analysis_dir / "dyad_state_counts.tsv"
            write_rows_tsv(
                state_counts_path,
                state_count_rows,
                [
                    "trajectory",
                    "ph",
                    "state",
                    "state_code",
                    "clean",
                    "nframes",
                    "matched_frames",
                    "clean_frames",
                    "fraction_of_matched",
                    "fraction_of_clean",
                ],
                float_fields=[
                    "ph",
                    "fraction_of_matched",
                    "fraction_of_clean",
                ],
            )
            state_outputs.append(state_counts_path)

            print(
                f"State join: {matched_count} coordinate frames matched exactly; "
                f"{outside_count} processed frames were outside the state-file window."
            )
        except RuntimeError as error:
            raise SystemExit(f"State-conditioned analysis failed: {error}") from error

    figures_dir = analysis_dir / "figures"
    for metric in METRIC_LABELS:
        plot_residue_metric(
            rows,
            metric,
            labels,
            residues,
            figures_dir,
            args.max_ph_lines,
        )

    plot_global_timeseries(
        analysis_dir,
        metadata,
        "rmsd.dat",
        r"Backbone RMSD ($\mathrm{\AA}$)",
        "backbone_rmsd_timeseries",
        figures_dir,
    )
    plot_global_timeseries(
        analysis_dir,
        metadata,
        "rg.dat",
        r"Radius of gyration ($\mathrm{\AA}$)",
        "radius_of_gyration_timeseries",
        figures_dir,
    )

    print(f"Wrote {analysis_dir / 'residue_summary.tsv'}")
    print(f"Wrote {analysis_dir / 'global_summary.tsv'}")
    for path in state_outputs:
        print(f"Wrote {path}")
    print(f"Figures: {figures_dir}")


if __name__ == "__main__":
    main()