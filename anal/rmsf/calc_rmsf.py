#!/usr/bin/env python3
"""Calculate AMBER per-residue RMSF with CPPTRAJ.

Modes:
  backbone  heavy backbone atoms N, CA, C, O
  sidechain heavy side-chain atoms (all non-H atoms except N, CA, C, O, OXT)
  total     all heavy atoms

Frames are fitted with heavy backbone atoms by default, regardless of the RMSF
mode. This removes global translation/rotation without fitting away side-chain
motion. Residue labels can be shifted and selected residues highlighted.
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import re
import shutil
import statistics
import subprocess
import sys
from pathlib import Path
from typing import List, Optional, Sequence, Tuple

DEFAULT_SELECTION_MASK = "!(:WAT,HOH,SOL,TIP3,TIP3P,TIP4P,OPC,SPC,SPCE)"
MODE_LABELS = {
    "backbone": "Backbone",
    "sidechain": "Side-chain",
    "total": "Total heavy-atom",
}


def fail(message: str, exit_code: int = 2) -> None:
    print(f"ERROR: {message}", file=sys.stderr)
    raise SystemExit(exit_code)


def warn(message: str) -> None:
    print(f"WARNING: {message}", file=sys.stderr)


def existing_file(value: str) -> Path:
    path = Path(value).expanduser()
    if not path.is_file():
        raise argparse.ArgumentTypeError(f"File not found: {value}")
    return path.resolve()


def positive_int(value: str) -> int:
    number = int(value)
    if number <= 0:
        raise argparse.ArgumentTypeError("Value must be greater than zero")
    return number


def validate_mask(mask: str, label: str) -> str:
    mask = mask.strip()
    if not mask:
        fail(f"{label} cannot be empty")
    if "\n" in mask or "\r" in mask:
        fail(f"{label} must be one line")
    return mask


def build_masks(args: argparse.Namespace) -> Tuple[str, str]:
    selection = validate_mask(args.selection_mask, "Selection mask")

    fit_mask = args.fit_mask or f"({selection})&@N,CA,C,O"
    fit_mask = validate_mask(fit_mask, "Fit mask")

    if args.rmsf_mask:
        rmsf_mask = args.rmsf_mask
    elif args.rmsf_mode == "backbone":
        rmsf_mask = f"({selection})&@N,CA,C,O"
    elif args.rmsf_mode == "sidechain":
        rmsf_mask = f"({selection})&(!@N,CA,C,O,OXT)&(!@H=)"
    else:
        rmsf_mask = f"({selection})&(!@H=)"

    return fit_mask, validate_mask(rmsf_mask, "RMSF mask")


def parse_highlight_residues(tokens: Optional[Sequence[str]]) -> List[int]:
    if not tokens:
        return []

    residues = set()
    for token in tokens:
        for item in token.split(","):
            item = item.strip()
            if not item:
                continue

            match = re.fullmatch(r"(-?\d+)\s*-\s*(-?\d+)", item)
            if match:
                start, stop = int(match.group(1)), int(match.group(2))
                if stop < start:
                    fail(f"Invalid residue range '{item}'")
                residues.update(range(start, stop + 1))
                continue

            try:
                residues.add(int(item))
            except ValueError:
                fail(
                    f"Invalid residue specification '{item}'. Use integers, "
                    "comma-separated lists, or ranges such as 7,14,20-25."
                )

    return sorted(residues)


def find_executable(name: str) -> str:
    expanded = os.path.expanduser(name)
    if os.path.sep in expanded:
        path = Path(expanded).resolve()
        if not path.is_file() or not os.access(str(path), os.X_OK):
            fail(f"CPPTRAJ executable is not usable: {path}")
        return str(path)

    located = shutil.which(expanded)
    if located is None:
        fail(
            f"Could not find '{name}' in PATH. Load AMBER/AmberTools or pass "
            "--cpptraj /full/path/to/cpptraj."
        )
    return located


def cpptraj_quote(path: Path) -> str:
    text = str(path.resolve())
    if any(character in text for character in ('\n', '\r', '"')):
        fail(f"Unsupported character in path: {text}")
    return f'"{text}"'


def run_command(command: Sequence[str]) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        list(command),
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )


def cpptraj_version(cpptraj: str) -> str:
    result = run_command([cpptraj, "--version"])
    text = result.stdout.strip()
    return text.splitlines()[0] if text else "unknown"


def parse_selected_numbers(output: str) -> Optional[List[int]]:
    match = re.search(r"Selected=\s*([0-9\s]+)", output)
    if not match:
        return None
    return [int(value) for value in re.findall(r"\d+", match.group(1))]


def selected_numbers(
    cpptraj: str, topology: Path, mask: str, option: str
) -> Optional[List[int]]:
    result = run_command([cpptraj, "-p", str(topology), option, mask])
    if result.returncode != 0:
        fail(f"CPPTRAJ rejected mask {mask!r}:\n{result.stdout.strip()}")
    return parse_selected_numbers(result.stdout)


def selected_atom_count(cpptraj: str, topology: Path, mask: str) -> Optional[int]:
    values = selected_numbers(cpptraj, topology, mask, "-ms")
    return None if values is None else len(values)


def selected_residue_count(cpptraj: str, topology: Path, mask: str) -> Optional[int]:
    values = selected_numbers(cpptraj, topology, mask, "-mr")
    return None if values is None else len(values)


def trajectory_frame_count(cpptraj: str, topology: Path, trajectory: Path) -> int:
    result = run_command([cpptraj, "-p", str(topology), "-y", str(trajectory), "-tl"])
    if result.returncode != 0:
        fail(f"Could not inspect trajectory {trajectory}:\n{result.stdout.strip()}")

    matches = re.findall(r"Frames:\s*(\d+)", result.stdout)
    if not matches:
        fail(f"Could not determine frame count for {trajectory}")
    return int(matches[-1])


def selected_frame_count(
    total_frames: int, start: int, stop: Optional[int], stride: int
) -> int:
    last = total_frames if stop is None else min(stop, total_frames)
    if start > last:
        return 0
    return ((last - start) // stride) + 1


def trajectory_lines(
    trajectories: Sequence[Path], start: int, stop: Optional[int], stride: int
) -> List[str]:
    stop_token = "last" if stop is None else str(stop)
    return [
        f"trajin {cpptraj_quote(trajectory)} {start} {stop_token} {stride}"
        for trajectory in trajectories
    ]


def build_cpptraj_input(
    topology: Path,
    trajectories: Sequence[Path],
    raw_output: Path,
    fit_mask: str,
    rmsf_mask: str,
    reference: Optional[Path],
    fit_reference: str,
    start: int,
    stop: Optional[int],
    stride: int,
    autoimage: bool,
    pdbres: bool,
) -> str:
    lines = [f"parm {cpptraj_quote(topology)}"]

    if reference is not None:
        lines.append(f"reference {cpptraj_quote(reference)} [RMSF_REF]")

    lines.extend(trajectory_lines(trajectories, start, stop, stride))

    if reference is None and fit_reference == "average":
        if autoimage:
            lines.append("autoimage")
        lines.append(f"rms InitialFit {fit_mask} first")
        lines.append("average crdset RMSF_AVERAGE")
        lines.append("run")
        lines.append("clear actions")

        if autoimage:
            lines.append("autoimage")
        lines.append(f"rms BackboneFit {fit_mask} ref RMSF_AVERAGE")
    else:
        if autoimage:
            lines.append("autoimage")
        if reference is not None:
            lines.append(f"rms BackboneFit {fit_mask} ref [RMSF_REF]")
        else:
            lines.append(f"rms BackboneFit {fit_mask} first")

    command = f"atomicfluct ResidueRMSF out {cpptraj_quote(raw_output)} {rmsf_mask} byres"
    if pdbres:
        command += " pdbres"
    lines.extend([command, "run", "quit"])
    return "\n".join(lines) + "\n"


def parse_cpptraj_rmsf(path: Path) -> Tuple[List[float], List[float]]:
    residues: List[float] = []
    values: List[float] = []

    if not path.is_file():
        fail(f"CPPTRAJ did not create the expected RMSF file: {path}")

    for line_number, line in enumerate(path.read_text(errors="replace").splitlines(), 1):
        stripped = line.strip()
        if not stripped or stripped.startswith(("#", "@")):
            continue
        fields = stripped.split()
        if len(fields) < 2:
            continue
        try:
            residue, rmsf = float(fields[0]), float(fields[1])
        except ValueError:
            continue
        if not math.isfinite(residue) or not math.isfinite(rmsf):
            fail(f"Non-finite RMSF data at {path}:{line_number}")
        if rmsf < 0:
            fail(f"Negative RMSF at {path}:{line_number}")
        residues.append(residue)
        values.append(rmsf)

    if not values:
        fail(f"No RMSF values could be parsed from {path}")
    return residues, values


def tidy_residue(value: float):
    rounded = round(value)
    return int(rounded) if abs(value - rounded) < 1.0e-8 else value


def metric_key(mode: str) -> str:
    return {
        "backbone": "backbone_rmsf_angstrom",
        "sidechain": "sidechain_rmsf_angstrom",
        "total": "total_heavy_atom_rmsf_angstrom",
    }[mode]


def write_csv(
    path: Path,
    source_residues: Sequence[float],
    plotted_residues: Sequence[float],
    values: Sequence[float],
    mode: str,
) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["source_residue", "plotted_residue", metric_key(mode)])
        for source, plotted, rmsf in zip(source_residues, plotted_residues, values):
            writer.writerow([tidy_residue(source), tidy_residue(plotted), f"{rmsf:.8f}"])


def make_plot(
    path: Path,
    residues: Sequence[float],
    values: Sequence[float],
    reference_label: str,
    mode: str,
    pdbres: bool,
    residue_offset: int,
    highlighted_residues: Sequence[int],
    dpi: int,
) -> None:
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        fail("matplotlib is required. Install it with: python -m pip install matplotlib")

    label = MODE_LABELS[mode]
    mean_rmsf = statistics.mean(values)
    if pdbres:
        x_label = "PDB residue number"
    elif residue_offset:
        x_label = "Corrected residue number"
    else:
        x_label = "Topology residue index"

    figure, axis = plt.subplots(figsize=(8.0, 4.8))
    axis.plot(residues, values, linewidth=1.1, label=f"{label} RMSF")
    axis.axhline(
        mean_rmsf,
        linewidth=1.0,
        linestyle="--",
        alpha=0.7,
        label=f"Mean = {mean_rmsf:.2f} Å",
    )

    lookup = {tidy_residue(r): v for r, v in zip(residues, values)}
    first_marker = True
    for residue in highlighted_residues:
        if residue not in lookup:
            continue
        rmsf = lookup[residue]
        axis.axvline(residue, linewidth=0.8, linestyle=":", alpha=0.45)
        axis.scatter(
            [residue],
            [rmsf],
            s=42,
            zorder=4,
            label="Highlighted residue(s)" if first_marker else None,
        )
        axis.annotate(
            str(residue),
            (residue, rmsf),
            xytext=(0, 8),
            textcoords="offset points",
            ha="center",
            va="bottom",
            fontsize=8,
        )
        first_marker = False

    axis.set_xlabel(x_label)
    axis.set_ylabel("RMSF (Å)")
    axis.set_title(f"{label} RMSF vs residue")
    axis.text(
        0.01,
        0.99,
        f"Backbone-fit reference: {reference_label}",
        transform=axis.transAxes,
        ha="left",
        va="top",
        fontsize=9,
    )
    axis.grid(True, alpha=0.25)
    axis.legend(frameon=False)
    figure.tight_layout()
    figure.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close(figure)


def write_summary(
    path: Path,
    args: argparse.Namespace,
    cpptraj_path: str,
    version: str,
    fit_mask: str,
    rmsf_mask: str,
    fit_atom_count: Optional[int],
    fit_residue_count: Optional[int],
    rmsf_atom_count: Optional[int],
    rmsf_residue_count: Optional[int],
    analyzed_frames: int,
    residues: Sequence[float],
    values: Sequence[float],
) -> None:
    reference_label = str(args.reference) if args.reference else args.fit_reference
    ranked = sorted(zip(residues, values), key=lambda item: item[1], reverse=True)
    top_n = min(args.top_residues, len(ranked))
    label = MODE_LABELS[args.rmsf_mode]

    lines = [
        f"AMBER {label.lower()} RMSF analysis",
        "",
        f"CPPTRAJ: {cpptraj_path}",
        f"CPPTRAJ version: {version}",
        f"Topology: {args.topology}",
        f"Trajectories: {', '.join(str(path) for path in args.trajectory)}",
        f"RMSF mode: {args.rmsf_mode}",
        f"Fit reference: {reference_label}",
        f"Selection mask: {args.selection_mask}",
        f"Backbone fit mask: {fit_mask}",
        f"RMSF mask: {rmsf_mask}",
        f"Fit atoms/residues: {fit_atom_count if fit_atom_count is not None else 'not determined'} / {fit_residue_count if fit_residue_count is not None else 'not determined'}",
        f"RMSF atoms/residues: {rmsf_atom_count if rmsf_atom_count is not None else 'not determined'} / {rmsf_residue_count if rmsf_residue_count is not None else 'not determined'}",
        f"Autoimage: {'yes' if not args.no_autoimage else 'no'}",
        f"Start/stop/stride per trajectory: {args.start}/{args.stop or 'last'}/{args.stride}",
        f"Analyzed frames: {analyzed_frames}",
        f"Residue numbering source: {'original PDB numbering' if args.pdbres else 'topology indices'}",
        f"Residue offset applied to CSV/plot labels: {args.residue_offset:+d}",
        f"Highlighted plotted residues: {', '.join(map(str, args.highlight_residues)) if args.highlight_residues else 'none'}",
        "By-residue aggregation: CPPTRAJ mass-weighted average of selected-atom RMSF",
        "",
        f"Mean RMSF (Å): {statistics.mean(values):.6f}",
        f"Median RMSF (Å): {statistics.median(values):.6f}",
        f"Minimum RMSF (Å): {min(values):.6f}",
        f"Maximum RMSF (Å): {max(values):.6f}",
        "",
        f"Top {top_n} RMSF residues:",
    ]
    for residue, rmsf in ranked[:top_n]:
        lines.append(f"  {tidy_residue(residue)}\t{rmsf:.6f} Å")

    if args.rmsf_mode == "sidechain":
        lines.extend(
            [
                "",
                "Note: residues with no selected side-chain heavy atoms (notably glycine)",
                "do not receive a side-chain RMSF value and may be absent from output.",
            ]
        )

    path.write_text("\n".join(lines) + "\n")


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Use CPPTRAJ to backbone-fit AMBER trajectories and calculate "
            "backbone, side-chain, or total heavy-atom RMSF by residue."
        )
    )
    parser.add_argument("-p", "--topology", required=True, type=existing_file)
    parser.add_argument(
        "-y",
        "--trajectory",
        required=True,
        nargs="+",
        type=existing_file,
        help="One or more trajectory files, in chronological order",
    )
    parser.add_argument(
        "-o", "--output-prefix", default="residue_rmsf", help="Default: residue_rmsf"
    )
    parser.add_argument(
        "--rmsf-mode",
        choices=("backbone", "sidechain", "total"),
        default="backbone",
        help="Atoms used for per-residue RMSF; default: backbone",
    )
    parser.add_argument(
        "--selection-mask",
        default=DEFAULT_SELECTION_MASK,
        help=(
            "Residues/solute included before atom filtering. For HEWL, prefer an "
            "explicit protein range such as ':1-129'."
        ),
    )
    parser.add_argument(
        "--fit-mask",
        help="Optional explicit fit mask; default is selection mask intersected with @N,CA,C,O",
    )
    parser.add_argument(
        "--rmsf-mask",
        help="Optional explicit RMSF mask; overrides --rmsf-mode",
    )
    parser.add_argument(
        "-r", "--reference", type=existing_file, help="Optional external fit reference"
    )
    parser.add_argument(
        "--fit-reference",
        choices=("average", "first"),
        default="average",
        help="Default: average",
    )
    parser.add_argument("--start", type=positive_int, default=1)
    parser.add_argument("--stop", type=positive_int)
    parser.add_argument("--stride", type=positive_int, default=1)
    parser.add_argument("--no-autoimage", action="store_true")
    parser.add_argument(
        "--pdbres", action="store_true", help="Use original PDB residue numbering"
    )
    parser.add_argument(
        "--residue-offset",
        type=int,
        default=0,
        help="Integer added to CPPTRAJ residue labels; default: 0",
    )
    parser.add_argument(
        "--highlight-residues",
        nargs="+",
        metavar="RESIDUE",
        help="Corrected/visible residues, e.g. 7 14 20-25 or 7,14,20-25",
    )
    parser.add_argument("--top-residues", type=positive_int, default=10)
    parser.add_argument("--dpi", type=positive_int, default=300)
    parser.add_argument("--cpptraj", default="cpptraj")
    return parser.parse_args()


def main() -> None:
    args = parse_arguments()
    args.highlight_residues = parse_highlight_residues(args.highlight_residues)
    fit_mask, rmsf_mask = build_masks(args)

    if args.stop is not None and args.stop < args.start:
        fail("--stop must be greater than or equal to --start")

    cpptraj = find_executable(args.cpptraj)
    version = cpptraj_version(cpptraj)

    output_prefix = Path(args.output_prefix).expanduser()
    if not output_prefix.is_absolute():
        output_prefix = (Path.cwd() / output_prefix).resolve()
    output_prefix.parent.mkdir(parents=True, exist_ok=True)

    input_path = Path(str(output_prefix) + ".cpptraj.in")
    log_path = Path(str(output_prefix) + ".cpptraj.log")
    raw_path = Path(str(output_prefix) + ".cpptraj.dat")
    csv_path = Path(str(output_prefix) + ".csv")
    plot_path = Path(str(output_prefix) + ".png")
    summary_path = Path(str(output_prefix) + ".summary.txt")

    fit_atom_count = selected_atom_count(cpptraj, args.topology, fit_mask)
    fit_residue_count = selected_residue_count(cpptraj, args.topology, fit_mask)
    rmsf_atom_count = selected_atom_count(cpptraj, args.topology, rmsf_mask)
    rmsf_residue_count = selected_residue_count(cpptraj, args.topology, rmsf_mask)

    if fit_atom_count == 0 or fit_residue_count == 0:
        fail(f"The fit mask selected no usable atoms/residues: {fit_mask}")
    if rmsf_atom_count == 0 or rmsf_residue_count == 0:
        fail(f"The RMSF mask selected no usable atoms/residues: {rmsf_mask}")

    print(f"Fit mask:  {fit_mask}")
    print(f"RMSF mask: {rmsf_mask}")
    if fit_atom_count is not None and fit_residue_count is not None:
        print(f"Fit selection: {fit_atom_count} atoms in {fit_residue_count} residues")
        expected = 4 * fit_residue_count
        if fit_atom_count != expected:
            warn(
                f"Fit mask selected {fit_atom_count} atoms across {fit_residue_count} "
                f"residues; four atoms per residue would be {expected}."
            )
    if rmsf_atom_count is not None and rmsf_residue_count is not None:
        print(f"RMSF selection: {rmsf_atom_count} atoms in {rmsf_residue_count} residues")

    analyzed_frames = 0
    for trajectory in args.trajectory:
        total = trajectory_frame_count(cpptraj, args.topology, trajectory)
        selected = selected_frame_count(total, args.start, args.stop, args.stride)
        if selected == 0:
            fail(
                f"No frames selected from {trajectory} with start={args.start}, "
                f"stop={args.stop or 'last'}, stride={args.stride}"
            )
        analyzed_frames += selected
        print(f"{trajectory.name}: {total} total, {selected} selected")

    if analyzed_frames < 2:
        fail("RMSF requires at least two analyzed frames")

    input_path.write_text(
        build_cpptraj_input(
            topology=args.topology,
            trajectories=args.trajectory,
            raw_output=raw_path,
            fit_mask=fit_mask,
            rmsf_mask=rmsf_mask,
            reference=args.reference,
            fit_reference=args.fit_reference,
            start=args.start,
            stop=args.stop,
            stride=args.stride,
            autoimage=not args.no_autoimage,
            pdbres=args.pdbres,
        )
    )

    result = run_command([cpptraj, "-i", str(input_path)])
    log_path.write_text(result.stdout)
    if result.returncode != 0:
        tail = "\n".join(result.stdout.splitlines()[-40:])
        fail(f"CPPTRAJ failed. Log: {log_path}\n\n{tail}", exit_code=1)

    residues, values = parse_cpptraj_rmsf(raw_path)
    if rmsf_residue_count is not None and len(values) != rmsf_residue_count:
        warn(
            f"Expected {rmsf_residue_count} residue rows but parsed {len(values)}. "
            f"Check {log_path}."
        )
    if len(set(residues)) != len(residues):
        warn("Duplicate residue labels found; this can occur with --pdbres in multi-chain systems.")

    plotted_residues = [residue + args.residue_offset for residue in residues]
    write_csv(csv_path, residues, plotted_residues, values, args.rmsf_mode)

    lookup = {tidy_residue(r): v for r, v in zip(plotted_residues, values)}
    missing = [r for r in args.highlight_residues if r not in lookup]
    if missing:
        warn("Highlighted residue(s) not present: " + ", ".join(map(str, missing)))
    matched = [r for r in args.highlight_residues if r in lookup]
    if matched:
        print("Highlighted residues:")
        for residue in matched:
            print(f"  {residue}: {lookup[residue]:.6f} Å")

    reference_label = (
        str(args.reference)
        if args.reference is not None
        else "trajectory average"
        if args.fit_reference == "average"
        else "first analyzed frame"
    )
    make_plot(
        plot_path,
        plotted_residues,
        values,
        reference_label,
        args.rmsf_mode,
        args.pdbres,
        args.residue_offset,
        args.highlight_residues,
        args.dpi,
    )
    write_summary(
        summary_path,
        args,
        cpptraj,
        version,
        fit_mask,
        rmsf_mask,
        fit_atom_count,
        fit_residue_count,
        rmsf_atom_count,
        rmsf_residue_count,
        analyzed_frames,
        plotted_residues,
        values,
    )

    if args.rmsf_mode == "sidechain":
        warn("Glycine has no side-chain heavy atoms, so glycine residues may be absent.")

    print("\nCreated:")
    for path in (plot_path, csv_path, summary_path, input_path, log_path, raw_path):
        print(f"  {path}")


if __name__ == "__main__":
    main()
