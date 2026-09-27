#!/usr/bin/env python3
"""Calculate AMBER trajectory RMSD with separate fitting and measurement masks.

Workflow in CPPTRAJ:
  1. Optionally autoimage each frame.
  2. Best-fit each frame to the reference using --fit-mask.
  3. Calculate no-fit RMSD on the already aligned coordinates using --rmsd-mask.

The default measurement mask is the fit mask, reproducing a conventional
single-mask RMSD calculation.
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

DEFAULT_FIT_MASK = (
    "(@N,CA,C,O)&!(:WAT,HOH,SOL,TIP3,TIP3P,TIP4P,OPC,SPC,SPCE)"
)


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


def positive_float(value: str) -> float:
    number = float(value)
    if number <= 0:
        raise argparse.ArgumentTypeError("Value must be greater than zero")
    return number


def positive_int(value: str) -> int:
    number = int(value)
    if number <= 0:
        raise argparse.ArgumentTypeError("Value must be greater than zero")
    return number


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
    if "\n" in text or "\r" in text or '"' in text:
        fail(f"Unsupported character in path: {text}")
    return f'"{text}"'


def validate_mask(mask: str, option_name: str) -> str:
    cleaned = mask.strip()
    if not cleaned:
        fail(f"{option_name} cannot be empty")
    if "\n" in cleaned or "\r" in cleaned:
        fail(f"{option_name} must be one line")
    return cleaned


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


def selected_atom_count(
    cpptraj: str, topology: Path, mask: str, option_name: str
) -> Optional[int]:
    result = run_command([cpptraj, "-p", str(topology), "-ms", mask])
    if result.returncode != 0:
        fail(
            f"CPPTRAJ rejected {option_name} or the topology:\n"
            f"{result.stdout.strip()}"
        )

    # Typical CPPTRAJ output contains a line such as:
    # Selected= 1 2 3 ...
    match = re.search(r"Selected=\s*([0-9\s]+)", result.stdout)
    if not match:
        return None
    return len(re.findall(r"\d+", match.group(1)))


def trajectory_frame_count(cpptraj: str, topology: Path, trajectory: Path) -> int:
    result = run_command(
        [cpptraj, "-p", str(topology), "-y", str(trajectory), "-tl"]
    )
    if result.returncode != 0:
        fail(f"Could not inspect trajectory {trajectory}:\n{result.stdout.strip()}")

    matches = re.findall(r"Frames:\s*(\d+)", result.stdout)
    if not matches:
        fail(f"Could not determine frame count for {trajectory}")
    return int(matches[-1])


def processed_local_frames(
    total_frames: int, start: int, stop: Optional[int], stride: int
) -> List[int]:
    last = total_frames if stop is None else min(stop, total_frames)
    if start > last:
        return []
    return list(range(start, last + 1, stride))


def parse_amber_timing(path: Path) -> Tuple[float, int]:
    """Read dt (ps/step) and ntwx (steps/saved frame) from mdin/mdout text."""
    text = path.read_text(errors="replace")
    text = re.sub(r"!.*$", "", text, flags=re.MULTILINE)
    number = r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eEdD][+-]?\d+)?"

    def last_value(key: str) -> Optional[str]:
        matches = re.findall(
            rf"\b{key}\s*=\s*({number})", text, flags=re.IGNORECASE
        )
        return matches[-1] if matches else None

    dt_text = last_value("dt")
    ntwx_text = last_value("ntwx")
    if dt_text is None or ntwx_text is None:
        fail(f"Could not find both dt and ntwx in {path}")

    dt_ps = float(dt_text.replace("D", "E").replace("d", "e"))
    ntwx = int(round(float(ntwx_text.replace("D", "E").replace("d", "e"))))
    if dt_ps <= 0 or ntwx <= 0:
        fail(f"Invalid timing values in {path}: dt={dt_ps}, ntwx={ntwx}")
    return dt_ps, ntwx


def reference_clause(reference: Optional[Path]) -> str:
    return "ref [RMSD_REF]" if reference is not None else "first"


def build_cpptraj_input(
    topology: Path,
    trajectories: Sequence[Path],
    raw_output: Path,
    fit_mask: str,
    rmsd_mask: str,
    reference: Optional[Path],
    start: int,
    stop: Optional[int],
    stride: int,
    autoimage: bool,
    mass_weighted: bool,
) -> str:
    lines = [f"parm {cpptraj_quote(topology)}"]

    if reference is not None:
        lines.append(f"reference {cpptraj_quote(reference)} [RMSD_REF]")

    stop_token = "last" if stop is None else str(stop)
    for trajectory in trajectories:
        lines.append(
            f"trajin {cpptraj_quote(trajectory)} {start} {stop_token} {stride}"
        )

    if autoimage:
        lines.append("autoimage")

    ref = reference_clause(reference)
    mass = " mass" if mass_weighted else ""

    # This first action must modify the coordinates; do not use nomod/nofit here.
    lines.append(f"rmsd FitToReference {fit_mask} {ref}{mass}")

    # Coordinates are now aligned by FitToReference. Calculate the requested
    # RMSD without a second fit, so motion is measured in the fit-mask frame.
    lines.append(
        f"rmsd MeasuredRMSD {rmsd_mask} {ref} nofit "
        f"out {cpptraj_quote(raw_output)}{mass}"
    )

    lines.extend(["run", "quit"])
    return "\n".join(lines) + "\n"


def parse_cpptraj_rmsd(path: Path) -> Tuple[List[int], List[float]]:
    frames: List[int] = []
    values: List[float] = []

    if not path.is_file():
        fail(f"CPPTRAJ did not create the expected RMSD file: {path}")

    for line_number, line in enumerate(path.read_text(errors="replace").splitlines(), 1):
        stripped = line.strip()
        if not stripped or stripped.startswith(("#", "@")):
            continue
        fields = stripped.split()
        if len(fields) < 2:
            continue
        try:
            frame = int(round(float(fields[0])))
            rmsd = float(fields[1])
        except ValueError:
            continue
        if not math.isfinite(rmsd):
            fail(f"Non-finite RMSD at {path}:{line_number}")
        frames.append(frame)
        values.append(rmsd)

    if not values:
        fail(f"No RMSD values could be parsed from {path}")
    return frames, values


def moving_average(values: Sequence[float], window: int) -> List[float]:
    if window <= 1:
        return list(values)

    output: List[float] = []
    running_sum = 0.0
    for index, value in enumerate(values):
        running_sum += value
        if index >= window:
            running_sum -= values[index - window]
        output.append(running_sum / min(index + 1, window))
    return output


def write_csv(
    path: Path,
    cpptraj_frames: Sequence[int],
    values: Sequence[float],
    times_ps: Optional[Sequence[float]],
    provenance: Sequence[Tuple[str, int]],
) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "processed_frame",
                "cpptraj_frame",
                "source_trajectory",
                "source_frame",
                "time_ps",
                "time_ns",
                "rmsd_angstrom",
            ]
        )

        for index, (cpptraj_frame, rmsd) in enumerate(
            zip(cpptraj_frames, values), 1
        ):
            if index - 1 < len(provenance):
                source_name, source_frame = provenance[index - 1]
            else:
                source_name, source_frame = "", ""

            if times_ps is None:
                time_ps = ""
                time_ns = ""
            else:
                time_ps = f"{times_ps[index - 1]:.8f}"
                time_ns = f"{times_ps[index - 1] / 1000.0:.8f}"

            writer.writerow(
                [
                    index,
                    cpptraj_frame,
                    source_name,
                    source_frame,
                    time_ps,
                    time_ns,
                    f"{rmsd:.8f}",
                ]
            )


def make_plot(
    path: Path,
    values: Sequence[float],
    times_ps: Optional[Sequence[float]],
    reference_label: str,
    smooth_window: int,
    dpi: int,
) -> None:
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        fail("matplotlib is required. Install it with: python -m pip install matplotlib")

    if times_ps is None:
        x_values = list(range(1, len(values) + 1))
        x_label = "Processed frame"
        title_suffix = "frame"
    else:
        x_values = [value / 1000.0 for value in times_ps]
        x_label = "Time (ns)"
        title_suffix = "time"

    figure, axis = plt.subplots(figsize=(8.0, 4.8))
    axis.plot(x_values, values, linewidth=1.0, label="Measured RMSD")

    if smooth_window > 1:
        axis.plot(
            x_values,
            moving_average(values, smooth_window),
            linewidth=1.8,
            label=f"Moving average ({smooth_window} frames)",
        )
        axis.legend(frameon=False)

    axis.set_xlabel(x_label)
    axis.set_ylabel("RMSD (Å)")
    axis.set_title(f"RMSD vs {title_suffix}")
    axis.text(
        0.01,
        0.99,
        f"Reference: {reference_label}",
        transform=axis.transAxes,
        ha="left",
        va="top",
        fontsize=9,
    )
    axis.grid(True, alpha=0.25)
    figure.tight_layout()
    figure.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close(figure)


def write_summary(
    path: Path,
    args: argparse.Namespace,
    cpptraj_path: str,
    version: str,
    fit_atom_count: Optional[int],
    rmsd_atom_count: Optional[int],
    values: Sequence[float],
    time_per_saved_frame_ps: Optional[float],
    jump_indices: Sequence[int],
) -> None:
    reference_label = str(args.reference) if args.reference else "first analyzed frame"
    lines = [
        "AMBER RMSD analysis with separate fit and measurement masks",
        "",
        f"CPPTRAJ: {cpptraj_path}",
        f"CPPTRAJ version: {version}",
        f"Topology: {args.topology}",
        f"Trajectories: {', '.join(str(path) for path in args.trajectory)}",
        f"Reference: {reference_label}",
        f"Fit mask: {args.fit_mask}",
        f"RMSD mask: {args.rmsd_mask}",
        f"Fit atoms: {fit_atom_count if fit_atom_count is not None else 'not determined'}",
        f"RMSD atoms: {rmsd_atom_count if rmsd_atom_count is not None else 'not determined'}",
        f"Autoimage: {'yes' if not args.no_autoimage else 'no'}",
        f"Mass-weighted fit and RMSD: {'yes' if args.mass_weighted else 'no'}",
        f"Start/stop/stride per trajectory: {args.start}/{args.stop or 'last'}/{args.stride}",
        f"Processed frames: {len(values)}",
        f"Time per saved input frame (ps): "
        f"{time_per_saved_frame_ps if time_per_saved_frame_ps is not None else 'not supplied'}",
        f"Time per analyzed frame (ps): "
        f"{time_per_saved_frame_ps * args.stride if time_per_saved_frame_ps is not None else 'not supplied'}",
        "",
        f"Mean RMSD (Å): {statistics.mean(values):.6f}",
        f"Median RMSD (Å): {statistics.median(values):.6f}",
        f"Minimum RMSD (Å): {min(values):.6f}",
        f"Maximum RMSD (Å): {max(values):.6f}",
        f"Final RMSD (Å): {values[-1]:.6f}",
        f"Jumps larger than {args.jump_threshold:.3f} Å: {len(jump_indices)}",
    ]
    path.write_text("\n".join(lines) + "\n")


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Fit an AMBER trajectory with one CPPTRAJ atom mask, then calculate "
            "RMSD for another mask on the aligned coordinates."
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
        "-o",
        "--output-prefix",
        default="rmsd",
        help="Output path prefix; default: rmsd",
    )
    parser.add_argument(
        "-r",
        "--reference",
        type=existing_file,
        help="Optional external reference; default: first analyzed trajectory frame",
    )
    parser.add_argument(
        "--fit-mask",
        default=DEFAULT_FIT_MASK,
        help=(
            "CPPTRAJ mask used to align each frame. Prefer an explicit protein "
            "range, e.g. ':23671-23928@N,CA,C,O'."
        ),
    )
    parser.add_argument(
        "--rmsd-mask",
        help=(
            "CPPTRAJ mask used for the reported RMSD after alignment. "
            "Default: the value of --fit-mask."
        ),
    )
    parser.add_argument(
        "--time-per-frame-ps",
        type=positive_float,
        help="Time between consecutive saved frames in each input trajectory, in ps",
    )
    parser.add_argument(
        "--mdin",
        type=existing_file,
        help="AMBER mdin/mdout file from which dt and ntwx will be read",
    )
    parser.add_argument("--start", type=positive_int, default=1)
    parser.add_argument("--stop", type=positive_int)
    parser.add_argument("--stride", type=positive_int, default=1)
    parser.add_argument("--no-autoimage", action="store_true")
    parser.add_argument(
        "--mass-weighted",
        action="store_true",
        help="Use mass weighting for both fitting and measured RMSD",
    )
    parser.add_argument(
        "--smooth-window",
        type=positive_int,
        default=1,
        help="Moving-average window in analyzed frames; raw RMSD is always plotted",
    )
    parser.add_argument(
        "--jump-threshold",
        type=positive_float,
        default=2.0,
        help="Warn when adjacent RMSD values differ by more than this many Å",
    )
    parser.add_argument("--dpi", type=positive_int, default=300)
    parser.add_argument("--cpptraj", default="cpptraj")
    return parser.parse_args()


def main() -> None:
    args = parse_arguments()
    args.fit_mask = validate_mask(args.fit_mask, "--fit-mask")
    if args.rmsd_mask is None:
        args.rmsd_mask = args.fit_mask
    args.rmsd_mask = validate_mask(args.rmsd_mask, "--rmsd-mask")

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

    fit_atom_count = selected_atom_count(
        cpptraj, args.topology, args.fit_mask, "--fit-mask"
    )
    rmsd_atom_count = selected_atom_count(
        cpptraj, args.topology, args.rmsd_mask, "--rmsd-mask"
    )
    if fit_atom_count == 0:
        fail(f"--fit-mask selected zero atoms: {args.fit_mask}")
    if rmsd_atom_count == 0:
        fail(f"--rmsd-mask selected zero atoms: {args.rmsd_mask}")

    print(
        "Fit selection: "
        f"{fit_atom_count if fit_atom_count is not None else 'count unavailable'} atoms"
    )
    print(
        "RMSD selection: "
        f"{rmsd_atom_count if rmsd_atom_count is not None else 'count unavailable'} atoms"
    )

    provenance: List[Tuple[str, int]] = []
    global_saved_indices: List[int] = []
    expected_frames = 0
    cumulative_saved_frames = 0

    for trajectory in args.trajectory:
        total = trajectory_frame_count(cpptraj, args.topology, trajectory)
        local_frames = processed_local_frames(total, args.start, args.stop, args.stride)
        if not local_frames:
            fail(
                f"No frames selected from {trajectory} with start={args.start}, "
                f"stop={args.stop or 'last'}, stride={args.stride}"
            )

        expected_frames += len(local_frames)
        provenance.extend((str(trajectory), frame) for frame in local_frames)
        global_saved_indices.extend(
            cumulative_saved_frames + frame - 1 for frame in local_frames
        )
        cumulative_saved_frames += total
        print(f"{trajectory.name}: {total} total, {len(local_frames)} selected")

    time_per_saved_frame_ps: Optional[float] = args.time_per_frame_ps
    if time_per_saved_frame_ps is None and args.mdin is not None:
        dt_ps, ntwx = parse_amber_timing(args.mdin)
        time_per_saved_frame_ps = dt_ps * ntwx
        print(
            f"Timing from {args.mdin.name}: dt={dt_ps:g} ps, ntwx={ntwx}, "
            f"saved-frame spacing={time_per_saved_frame_ps:g} ps"
        )
    elif time_per_saved_frame_ps is not None and args.mdin is not None:
        warn("--time-per-frame-ps overrides timing parsed from --mdin")

    input_text = build_cpptraj_input(
        topology=args.topology,
        trajectories=args.trajectory,
        raw_output=raw_path,
        fit_mask=args.fit_mask,
        rmsd_mask=args.rmsd_mask,
        reference=args.reference,
        start=args.start,
        stop=args.stop,
        stride=args.stride,
        autoimage=not args.no_autoimage,
        mass_weighted=args.mass_weighted,
    )
    input_path.write_text(input_text)

    result = run_command([cpptraj, "-i", str(input_path)])
    log_path.write_text(result.stdout)
    if result.returncode != 0:
        tail = "\n".join(result.stdout.splitlines()[-30:])
        fail(f"CPPTRAJ failed. Log: {log_path}\n\n{tail}", exit_code=1)

    cpptraj_frames, values = parse_cpptraj_rmsd(raw_path)
    if len(values) != expected_frames:
        warn(
            f"Expected {expected_frames} processed frames but CPPTRAJ returned "
            f"{len(values)}. Check {log_path}."
        )

    if args.reference is None and abs(values[0]) > 1.0e-3:
        warn(
            f"First-frame RMSD is {values[0]:.6f} Å rather than approximately zero. "
            "Check the CPPTRAJ log and masks."
        )

    jump_indices = [
        index
        for index in range(1, len(values))
        if abs(values[index] - values[index - 1]) > args.jump_threshold
    ]
    if jump_indices:
        largest = max(
            abs(values[index] - values[index - 1]) for index in jump_indices
        )
        warn(
            f"Found {len(jump_indices)} adjacent RMSD jump(s) above "
            f"{args.jump_threshold:g} Å; largest={largest:.3f} Å."
        )

    times_ps: Optional[List[float]]
    if time_per_saved_frame_ps is None:
        times_ps = None
        warn(
            "No time spacing supplied. The plot uses processed frame number. "
            "Use --mdin or --time-per-frame-ps for a time axis."
        )
    elif len(global_saved_indices) != len(values):
        warn(
            "Could not map every RMSD value to its source frame; using a "
            "uniform analyzed-frame time axis."
        )
        analyzed_spacing_ps = time_per_saved_frame_ps * args.stride
        times_ps = [index * analyzed_spacing_ps for index in range(len(values))]
    else:
        first_saved_index = global_saved_indices[0]
        times_ps = [
            (index - first_saved_index) * time_per_saved_frame_ps
            for index in global_saved_indices
        ]

    write_csv(csv_path, cpptraj_frames, values, times_ps, provenance)
    reference_label = str(args.reference) if args.reference else "first analyzed frame"
    make_plot(
        plot_path,
        values,
        times_ps,
        reference_label,
        args.smooth_window,
        args.dpi,
    )
    write_summary(
        summary_path,
        args,
        cpptraj,
        version,
        fit_atom_count,
        rmsd_atom_count,
        values,
        time_per_saved_frame_ps,
        jump_indices,
    )

    print("\nCreated:")
    for path in (plot_path, csv_path, summary_path, input_path, log_path, raw_path):
        print(f"  {path}")


if __name__ == "__main__":
    main()
