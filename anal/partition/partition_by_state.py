#!/usr/bin/env python3
"""
partition_by_state.py

Generic post-processing layer for state-conditioning frame-wise observables.

The program does not calculate structural observables and does not interpret
CpHMD lambda coordinates. It joins two already-defined data products:

  1. a frame-wise observation table (for example, calc_distance.py output), and
  2. a frame/state assignment table.

For the current HEWL coupled dyad, dyad_frame_assignments.tsv can be used
directly with the default join key:

    trajectory,trajectory_frame

Shared frame/provenance metadata are validated rather than duplicated in the
joined output. In particular:

  - analysis_row, trajectory_frame, and step must agree when populated in both
    input tables;
  - pH/ph values are treated as equivalent metadata, validated numerically,
    and canonicalized to a single output column named "pH";
  - state-side duplicate frame identifiers are not copied into the joined TSV.

All state rows are retained during the exact frame join. State-quality filters
(clean/matched by default) are applied only when calculating summaries, so
mixed and outside-window frames remain auditable and retain their positions on
the original trajectory timeline.

For block SEM, contiguous blocks are constructed on the complete joined
trajectory timeline first. Eligible frames belonging to each requested state
are then selected within each block. Intermittent state visits are therefore
not concatenated before blocking.

For future independent residue protonation, a long-format state table can be
joined with, for example:

    --join-on trajectory,trajectory_frame,resid

without changing the partitioning/statistics code.
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

import numpy as np


FRAME_ID_COLUMNS = {
    "trajectory",
    "analysis_row",
    "trajectory_frame",
    "step",
    "absolute_time_ns",
}


def read_tsv(path: Path) -> list[dict]:
    with path.open(newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
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
                if isinstance(value, (float, np.floating)):
                    out[key] = (
                        ""
                        if not np.isfinite(value)
                        else f"{float(value):.10g}"
                    )
            writer.writerow(out)


def truthy(value: str | None) -> bool:
    if value is None:
        return False
    return str(value).strip().lower() in {"1", "true", "yes", "y"}


def parse_float(value: str | None) -> float:
    if value is None or not str(value).strip():
        return math.nan

    try:
        return float(value)
    except ValueError:
        return math.nan


def parse_int_exact(value: str | None, context: str) -> int:
    if value is None or not str(value).strip():
        raise RuntimeError(f"Missing integer {context}")

    try:
        x = float(value)
    except ValueError as exc:
        raise RuntimeError(
            f"Non-numeric integer {context}: {value!r}"
        ) from exc

    if not math.isfinite(x):
        raise RuntimeError(f"Non-finite integer {context}: {value!r}")

    i = int(round(x))
    if abs(x - i) > 1.0e-8:
        raise RuntimeError(f"Non-integral {context}: {value!r}")

    return i


def nonempty(value: str | None) -> bool:
    return value is not None and bool(str(value).strip())


def first_nonempty(
    row: dict,
    names: tuple[str, ...],
) -> tuple[str | None, str | None]:
    for name in names:
        if name in row and nonempty(row.get(name)):
            return name, str(row[name]).strip()

    return None, None


def validate_equal_int_metadata(
    obs_row: dict,
    state_row: dict,
    column: str,
    key: tuple,
) -> None:
    """
    Validate a shared integer-valued metadata field when populated in both rows.

    Missing metadata on one side is allowed because some observation/state
    formats may not carry every provenance field.
    """
    obs_value = obs_row.get(column)
    state_value = state_row.get(column)

    if not (nonempty(obs_value) and nonempty(state_value)):
        return

    obs_int = parse_int_exact(obs_value, f"observation.{column}")
    state_int = parse_int_exact(state_value, f"state.{column}")

    if obs_int != state_int:
        raise RuntimeError(
            f"{column} mismatch for join key {key}: "
            f"observations={obs_int}, states={state_int}"
        )


def canonical_ph(
    obs_row: dict,
    state_row: dict,
    key: tuple,
) -> str:
    """
    Validate pH/ph metadata and return one canonical string representation.

    "pH" and "ph" are treated as equivalent input column names. If both input
    tables provide a value, they must agree numerically.
    """
    _, obs_value = first_nonempty(obs_row, ("pH", "ph"))
    _, state_value = first_nonempty(state_row, ("pH", "ph"))

    if obs_value is not None and state_value is not None:
        try:
            obs_ph = float(obs_value)
            state_ph = float(state_value)
        except ValueError as exc:
            raise RuntimeError(
                f"Non-numeric pH metadata for join key {key}: "
                f"observations={obs_value!r}, states={state_value!r}"
            ) from exc

        if not (math.isfinite(obs_ph) and math.isfinite(state_ph)):
            raise RuntimeError(
                f"Non-finite pH metadata for join key {key}: "
                f"observations={obs_value!r}, states={state_value!r}"
            )

        if not math.isclose(
            obs_ph,
            state_ph,
            rel_tol=0.0,
            abs_tol=1.0e-8,
        ):
            raise RuntimeError(
                f"pH mismatch for join key {key}: "
                f"observations={obs_ph:g}, states={state_ph:g}"
            )

    value = obs_value if obs_value is not None else state_value
    return "" if value is None else value


def parse_csv_names(text: str) -> list[str]:
    values = [x.strip() for x in text.split(",") if x.strip()]
    if not values:
        raise ValueError(
            "Expected at least one comma-separated column name."
        )

    if len(values) != len(set(values)):
        raise ValueError(
            f"Duplicate column name in comma-separated list: {text!r}"
        )

    return values


def make_key(
    row: dict,
    columns: list[str],
    table_name: str,
) -> tuple:
    values = []

    for column in columns:
        if column not in row:
            raise RuntimeError(
                f"{table_name} is missing join column {column!r}"
            )

        value = row[column]

        if column in {
            "analysis_row",
            "trajectory_frame",
            "step",
            "resid",
        }:
            if value is None or not str(value).strip():
                values.append("")
            else:
                values.append(
                    parse_int_exact(
                        value,
                        f"{table_name}.{column}",
                    )
                )
        else:
            values.append(str(value).strip())

    return tuple(values)


def summarize_state_masked(
    values: np.ndarray,
    state_labels: np.ndarray,
    eligible_mask: np.ndarray,
    target_state: str,
    timeline_rows: np.ndarray,
    nblocks: int,
) -> dict:
    """
    Summarize one state while preserving the complete original trajectory timeline.

    Mean, SD, and quantiles use frames that are both:
      1. eligible under clean/matched filtering, and
      2. assigned to target_state.

    For block SEM, contiguous blocks are formed over ALL joined trajectory
    frames first. Only then are eligible target-state frames selected within
    each block. Mixed, unmatched, and outside-window frames therefore retain
    their temporal positions and never cause artificial timeline compression.

    block_sem is a blocking diagnostic and should not be interpreted as proof
    that the retained state samples are statistically independent.
    """
    values = np.asarray(values, dtype=float)
    state_labels = np.asarray(state_labels, dtype=object)
    eligible_mask = np.asarray(eligible_mask, dtype=bool)
    timeline_rows = np.asarray(timeline_rows, dtype=int)

    if not (
        len(values)
        == len(state_labels)
        == len(eligible_mask)
        == len(timeline_rows)
    ):
        raise RuntimeError(
            "Internal length mismatch in state summary."
        )

    target_mask = (
        eligible_mask
        & (state_labels == target_state)
        & np.isfinite(values)
    )

    x = values[target_mask]
    n = int(x.size)

    if n == 0:
        return {
            "nframes": 0,
            "mean": np.nan,
            "sd": np.nan,
            "block_sem": np.nan,
            "nblocks": 0,
            "median": np.nan,
            "q10": np.nan,
            "q90": np.nan,
            "min": np.nan,
            "max": np.nan,
        }

    # Construct blocks on the COMPLETE joined timeline.
    order = np.argsort(timeline_rows)
    values_ordered = values[order]
    states_ordered = state_labels[order]
    eligible_ordered = eligible_mask[order]

    nb = min(max(1, nblocks), len(values_ordered))
    index_blocks = [
        block
        for block in np.array_split(
            np.arange(len(values_ordered)),
            nb,
        )
        if len(block)
    ]

    block_means = []

    for indices in index_blocks:
        block_mask = (
            eligible_ordered[indices]
            & (states_ordered[indices] == target_state)
            & np.isfinite(values_ordered[indices])
        )

        if not np.any(block_mask):
            continue

        block_values = values_ordered[indices][block_mask]
        if len(block_values):
            block_means.append(
                float(np.mean(block_values))
            )

    if len(block_means) > 1:
        block_sem = float(
            np.std(
                np.asarray(block_means),
                ddof=1,
            )
            / np.sqrt(len(block_means))
        )
    else:
        block_sem = np.nan

    return {
        "nframes": n,
        "mean": float(np.mean(x)),
        "sd": (
            float(np.std(x, ddof=1))
            if n > 1
            else np.nan
        ),
        "block_sem": block_sem,
        "nblocks": len(block_means),
        "median": float(np.median(x)),
        "q10": float(np.quantile(x, 0.10)),
        "q90": float(np.quantile(x, 0.90)),
        "min": float(np.min(x)),
        "max": float(np.max(x)),
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Join frame-wise observables to state assignments "
            "and summarize by state."
        )
    )

    parser.add_argument(
        "--observations",
        type=Path,
        required=True,
        help="Frame-wise observable TSV.",
    )
    parser.add_argument(
        "--states",
        type=Path,
        required=True,
        help="Authoritative frame/state assignment TSV.",
    )
    parser.add_argument(
        "--join-on",
        default="trajectory,trajectory_frame",
        help=(
            "Comma-separated exact join key. Current coupled-dyad "
            "default: trajectory,trajectory_frame. Future per-residue "
            "states can use trajectory,trajectory_frame,resid."
        ),
    )
    parser.add_argument(
        "--value-column",
        default="value",
        help="Observation column containing the numeric value.",
    )
    parser.add_argument(
        "--state-column",
        default="state_label",
        help="State-table column containing the state label.",
    )
    parser.add_argument(
        "--clean-column",
        default="clean",
        help="Optional state-table clean-state flag column.",
    )
    parser.add_argument(
        "--matched-column",
        default="matched",
        help="Optional state-table exact-match flag column.",
    )
    parser.add_argument(
        "--include-unclean",
        action="store_true",
        help=(
            "Include rows marked unclean in state-conditioned summaries. "
            "They remain present in the joined audit table regardless."
        ),
    )
    parser.add_argument(
        "--include-unmatched",
        action="store_true",
        help=(
            "Include rows marked unmatched in state-conditioned summaries. "
            "They remain present in the joined audit table regardless."
        ),
    )
    parser.add_argument(
        "--group-by",
        default=None,
        help=(
            "Comma-separated descriptor columns for summaries. "
            "If omitted, observation columns except frame IDs and "
            "the value column are used."
        ),
    )
    parser.add_argument(
        "--blocks",
        type=int,
        default=5,
        help=(
            "Number of contiguous blocks constructed on each complete "
            "trajectory timeline for block-SEM diagnostics."
        ),
    )
    parser.add_argument(
        "--joined-output",
        type=Path,
        default=Path("observations_with_states.tsv"),
        help="Per-frame joined audit TSV.",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=Path("state_partition_summary.tsv"),
        help="State-conditioned summary TSV.",
    )

    return parser.parse_args()


def main() -> None:
    args = parse_args()

    if args.blocks < 1:
        raise ValueError("--blocks must be >= 1")

    observations = read_tsv(args.observations)
    states = read_tsv(args.states)

    if not observations:
        raise RuntimeError(
            f"Observation table is empty: {args.observations}"
        )

    if not states:
        raise RuntimeError(
            f"State table is empty: {args.states}"
        )

    observation_columns = list(observations[0].keys())
    state_columns = list(states[0].keys())

    if args.value_column not in observation_columns:
        raise RuntimeError(
            f"Observation table lacks value column "
            f"{args.value_column!r}"
        )

    if args.state_column not in state_columns:
        raise RuntimeError(
            f"State table lacks state column "
            f"{args.state_column!r}"
        )

    if args.state_column in observation_columns:
        raise RuntimeError(
            f"Observation table already contains the state column "
            f"{args.state_column!r}; state annotations must come "
            f"from --states."
        )

    join_columns = parse_csv_names(args.join_on)

    # Build the lookup from ALL state rows. clean/matched filtering is a
    # summary decision, not a frame-identity decision.
    state_lookup: dict[tuple, dict] = {}

    for row in states:
        key = make_key(
            row,
            join_columns,
            "state table",
        )

        if key in state_lookup:
            raise RuntimeError(
                f"State join key is not unique: {key}. "
                f"Use a more specific --join-on key."
            )

        state_lookup[key] = row

    if not state_lookup:
        raise RuntimeError(
            "State table contains no joinable rows."
        )

    # Shared frame/provenance metadata are validated rather than copied
    # twice into the joined output.
    canonical_metadata = {
        "trajectory",
        "analysis_row",
        "trajectory_frame",
        "step",
        "pH",
        "ph",
    }

    state_columns_to_copy = [
        column
        for column in state_columns
        if (
            column not in join_columns
            and column not in canonical_metadata
        )
    ]

    joined: list[dict] = []
    missing: list[tuple] = []

    for observation in observations:
        key = make_key(
            observation,
            join_columns,
            "observation table",
        )

        state = state_lookup.get(key)

        if state is None:
            missing.append(key)
            continue

        # Validate shared frame/provenance metadata before merging.
        for column in (
            "analysis_row",
            "trajectory_frame",
            "step",
        ):
            validate_equal_int_metadata(
                observation,
                state,
                column,
                key,
            )

        ph_value = canonical_ph(
            observation,
            state,
            key,
        )

        # Start from observation columns, but canonicalize pH/ph to "pH".
        out = {
            column: value
            for column, value in observation.items()
            if column not in {"pH", "ph"}
        }

        # Keep canonical pH near trajectory in the joined output.
        if (
            ph_value != ""
            or "pH" in observation
            or "ph" in observation
            or "pH" in state
            or "ph" in state
        ):
            reordered: dict = {}

            for column, value in out.items():
                reordered[column] = value

                if column == "trajectory":
                    reordered["pH"] = ph_value

            if "pH" not in reordered:
                reordered["pH"] = ph_value

            out = reordered

        # Append only genuinely state-specific information.
        for column in state_columns_to_copy:
            destination = (
                column
                if column not in out
                else f"state_{column}"
            )
            out[destination] = state.get(column, "")

        joined.append(out)

    if missing:
        examples = "; ".join(
            map(str, missing[:10])
        )
        raise RuntimeError(
            f"{len(missing)} observation rows have no state "
            f"assignment. Examples: {examples}"
        )

    if not joined:
        raise RuntimeError(
            "No observation rows were joined to states."
        )

    joined_fields = list(joined[0].keys())
    write_tsv(
        args.joined_output,
        joined,
        joined_fields,
    )

    # Apply state-quality filters only for statistical summaries.
    eligible: list[dict] = []
    excluded_unclean = 0
    excluded_unmatched = 0

    for row in joined:
        matched_column = (
            args.matched_column
            if args.matched_column in row
            else f"state_{args.matched_column}"
        )
        clean_column = (
            args.clean_column
            if args.clean_column in row
            else f"state_{args.clean_column}"
        )

        if (
            not args.include_unmatched
            and matched_column in row
            and not truthy(row.get(matched_column))
        ):
            excluded_unmatched += 1
            continue

        if (
            not args.include_unclean
            and clean_column in row
            and not truthy(row.get(clean_column))
        ):
            excluded_unclean += 1
            continue

        eligible.append(row)

    if not eligible:
        raise RuntimeError(
            "All joined observation rows were excluded from "
            "summaries by clean/matched filtering."
        )

    if args.group_by:
        group_columns = parse_csv_names(
            args.group_by
        )
    else:
        excluded = set(FRAME_ID_COLUMNS)
        excluded.add(args.value_column)

        # Exclude copied state-specific columns from the default observable
        # descriptor key.
        excluded.update(
            (
                column
                if column not in observation_columns
                else f"state_{column}"
            )
            for column in state_columns_to_copy
        )

        group_columns = [
            column
            for column in observation_columns
            if (
                column not in excluded
                and column != "ph"
            )
        ]

        # If the observation table used "ph", the joined table now uses "pH".
        if (
            "ph" in observation_columns
            and "pH" not in group_columns
        ):
            group_columns.append("pH")

    # Validate requested group columns against the joined schema.
    missing_group_columns = [
        column
        for column in group_columns
        if column not in joined_fields
    ]

    if missing_group_columns:
        raise RuntimeError(
            "Joined table lacks requested --group-by column(s): "
            + ", ".join(repr(x) for x in missing_group_columns)
        )

    # Always summarize each original trajectory independently. This
    # preserves a well-defined contiguous timeline and avoids accidental
    # blocking across pH slots or replica trajectories.
    if "trajectory" not in group_columns:
        group_columns = ["trajectory"] + group_columns

    # Group using ALL joined rows so original temporal spacing is retained.
    groups: dict[tuple, list[dict]] = {}

    for row in joined:
        key = tuple(
            row.get(column, "")
            for column in group_columns
        )
        groups.setdefault(key, []).append(row)

    summary_rows: list[dict] = []

    for key, rows in groups.items():
        values = np.asarray(
            [
                parse_float(
                    row.get(args.value_column)
                )
                for row in rows
            ],
            dtype=float,
        )

        state_labels = np.asarray(
            [
                str(
                    row.get(
                        args.state_column,
                        "",
                    )
                ).strip()
                for row in rows
            ],
            dtype=object,
        )

        timeline_rows = np.asarray(
            [
                parse_int_exact(
                    row.get("analysis_row")
                    or row.get("trajectory_frame"),
                    "analysis_row/trajectory_frame",
                )
                for row in rows
            ],
            dtype=int,
        )

        eligible_mask = []

        for row in rows:
            matched_column = (
                args.matched_column
                if args.matched_column in row
                else f"state_{args.matched_column}"
            )
            clean_column = (
                args.clean_column
                if args.clean_column in row
                else f"state_{args.clean_column}"
            )

            keep = True

            if (
                not args.include_unmatched
                and matched_column in row
                and not truthy(
                    row.get(matched_column)
                )
            ):
                keep = False

            if (
                not args.include_unclean
                and clean_column in row
                and not truthy(
                    row.get(clean_column)
                )
            ):
                keep = False

            eligible_mask.append(keep)

        eligible_mask = np.asarray(
            eligible_mask,
            dtype=bool,
        )

        unique_states = sorted(
            {
                state
                for state, keep in zip(
                    state_labels,
                    eligible_mask,
                )
                if keep and state
            }
        )

        descriptor = dict(
            zip(
                group_columns,
                key,
            )
        )

        for state in unique_states:
            statistics = summarize_state_masked(
                values=values,
                state_labels=state_labels,
                eligible_mask=eligible_mask,
                target_state=state,
                timeline_rows=timeline_rows,
                nblocks=args.blocks,
            )

            summary_rows.append(
                {
                    **descriptor,
                    "state": state,
                    **statistics,
                }
            )

    fields = (
        group_columns
        + [
            "state",
            "nframes",
            "mean",
            "sd",
            "block_sem",
            "nblocks",
            "median",
            "q10",
            "q90",
            "min",
            "max",
        ]
    )

    write_tsv(
        args.output,
        summary_rows,
        fields,
    )

    print(
        f"Joined rows              : {len(joined)}"
    )
    print(
        f"Rows used in summaries   : {len(eligible)}"
    )
    print(
        f"Excluded unmatched rows  : "
        f"{excluded_unmatched}"
    )
    print(
        f"Excluded unclean rows    : "
        f"{excluded_unclean}"
    )
    print(
        f"Wrote {args.joined_output}"
    )
    print(
        f"Wrote {args.output}"
    )


if __name__ == "__main__":
    main()
