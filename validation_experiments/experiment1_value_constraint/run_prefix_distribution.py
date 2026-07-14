#!/usr/bin/env python3
"""Experiment 1c: random distribution of all 3-nibble prefixes.

This diagnostic keeps A={12}, B={7}, and prefix length s=3, but increases the
middle core to r_dist=10.  It samples random keys and base states, records the
observed 3-nibble output-difference prefix, and counts all 4096 possible
prefixes.
"""

from __future__ import annotations

import argparse
import math
import random
import sys
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from experiment_utils import (  # noqa: E402
    format_power_of_two,
    format_prefix,
    format_probability,
    log2_or_none,
    write_json,
)
from skinny64 import difference_prefix, random_state  # noqa: E402


DEFAULT_A_CELL = 12
DEFAULT_B_CELL = 7
DEFAULT_R_DIST = 10
DEFAULT_PREFIX_LEN = 3
DEFAULT_TRIALS = 1 << 20
DEFAULT_SEED = 20260714
DEFAULT_TOP = 12
DEFAULT_HISTOGRAM_WIDTH = 16


def prefix_to_code(prefix: tuple[int, ...]) -> int:
    code = 0
    for index, value in enumerate(prefix):
        code |= value << (4 * index)
    return code


def code_to_prefix(code: int, prefix_len: int) -> tuple[int, ...]:
    return tuple((code >> (4 * index)) & 0xF for index in range(prefix_len))


def make_bucket(
    code: int,
    hits: int,
    trials: int,
    prefix_len: int,
) -> dict[str, object]:
    probability = hits / trials
    return {
        "prefix": list(code_to_prefix(code, prefix_len)),
        "hits": hits,
        "observed_probability": probability,
        "observed_probability_log2": log2_or_none(probability),
        "observed_probability_text": format_probability(probability, trials),
    }


def percentile(sorted_values: list[int], q: float) -> int:
    if not sorted_values:
        raise ValueError("cannot take percentile of an empty list")
    index = round(q * (len(sorted_values) - 1))
    return sorted_values[index]


def range_label(value: int, width: int) -> str:
    start = (value // width) * width
    end = start + width - 1
    return f"{start}-{end}"


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Sample the hit-count distribution of all 3-nibble prefixes."
    )
    parser.add_argument("--rounds", type=int, default=DEFAULT_R_DIST)
    parser.add_argument("--active-cell", type=int, default=DEFAULT_A_CELL)
    parser.add_argument("--output-cell", type=int, default=DEFAULT_B_CELL)
    parser.add_argument("--prefix-len", type=int, default=DEFAULT_PREFIX_LEN)
    parser.add_argument("--trials", type=int, default=DEFAULT_TRIALS)
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    parser.add_argument("--top", type=int, default=DEFAULT_TOP)
    parser.add_argument("--histogram-width", type=int, default=DEFAULT_HISTOGRAM_WIDTH)
    parser.add_argument(
        "--output-prefix",
        default="prefix_distribution_r10",
        help="Base filename for JSON and summary outputs.",
    )
    args = parser.parse_args()

    if args.prefix_len <= 0:
        raise ValueError("--prefix-len must be positive")
    if args.active_cell < 0 or args.active_cell >= 16:
        raise ValueError("--active-cell must be in [0, 15]")
    if args.output_cell < 0 or args.output_cell >= 16:
        raise ValueError("--output-cell must be in [0, 15]")
    if args.trials <= 0:
        raise ValueError("--trials must be positive")
    if args.histogram_width <= 0:
        raise ValueError("--histogram-width must be positive")

    prefix_count = 1 << (4 * args.prefix_len)
    counts = [0] * prefix_count
    rng = random.Random(args.seed)

    for _ in range(args.trials):
        master_key = random_state(rng)
        base_state = random_state(rng)
        prefix = difference_prefix(
            master_key,
            base_state,
            args.active_cell,
            args.output_cell,
            args.rounds,
            args.prefix_len,
        )
        counts[prefix_to_code(prefix)] += 1

    expected_probability = 1 / prefix_count
    expected_hits = args.trials * expected_probability
    expected_stddev = math.sqrt(args.trials * expected_probability * (1 - expected_probability))
    hit_values = sorted(counts)
    non_empty_count = sum(1 for hits in counts if hits > 0)
    empty_count = prefix_count - non_empty_count
    empirical_mean = sum(counts) / prefix_count
    empirical_stddev = math.sqrt(
        sum((hits - empirical_mean) ** 2 for hits in counts) / prefix_count
    )
    chi_square = sum((hits - expected_hits) ** 2 / expected_hits for hits in counts)
    reduced_chi_square = chi_square / (prefix_count - 1)

    ranked_codes = sorted(
        range(prefix_count),
        key=lambda code: (counts[code], code),
        reverse=True,
    )
    top_buckets = [
        make_bucket(code, counts[code], args.trials, args.prefix_len)
        for code in ranked_codes[: args.top]
    ]
    bottom_buckets = [
        make_bucket(code, counts[code], args.trials, args.prefix_len)
        for code in reversed(ranked_codes[-args.top :])
    ]

    max_abs_z_code = max(
        range(prefix_count),
        key=lambda code: abs(counts[code] - expected_hits) / expected_stddev,
    )
    selected_prefixes = {}
    for prefix in ((0, 0, 0), (7, 11, 12), (5, 1, 4)):
        if len(prefix) == args.prefix_len:
            code = prefix_to_code(prefix)
            selected_prefixes[format_prefix(prefix)] = make_bucket(
                code, counts[code], args.trials, args.prefix_len
            )

    hit_count_histogram = Counter(counts)
    range_histogram = Counter(
        range_label(hits, args.histogram_width) for hits in counts
    )
    all_buckets = [
        make_bucket(code, counts[code], args.trials, args.prefix_len)
        for code in range(prefix_count)
    ]

    result = {
        "experiment": "experiment1_random_prefix_distribution",
        "cipher": "SKINNY-64-64",
        "r_dist": args.rounds,
        "active_input_cell_A": args.active_cell,
        "output_cell_B": args.output_cell,
        "prefix_length": args.prefix_len,
        "trials": args.trials,
        "trials_text": format_power_of_two(args.trials),
        "seed": args.seed,
        "prefix_count": prefix_count,
        "expected_probability_per_prefix": expected_probability,
        "expected_probability_per_prefix_text": format_power_of_two(expected_probability),
        "expected_hits_per_prefix": expected_hits,
        "expected_hits_per_prefix_text": format_power_of_two(expected_hits),
        "expected_stddev_per_prefix": expected_stddev,
        "non_empty_prefix_count": non_empty_count,
        "empty_prefix_count": empty_count,
        "minimum_hits": hit_values[0],
        "maximum_hits": hit_values[-1],
        "empirical_mean_hits": empirical_mean,
        "empirical_stddev_hits": empirical_stddev,
        "chi_square_against_uniform": chi_square,
        "reduced_chi_square_against_uniform": reduced_chi_square,
        "max_abs_z_prefix": make_bucket(
            max_abs_z_code,
            counts[max_abs_z_code],
            args.trials,
            args.prefix_len,
        ),
        "max_abs_z": abs(counts[max_abs_z_code] - expected_hits) / expected_stddev,
        "hit_count_quantiles": {
            "0%": percentile(hit_values, 0.00),
            "1%": percentile(hit_values, 0.01),
            "5%": percentile(hit_values, 0.05),
            "25%": percentile(hit_values, 0.25),
            "50%": percentile(hit_values, 0.50),
            "75%": percentile(hit_values, 0.75),
            "95%": percentile(hit_values, 0.95),
            "99%": percentile(hit_values, 0.99),
            "100%": percentile(hit_values, 1.00),
        },
        "selected_prefixes": selected_prefixes,
        "top_buckets_by_hits": top_buckets,
        "bottom_buckets_by_hits": bottom_buckets,
        "hit_count_histogram": {
            str(hits): count for hits, count in sorted(hit_count_histogram.items())
        },
        "range_histogram_width": args.histogram_width,
        "hit_count_range_histogram": {
            label: range_histogram[label]
            for label in sorted(
                range_histogram,
                key=lambda text: int(text.split("-", 1)[0]),
            )
        },
        "buckets": all_buckets,
    }

    lines = [
        "Experiment 1c: random prefix distribution",
        (
            "cipher: SKINNY-64-64, "
            f"r_dist={args.rounds}, A={{{args.active_cell}}}, "
            f"B={{{args.output_cell}}}, s={args.prefix_len}"
        ),
        f"trials: {args.trials} = {format_power_of_two(args.trials)}",
        f"prefix buckets: {prefix_count}",
        (
            "expected hits per prefix under uniformity: "
            f"{expected_hits:.2f} = {format_power_of_two(expected_hits)}"
        ),
        f"expected probability per prefix: {format_power_of_two(expected_probability)}",
        "",
        "aggregate distribution:",
        f"  non-empty prefixes: {non_empty_count}/{prefix_count}",
        f"  empty prefixes: {empty_count}/{prefix_count}",
        f"  min hits: {hit_values[0]}",
        f"  max hits: {hit_values[-1]}",
        f"  empirical mean hits: {empirical_mean:.2f}",
        f"  empirical stddev hits: {empirical_stddev:.2f}",
        f"  expected per-prefix stddev: {expected_stddev:.2f}",
        f"  reduced chi-square vs uniform: {reduced_chi_square:.4f}",
        (
            "  largest absolute z-score: "
            f"{result['max_abs_z']:.3f} at "
            f"{format_prefix(tuple(result['max_abs_z_prefix']['prefix']))}"
        ),
        "",
        "hit-count quantiles:",
    ]
    for label, value in result["hit_count_quantiles"].items():
        lines.append(f"  {label}: {value}")

    lines.extend(["", f"hit-count range histogram (width={args.histogram_width}):"])
    for label, count in result["hit_count_range_histogram"].items():
        lines.append(f"  {label}: {count} prefixes")

    lines.extend(["", "selected prefixes:"])
    for prefix_text, bucket in selected_prefixes.items():
        lines.append(
            f"  {prefix_text}: hits={bucket['hits']}, "
            f"probability={bucket['observed_probability_text']}"
        )

    lines.extend(["", f"top {args.top} prefixes by hits:"])
    for bucket in top_buckets:
        lines.append(
            f"  {format_prefix(tuple(bucket['prefix']))}: "
            f"hits={bucket['hits']}, probability={bucket['observed_probability_text']}"
        )

    lines.extend(["", f"bottom {args.top} prefixes by hits:"])
    for bucket in bottom_buckets:
        lines.append(
            f"  {format_prefix(tuple(bucket['prefix']))}: "
            f"hits={bucket['hits']}, probability={bucket['observed_probability_text']}"
        )

    summary = "\n".join(lines).rstrip() + "\n"
    out_dir = Path(__file__).resolve().parent
    write_json(out_dir / f"{args.output_prefix}.json", result)
    (out_dir / f"{args.output_prefix}_summary.txt").write_text(summary, encoding="utf-8")
    print(summary)


if __name__ == "__main__":
    main()
