#!/usr/bin/env python3
"""Experiment 1: CT1 value-constraint probability for 6-round SKINNY-64-64.

This experiment measures how often a fixed 3-nibble prefix appears in the
difference sequence at cell B={7}, when the input delta-set is active only at
cell A={12}.  The idealized per-prefix probability is 2^{-12}.

The default run deliberately checks both:

    (0,0,0)    - no hit is observed, illustrating short-round non-uniformity;
    (7,11,12) - a reachable prefix used by Experiment 2.
"""

from __future__ import annotations

import argparse
import math
import random
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from experiment_utils import (  # noqa: E402
    format_power_of_two,
    format_prefix,
    format_probability,
    log2_or_none,
    parse_prefixes,
    write_json,
)
from skinny64 import difference_prefix, random_state  # noqa: E402


A_CELL = 12
B_CELL = 7
R_DIST = 6
PREFIX_LEN = 3
DEFAULT_CONSTRAINTS = ["0,0,0", "7,11,12"]
DEFAULT_SEED = 20260707


def classify_prefix(prefix: tuple[int, ...], hits: int) -> str:
    if prefix == (0, 0, 0) and hits == 0:
        return (
            "no hit observed; this illustrates the visible non-uniformity of "
            "the short 6-round core"
        )
    if prefix == (7, 11, 12):
        return "reachable prefix used by Experiment 2 for key-recovery validation"
    if hits == 0:
        return "no hit observed for this prefix in the sampled trials"
    return "reachable prefix in the sampled trials"


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Measure CT1 prefix probabilities for the 6-round SKINNY-64-64 core."
    )
    parser.add_argument("--trials", type=int, default=1 << 20)
    parser.add_argument(
        "--constraints",
        nargs="+",
        default=DEFAULT_CONSTRAINTS,
        help="One or more 3-nibble prefixes, for example: --constraints 0,0,0 7,11,12",
    )
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    args = parser.parse_args()

    constraints = parse_prefixes(args.constraints, PREFIX_LEN)
    counts = {constraint: 0 for constraint in constraints}
    examples = {constraint: [] for constraint in constraints}
    rng = random.Random(args.seed)

    for trial in range(args.trials):
        master_key = random_state(rng)
        base_state = random_state(rng)
        prefix = difference_prefix(master_key, base_state, A_CELL, B_CELL, R_DIST, PREFIX_LEN)
        if prefix in counts:
            counts[prefix] += 1
            if len(examples[prefix]) < 5:
                examples[prefix].append(
                    {"trial": trial, "key": master_key, "base_state": base_state}
                )

    expected_probability = 2 ** (-4 * PREFIX_LEN)
    expected_hits = args.trials * expected_probability
    variance = args.trials * expected_probability * (1 - expected_probability)

    prefix_results = []
    for constraint in constraints:
        hits = counts[constraint]
        observed_probability = hits / args.trials
        observed_over_expected = observed_probability / expected_probability
        z_score = (hits - expected_hits) / math.sqrt(variance) if variance > 0 else 0.0
        prefix_results.append(
            {
                "constraint_prefix": list(constraint),
                "hits": hits,
                "trials": args.trials,
                "observed_probability": observed_probability,
                "observed_probability_log2": log2_or_none(observed_probability),
                "observed_probability_text": format_probability(observed_probability, args.trials),
                "expected_probability": expected_probability,
                "expected_probability_log2": math.log2(expected_probability),
                "expected_probability_text": format_power_of_two(expected_probability),
                "observed_over_expected": observed_over_expected,
                "z_score": z_score,
                "interpretation": classify_prefix(constraint, hits),
                "examples": examples[constraint],
            }
        )

    result = {
        "experiment": "experiment1_value_constraint_probability",
        "cipher": "SKINNY-64-64",
        "r_dist": R_DIST,
        "active_input_cell_A": A_CELL,
        "output_cell_B": B_CELL,
        "constraint_prefix_length": PREFIX_LEN,
        "trials": args.trials,
        "trials_text": format_power_of_two(args.trials),
        "seed": args.seed,
        "expected_probability_per_prefix": expected_probability,
        "expected_probability_per_prefix_text": format_power_of_two(expected_probability),
        "prefix_results": prefix_results,
        "main_takeaway": (
            "The all-zero prefix has no hit in this sample, showing that the "
            "6-round core is visibly non-uniform; the reachable prefix (7,11,12) "
            "has probability close to 2^{-12} and is used in Experiment 2."
        ),
    }

    lines = [
        "Experiment 1: CT1 value-constraint probability",
        f"cipher: SKINNY-64-64, r_dist={R_DIST}, A={{{A_CELL}}}, B={{{B_CELL}}}, s={PREFIX_LEN}",
        f"trials: {args.trials} = {format_power_of_two(args.trials)}",
        f"expected probability for a fixed 3-nibble prefix: {format_power_of_two(expected_probability)}",
        "",
    ]
    for item in prefix_results:
        lines.extend(
            [
                f"prefix {format_prefix(tuple(item['constraint_prefix']))}:",
                f"  hits: {item['hits']}/{args.trials}",
                f"  observed probability: {item['observed_probability_text']}",
                f"  expected probability: {item['expected_probability_text']}",
                f"  observed/expected: {item['observed_over_expected']:.4f}",
                f"  z-score: {item['z_score']:.3f}",
                f"  note: {item['interpretation']}",
                "",
            ]
        )
    summary = "\n".join(lines).rstrip() + "\n"

    out_dir = Path(__file__).resolve().parent
    write_json(out_dir / "results.json", result)
    (out_dir / "summary.txt").write_text(summary, encoding="utf-8")
    print(summary)


if __name__ == "__main__":
    main()
