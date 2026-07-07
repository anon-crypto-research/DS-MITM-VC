#!/usr/bin/env python3
"""Experiment 2: reduced 0+6+2 SKINNY-64-64 filtering validation.

The MILP output in ``expiriment/results/output`` gives the reduced case:

    r_in = 0, r_dist = 6, r_out = 2
    A = {12}, B = {7}, Val_Con = 3, Key = Kg = 1.

This script fixes the reachable CT1 prefix (7,11,12), runs the corresponding
8-round filtering procedure, and checks whether the model-indicated online
master-key nibble is recovered.
"""

from __future__ import annotations

import argparse
import random
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from experiment_utils import (  # noqa: E402
    format_power_of_two,
    format_prefix,
    log2_or_none,
    parse_prefix,
    write_json,
)
from skinny64 import (  # noqa: E402
    check_inverse,
    decrypt_to_round,
    delta_states,
    difference_prefix,
    difference_sequence,
    encrypt_with_trace,
    random_state,
)


A_CELL = 12
B_CELL = 7
R_DIST = 6
R_OUT = 2
TOTAL_ROUNDS = R_DIST + R_OUT
PREFIX_LEN = 3
DEFAULT_CONSTRAINT = "7,11,12"
DEFAULT_SEED = 20260707
GUESSED_MASTER_CELL = 15
MODEL_DATA_EXPONENT = 16
MODEL_ONLINE_EXPONENT = 16


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Run the reduced 0+6+2 SKINNY-64-64 filtering validation."
    )
    parser.add_argument("--valid-trials", type=int, default=64)
    parser.add_argument("--max-attempts", type=int, default=1 << 22)
    parser.add_argument(
        "--constraint",
        default=DEFAULT_CONSTRAINT,
        help="Fixed reachable 3-nibble CT1 prefix. Default: 7,11,12.",
    )
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    parser.add_argument("--self-test-trials", type=int, default=100)
    args = parser.parse_args()

    check_inverse(TOTAL_ROUNDS, trials=args.self_test_trials)

    constraint = parse_prefix(args.constraint, PREFIX_LEN)
    rng = random.Random(args.seed)

    attempts = 0
    valid = 0
    correct_retained = 0
    wrong_survivors = 0
    survivor_histogram: dict[str, int] = {}
    examples = []

    while valid < args.valid_trials and attempts < args.max_attempts:
        attempts += 1
        master_key = random_state(rng)
        base_state = random_state(rng)
        if difference_prefix(master_key, base_state, A_CELL, B_CELL, R_DIST, PREFIX_LEN) != constraint:
            continue

        plaintexts = delta_states(base_state, A_CELL)
        traces = [encrypt_with_trace(plaintext, master_key, TOTAL_ROUNDS) for plaintext in plaintexts]
        true_x6 = [trace[R_DIST] for trace in traces]
        ciphertexts = [trace[TOTAL_ROUNDS] for trace in traces]
        true_sequence = difference_sequence(true_x6, B_CELL)

        if true_sequence[:PREFIX_LEN] != constraint:
            continue

        valid += 1
        h_suffix = true_sequence[PREFIX_LEN:]
        true_guess = master_key[GUESSED_MASTER_CELL]
        survivors = []

        for guess in range(16):
            guessed_key = master_key[:]
            guessed_key[GUESSED_MASTER_CELL] = guess
            guessed_x6 = [
                decrypt_to_round(ciphertext, guessed_key, R_DIST, TOTAL_ROUNDS)
                for ciphertext in ciphertexts
            ]
            guessed_sequence = difference_sequence(guessed_x6, B_CELL)
            if guessed_sequence[:PREFIX_LEN] == constraint and guessed_sequence[PREFIX_LEN:] == h_suffix:
                survivors.append(guess)

        if true_guess in survivors:
            correct_retained += 1
        wrong_survivors += sum(1 for guess in survivors if guess != true_guess)
        survivor_histogram[str(len(survivors))] = survivor_histogram.get(str(len(survivors)), 0) + 1

        if len(examples) < 5:
            examples.append(
                {
                    "attempt": attempts,
                    "true_guess": true_guess,
                    "survivors": survivors,
                    "key": master_key,
                    "base_state": base_state,
                    "sequence_prefix": list(true_sequence[:PREFIX_LEN]),
                    "sequence_suffix": list(h_suffix),
                }
            )

    if valid == 0:
        raise RuntimeError("No constrained delta-set found; increase --max-attempts.")

    empirical_constraint_probability = valid / attempts
    expected_constraint_probability = 2 ** (-4 * PREFIX_LEN)
    attempts_per_valid = attempts / valid
    average_data_search_cost = 16 * attempts_per_valid

    result = {
        "experiment": "experiment2_reduced_attack_filtering",
        "cipher": "SKINNY-64-64",
        "round_split": "0+6+2",
        "active_input_cell_A": A_CELL,
        "output_cell_B": B_CELL,
        "constraint_prefix": list(constraint),
        "constraint_prefix_text": format_prefix(constraint),
        "guessed_master_cell": GUESSED_MASTER_CELL,
        "requested_valid_trials": args.valid_trials,
        "valid_trials": valid,
        "attempts": attempts,
        "seed": args.seed,
        "empirical_constraint_probability": empirical_constraint_probability,
        "empirical_constraint_probability_log2": log2_or_none(empirical_constraint_probability),
        "empirical_constraint_probability_text": format_power_of_two(empirical_constraint_probability),
        "expected_constraint_probability": expected_constraint_probability,
        "expected_constraint_probability_log2": -4 * PREFIX_LEN,
        "expected_constraint_probability_text": format_power_of_two(expected_constraint_probability),
        "attempts_per_valid_delta_set": attempts_per_valid,
        "attempts_per_valid_delta_set_log2": log2_or_none(attempts_per_valid),
        "attempts_per_valid_delta_set_text": format_power_of_two(attempts_per_valid),
        "plaintexts_per_delta_set": 16,
        "plaintexts_per_delta_set_text": "2^{4}",
        "average_data_search_cost": average_data_search_cost,
        "average_data_search_cost_log2": log2_or_none(average_data_search_cost),
        "average_data_search_cost_text": format_power_of_two(average_data_search_cost),
        "model_prediction": {
            "data_complexity_text": f"2^{{{MODEL_DATA_EXPONENT}}}",
            "online_time_complexity_text": f"2^{{{MODEL_ONLINE_EXPONENT}}}",
        },
        "correct_retained": correct_retained,
        "correct_retention_rate": correct_retained / valid,
        "wrong_survivors": wrong_survivors,
        "survivor_histogram": survivor_histogram,
        "examples": examples,
    }

    summary = (
        "Experiment 2: reduced 0+6+2 filtering validation\n"
        f"cipher: SKINNY-64-64, A={{{A_CELL}}}, B={{{B_CELL}}}, s={PREFIX_LEN}\n"
        f"fixed reachable CT1 prefix: {format_prefix(constraint)}\n"
        f"model prediction: Data=2^{{{MODEL_DATA_EXPONENT}}}, Online=2^{{{MODEL_ONLINE_EXPONENT}}}\n"
        f"valid constrained delta-sets: {valid} found in {attempts} attempts\n"
        f"empirical constraint probability: {format_power_of_two(empirical_constraint_probability)} "
        f"({empirical_constraint_probability:.8e})\n"
        f"expected CT1 probability: {format_power_of_two(expected_constraint_probability)} "
        f"({expected_constraint_probability:.8e})\n"
        f"attempts per valid delta-set: {attempts_per_valid:.2f} = "
        f"{format_power_of_two(attempts_per_valid)}\n"
        f"average data/online-search cost: {average_data_search_cost:.2f} = "
        f"{format_power_of_two(average_data_search_cost)}\n"
        f"correct key-nibble retained: {correct_retained}/{valid}\n"
        f"wrong surviving key-nibble guesses: {wrong_survivors}\n"
        f"survivor histogram: {survivor_histogram}\n"
    )

    out_dir = Path(__file__).resolve().parent
    write_json(out_dir / "results.json", result)
    (out_dir / "summary.txt").write_text(summary, encoding="utf-8")
    print(summary)


if __name__ == "__main__":
    main()
