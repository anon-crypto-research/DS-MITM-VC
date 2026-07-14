#!/usr/bin/env python3
"""Experiment 2: reduced 0+6+2 SKINNY-64-64 key-recovery validation.

The MILP output in ``validation_experiments/results/output`` gives the reduced case:

    r_in = 0, r_dist = 6, r_out = 2
    A = {12}, B = {7}, Val_Con = 3, Key = Kg = 1.

This script simulates the reduced key-recovery loop: an oracle samples a master
key, the attacker asks for 16 ciphertexts from one A={12} delta-set at a time,
enumerates the model-indicated key nibble, and checks the resulting 6-round
ordered sequence against the offline table for a fixed CT1 prefix.
"""

from __future__ import annotations

import argparse
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
    log2_or_none,
    parse_prefix,
    write_json,
)
from skinny64 import (  # noqa: E402
    INV_SBOX,
    ROUND_CONSTANTS,
    check_inverse,
    delta_states,
    difference_sequence,
    encrypt_rounds,
    random_state,
    round_keys,
)

from experiment1_value_constraint.run_sequence_distribution import (  # noqa: E402
    IDENTITY_SEQUENCE_CODE,
    R_DIST as TABLE_R_DIST,
    apply_transform,
    build_diff_permutations,
    build_transform_tables,
    decode_prefix,
)


A_CELL = 12
B_CELL = 7
R_DIST = 6
R_OUT = 2
TOTAL_ROUNDS = R_DIST + R_OUT
PREFIX_LEN = 3
DEFAULT_CONSTRAINTS = ["5,1,4", "7,11,12"]
DEFAULT_SEED = 20260707
GUESSED_MASTER_CELL = 15
MODEL_DATA_EXPONENT = 16
MODEL_ONLINE_EXPONENT = 16


def pack_sequence(sequence: tuple[int, ...]) -> int:
    code = 0
    for index, value in enumerate(sequence):
        code |= value << (4 * index)
    return code


def build_offline_table(
    constraint: tuple[int, ...],
    transform_tables: list[tuple[list[int], list[int], list[int], list[int]]],
) -> set[int]:
    if TABLE_R_DIST != R_DIST:
        raise ValueError(f"offline table rounds mismatch: {TABLE_R_DIST} != {R_DIST}")

    offline_table: set[int] = set()
    stack = [(0, IDENTITY_SEQUENCE_CODE)]
    while stack:
        round_index, sequence_code = stack.pop()
        if round_index == R_DIST:
            if decode_prefix(sequence_code, PREFIX_LEN) == constraint:
                offline_table.add(sequence_code)
            continue
        for tables in transform_tables:
            stack.append((round_index + 1, apply_transform(sequence_code, tables)))
    return offline_table


def round6_cell7_sequence(ciphertexts: list[list[int]], round6_key_7: int, round7_key_4: int) -> tuple[int, ...]:
    rc7_cell4 = (ROUND_CONSTANTS[7] >> 4) & 0x3
    values = []
    for ciphertext in ciphertexts:
        q7_5 = ciphertext[5] ^ ciphertext[9] ^ ciphertext[13]
        q7_10 = ciphertext[6] ^ ciphertext[14]
        q7_15 = ciphertext[3] ^ ciphertext[15]
        x7_4 = INV_SBOX[q7_5 ^ round7_key_4 ^ rc7_cell4]
        x7_8 = INV_SBOX[q7_10 ^ 0x2]
        x7_12 = INV_SBOX[q7_15]
        q6_4 = x7_4 ^ x7_8 ^ x7_12
        values.append(INV_SBOX[q6_4 ^ round6_key_7])
    return difference_sequence([[value] for value in values], 0)


def parse_constraint_list(args: argparse.Namespace) -> list[tuple[int, ...]]:
    if args.constraints is not None:
        constraint_texts = args.constraints
    elif args.constraint is not None:
        constraint_texts = [args.constraint]
    else:
        constraint_texts = DEFAULT_CONSTRAINTS

    constraints = []
    seen = set()
    for text in constraint_texts:
        constraint = parse_prefix(text, PREFIX_LEN)
        if constraint in seen:
            continue
        constraints.append(constraint)
        seen.add(constraint)
    return constraints


def run_attack_for_constraint(
    constraint: tuple[int, ...],
    args: argparse.Namespace,
    transform_tables: list[tuple[list[int], list[int], list[int], list[int]]],
) -> dict[str, object]:
    offline_table = build_offline_table(constraint, transform_tables)
    rng = random.Random(args.seed)
    total_delta_set_attempts = 0
    recovered_trials = 0
    correct_retained = 0
    wrong_submissions = 0
    ambiguous_matches = 0
    candidate_histogram: Counter[int] = Counter()
    attempts_histogram: Counter[int] = Counter()
    examples = []

    for trial in range(args.valid_trials):
        master_key = random_state(rng)
        known_key_template = master_key[:]
        true_guess = master_key[GUESSED_MASTER_CELL]
        guessed_round_keys = []
        for guess in range(16):
            guessed_key = known_key_template[:]
            guessed_key[GUESSED_MASTER_CELL] = guess
            guessed_round_keys.append(round_keys(guessed_key, TOTAL_ROUNDS))
        trial_attempts = 0
        submitted_guess = None
        submitted_sequence = None
        submitted_base_state = None
        submitted_candidates: list[int] = []

        while submitted_guess is None and trial_attempts < args.max_attempts:
            trial_attempts += 1
            total_delta_set_attempts += 1

            base_state = random_state(rng)
            plaintexts = delta_states(base_state, A_CELL)
            ciphertexts = [encrypt_rounds(plaintext, master_key, TOTAL_ROUNDS) for plaintext in plaintexts]

            candidates: list[tuple[int, tuple[int, ...]]] = []
            for guess in range(16):
                guessed_sequence = round6_cell7_sequence(
                    ciphertexts,
                    guessed_round_keys[guess][6][7],
                    guessed_round_keys[guess][7][4],
                )
                if guessed_sequence[:PREFIX_LEN] != constraint:
                    continue
                if pack_sequence(guessed_sequence) not in offline_table:
                    continue
                candidates.append((guess, guessed_sequence))
                if args.stop_on_first_candidate:
                    break

            if candidates:
                submitted_guess = candidates[0][0]
                submitted_sequence = candidates[0][1]
                submitted_base_state = base_state
                submitted_candidates = [guess for guess, _sequence in candidates]

        if submitted_guess is None:
            raise RuntimeError(
                f"No table-matching candidate found in trial {trial}; "
                "increase --max-attempts."
            )

        recovered_trials += 1
        attempts_histogram[trial_attempts] += 1
        candidate_histogram[len(submitted_candidates)] += 1
        if len(submitted_candidates) > 1:
            ambiguous_matches += 1

        if submitted_guess == true_guess:
            correct_retained += 1
        else:
            wrong_submissions += 1

        if len(examples) < 5:
            examples.append(
                {
                    "trial": trial,
                    "attempts": trial_attempts,
                    "true_guess": true_guess,
                    "submitted_guess": submitted_guess,
                    "submitted_candidates": submitted_candidates,
                    "submitted_sequence_prefix": list(submitted_sequence[:PREFIX_LEN]),
                    "submitted_sequence_suffix": list(submitted_sequence[PREFIX_LEN:]),
                    "key": master_key,
                    "base_state": submitted_base_state,
                }
            )

    empirical_constraint_probability = recovered_trials / total_delta_set_attempts
    expected_constraint_probability = 2 ** (-4 * PREFIX_LEN)
    attempts_per_valid = total_delta_set_attempts / recovered_trials
    average_data_search_cost = 16 * attempts_per_valid

    return {
        "experiment": "experiment2_reduced_key_recovery_attack",
        "cipher": "SKINNY-64-64",
        "round_split": "0+6+2",
        "active_input_cell_A": A_CELL,
        "output_cell_B": B_CELL,
        "constraint_prefix": list(constraint),
        "constraint_prefix_text": format_prefix(constraint),
        "offline_table_sequences": len(offline_table),
        "guessed_master_cell": GUESSED_MASTER_CELL,
        "known_key_nibbles_except_guessed_cell": True,
        "requested_attack_trials": args.valid_trials,
        "recovered_trials": recovered_trials,
        "delta_set_attempts": total_delta_set_attempts,
        "oracle_plaintext_queries": 16 * total_delta_set_attempts,
        "key_nibble_guess_tests": 16 * total_delta_set_attempts,
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
        "average_key_nibble_guess_tests": 16 * attempts_per_valid,
        "average_key_nibble_guess_tests_text": format_power_of_two(16 * attempts_per_valid),
        "model_prediction": {
            "data_complexity_text": f"2^{{{MODEL_DATA_EXPONENT}}}",
            "online_time_complexity_text": f"2^{{{MODEL_ONLINE_EXPONENT}}}",
        },
        "correct_retained": correct_retained,
        "correct_retention_rate": correct_retained / recovered_trials,
        "wrong_submissions": wrong_submissions,
        "ambiguous_matches": ambiguous_matches,
        "candidate_histogram": {str(size): count for size, count in sorted(candidate_histogram.items())},
        "attempts_histogram": {
            str(attempts): count for attempts, count in sorted(attempts_histogram.items())
        },
        "examples": examples,
    }


def format_prefix_result(item: dict[str, object]) -> list[str]:
    return [
        f"prefix {item['constraint_prefix_text']}:",
        f"  offline table sequences for prefix: {item['offline_table_sequences']}",
        f"  attack trials: {item['recovered_trials']}/{item['requested_attack_trials']}",
        f"  delta-set attempts: {item['delta_set_attempts']}",
        f"  oracle plaintext queries: {item['oracle_plaintext_queries']}",
        "  empirical table-hit probability per delta-set: "
        f"{item['empirical_constraint_probability_text']} "
        f"({item['empirical_constraint_probability']:.8e})",
        f"  expected CT1 probability: {item['expected_constraint_probability_text']} "
        f"({item['expected_constraint_probability']:.8e})",
        "  average delta-set attempts per recovered nibble: "
        f"{item['attempts_per_valid_delta_set']:.2f} = "
        f"{item['attempts_per_valid_delta_set_text']}",
        f"  average data cost: {item['average_data_search_cost']:.2f} = "
        f"{item['average_data_search_cost_text']}",
        f"  average key-nibble guess tests: {item['average_key_nibble_guess_tests']:.2f} = "
        f"{item['average_key_nibble_guess_tests_text']}",
        f"  correct submitted key nibble: {item['correct_retained']}/{item['recovered_trials']}",
        f"  wrong submitted key nibble: {item['wrong_submissions']}",
        f"  ambiguous matching delta-sets: {item['ambiguous_matches']}",
        f"  candidate histogram: {item['candidate_histogram']}",
    ]


def build_summary(prefix_results: list[dict[str, object]]) -> str:
    lines = [
        "Experiment 2: reduced 0+6+2 key-recovery attack validation",
        f"cipher: SKINNY-64-64, A={{{A_CELL}}}, B={{{B_CELL}}}, s={PREFIX_LEN}",
        f"model prediction: Data=2^{{{MODEL_DATA_EXPONENT}}}, Online=2^{{{MODEL_ONLINE_EXPONENT}}}",
        "",
    ]
    for index, item in enumerate(prefix_results):
        if index:
            lines.append("")
        lines.extend(format_prefix_result(item))

    if len(prefix_results) > 1:
        lines.extend(["", "comparison:"])
        for item in prefix_results:
            lines.append(
                f"  {item['constraint_prefix_text']}: data={item['average_data_search_cost_text']}, "
                f"table-hit={item['empirical_constraint_probability_text']}, "
                f"correct={item['correct_retained']}/{item['recovered_trials']}"
            )
    return "\n".join(lines).rstrip() + "\n"


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Run the reduced 0+6+2 SKINNY-64-64 key-recovery validation."
    )
    parser.add_argument("--valid-trials", type=int, default=64)
    parser.add_argument("--max-attempts", type=int, default=1 << 22)
    parser.add_argument(
        "--constraint",
        default=None,
        help="Run one fixed 3-nibble CT1 prefix, for example 5,1,4.",
    )
    parser.add_argument(
        "--constraints",
        nargs="+",
        default=None,
        help="Run multiple fixed 3-nibble CT1 prefixes. Default: 5,1,4 7,11,12.",
    )
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    parser.add_argument("--self-test-trials", type=int, default=100)
    parser.add_argument(
        "--stop-on-first-candidate",
        action=argparse.BooleanOptionalAction,
        default=False,
        help=(
            "Submit the first table-matching key-nibble candidate to the oracle. "
            "The default is --no-stop-on-first-candidate, which checks all 16 "
            "guesses for each queried delta-set."
        ),
    )
    args = parser.parse_args()

    check_inverse(TOTAL_ROUNDS, trials=args.self_test_trials)

    constraints = parse_constraint_list(args)
    transform_tables = build_transform_tables(build_diff_permutations())
    prefix_results = [
        run_attack_for_constraint(constraint, args, transform_tables)
        for constraint in constraints
    ]

    result = {
        "experiment": "experiment2_reduced_key_recovery_attack_multi_prefix",
        "cipher": "SKINNY-64-64",
        "round_split": "0+6+2",
        "active_input_cell_A": A_CELL,
        "output_cell_B": B_CELL,
        "guessed_master_cell": GUESSED_MASTER_CELL,
        "known_key_nibbles_except_guessed_cell": True,
        "seed": args.seed,
        "requested_attack_trials_per_prefix": args.valid_trials,
        "model_prediction": {
            "data_complexity_text": f"2^{{{MODEL_DATA_EXPONENT}}}",
            "online_time_complexity_text": f"2^{{{MODEL_ONLINE_EXPONENT}}}",
        },
        "prefix_results": prefix_results,
    }

    summary = build_summary(prefix_results)

    out_dir = Path(__file__).resolve().parent
    write_json(out_dir / "results.json", result)
    (out_dir / "summary.txt").write_text(summary, encoding="utf-8")
    print(summary)


if __name__ == "__main__":
    main()
