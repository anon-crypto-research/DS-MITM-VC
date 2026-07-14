#!/usr/bin/env python3
"""Exhaustive ordered-output-sequence distribution for Experiment 1.

The MILP output fixes the value-relevant Z path

    Z[0][12] -> Z[1][3] -> ... -> Z[5][3] -> Z[6][7].

For each assignment of the six pre-output Z nibbles, this script computes the
ordered sequence of the 15 nonzero output differences at B=7 and buckets the
full sequence by its first three differences.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from experiment_utils import format_power_of_two, format_prefix  # noqa: E402
from skinny64 import SBOX, SR, mix_columns  # noqa: E402


A_CELL = 12
B_CELL = 7
R_DIST = 6
PREFIX_LEN = 3
DEFAULT_OUTPUT_PATH = ROOT / "results" / "output"

NIBBLES = 16
NONZERO_NIBBLES = tuple(range(1, NIBBLES))
FULL_SEQUENCE_LEN = len(NONZERO_NIBBLES)
IDENTITY_SEQUENCE_CODE = sum(value << (4 * index) for index, value in enumerate(NONZERO_NIBBLES))


@dataclass
class BucketStats:
    total_assignments: int = 0
    distinct_sequences: int = 0
    multiplicity_histogram: Counter[int] = field(default_factory=Counter)
    top_sequences: list[tuple[int, int]] = field(default_factory=list)

    def add(self, sequence_code: int, multiplicity: int, top_limit: int) -> None:
        self.total_assignments += multiplicity
        self.distinct_sequences += 1
        self.multiplicity_histogram[multiplicity] += 1
        if top_limit <= 0:
            return
        self.top_sequences.append((multiplicity_sort_key(multiplicity, sequence_code), sequence_code))
        self.top_sequences.sort(reverse=True)
        if len(self.top_sequences) > top_limit:
            self.top_sequences.pop()


def multiplicity_sort_key(multiplicity: int, sequence_code: int) -> int:
    return (multiplicity << 64) - sequence_code


def decode_sequence(code: int, length: int = FULL_SEQUENCE_LEN) -> list[int]:
    return [(code >> (4 * index)) & 0xF for index in range(length)]


def decode_prefix(prefix_code: int, length: int = PREFIX_LEN) -> tuple[int, ...]:
    return tuple((prefix_code >> (4 * index)) & 0xF for index in range(length))


def parse_z_positions(output_path: Path) -> list[int]:
    lines = output_path.read_text(encoding="utf-8").splitlines()
    in_z_section = False
    positions: list[int] = []
    index = 0

    while index < len(lines):
        line = lines[index].strip()
        if line == "---------- Var Z ----------":
            in_z_section = True
            index += 1
            continue
        if in_z_section and line.startswith("----------"):
            break
        if in_z_section and line.startswith("Z["):
            round_index = int(line.split("[", 1)[1].split("]", 1)[0].strip())
            matrix = []
            for row_offset in range(1, 5):
                matrix.extend(int(value) for value in lines[index + row_offset].split())
            active_positions = [cell for cell, value in enumerate(matrix) if value == 1]
            if len(active_positions) != 1:
                raise ValueError(
                    f"Z[{round_index}] must contain exactly one active cell; "
                    f"got {active_positions}"
                )
            positions.append(active_positions[0])
            index += 5
            continue
        index += 1

    if len(positions) != R_DIST + 1:
        raise ValueError(f"expected {R_DIST + 1} Z rounds, got {len(positions)}")
    if positions[0] != A_CELL or positions[-1] != B_CELL:
        raise ValueError(f"unexpected Z path endpoints: {positions[0]} -> {positions[-1]}")
    return positions


def linear_image_positions(source_position: int) -> set[int]:
    state = [0] * NIBBLES
    state[source_position] = 1
    shifted = [state[SR[index]] for index in range(NIBBLES)]
    mixed = mix_columns(shifted)
    return {index for index, value in enumerate(mixed) if value != 0}


def validate_z_path(z_positions: list[int]) -> list[dict[str, object]]:
    checks = []
    for round_index, (source, target) in enumerate(zip(z_positions, z_positions[1:])):
        image = sorted(linear_image_positions(source))
        if target not in image:
            raise ValueError(
                f"Z path is inconsistent at round {round_index}: "
                f"{source} maps to {image}, not {target}"
            )
        checks.append(
            {
                "round": round_index,
                "source": source,
                "target": target,
                "linear_image": image,
            }
        )
    return checks


def build_diff_permutations() -> list[list[int]]:
    permutations = []
    for base_value in range(NIBBLES):
        permutation = [0]
        permutation.extend(SBOX[base_value] ^ SBOX[base_value ^ diff] for diff in NONZERO_NIBBLES)
        if sorted(permutation[1:]) != list(NONZERO_NIBBLES):
            raise ValueError(f"S-box difference map for base value {base_value} is not a permutation")
        permutations.append(permutation)
    return permutations


def build_transform_tables(permutations: list[list[int]]) -> list[tuple[list[int], list[int], list[int], list[int]]]:
    tables = []
    for permutation in permutations:
        chunks = []
        for chunk_index, chunk_nibbles in enumerate((4, 4, 4, 3)):
            table_size = 1 << (4 * chunk_nibbles)
            base_shift = 16 * chunk_index
            table = [0] * table_size
            for chunk in range(table_size):
                transformed = 0
                for local_index in range(chunk_nibbles):
                    diff = (chunk >> (4 * local_index)) & 0xF
                    if diff:
                        transformed |= permutation[diff] << (base_shift + 4 * local_index)
                table[chunk] = transformed
            chunks.append(table)
        tables.append((chunks[0], chunks[1], chunks[2], chunks[3]))
    return tables


def apply_transform(code: int, tables: tuple[list[int], list[int], list[int], list[int]]) -> int:
    return (
        tables[0][code & 0xFFFF]
        | tables[1][(code >> 16) & 0xFFFF]
        | tables[2][(code >> 32) & 0xFFFF]
        | tables[3][(code >> 48) & 0xFFF]
    )


def enumerate_sequence_codes(rounds: int, verbose: bool) -> list[int]:
    transform_tables = build_transform_tables(build_diff_permutations())
    codes = [IDENTITY_SEQUENCE_CODE]
    for round_index in range(rounds):
        next_codes: list[int] = []
        next_codes_append = next_codes.append
        for code in codes:
            for tables in transform_tables:
                next_codes_append(apply_transform(code, tables))
        codes = next_codes
        if verbose:
            message = f"round {round_index + 1}: assignments={len(codes)}"
            if len(codes) <= (1 << 20):
                message += f" distinct_sequences={len(set(codes))}"
            print(message, flush=True)
    return codes


def iter_runs(sorted_codes: Iterable[int]) -> Iterable[tuple[int, int]]:
    iterator = iter(sorted_codes)
    try:
        current = next(iterator)
    except StopIteration:
        return
    count = 1
    for code in iterator:
        if code == current:
            count += 1
        else:
            yield current, count
            current = code
            count = 1
    yield current, count


def sequence_to_hex(sequence: list[int]) -> str:
    return "".join(format(value, "x") for value in sequence)


def make_jsonable_bucket(prefix_code: int, stats: BucketStats) -> dict[str, object]:
    top_sequences = []
    for packed_key, sequence_code in sorted(stats.top_sequences, reverse=True):
        multiplicity = packed_key >> 64
        sequence = decode_sequence(sequence_code)
        top_sequences.append(
            {
                "sequence": sequence,
                "sequence_hex": sequence_to_hex(sequence),
                "multiplicity": multiplicity,
            }
        )
    return {
        "prefix": list(decode_prefix(prefix_code)),
        "prefix_text": format_prefix(decode_prefix(prefix_code)),
        "total_assignments": stats.total_assignments,
        "distinct_sequences": stats.distinct_sequences,
        "duplicate_assignments": stats.total_assignments - stats.distinct_sequences,
        "multiplicity_histogram": {
            str(multiplicity): count
            for multiplicity, count in sorted(stats.multiplicity_histogram.items())
        },
        "top_sequences": top_sequences,
    }


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Exhaustively bucket 6-round SKINNY ordered output sequences."
    )
    parser.add_argument("--milp-output", type=Path, default=DEFAULT_OUTPUT_PATH)
    parser.add_argument("--top-buckets", type=int, default=12)
    parser.add_argument("--top-sequences-per-bucket", type=int, default=5)
    parser.add_argument("--quiet", action="store_true")
    args = parser.parse_args()

    z_positions = parse_z_positions(args.milp_output)
    z_path_checks = validate_z_path(z_positions)

    codes = enumerate_sequence_codes(R_DIST, verbose=not args.quiet)
    total_assignments = len(codes)
    codes.sort()

    buckets: dict[int, BucketStats] = {}
    multiplicity_histogram: Counter[int] = Counter()
    for sequence_code, multiplicity in iter_runs(codes):
        prefix_code = sequence_code & ((1 << (4 * PREFIX_LEN)) - 1)
        bucket = buckets.setdefault(prefix_code, BucketStats())
        bucket.add(sequence_code, multiplicity, args.top_sequences_per_bucket)
        multiplicity_histogram[multiplicity] += 1

    distinct_sequences = sum(bucket.distinct_sequences for bucket in buckets.values())
    reachable_bucket_count = len(buckets)
    all_bucket_count = 1 << (4 * PREFIX_LEN)
    empty_bucket_count = all_bucket_count - reachable_bucket_count
    ideal_per_bucket = total_assignments / all_bucket_count
    permutation_prefix_bucket_count = 15 * 14 * 13
    ideal_per_permutation_prefix_bucket = total_assignments / permutation_prefix_bucket_count
    bucket_total_histogram = Counter(
        buckets.get(prefix_code, BucketStats()).total_assignments
        for prefix_code in range(all_bucket_count)
    )

    bucket_items = sorted(
        buckets.items(),
        key=lambda item: (
            item[1].total_assignments,
            item[1].distinct_sequences,
            -item[0],
        ),
        reverse=True,
    )
    top_buckets = [
        make_jsonable_bucket(prefix_code, stats)
        for prefix_code, stats in bucket_items[: args.top_buckets]
    ]

    selected_prefixes = {}
    for prefix in ((0, 0, 0), (7, 11, 12)):
        prefix_code = sum(value << (4 * index) for index, value in enumerate(prefix))
        stats = buckets.get(prefix_code, BucketStats())
        selected_prefixes[format_prefix(prefix)] = make_jsonable_bucket(prefix_code, stats)
    all_buckets = [
        make_jsonable_bucket(prefix_code, buckets.get(prefix_code, BucketStats()))
        for prefix_code in range(all_bucket_count)
    ]

    result = {
        "experiment": "experiment1_ordered_output_sequence_distribution",
        "cipher": "SKINNY-64-64",
        "r_dist": R_DIST,
        "active_input_cell_A": A_CELL,
        "output_cell_B": B_CELL,
        "z_positions": z_positions,
        "enumerated_z_parameters": [
            {"round": round_index, "cell": cell}
            for round_index, cell in enumerate(z_positions[:-1])
        ],
        "z_path_checks": z_path_checks,
        "ordered_input_differences": list(NONZERO_NIBBLES),
        "sequence_length": FULL_SEQUENCE_LEN,
        "bucket_prefix_length": PREFIX_LEN,
        "total_assignments": total_assignments,
        "total_assignments_text": format_power_of_two(total_assignments),
        "distinct_sequences": distinct_sequences,
        "duplicate_assignments": total_assignments - distinct_sequences,
        "reachable_bucket_count": reachable_bucket_count,
        "empty_bucket_count": empty_bucket_count,
        "all_bucket_count": all_bucket_count,
        "ideal_assignments_per_12bit_prefix_bucket": ideal_per_bucket,
        "ideal_assignments_per_12bit_prefix_bucket_text": format_power_of_two(ideal_per_bucket),
        "permutation_prefix_bucket_count": permutation_prefix_bucket_count,
        "ideal_assignments_per_nonzero_distinct_prefix_bucket": ideal_per_permutation_prefix_bucket,
        "ideal_assignments_per_nonzero_distinct_prefix_bucket_text": format_power_of_two(
            ideal_per_permutation_prefix_bucket
        ),
        "bucket_total_assignment_histogram": {
            str(total): count for total, count in sorted(bucket_total_histogram.items())
        },
        "sequence_multiplicity_histogram": {
            str(multiplicity): count
            for multiplicity, count in sorted(multiplicity_histogram.items())
        },
        "selected_prefixes": selected_prefixes,
        "top_buckets_by_total_assignments": top_buckets,
        "buckets": all_buckets,
    }

    lines = [
        "Experiment 1b: ordered output sequence distribution",
        f"cipher: SKINNY-64-64, r_dist={R_DIST}, A={{{A_CELL}}}, B={{{B_CELL}}}",
        "Z path: " + " -> ".join(f"Z[{round_index}][{cell}]" for round_index, cell in enumerate(z_positions)),
        "enumerated Z parameters: "
        + ", ".join(f"Z[{round_index}][{cell}]" for round_index, cell in enumerate(z_positions[:-1])),
        f"assignments: {total_assignments} = {format_power_of_two(total_assignments)}",
        f"distinct full sequences: {distinct_sequences}",
        f"duplicate assignments: {total_assignments - distinct_sequences}",
        f"non-empty prefix buckets: {reachable_bucket_count}/{all_bucket_count}",
        f"empty prefix buckets: {empty_bucket_count}/{all_bucket_count}",
        f"ideal count per 12-bit 3-nibble bucket: {ideal_per_bucket:.2f} = "
        f"{format_power_of_two(ideal_per_bucket)}",
        f"ideal count per nonzero-distinct 3-nibble bucket: "
        f"{ideal_per_permutation_prefix_bucket:.2f} = "
        f"{format_power_of_two(ideal_per_permutation_prefix_bucket)}",
        "",
        "selected buckets:",
    ]
    for prefix_text, item in selected_prefixes.items():
        lines.extend(
            [
                f"  prefix {prefix_text}:",
                f"    total assignments: {item['total_assignments']}",
                f"    distinct sequences: {item['distinct_sequences']}",
                f"    duplicate assignments: {item['duplicate_assignments']}",
                f"    multiplicity histogram: {item['multiplicity_histogram']}",
            ]
        )
    lines.extend(["", f"top {args.top_buckets} buckets by total assignments:"])
    for item in top_buckets:
        probability = item["total_assignments"] / total_assignments
        lines.append(
            f"  {item['prefix_text']}: total={item['total_assignments']}, "
            f"distinct={item['distinct_sequences']}, probability={format_power_of_two(probability)}"
        )
    summary = "\n".join(lines).rstrip() + "\n"

    out_dir = Path(__file__).resolve().parent
    (out_dir / "sequence_distribution.json").write_text(
        json.dumps(result, indent=2) + "\n", encoding="utf-8"
    )
    (out_dir / "sequence_distribution_summary.txt").write_text(summary, encoding="utf-8")
    print(summary)


if __name__ == "__main__":
    main()
