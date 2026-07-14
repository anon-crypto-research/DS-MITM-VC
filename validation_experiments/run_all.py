#!/usr/bin/env python3
"""Regenerate the reduced SKINNY-64-64 validation results."""

from __future__ import annotations

import shlex
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parent


def run(label: str, command: list[str]) -> None:
    printable = " ".join(shlex.quote(part) for part in command)
    print(f"\n== {label} ==", flush=True)
    print(f"$ {printable}", flush=True)
    subprocess.run(command, cwd=ROOT, check=True)


def main() -> None:
    run(
        "Experiment 1: CT1 value-constraint probability",
        [
            sys.executable,
            "experiment1_value_constraint/run_experiment.py",
            "--constraints",
            "0,0,0",
            "7,11,12",
            "5,1,4",
            "--trials",
            str(1 << 20),
        ],
    )
    run(
        "Experiment 1b: exhaustive ordered sequence distribution",
        [
            sys.executable,
            "experiment1_value_constraint/run_sequence_distribution.py",
            "--quiet",
        ],
    )
    run(
        "Experiment 1c: 10-round random prefix distribution",
        [
            sys.executable,
            "experiment1_value_constraint/run_prefix_distribution.py",
            "--rounds",
            "10",
            "--trials",
            str(1 << 20),
        ],
    )
    run(
        "Experiment 2: reduced 0+6+2 key-recovery attack validation",
        [
            sys.executable,
            "experiment2_reduced_attack/run_experiment.py",
            "--constraints",
            "5,1,4",
            "7,11,12",
            "--valid-trials",
            "64",
            "--max-attempts",
            str(1 << 22),
        ],
    )


if __name__ == "__main__":
    main()
