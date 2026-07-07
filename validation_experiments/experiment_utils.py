#!/usr/bin/env python3
"""Shared helpers for the reduced SKINNY validation experiments."""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any


Prefix = tuple[int, ...]


def parse_prefix(text: str, expected_len: int) -> Prefix:
    parts = text.split(",")
    if len(parts) != expected_len:
        raise ValueError(f"prefix must contain {expected_len} comma-separated nibbles")
    values = tuple(int(part, 0) for part in parts)
    if any(value < 0 or value > 15 for value in values):
        raise ValueError("prefix nibbles must be in [0, 15]")
    return values


def parse_prefixes(texts: list[str], expected_len: int) -> list[Prefix]:
    prefixes = []
    seen = set()
    for text in texts:
        prefix = parse_prefix(text, expected_len)
        if prefix not in seen:
            prefixes.append(prefix)
            seen.add(prefix)
    return prefixes


def format_prefix(prefix: Prefix) -> str:
    return "(" + ", ".join(str(value) for value in prefix) + ")"


def log2_or_none(value: float) -> float | None:
    if value <= 0:
        return None
    return math.log2(value)


def format_power_of_two(value: float, digits: int = 2) -> str:
    if value == 0:
        return "0"
    if value < 0:
        raise ValueError("cannot format a negative value as a power of two")
    exponent = math.log2(value)
    nearest = round(exponent)
    if abs(exponent - nearest) < 1e-12:
        return f"2^{{{nearest}}}"
    return f"2^{{{exponent:.{digits}f}}}"


def format_probability(value: float, trials: int | None = None) -> str:
    if value > 0:
        return format_power_of_two(value)
    if trials is None:
        return "0"
    return f"0 (no hit; one-hit resolution {format_power_of_two(1 / trials)})"


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
