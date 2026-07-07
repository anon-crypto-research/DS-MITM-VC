#!/usr/bin/env python3
"""Minimal SKINNY-64-64 implementation used by the validation experiments.

The state and master key are represented as lists of 16 nibbles.  This module
is intentionally dependency-free so that the reduced-round experiments can be
run without Gurobi or any third-party Python package.
"""

from __future__ import annotations

import random


State = list[int]

SBOX = [0xC, 0x6, 0x9, 0x0, 0x1, 0xA, 0x2, 0xB, 0x3, 0x8, 0x5, 0xD, 0x4, 0xE, 0x7, 0xF]
INV_SBOX = [0] * 16
for index, value in enumerate(SBOX):
    INV_SBOX[value] = index

SR = [0, 1, 2, 3, 7, 4, 5, 6, 10, 11, 8, 9, 13, 14, 15, 12]
PT = [9, 15, 8, 13, 10, 14, 12, 11, 0, 1, 2, 3, 4, 5, 6, 7]

ROUND_CONSTANTS = [
    0x01, 0x03, 0x07, 0x0F, 0x1F, 0x3E, 0x3D, 0x3B,
    0x37, 0x2F, 0x1E, 0x3C, 0x39, 0x33, 0x27, 0x0E,
    0x1D, 0x3A, 0x35, 0x2B, 0x16, 0x2C, 0x18, 0x30,
    0x21, 0x02, 0x05, 0x0B, 0x17, 0x2E, 0x1C, 0x38,
    0x31, 0x23, 0x06, 0x0D, 0x1B, 0x36, 0x2D, 0x1A,
    0x34, 0x29, 0x12, 0x24, 0x08, 0x11, 0x22, 0x04,
    0x09, 0x13, 0x26, 0x0C, 0x19, 0x32, 0x25, 0x0A,
    0x15, 0x2A, 0x14, 0x28, 0x10, 0x20,
]


def mix_columns(state: State) -> State:
    out = [0] * 16
    for column in range(4):
        a0 = state[column]
        a1 = state[column + 4]
        a2 = state[column + 8]
        a3 = state[column + 12]
        out[column] = a0 ^ a2 ^ a3
        out[column + 4] = a0
        out[column + 8] = a1 ^ a2
        out[column + 12] = a0 ^ a2
    return out


def inv_mix_columns(state: State) -> State:
    out = [0] * 16
    for column in range(4):
        b0 = state[column]
        b1 = state[column + 4]
        b2 = state[column + 8]
        b3 = state[column + 12]
        out[column] = b1
        out[column + 4] = b1 ^ b2 ^ b3
        out[column + 8] = b1 ^ b3
        out[column + 12] = b0 ^ b3
    return out


def round_keys(master_key: State, rounds: int) -> list[State]:
    tk = master_key[:]
    keys = []
    for _ in range(rounds):
        keys.append(tk[:])
        tk = [tk[PT[index]] for index in range(16)]
    return keys


def skinny_round(state: State, round_key: State, round_number: int) -> State:
    state = [SBOX[value] for value in state]
    rc = ROUND_CONSTANTS[round_number]
    state[0] ^= rc & 0xF
    state[4] ^= (rc >> 4) & 0x3
    state[8] ^= 0x2
    for index in range(8):
        state[index] ^= round_key[index]
    state = [state[SR[index]] for index in range(16)]
    return mix_columns(state)


def inv_skinny_round(state: State, round_key: State, round_number: int) -> State:
    state = inv_mix_columns(state)
    pre_shift_rows = [0] * 16
    for index in range(16):
        pre_shift_rows[SR[index]] = state[index]
    state = pre_shift_rows
    for index in range(8):
        state[index] ^= round_key[index]
    rc = ROUND_CONSTANTS[round_number]
    state[0] ^= rc & 0xF
    state[4] ^= (rc >> 4) & 0x3
    state[8] ^= 0x2
    return [INV_SBOX[value] for value in state]


def encrypt_rounds(state: State, master_key: State, rounds: int) -> State:
    keys = round_keys(master_key, rounds)
    state = state[:]
    for round_number in range(rounds):
        state = skinny_round(state, keys[round_number], round_number)
    return state


def encrypt_with_trace(state: State, master_key: State, rounds: int) -> list[State]:
    keys = round_keys(master_key, rounds)
    state = state[:]
    trace = [state[:]]
    for round_number in range(rounds):
        state = skinny_round(state, keys[round_number], round_number)
        trace.append(state[:])
    return trace


def decrypt_to_round(ciphertext: State, master_key: State, start_round: int, total_rounds: int) -> State:
    keys = round_keys(master_key, total_rounds)
    state = ciphertext[:]
    for round_number in range(total_rounds - 1, start_round - 1, -1):
        state = inv_skinny_round(state, keys[round_number], round_number)
    return state


def random_state(rng: random.Random) -> State:
    return [rng.randrange(16) for _ in range(16)]


def delta_states(base_state: State, active_cell: int, count: int = 16) -> list[State]:
    states = []
    for delta in range(count):
        state = base_state[:]
        state[active_cell] ^= delta
        states.append(state)
    return states


def difference_sequence(states: list[State], cell: int) -> tuple[int, ...]:
    return tuple(states[0][cell] ^ states[index][cell] for index in range(1, len(states)))


def difference_prefix(
    master_key: State,
    base_state: State,
    active_cell: int,
    output_cell: int,
    rounds: int,
    prefix_len: int,
) -> tuple[int, ...]:
    outputs = []
    for state in delta_states(base_state, active_cell, prefix_len + 1):
        outputs.append(encrypt_rounds(state, master_key, rounds))
    return difference_sequence(outputs, output_cell)


def check_inverse(rounds: int, trials: int = 100, seed: int = 1) -> None:
    rng = random.Random(seed)
    for _ in range(trials):
        key = random_state(rng)
        state = random_state(rng)
        trace = encrypt_with_trace(state, key, rounds)
        back = trace[rounds]
        keys = round_keys(key, rounds)
        for round_number in range(rounds - 1, -1, -1):
            back = inv_skinny_round(back, keys[round_number], round_number)
        assert back == state
