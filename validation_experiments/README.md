# Small-instance SKINNY-64-64 validation experiments

This directory is a self-contained validation package for a reduced
SKINNY-64-64 instance.  It is intended for code release and reviewer inspection:
the scripts translate one MILP-selected small instance into executable
experiments, then record the observed value-constraint and key-recovery
behavior.

The reduced instance comes from `results/output`:

- round split: `r_in = 0`, `r_dist = 6`, `r_out = 2`;
- active input cell: `A = {12}`;
- observed output cell: `B = {7}`;
- CT1 value constraints: `Val_Con = 3`;
- online key recovery target: one master-key nibble, `Key = Kg = 1`;
- MILP prediction for the reduced instance: `Data = 2^{16}`,
  `Online = 2^{16}`.

The experiments use only the Python standard library.  Re-running
`milp_skinny.py` requires the original MILP/Gurobi environment, but the
validation scripts do not.

## Directory layout

| Path | Purpose |
| --- | --- |
| `skinny64.py` | Minimal SKINNY-64-64 implementation used by the experiments. |
| `experiment_utils.py` | Shared formatting, parsing, and JSON helpers. |
| `milp_skinny.py` | Reduced copy of the MILP model that produced `results/output`; not needed for re-running the validation. |
| `results/output` | MILP-selected small instance and variable dump. |
| `experiment1_value_constraint/` | Prefix-probability, all-prefix random distribution, and exhaustive sequence-distribution experiments. |
| `experiment2_reduced_attack/` | Oracle-style reduced key-recovery validation. |
| `run_all.py` | Regenerates all experiment summaries and JSON result files. |

## Reproduce the results

From this directory:

```bash
python3 run_all.py
```

This regenerates:

- `experiment1_value_constraint/results.json`;
- `experiment1_value_constraint/summary.txt`;
- `experiment1_value_constraint/sequence_distribution.json`;
- `experiment1_value_constraint/sequence_distribution_summary.txt`;
- `experiment1_value_constraint/prefix_distribution_r10.json`;
- `experiment1_value_constraint/prefix_distribution_r10_summary.txt`;
- `experiment2_reduced_attack/results.json`;
- `experiment2_reduced_attack/summary.txt`.

Individual experiments can also be run directly:

```bash
python3 experiment1_value_constraint/run_experiment.py
python3 experiment1_value_constraint/run_sequence_distribution.py
python3 experiment1_value_constraint/run_prefix_distribution.py
python3 experiment2_reduced_attack/run_experiment.py
```

The scripts overwrite their local result files.

## Experiment 1: value-constraint behavior

`experiment1_value_constraint/run_experiment.py` samples random master keys and
base states, then measures how often the 6-round output-difference prefix at
cell `B=7` equals a fixed 3-nibble value.

With the default seed and `2^{20}` trials:

| Prefix | Hits | Observed probability | Comment |
| --- | ---: | ---: | --- |
| `(0, 0, 0)` | 0 | no hit at `2^{-20}` resolution | Demonstrates short-round non-uniformity. |
| `(7, 11, 12)` | 277 | `2^{-11.89}` | Original reachable reference prefix. |
| `(5, 1, 4)` | 523 | `2^{-10.97}` | Largest prefix bucket in the exhaustive distribution. |

`experiment1_value_constraint/run_sequence_distribution.py` exhaustively
enumerates the six value-relevant Z nibbles selected by the MILP:

```text
Z[0][12], Z[1][3], Z[2][3], Z[3][3], Z[4][3], Z[5][3]
```

For each of the `16^6 = 2^{24}` assignments, it computes the full 15-nibble
ordered output-difference sequence at `B=7` and buckets it by the first three
differences.

Key exhaustive-distribution results:

| Quantity | Value |
| --- | ---: |
| Total assignments | `2^{24}` |
| Distinct full sequences | 16691380 |
| Non-empty 3-nibble prefix buckets | 2730/4096 |
| Empty 3-nibble prefix buckets | 1366/4096 |
| `(5, 1, 4)` bucket size | 8060 assignments, 8048 distinct sequences |
| `(7, 11, 12)` bucket size | 4728 assignments, 4716 distinct sequences |
| `(0, 0, 0)` bucket size | 0 |

The idealized value-constraint estimate is `2^{-12}` for a fixed 3-nibble
prefix.  The exhaustive result shows why this small instance is not exactly
ideal: only `15 * 14 * 13 = 2730` of the 4096 possible prefixes are reachable,
and the reachable buckets are not uniform.

As an additional diagnostic, `experiment1_value_constraint/run_prefix_distribution.py`
keeps `A={12}`, `B={7}`, and prefix length 3, but increases the middle core to
`r_dist=10` and samples all 4096 prefix buckets.  With `2^{20}` random trials,
all buckets are non-empty and the hit-count distribution is close to the
uniform expectation of 256 hits per prefix:

| Quantity | Value |
| --- | ---: |
| Non-empty prefix buckets | 4096/4096 |
| Min / max hits | 195 / 317 |
| Mean / stddev hits | 256.00 / 16.09 |
| Reduced chi-square vs uniform | 1.0121 |
| `(0, 0, 0)` hits | 254 = about `2^{-12.01}` |
| `(7, 11, 12)` hits | 268 = about `2^{-11.93}` |
| `(5, 1, 4)` hits | 246 = about `2^{-12.06}` |

## Experiment 2: reduced key recovery

`experiment2_reduced_attack/run_experiment.py` simulates the reduced
key-recovery loop for two prefixes, `(5, 1, 4)` and `(7, 11, 12)`.

For each prefix, the script first builds the offline table by enumerating the
six Z nibbles from Experiment 1.  It then repeats 64 oracle trials:

1. The oracle samples a random master key.
2. The attacker queries a random `A={12}` delta-set and receives 16 ciphertexts.
3. The attacker enumerates all 16 guesses for the model-indicated master-key
   nibble 15.
4. For each guess, the attacker computes the 6-round sequence at `B=7`.
5. If the sequence has the selected prefix and is present in the offline table,
   the attacker submits that nibble to the oracle.

Default-seed results:

| Prefix | Offline table sequences | Correct submitted nibble | Wrong submissions | Average data cost |
| --- | ---: | ---: | ---: | ---: |
| `(5, 1, 4)` | 8048 | 64/64 | 0 | `2^{14.85}` |
| `(7, 11, 12)` | 4716 | 64/64 | 0 | `2^{15.66}` |

The data cost is controlled by the table-hit probability:

```text
average data cost = 16 / Pr[one queried delta-set produces a table hit].
```

The larger `(5, 1, 4)` bucket therefore gives lower observed data cost than
`(7, 11, 12)`.

## Interpretation and scope

These experiments validate the executable pipeline on a deliberately small
instance: MILP output, value-constraint table construction, and reduced
key-nibble recovery are consistent with each other.  They should not be read as
evidence that 6-round SKINNY-64-64 behaves like an ideal permutation.  On the
contrary, the all-zero prefix and the unequal bucket sizes explicitly show
short-round non-uniformity.

The key-recovery experiment is also reduced in scope: all key material outside
the model-indicated online nibble is treated as fixed by the surrounding small
instance.  The experiment validates recovery of that remaining nibble, not a
full 64-bit key recovery attack.
