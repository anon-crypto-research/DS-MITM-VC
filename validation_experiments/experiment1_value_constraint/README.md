# Experiment 1: Value-constraint and sequence distribution

This directory contains the value-constraint validation for the reduced
SKINNY-64-64 small instance selected in `../results/output`.

Parameters:

- `r_dist = 6`;
- active input cell `A = {12}`;
- observed output cell `B = {7}`;
- prefix length `s = 3` nibbles;
- idealized probability for a fixed 3-nibble prefix: `2^{-12}`.

There are two complementary scripts:

- `run_experiment.py` performs random sampling of full keys and base states;
- `run_sequence_distribution.py` exhaustively enumerates the six
  value-relevant Z nibbles from the MILP output.

## Random prefix sampling

Run from this directory:

```bash
python3 run_experiment.py
```

or from the parent `validation_experiments` directory:

```bash
python3 experiment1_value_constraint/run_experiment.py
```

Default command:

```bash
python3 run_experiment.py --constraints 0,0,0 7,11,12 5,1,4 --trials 1048576
```

Default-seed results:

| Prefix | Hits in `2^{20}` trials | Observed probability | Interpretation |
| --- | ---: | ---: | --- |
| `(0, 0, 0)` | 0 | no hit at `2^{-20}` resolution | Empty/unreachable prefix in the exhaustive distribution. |
| `(7, 11, 12)` | 277 | `2^{-11.89}` | Original reachable reference prefix. |
| `(5, 1, 4)` | 523 | `2^{-10.97}` | Largest prefix bucket in the exhaustive distribution. |

Outputs:

- `results.json`: machine-readable sampling statistics and examples;
- `summary.txt`: compact text summary.

## Exhaustive sequence distribution

Run:

```bash
python3 run_sequence_distribution.py
```

This script parses the Z path from `../results/output`:

```text
Z[0][12] -> Z[1][3] -> Z[2][3] -> Z[3][3] -> Z[4][3] -> Z[5][3] -> Z[6][7]
```

It enumerates the six pre-output Z nibbles:

```text
Z[0][12], Z[1][3], Z[2][3], Z[3][3], Z[4][3], Z[5][3]
```

For each of the `16^6 = 2^{24}` assignments, it computes the full ordered
sequence of 15 output differences at `B=7`, then buckets the sequence by its
first three differences.

Default-seed/exhaustive results:

| Quantity | Value |
| --- | ---: |
| Total assignments | `2^{24}` |
| Distinct full sequences | 16691380 |
| Duplicate assignments | 85836 |
| Non-empty prefix buckets | 2730/4096 |
| Empty prefix buckets | 1366/4096 |
| Ideal count per 12-bit prefix bucket | `2^{12}` |
| Ideal count per reachable nonzero-distinct prefix bucket | `2^{12.59}` |

Selected buckets:

| Prefix | Assignments | Distinct sequences | Probability |
| --- | ---: | ---: | ---: |
| `(0, 0, 0)` | 0 | 0 | 0 |
| `(7, 11, 12)` | 4728 | 4716 | about `2^{-11.79}` |
| `(5, 1, 4)` | 8060 | 8048 | `2^{-11.02}` |

Outputs:

- `sequence_distribution.json`: all 4096 prefix buckets, including empty
  buckets, multiplicity histograms, and representative top sequences;
- `sequence_distribution_summary.txt`: compact text summary.

## Interpretation

The idealized `2^{-12}` estimate assumes that the first three output differences
behave like three independent uniform 4-bit values.  This small 6-round core
does not satisfy that assumption.  The 15 nonzero input differences pass through
S-box differential maps as nonzero-difference permutations, so prefixes
containing zero or repeated values are structurally excluded.  This explains the
1366 empty buckets.

The non-empty buckets are also uneven.  The larger `(5, 1, 4)` bucket produces a
higher random hit rate and lower Experiment 2 data cost than `(7, 11, 12)`.
