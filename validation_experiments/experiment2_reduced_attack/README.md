# Experiment 2: Reduced key-recovery validation

This directory contains the reduced `0+6+2` SKINNY-64-64 key-recovery
validation for the small instance selected in `../results/output`.

Parameters:

- `r_in = 0`;
- `r_dist = 6`;
- `r_out = 2`;
- active input cell `A = {12}`;
- observed output cell `B = {7}`;
- CT1 value constraints: `Val_Con = 3`;
- online guessed key material: one master-key nibble, `Key = Kg = 1`;
- MILP prediction: `Data = 2^{16}`, `Online = 2^{16}`.

The experiment uses the exhaustive sequence-distribution logic from Experiment
1 to build offline tables for fixed 3-nibble prefixes.  It then simulates an
oracle-style reduced key-recovery loop.

## Attack simulation

For each prefix and each trial:

1. The oracle samples one random master key.
2. The attacker chooses a random `A={12}` delta-set and receives the
   corresponding 16 ciphertexts.
3. The attacker enumerates all 16 guesses for master-key nibble 15.
4. For each guess, the attacker computes the 6-round ordered output sequence at
   `B=7`.
5. If the sequence has the selected prefix and appears in the offline table, the
   attacker submits that key-nibble candidate.
6. The oracle verifies whether the submitted nibble is correct.

All key material outside the model-indicated online nibble is treated as fixed
by the surrounding reduced setting.  This validates recovery of that remaining
nibble, not a full 64-bit key recovery attack.

## Run

Run from this directory:

```bash
python3 run_experiment.py
```

or from the parent `validation_experiments` directory:

```bash
python3 experiment2_reduced_attack/run_experiment.py
```

Default command:

```bash
python3 run_experiment.py --constraints 5,1,4 7,11,12 --valid-trials 64
```

Outputs:

- `results.json`: machine-readable attack statistics for both prefixes;
- `summary.txt`: compact text summary.

## Results

Default-seed results:

| Prefix | Offline table sequences | Table-hit probability | Average data cost | Correct submitted nibble | Wrong submissions |
| --- | ---: | ---: | ---: | ---: | ---: |
| `(5, 1, 4)` | 8048 | `2^{-10.85}` | `2^{14.85}` | 64/64 | 0 |
| `(7, 11, 12)` | 4716 | `2^{-11.66}` | `2^{15.66}` | 64/64 | 0 |

The data cost is:

```text
average data cost = 16 / Pr[one queried delta-set produces a table hit].
```

The `(5, 1, 4)` table is larger than the `(7, 11, 12)` table, so the attacker
finds a table hit faster and needs fewer oracle plaintext queries on average.

## Relation to Experiment 1

The offline table is built by exhaustive enumeration of the six Z nibbles from
Experiment 1.  A table built from random hits is a sampled subset of this
enumerated table: every real random hit induces concrete values for those six Z
nibbles, and those values are included in the `16^6` enumeration.

This is why Experiment 2 uses the enumerated table rather than a random-sample
table.  The random-sampling result estimates how frequently a prefix appears;
the enumerated table provides all full sequences for that prefix in the reduced
model.
