# Reduced SKINNY-64-64 validation experiments

This directory contains two small executable experiments added for the ToSC
review response.  They validate, on a practical reduced-round instance, that the
MILP output can be translated into an executable filtering procedure and that
the observed cost is consistent with the predicted exponents.

The reduced case is the 0+6+2 SKINNY-64-64 setting recorded in
`results/output`:

- `r_in = 0`, `r_dist = 6`, `r_out = 2`
- active input cell `A = {12}`
- output cell `B = {7}`
- three CT1 value constraints, with expected probability `2^{-12}`
- one online-guessed master-key nibble
- model prediction: `Data = 2^{16}`, `Online = 2^{16}`

The validation scripts use only the Python standard library.  Re-running
`milp_skinny.py` requires the original MILP/Gurobi environment, but these
experiments do not.

## Run

From this directory:

```bash
python3 run_all.py
```

or run each experiment separately:

```bash
python3 experiment1_value_constraint/run_experiment.py
python3 experiment2_reduced_attack/run_experiment.py
```

The scripts overwrite their local `results.json` and `summary.txt` files.

## Experiment 1: value-constraint probability

`experiment1_value_constraint` checks how often the first three differences at
cell `B=7` after the 6-round core equal a fixed CT1 prefix.

By default it tests two prefixes:

| Prefix | Role | Observed result |
| --- | --- | --- |
| `(0, 0, 0)` | sanity check for non-uniformity | no hit in `2^{20}` trials |
| `(7, 11, 12)` | reachable prefix used by Experiment 2 | `277/2^{20} = 2^{-11.89}` |

The idealized estimate for a fixed 3-nibble prefix is `2^{-12}`.  The all-zero
prefix is intentionally included to show that 6-round SKINNY-64 is not behaving
as an ideal permutation in this short setting.  The reachable prefix
`(7, 11, 12)` has probability close to the expected scale and is therefore used
for the executable key-recovery validation.

## Experiment 2: reduced key-recovery validation

`experiment2_reduced_attack` fixes the reachable prefix `(7, 11, 12)` and runs
the induced 8-round filtering procedure.  For each valid delta-set, it enumerates
the model-indicated online master-key nibble and checks which guesses reproduce
the table sequence after decrypting the last two rounds.

With the default seed, 64 valid constrained delta-sets are found in 267773
attempts.  In all 64 cases the correct key nibble is retained and no wrong
key-nibble guess survives.  The average data/online-search cost is
`2^{16.03}`, close to the model prediction `2^{16}`.

## Interpretation

These experiments are a sanity check for the implementation and the complexity
estimate on a fully executable reduced instance.  They are not intended to prove
that every 6-round SKINNY-64 prefix is uniformly distributed.  On the contrary,
Experiment 1 explicitly shows visible non-uniformity for the short 6-round core.
For the longer distinguishers used in the paper, the output behavior is expected
to be closer to random, making the independence heuristic more reliable.
