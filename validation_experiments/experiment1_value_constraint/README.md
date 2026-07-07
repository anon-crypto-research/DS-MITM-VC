# Experiment 1: Value-Constraint Probability

This experiment checks the probability of the three CT1 value constraints in
the reduced SKINNY-64-64 setting used for the executable validation.

Parameters:

- `r_dist = 6`
- active input cell `A = {12}`
- output cell `B = {7}`
- CT1 prefix length `s = 3`
- idealized probability for a fixed prefix: `2^{-12}`

Run from this directory:

```bash
python3 run_experiment.py
```

or from the parent `expiriment` directory:

```bash
python3 experiment1_value_constraint/run_experiment.py
```

The default run tests two prefixes:

```bash
python3 run_experiment.py --constraints 0,0,0 7,11,12 --trials 1048576
```

Expected output with the default seed:

- prefix `(0, 0, 0)`: `0` hits in `2^{20}` trials;
- prefix `(7, 11, 12)`: `277` hits in `2^{20}` trials, i.e. `2^{-11.89}`.

The all-zero prefix is included on purpose.  It shows that the 6-round core is
not uniformly distributed in this very short setting.  The reachable prefix
`(7, 11, 12)` is close to the expected `2^{-12}` scale and is used by
Experiment 2 for the key-recovery validation.

Outputs:

- `results.json`: machine-readable statistics, including log2 probabilities;
- `summary.txt`: short text summary for the rebuttal and release notes.
