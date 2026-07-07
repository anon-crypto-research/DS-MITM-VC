# Experiment 2: Reduced 0+6+2 Filtering Validation

This experiment validates the executable filtering chain for the reduced
SKINNY-64-64 case in `../results/output`:

- `r_in = 0`
- `r_dist = 6`
- `r_out = 2`
- `A = {12}`
- `B = {7}`
- `Val_Con = 3`
- `Key = Kg = 1`
- model prediction: `Data = 2^{16}`, `Online = 2^{16}`

The script fixes the reachable CT1 prefix `(7, 11, 12)`, which is checked in
Experiment 1.  It then searches for random delta-sets whose true 6-round
internal sequence satisfies this prefix.  For each constrained delta-set, it
varies the model-indicated master-key nibble 15, decrypts the last two rounds,
and checks whether the table sequence is retained.

Run from this directory:

```bash
python3 run_experiment.py
```

or from the parent `expiriment` directory:

```bash
python3 experiment2_reduced_attack/run_experiment.py
```

Expected output with the default seed:

- 64 valid constrained delta-sets are found in 267773 attempts;
- empirical CT1 probability is `2^{-12.03}`, close to `2^{-12}`;
- the correct key nibble is retained in all 64 cases;
- no wrong key-nibble guess survives;
- average data/online-search cost is `2^{16.03}`, close to the model prediction
  `2^{16}`.

Outputs:

- `results.json`: machine-readable statistics, including log2 costs;
- `summary.txt`: short text summary for the rebuttal and release notes.

Scope: this is a reduced-round filtering validation.  It checks that the
correct key-nibble candidate is retained by the value-constraint/table-sequence
filter; it does not enumerate the full 64-bit key space.
