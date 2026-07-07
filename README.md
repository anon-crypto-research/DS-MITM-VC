# DS-MITM-VC

This repository contains the MILP models used to search DS-MITM attacks with value constraints.
The Python scripts require Gurobi and its Python package:

```sh
python3 -m pip install gurobipy
```

A valid Gurobi license is also required. See the Gurobi documentation for license
setup if `gurobipy` is installed but the solver cannot start.

## Running the models

Run a model from the repository root, for example:

```sh
python3 SIMON/milp_simon.py
python3 SKINNY/milp_skinny.py
python3 TWINE/milp_twine.py
python3 Rijndael/milp_rijndael.py
python3 Rijndael/milp_rijndael_td.py
python3 Vistrutah/milp_vistrutah.py
```

Each script writes its output to the `results` directory under the corresponding
cipher folder. In this repository the saved files are named by the attack
parameters, for example `Rijndael/results/9R-256-256` and
`SKINNY/results/22R-64-192`; in older scripts the same kind of file may be
called `output`.

## Reading the result files

The files under `*/results/` are direct solver dumps. The first lines summarize
the selected MILP solution, and the later matrices print the binary variables
that define the distinguisher, the outer rounds, the value constraints, and the
key material.

Common fields in the summary lines are:

- `Model Status`: Gurobi status code. Status `2` means that Gurobi proved
  optimality for the configured model.
- `Min_Obj`: value of the primary optimized objective in the current run.
- `Obj` or `Time`: the overall time objective reported by the model. When
  `Offline` and `Online` are printed next to it, they are the precomputation and
  online key-recovery time exponents used to derive the overall time.
- `Data`: data complexity exponent.
- `Memory`: memory complexity exponent.
- `Deg`: number of offline parameters in the DS-MITM distinguisher, i.e. the
  size of the printed `Z` overlap set before reductions such as value
  constraints or key-dependent sieving. In SKINNY output, a line such as
  `Deg = 69 - 17 - 2 = 50` shows the raw distinguisher parameters, the
  non-full-key-addition reduction, the key-sieve reduction, and the resulting
  effective offline parameter count.
- `Plaintext Difference`, `Plain_diff`: dimension of the plaintext structure 
  used for data generation.
- `Key`, `K_on`, `Key_guess`, `Key_space`: independent online key material after
  key-schedule relations are taken into account.
- `Val_Con` or `Val_con`: number of value constraints imposed on the output
  sequence.
- `Kg`, `Kg1`, `Kg2`, `K_g1`, `K_g2`: key material guessed for the table lookup.
  The `1` and `2` suffixes are the two guessing stages in a two-stage lookup.
- `Kbg1`, `Kbg2`: key-bridged sizes of the corresponding guessed key sets.
- `St`, `S_t`, `Sdt`, `Svt`: lookup-related state terms used in the table-size
  and online-time constraints. In SIMON, `Sdt` and `Svt` separate difference and
  value terms. Suffix `1`, as in `St1`, `S_t1`, `Sdt1`, or `Svt1`, denotes the
  part already covered by the first lookup stage.
- `Kt` or `K_t`: key material supplied by the auxiliary lookup table.
- `K_off`, `K_cup`, `Key_sieve`: key-dependent-sieve and weak-key-partition
  counts. `K_off` is the offline key material, `K_cup` is the bridged union of
  offline and online key material, and `Key_sieve` is the resulting reduction
  applied to the offline table.
- `A` and `B`: active input and output cells of the distinguisher. These are the
  main cell lists transferred from the MILP output into the paper tables.

The matrix sections should be read as binary state or key-cell selections. A
`1` means that the corresponding cell is active, known, required, guessed, or
stored, depending on the variable family; a `0` means it is not selected.

- `X`: forward difference propagation through the distinguisher.
- `Y`: backward determination propagation from the output side of the
  distinguisher.
- `Z`: overlap of `X` and `Y`; these cells form the offline parameter set
  `Boff`.
- `M`: plaintext-side outer-round propagation used to generate the required
  input structure.
- `W`: ciphertext-side outer-round propagation used to recover the selected
  output difference sequence.
- `D`, `D1`, `D2`: cells known after online key guesses for the table lookup.
  `D1` and `D2` are the two stages of a two-stage lookup.
- `H`: residual information supplied by, or propagated through, the auxiliary
  lookup table.
- `Key D`, `Key D1`, `Key D2`, `Key H`, `Key V`, `Round Key Guess`: selected
  key cells associated with the corresponding state propagation.

Labels such as `_linear` and `_non_linear` mark the state before and after the
linear and nonlinear parts of a round. Arrows such as `(ARK, SB, SR)`, `(MC)`,
or `(ML)` show which cipher operation connects the printed states.

## Solver output and manual work

The solver output determines the MILP-selected objects: active cell lists `A`
and `B`, the offline parameter set `Boff` through `Z`, the value-constraint
count, lookup-related cells, guessed or table-supplied key cells, and the
reported complexity exponents. The human-written part of the paper translates
these binary selections into attack algorithms, table layouts, equations,
round-by-round descriptions, and final complexity accounting with cipher-specific
normalization and small constants.

## Objective order

The scripts use Gurobi multi-objective minimization. For the paper-level
complexities, the intended optimization order is:

1. Minimize time.
2. Subject to the best time, minimize data.
3. Subject to the best time and data, minimize memory.

## Structure of the MILP scripts

The repository does not have a single file named `milp.py`; instead each cipher
has its own script, such as `Rijndael/milp_rijndael.py`,
`SKINNY/milp_skinny.py`, `SIMON/milp_simon.py`, `TWINE/milp_twine.py`, and
`Vistrutah/milp_vistrutah.py`. They follow the same overall structure:

- `Reset_model(...)`: creates a fresh Gurobi model, initializes global variable
  containers, and creates the objective and complexity expressions.
- `Evaluate_expr(...)`: evaluates a Gurobi linear expression under the current
  solution, mainly for printing derived counts.
- Cipher-operation helpers, such as `Build_shiftrow`, `Build_mixcolumn`,
  `Build_permutation`, `Build_F`, `SR`, `calc`, and key-index helpers: encode
  the propagation rules and key-schedule indexing for the target cipher.
- `Build_distinguisher(...)`: builds the offline DS-MITM distinguisher using the
  `X`, `Y`, and `Z` variables, applies fixed or extra constraints on `A` and
  `B`, and accumulates the offline parameter count `Deg`.
- `Build_key_recovery(...)`: builds the plaintext-side and ciphertext-side
  outer-round constraints, represented by variables such as `M` and `W`, and
  counts the online key material and data-generation requirements.
- Key-counting helpers, such as `Build_key_bridging(...)`,
  `Build_key_counting(...)`, and `Extract_involved_key(...)`: convert selected
  round-key cells into independent key counts after key-schedule relations.
- `Build_value_constraint(...)` or `Build_value_constraints(...)`: adds the
  value-constraint and table-lookup constraints, including variables such as
  `D`, `D1`, `D2`, and `H`, and computes `Val_Con`, `Kg`, `St`, and `Kt`.
- `Build_key_dependent_sieve(...)`: when present, accounts for the
  key-dependent sieve and weak-key partition through `K_off`, `K_cup`, and
  `Key_sieve`.
- `Set_objective(...)`: connects the counted quantities to the time, data, and
  memory complexity expressions and registers the Gurobi objectives.
- Printing helpers, such as `print_block`, `print_flow`, `Print_state`, and
  `Print_result`: format the selected solution into the result files.
- `Start_solver(...)`: configures Gurobi, runs optimization, and writes the
  result dump.
- `Search_attack(...)`, `Search_ds_mitm_attack(...)`, `Parse_args()`, and
  `main()`: provide the top-level search workflow. The usual order is
  model reset, distinguisher construction, key-recovery construction,
  value-constraint construction, optional key-dependent sieve, objective setup,
  optimization, and result printing.

## Reproducibility notes

To reproduce the paper results directly, some variables corresponding to 
optimal solutions are fixed in the code. These constraints reduce solving time.
Removing them gives a more open search model, but the solver may take much
longer to finish.

Output files may be overwritten by a new run, so keep a copy of previous results
when comparing different parameter choices.
