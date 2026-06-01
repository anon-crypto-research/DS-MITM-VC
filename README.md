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
cipher folder, usually as `results/output`.

## Reproducibility notes

To reproduce the paper results directly, some variables corresponding to 
optimal solutions are fixed in the code. These constraints reduce solving time.
Removing them gives a more open search model, but the solver may take much
longer to finish.

Output files may be overwritten by a new run, so keep a copy of previous results
when comparing different parameter choices.
