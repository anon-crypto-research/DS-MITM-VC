import argparse
from pathlib import Path

from gurobipy import Model, GRB, quicksum, LinExpr, QuadExpr, Var

DEFAULT_RESULTS_DIR = Path(__file__).resolve().parent / "results"
DEFAULT_OUTPUT_PATH = DEFAULT_RESULTS_DIR / "output"
FEASIBLE_STATUSES = [GRB.OPTIMAL, GRB.TIME_LIMIT, GRB.INTERRUPTED, GRB.SUBOPTIMAL]

# SIMON32/64, 17 rounds
WORD_SIZE = 16
KEY_WORDS = 4
R_IN = 2
R_DIST = 8
R_OUT = 7
FIXED_A = [9, 11, 17, 19, 24, 25, 26]
FIXED_B = [14]

# SIMON64/128, 20 rounds
# WORD_SIZE = 32
# KEY_WORDS = 4
# R_IN = 3
# R_DIST = 10
# R_OUT = 7
# FIXED_A = [13, 27, 37, 43, 44, 51, 57, 58]
# FIXED_B = [15]


def Reset_model(N, M, r_in, r_dist, r_out):
    assert r_in >= 1 and r_dist >= 4 and r_out >= 2, "Invalid parameters!"

    global SIMON
    SIMON = Model(f"SIMON{N}/{M*N}")

    global n, m
    n = N  # word size = n bits
    m = M  # key size = m * n bits

    global state_X, state_Y, state_Z1, state_Z2, state_Z
    state_X = {}  # VAR X
    state_Y = {}  # VAR Y
    state_Z1, state_Z2, state_Z = {}, {}, {}  # VAR Z

    global state_DM, state_VM, Key_M
    global state_DW, state_VW, Key_W
    # Active, Value needed, Key bits involved in E_0
    state_DM, state_VM, Key_M = {}, {}, {}
    # Difference needed, Value needed, Key bits involved in E_2
    state_DW, state_VW, Key_W = {}, {}, {}

    global state_D1V, state_D1D, Key_D1, state_D2V, state_D2D, Key_D2
    global state_HD, state_HV, state_RHD, state_RHV, Key_H
    state_D1V, state_D1D, Key_D1 = {}, {}, {}
    state_D2V, state_D2D, Key_D2 = {}, {}, {}
    state_HD, state_HV = {}, {}
    state_RHD, state_RHV, Key_H = {}, {}, {}

    global Key_guess, Key_bridge, Key_space, Key_M_count, Key_W_count
    Key_guess = LinExpr()
    Key_bridge = LinExpr()
    Key_space = LinExpr()
    Key_M_count = SIMON.addVar(vtype=GRB.INTEGER, name="Key_M_count")
    Key_W_count = SIMON.addVar(vtype=GRB.INTEGER, name="Key_W_count")

    global Kg_space, Kg1_space, Kg2_space
    global Kbg1_space, Kbg2_space, Kt, St
    global Svt, Sdt, Svt1, Sdt1
    Kg_space = LinExpr()
    Kg1_space = LinExpr()
    Kg2_space = LinExpr()
    Kbg1_space = LinExpr()
    Kbg2_space = LinExpr()
    Kt = LinExpr()
    St = SIMON.addVar(vtype=GRB.INTEGER, name="St")
    Svt = LinExpr()
    Sdt = LinExpr()
    Svt1 = LinExpr()
    Sdt1 = LinExpr()

    global Deg, Size_A, Size_B, Dist_range
    Deg = LinExpr()
    Size_A = LinExpr()
    Size_B = LinExpr()
    Dist_range = LinExpr()

    global Obj_offline, Obj_online, Obj_time, Obj_data, Plain_diff, Val_Con
    Obj_offline = SIMON.addVar(vtype=GRB.INTEGER, name="Obj_offline")
    Obj_online = SIMON.addVar(vtype=GRB.INTEGER, name="Obj_online")
    Obj_time = SIMON.addVar(vtype=GRB.INTEGER, name="Obj_time")
    Obj_data = SIMON.addVar(vtype=GRB.INTEGER, name="Obj_data")
    Plain_diff = SIMON.addVar(vtype=GRB.INTEGER, name="Plain_diff")
    Val_Con = SIMON.addVar(vtype=GRB.INTEGER, name="Val_Con")


def Evaluate_expr(expr):
    val = 0
    if isinstance(expr, QuadExpr):
        for i in range(expr.size()):
            val += (
                expr.getCoeff(i) * round(expr.getVar1(i).Xn) * round(expr.getVar2(i).Xn)
            )
        val += Evaluate_expr(expr.getLinExpr())
    elif isinstance(expr, LinExpr):
        for i in range(expr.size()):
            val += expr.getCoeff(i) * round(expr.getVar(i).Xn)
        val += expr.getConstant()
    elif isinstance(expr, Var):
        val = round(expr.Xn)
    elif isinstance(expr, (int, float)):
        val = expr
    return val


# ================= Offline Phase =================
def Build_distinguisher(n, r_in, r_dist):
    global Dist_range

    # VAR X - forward differential
    state_X[r_in] = SIMON.addVars(2 * n, vtype=GRB.BINARY, name=f"state_X_{r_in}")
    for rd in range(r_in + 1, r_in + r_dist + 1):
        state_X[rd] = SIMON.addVars(2 * n, vtype=GRB.BINARY, name=f"state_X_{rd}")
        for i in range(n):
            SIMON.addConstr(state_X[rd][i + n] == state_X[rd - 1][i])
            SIMON.addGenConstrOr(
                state_X[rd][i],
                [
                    state_X[rd - 1][i + n],
                    state_X[rd - 1][(i + 1) % n],
                    state_X[rd - 1][(i + 8) % n],
                    state_X[rd - 1][(i + 2) % n],
                ],
            )

    # VAR Y - backward determination
    state_Y[r_in + r_dist] = SIMON.addVars(
        2 * n, vtype=GRB.BINARY, name=f"state_Y_{r_in + r_dist}"
    )
    for rd in range(r_in + r_dist - 1, r_in - 1, -1):
        state_Y[rd] = SIMON.addVars(2 * n, vtype=GRB.BINARY, name=f"state_Y_{rd}")
        for i in range(n):
            SIMON.addConstr(state_Y[rd][i + n] == state_Y[rd + 1][i])
            SIMON.addGenConstrOr(
                state_Y[rd][i],
                [
                    state_Y[rd + 1][i + n],
                    state_Y[rd + 1][(i + n - 1) % n],
                    state_Y[rd + 1][(i + n - 8) % n],
                    state_Y[rd + 1][(i + n - 2) % n],
                ],
            )

    # VAR Z
    for rd in range(r_in, r_in + r_dist):
        state_Z1[rd] = SIMON.addVars(n, vtype=GRB.BINARY, name=f"state_Z1_{rd}")
        state_Z2[rd] = SIMON.addVars(n, vtype=GRB.BINARY, name=f"state_Z2_{rd}")
        state_Z[rd] = SIMON.addVars(n, vtype=GRB.BINARY, name=f"state_Z_{rd}")
        for i in range(n):
            SIMON.addGenConstrAnd(
                state_Z1[rd][i], [state_X[rd][(i + 8) % n], state_Y[rd + 1][i]]
            )
            SIMON.addGenConstrAnd(
                state_Z2[rd][i], [state_X[rd][(i + 1) % n], state_Y[rd + 1][i]]
            )
            SIMON.addGenConstrOr(
                state_Z[rd][i],
                [state_Z1[rd][(i + n - 1) % n], state_Z2[rd][(i + n - 8) % n]],
            )
        Deg.add(quicksum(state_Z[rd][i] for i in range(n)))

    # Validity of the distinguisher
    A = quicksum(state_X[r_in][i] for i in range(2 * n))
    B = quicksum(state_Y[r_in + r_dist][i] for i in range(2 * n))
    Size_A.add(A)
    Size_B.add(B)
    numA = SIMON.addVars(2 * n + 1, vtype=GRB.BINARY)
    lhs = quicksum(numA[i] for i in range(2 * n + 1))
    SIMON.addConstr(lhs == 1)
    SIMON.addConstr(A == quicksum(numA[i] * i for i in range(2 * n + 1)))
    SIMON.addConstr(
        Deg + B + 1 <= quicksum(2**i * numA[i] * B for i in range(2 * n + 1))
    )
    Dist_range = quicksum((2**i - 1) * numA[i] * B for i in range(2 * n + 1)) - Deg

    # Nontrivial
    SIMON.addConstr(A >= 1)
    SIMON.addConstr(B >= 1)


# ================= Online Phase =================
def Build_key_recovery(r_in, r_dist, r_out):
    # E_0
    state_DM[r_in] = SIMON.addVars(2 * n, vtype=GRB.BINARY, name=f"state_DM_{r_in}")
    SIMON.addConstrs(state_DM[r_in][i] == state_X[r_in][i] for i in range(2 * n))
    state_VM[r_in] = SIMON.addVars(2 * n, vtype=GRB.BINARY, name=f"state_VM_{r_in}")
    SIMON.addConstrs(state_VM[r_in][i] == 0 for i in range(2 * n))
    for rd in range(r_in - 1, -1, -1):
        state_DM[rd] = SIMON.addVars(2 * n, vtype=GRB.BINARY, name=f"state_DM_{rd}")
        state_VM[rd] = SIMON.addVars(2 * n, vtype=GRB.BINARY, name=f"state_VM_{rd}")
        Key_M[rd] = SIMON.addVars(n, vtype=GRB.BINARY, name=f"Key_M_{rd}")
        for i in range(n):
            SIMON.addConstr(state_DM[rd][i] == state_DM[rd + 1][i + n])
            SIMON.addGenConstrOr(
                state_DM[rd][i + n],
                [
                    state_DM[rd + 1][i],
                    state_DM[rd][(i + 1) % n],
                    state_DM[rd][(i + 8) % n],
                    state_DM[rd][(i + 2) % n],
                ],
            )
            SIMON.addConstr(Key_M[rd][i] == state_VM[rd + 1][i])
            SIMON.addConstr(state_VM[rd][i + n] == state_VM[rd + 1][i])
            SIMON.addGenConstrOr(
                state_VM[rd][i],
                [
                    state_VM[rd + 1][i + n],
                    state_DM[rd][(i + 7) % n],
                    state_DM[rd][(i + n - 7) % n],
                    state_VM[rd + 1][(i + n - 2) % n],
                    state_VM[rd + 1][(i + n - 1) % n],
                    state_VM[rd + 1][(i + n - 8) % n],
                ],
            )
    SIMON.addConstr(Plain_diff == quicksum(state_DM[0][i] for i in range(2 * n)))

    # E_2
    state_DW[r_in + r_dist] = SIMON.addVars(
        2 * n, vtype=GRB.BINARY, name=f"state_DW_{r_in + r_dist}"
    )
    SIMON.addConstrs(
        state_DW[r_in + r_dist][i] == state_Y[r_in + r_dist][i] for i in range(2 * n)
    )
    state_VW[r_in + r_dist] = SIMON.addVars(
        2 * n, vtype=GRB.BINARY, name=f"state_VW_{r_in + r_dist}"
    )
    SIMON.addConstrs(state_VW[r_in + r_dist][i] == 0 for i in range(2 * n))
    for rd in range(r_in + r_dist + 1, r_in + r_dist + r_out + 1):
        state_DW[rd] = SIMON.addVars(2 * n, vtype=GRB.BINARY, name=f"state_DW_{rd}")
        state_VW[rd] = SIMON.addVars(2 * n, vtype=GRB.BINARY, name=f"state_VW_{rd}")
        Key_W[rd - 1] = SIMON.addVars(n, vtype=GRB.BINARY, name=f"Key_W_{rd-1}")
        for i in range(n):
            SIMON.addConstr(state_DW[rd][i] == state_DW[rd - 1][i + n])
            SIMON.addGenConstrOr(
                state_DW[rd][i + n],
                [
                    state_DW[rd - 1][i],
                    state_DW[rd][(i + n - 1) % n],
                    state_DW[rd][(i + n - 8) % n],
                    state_DW[rd][(i + n - 2) % n],
                ],
            )
            SIMON.addConstr(Key_W[rd - 1][i] == state_VW[rd - 1][i + n])
            SIMON.addConstr(state_VW[rd][i] == state_VW[rd - 1][i + n])
            SIMON.addGenConstrOr(
                state_VW[rd][i + n],
                [
                    state_VW[rd - 1][i],
                    state_VW[rd][(i + n - 1) % n],
                    state_VW[rd][(i + n - 8) % n],
                    state_VW[rd][(i + n - 2) % n],
                    state_DW[rd][(i + n - 1) % n],
                    state_DW[rd][(i + n - 8) % n],
                ],
            )


def Build_key_counting():
    SIMON.addConstr(
        Key_M_count == quicksum(Key_M[rd][i] for rd in Key_M for i in range(n))
    )
    SIMON.addConstr(
        Key_W_count == quicksum(Key_W[rd][i] for rd in Key_W for i in range(n))
    )
    Key_guess.add(Key_M_count)
    Key_guess.add(Key_W_count)
    Key_space.add(Key_guess)


def _count_round_keys(round_keys):
    return quicksum(round_keys[rd][i] for rd in round_keys for i in range(n))


def Build_backward_known_dv(r_in, r_dist, r_out, num, state_DV, state_DD, Key_D):
    """Step 1: values/differences derived after guessing pre-lookup keys."""
    boundary = r_in + r_dist
    r_total = boundary + r_out
    state_DV[r_total] = SIMON.addVars(
        2 * n, vtype=GRB.BINARY, name=f"state_D{num}V_{r_total}"
    )
    state_DD[r_total] = SIMON.addVars(
        2 * n, vtype=GRB.BINARY, name=f"state_D{num}D_{r_total}"
    )
    SIMON.addConstrs(state_DV[r_total][i] == 1 for i in range(2 * n))
    SIMON.addConstrs(state_DD[r_total][i] == 1 for i in range(2 * n))

    for rd in range(r_total - 1, boundary - 1, -1):
        state_DV[rd] = SIMON.addVars(
            2 * n, vtype=GRB.BINARY, name=f"state_D{num}V_{rd}"
        )
        state_DD[rd] = SIMON.addVars(
            2 * n, vtype=GRB.BINARY, name=f"state_D{num}D_{rd}"
        )
        Key_D[rd] = SIMON.addVars(n, vtype=GRB.BINARY, name=f"Key_D{num}_{rd}")
        for i in range(n):
            # L_rd = R_{rd+1}
            SIMON.addConstr(state_DV[rd][i] == state_DV[rd + 1][i + n])
            SIMON.addGenConstrAnd(
                state_DV[rd][i + n],
                [
                    state_DV[rd + 1][i],
                    state_DV[rd][(i + 1) % n],
                    state_DV[rd][(i + 8) % n],
                    state_DV[rd][(i + 2) % n],
                    Key_D[rd][i],
                ],
            )
        for i in range(n):
            SIMON.addConstr(state_DD[rd][i] == state_DD[rd + 1][i + n])
        for i in range(n):
            SIMON.addGenConstrAnd(
                state_DD[rd][i + n],
                [
                    state_DD[rd + 1][i],
                    state_DD[rd][(i + 1) % n],
                    state_DD[rd][(i + 8) % n],
                    state_DD[rd][(i + 2) % n],
                    state_DV[rd][(i + 1) % n],
                    state_DV[rd][(i + 8) % n],
                ],
            )


def _build_not_d_and_h(residual, h_state, d1_state, d2_state, rd):
    """residual = H and not (D1 or D2), bitwise."""
    for i in range(2 * n):
        SIMON.addConstr(residual[rd][i] <= h_state[rd][i])
        SIMON.addConstr(residual[rd][i] <= 1 - d1_state[rd][i])
        SIMON.addConstr(residual[rd][i] <= 1 - d2_state[rd][i])
        SIMON.addConstr(
            residual[rd][i] >= h_state[rd][i] - d1_state[rd][i] - d2_state[rd][i]
        )


def _add_table_size_terms(
    hv_state, d1v_state, d2v_state, hd_state, d1d_state, d2d_state, rd
):
    for i in range(2 * n):
        dv_union = SIMON.addVar(vtype=GRB.BINARY, name=f"state_DVU_{rd}_{i}")
        SIMON.addGenConstrOr(dv_union, [d1v_state[rd][i], d2v_state[rd][i]])

        svt_bit = SIMON.addVar(vtype=GRB.BINARY, name=f"state_Svt_{rd}_{i}")
        SIMON.addGenConstrAnd(svt_bit, [hv_state[rd][i], dv_union])
        Svt.add(svt_bit)

        svt1_bit = SIMON.addVar(vtype=GRB.BINARY, name=f"state_Svt1_{rd}_{i}")
        SIMON.addGenConstrAnd(svt1_bit, [hv_state[rd][i], d1v_state[rd][i]])
        Svt1.add(svt1_bit)

        dd_union = SIMON.addVar(vtype=GRB.BINARY, name=f"state_DDU_{rd}_{i}")
        SIMON.addGenConstrOr(dd_union, [d1d_state[rd][i], d2d_state[rd][i]])

        sdt_bit = SIMON.addVar(vtype=GRB.BINARY, name=f"state_Sdt_{rd}_{i}")
        SIMON.addGenConstrAnd(sdt_bit, [hd_state[rd][i], dd_union])
        Sdt.add(sdt_bit)

        sdt1_bit = SIMON.addVar(vtype=GRB.BINARY, name=f"state_Sdt1_{rd}_{i}")
        SIMON.addGenConstrAnd(sdt1_bit, [hd_state[rd][i], d1d_state[rd][i]])
        Sdt1.add(sdt1_bit)


def Build_forward_table_hd_hv(
    r_in, r_dist, r_out, state_D1V, state_D2V, state_D1D, state_D2D
):
    """Step 2: table-side propagation for SIMON."""
    boundary = r_in + r_dist
    r_total = boundary + r_out

    state_HD[boundary] = SIMON.addVars(
        2 * n, vtype=GRB.BINARY, name=f"state_HD_{boundary}"
    )
    state_HV[boundary] = SIMON.addVars(
        2 * n, vtype=GRB.BINARY, name=f"state_HV_{boundary}"
    )
    SIMON.addConstrs(
        state_HD[boundary][i] == state_DW[boundary][i] for i in range(2 * n)
    )
    SIMON.addConstrs(
        state_HV[boundary][i] == state_VW[boundary][i] for i in range(2 * n)
    )

    for rd in range(boundary, r_total):
        state_RHD[rd] = SIMON.addVars(2 * n, vtype=GRB.BINARY, name=f"state_RHD_{rd}")
        state_RHV[rd] = SIMON.addVars(2 * n, vtype=GRB.BINARY, name=f"state_RHV_{rd}")
        _build_not_d_and_h(state_RHD, state_HD, state_D1D, state_D2D, rd)
        _build_not_d_and_h(state_RHV, state_HV, state_D1V, state_D2V, rd)

        _add_table_size_terms(
            state_HV, state_D1V, state_D2V, state_HD, state_D1D, state_D2D, rd
        )

        state_HD[rd + 1] = SIMON.addVars(
            2 * n, vtype=GRB.BINARY, name=f"state_HD_{rd + 1}"
        )
        state_HV[rd + 1] = SIMON.addVars(
            2 * n, vtype=GRB.BINARY, name=f"state_HV_{rd + 1}"
        )
        Key_H[rd] = SIMON.addVars(n, vtype=GRB.BINARY, name=f"Key_H_{rd}")

        for i in range(n):
            SIMON.addConstr(state_HD[rd + 1][i] == state_RHD[rd][i + n])
            SIMON.addGenConstrOr(
                state_HD[rd + 1][i + n],
                [
                    state_RHD[rd][i],
                    state_HD[rd + 1][(i + n - 1) % n],
                    state_HD[rd + 1][(i + n - 8) % n],
                    state_HD[rd + 1][(i + n - 2) % n],
                ],
            )

            SIMON.addConstr(state_HV[rd + 1][i] == state_RHV[rd][i + n])
            SIMON.addGenConstrOr(
                state_HV[rd + 1][i + n],
                [
                    state_RHV[rd][i],
                    state_HV[rd + 1][(i + n - 1) % n],
                    state_HV[rd + 1][(i + n - 8) % n],
                    state_HV[rd + 1][(i + n - 2) % n],
                    state_HD[rd + 1][(i + n - 1) % n],
                    state_HD[rd + 1][(i + n - 8) % n],
                ],
            )

            SIMON.addConstr(Key_H[rd][i] == state_RHV[rd][i + n])
            Kt.add(Key_H[rd][i])


def Build_value_constraints(
    r_in,
    r_dist,
    r_out,
    min_value_constraints=None,
    max_value_constraints=None,
):
    """Step 3: connect Kg / Kt / St and add value-only lookup constraints."""
    if min_value_constraints is not None:
        SIMON.addConstr(Val_Con >= min_value_constraints)
    if max_value_constraints is not None:
        SIMON.addConstr(Val_Con <= max_value_constraints)

    Build_backward_known_dv(r_in, r_dist, r_out, 1, state_D1V, state_D1D, Key_D1)
    Build_backward_known_dv(r_in, r_dist, r_out, 2, state_D2V, state_D2D, Key_D2)
    Build_forward_table_hd_hv(
        r_in, r_dist, r_out, state_D1V, state_D2V, state_D1D, state_D2D
    )

    Kg1_space.add(_count_round_keys(Key_D1))
    Kg2_space.add(_count_round_keys(Key_D2))
    Kg_space.add(Kg1_space)
    Kg_space.add(Kg2_space)

    Kbg1_space.add(Key_M_count)
    Kbg1_space.add(_count_round_keys(Key_D1))
    Kbg2_space.add(Key_M_count)
    Kbg2_space.add(_count_round_keys(Key_D2))

    SIMON.addConstr(St == Val_Con * Sdt + Svt + Kt)
    SIMON.addConstr(St <= Deg)


# Objective function
def Set_objective(use_value_constraints=False):
    if use_value_constraints:
        SIMON.addConstr(Obj_offline == Deg - Val_Con + Size_A)
    else:
        SIMON.addConstr(Val_Con == 0)
        SIMON.addConstr(Obj_offline == Deg + Size_A)
    SIMON.addConstr(Obj_online >= Key_space + Size_A)
    SIMON.addConstr(Obj_online >= m * n - Dist_range)
    if use_value_constraints:
        SIMON.addConstr(Obj_online >= Val_Con * Sdt + Svt + Kt - Val_Con + 3)
        SIMON.addConstr(
            Obj_online >= Kbg1_space + (Sdt - Sdt1) * Val_Con + (Svt - Svt1) + Kt + 3
        )
        SIMON.addConstr(Obj_online >= Val_Con + Kbg1_space + 3)
        SIMON.addConstr(Obj_online >= Val_Con + Kbg2_space + 3)

    SIMON.addConstr(Obj_time >= Obj_online)
    SIMON.addConstr(2 * Obj_time >= Obj_offline + Obj_online)

    SIMON.addConstr(Obj_data >= Plain_diff)
    SIMON.addConstr(2 * Obj_data >= 2 * Size_A + Obj_offline - Obj_online)

    if use_value_constraints:
        SIMON.setObjectiveN(Obj_time, index=0, priority=6, name="Obj_time", weight=1.0)
        SIMON.setObjectiveN(Obj_data, index=1, priority=5, name="Obj_data", weight=1.0)
        SIMON.setObjectiveN(St, index=2, priority=4, name="Obj_St", weight=1.0)
        SIMON.setObjectiveN(Kt, index=3, priority=3, name="Obj_Kt", weight=1.0)
        SIMON.setObjectiveN(Kg1_space, index=4, priority=2, name="Obj_Kg1", weight=1.0)
        SIMON.setObjectiveN(Kg2_space, index=5, priority=1, name="Obj_Kg2", weight=1.0)
    else:
        SIMON.setObjectiveN(Obj_time, index=0, priority=2, name="Obj_time", weight=1.0)
        SIMON.setObjectiveN(Obj_data, index=1, priority=1, name="Obj_data", weight=1.0)
    SIMON.ModelSense = GRB.MINIMIZE


def Set_objective_time_limit(obj_time_limit):
    if obj_time_limit is None:
        return
    assert obj_time_limit > 0, "Objective time limit must be positive!"

    SIMON.update()
    for obj_index in range(SIMON.NumObj):
        SIMON.getMultiobjEnv(obj_index).setParam("TimeLimit", obj_time_limit)


def Add_extra_constraints(r_in, r_dist, fixed_a=None, fixed_b=None):
    """Fix active bit positions of A=state_X[r_in] and B=state_Y[r_in+r_dist]."""
    if fixed_a is not None:
        fixed_a = set(fixed_a)
        assert all(0 <= i < 2 * n for i in fixed_a), "Invalid A bit index!"
        SIMON.addConstrs(
            state_X[r_in][i] == (1 if i in fixed_a else 0) for i in range(2 * n)
        )

    if fixed_b is not None:
        fixed_b = set(fixed_b)
        assert all(0 <= i < 2 * n for i in fixed_b), "Invalid B bit index!"
        SIMON.addConstrs(
            state_Y[r_in + r_dist][i] == (1 if i in fixed_b else 0)
            for i in range(2 * n)
        )


def _values(var_dict, width):
    return [round(var_dict[i].Xn) for i in range(width)]


def _active_indices(var_dict, width):
    return [i for i in range(width) if round(var_dict[i].Xn) == 1]


def Print_state(name, states, width):
    print(f"---------- Var {name} ----------")
    for rd in sorted(states):
        vals = _values(states[rd], width)
        print(f"{name}[{rd}]")
        if width == 2 * n:
            print("L:", " ".join(map(str, vals[:n])))
            print("R:", " ".join(map(str, vals[n:])))
        else:
            print(" ".join(map(str, vals)))
        print("")


def Print_result(r_in, r_dist, r_out, output_path=None):
    lines = []

    def emit(text=""):
        lines.append(text)

    emit(f"Model Status: {SIMON.Status}")
    if SIMON.SolCount:
        emit("Min_Obj: %g" % SIMON.ObjVal)
        for k in range(SIMON.SolCount):
            SIMON.Params.SolutionNumber = k
            emit(
                "******** Sol no.{}    Time = {} (Offline = {}, Online = {})    Data = {} ********".format(
                    k + 1,
                    round(Obj_time.Xn),
                    round(Obj_offline.Xn),
                    round(Obj_online.Xn),
                    round(Obj_data.Xn),
                )
            )
            emit(
                "******** Deg = {}    Size_A = {}    Plain_diff = {} ********".format(
                    round(Evaluate_expr(Deg)),
                    round(Evaluate_expr(Size_A)),
                    round(Plain_diff.Xn),
                )
            )
            emit(
                "******** Size_B = {}    Dist_range = {} ********".format(
                    round(Evaluate_expr(Size_B)),
                    round(Evaluate_expr(Dist_range)),
                )
            )
            emit(
                "******** Key_guess = {}    Key_bridge = {}    Key_space = {}    Key_M_count = {}    Key_W_count = {} ********".format(
                    round(Evaluate_expr(Key_guess)),
                    round(Evaluate_expr(Key_bridge)),
                    round(Evaluate_expr(Key_space)),
                    round(Key_M_count.Xn),
                    round(Key_W_count.Xn),
                )
            )
            emit(
                "******** Val_Con = {}    Kg = {}    Kg1 = {}    Kg2 = {}    Kbg1 = {}    Kbg2 = {}    Kt = {}    St = {} ********".format(
                    round(Val_Con.Xn),
                    round(Evaluate_expr(Kg_space)),
                    round(Evaluate_expr(Kg1_space)),
                    round(Evaluate_expr(Kg2_space)),
                    round(Evaluate_expr(Kbg1_space)),
                    round(Evaluate_expr(Kbg2_space)),
                    round(Evaluate_expr(Kt)),
                    round(Evaluate_expr(St)),
                )
            )
            emit(
                "******** Svt = {}    Sdt = {}    Svt1 = {}    Sdt1 = {} ********".format(
                    round(Evaluate_expr(Svt)),
                    round(Evaluate_expr(Sdt)),
                    round(Evaluate_expr(Svt1)),
                    round(Evaluate_expr(Sdt1)),
                )
            )
            emit(
                "A = [ {} ]    B = [ {} ]".format(
                    " ".join(map(str, _active_indices(state_X[r_in], 2 * n))),
                    " ".join(map(str, _active_indices(state_Y[r_in + r_dist], 2 * n))),
                )
            )

    text = "\n".join(lines)
    print(text)

    if output_path is not None:
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with output_path.open("w", encoding="utf-8") as file:
            print(text, file=file)
            if SIMON.SolCount:
                Print_state_to_file(file, "Z", state_Z, n)
                Print_state_to_file(file, "X", state_X, 2 * n)
                Print_state_to_file(file, "Y", state_Y, 2 * n)
                Print_state_to_file(file, "DM", state_DM, 2 * n)
                Print_state_to_file(file, "VM", state_VM, 2 * n)
                Print_state_to_file(file, "DW", state_DW, 2 * n)
                Print_state_to_file(file, "VW", state_VW, 2 * n)
                Print_state_to_file(file, "Key_M", Key_M, n)
                Print_state_to_file(file, "Key_W", Key_W, n)
                if state_D1V:
                    Print_state_to_file(file, "D1V", state_D1V, 2 * n)
                    Print_state_to_file(file, "D1D", state_D1D, 2 * n)
                    Print_state_to_file(file, "D2V", state_D2V, 2 * n)
                    Print_state_to_file(file, "D2D", state_D2D, 2 * n)
                    Print_state_to_file(file, "HD", state_HD, 2 * n)
                    Print_state_to_file(file, "HV", state_HV, 2 * n)
                    Print_state_to_file(file, "RHD", state_RHD, 2 * n)
                    Print_state_to_file(file, "RHV", state_RHV, 2 * n)
                    Print_state_to_file(file, "Key_D1", Key_D1, n)
                    Print_state_to_file(file, "Key_D2", Key_D2, n)
                    Print_state_to_file(file, "Key_H", Key_H, n)


def Print_state_to_file(file, name, states, width):
    print(f"---------- Var {name} ----------", file=file)
    for rd in sorted(states):
        vals = _values(states[rd], width)
        print(f"{name}[{rd}]", file=file)
        if width == 2 * n:
            print("L:", " ".join(map(str, vals[:n])), file=file)
            print("R:", " ".join(map(str, vals[n:])), file=file)
        else:
            print(" ".join(map(str, vals)), file=file)
        print("", file=file)


def Start_solver(r_in, r_dist, r_out, output_path=None):
    SIMON.optimize()
    Print_result(r_in, r_dist, r_out, output_path)


def Search_attack(
    N,
    M,
    r_in,
    r_dist,
    r_out,
    # threads=0,
    time_limit=None,
    obj_time_limit=None,
    pool_solutions=1,
    output_path=None,
    log_path=None,
    verbose_gurobi=False,
    fixed_a=None,
    fixed_b=None,
    max_time_exp=None,
    max_data_exp=None,
    use_value_constraints=False,
    min_value_constraints=None,
    max_value_constraints=None,
):
    Reset_model(N, M, r_in, r_dist, r_out)

    SIMON.Params.OutputFlag = 1 if verbose_gurobi else 0
    SIMON.Params.PoolSolutions = pool_solutions
    # if threads > 0:
    #     SIMON.Params.Threads = threads
    if time_limit is not None:
        SIMON.Params.TimeLimit = time_limit
    if log_path is not None:
        log_path = Path(log_path)
        log_path.parent.mkdir(parents=True, exist_ok=True)
        SIMON.Params.LogFile = str(log_path)

    SIMON.message("=== Solving Target: " + SIMON.ModelName + " ===")
    SIMON.message(f"r_in: {r_in}, r_dist: {r_dist}, r_out: {r_out}")

    Build_distinguisher(n, r_in, r_dist)
    # fixed_a = [17, 19, 21, 23, 25, 31]
    # fixed_b = [16]
    Add_extra_constraints(r_in, r_dist, fixed_a, fixed_b)
    Build_key_recovery(r_in, r_dist, r_out)
    Build_key_counting()
    if use_value_constraints:
        Build_value_constraints(
            r_in,
            r_dist,
            r_out,
            min_value_constraints=min_value_constraints,
            max_value_constraints=max_value_constraints,
        )
    Set_objective(use_value_constraints=use_value_constraints)
    Set_objective_time_limit(obj_time_limit)
    if max_time_exp is not None:
        SIMON.addConstr(Obj_time <= max_time_exp)
    if max_data_exp is not None:
        SIMON.addConstr(Obj_data <= max_data_exp)
    Start_solver(r_in, r_dist, r_out, output_path)
    return SIMON.Status


def Parse_args():
    parser = argparse.ArgumentParser(
        description="Search SIMON DS-MITM distinguishers and key-recovery attacks."
    )
    parser.add_argument("-N", "--word-size", type=int, default=16)
    parser.add_argument("-M", "--key-words", type=int, default=4)
    parser.add_argument("--r-in", type=int, default=2)
    parser.add_argument("--r-dist", type=int, default=4)
    parser.add_argument("--r-out", type=int, default=1)
    # parser.add_argument("--threads", type=int, default=0)
    parser.add_argument("--time-limit", type=float, default=None)
    parser.add_argument("--obj-time-limit", type=float, default=None)
    parser.add_argument("--pool-solutions", type=int, default=1)
    parser.add_argument("--verbose-gurobi", action="store_true")
    parser.add_argument("--output", type=Path, default=None)
    parser.add_argument("--log", type=Path, default=None)
    parser.add_argument("--fix-a", type=int, nargs="*", default=None)
    parser.add_argument("--fix-b", type=int, nargs="*", default=None)
    parser.add_argument("--max-time-exp", type=int, default=None)
    parser.add_argument("--max-data-exp", type=int, default=None)
    parser.add_argument("--value-constraints", action="store_true")
    parser.add_argument("--min-value-constraints", type=int, default=None)
    parser.add_argument("--max-value-constraints", type=int, default=None)
    return parser.parse_args()


def main():
    status = Search_attack(
        WORD_SIZE,
        KEY_WORDS,
        R_IN,
        R_DIST,
        R_OUT,
        output_path=DEFAULT_OUTPUT_PATH,
        verbose_gurobi=True,
        fixed_a=FIXED_A,
        fixed_b=FIXED_B,
        use_value_constraints=True,
    )
    return 0 if status in FEASIBLE_STATUSES else 1


if __name__ == "__main__":
    raise SystemExit(main())
