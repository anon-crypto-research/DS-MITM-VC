# Starts from Round 1

from gurobipy import Model, GRB, quicksum, LinExpr, QuadExpr, Var
from pathlib import Path
import sys

MODULE_DIR = Path(__file__).resolve().parent
DEFAULT_OUTPUT_PATH = MODULE_DIR / "results" / "output"

rho0 = [9, 7, 13, 14, 0, 10, 3, 5, 1, 2, 15, 4, 6, 12, 11, 8]
rho1 = [12, 8, 1, 9, 15, 4, 0, 3, 14, 10, 6, 7, 2, 5, 13, 11]


def Reset_model(r_in, r_dist, r_out):
    assert r_in >= 1 and r_dist >= 4 and r_out >= 1, "Invalid parameters!"

    global Vistrutah
    Vistrutah = Model("Vistrutah-256")

    # Odd rounds: C[r-1] -- (ARK, SB, SR) -> A[r] -- (MC) -> B[r] -- (==) -> C[r]
    # Even rounds: C[r-1] -- (ARK, SB, SR) -> A[r] -- (MC) -> B[r] -- (ML) -> C[r]
    # Last round: C[r-1] -- (ARK, SB, SR) -> A[r] -- (ARK) -> B[r](without MC)
    # Last round: C[r-1] -- (ARK, SB, SR) -> A[r] -- (MC, ARK) -> B[r] (with MC)

    global state_XA, state_XB, state_XC, state_YA, state_YB, state_YC, state_Z
    state_XA, state_XB, state_XC = {}, {}, {}  # VAR X
    state_YA, state_YB, state_YC = {}, {}, {}  # VAR Y
    state_Z = {}  # VAR Z

    global state_MA, state_MB, state_MC, state_WA, state_WB, state_WC
    state_MA, state_MB, state_MC = {}, {}, {}  # VAR M
    state_WA, state_WB, state_WC = {}, {}, {}  # VAR W

    global state_D1A, state_D1B, state_D1C, key_D1
    global state_D2A, state_D2B, state_D2C, key_D2
    global state_HA, state_HB, state_HC, key_H
    state_D1A, state_D1B, state_D1C, key_D1 = {}, {}, {}, {}  # VAR D1
    state_D2A, state_D2B, state_D2C, key_D2 = {}, {}, {}, {}  # VAR D2
    state_HA, state_HB, state_HC, key_H = {}, {}, {}, {}  # VAR H

    global Deg, K_on, K_g1, K_g2, K_t, S_t, S_t1
    Deg = LinExpr()
    K_on = LinExpr()
    K_g1 = LinExpr()
    K_g2 = LinExpr()
    K_t = LinExpr()
    S_t = LinExpr()
    S_t1 = LinExpr()

    global Switch_diff_seq, Switch_H
    Switch_diff_seq = Vistrutah.addVar(vtype=GRB.BINARY, name="Switch_diff_seq")
    Switch_H = Vistrutah.addVar(vtype=GRB.BINARY, name="Switch_H")

    global Obj_offline, Obj_online, Obj_data, Obj_memory, Plain_diff, Obj_time, Val_con
    Obj_offline = Vistrutah.addVar(vtype=GRB.INTEGER, name="Obj_offline")
    Obj_online = Vistrutah.addVar(vtype=GRB.INTEGER, name="Obj_online")
    Obj_data = Vistrutah.addVar(vtype=GRB.INTEGER, name="Obj_data")
    Obj_memory = Vistrutah.addVar(vtype=GRB.INTEGER, name="Obj_memory")
    Plain_diff = Vistrutah.addVar(vtype=GRB.INTEGER, name="Plain_diff")
    Obj_time = Vistrutah.addVar(vtype=GRB.INTEGER, name="Obj_time")
    Val_con = Vistrutah.addVar(vtype=GRB.INTEGER, name="Val_con")


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


# var1 = Or(var2)
def Build_mixcolumn_or(var1, var2):
    for i in range(32):
        Vistrutah.addGenConstrOr(
            var1[i],
            [
                var2[4 * (i // 4)],
                var2[4 * (i // 4) + 1],
                var2[4 * (i // 4) + 2],
                var2[4 * (i // 4) + 3],
            ],
        )


# var1 = And(var2)
def Build_mixcolumn_and(var1, var2):
    for i in range(32):
        Vistrutah.addGenConstrAnd(
            var1[i],
            [
                var2[4 * (i // 4)],
                var2[4 * (i // 4) + 1],
                var2[4 * (i // 4) + 2],
                var2[4 * (i // 4) + 3],
            ],
        )


# var1 = Mixing Layer(var2)
def Build_mixing_layer(rd, var1, var2):
    if rd % 2 == 0:  # Even round
        Vistrutah.addConstrs(var1[i] == var2[i * 2] for i in range(16))
        Vistrutah.addConstrs(var1[i + 16] == var2[i * 2 + 1] for i in range(16))
    else:  # Odd round
        Vistrutah.addConstrs(var1[i] == var2[i] for i in range(32))


def SR(i):
    return (i + 4 * (i % 4)) % 16


# var1 = ShiftRow(var2)
def Build_shiftrow(var1, var2):
    for i in range(16):
        Vistrutah.addConstr(var1[i] == var2[SR(i)])
        Vistrutah.addConstr(var1[i + 16] == var2[SR(i) + 16])


# ================= Offline Phase =================
def Build_distinguisher(r_in, r_dist):
    # VAR X - forward differential
    state_XA[r_in] = Vistrutah.addVars(32, vtype=GRB.BINARY, name=f"state_XA_{r_in}")
    for rd in range(r_in + 1, r_in + r_dist + 1):
        state_XB[rd - 1] = Vistrutah.addVars(
            32, vtype=GRB.BINARY, name=f"state_XB_{rd-1}"
        )
        state_XC[rd - 1] = Vistrutah.addVars(
            32, vtype=GRB.BINARY, name=f"state_XC_{rd-1}"
        )
        state_XA[rd] = Vistrutah.addVars(32, vtype=GRB.BINARY, name=f"state_XA_{rd}")
        # MC
        if rd == r_in + 1:  # GDS-MITM
            dummy_MC_st = Vistrutah.addVars(8, vtype=GRB.BINARY, name="dummy_MC_st")
            for i in range(8):
                Vistrutah.addConstr(
                    quicksum(state_XB[rd - 1][j + 4 * i] for j in range(4))
                    + quicksum(state_XA[rd - 1][j + 4 * i] for j in range(4))
                    >= 5 * dummy_MC_st[i]
                )
                for j in range(4):
                    Vistrutah.addConstr(dummy_MC_st[i] >= state_XB[rd - 1][j + 4 * i])
                    Vistrutah.addConstr(dummy_MC_st[i] >= state_XA[rd - 1][j + 4 * i])
        else:
            Build_mixcolumn_or(state_XB[rd - 1], state_XA[rd - 1])
        # ML
        Build_mixing_layer(rd - 1, state_XC[rd - 1], state_XB[rd - 1])
        # ARK, SB, SR
        Build_shiftrow(state_XA[rd], state_XC[rd - 1])

    # VAR Y - backward determination
    state_YB[r_in + r_dist] = Vistrutah.addVars(
        32, vtype=GRB.BINARY, name=f"state_YB_{r_in+r_dist}"
    )
    for rd in range(r_in + r_dist, r_in, -1):
        state_YA[rd] = Vistrutah.addVars(32, vtype=GRB.BINARY, name=f"state_YA_{rd}")
        state_YC[rd - 1] = Vistrutah.addVars(
            32, vtype=GRB.BINARY, name=f"state_YC_{rd-1}"
        )
        state_YB[rd - 1] = Vistrutah.addVars(
            32, vtype=GRB.BINARY, name=f"state_YB_{rd-1}"
        )
        # MC
        if rd == r_in + r_dist:  # GDS-MITM
            dummy_MC_ed = Vistrutah.addVars(8, vtype=GRB.BINARY, name="dummy_MC_ed")
            for i in range(8):
                Vistrutah.addConstr(
                    quicksum(state_YB[rd][j + 4 * i] for j in range(4))
                    + quicksum(state_YA[rd][j + 4 * i] for j in range(4))
                    >= 5 * dummy_MC_ed[i]
                )
                for j in range(4):
                    Vistrutah.addConstr(dummy_MC_ed[i] >= state_YB[rd][j + 4 * i])
                    Vistrutah.addConstr(dummy_MC_ed[i] >= state_YA[rd][j + 4 * i])
        else:
            Build_mixcolumn_or(state_YA[rd], state_YB[rd])
        # ARK, SB, SR
        Build_shiftrow(state_YA[rd], state_YC[rd - 1])
        # ML
        Build_mixing_layer(rd - 1, state_YC[rd - 1], state_YB[rd - 1])

    # VAR Z
    for rd in range(r_in + 1, r_in + r_dist + 1):
        state_Z[rd] = Vistrutah.addVars(32, vtype=GRB.BINARY, name=f"state_Z_{rd}")
        for i in range(32):
            Vistrutah.addGenConstrAnd(
                state_Z[rd][i], [state_XC[rd - 1][i], state_YC[rd - 1][i]]
            )
            # Add to Deg
            Deg.add(state_Z[rd][i])

    # Nontrivial
    Vistrutah.addConstr(quicksum(state_XA[r_in][i] for i in range(32)) >= 1)
    Vistrutah.addConstr(quicksum(state_YB[r_in + r_dist][i] for i in range(32)) >= 1)


# ================= Key Bridging =================
def Build_key_bridging(r_in, r_dist, r_out, k_in, sp, u_out, K_res, switch_last_mc):
    global k_target, k_initial, u_initial
    k_target = Vistrutah.addVars(32, vtype=GRB.BINARY, name="k_target")
    k_initial = Vistrutah.addVars(32, vtype=GRB.BINARY, name="k_initial")
    K_res.add(quicksum(k_initial[i] for i in range(32)))
    u_initial = {}

    k_even_pos = list(range(16, 32)) + list(range(0, 16))
    r_total = r_in + r_dist + r_out

    def get_ml_idx(r, idx):
        if r % 2 != 0 or r == r_total:
            return idx
        # Inverse of ML
        if idx % 2 == 0:
            return idx // 2
        else:
            return 16 + idx // 2

    for rd in range(r_total + 1):

        def mk_idx(i):
            return k_even_pos[i] if rd % 2 == 0 else i

        # for i in range(4):
        #     for j in range(4):
        #         print(mk_idx(4 * j + i), end=" ")
        #     print("|", end=" ")
        #     for j in range(4, 8):
        #         print(mk_idx(4 * j + i), end=" ")
        #     print("")

        # Initial rounds: [0, r_in-1], or last round with switch_last_mc = 0 (without MC)
        if rd < r_in or (rd == r_total and switch_last_mc == 0):
            ref_state = k_in if rd < r_in else u_out
            Vistrutah.addConstrs(
                k_target[mk_idx(i)] >= ref_state[rd][i] for i in range(32)
            )

        # Round r_in + r_dist (meet-in-the-middle connection)
        elif rd == r_in + r_dist:
            for i in range(32):
                Vistrutah.addGenConstrIndicator(
                    Switch_diff_seq, False, k_target[mk_idx(i)] >= sp[i]
                )

        # Final rounds: [r_in + r_dist + 1, r_total - 1], or last round with switch_last_mc = 1 (with MC)
        elif (r_in + r_dist < rd < r_total) or (rd == r_total and switch_last_mc == 1):
            u_initial[rd] = Vistrutah.addVars(
                32, vtype=GRB.BINARY, name=f"u_initial_{rd}"
            )
            K_res.add(quicksum(u_initial[rd][i] for i in range(32)))

            for i in range(8):
                temp_v = Vistrutah.addVar(vtype=GRB.BINARY, name=f"temp_v_{rd}_{i}")
                k_idxs = [mk_idx(get_ml_idx(rd, 4 * i + j)) for j in range(4)]
                expr = quicksum(u_initial[rd][4 * i + j] for j in range(4)) + quicksum(
                    k_initial[idx] for idx in k_idxs
                )
                Vistrutah.addConstr(expr >= 4 * temp_v)
                Vistrutah.addConstr(expr <= 3 + 5 * temp_v)
                Vistrutah.addConstrs(
                    temp_v + u_initial[rd][4 * i + j] >= u_out[rd][4 * i + j]
                    for j in range(4)
                )
                Vistrutah.addConstrs(
                    temp_v + k_initial[idx] >= k_target[idx] for idx in k_idxs
                )

        if rd % 2 == 0:  # Even round key schedule update
            new_k0 = [k_even_pos[i] for i in range(16)]
            new_k1 = [k_even_pos[16 + i] for i in range(16)]
            for i in range(16):
                k_even_pos[i] = new_k0[rho0[i]]
                k_even_pos[16 + i] = new_k1[rho1[i]]


# ================= Online Phase =================
def Build_key_recovery(r_in, r_dist, r_out, switch_last_mc):
    # VAR M - backward determination
    state_MA[r_in] = Vistrutah.addVars(32, vtype=GRB.BINARY, name=f"state_MA_{r_in}")
    state_MC[r_in - 1] = Vistrutah.addVars(
        32, vtype=GRB.BINARY, name=f"state_MC_{r_in-1}"
    )
    Vistrutah.addConstrs(state_MA[r_in][i] == state_XA[r_in][i] for i in range(32))
    # ARK, SB, SR
    Build_shiftrow(state_MA[r_in], state_MC[r_in - 1])

    for rd in range(r_in - 1, 0, -1):
        state_MB[rd] = Vistrutah.addVars(32, vtype=GRB.BINARY, name=f"state_MB_{rd}")
        state_MA[rd] = Vistrutah.addVars(32, vtype=GRB.BINARY, name=f"state_MA_{rd}")
        state_MC[rd - 1] = Vistrutah.addVars(
            32, vtype=GRB.BINARY, name=f"state_MC_{rd-1}"
        )
        # ML
        Build_mixing_layer(rd, state_MC[rd], state_MB[rd])
        # MC
        Build_mixcolumn_or(state_MA[rd], state_MB[rd])
        # ARK, SB, SR
        Build_shiftrow(state_MA[rd], state_MC[rd - 1])

    # Calculate the number of active bytes in the plaintext
    Vistrutah.addConstr(Plain_diff == quicksum(state_MC[0][i] for i in range(32)))

    # VAR W - forward determination
    state_WB[r_in + r_dist] = Vistrutah.addVars(
        32, vtype=GRB.BINARY, name=f"state_WB_{r_in + r_dist}"
    )
    Vistrutah.addConstrs(
        state_WB[r_in + r_dist][i] == state_YB[r_in + r_dist][i] for i in range(32)
    )
    for rd in range(r_in + r_dist, r_in + r_dist + r_out):
        # ML
        state_WC[rd] = Vistrutah.addVars(32, vtype=GRB.BINARY, name=f"state_WC_{rd}")
        Build_mixing_layer(rd, state_WC[rd], state_WB[rd])
        # ARK, SB, SR
        state_WA[rd + 1] = Vistrutah.addVars(
            32, vtype=GRB.BINARY, name=f"state_WA_{rd+1}"
        )
        Build_shiftrow(state_WA[rd + 1], state_WC[rd])
        # MC
        state_WB[rd + 1] = Vistrutah.addVars(
            32, vtype=GRB.BINARY, name=f"state_WB_{rd+1}"
        )
        if rd + 1 == r_in + r_dist + r_out and switch_last_mc == 0:
            Vistrutah.addConstrs(
                state_WB[rd + 1][i] == state_WA[rd + 1][i] for i in range(32)
            )
        else:
            Build_mixcolumn_or(state_WB[rd + 1], state_WA[rd + 1])

    Build_key_bridging(
        r_in,
        r_dist,
        r_out,
        state_MC,
        state_WC[r_in + r_dist],
        state_WA,
        K_on,
        switch_last_mc,
    )


def Build_backward_guess(r_in, r_dist, r_out, num, state_DA, state_DB, state_DC, key_D):
    r_total = r_in + r_dist + r_out
    state_DB[r_total] = Vistrutah.addVars(
        32, vtype=GRB.BINARY, name=f"state_D{num}B_{r_total}"
    )
    Vistrutah.addConstrs(state_DB[r_total][i] == 1 for i in range(32))
    state_DA[r_total] = Vistrutah.addVars(
        32, vtype=GRB.BINARY, name=f"state_D{num}A_{r_total}"
    )
    key_D[r_total] = Vistrutah.addVars(
        32, vtype=GRB.BINARY, name=f"key_D{num}_{r_total}"
    )
    for i in range(32):
        Vistrutah.addGenConstrAnd(
            state_DA[r_total][i], [state_DB[r_total][i], key_D[r_total][i]]
        )
    for rd in range(r_total - 1, r_in + r_dist, -1):
        state_DC[rd] = Vistrutah.addVars(
            32, vtype=GRB.BINARY, name=f"state_D{num}C_{rd}"
        )
        state_DB[rd] = Vistrutah.addVars(
            32, vtype=GRB.BINARY, name=f"state_D{num}B_{rd}"
        )
        state_DA[rd] = Vistrutah.addVars(
            32, vtype=GRB.BINARY, name=f"state_D{num}A_{rd}"
        )
        key_D[rd] = Vistrutah.addVars(32, vtype=GRB.BINARY, name=f"key_D{num}_{rd}")
        # SR^-1, SB^-1
        Build_shiftrow(state_DA[rd + 1], state_DC[rd])
        # ML^-1
        Build_mixing_layer(rd, state_DC[rd], state_DB[rd])
        # MC^-1, ARK
        temp_MC = Vistrutah.addVars(32, vtype=GRB.BINARY)
        Build_mixcolumn_and(temp_MC, state_DB[rd])
        for i in range(32):
            Vistrutah.addGenConstrAnd(state_DA[rd][i], [temp_MC[i], key_D[rd][i]])
    state_DC[r_in + r_dist] = Vistrutah.addVars(
        32, vtype=GRB.BINARY, name=f"state_D{num}C_{r_in + r_dist}"
    )
    state_DB[r_in + r_dist] = Vistrutah.addVars(
        32, vtype=GRB.BINARY, name=f"state_D{num}B_{r_in + r_dist}"
    )
    key_D[r_in + r_dist] = Vistrutah.addVars(
        32, vtype=GRB.BINARY, name=f"key_D{num}_{r_in + r_dist}"
    )
    temp_SR = Vistrutah.addVars(32, vtype=GRB.BINARY)
    Build_shiftrow(state_DA[r_in + r_dist + 1], temp_SR)
    for i in range(32):
        Vistrutah.addGenConstrAnd(
            state_DC[r_in + r_dist][i], [temp_SR[i], key_D[r_in + r_dist][i]]
        )
    Build_mixing_layer(r_in + r_dist, state_DC[r_in + r_dist], state_DB[r_in + r_dist])


def Build_not_D_and_H(res, D1, D2, H):
    for i in range(32):
        Vistrutah.addConstr(res[i] <= 1 - D1[i])
        Vistrutah.addConstr(res[i] <= 1 - D2[i])
        Vistrutah.addConstr(res[i] <= H[i])
        Vistrutah.addConstr(res[i] >= H[i] - D1[i] - D2[i])


def Calc_table_size(res, res1, H, D1, D2):
    for i in range(32):
        v1 = Vistrutah.addVar(vtype=GRB.BINARY)
        Vistrutah.addGenConstrOr(v1, [D1[i], D2[i]])
        v2 = Vistrutah.addVar(vtype=GRB.BINARY)
        Vistrutah.addGenConstrAnd(v2, [H[i], v1])
        res.add(v2)
        v3 = Vistrutah.addVar(vtype=GRB.BINARY)
        Vistrutah.addGenConstrAnd(v3, [D1[i], H[i]])
        res1.add(v3)


def Build_value_constraint(r_in, r_dist, r_out, switch_last_mc):
    # Limited by the table lookup method of key recovery
    # VAR D1 - guessed part one
    Build_backward_guess(
        r_in, r_dist, r_out, 1, state_D1A, state_D1B, state_D1C, key_D1
    )
    Build_key_bridging(
        r_in,
        r_dist,
        r_out,
        state_MC,
        key_D1[r_in + r_dist],
        key_D1,
        K_g1,
        switch_last_mc,
    )
    # VAR D2 - guessed part two
    Build_backward_guess(
        r_in, r_dist, r_out, 2, state_D2A, state_D2B, state_D2C, key_D2
    )
    Build_key_bridging(
        r_in,
        r_dist,
        r_out,
        state_MC,
        key_D2[r_in + r_dist],
        key_D2,
        K_g2,
        switch_last_mc,
    )
    # VAR H - table part
    state_HB[r_in + r_dist] = Vistrutah.addVars(
        32, vtype=GRB.BINARY, name=f"state_HB_{r_in + r_dist}"
    )
    for i in range(32):
        Vistrutah.addGenConstrIndicator(
            Switch_H, False, state_HB[r_in + r_dist][i] == 0
        )
        Vistrutah.addGenConstrIndicator(
            Switch_H,
            False,
            state_D1B[r_in + r_dist][i] + state_D2B[r_in + r_dist][i]
            >= state_YB[r_in + r_dist][i],
        )
        Vistrutah.addGenConstrIndicator(
            Switch_H, True, state_HB[r_in + r_dist][i] == state_YB[r_in + r_dist][i]
        )
    Calc_table_size(
        S_t,
        S_t1,
        state_HB[r_in + r_dist],
        state_D1B[r_in + r_dist],
        state_D2B[r_in + r_dist],
    )
    # ML
    temp = Vistrutah.addVars(32, vtype=GRB.BINARY)
    Build_not_D_and_H(
        temp,
        state_D1B[r_in + r_dist],
        state_D2B[r_in + r_dist],
        state_HB[r_in + r_dist],
    )
    state_HC[r_in + r_dist] = Vistrutah.addVars(
        32, vtype=GRB.BINARY, name=f"state_HC_{r_in + r_dist}"
    )
    Build_mixing_layer(r_in + r_dist, state_HC[r_in + r_dist], temp)
    Calc_table_size(
        S_t,
        S_t1,
        state_HC[r_in + r_dist],
        state_D1C[r_in + r_dist],
        state_D2C[r_in + r_dist],
    )
    # ARK
    temp = Vistrutah.addVars(32, vtype=GRB.BINARY)
    Build_not_D_and_H(
        temp,
        state_D1C[r_in + r_dist],
        state_D2C[r_in + r_dist],
        state_HC[r_in + r_dist],
    )
    key_H[r_in + r_dist] = Vistrutah.addVars(
        32, vtype=GRB.BINARY, name=f"key_H_{r_in + r_dist}"
    )
    for i in range(32):
        Vistrutah.addConstr(key_H[r_in + r_dist][i] <= temp[i])
        Vistrutah.addConstr(key_H[r_in + r_dist][i] <= 1 - Switch_diff_seq)
        Vistrutah.addConstr(key_H[r_in + r_dist][i] >= temp[i] - Switch_diff_seq)
    # SB, SR
    state_HA[r_in + r_dist + 1] = Vistrutah.addVars(
        32, vtype=GRB.BINARY, name=f"state_HA_{r_in + r_dist + 1}"
    )
    Build_shiftrow(state_HA[r_in + r_dist + 1], temp)
    Calc_table_size(
        S_t,
        S_t1,
        state_HA[r_in + r_dist + 1],
        state_D1A[r_in + r_dist + 1],
        state_D2A[r_in + r_dist + 1],
    )

    r_total = r_in + r_dist + r_out
    for rd in range(r_in + r_dist + 1, r_total):
        # ARK
        temp = Vistrutah.addVars(32, vtype=GRB.BINARY)
        Build_not_D_and_H(temp, state_D1A[rd], state_D2A[rd], state_HA[rd])
        key_H[rd] = Vistrutah.addVars(32, vtype=GRB.BINARY, name=f"key_H_{rd}")
        Vistrutah.addConstrs(key_H[rd][i] == temp[i] for i in range(32))
        # MC
        state_HB[rd] = Vistrutah.addVars(32, vtype=GRB.BINARY, name=f"state_HB_{rd}")
        Build_mixcolumn_or(state_HB[rd], temp)
        Calc_table_size(S_t, S_t1, state_HB[rd], state_D1B[rd], state_D2B[rd])
        # ML
        temp = Vistrutah.addVars(32, vtype=GRB.BINARY)
        Build_not_D_and_H(temp, state_D1B[rd], state_D2B[rd], state_HB[rd])
        state_HC[rd] = Vistrutah.addVars(32, vtype=GRB.BINARY, name=f"state_HC_{rd}")
        Build_mixing_layer(rd, state_HC[rd], temp)
        Calc_table_size(S_t, S_t1, state_HC[rd], state_D1C[rd], state_D2C[rd])
        # SB, SR
        temp = Vistrutah.addVars(32, vtype=GRB.BINARY)
        Build_not_D_and_H(temp, state_D1C[rd], state_D2C[rd], state_HC[rd])
        state_HA[rd + 1] = Vistrutah.addVars(
            32, vtype=GRB.BINARY, name=f"state_HA_{rd+1}"
        )
        Build_shiftrow(state_HA[rd + 1], temp)
        Calc_table_size(
            S_t, S_t1, state_HA[rd + 1], state_D1A[rd + 1], state_D2A[rd + 1]
        )
    temp = Vistrutah.addVars(32, vtype=GRB.BINARY)
    Build_not_D_and_H(temp, state_D1A[r_total], state_D2A[r_total], state_HA[r_total])
    key_H[r_total] = Vistrutah.addVars(32, vtype=GRB.BINARY, name=f"key_H_{r_total}")
    Vistrutah.addConstrs(key_H[r_total][i] == temp[i] for i in range(32))
    # Calculate K_t
    for rd in range(r_in + r_dist, r_total + 1):
        K_t.add(quicksum(key_H[rd][i] for i in range(32)))
    # Memory Limit
    Vistrutah.addConstr((Val_con + Switch_diff_seq) * S_t + K_t <= Deg)


# Objective function
def Set_objective():
    # Obj_offline
    Vistrutah.addConstr(Obj_offline == 8 * (Deg - Val_con))
    # Obj_Online
    Vistrutah.addConstr(Obj_online >= 8 * (Val_con + K_g1))
    Vistrutah.addConstr(
        Obj_online >= 8 * (K_g1 + (Val_con + Switch_diff_seq) * (S_t - S_t1) + K_t)
    )
    Vistrutah.addConstr(Obj_online >= 8 * (Val_con + K_g2))
    Vistrutah.addConstr(Obj_online >= 8 * K_on)
    # Obj_time
    Vistrutah.addConstr(Obj_time >= Obj_online)
    Vistrutah.addConstr(2 * Obj_time >= Obj_offline + Obj_online)
    Vistrutah.addConstr(
        Obj_time >= 8 * ((Val_con + Switch_diff_seq) * S_t + K_t - Val_con)
    )

    # Obj_Data: Consider DTM trade-off
    Vistrutah.addConstr(Obj_data >= 8 * Plain_diff)
    Vistrutah.addConstr(2 * Obj_data >= 16 + Obj_offline - Obj_online)
    # Obj_Memory
    Vistrutah.addGenConstrMin(Obj_memory, [Obj_offline, Obj_time])

    # Objective
    Vistrutah.setObjectiveN(Obj_time, index=0, priority=4, name="Obj_time", weight=1.0)
    Vistrutah.setObjectiveN(Obj_data, index=1, priority=3, name="Obj_data", weight=1.0)
    Vistrutah.setObjectiveN(K_on, index=2, priority=2, name="Obj_K_on", weight=1.0)
    Vistrutah.setObjectiveN(
        Obj_memory, index=3, priority=1, name="Obj_memory", weight=1.0
    )
    Vistrutah.ModelSense = GRB.MINIMIZE


def print_block(vals):
    for i in range(4):
        for j in range(4):
            print(vals[4 * j + i], end=" ")
        print("|", end=" ")
        for j in range(4, 8):
            print(vals[4 * j + i], end=" ")
        print("")


def print_var(name, dict):
    print(f"---------- Var {name} ----------")
    for k in sorted(dict):
        print(f"{name}[{k}]")
        print_block([round(dict[k][i].Xn) for i in range(32)])
        print("")


def print_flow(name, dict_A, dict_B, dict_C):
    print(f"------- Var {name} -------")
    timeline = (
        [(r, "A", dict_A[r]) for r in dict_A]
        + [(r, "B", dict_B[r]) for r in dict_B]
        + [(r, "C", dict_C[r]) for r in dict_C]
    )
    timeline.sort(key=lambda x: (x[0], x[1]))
    for idx, (r, stage, var_dict) in enumerate(timeline):
        print(f"{name}{stage}[{r}]")
        print_block([round(var_dict[i].Xn) for i in range(32)])
        if idx < len(timeline) - 1:
            if stage == "A":
                print("-- (MC) ->")
            elif stage == "B":
                print("-- (ML) ->")
            elif stage == "C":
                print("-- (ARK, SB, SR) ->")
    print("")


def Start_solver(r_in, r_dist, r_out, switch_last_mc, Print_result):
    Vistrutah.optimize()
    if not Print_result:
        return

    DEFAULT_OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    file = DEFAULT_OUTPUT_PATH.open("w", encoding="utf-8")
    sys.stdout = file
    print(f"\n================= Switch_last_mc = {switch_last_mc} =================")
    print("Model Status:", Vistrutah.Status)
    if Vistrutah.Status in [2, 9, 11]:
        print("Min_Obj: %g" % Vistrutah.ObjVal)
        for k in range(Vistrutah.SolCount):
            Vistrutah.Params.SolutionNumber = k
            print(
                "******** Complexities: Time = {} (Offline = {}, Online = {})    Data = {}    Memory = {} ********".format(
                    round(Obj_time.Xn),
                    round(Obj_offline.Xn),
                    round(Obj_online.Xn),
                    round(Obj_data.Xn),
                    round(Obj_memory.Xn),
                )
            )
            print(
                "******** Deg = {}    K_on = {}    K_g1 = {}    K_g2 = {}    Val_con = {}    Plaintext Difference = {} ********".format(
                    round(Evaluate_expr(Deg)),
                    round(Evaluate_expr(K_on)),
                    round(Evaluate_expr(K_g1)),
                    round(Evaluate_expr(K_g2)),
                    round(Evaluate_expr(Val_con)),
                    round(Evaluate_expr(Plain_diff)),
                )
            )
            print(
                "******** K_t = {}    S_t = {}    S_t1 = {} ********".format(
                    round(Evaluate_expr(K_t)),
                    round(Evaluate_expr(S_t)),
                    round(Evaluate_expr(S_t1)),
                )
            )
            print(
                "******** Switch_diff_seq = {}    Switch_H = {} ********".format(
                    round(Switch_diff_seq.Xn), round(Switch_H.Xn)
                )
            )

            strr = (
                "A = [ "
                + " ".join(
                    [str(i) for i in range(32) if round(state_XC[r_in][i].Xn) == 1]
                )
                + " ]    "
            )
            strr += (
                "B = [ "
                + " ".join(
                    [
                        str(i)
                        for i in range(32)
                        if round(state_YC[r_in + r_dist - 1][i].Xn) == 1
                    ]
                )
                + " ]"
            )
            print(strr)

            print_var("Z", state_Z)
            print_flow("X", state_XA, state_XB, state_XC)
            print_flow("Y", state_YA, state_YB, state_YC)
            print_flow("M", state_MA, state_MB, state_MC)
            print_flow("W", state_WA, state_WB, state_WC)

            print_flow("D1", state_D1A, state_D1B, state_D1C)
            print_var("key_D1", key_D1)

            print_flow("D2", state_D2A, state_D2B, state_D2C)
            print_var("key_D2", key_D2)

            print_flow("H", state_HA, state_HB, state_HC)
            print_var("key_H", key_H)

            # if u_initial:
            #     print_var("u_initial", u_initial)

            # print("---------- Var k_target ----------")
            # print_block([round(k_target[i].Xn) for i in range(32)])
            # print("\n---------- Var k_initial ----------")
            # print_block([round(k_initial[i].Xn) for i in range(32)])
            # print("")
    sys.stdout = sys.__stdout__
    file.close()


def _fix_active_indices(var_dict, rd, active_indices, name):
    active_set = set(active_indices)
    if len(active_set) != len(active_indices):
        raise ValueError(f"Duplicate indices in {name}: {active_indices}")
    invalid_indices = sorted(i for i in active_set if i < 0 or i >= 32)
    if invalid_indices:
        raise ValueError(f"Invalid indices in {name}: {invalid_indices}")

    for i in range(32):
        Vistrutah.addConstr(
            var_dict[rd][i] == int(i in active_set), name=f"fix_{name}_{rd}_{i}"
        )


def Fix_A_B_extra(A_indices, B_indices):
    def add_constraints(r_in, r_dist, r_out, switch_last_mc):
        _fix_active_indices(state_XC, r_in, A_indices, "A")
        _fix_active_indices(state_YC, r_in + r_dist - 1, B_indices, "B")

    return add_constraints


def Search_attack(r_in, r_dist, r_out, switch_last_mc, extra_constr):
    Reset_model(r_in, r_dist, r_out)
    # Vistrutah.params.OutputFlag = 0
    Vistrutah.Params.PoolSearchMode = 0
    Vistrutah.Params.PoolSolutions = 1
    # Vistrutah.Params.Threads = 192
    Vistrutah.message("=== Solving Target: " + Vistrutah.ModelName + " ===")
    Vistrutah.message(f"=== Parameters [Switch MC={switch_last_mc}] ===")
    Vistrutah.message(f"r_in: {r_in}, r_dist: {r_dist}, r_out: {r_out}")
    Build_distinguisher(r_in, r_dist)
    for add_constraints in extra_constr:
        add_constraints(r_in, r_dist, r_out, switch_last_mc)
    Build_key_recovery(r_in, r_dist, r_out, switch_last_mc)
    Build_value_constraint(r_in, r_dist, r_out, switch_last_mc)
    Set_objective()
    Start_solver(r_in, r_dist, r_out, switch_last_mc, True)


if __name__ == "__main__":
    r_in = 3
    r_dist = 4
    r_out = 3
    switch_last_mc = 0  # 0: without MC in the last round, 1: with MC in the last round
    extra_constr = [Fix_A_B_extra([0, 1, 2], [17, 22, 28])]

    DEFAULT_OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    DEFAULT_OUTPUT_PATH.write_text("", encoding="utf-8")
    Search_attack(r_in, r_dist, r_out, switch_last_mc, extra_constr)
