# MC involved in last round
# Consider difference sequence first
# r_dist >= 2
# Nk fixed to 8

from gurobipy import Model, GRB, quicksum, LinExpr, QuadExpr, Var
from pathlib import Path
import sys

MODULE_DIR = Path(__file__).resolve().parent
DEFAULT_OUTPUT_PATH = MODULE_DIR / "results" / "output"


def Reset_model(nb, nk):
    assert nb in [4, 5, 6, 7, 8] and nk == 8
    global Nb, Nk, Rijndael
    Nb = nb
    Nk = nk
    Rijndael = Model("Rijndael-" + str(32 * Nb) + "-" + str(32 * Nk))
    print("=== Solving Target: " + Rijndael.ModelName + " ===")

    # l -(ARK,SB,SR)-> nl -(MC)-> l
    global state_X_l, state_X_nl, state_Y_l, state_Y_nl, state_Z
    state_X_l = {}  # VAR X: parameters
    state_X_nl = {}
    state_Y_l = {}  # VAR Y: parameters
    state_Y_nl = {}
    state_Z = {}  # VAR Z: parameters

    global state_W_l, state_W_nl, state_M_l, state_M_nl
    state_W_l = {}
    state_W_nl = {}  # [r_in + r_dist + 1, r_in + r_dist + r_out + 1)
    state_M_l = {}  # [0, r_in + 1)
    state_M_nl = {}

    global state_D1_l, state_D1_nl, key_D1, state_D2_l, state_D2_nl, key_D2
    global state_H_l, state_H_nl, key_H
    state_D1_l = {}
    state_D1_nl = {}
    key_D1 = {}
    state_D2_l = {}
    state_D2_nl = {}
    key_D2 = {}
    state_H_l = {}
    state_H_nl = {}
    key_H = {}

    global key_dist, ukey_dist
    key_dist = {}  # [r_in + 1, r_in + r_dist)
    ukey_dist = {}  # [r_in + 1, r_in + r_dist)

    global Deg, Plain_diff, Key, Key_sieve, K_off, K_cup, Kg1, Kg2, Kt, St, St1
    Deg = LinExpr()
    Plain_diff = LinExpr()
    Key = LinExpr()
    Key_sieve = LinExpr()
    K_off = LinExpr()
    K_cup = LinExpr()
    Kg1 = LinExpr()
    Kg2 = LinExpr()
    Kt = LinExpr()
    St = LinExpr()
    St1 = LinExpr()

    global H_off, Val_Con, Obj_Data, Obj_Memory, Obj_Offline, Obj_Online, Obj
    H_off = Rijndael.addVar(vtype=GRB.BINARY, name="H_off")
    Val_Con = Rijndael.addVar(vtype=GRB.INTEGER, name="Val_Con")
    Obj_Data = Rijndael.addVar(vtype=GRB.INTEGER, name="Obj_Data")
    Obj_Memory = Rijndael.addVar(vtype=GRB.INTEGER, name="Obj_Memory")
    Obj_Offline = Rijndael.addVar(vtype=GRB.INTEGER, name="Obj_Offline")
    Obj_Online = Rijndael.addVar(vtype=GRB.INTEGER, name="Obj_Online")
    Obj = Rijndael.addVar(vtype=GRB.INTEGER, name="Obj")


def Evaluate_expr(expr):
    val = 0
    if isinstance(expr, QuadExpr):
        for i in range(expr.size()):
            coeff = expr.getCoeff(i)
            var1 = expr.getVar1(i)
            var2 = expr.getVar2(i)
            val += coeff * round(var1.Xn) * round(var2.Xn)
        val += Evaluate_expr(expr.getLinExpr())
    elif isinstance(expr, LinExpr):
        for i in range(expr.size()):
            coeff = expr.getCoeff(i)
            var = expr.getVar(i)
            val += coeff * round(var.Xn)
        val += expr.getConstant()
    elif isinstance(expr, Var):
        val = round(expr.Xn)
    elif isinstance(expr, (int, float)):
        val = expr
    return val


# var1 = And(var2)
def Build_mixcolumn_and(var1, var2):
    for i in range(4 * Nb):
        Rijndael.addGenConstrAnd(
            var1[i],
            [
                var2[4 * (i // 4)],
                var2[4 * (i // 4) + 1],
                var2[4 * (i // 4) + 2],
                var2[4 * (i // 4) + 3],
            ],
        )


# var1 = Or(var2)
def Build_mixcolumn_or(var1, var2):
    for i in range(4 * Nb):
        Rijndael.addGenConstrOr(
            var1[i],
            [
                var2[4 * (i // 4)],
                var2[4 * (i // 4) + 1],
                var2[4 * (i // 4) + 2],
                var2[4 * (i // 4) + 3],
            ],
        )


def SR(i):
    if Nb <= 6:
        j = (i + 4 * (i % 4)) % (4 * Nb)
    elif Nb == 7:
        if i % 4 == 0:
            j = i
        elif i % 4 == 1:
            j = (i + 4) % 28
        elif i % 4 == 2:
            j = (i + 8) % 28
        else:  # i % 4 == 3
            j = (i + 16) % 28
    else:  # Nb == 8
        if i % 4 == 0:
            j = i
        elif i % 4 == 1:
            j = (i + 4) % 32
        elif i % 4 == 2:
            j = (i + 12) % 32
        else:  # i % 4 == 3
            j = (i + 16) % 32
    return j


# var_nl = ShiftRow(var_l)
def Build_shiftrow(var_nl, var_l):
    for i in range(4 * Nb):
        Rijndael.addConstr(var_nl[i] == var_l[SR(i)])


# Involved key bytes in the distinguisher part
def Extract_involved_key(r_dist):
    for rd in range(2, r_dist):
        key_dist[rd] = Rijndael.addVars(
            4 * Nb, vtype=GRB.BINARY, name="key_dist_" + str(rd)
        )
        ukey_dist[rd] = Rijndael.addVars(
            4 * Nb, vtype=GRB.BINARY, name="ukey_dist_" + str(rd)
        )
        for i in range(Nb):
            Z_all_one = Rijndael.addVar(vtype=GRB.BINARY)
            X_all_one = Rijndael.addVar(vtype=GRB.BINARY)
            X_zero = Rijndael.addVar(vtype=GRB.BINARY)
            j0 = 4 * i
            j1 = 4 * i + 1
            j2 = 4 * i + 2
            j3 = 4 * i + 3
            Rijndael.addGenConstrAnd(
                X_all_one,
                [state_Z[rd][j0], state_Z[rd][j1], state_Z[rd][j2], state_Z[rd][j3]],
            )
            Rijndael.addConstr(X_zero == 1 - X_all_one)
            Rijndael.addGenConstrAnd(
                Z_all_one,
                [
                    state_Z[rd - 1][SR(j0)],
                    state_Z[rd - 1][SR(j1)],
                    state_Z[rd - 1][SR(j2)],
                    state_Z[rd - 1][SR(j3)],
                    X_zero,
                ],
            )
            for j in range(4):
                Rijndael.addGenConstrAnd(
                    key_dist[rd][j + 4 * i], [Z_all_one, state_Z[rd][j + 4 * i]]
                )
                Rijndael.addGenConstrAnd(
                    ukey_dist[rd][j + 4 * i],
                    [X_all_one, state_Z[rd - 1][SR(j + 4 * i)]],
                )
        for i in range(4 * Nb):
            # Add to Key_sieve
            Key_sieve.add(key_dist[rd][i])
            Key_sieve.add(ukey_dist[rd][i])


# Offline phase
def Build_distinguisher(r_dist):
    # VAR X - forward differential
    state_X_nl[0] = Rijndael.addVars(4 * Nb, vtype=GRB.BINARY, name="state_X_nl_0")
    for rd in range(1, r_dist):
        state_X_l[rd] = Rijndael.addVars(
            4 * Nb, vtype=GRB.BINARY, name="state_X_l_" + str(rd)
        )
        state_X_nl[rd] = Rijndael.addVars(
            4 * Nb, vtype=GRB.BINARY, name="state_X_nl_" + str(rd)
        )
        # MC
        if rd == 1:  # GDS-MITM
            dummy_MC_st = Rijndael.addVars(Nb, vtype=GRB.BINARY, name="dummy_MC_st")
            for i in range(Nb):
                Rijndael.addConstr(
                    quicksum(state_X_l[rd][j + 4 * i] for j in range(4))
                    + quicksum(state_X_nl[rd - 1][j + 4 * i] for j in range(4))
                    >= 5 * dummy_MC_st[i]
                )
                for j in range(4):
                    Rijndael.addConstr(dummy_MC_st[i] >= state_X_l[rd][j + 4 * i])
                    Rijndael.addConstr(dummy_MC_st[i] >= state_X_nl[rd - 1][j + 4 * i])
        else:
            Build_mixcolumn_or(state_X_l[rd], state_X_nl[rd - 1])
        # SR
        Build_shiftrow(state_X_nl[rd], state_X_l[rd])

    # VAR Y - backward determination
    state_Y_l[r_dist] = Rijndael.addVars(
        4 * Nb, vtype=GRB.BINARY, name="state_Y_l_" + str(r_dist)
    )
    for rd in range(r_dist - 1, 0, -1):
        state_Y_nl[rd] = Rijndael.addVars(
            4 * Nb, vtype=GRB.BINARY, name="state_Y_nl_" + str(rd)
        )
        state_Y_l[rd] = Rijndael.addVars(
            4 * Nb, vtype=GRB.BINARY, name="state_Y_l_" + str(rd)
        )
        # MC
        if rd == r_dist - 1:  # GDS-MITM
            dummy_MC_ed = Rijndael.addVars(Nb, vtype=GRB.BINARY, name="dummy_MC_ed")
            for i in range(Nb):
                Rijndael.addConstr(
                    quicksum(state_Y_l[rd + 1][j + 4 * i] for j in range(4))
                    + quicksum(state_Y_nl[rd][j + 4 * i] for j in range(4))
                    >= 5 * dummy_MC_ed[i]
                )
                for j in range(4):
                    Rijndael.addConstr(dummy_MC_ed[i] >= state_Y_l[rd + 1][j + 4 * i])
                    Rijndael.addConstr(dummy_MC_ed[i] >= state_Y_nl[rd][j + 4 * i])
        else:
            Build_mixcolumn_or(state_Y_nl[rd], state_Y_l[rd + 1])
        # SR
        Build_shiftrow(state_Y_nl[rd], state_Y_l[rd])

    # VAR Z
    for rd in range(1, r_dist):
        state_Z[rd] = Rijndael.addVars(
            4 * Nb, vtype=GRB.BINARY, name="state_Z_" + str(rd)
        )
        for i in range(4 * Nb):
            Rijndael.addGenConstrAnd(
                state_Z[rd][i], [state_X_l[rd][i], state_Y_l[rd][i]]
            )
            # Add to Deg
            Deg.add(state_Z[rd][i])

    # Nontrivial
    Rijndael.addConstr(quicksum(state_X_nl[0][i] for i in range(4 * Nb)) >= 1)
    Rijndael.addConstr(quicksum(state_Y_l[r_dist][i] for i in range(4 * Nb)) >= 1)

    # Extract involved key bytes in the distinguisher part
    Extract_involved_key(r_dist)


# Online phase
def Build_key_recovery(beta, r_in, r_dist, r_out):
    # VAR M - backward determination
    state_M_nl[r_in] = Rijndael.addVars(
        4 * Nb, vtype=GRB.BINARY, name="state_M_nl_" + str(r_in)
    )
    state_M_l[r_in] = Rijndael.addVars(
        4 * Nb, vtype=GRB.BINARY, name="state_M_l_" + str(r_in)
    )
    Rijndael.addConstrs(state_M_nl[r_in][i] == state_X_nl[0][i] for i in range(4 * Nb))
    # SR
    Build_shiftrow(state_M_nl[r_in], state_M_l[r_in])
    for rd in range(r_in - 1, -1, -1):
        state_M_nl[rd] = Rijndael.addVars(
            4 * Nb, vtype=GRB.BINARY, name="state_M_nl_" + str(rd)
        )
        state_M_l[rd] = Rijndael.addVars(
            4 * Nb, vtype=GRB.BINARY, name="state_M_l_" + str(rd)
        )
        # MC
        Build_mixcolumn_or(state_M_nl[rd], state_M_l[rd + 1])
        # SR
        Build_shiftrow(state_M_nl[rd], state_M_l[rd])

    # Calculate the number of active bytes in the plaintext
    Plain_diff.add(quicksum(state_M_l[0][i] for i in range(4 * Nb)))

    # VAR W - forward determination
    state_W_l[0] = Rijndael.addVars(4 * Nb, vtype=GRB.BINARY, name="state_W_l_0")
    Rijndael.addConstrs(state_W_l[0][i] == state_Y_l[r_dist][i] for i in range(4 * Nb))
    for rd in range(r_out):
        state_W_nl[rd] = Rijndael.addVars(
            4 * Nb, vtype=GRB.BINARY, name="state_W_nl_" + str(rd)
        )
        state_W_l[rd + 1] = Rijndael.addVars(
            4 * Nb, vtype=GRB.BINARY, name="state_W_l_" + str(rd + 1)
        )
        # SR
        Build_shiftrow(state_W_nl[rd], state_W_l[rd])
        # MC
        # if rd < r_out - 1:
        Build_mixcolumn_or(state_W_l[rd + 1], state_W_nl[rd])
        # else:
        #     Rijndael.addConstrs(state_W_l[rd + 1][i] == state_W_nl[rd][i] for i in range(4 * Nb))

    # for rd in range(r_in + 1):
    #     # Multiset method -> r_in
    #     # Sequence method -> r_in + 1
    #     for i in range(4 * Nb):
    #         Key.add(state_M_l[rd][i])

    # for rd in range(r_out):
    #     for i in range(4 * Nb):
    #         Key.add(state_W_nl[rd][i])

    Build_key_bridging(beta, r_in, r_dist, r_out, state_M_l, {}, {}, state_W_nl, Key)


def calc(x, y):
    return 4 * x + y


def Get_idx_k(m, x, y, typ):
    # typ0: len = 1 & (x-8,x-1) -> x
    a, b = -1, -1
    if typ == 0:
        if x - 8 < 0:  # OUT
            return -1, -1
        a = calc(x - 8, y)
        if x % 8 == 0:  # column 0
            b = calc(x - 1, (y + 1) % 4)
        else:  # column 1,2,3,4,5,6,7
            b = calc(x - 1, y)
    # typ1: len = 1 & (x+7,x+8) -> x
    if typ == 1:
        if x + 8 >= m:  # OUT
            return -1, -1
        a = calc(x + 8, y)
        if x % 8 == 0:  # column 0
            b = calc(x + 7, (y + 1) % 4)
        else:  # column 1,2,3,4,5,6,7
            b = calc(x + 7, y)
    # typ2: len = 1 & (x-7,x+1) -> x
    if typ == 2:
        if x + 1 >= m or x - 7 < 0:  # OUT
            return -1, -1
        if x % 8 == 7:  # column 7
            a = calc(x - 7, (y + 3) % 4)
            b = calc(x + 1, (y + 3) % 4)
        else:  # column 0,1,2,3,4,5,6
            a = calc(x - 7, y)
            b = calc(x + 1, y)
    # typ3: len = 2 & (x-16,x-2) -> x
    if typ == 3:
        if x - 16 < 0 or x % 8 == 0 or x % 8 == 4:  # OUT or column 0 or column 4
            return -1, -1
        a = calc(x - 16, y)
        if x % 8 == 1:  # column 1
            b = calc(x - 2, (y + 1) % 4)
        else:  # column 2,3,5,6,7
            b = calc(x - 2, y)
    # typ4: len = 2 & (x+14,x+16) -> x
    if typ == 4:
        if x + 16 >= m or x % 8 == 0 or x % 8 == 4:  # OUT or column 0 or column 4
            return -1, -1
        a = calc(x + 16, y)
        if x % 8 == 1:  # column 1
            b = calc(x + 14, (y + 1) % 4)
        else:  # column 2,3,5,6,7
            b = calc(x + 14, y)
    # typ5: len = 2 & (x-14,x+2) -> x
    if typ == 5:
        if (
            x + 2 >= m or x - 14 < 0 or x % 8 == 6 or x % 8 == 2
        ):  # OUT or column 2 or column 6
            return -1, -1
        if x % 8 == 7:  # column 7
            a = calc(x - 14, (y + 3) % 4)
            b = calc(x + 2, (y + 3) % 4)
        else:  # column 0,1,3,4,5
            a = calc(x - 14, y)
            b = calc(x + 2, y)
    # typ6: len = 4 & (x-32,x-4) -> x
    if typ == 6:
        if x - 32 < 0 or (x % 8 != 3 and x % 8 != 7):  # OUT or column 0,1,2,4,5,6
            return -1, -1
        if x % 8 == 3:  # column 3
            a = calc(x - 32, y)
            b = calc(x - 4, (y + 1) % 4)
        else:  # column 7
            a = calc(x - 32, y)
            b = calc(x - 4, y)
    # typ7: len = 4 & (x+28,x+32) -> x
    if typ == 7:
        if x + 32 >= m or (x % 8 != 3 and x % 8 != 7):  # OUT or column 0,1,2,4,5,6
            return -1, -1
        if x % 8 == 3:  # column 3
            a = calc(x + 28, (y + 1) % 4)
            b = calc(x + 32, y)
        else:  # column 7
            a = calc(x + 28, y)
            b = calc(x + 32, y)
    # typ8: len = 4 & (x-28,x+4) -> x
    if typ == 8:
        if (
            x + 4 >= m or x - 28 < 0 or (x % 8 != 3 and x % 8 != 7)
        ):  # OUT or column 0,1,2,4,5,6
            return -1, -1
        if x % 8 == 7:  # column 7
            a = calc(x - 28, (y + 3) % 4)
            b = calc(x + 4, (y + 3) % 4)
        else:  # column 3
            a = calc(x - 28, y)
            b = calc(x + 4, y)
    return a, b


def Get_idx_u(m, x, y, typ):
    # All relations of u with rotation are invalid
    assert typ <= 5
    a, b = -1, -1
    # typ0: len = 1 & (x-8,x-1) -> x
    if typ == 0:
        if x - 8 < 0 or x % 8 == 0 or x % 8 == 4:  # OUT or column 0 or column 4
            return -1, -1
        a = calc(x - 8, y)
        b = calc(x - 1, y)
    # typ1: len = 1 & (x+7,x+8) -> x
    if typ == 1:
        if x + 8 >= m or x % 8 == 0 or x % 8 == 4:  # OUT or column 0 or column 4
            return -1, -1
        a = calc(x + 7, y)
        b = calc(x + 8, y)
    # typ2: len = 1 & (x-7,x+1) -> x
    if typ == 2:
        if (
            x + 1 >= m or x - 7 < 0 or x % 8 == 7 or x % 8 == 3
        ):  # OUT or column 7 or column 3
            return -1, -1
        a = calc(x - 7, y)
        b = calc(x + 1, y)
    # typ3: len = 2 & (x-16,x-2) -> x
    if typ == 3:
        if (
            x - 16 < 0 or x % 8 <= 1 or x % 8 == 4 or x % 8 == 5
        ):  # OUT or column 0,1,4,5
            return -1, -1
        a = calc(x - 16, y)
        b = calc(x - 2, y)
    # typ4: len = 2 & (x+14,x+16) -> x
    if typ == 4:
        if (
            x + 16 >= m or x % 8 <= 1 or x % 8 == 4 or x % 8 == 5
        ):  # OUT or column 0,1,4,5
            return -1, -1
        a = calc(x + 14, y)
        b = calc(x + 16, y)
    # typ5: len = 2 & (x-14,x+2) -> x
    if typ == 5:
        if (
            x + 2 >= m or x - 14 < 0 or x % 8 >= 6 or x % 8 == 2 or x % 8 == 3
        ):  # OUT or column 2,3,6,7
            return -1, -1
        a = calc(x - 14, y)
        b = calc(x + 2, y)
    return a, b


def Build_key_bridging(beta, r_in, r_dist, r_out, k_in, k_dist, uk_dist, uk_out, K_res):
    # beta: How many steps need to deduce a variable at most
    involved_k = {}
    involved_u = {}
    kb_state_k = {}
    kb_path_k = {}
    kb_state_u = {}
    kb_path_u = {}
    m = Nb * (r_in + r_dist + r_out + 1)
    for i in range(m):
        involved_k[i] = Rijndael.addVars(4, vtype=GRB.BINARY)
        involved_u[i] = Rijndael.addVars(4, vtype=GRB.BINARY)

    for rd in range(r_in + 1):
        for i in range(Nb):
            if k_in.get(rd) is None:
                Rijndael.addConstrs(involved_k[Nb * rd + i][j] == 0 for j in range(4))
            else:
                Rijndael.addConstrs(
                    involved_k[Nb * rd + i][j] == k_in[rd][4 * i + j] for j in range(4)
                )
            Rijndael.addConstrs(involved_u[Nb * rd + i][j] == 0 for j in range(4))

    for rd in range(1, r_dist):
        for i in range(Nb):
            idx = Nb * (rd + r_in) + i
            if k_dist.get(rd) is None:
                Rijndael.addConstrs(involved_k[idx][j] == 0 for j in range(4))
            else:
                Rijndael.addConstrs(
                    involved_k[idx][j] == k_dist[rd][4 * i + j] for j in range(4)
                )
            if uk_dist.get(rd) is None:
                Rijndael.addConstrs(involved_u[idx][j] == 0 for j in range(4))
            else:
                Rijndael.addConstrs(
                    involved_u[idx][j] == uk_dist[rd][4 * i + j] for j in range(4)
                )

    # Last round: MC omitted - k; MC exists - u
    for rd in range(r_out):
        for i in range(Nb):
            idx = Nb * (rd + r_in + r_dist + 1) + i
            # if rd == r_out - 1:
            #     Rijndael.addConstrs(
            #         involved_k[idx][j] == state_W_nl[rd][4 * i + j] for j in range(4)
            #     )
            #     Rijndael.addConstrs(involved_u[idx][j] == 0 for j in range(4))
            # else:
            Rijndael.addConstrs(involved_k[idx][j] == 0 for j in range(4))
            if uk_out.get(rd) is None:
                Rijndael.addConstrs(involved_u[idx][j] == 0 for j in range(4))
            else:
                Rijndael.addConstrs(
                    involved_u[idx][j] == uk_out[rd][4 * i + j] for j in range(4)
                )

    n = 4 * m
    kb_state_k[0] = Rijndael.addVars(n, vtype=GRB.BINARY)
    kb_state_u[0] = Rijndael.addVars(n, vtype=GRB.BINARY)
    for i in range(beta):  # For each step
        kb_state_k[i + 1] = Rijndael.addVars(n, vtype=GRB.BINARY)
        kb_path_k[i + 1] = {}
        for j in range(n):  # For each variable in k
            # Model of k
            kb_path_k[i + 1][j] = Rijndael.addVars(10, vtype=GRB.BINARY)
            col = j >> 2
            row = j % 4
            # For 10 relation
            for k in range(9):
                a, b = Get_idx_k(m, col, row, k)
                if a != -1 and b != -1:
                    Rijndael.addGenConstrAnd(
                        kb_path_k[i + 1][j][k], [kb_state_k[i][a], kb_state_k[i][b]]
                    )
                else:
                    Rijndael.addConstr(kb_path_k[i + 1][j][k] == 0)
            # Relation: k = MC(u)
            expr_sum = quicksum(
                kb_state_k[i][calc(col, r)] for r in range(4)
            ) + quicksum(kb_state_u[i][calc(col, r)] for r in range(4))
            Rijndael.addGenConstrIndicator(
                kb_path_k[i + 1][j][9], True, expr_sum, GRB.GREATER_EQUAL, 4
            )
            Rijndael.addGenConstrIndicator(
                kb_path_k[i + 1][j][9], False, expr_sum, GRB.LESS_EQUAL, 3
            )
            Rijndael.addConstr(kb_state_k[i + 1][j] >= kb_state_k[i][j])  # already 1
            Rijndael.addConstrs(
                kb_state_k[i + 1][j] >= kb_path_k[i + 1][j][k] for k in range(10)
            )
            Rijndael.addConstr(
                kb_state_k[i + 1][j]
                <= kb_state_k[i][j]
                + quicksum(kb_path_k[i + 1][j][k] for k in range(10))
            )

        # Model of u
        kb_state_u[i + 1] = Rijndael.addVars(n, vtype=GRB.BINARY)
        kb_path_u[i + 1] = {}
        for j in range(n):  # For each variable in u
            kb_path_u[i + 1][j] = Rijndael.addVars(7, vtype=GRB.BINARY)
            col = j >> 2
            row = j % 4
            # For 7 relation
            # len = 4 does not apply to u
            for k in range(6):
                a, b = Get_idx_u(m, col, row, k)
                if a != -1 and b != -1:
                    Rijndael.addGenConstrAnd(
                        kb_path_u[i + 1][j][k], [kb_state_u[i][a], kb_state_u[i][b]]
                    )
                else:
                    Rijndael.addConstr(kb_path_u[i + 1][j][k] == 0)
            # Relation: k = MC(u)
            expr_sum = quicksum(
                kb_state_k[i][calc(col, r)] for r in range(4)
            ) + quicksum(kb_state_u[i][calc(col, r)] for r in range(4))
            Rijndael.addGenConstrIndicator(
                kb_path_u[i + 1][j][6], True, expr_sum, GRB.GREATER_EQUAL, 4
            )
            Rijndael.addGenConstrIndicator(
                kb_path_u[i + 1][j][6], False, expr_sum, GRB.LESS_EQUAL, 3
            )
            Rijndael.addConstr(kb_state_u[i + 1][j] >= kb_state_u[i][j])  # already 1
            Rijndael.addConstrs(
                kb_state_u[i + 1][j] >= kb_path_u[i + 1][j][k] for k in range(7)
            )
            Rijndael.addConstr(
                kb_state_u[i + 1][j]
                <= kb_state_u[i][j] + quicksum(kb_path_u[i + 1][j][k] for k in range(7))
            )

    for i in range(m):
        Rijndael.addConstrs(
            kb_state_k[beta][4 * i + j] >= involved_k[i][j] for j in range(4)
        )
        Rijndael.addConstrs(
            kb_state_u[beta][4 * i + j] >= involved_u[i][j] for j in range(4)
        )
    K_res.add(
        quicksum(kb_state_k[0][i] for i in range(n))
        + quicksum(kb_state_u[0][i] for i in range(n))
    )


def Build_value_constraint(beta, r_in, r_dist, r_out):
    Rijndael.addConstr(Val_Con >= 0)

    # Limited by the table lookup method of distinguisher

    # Limited by the table lookup method of key recovery
    # VAR D1 - guessed part one
    state_D1_l[r_out] = Rijndael.addVars(
        4 * Nb, vtype=GRB.BINARY, name="state_D1_l_" + str(r_out)
    )
    Rijndael.addConstrs(state_D1_l[r_out][i] == 1 for i in range(4 * Nb))
    for rd in range(r_out - 1, -1, -1):
        state_D1_nl[rd] = Rijndael.addVars(
            4 * Nb, vtype=GRB.BINARY, name="state_D1_nl_" + str(rd)
        )
        state_D1_l[rd] = Rijndael.addVars(
            4 * Nb, vtype=GRB.BINARY, name="state_D1_l_" + str(rd)
        )
        # MC
        Build_mixcolumn_and(state_D1_nl[rd], state_D1_l[rd + 1])
        # SR
        key_D1[rd] = Rijndael.addVars(
            4 * Nb, vtype=GRB.BINARY, name="key_D1_" + str(rd)
        )
        for i in range(4 * Nb):
            Rijndael.addGenConstrAnd(
                state_D1_l[rd][SR(i)], [key_D1[rd][i], state_D1_nl[rd][i]]
            )
    # VAR D2 - guessed part two
    state_D2_l[r_out] = Rijndael.addVars(
        4 * Nb, vtype=GRB.BINARY, name="state_D2_l_" + str(r_out)
    )
    Rijndael.addConstrs(state_D2_l[r_out][i] == 1 for i in range(4 * Nb))
    for rd in range(r_out - 1, -1, -1):
        state_D2_nl[rd] = Rijndael.addVars(
            4 * Nb, vtype=GRB.BINARY, name="state_D2_nl_" + str(rd)
        )
        state_D2_l[rd] = Rijndael.addVars(
            4 * Nb, vtype=GRB.BINARY, name="state_D2_l_" + str(rd)
        )
        # MC
        Build_mixcolumn_and(state_D2_nl[rd], state_D2_l[rd + 1])
        # SR
        key_D2[rd] = Rijndael.addVars(
            4 * Nb, vtype=GRB.BINARY, name="key_D2_" + str(rd)
        )
        for i in range(4 * Nb):
            Rijndael.addGenConstrAnd(
                state_D2_l[rd][SR(i)], [key_D2[rd][i], state_D2_nl[rd][i]]
            )
    # Calculate Kg1, Kg2
    Build_key_bridging(beta, r_in, r_dist, r_out, state_M_l, {}, {}, key_D1, Kg1)
    Build_key_bridging(beta, r_in, r_dist, r_out, state_M_l, {}, {}, key_D2, Kg2)
    # VAR H - table part
    state_H_l[0] = Rijndael.addVars(4 * Nb, vtype=GRB.BINARY, name="state_H_l_0")
    for i in range(4 * Nb):
        Rijndael.addGenConstrIndicator(H_off, True, state_H_l[0][i] == 0)
        Rijndael.addGenConstrIndicator(
            H_off, True, state_D1_l[0][i] + state_D2_l[0][i] >= state_W_l[0][i]
        )
        Rijndael.addGenConstrIndicator(H_off, False, state_H_l[0][i] == state_W_l[0][i])
    for rd in range(r_out):
        state_H_nl[rd] = Rijndael.addVars(
            4 * Nb, vtype=GRB.BINARY, name="state_H_nl_" + str(rd)
        )
        state_H_l[rd + 1] = Rijndael.addVars(
            4 * Nb, vtype=GRB.BINARY, name="state_H_l_" + str(rd + 1)
        )
        key_H[rd] = Rijndael.addVars(4 * Nb, vtype=GRB.BINARY, name="key_H_" + str(rd))
        # SR Logic: H = not(D1 or D2) and W
        for i in range(4 * Nb):
            Rijndael.addConstr(state_H_nl[rd][i] <= 1 - state_D1_l[rd][SR(i)])
            Rijndael.addConstr(state_H_nl[rd][i] <= 1 - state_D2_l[rd][SR(i)])
            Rijndael.addConstr(state_H_nl[rd][i] <= state_H_l[rd][SR(i)])
            Rijndael.addConstr(
                state_H_nl[rd][i]
                >= state_H_l[rd][SR(i)] - state_D1_l[rd][SR(i)] - state_D2_l[rd][SR(i)]
            )
            Rijndael.addConstr(key_H[rd][i] == state_H_nl[rd][i])
            Kt.add(key_H[rd][i])
        # MC
        tmp_nl = Rijndael.addVars(4 * Nb, vtype=GRB.BINARY)
        for i in range(4 * Nb):
            Rijndael.addConstr(tmp_nl[i] <= 1 - state_D1_nl[rd][i])
            Rijndael.addConstr(tmp_nl[i] <= 1 - state_D2_nl[rd][i])
            Rijndael.addConstr(tmp_nl[i] <= state_H_nl[rd][i])
            Rijndael.addConstr(
                tmp_nl[i] >= state_H_nl[rd][i] - state_D1_nl[rd][i] - state_D2_nl[rd][i]
            )
        Build_mixcolumn_or(state_H_l[rd + 1], tmp_nl)

    # Calculate St, St1
    for i in range(4 * Nb):
        tmp_a_init = Rijndael.addVar(vtype=GRB.BINARY)
        Rijndael.addGenConstrOr(tmp_a_init, [state_D1_l[0][i], state_D2_l[0][i]])
        res_a_init = Rijndael.addVar(vtype=GRB.BINARY)
        Rijndael.addGenConstrAnd(res_a_init, [state_H_l[0][i], tmp_a_init])
        St.add(res_a_init)

        res_b_init = Rijndael.addVar(vtype=GRB.BINARY)
        Rijndael.addGenConstrAnd(res_b_init, [state_H_l[0][i], state_D1_l[0][i]])
        St1.add(res_b_init)
    for rd in range(r_out):
        for i in range(4 * Nb):
            tmp_a = Rijndael.addVar(vtype=GRB.BINARY)
            Rijndael.addGenConstrOr(
                tmp_a, [state_D1_l[rd + 1][i], state_D2_l[rd + 1][i]]
            )
            res_a = Rijndael.addVar(vtype=GRB.BINARY)
            Rijndael.addGenConstrAnd(res_a, [state_H_l[rd + 1][i], tmp_a])
            St.add(res_a)
            tmp_b = Rijndael.addVar(vtype=GRB.BINARY)
            Rijndael.addGenConstrOr(tmp_b, [state_D1_nl[rd][i], state_D2_nl[rd][i]])
            res_b = Rijndael.addVar(vtype=GRB.BINARY)
            Rijndael.addGenConstrAnd(res_b, [state_H_nl[rd][i], tmp_b])
            St.add(res_b)

    for rd in range(r_out):
        for i in range(4 * Nb):
            res_a = Rijndael.addVar(vtype=GRB.BINARY)
            Rijndael.addGenConstrAnd(
                res_a, [state_H_l[rd + 1][i], state_D1_l[rd + 1][i]]
            )
            St1.add(res_a)
            res_b = Rijndael.addVar(vtype=GRB.BINARY)
            Rijndael.addGenConstrAnd(res_b, [state_H_nl[rd][i], state_D1_nl[rd][i]])
            St1.add(res_b)

    # Memory Limit
    Rijndael.addConstr((Val_Con + 1) * St + Kt <= Deg)


def Build_key_dependent_sieve(beta, r_in, r_dist, r_out):
    Build_key_bridging(beta, r_in, r_dist, r_out, {}, key_dist, ukey_dist, {}, K_off)
    Build_key_bridging(
        beta, r_in, r_dist, r_out, state_M_l, key_dist, ukey_dist, state_W_nl, K_cup
    )
    Key_sieve.add(-K_off)


# Objective function
def Set_objective():
    # Obj_Offline
    Rijndael.addConstr(Obj_Offline == 8 * (Deg - Val_Con - Key_sieve))
    # Obj_Online
    Rijndael.addConstr(Obj_Online >= 8 * ((Val_Con + 1) * St + Kt - Val_Con) - 6)
    Rijndael.addConstr(Obj_Online >= 8 * (Val_Con + Kg1))
    Rijndael.addConstr(Obj_Online >= 8 * (Kg1 + (Val_Con + 1) * (St - St1) + Kt))
    Rijndael.addConstr(Obj_Online >= 8 * (Val_Con + Kg2))
    Rijndael.addConstr(Obj_Online >= 8 * Key)
    # Obj_Data: Consider DTM trade-off
    Rijndael.addConstr(Obj_Data >= 8 * Plain_diff)
    Rijndael.addConstr(2 * Obj_Data >= 16 * (1 + Val_Con) + (Obj_Offline - Obj_Online))
    # Obj_Memory
    Rijndael.addConstr(Obj_Memory == Obj_Offline - 8 * (K_off + Key - K_cup))
    # Obj
    Rijndael.addConstr(Obj >= Obj_Online)
    Rijndael.addConstr(2 * Obj >= Obj_Offline + Obj_Online)

    # Objective
    Rijndael.setObjectiveN(Obj, index=0, priority=4, name="Obj_Complexity", weight=1.0)
    Rijndael.setObjectiveN(Obj_Data, index=1, priority=3, name="Obj_Data", weight=1.0)
    Rijndael.setObjectiveN(Key, index=2, priority=2, name="Key", weight=1.0)
    Rijndael.setObjectiveN(
        Obj_Memory, index=3, priority=1, name="Obj_Memory", weight=1.0
    )
    # Rijndael.setObjectiveN(St, index=3, priority=4, name="Obj_St", weight=1.0)
    # Rijndael.setObjectiveN(Kt, index=4, priority=3, name="Obj_Kt", weight=1.0)
    # Rijndael.setObjectiveN(Kg1, index=5, priority=2, name="Obj_Kg1", weight=1.0)
    # Rijndael.setObjectiveN(Kg2, index=6, priority=1, name="Obj_Kg2", weight=1.0)
    Rijndael.ModelSense = GRB.MINIMIZE


def print_block(vals):
    for i in range(4):
        for j in range(Nb):
            print(vals[4 * j + i], end=" ")
        print("")


def print_flow(name, keys, l_dict, nl_dict, offset=0):
    print(f"---------- Var {name} ----------")
    keys = list(keys)
    for i, rd in enumerate(keys):
        disp_idx = rd + offset
        if l_dict and rd in l_dict:
            print(f"{name}_linear[{disp_idx}]")
            vals = [round(l_dict[rd][k].Xn) for k in range(4 * Nb)]
            print_block(vals)
            if nl_dict and rd in nl_dict:
                print("-- (ARK, SB, SR) ->")

        if nl_dict and rd in nl_dict:
            print(f"{name}_non_linear[{disp_idx}]")
            vals = [round(nl_dict[rd][k].Xn) for k in range(4 * Nb)]
            print_block(vals)
            if i < len(keys) - 1:
                print("-- (MC) ->")


def Print_key_out(r_in, r_dist, r_out, key_str, kv):
    print(f"---------- Key {key_str} ----------")
    offset = r_dist + r_in
    for rd in range(r_out):
        print(f"Key {key_str}[{rd + offset}]")
        vals = [round(kv[rd][k].Xn) for k in range(4 * Nb)]
        print_block(vals)


def Start_solver(r_in, r_dist, r_out, Print_result):
    # AES.setParam("OutputFlag", 0)
    Rijndael.Params.PoolSearchMode = 0
    Rijndael.Params.PoolSolutions = 1
    # Rijndael.Params.PoolGap = 2.0
    # Rijndael.Params.Threads = 128
    # Rijndael.Params.TimeLimit = 100
    # Rijndael.setParam("IntFeasTol", 1e-9)
    Rijndael.optimize()
    if Print_result == False:
        return

    DEFAULT_OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    file = DEFAULT_OUTPUT_PATH.open("w", encoding="utf-8")
    sys.stdout = file
    print("Model Status:", Rijndael.Status)
    if Rijndael.Status == 2 or Rijndael.Status == 9 or Rijndael.Status == 11:
        print("Min_Obj: %g" % Rijndael.ObjVal)

        # All solutions
        for k in range(Rijndael.SolCount):
            Rijndael.Params.SolutionNumber = k
            print(
                "******** Sol no.{}    Obj = {}    Offline = {}    Online = {}    Data = {}    Memory = {} ********".format(
                    k + 1,
                    round(Obj.Xn),
                    round(Obj_Offline.Xn),
                    round(Obj_Online.Xn),
                    round(Obj_Data.Xn),
                    round(Obj_Memory.Xn),
                )
            )
            print(
                "******** Deg = {}    Key = {}    Plaintext Difference = {} ********".format(
                    round(Evaluate_expr(Deg)),
                    round(Evaluate_expr(Key)),
                    round(Evaluate_expr(Plain_diff)),
                )
            )
            print(
                "******** Val_Con = {}    Kg1 = {}    Kg2 = {}    St = {}    Kt = {}    St1 = {} ********".format(
                    round(Val_Con.Xn),
                    round(Evaluate_expr(Kg1)),
                    round(Evaluate_expr(Kg2)),
                    round(Evaluate_expr(St)),
                    round(Evaluate_expr(Kt)),
                    round(Evaluate_expr(St1)),
                )
            )
            print(
                "******** K_off = {}    K_cup = {}    Key_sieve = {} ********".format(
                    round(Evaluate_expr(K_off)),
                    round(Evaluate_expr(K_cup)),
                    round(Evaluate_expr(Key_sieve)),
                )
            )
            strr = "A = [ "
            for i in range(4 * Nb):
                if round(state_Z[1][i].Xn) == 1:
                    strr += str(i) + " "
            strr += "]    B = [ "
            for i in range(4 * Nb):
                if round(state_Z[r_dist - 1][i].Xn) == 1:
                    strr += str(i) + " "
            strr += "]"
            print(strr)

            print_flow("Z", range(1, r_dist), state_Z, None, offset=r_in)
            print_flow("X", range(r_dist), state_X_l, state_X_nl, offset=r_in)
            print_flow("Y", range(1, r_dist + 1), state_Y_l, state_Y_nl, offset=r_in)

            print_flow("M", range(r_in + 1), state_M_l, state_M_nl)
            print_flow(
                "W", range(r_out + 1), state_W_l, state_W_nl, offset=r_in + r_dist
            )

            print_flow(
                "D1", range(r_out + 1), state_D1_l, state_D1_nl, offset=r_in + r_dist
            )
            print_flow(
                "D2", range(r_out + 1), state_D2_l, state_D2_nl, offset=r_in + r_dist
            )
            print_flow(
                "H", range(r_out + 1), state_H_l, state_H_nl, offset=r_in + r_dist
            )
            Print_key_out(r_in, r_dist, r_out, "D1", key_D1)
            Print_key_out(r_in, r_dist, r_out, "D2", key_D2)
            Print_key_out(r_in, r_dist, r_out, "H", key_H)

            print_flow("k_dist", range(2, r_dist), key_dist, None, offset=r_in)
            print_flow("u_dist", range(2, r_dist), ukey_dist, None, offset=r_in)

    file.close()
    sys.stdout = sys.__stdout__


def Search_attack(nb, nk, r_dist, r_in, r_out, extra_constr, beta):
    Reset_model(nb, nk)
    print("=== Parameters ===")
    print("r_dist:", r_dist)
    print("r_in:", r_in)
    print("r_out:", r_out)
    print("extra_constr:", extra_constr)
    print("beta:", beta)
    print("==================")
    Build_distinguisher(r_dist)
    Build_key_recovery(beta, r_in, r_dist, r_out)
    # Rijndael.addConstr(quicksum(state_Z[r_dist - 1][i] for i in range(4 * Nb)) == 3)
    # Rijndael.addConstr(Obj >= 208)
    # Rijndael.addConstr(Obj_Data >= 64)
    # Rijndael.addConstr(Key >= 26)

    if len(extra_constr) > 0:
        st0 = extra_constr[0]
        st1 = extra_constr[1]
        ed0 = extra_constr[2]
        ed1 = extra_constr[3]
        if len(st0) > 0:
            Rijndael.addConstrs(state_X_nl[0][i] == st0[i] for i in range(4 * Nb))
        if len(st1) > 0:
            Rijndael.addConstrs(state_X_l[1][i] == st1[i] for i in range(4 * Nb))
        if len(ed0) > 0:
            Rijndael.addConstrs(
                state_Y_nl[r_dist - 1][i] == ed0[i] for i in range(4 * Nb)
            )
        if len(ed1) > 0:
            Rijndael.addConstrs(state_Y_l[r_dist][i] == ed1[i] for i in range(4 * Nb))
    Build_value_constraint(beta, r_in, r_dist, r_out)
    Build_key_dependent_sieve(beta, r_in, r_dist, r_out)
    Set_objective()
    Start_solver(r_in, r_dist, r_out, True)
    if Rijndael.Status == 2 or Rijndael.Status == 9 or Rijndael.Status == 11:
        R = r_in + r_dist + r_out + 1
        Guessed_k = [[0] * (4 * Nb) for _ in range(R)]
        Guessed_u = [[0] * (4 * Nb) for _ in range(R)]
        for rd in range(r_in + 1):
            for i in range(4 * Nb):
                if round(state_M_l[rd][i].Xn) == 1:
                    Guessed_k[rd][i] = 1
        for rd in range(r_out):
            for i in range(4 * Nb):
                if round(state_W_nl[rd][i].Xn) == 1:
                    Guessed_u[rd + r_in + 1 + r_dist][i] = 1
        return (
            round(Evaluate_expr(Deg)),
            round(Evaluate_expr(Plain_diff)),
            Guessed_k,
            Guessed_u,
        )
    return -1, -1, [], []


if __name__ == "__main__":
    # fmt: off
    # st0 = [0, 0, 0, 0,
    #        0, 0, 0, 0,
    #        0, 0, 0, 0,
    #        0, 0, 0, 0,
    #        0, 1, 0, 1,
    #        0, 0, 0, 0,
    #        0, 0, 0, 0,
    #        0, 0, 0, 0]
    # st1 = [0, 0, 0, 0,
    #        0, 0, 0, 0,
    #        0, 0, 0, 0,
    #        0, 0, 0, 0,
    #        1, 0, 1, 1,
    #        0, 0, 0, 0,
    #        0, 0, 0, 0,
    #        0, 0, 0, 0]
    # ed0 = [0, 0, 0, 0,
    #        0, 0, 0, 0,
    #        0, 0, 0, 0,
    #        0, 0, 0, 0,
    #        0, 0, 0, 0,
    #        0, 0, 0, 0,
    #        1, 0, 0, 0,
    #        0, 0, 0, 0]
    # ed1 = [0, 0, 0, 0,
    #        0, 0, 0, 0,
    #        0, 0, 0, 0,
    #        0, 0, 0, 0,
    #        0, 0, 0, 0,
    #        0, 0, 0, 0,
    #        1, 1, 1, 1,
    #        0, 0, 0, 0]
    # fmt: on

    # nb,nk,r_dist,r_in,r_out,extra_constr,beta
    # deg, data, Gk, Gu = Search_attack(8, 8, 6, 1, 2, [st0, st1, ed0, ed1], 1)

    # fmt: off
    st0 = [0, 0, 0, 0,
           0, 0, 0, 1,
           0, 0, 0, 0,
           0, 0, 0, 0,
           0, 0, 0, 0]
    st1 = [0, 0, 0, 0,
           1, 1, 1, 1,
           0, 0, 0, 0,
           0, 0, 0, 0,
           0, 0, 0, 0]
    ed0 = [0, 0, 0, 0,
           0, 0, 0, 0,
           0, 0, 0, 0,
           0, 0, 0, 0,
           0, 0, 0, 1]
    ed1 = [0, 0, 0, 0,
           0, 0, 0, 0,
           0, 0, 0, 0,
           0, 0, 0, 0,
           1, 1, 1, 1]
    # fmt: on
    deg, data, Gk, Gu = Search_attack(5, 8, 6, 1, 2, [st0, st1, ed0, ed1], 1)

    print(Gk)
    print(Gu)
