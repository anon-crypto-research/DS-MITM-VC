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
    global state_X_l, state_X_nl, state_Y_l, state_Y_nl
    global state_Z, state_Z_l, state_Z_nl
    state_X_l = {}  # VAR X: parameters
    state_X_nl = {}
    state_Y_l = {}  # VAR Y: parameters
    state_Y_nl = {}
    state_Z = {}
    state_Z_l = {}  # VAR Z: linear state
    state_Z_nl = {}  # VAR Z: non-linear state

    global state_W_l, state_W_nl, state_M_l, state_M_nl
    state_W_l = {}
    state_W_nl = {}  # [r_in + r_dist + 1, r_in + r_dist + r_out + 1)
    state_M_l = {}  # [0, r_in + 1)
    state_M_nl = {}

    global key_dist, ukey_dist
    key_dist = {}  # [r_in + 1, r_in + r_dist)
    ukey_dist = {}  # [r_in + 1, r_in + r_dist)

    global Deg, Plain_diff, Key, Key_sieve, K_off, K_cup
    Deg = LinExpr()
    Plain_diff = LinExpr()
    Key = LinExpr()
    Key_sieve = LinExpr()
    K_off = LinExpr()
    K_cup = LinExpr()

    global Num_structures, Pr_td
    Num_structures = Rijndael.addVar(vtype=GRB.INTEGER, name="Num_structures")
    Pr_td = Rijndael.addVar(vtype=GRB.INTEGER, name="Pr_td")

    global Obj_Data, Obj_Memory, Obj_Offline, Obj_Online, Obj
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

    # Nontrivial
    Rijndael.addConstr(quicksum(state_X_nl[0][i] for i in range(4 * Nb)) >= 1)
    Rijndael.addConstr(quicksum(state_Y_l[r_dist][i] for i in range(4 * Nb)) >= 1)


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


def Build_key_dependent_sieve(beta, r_in, r_dist, r_out):
    Build_key_bridging(beta, r_in, r_dist, r_out, {}, key_dist, ukey_dist, {}, K_off)
    Build_key_bridging(
        beta, r_in, r_dist, r_out, state_M_l, key_dist, ukey_dist, state_W_nl, K_cup
    )
    Key_sieve.add(-K_off)


def _add_and_with_default_one(target, operands):
    operands = [operand for operand in operands if operand is not None]
    if len(operands) == 0:
        Rijndael.addConstr(target == 1)
    elif len(operands) == 1:
        Rijndael.addConstr(target == operands[0])
    else:
        Rijndael.addGenConstrAnd(target, operands)


def _add_nonzero_column_zero_count(state, rd, j, expr, name):
    col_state = [state[rd][4 * j + i] for i in range(4)]
    col_active = Rijndael.addVar(vtype=GRB.BINARY, name=f"{name}_{rd}_{j}")
    Rijndael.addGenConstrOr(col_active, col_state)
    expr.add(4 * col_active)
    expr.add(-quicksum(col_state))


def Build_differential_enumeration(r_in, r_dist, r_out):
    global state_Z

    for rd in range(1, r_dist + 1):
        state_Z_l[rd] = Rijndael.addVars(
            4 * Nb, vtype=GRB.BINARY, name="state_Z_l_" + str(rd)
        )
        for i in range(4 * Nb):
            _add_and_with_default_one(
                state_Z_l[rd][i],
                [
                    state_X_l.get(rd, {}).get(i),
                    state_Y_l.get(rd, {}).get(i),
                ],
            )

    for rd in range(r_dist):
        state_Z_nl[rd] = Rijndael.addVars(
            4 * Nb, vtype=GRB.BINARY, name="state_Z_nl_" + str(rd)
        )
        for i in range(4 * Nb):
            _add_and_with_default_one(
                state_Z_nl[rd][i],
                [
                    state_X_nl.get(rd, {}).get(i),
                    state_Y_nl.get(rd, {}).get(i),
                ],
            )

    state_Z = state_Z_l

    for rd in range(r_dist):
        for j in range(Nb):
            col_state = [state_Z_nl[rd][4 * j + i] for i in range(4)] + [
                state_Z_l[rd + 1][4 * j + i] for i in range(4)
            ]
            col_active = Rijndael.addVar(
                vtype=GRB.BINARY, name=f"deg_col_active_{rd}_{j}"
            )
            Rijndael.addGenConstrOr(col_active, col_state)
            Deg.add(quicksum(col_state))
            Deg.add(-4 * col_active)

    pr_td_expr = LinExpr()
    for rd in range(1, r_in + 1):
        for j in range(Nb):
            _add_nonzero_column_zero_count(
                state_M_l, rd, j, pr_td_expr, "pr_td_m_col_active"
            )
    for rd in range(1, r_dist + 1):
        for j in range(Nb):
            _add_nonzero_column_zero_count(
                state_Z_l, rd, j, pr_td_expr, "pr_td_z_col_active"
            )
    Rijndael.addConstr(Pr_td == pr_td_expr)
    Rijndael.addConstr(Num_structures + 16 * Plain_diff - 1 == 8 * Pr_td)

    # Extract involved key bytes in the distinguisher part
    Extract_involved_key(r_dist)


# Objective function
def Set_objective():
    # Obj_Offline
    Rijndael.addConstr(Obj_Offline == 8 * (Deg - Key_sieve))
    # Obj_Online
    Rijndael.addConstr(Obj_Online >= 8 * Key)
    Rijndael.addConstr(Obj_Online >= Obj_Data)
    # Obj_Data: Consider DTM trade-off
    Rijndael.addConstr(Obj_Data >= 8 * Plain_diff)
    Rijndael.addConstr(Obj_Data >= Num_structures + 8 * Plain_diff)
    Rijndael.addConstr(2 * Obj_Data >= 16 + (Obj_Offline - Obj_Online))
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
    Rijndael.getMultiobjEnv(0).setParam("TimeLimit", 600)
    Rijndael.getMultiobjEnv(1).setParam("TimeLimit", 180)
    Rijndael.getMultiobjEnv(2).setParam("TimeLimit", 240)
    Rijndael.getMultiobjEnv(3).setParam("TimeLimit", 240)
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
    if Rijndael.Status in [2, 9, 11, 13]:
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
                "******** K_off = {}    K_cup = {}    Key_sieve = {} ********".format(
                    round(Evaluate_expr(K_off)),
                    round(Evaluate_expr(K_cup)),
                    round(Evaluate_expr(Key_sieve)),
                )
            )
            print(
                "******** Num_structures = {}    Pr_td = {} ********".format(
                    round(Num_structures.Xn),
                    round(Pr_td.Xn),
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

            print_flow("Z", range(r_dist + 1), state_Z_l, state_Z_nl, offset=r_in)
            print_flow("X", range(r_dist), state_X_l, state_X_nl, offset=r_in)
            print_flow("Y", range(1, r_dist + 1), state_Y_l, state_Y_nl, offset=r_in)

            print_flow("M", range(r_in + 1), state_M_l, state_M_nl)
            print_flow(
                "W", range(r_out + 1), state_W_l, state_W_nl, offset=r_in + r_dist
            )

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
    Build_differential_enumeration(r_in, r_dist, r_out)
    Build_key_dependent_sieve(beta, r_in, r_dist, r_out)
    Set_objective()
    Start_solver(r_in, r_dist, r_out, True)
    if Rijndael.Status in [2, 9, 11, 13]:
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
    st0 = [0, 0, 0, 0,
           0, 0, 0, 0,
           0, 0, 0, 0,
           1, 1, 0, 0,
           0, 0, 0, 0,
           0, 0, 0, 0,
           0, 0, 0, 0,
           0, 0, 0, 0]
    st1 = [0, 0, 0, 0,
           0, 0, 0, 0,
           0, 0, 0, 0,
           1, 1, 0, 1,
           0, 0, 0, 0,
           0, 0, 0, 0,
           0, 0, 0, 0,
           0, 0, 0, 0]
    ed0 = [0, 0, 0, 0,
           1, 1, 1, 1,
           0, 0, 0, 0,
           0, 0, 0, 0,
           0, 0, 0, 0,
           0, 0, 0, 0,
           0, 0, 0, 0,
           0, 0, 0, 0]
    ed1 = [0, 0, 0, 0,
           0, 0, 1, 0,
           0, 0, 0, 0,
           0, 0, 0, 0,
           0, 0, 0, 0,
           0, 0, 0, 0,
           0, 0, 0, 0,
           0, 0, 0, 0]
    # fmt: on

    # nb,nk,r_dist,r_in,r_out,extra_constr,beta
    deg, data, Gk, Gu = Search_attack(8, 8, 6, 1, 3, [st0, st1, ed0, ed1], 1)
    # deg, data, Gk, Gu = Search_attack(8, 8, 6, 1, 3, [], 1)
    # print(deg, data)
    print(Gk)
    print(Gu)
