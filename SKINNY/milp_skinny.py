from typing import Any, Dict

from gurobipy import GRB, LinExpr, Model, QuadExpr, Var, quicksum
from pathlib import Path
import sys

MODULE_DIR = Path(__file__).resolve().parent
DEFAULT_OUTPUT_PATH = MODULE_DIR / "results" / "output"

VarMap = Dict[int, Any]

SKINNY: Model

state_X_l: VarMap = {}
state_X_nl: VarMap = {}
state_Y_l: VarMap = {}
state_Y_nl: VarMap = {}
state_Z: VarMap = {}

state_V_sb: VarMap = {}
state_V: VarMap = {}
key_V: VarMap = {}

con_1: VarMap = {}
con_2: VarMap = {}
con_3: VarMap = {}
con_4: VarMap = {}

state_W_l: VarMap = {}
state_W_nl: VarMap = {}
state_M_l: VarMap = {}
state_M_nl: VarMap = {}
state_O_l: VarMap = {}
state_O_nl: VarMap = {}

state_D_l: VarMap = {}
state_D_nl: VarMap = {}
key_D: VarMap = {}
state_H_l: VarMap = {}
state_H_nl: VarMap = {}
key_H: VarMap = {}

key_bri: Dict[int, LinExpr] = {}
key_dep: Dict[int, LinExpr] = {}
key_bri_g: Dict[int, LinExpr] = {}

Deg: LinExpr
Con: LinExpr
Key_sieve: LinExpr
Key: LinExpr
Kg: LinExpr
St: LinExpr
Kt: LinExpr
Plain_diff: LinExpr
Start: LinExpr
Val_Con: Var
Obj_Data: Var
Obj_Offline: Var
Obj_Online: Var
Obj: Var


SR = [0, 1, 2, 3, 7, 4, 5, 6, 10, 11, 8, 9, 13, 14, 15, 12]
PT = [9, 15, 8, 13, 10, 14, 12, 11, 0, 1, 2, 3, 4, 5, 6, 7]


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


def Reset_model(r_in, r_dist, r_out):
    # l -(SB,AK,SR)-> nl -(MC)-> l
    global SKINNY
    SKINNY = Model("SKINNY")

    global state_X_l, state_X_nl, state_Y_l, state_Y_nl, state_Z
    state_X_l = {}  # VAR X
    state_X_nl = {}
    state_Y_l = {}  # VAR Y
    state_Y_nl = {}
    state_Z = {}  # VAR Z

    global state_V_sb, state_V, key_V
    state_V_sb = {}  # VAR V
    state_V = {}
    key_V = {}  # key for VAR V

    global con_1, con_2, con_3, con_4
    con_1 = {}  # cipher-specific constraints (non-full key addition)
    con_2 = {}
    con_3 = {}
    con_4 = {}

    global state_W_l, state_W_nl, state_M_l, state_M_nl, state_O_l, state_O_nl
    state_W_l = {}  # VAR W
    state_W_nl = {}
    state_M_l = {}  # VAR M
    state_M_nl = {}
    state_O_l = {}  # VAR O
    state_O_nl = {}

    global state_D_l, state_D_nl, key_D, state_H_l, state_H_nl, key_H
    state_D_l = {}  # VAR D
    state_D_nl = {}
    key_D = {}
    state_H_l = {}  # VAR H
    state_H_nl = {}
    key_H = {}

    global key_bri, key_dep, key_bri_g
    key_bri = {}
    key_dep = {}
    key_bri_g = {}

    global Deg, Con, Key_sieve, Key, Kg, St, Kt, Plain_diff, Start
    global Val_Con, Obj_Data, Obj_Offline, Obj_Online, Obj
    Deg = LinExpr()
    Con = LinExpr()
    Key_sieve = LinExpr()
    Key = LinExpr()
    Kg = LinExpr()
    St = LinExpr()
    Kt = LinExpr()
    Plain_diff = LinExpr()
    Start = LinExpr()
    Val_Con = SKINNY.addVar(vtype=GRB.INTEGER, name="Val_Con")
    Obj_Data = SKINNY.addVar(vtype=GRB.INTEGER, name="Obj_Data")
    Obj_Offline = SKINNY.addVar(vtype=GRB.INTEGER, name="Obj_Offline")
    Obj_Online = SKINNY.addVar(vtype=GRB.INTEGER, name="Obj_Online")
    Obj = SKINNY.addVar(vtype=GRB.INTEGER, name="Obj")


def Build_shiftrow(a, ra, b, rb):
    # a - SR -> b
    SKINNY.addConstrs(a[ra][SR[i]] == b[rb][i] for i in range(16))


def Build_mixcolumn(a, ra, b, rb, typ):
    assert (
        typ == "forward_differential"
        or typ == "backward_differential"
        or typ == "forward_determination"
        or typ == "backward_determination"
    )
    # a - MC -> b
    if typ == "forward_differential":
        # MC - COL1
        SKINNY.addConstrs(b[rb][i] >= a[ra][i] for i in range(4))
        SKINNY.addConstrs(b[rb][i] >= a[ra][i + 8] for i in range(4))
        SKINNY.addConstrs(b[rb][i] >= a[ra][i + 12] for i in range(4))
        SKINNY.addConstrs(
            b[rb][i] <= a[ra][i] + a[ra][i + 8] + a[ra][i + 12] for i in range(4)
        )
        # MC - COL2
        SKINNY.addConstrs(b[rb][i + 4] == a[ra][i] for i in range(4))
        # MC - COL3
        SKINNY.addConstrs(b[rb][i + 8] >= a[ra][i + 4] for i in range(4))
        SKINNY.addConstrs(b[rb][i + 8] >= a[ra][i + 8] for i in range(4))
        SKINNY.addConstrs(b[rb][i + 8] <= a[ra][i + 4] + a[ra][i + 8] for i in range(4))
        # MC - COL4
        SKINNY.addConstrs(b[rb][i + 12] >= a[ra][i] for i in range(4))
        SKINNY.addConstrs(b[rb][i + 12] >= a[ra][i + 8] for i in range(4))
        SKINNY.addConstrs(b[rb][i + 12] <= a[ra][i] + a[ra][i + 8] for i in range(4))
    if typ == "backward_determination":
        # MC - COL1
        SKINNY.addConstrs(a[ra][i] >= b[rb][i] for i in range(4))
        SKINNY.addConstrs(a[ra][i] >= b[rb][i + 4] for i in range(4))
        SKINNY.addConstrs(a[ra][i] >= b[rb][i + 12] for i in range(4))
        SKINNY.addConstrs(
            a[ra][i] <= b[rb][i] + b[rb][i + 4] + b[rb][i + 12] for i in range(4)
        )
        # MC - COL2
        SKINNY.addConstrs(a[ra][i + 4] == b[rb][i + 8] for i in range(4))
        # MC - COL3
        SKINNY.addConstrs(a[ra][i + 8] >= b[rb][i] for i in range(4))
        SKINNY.addConstrs(a[ra][i + 8] >= b[rb][i + 8] for i in range(4))
        SKINNY.addConstrs(a[ra][i + 8] >= b[rb][i + 12] for i in range(4))
        SKINNY.addConstrs(
            a[ra][i + 8] <= b[rb][i] + b[rb][i + 8] + b[rb][i + 12] for i in range(4)
        )
        # MC - COL4
        SKINNY.addConstrs(a[ra][i + 12] == b[rb][i] for i in range(4))
    if typ == "forward_determination":
        # MC - COL1
        SKINNY.addConstrs(b[rb][i] == a[ra][i + 12] for i in range(4))
        # MC - COL2
        SKINNY.addConstrs(b[rb][i + 4] >= a[ra][i] for i in range(4))
        SKINNY.addConstrs(b[rb][i + 4] >= a[ra][i + 4] for i in range(4))
        SKINNY.addConstrs(b[rb][i + 4] >= a[ra][i + 8] for i in range(4))
        SKINNY.addConstrs(
            b[rb][i + 4] <= a[ra][i] + a[ra][i + 4] + a[ra][i + 8] for i in range(4)
        )
        # MC - COL3
        SKINNY.addConstrs(b[rb][i + 8] == a[ra][i + 4] for i in range(4))
        # MC - COL4
        SKINNY.addConstrs(b[rb][i + 12] >= a[ra][i + 4] for i in range(4))
        SKINNY.addConstrs(b[rb][i + 12] >= a[ra][i + 8] for i in range(4))
        SKINNY.addConstrs(b[rb][i + 12] >= a[ra][i + 12] for i in range(4))
        SKINNY.addConstrs(
            b[rb][i + 12] <= a[ra][i + 4] + a[ra][i + 8] + a[ra][i + 12]
            for i in range(4)
        )
    if typ == "backward_differential":
        # MC - COL1
        SKINNY.addConstrs(a[ra][i] == b[rb][i + 4] for i in range(4))
        # MC - COL2
        SKINNY.addConstrs(a[ra][i + 4] >= b[rb][i + 4] for i in range(4))
        SKINNY.addConstrs(a[ra][i + 4] >= b[rb][i + 8] for i in range(4))
        SKINNY.addConstrs(a[ra][i + 4] >= b[rb][i + 12] for i in range(4))
        SKINNY.addConstrs(
            a[ra][i + 4] <= b[rb][i + 4] + b[rb][i + 8] + b[rb][i + 12]
            for i in range(4)
        )
        # MC - COL3
        SKINNY.addConstrs(a[ra][i + 8] >= b[rb][i + 4] for i in range(4))
        SKINNY.addConstrs(a[ra][i + 8] >= b[rb][i + 12] for i in range(4))
        SKINNY.addConstrs(
            a[ra][i + 8] <= b[rb][i + 4] + b[rb][i + 12] for i in range(4)
        )
        # MC - COL4
        SKINNY.addConstrs(a[ra][i + 12] >= b[rb][i] for i in range(4))
        SKINNY.addConstrs(a[ra][i + 12] >= b[rb][i + 12] for i in range(4))
        SKINNY.addConstrs(a[ra][i + 12] <= b[rb][i] + b[rb][i + 12] for i in range(4))


# Non-full key additions
def Build_non_full_key_addition(r_dist):
    for rd in range(r_dist):
        con_1[rd] = SKINNY.addVars(4, vtype=GRB.BINARY, name="con_1_" + str(rd))
        con_2[rd] = SKINNY.addVars(4, vtype=GRB.BINARY, name="con_2_" + str(rd))
        con_3[rd] = SKINNY.addVars(4, vtype=GRB.BINARY, name="con_3_" + str(rd))
        con_4[rd] = SKINNY.addVars(4, vtype=GRB.BINARY, name="con_4_" + str(rd))
        for i in range(4):
            SKINNY.addConstr(
                state_Z[rd + 1][i + 4]
                + state_Z[rd + 1][i + 12]
                + state_Z[rd][SR[i + 8]]
                >= 3 * con_1[rd][i]
            )
            SKINNY.addConstr(
                state_Z[rd + 1][i + 4]
                + state_Z[rd + 1][i + 12]
                + state_Z[rd][SR[i + 8]]
                <= 2 + con_1[rd][i]
            )
            SKINNY.addConstr(
                state_Z[rd + 1][i] + state_Z[rd + 1][i + 12] + state_Z[rd][SR[i + 12]]
                >= 3 * con_2[rd][i]
            )
            SKINNY.addConstr(
                state_Z[rd + 1][i] + state_Z[rd + 1][i + 12] + state_Z[rd][SR[i + 12]]
                <= 2 + con_2[rd][i]
            )
            SKINNY.addConstr(
                state_Z[rd + 1][i]
                + state_Z[rd + 1][i + 4]
                + state_Z[rd][SR[i + 8]]
                + state_Z[rd][SR[i + 12]]
                >= 4 * con_3[rd][i]
            )
            SKINNY.addConstr(
                state_Z[rd + 1][i]
                + state_Z[rd + 1][i + 4]
                + state_Z[rd][SR[i + 8]]
                + state_Z[rd][SR[i + 12]]
                <= 3 + con_3[rd][i]
            )
            SKINNY.addGenConstrIndicator(
                con_4[rd][i],
                True,
                con_1[rd][i] + con_2[rd][i] + con_3[rd][i],
                GRB.GREATER_EQUAL,
                3,
            )
            SKINNY.addGenConstrIndicator(
                con_4[rd][i],
                False,
                con_1[rd][i] + con_2[rd][i] + con_3[rd][i],
                GRB.LESS_EQUAL,
                1,
            )
            if rd < r_dist - 1:
                Con.add(con_1[rd][i])
                Con.add(con_2[rd][i])
                Con.add(con_3[rd][i])
                Con.add(-con_4[rd][i])


# Key-dependent-sieve technique
def Build_key_dependent_sieve(attack_type, r_in, r_dist):
    # VAR V
    for rd in range(r_dist - 1):
        state_V_sb[rd] = SKINNY.addVars(
            16, vtype=GRB.BINARY, name="state_V_sb_" + str(rd)
        )
        state_V[rd + 1] = SKINNY.addVars(
            16, vtype=GRB.BINARY, name="state_V_" + str(rd + 1)
        )
        dummy_V = SKINNY.addVars(28, vtype=GRB.BINARY, name="dummy_V_" + str(rd + 1))
        key_V[rd] = SKINNY.addVars(8, vtype=GRB.BINARY, name="key_V_" + str(rd))
        for i in range(4):
            # Row 1 for V_sb
            SKINNY.addConstr(state_V_sb[rd][i] == state_V[rd + 1][i + 4])
            # Row 2 for V_sb
            SKINNY.addGenConstrAnd(
                dummy_V[i + 20],
                [state_V[rd + 1][i + 4], state_V[rd + 1][i + 12]],
            )
            SKINNY.addGenConstrOr(
                dummy_V[i + 24],
                [state_V_sb[rd][i + 8], dummy_V[i + 20]],
            )
            SKINNY.addGenConstrAnd(
                state_V_sb[rd][i + 4],
                [state_V[rd + 1][i + 8], dummy_V[i + 24]],
            )
            # Row 3 for V_sb
            SKINNY.addConstr(state_V_sb[rd][i + 8] == state_Z[rd][SR[i + 8]])
            # Row 4 for V_sb
            SKINNY.addConstr(state_V_sb[rd][i + 12] == state_Z[rd][SR[i + 12]])
            # Row 1 for V
            SKINNY.addGenConstrAnd(
                dummy_V[i], [state_V_sb[rd][i + 12], state_V[rd + 1][i + 12]]
            )
            SKINNY.addGenConstrOr(state_V[rd + 1][i], [dummy_V[i], state_Z[rd + 1][i]])
            # Row 2 for V
            SKINNY.addGenConstrAnd(
                dummy_V[i + 4], [state_V_sb[rd][i + 8], state_V[rd + 1][i + 12]]
            )
            SKINNY.addGenConstrOr(
                state_V[rd + 1][i + 4], [dummy_V[i + 4], state_Z[rd + 1][i + 4]]
            )
            # Row 3 for V
            SKINNY.addConstr(state_V[rd + 1][i + 8] == state_Z[rd + 1][i + 8])
            # Row 4 for V
            SKINNY.addGenConstrAnd(
                dummy_V[i + 8], [state_V_sb[rd][i + 12], state_V[rd + 1][i]]
            )
            SKINNY.addGenConstrAnd(
                dummy_V[i + 12], [state_V_sb[rd][i + 8], state_V[rd + 1][i + 4]]
            )
            SKINNY.addGenConstrOr(
                dummy_V[i + 16],
                [dummy_V[i + 8], dummy_V[i + 12]],
            )
            SKINNY.addGenConstrOr(
                state_V[rd + 1][i + 12],
                [dummy_V[i + 16], state_Z[rd + 1][i + 12]],
            )
        for i in range(8):
            SKINNY.addGenConstrAnd(
                key_V[rd][SR[i]], [state_V_sb[rd][i], state_Z[rd][SR[i]]]
            )

    for i in range(16):
        key_dep[i] = LinExpr()

    tmp = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]
    for rd in range(r_in + r_dist):
        if rd >= r_in and rd < r_in + r_dist - 1:
            for i in range(8):
                key_dep[tmp[i]].add(key_V[rd - r_in][i])
        for i in range(16):
            tmp[i] = PT[tmp[i]]

    big_M = 50
    dummy_K = SKINNY.addVars(16, vtype=GRB.INTEGER, name="dummy_K")
    dummy_B = SKINNY.addVars(16, vtype=GRB.BINARY, name="dummy_B")
    for i in range(16):
        if attack_type == "k = 3n":
            SKINNY.addGenConstrIndicator(
                dummy_B[i], True, key_dep[i], GRB.GREATER_EQUAL, 4
            )
            SKINNY.addGenConstrIndicator(
                dummy_B[i], False, key_dep[i], GRB.LESS_EQUAL, 3
            )
            SKINNY.addConstr(dummy_K[i] >= key_dep[i] - 3)
            SKINNY.addConstr(dummy_K[i] <= key_dep[i] - 3 + big_M * (1 - dummy_B[i]))
        if attack_type == "k = 2n":
            SKINNY.addGenConstrIndicator(
                dummy_B[i], True, key_dep[i], GRB.GREATER_EQUAL, 3
            )
            SKINNY.addGenConstrIndicator(
                dummy_B[i], False, key_dep[i], GRB.LESS_EQUAL, 2
            )
            SKINNY.addConstr(dummy_K[i] >= key_dep[i] - 2)
            SKINNY.addConstr(dummy_K[i] <= key_dep[i] - 2 + big_M * (1 - dummy_B[i]))
        if attack_type == "k = n":
            SKINNY.addGenConstrIndicator(
                dummy_B[i], True, key_dep[i], GRB.GREATER_EQUAL, 2
            )
            SKINNY.addGenConstrIndicator(
                dummy_B[i], False, key_dep[i], GRB.LESS_EQUAL, 1
            )
            SKINNY.addConstr(dummy_K[i] >= key_dep[i] - 1)
            SKINNY.addConstr(dummy_K[i] <= key_dep[i] - 1 + big_M * (1 - dummy_B[i]))
        SKINNY.addConstr(dummy_K[i] >= 0)
        SKINNY.addConstr(dummy_K[i] <= big_M * dummy_B[i])
        Key_sieve.add(dummy_K[i])


# Offline phase
def Build_distinguisher(attack_type, block_size, r_in, r_dist, ex_constr):
    # VAR X - forward differential
    state_X_l[0] = SKINNY.addVars(16, vtype=GRB.BINARY, name="state_X_l_0")
    for rd in range(r_dist):
        state_X_nl[rd] = SKINNY.addVars(
            16, vtype=GRB.BINARY, name="state_X_nl_" + str(rd)
        )
        state_X_l[rd + 1] = SKINNY.addVars(
            16, vtype=GRB.BINARY, name="state_X_l_" + str(rd + 1)
        )
        # SR
        Build_shiftrow(state_X_l, rd, state_X_nl, rd)
        # MC
        Build_mixcolumn(state_X_nl, rd, state_X_l, rd + 1, "forward_differential")

    # VAR Y - backward determination
    state_Y_l[r_dist] = SKINNY.addVars(
        16, vtype=GRB.BINARY, name="state_Y_l_" + str(r_dist)
    )
    for rd in range(r_dist - 1, -1, -1):
        state_Y_nl[rd] = SKINNY.addVars(
            16, vtype=GRB.BINARY, name="state_Y_nl_" + str(rd)
        )
        state_Y_l[rd] = SKINNY.addVars(
            16, vtype=GRB.BINARY, name="state_Y_l_" + str(rd)
        )
        # MC
        Build_mixcolumn(state_Y_nl, rd, state_Y_l, rd + 1, "backward_determination")
        # SR
        Build_shiftrow(state_Y_l, rd, state_Y_nl, rd)

    # VAR Z
    for rd in range(r_dist + 1):
        state_Z[rd] = SKINNY.addVars(16, vtype=GRB.BINARY, name="state_Z_" + str(rd))
        SKINNY.addConstrs(state_Z[rd][i] <= state_X_l[rd][i] for i in range(16))
        SKINNY.addConstrs(state_Z[rd][i] <= state_Y_l[rd][i] for i in range(16))
        SKINNY.addConstrs(
            state_X_l[rd][i] + state_Y_l[rd][i] - state_Z[rd][i] <= 1 for i in range(16)
        )
        if rd < r_dist and rd >= 0:
            for i in range(16):
                Deg.add(state_Z[rd][i])

    # Size of delta-set
    for i in range(16):
        Start.add(state_Z[0][i])

    if block_size == 64 and attack_type != "k = n":
        SKINNY.addConstr(Start >= 2)

    # Extra constraints of start and end
    if len(ex_constr) != 0:
        A = ex_constr[0]
        B = ex_constr[1]
        for i in range(16):
            SKINNY.addConstr(state_Z[0][i] == A[i])
            SKINNY.addConstr(state_Z[r_dist][i] == B[i])

    # Nontrivial
    SKINNY.addConstr(quicksum(state_X_l[0][i] for i in range(16)) >= 1)
    SKINNY.addConstr(quicksum(state_Y_l[r_dist][i] for i in range(16)) >= 1)

    # Non-full key addition
    Build_non_full_key_addition(r_dist)
    Build_key_dependent_sieve(attack_type, r_in, r_dist)


# Key-bridging technique
def Build_key_bridging(attack_type, r_in, r_dist, r_out, begin_out, kb, var_out, K_tmp):
    for i in range(16):
        kb[i] = LinExpr()

    tmp = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]
    for rd in range(r_in + r_dist + r_out):
        if rd < r_in - 1:
            for i in range(8):
                kb[tmp[SR[i]]].add(state_O_nl[rd][i])
        elif rd >= begin_out:
            for i in range(8):
                kb[tmp[i]].add(var_out[rd - r_in - r_dist][i])
        for i in range(16):
            tmp[i] = PT[tmp[i]]

    assert attack_type == "k = 3n" or attack_type == "k = 2n" or attack_type == "k = n"
    if attack_type == "k = 3n":
        alpha = SKINNY.addVars(16, vtype=GRB.BINARY)
        beta = SKINNY.addVars(16, vtype=GRB.BINARY)
        gamma = SKINNY.addVars(16, vtype=GRB.BINARY)
        for i in range(16):
            SKINNY.addGenConstrIndicator(alpha[i], True, kb[i], GRB.GREATER_EQUAL, 3)
            SKINNY.addGenConstrIndicator(alpha[i], False, kb[i], GRB.LESS_EQUAL, 2)
            SKINNY.addGenConstrIndicator(beta[i], True, kb[i], GRB.GREATER_EQUAL, 2)
            SKINNY.addGenConstrIndicator(beta[i], False, kb[i], GRB.LESS_EQUAL, 1)
            SKINNY.addGenConstrIndicator(gamma[i], True, kb[i], GRB.GREATER_EQUAL, 1)
            SKINNY.addGenConstrIndicator(gamma[i], False, kb[i], GRB.EQUAL, 0)
            K_tmp.add(alpha[i])
            K_tmp.add(beta[i])
            K_tmp.add(gamma[i])
    if attack_type == "k = 2n":
        alpha = SKINNY.addVars(16, vtype=GRB.BINARY)
        beta = SKINNY.addVars(16, vtype=GRB.BINARY)
        for i in range(16):
            SKINNY.addGenConstrIndicator(alpha[i], True, kb[i], GRB.GREATER_EQUAL, 2)
            SKINNY.addGenConstrIndicator(alpha[i], False, kb[i], GRB.LESS_EQUAL, 1)
            SKINNY.addGenConstrIndicator(beta[i], True, kb[i], GRB.GREATER_EQUAL, 1)
            SKINNY.addGenConstrIndicator(beta[i], False, kb[i], GRB.EQUAL, 0)
            K_tmp.add(alpha[i])
            K_tmp.add(beta[i])
    if attack_type == "k = n":
        alpha = SKINNY.addVars(16, vtype=GRB.BINARY)
        for i in range(16):
            SKINNY.addGenConstrIndicator(alpha[i], True, kb[i], GRB.GREATER_EQUAL, 1)
            SKINNY.addGenConstrIndicator(alpha[i], False, kb[i], GRB.EQUAL, 0)
            K_tmp.add(alpha[i])


# Online phase
def Build_key_recovery(attack_type, r_in, r_dist, r_out):
    # VAR W - forward determination
    state_W_l[0] = SKINNY.addVars(16, vtype=GRB.BINARY, name="state_W_l_0")
    SKINNY.addConstrs(state_W_l[0][i] == state_Z[r_dist][i] for i in range(16))
    for rd in range(r_out):
        state_W_nl[rd] = SKINNY.addVars(
            16, vtype=GRB.BINARY, name="state_W_nl_" + str(rd)
        )
        state_W_l[rd + 1] = SKINNY.addVars(
            16, vtype=GRB.BINARY, name="state_W_l_" + str(rd + 1)
        )
        # SR
        Build_shiftrow(state_W_l, rd, state_W_nl, rd)
        # MC
        Build_mixcolumn(state_W_nl, rd, state_W_l, rd + 1, "forward_determination")

    # VAR M - backward differential
    state_M_l[r_in] = SKINNY.addVars(
        16, vtype=GRB.BINARY, name="state_M_l_" + str(r_in)
    )
    SKINNY.addConstrs(state_M_l[r_in][i] == state_Z[0][i] for i in range(16))
    for rd in range(r_in - 1, -1, -1):
        state_M_nl[rd] = SKINNY.addVars(
            16, vtype=GRB.BINARY, name="state_M_nl_" + str(rd)
        )
        state_M_l[rd] = SKINNY.addVars(
            16, vtype=GRB.BINARY, name="state_M_l_" + str(rd)
        )
        # MC
        Build_mixcolumn(state_M_nl, rd, state_M_l, rd + 1, "backward_differential")
        # SR
        Build_shiftrow(state_M_l, rd, state_M_nl, rd)

    # Data
    for i in range(16):
        Plain_diff.add(state_M_l[0][i])

    # VAR O - backward determination
    if r_in > 0:
        # state_O_nl[r_in - 1] = SKINNY.addVars(
        #     16, vtype=GRB.BINARY, name="state_O_nl_" + str(r_in - 1)
        # )
        state_O_l[r_in - 1] = SKINNY.addVars(
            16, vtype=GRB.BINARY, name="state_O_l_" + str(r_in - 1)
        )
        # SKINNY.addConstrs(
        #     state_O_nl[r_in - 1][i] == state_M_nl[r_in - 1][i] for i in range(16)
        # )
        SKINNY.addConstrs(
            state_O_l[r_in - 1][i] == state_M_l[r_in - 1][i] for i in range(16)
        )
    for rd in range(r_in - 2, -1, -1):
        state_O_nl[rd] = SKINNY.addVars(
            16, vtype=GRB.BINARY, name="state_O_nl_" + str(rd)
        )
        state_O_l[rd] = SKINNY.addVars(
            16, vtype=GRB.BINARY, name="state_O_l_" + str(rd)
        )
        # MC
        Build_mixcolumn(state_O_nl, rd, state_O_l, rd + 1, "backward_determination")
        # SR & VAR M
        SKINNY.addConstrs(state_O_l[rd][i] >= state_M_l[rd][i] for i in range(16))
        SKINNY.addConstrs(state_O_l[rd][SR[i]] >= state_O_nl[rd][i] for i in range(16))
        SKINNY.addConstrs(
            state_O_l[rd][SR[i]] <= state_M_l[rd][SR[i]] + state_O_nl[rd][i]
            for i in range(16)
        )

    # Key-Bridging
    Build_key_bridging(
        attack_type, r_in, r_dist, r_out, r_in + r_dist + 1, key_bri, state_W_l, Key
    )


# Value Constraints
def Build_value_constraint(attack_type, r_in, r_dist, r_out):

    SKINNY.addConstr(Val_Con >= 0)

    # Limited by the table lookup method of distinguisher

    # Limited by the table lookup method of key recovery
    # VAR D - guessed part
    state_D_l[r_out] = SKINNY.addVars(
        16, vtype=GRB.BINARY, name="state_D_l_" + str(r_out)
    )
    SKINNY.addConstrs(state_D_l[r_out][i] == 1 for i in range(16))
    for rd in range(r_out - 1, -1, -1):
        state_D_nl[rd] = SKINNY.addVars(
            16, vtype=GRB.BINARY, name="state_D_nl_" + str(rd)
        )
        state_D_l[rd] = SKINNY.addVars(
            16, vtype=GRB.BINARY, name="state_D_l_" + str(rd)
        )
        # MC
        for i in range(4):
            SKINNY.addConstr(state_D_nl[rd][i] == state_D_l[rd + 1][i + 4])
            SKINNY.addGenConstrAnd(
                state_D_nl[rd][i + 4],
                [
                    state_D_l[rd + 1][i + 4],
                    state_D_l[rd + 1][i + 8],
                    state_D_l[rd + 1][i + 12],
                ],
            )
            SKINNY.addGenConstrAnd(
                state_D_nl[rd][i + 8],
                [
                    state_D_l[rd + 1][i + 4],
                    state_D_l[rd + 1][i + 12],
                ],
            )
            SKINNY.addGenConstrAnd(
                state_D_nl[rd][i + 12],
                [
                    state_D_l[rd + 1][i],
                    state_D_l[rd + 1][i + 12],
                ],
            )
        # SR
        key_D[rd] = SKINNY.addVars(8, vtype=GRB.BINARY, name="key_D_" + str(rd))
        for i in range(8):
            SKINNY.addGenConstrAnd(
                state_D_l[rd][SR[i]], [state_D_nl[rd][i], key_D[rd][SR[i]]]
            )
        SKINNY.addConstrs(
            state_D_l[rd][SR[i]] == state_D_nl[rd][i] for i in range(8, 16)
        )

    # Calculate Kg
    Build_key_bridging(
        attack_type, r_in, r_dist, r_out, r_in + r_dist + 1, key_bri_g, key_D, Kg
    )

    # Calculate St, Kt
    # VAR H - table part
    state_H_l[0] = SKINNY.addVars(16, vtype=GRB.BINARY, name="state_H_l_0")
    SKINNY.addConstrs(state_H_l[0][i] == state_Z[r_dist][i] for i in range(16))
    for rd in range(r_out):
        state_H_nl[rd] = SKINNY.addVars(
            16, vtype=GRB.BINARY, name="state_H_nl_" + str(rd)
        )
        state_H_l[rd + 1] = SKINNY.addVars(
            16, vtype=GRB.BINARY, name="state_H_l_" + str(rd + 1)
        )
        key_H[rd] = SKINNY.addVars(8, vtype=GRB.BINARY, name="key_H_" + str(rd))
        # SR
        for i in range(16):
            SKINNY.addConstr(state_H_nl[rd][i] <= 1 - state_D_l[rd][SR[i]])
            SKINNY.addConstr(state_H_nl[rd][i] <= state_H_l[rd][SR[i]])
            SKINNY.addConstr(
                state_H_nl[rd][i] >= state_H_l[rd][SR[i]] - state_D_l[rd][SR[i]]
            )
            if i < 8:
                SKINNY.addConstr(key_H[rd][SR[i]] == state_H_nl[rd][i])
                Kt.add(key_H[rd][i])
        # MC
        tmp_nl = SKINNY.addVars(16, vtype=GRB.BINARY)
        for i in range(16):
            SKINNY.addConstr(tmp_nl[i] <= 1 - state_D_nl[rd][i])
            SKINNY.addConstr(tmp_nl[i] <= state_H_nl[rd][i])
            SKINNY.addConstr(tmp_nl[i] >= state_H_nl[rd][i] - state_D_nl[rd][i])
        # MC - COL1
        SKINNY.addConstrs(state_H_l[rd + 1][i] == tmp_nl[i + 12] for i in range(4))
        # MC - COL2
        SKINNY.addConstrs(state_H_l[rd + 1][i + 4] >= tmp_nl[i] for i in range(4))
        SKINNY.addConstrs(state_H_l[rd + 1][i + 4] >= tmp_nl[i + 4] for i in range(4))
        SKINNY.addConstrs(state_H_l[rd + 1][i + 4] >= tmp_nl[i + 8] for i in range(4))
        SKINNY.addConstrs(
            state_H_l[rd + 1][i + 4] <= tmp_nl[i] + tmp_nl[i + 4] + tmp_nl[i + 8]
            for i in range(4)
        )
        # MC - COL3
        SKINNY.addConstrs(state_H_l[rd + 1][i + 8] == tmp_nl[i + 4] for i in range(4))
        # MC - COL4
        SKINNY.addConstrs(state_H_l[rd + 1][i + 12] >= tmp_nl[i + 4] for i in range(4))
        SKINNY.addConstrs(state_H_l[rd + 1][i + 12] >= tmp_nl[i + 8] for i in range(4))
        SKINNY.addConstrs(state_H_l[rd + 1][i + 12] >= tmp_nl[i + 12] for i in range(4))
        SKINNY.addConstrs(
            state_H_l[rd + 1][i + 12] <= tmp_nl[i + 4] + tmp_nl[i + 8] + tmp_nl[i + 12]
            for i in range(4)
        )

    for rd in range(r_out + 1):
        for i in range(16):
            res_a = SKINNY.addVar(vtype=GRB.BINARY)
            SKINNY.addGenConstrAnd(res_a, [state_H_l[rd][i], state_D_l[rd][i]])
            St.add(res_a)
            if rd < r_out:
                res_b = SKINNY.addVar(vtype=GRB.BINARY)
                SKINNY.addGenConstrAnd(res_b, [state_H_nl[rd][i], state_D_nl[rd][i]])
                St.add(res_b)

    # Memory Limit
    SKINNY.addConstr((Val_Con + 1) * St + Kt <= Deg - Con - Key_sieve)


# Objective function
def Set_objective(attack_type, ns, mx_obj, mx_off, mx_on, mx_data):
    if mx_obj != 0:
        SKINNY.addConstr(Obj <= mx_obj)
    if mx_on != 0:
        SKINNY.addConstr(Obj_Online <= mx_on)
    # else:
    #     if attack_type == "k = 3n":
    #         SKINNY.addConstr(Obj_Online <= 47)
    #     if attack_type == "k = 2n":
    #         SKINNY.addConstr(Obj_Online <= 31)
    #     if attack_type == "k = n":
    #         SKINNY.addConstr(Obj_Online <= 15)
    if mx_off != 0:
        SKINNY.addConstr(Obj_Offline <= mx_off)
    if mx_data != 0:
        SKINNY.addConstr(Obj_Data <= mx_data)
    else:
        SKINNY.addConstr(Obj_Data <= 16 * ns - 1)

    # Obj_Offline
    SKINNY.addConstr(Obj_Offline == 4 + ns * (Deg - Con - Key_sieve - Val_Con))
    # Obj_Online
    SKINNY.addConstr(Obj_Online >= ns * ((Val_Con + 1) * St + Kt - Val_Con) - 3)
    SKINNY.addConstr(Obj_Online >= ns * (Val_Con + Kg))
    SKINNY.addConstr(Obj_Online >= 4 + ns * Key)
    # Obj_Data: Consider DTM trade-off
    SKINNY.addConstr(Obj_Data >= ns * Plain_diff)
    SKINNY.addConstr(
        2 * Obj_Data >= 2 * ns * (Start + Val_Con) + (Obj_Offline - Obj_Online)
    )
    # Obj
    SKINNY.addConstr(Obj >= Obj_Online)
    SKINNY.addConstr(2 * Obj >= Obj_Offline + Obj_Online)
    # Objective
    SKINNY.setObjectiveN(Obj, index=0, priority=6, name="Obj_Complexity", weight=1.0)
    SKINNY.setObjectiveN(Obj_Data, index=1, priority=5, name="Obj_Data", weight=1.0)
    SKINNY.setObjectiveN(Obj_Offline, index=2, priority=4, name="Obj_Deg", weight=1.0)
    SKINNY.setObjectiveN(Kg, index=3, priority=3, name="Obj_Kg", weight=1.0)
    SKINNY.setObjectiveN(St, index=4, priority=2, name="Obj_St", weight=1.0)
    SKINNY.setObjectiveN(Kt, index=5, priority=1, name="Obj_Kt", weight=1.0)
    SKINNY.ModelSense = GRB.MINIMIZE
    # SKINNY.setObjective(Obj, GRB.MINIMIZE)


# Start solver
def Start_solver():
    # SKINNY.Params.OutputFlag = 0
    # SKINNY.Params.Threads = 192
    SKINNY.Params.PoolSearchMode = 2
    SKINNY.Params.PoolSolutions = 1
    SKINNY.Params.PoolGap = 2.0
    # SKINNY.Params.TimeLimit = 200
    # SKINNY.setParam("IntFeasTol", 1e-9)
    SKINNY.optimize()


def Print_result(r_in, r_dist, r_out):
    DEFAULT_OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    file = DEFAULT_OUTPUT_PATH.open("w", encoding="utf-8")
    sys.stdout = file
    print("Model Status:", SKINNY.Status)
    if SKINNY.Status == 2 or SKINNY.Status == 9:
        print("Min_Obj: %g" % SKINNY.ObjVal)

        # All solutions
        for k in range(SKINNY.SolCount):
            SKINNY.Params.SolutionNumber = k
            print(
                "******** Sol no.{}    Obj = {}    Offline = {}    Online = {}    Data = {} ********".format(
                    k + 1,
                    round(Obj.Xn),
                    round(Obj_Offline.Xn),
                    round(Obj_Online.Xn),
                    round(Obj_Data.Xn),
                )
            )
            print(
                "******** Deg = {} - {} - {} = {}    Key = {}    Plaintext Difference = {} ********".format(
                    round(Evaluate_expr(Deg)),
                    round(Evaluate_expr(Con)),
                    round(Evaluate_expr(Key_sieve)),
                    round(
                        Evaluate_expr(Deg)
                        - Evaluate_expr(Con)
                        - Evaluate_expr(Key_sieve)
                    ),
                    round(Evaluate_expr(Key)),
                    round(Evaluate_expr(Plain_diff)),
                )
            )
            print(
                "******** Val_Con = {}    Kg = {}    St = {}    Kt = {} ********".format(
                    round(Val_Con.Xn),
                    round(Evaluate_expr(Kg)),
                    round(Evaluate_expr(St)),
                    round(Evaluate_expr(Kt)),
                )
            )
            strr = "A = [ "
            for i in range(16):
                if round(state_Z[0][i].Xn) == 1:
                    strr += str(i) + " "
            strr += "]    B = [ "
            for i in range(16):
                if round(state_Z[r_dist][i].Xn) == 1:
                    strr += str(i) + " "
            strr += "]"
            print(strr)

            print("---------- Var X ----------")
            for rd in range(r_dist + 1):
                print("X_linear[", rd, "]")
                for i in range(4):
                    print(
                        round(state_X_l[rd][4 * i].Xn),
                        round(state_X_l[rd][4 * i + 1].Xn),
                        round(state_X_l[rd][4 * i + 2].Xn),
                        round(state_X_l[rd][4 * i + 3].Xn),
                    )
                if rd < r_dist:
                    print("X_non_linear[", rd, "]")
                    for i in range(4):
                        print(
                            round(state_X_nl[rd][4 * i].Xn),
                            round(state_X_nl[rd][4 * i + 1].Xn),
                            round(state_X_nl[rd][4 * i + 2].Xn),
                            round(state_X_nl[rd][4 * i + 3].Xn),
                        )

            print("---------- Var Y ----------")
            for rd in range(r_dist + 1):
                print("Y_linear[", rd, "]")
                for i in range(4):
                    print(
                        round(state_Y_l[rd][4 * i].Xn),
                        round(state_Y_l[rd][4 * i + 1].Xn),
                        round(state_Y_l[rd][4 * i + 2].Xn),
                        round(state_Y_l[rd][4 * i + 3].Xn),
                    )
                if rd < r_dist:
                    print("Y_non_linear[", rd, "]")
                    for i in range(4):
                        print(
                            round(state_Y_nl[rd][4 * i].Xn),
                            round(state_Y_nl[rd][4 * i + 1].Xn),
                            round(state_Y_nl[rd][4 * i + 2].Xn),
                            round(state_Y_nl[rd][4 * i + 3].Xn),
                        )

            print("---------- Var Z ----------")
            for rd in range(r_dist + 1):
                print("Z[", rd, "]")
                for i in range(4):
                    print(
                        round(state_Z[rd][4 * i].Xn),
                        round(state_Z[rd][4 * i + 1].Xn),
                        round(state_Z[rd][4 * i + 2].Xn),
                        round(state_Z[rd][4 * i + 3].Xn),
                    )

            # print("---------- Var V ----------")
            # for rd in range(r_dist):
            #     print("V_sb[", rd, "]")
            #     for i in range(4):
            #         print(
            #             round(state_V_sb[rd][4 * i].Xn),
            #             round(state_V_sb[rd][4 * i + 1].Xn),
            #             round(state_V_sb[rd][4 * i + 2].Xn),
            #             round(state_V_sb[rd][4 * i + 3].Xn),
            #         )
            #     print("V[", rd + 1, "]")
            #     for i in range(4):
            #         print(
            #             round(state_V[rd + 1][4 * i].Xn),
            #             round(state_V[rd + 1][4 * i + 1].Xn),
            #             round(state_V[rd + 1][4 * i + 2].Xn),
            #             round(state_V[rd + 1][4 * i + 3].Xn),
            #         )

            print("---------- Key V ----------")
            for rd in range(r_dist - 1):
                print("Key_V[", rd, "]")
                for i in range(2):
                    print(
                        round(key_V[rd][4 * i].Xn),
                        round(key_V[rd][4 * i + 1].Xn),
                        round(key_V[rd][4 * i + 2].Xn),
                        round(key_V[rd][4 * i + 3].Xn),
                    )

            print("---------- Var M ----------")
            for rd in range(r_in + 1):
                print("M_linear[", rd, "]")
                for i in range(4):
                    print(
                        round(state_M_l[rd][4 * i].Xn),
                        round(state_M_l[rd][4 * i + 1].Xn),
                        round(state_M_l[rd][4 * i + 2].Xn),
                        round(state_M_l[rd][4 * i + 3].Xn),
                    )
                if rd < r_in:
                    print("M_non_linear[", rd, "]")
                    for i in range(4):
                        print(
                            round(state_M_nl[rd][4 * i].Xn),
                            round(state_M_nl[rd][4 * i + 1].Xn),
                            round(state_M_nl[rd][4 * i + 2].Xn),
                            round(state_M_nl[rd][4 * i + 3].Xn),
                        )

            print("---------- Var O ----------")
            for rd in range(r_in):
                print("O_linear[", rd, "]")
                for i in range(4):
                    print(
                        round(state_O_l[rd][4 * i].Xn),
                        round(state_O_l[rd][4 * i + 1].Xn),
                        round(state_O_l[rd][4 * i + 2].Xn),
                        round(state_O_l[rd][4 * i + 3].Xn),
                    )
                if rd == r_in - 1:
                    break
                print("O_non_linear[", rd, "]")
                for i in range(4):
                    print(
                        round(state_O_nl[rd][4 * i].Xn),
                        round(state_O_nl[rd][4 * i + 1].Xn),
                        round(state_O_nl[rd][4 * i + 2].Xn),
                        round(state_O_nl[rd][4 * i + 3].Xn),
                    )

            print("---------- Var W ----------")
            for rd in range(r_out + 1):
                print("W_linear[", rd, "]")
                for i in range(4):
                    print(
                        round(state_W_l[rd][4 * i].Xn),
                        round(state_W_l[rd][4 * i + 1].Xn),
                        round(state_W_l[rd][4 * i + 2].Xn),
                        round(state_W_l[rd][4 * i + 3].Xn),
                    )
                if rd < r_out:
                    print("W_non_linear[", rd, "]")
                    for i in range(4):
                        print(
                            round(state_W_nl[rd][4 * i].Xn),
                            round(state_W_nl[rd][4 * i + 1].Xn),
                            round(state_W_nl[rd][4 * i + 2].Xn),
                            round(state_W_nl[rd][4 * i + 3].Xn),
                        )

            print("---------- Var D ----------")
            for rd in range(r_out + 1):
                print("D_linear[", rd, "]")
                for i in range(4):
                    print(
                        round(state_D_l[rd][4 * i].Xn),
                        round(state_D_l[rd][4 * i + 1].Xn),
                        round(state_D_l[rd][4 * i + 2].Xn),
                        round(state_D_l[rd][4 * i + 3].Xn),
                    )
                if rd < r_out:
                    print("D_non_linear[", rd, "]")
                    for i in range(4):
                        print(
                            round(state_D_nl[rd][4 * i].Xn),
                            round(state_D_nl[rd][4 * i + 1].Xn),
                            round(state_D_nl[rd][4 * i + 2].Xn),
                            round(state_D_nl[rd][4 * i + 3].Xn),
                        )

            print("---------- Key D ----------")
            for rd in range(r_out):
                print("Key_D[", rd, "]")
                for i in range(2):
                    print(
                        round(key_D[rd][4 * i].Xn),
                        round(key_D[rd][4 * i + 1].Xn),
                        round(key_D[rd][4 * i + 2].Xn),
                        round(key_D[rd][4 * i + 3].Xn),
                    )

            print("---------- Var H ----------")
            for rd in range(r_out + 1):
                print("H_linear[", rd, "]")
                for i in range(4):
                    print(
                        round(state_H_l[rd][4 * i].Xn),
                        round(state_H_l[rd][4 * i + 1].Xn),
                        round(state_H_l[rd][4 * i + 2].Xn),
                        round(state_H_l[rd][4 * i + 3].Xn),
                    )
                if rd < r_out:
                    print("H_non_linear[", rd, "]")
                    for i in range(4):
                        print(
                            round(state_H_nl[rd][4 * i].Xn),
                            round(state_H_nl[rd][4 * i + 1].Xn),
                            round(state_H_nl[rd][4 * i + 2].Xn),
                            round(state_H_nl[rd][4 * i + 3].Xn),
                        )

            print("---------- Key H ----------")
            for rd in range(r_out):
                print("Key_H[", rd, "]")
                for i in range(2):
                    print(
                        round(key_H[rd][4 * i].Xn),
                        round(key_H[rd][4 * i + 1].Xn),
                        round(key_H[rd][4 * i + 2].Xn),
                        round(key_H[rd][4 * i + 3].Xn),
                    )

            print("---------- Con ----------")
            for rd in range(r_dist):
                print(
                    "con[",
                    rd,
                    "] :",
                    round(
                        con_1[rd][0].Xn
                        + con_2[rd][0].Xn
                        + con_3[rd][0].Xn
                        - con_4[rd][0].Xn
                    ),
                    round(
                        con_1[rd][1].Xn
                        + con_2[rd][1].Xn
                        + con_3[rd][1].Xn
                        - con_4[rd][1].Xn
                    ),
                    round(
                        con_1[rd][2].Xn
                        + con_2[rd][2].Xn
                        + con_3[rd][2].Xn
                        - con_4[rd][2].Xn
                    ),
                    round(
                        con_1[rd][3].Xn
                        + con_2[rd][3].Xn
                        + con_3[rd][3].Xn
                        - con_4[rd][3].Xn
                    ),
                )

            print("---------- Key-Bridge ----------")
            strr = ""
            for i in range(16):
                strr += str(round(Evaluate_expr(key_bri[i]))) + " "
            print(strr)

            print("---------- Key-Dependent-Seive ----------")
            strr = ""
            for i in range(16):
                strr += str(round(Evaluate_expr(key_dep[i]))) + " "
            print(strr)

    file.close()
    sys.stdout = sys.__stdout__


def Search_ds_mitm_attack(
    attack_type,
    block_size,
    r_dist,
    r_in,
    r_out,
    mx_obj,
    mx_off,
    mx_on,
    mx_data,
    ex_constr,
):
    Reset_model(r_in, r_dist, r_out)
    Build_distinguisher(attack_type, block_size, r_in, r_dist, ex_constr)
    Build_key_recovery(attack_type, r_in, r_dist, r_out)
    Build_value_constraint(attack_type, r_in, r_dist, r_out)
    Set_objective(attack_type, block_size // 16, mx_obj, mx_off, mx_on, mx_data)
    Start_solver()
    Print_result(r_in, r_dist, r_out)
    if SKINNY.Status == 2:
        print("Min_Obj: %g" % SKINNY.ObjVal)
        print("# Value Constraint = {}".format(round(Val_Con.Xn)))
        return True
    return False


if __name__ == "__main__":
    # attack_type, block_size, r_dist, r_in, r_out, mx_obj, mx_off, mx_on, mx_data, [A, B]

    A = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0]
    # B = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    B1 = [1, 1, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 1, 0]
    B2 = [0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1]
    Search_ds_mitm_attack(
        attack_type="k = 3n",
        block_size=64,
        r_dist=9,
        r_in=3,
        r_out=10,
        mx_obj=0,
        mx_off=0,
        mx_on=0,
        mx_data=0,
        ex_constr=[],
    )

    # Exploiting Non-Full Key Additions: Full-Fledged Automatic Demirci-Selcuk Meet-in-the-Middle Cryptanalysis of SKINNY
    # Search_ds_mitm_attack(9, 4, 10, 0, 0, 0, [])

    # Automatic DS Meet-in-the-Middle Attack on SKINNY with Key-Bridging
    # Search_ds_mitm_attack(10, 3, 9, 0, 0, 0, [])
