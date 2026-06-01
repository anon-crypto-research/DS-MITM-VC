from gurobipy import *
from pathlib import Path
import time
import sys

MODULE_DIR = Path(__file__).resolve().parent
DEFAULT_OUTPUT_PATH = MODULE_DIR / "results" / "output"


PIE = [5, 0, 1, 4, 7, 12, 3, 8, 13, 6, 9, 2, 15, 10, 11, 14]
ROT1_80 = [4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 1, 2, 3, 0]
ROT2_80 = [19, 16, 17, 18, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]
ROT1_128 = [
    4,
    5,
    6,
    7,
    8,
    9,
    10,
    11,
    12,
    13,
    14,
    15,
    16,
    17,
    18,
    19,
    20,
    21,
    22,
    23,
    24,
    25,
    26,
    27,
    28,
    29,
    30,
    31,
    1,
    2,
    3,
    0,
]
ROT2_128 = [
    31,
    28,
    29,
    30,
    0,
    1,
    2,
    3,
    4,
    5,
    6,
    7,
    8,
    9,
    10,
    11,
    12,
    13,
    14,
    15,
    16,
    17,
    18,
    19,
    20,
    21,
    22,
    23,
    24,
    25,
    26,
    27,
]
IDX_80 = [1, 3, 4, 6, 13, 14, 15, 16]
IDX_128 = [2, 3, 12, 15, 17, 18, 28, 31]


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


def Reset_model():
    # l -(ARK, Sbox, XOR)-> nl -(Permutation)-> l
    global TWINE
    TWINE = Model("TWINE")

    global state_X_l, state_X_nl, state_Y_l, state_Y_nl
    global state_Z_l, state_Z_nl, state_Zg
    state_X_l = {}  # VAR X
    state_X_nl = {}
    state_Y_l = {}  # VAR Y
    state_Y_nl = {}
    state_Z_l = {}  # VAR Z
    state_Z_nl = {}
    state_Zg = {}

    global state_W_l, state_W_nl, state_M_l, state_M_nl
    global state_O_l, state_O_nl, key_V
    state_W_l = {}  # VAR W
    state_W_nl = {}
    state_M_l = {}  # VAR M
    state_M_nl = {}
    state_O_l = {}  # VAR O
    state_O_nl = {}
    key_V = {}

    global kb_state, kb_path
    kb_state = {}
    kb_path = {}

    global Deg, Start, End, Data, Key, Obj_Online, Val_Con, Obj
    Deg = LinExpr()
    Start = LinExpr()
    End = LinExpr()
    Data = LinExpr()
    Key = LinExpr()
    Obj_Online = TWINE.addVar(vtype=GRB.INTEGER, name="Obj_Online")
    Val_Con = TWINE.addVar(vtype=GRB.INTEGER, name="Val_Con")
    Obj = TWINE.addVar(vtype=GRB.INTEGER, name="Obj")


def Build_permutation(a, b):
    # a -(Permutation)-> b
    TWINE.addConstrs(b[PIE[i]] == a[i] for i in range(16))


def Build_F(a, b, typ):
    # b_i = a_i
    # b_{i+1} = a_i xor a_{i+1}
    assert typ == "determination" or typ == "differential"
    if typ == "differential":  # a active -> b active (diffusion)
        for i in range(8):
            TWINE.addGenConstrOr(b[2 * i + 1], [a[2 * i], a[2 * i + 1]])
            TWINE.addConstr(b[2 * i] == a[2 * i])
    if typ == "determination":  # To determine a, we need the value of b
        for i in range(8):
            TWINE.addGenConstrOr(b[2 * i], [a[2 * i], a[2 * i + 1]])
            TWINE.addConstr(b[2 * i + 1] == a[2 * i + 1])


# Key-dependent-sieve technique
# def Build_key_dependent_sieve():


# Offline phase
def Build_distinguisher(r_dist, ex_constr):
    # VAR X - differential
    state_X_l[0] = TWINE.addVars(16, vtype=GRB.BINARY, name="state_X_l_0")
    for rd in range(r_dist):
        state_X_nl[rd] = TWINE.addVars(
            16, vtype=GRB.BINARY, name="state_X_nl_" + str(rd)
        )
        state_X_l[rd + 1] = TWINE.addVars(
            16, vtype=GRB.BINARY, name="state_X_l_" + str(rd + 1)
        )
        # F
        Build_F(state_X_l[rd], state_X_nl[rd], "differential")
        # Permutation
        Build_permutation(state_X_nl[rd], state_X_l[rd + 1])

    # VAR Y - determination
    state_Y_l[r_dist] = TWINE.addVars(
        16, vtype=GRB.BINARY, name="state_Y_l_" + str(r_dist)
    )
    for rd in range(r_dist - 1, -1, -1):
        state_Y_nl[rd] = TWINE.addVars(
            16, vtype=GRB.BINARY, name="state_Y_nl_" + str(rd)
        )
        state_Y_l[rd] = TWINE.addVars(16, vtype=GRB.BINARY, name="state_Y_l_" + str(rd))
        # Permutation
        Build_permutation(state_Y_nl[rd], state_Y_l[rd + 1])
        # F
        Build_F(state_Y_nl[rd], state_Y_l[rd], "determination")

    # VAR Z
    for rd in range(r_dist + 1):
        state_Z_l[rd] = TWINE.addVars(16, vtype=GRB.BINARY, name="state_Z_l_" + str(rd))
        for i in range(16):
            TWINE.addGenConstrAnd(
                state_Z_l[rd][i], [state_X_l[rd][i], state_Y_l[rd][i]]
            )

    for rd in range(r_dist):
        state_Z_nl[rd] = TWINE.addVars(
            16, vtype=GRB.BINARY, name="state_Z_nl_" + str(rd)
        )
        for i in range(16):
            TWINE.addGenConstrAnd(
                state_Z_nl[rd][i], [state_X_nl[rd][i], state_Y_nl[rd][i]]
            )
        state_Zg[rd] = TWINE.addVars(8, vtype=GRB.BINARY, name="state_Zg_" + str(rd))
        for i in range(8):
            TWINE.addGenConstrAnd(
                state_Zg[rd][i], [state_Z_l[rd][2 * i], state_Z_nl[rd][2 * i + 1]]
            )
            Deg.add(state_Zg[rd][i])

    # Size of delta-set
    for i in range(16):
        Start.add(state_Z_l[0][i])
        End.add(state_Z_l[r_dist][i])

    # Nontrivial and effective
    TWINE.addConstr(Start >= 1)
    TWINE.addConstr(End >= 1)
    # TWINE.addConstr(Start + End >= 3)

    # Extra constraint at the head and tail
    if len(ex_constr) != 0:
        A = ex_constr[0]
        B = ex_constr[1]
        if any(A):
            TWINE.addConstrs(state_X_l[0][i] == A[i] for i in range(16))
        if any(B):
            TWINE.addConstrs(state_Y_l[r_dist][i] == B[i] for i in range(16))


def calc_80(rd, idx):
    return 20 * rd + idx


def calc_128(rd, idx):
    return 32 * rd + idx


# Key-bridging technique - TWINE-80
def Build_key_bridging_80(r_in, r_dist, r_out):
    R = r_in + r_dist + r_out
    n = 20 * R
    beta = 10  # Max derivation step
    kb_state[0] = TWINE.addVars(n, vtype=GRB.BINARY, name="kb_state_0")
    for i in range(beta):  # For each step
        kb_state[i + 1] = TWINE.addVars(
            n, vtype=GRB.BINARY, name="kb_state_" + str(i + 1)
        )
        kb_path[i + 1] = {}
        for j in range(n):  # For each variable in k
            kb_path[i + 1][j] = TWINE.addVars(
                2, vtype=GRB.BINARY, name="kb_path_" + str(i + 1) + "_" + str(j)
            )
            rd = j // 20
            idx = j % 20
            # relation 0: from last row
            if rd > 0:
                if idx == 16:
                    TWINE.addGenConstrAnd(
                        kb_path[i + 1][j][0],
                        [
                            kb_state[i][calc_80(rd - 1, 0)],
                            kb_state[i][calc_80(rd - 1, 1)],
                        ],
                    )
                elif idx == 0:
                    TWINE.addGenConstrAnd(
                        kb_path[i + 1][j][0],
                        [
                            kb_state[i][calc_80(rd - 1, 4)],
                            kb_state[i][calc_80(rd - 1, 16)],
                        ],
                    )
                else:
                    TWINE.addConstr(
                        kb_path[i + 1][j][0]
                        == kb_state[i][calc_80(rd - 1, ROT1_80[idx])]
                    )
            else:
                TWINE.addConstr(kb_path[i + 1][j][0] == 0)
            # relation 1: from next row
            if rd < R - 1:
                if idx == 1:
                    TWINE.addGenConstrAnd(
                        kb_path[i + 1][j][1],
                        [
                            kb_state[i][calc_80(rd + 1, 16)],
                            kb_state[i][calc_80(rd + 1, 19)],
                        ],
                    )
                elif idx == 4:
                    TWINE.addGenConstrAnd(
                        kb_path[i + 1][j][1],
                        [
                            kb_state[i][calc_80(rd + 1, 0)],
                            kb_state[i][calc_80(rd + 1, 12)],
                        ],
                    )
                else:
                    TWINE.addConstr(
                        kb_path[i + 1][j][1]
                        == kb_state[i][calc_80(rd + 1, ROT2_80[idx])]
                    )
            else:
                TWINE.addConstr(kb_path[i + 1][j][1] == 0)
            TWINE.addGenConstrOr(
                kb_state[i + 1][j],
                [kb_state[i][j], kb_path[i + 1][j][0], kb_path[i + 1][j][1]],
            )
    for rd in range(R):
        if rd < r_in:
            for i in range(8):
                TWINE.addConstr(kb_state[beta][calc_80(rd, IDX_80[i])] >= key_V[rd][i])
        if rd >= r_in + r_dist:
            for i in range(8):
                TWINE.addConstr(
                    kb_state[beta][calc_80(rd, IDX_80[i])]
                    >= state_W_l[rd - r_in - r_dist][2 * i + 1]
                )
    Key.add(quicksum(kb_state[0][i] for i in range(n)))


# Key-bridging technique - TWINE-128
def Build_key_bridging_128(r_in, r_dist, r_out):
    R = r_in + r_dist + r_out
    n = 32 * R
    beta = 0  # Max derivation step
    kb_state[0] = TWINE.addVars(n, vtype=GRB.BINARY, name="kb_state_0")
    for i in range(beta):  # For each step
        kb_state[i + 1] = TWINE.addVars(
            n, vtype=GRB.BINARY, name="kb_state_" + str(i + 1)
        )
        kb_path[i + 1] = {}
        for j in range(n):  # For each variable in k
            kb_path[i + 1][j] = TWINE.addVars(
                2, vtype=GRB.BINARY, name="kb_path_" + str(i + 1) + "_" + str(j)
            )
            rd = j // 32
            idx = j % 32
            # relation 0: from last row
            if rd > 0:
                if idx == 28:
                    TWINE.addGenConstrAnd(
                        kb_path[i + 1][j][0],
                        [
                            kb_state[i][calc_128(rd - 1, 0)],
                            kb_state[i][calc_128(rd - 1, 1)],
                        ],
                    )
                elif idx == 0:
                    TWINE.addGenConstrAnd(
                        kb_path[i + 1][j][0],
                        [
                            kb_state[i][calc_128(rd - 1, 4)],
                            kb_state[i][calc_128(rd - 1, 16)],
                        ],
                    )
                elif idx == 19:
                    TWINE.addGenConstrAnd(
                        kb_path[i + 1][j][0],
                        [
                            kb_state[i][calc_128(rd - 1, 23)],
                            kb_state[i][calc_128(rd - 1, 30)],
                        ],
                    )
                else:
                    TWINE.addConstr(
                        kb_path[i + 1][j][0]
                        == kb_state[i][calc_128(rd - 1, ROT1_128[idx])]
                    )
            else:
                TWINE.addConstr(kb_path[i + 1][j][0] == 0)
            # relation 1: from next row
            if rd < R - 1:
                if idx == 1:
                    TWINE.addGenConstrAnd(
                        kb_path[i + 1][j][1],
                        [
                            kb_state[i][calc_128(rd + 1, 28)],
                            kb_state[i][calc_128(rd + 1, 31)],
                        ],
                    )
                elif idx == 4:
                    TWINE.addGenConstrAnd(
                        kb_path[i + 1][j][1],
                        [
                            kb_state[i][calc_128(rd + 1, 0)],
                            kb_state[i][calc_128(rd + 1, 12)],
                        ],
                    )
                elif idx == 23:
                    TWINE.addGenConstrAnd(
                        kb_path[i + 1][j][1],
                        [
                            kb_state[i][calc_128(rd + 1, 19)],
                            kb_state[i][calc_128(rd + 1, 26)],
                        ],
                    )
                else:
                    TWINE.addConstr(
                        kb_path[i + 1][j][1]
                        == kb_state[i][calc_128(rd + 1, ROT2_128[idx])]
                    )
            else:
                TWINE.addConstr(kb_path[i + 1][j][1] == 0)
            TWINE.addGenConstrOr(
                kb_state[i + 1][j],
                [kb_state[i][j], kb_path[i + 1][j][0], kb_path[i + 1][j][1]],
            )
    for rd in range(R):
        if rd < r_in:
            for i in range(8):
                TWINE.addConstr(
                    kb_state[beta][calc_128(rd, IDX_128[i])] >= key_V[rd][i]
                )
        if rd >= r_in + r_dist:
            for i in range(8):
                TWINE.addConstr(
                    kb_state[beta][calc_128(rd, IDX_128[i])]
                    >= state_W_l[rd - r_in - r_dist][2 * i + 1]
                )
    Key.add(quicksum(kb_state[0][i] for i in range(n)))


# Online phase
def Build_key_recovery(attack_type, r_in, r_dist, r_out):
    assert attack_type == "TWINE-80" or attack_type == "TWINE-128"
    # VAR W - determination
    state_W_l[0] = TWINE.addVars(16, vtype=GRB.BINARY, name="state_W_l_0")
    TWINE.addConstrs(state_W_l[0][i] == state_Z_l[r_dist][i] for i in range(16))
    for rd in range(r_out):
        state_W_nl[rd] = TWINE.addVars(
            16, vtype=GRB.BINARY, name="state_W_nl_" + str(rd)
        )
        state_W_l[rd + 1] = TWINE.addVars(
            16, vtype=GRB.BINARY, name="state_W_l_" + str(rd + 1)
        )
        # F
        Build_F(state_W_l[rd], state_W_nl[rd], "determination")
        # Permutation
        Build_permutation(state_W_nl[rd], state_W_l[rd + 1])

    # VAR M - differential
    state_M_l[r_in] = TWINE.addVars(16, vtype=GRB.BINARY, name="state_M_l_" + str(r_in))
    TWINE.addConstrs(state_M_l[r_in][i] == state_Z_l[0][i] for i in range(16))
    for rd in range(r_in - 1, -1, -1):
        state_M_nl[rd] = TWINE.addVars(
            16, vtype=GRB.BINARY, name="state_M_nl_" + str(rd)
        )
        state_M_l[rd] = TWINE.addVars(16, vtype=GRB.BINARY, name="state_M_l_" + str(rd))
        # Permutation
        Build_permutation(state_M_nl[rd], state_M_l[rd + 1])
        # F
        Build_F(state_M_nl[rd], state_M_l[rd], "differential")

    # VAR O - determination
    state_O_l[r_in] = TWINE.addVars(16, vtype=GRB.BINARY, name="state_O_l_" + str(r_in))
    TWINE.addConstrs(state_O_l[r_in][i] == 0 for i in range(16))
    for rd in range(r_in - 1, -1, -1):
        state_O_nl[rd] = TWINE.addVars(
            16, vtype=GRB.BINARY, name="state_O_nl_" + str(rd)
        )
        state_O_l[rd] = TWINE.addVars(16, vtype=GRB.BINARY, name="state_O_l_" + str(rd))
        # Permutation
        Build_permutation(state_O_nl[rd], state_O_l[rd + 1])
        # F
        for i in range(8):
            TWINE.addGenConstrOr(
                state_O_l[rd][2 * i],
                [
                    state_O_nl[rd][2 * i],
                    state_O_nl[rd][2 * i + 1],
                    state_M_nl[rd][2 * i],
                ],
            )
            TWINE.addConstr(state_O_l[rd][2 * i + 1] == state_O_nl[rd][2 * i + 1])

    # Data
    for i in range(16):
        Data.add(state_M_l[0][i])

    # Kin
    for rd in range(r_in):
        key_V[rd] = TWINE.addVars(8, vtype=GRB.BINARY, name="key_V_" + str(rd))
        for i in range(8):
            TWINE.addGenConstrOr(
                key_V[rd][i], [state_M_l[rd][2 * i], state_O_l[rd][2 * i + 1]]
            )

    # Key-Bridging
    if attack_type == "TWINE-80":
        Build_key_bridging_80(r_in, r_dist, r_out)
    if attack_type == "TWINE-128":
        Build_key_bridging_128(r_in, r_dist, r_out)


# Value Constraints
# def Build_value_constraint():


# Objective function
def Set_objective(mx_obj, mx_off, mx_on, mx_data):
    if mx_obj != 0:
        TWINE.addConstr(Obj <= mx_obj)
    if mx_on != 0:
        TWINE.addConstr(Obj_Online <= mx_on)
    if mx_off != 0:
        TWINE.addConstr(Deg <= mx_off)
    if mx_data != 0:
        TWINE.addConstr(Data <= mx_data)
    else:
        TWINE.addConstr(Data <= 15)

    # ----------
    TWINE.addConstr(Val_Con == 0)
    # TWINE.addConstr(Key >= 31)
    # TWINE.addConstr(Data >= 8)
    # TWINE.addConstr(Deg >= 14)
    # Data <= 15
    # TWINE.addConstr(Deg + Val_Con - Obj_Online <= 28)
    # ----------

    TWINE.addConstr(Obj_Online >= Start + Key)
    TWINE.addConstr(Obj >= Obj_Online)
    TWINE.addConstr(2 * Obj >= (Start + Deg - Val_Con) + Obj_Online)
    TWINE.setObjectiveN(Obj, index=0, priority=4, name="Obj_Complexity", weight=1.0)
    TWINE.setObjectiveN(Key, index=1, priority=3, name="Obj_Online", weight=1.0)
    TWINE.setObjectiveN(Data, index=2, priority=2, name="Obj_Data", weight=1.0)
    TWINE.setObjectiveN(
        Deg - Val_Con, index=3, priority=1, name="Obj_Offline", weight=1.0
    )
    TWINE.ModelSense = GRB.MINIMIZE
    # SKINNY.setObjective(Obj, GRB.MINIMIZE)
    # SKINNY.setObjective(Val_Con, GRB.MINIMIZE)


# Start solver
def Start_solver():
    # TWINE.Params.OutputFlag = 0
    # TWINE.Params.PoolSearchMode = 2
    # TWINE.Params.Threads = 192
    TWINE.Params.PoolSolutions = 1
    # TWINE.Params.PoolGap = 2.0
    # TWINE.Params.TimeLimit = 200
    # TWINE.setParam("IntFeasTol", 1e-9)
    TWINE.optimize()


def Ouput_var(var_name, rd, var):
    print(var_name + "[", rd, "]")
    for i in range(8):
        print(
            "(" + str(round(var[2 * i].Xn)) + " " + str(round(var[2 * i + 1].Xn)) + ")",
            end=" ",
        )
    print("\n")


def Print_result(r_in, r_dist, r_out):
    DEFAULT_OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    file = DEFAULT_OUTPUT_PATH.open("w", encoding="utf-8")
    sys.stdout = file
    print("Model Status:", TWINE.Status)
    if TWINE.Status != 3:
        print("Min_Obj: %g" % TWINE.ObjVal)

        # All solutions
        for k in range(TWINE.SolCount):
            TWINE.Params.SolutionNumber = k
            print(
                "******** Sol no.{}    Obj = {}    Deg = {}    Key = {}    Data = {}    Val_Con = {} ********".format(
                    k + 1,
                    round(Obj.Xn),
                    round(Evaluate_expr(Deg)),
                    round(Evaluate_expr(Key)),
                    round(Evaluate_expr(Data)),
                    round(Evaluate_expr(Val_Con)),
                )
            )
            print("A = [", end=" ")
            for i in range(16):
                if round(state_Z_l[0][i].Xn) == 1:
                    print(i, end=" ")
            print("]    B = [", end=" ")
            for i in range(16):
                if round(state_Z_l[r_dist][i].Xn) == 1:
                    print(i, end=" ")
            print("]")

            print("---------- Var X ----------")
            for rd in range(r_dist + 1):
                Ouput_var("X_linear", rd, state_X_l[rd])
                if rd < r_dist:
                    Ouput_var("X_non_linear", rd, state_X_nl[rd])

            print("---------- Var Y ----------")
            for rd in range(r_dist + 1):
                Ouput_var("Y_linear", rd, state_Y_l[rd])
                if rd < r_dist:
                    Ouput_var("Y_non_linear", rd, state_Y_nl[rd])

            print("---------- Var Z ----------")
            for rd in range(r_dist + 1):
                Ouput_var("Z_linear", rd, state_Z_l[rd])
                if rd < r_dist:
                    Ouput_var("Z_non_linear", rd, state_Z_nl[rd])

            print("---------- Var Zg ----------")
            for rd in range(r_dist):
                print("Zg[", rd, "]")
                for i in range(8):
                    print(round(state_Zg[rd][i].Xn), end=" ")
                print("\n")

            print("---------- Var M ----------")
            for rd in range(r_in + 1):
                Ouput_var("M_linear", rd, state_M_l[rd])
                if rd < r_in:
                    Ouput_var("M_non_linear", rd, state_M_nl[rd])

            print("---------- Var O ----------")
            for rd in range(r_in + 1):
                Ouput_var("O_linear", rd, state_O_l[rd])
                if rd < r_in:
                    Ouput_var("O_non_linear", rd, state_O_nl[rd])

            print("---------- Var W ----------")
            for rd in range(r_out + 1):
                Ouput_var("W_linear", rd, state_W_l[rd])
                if rd < r_out:
                    Ouput_var("W_non_linear", rd, state_W_nl[rd])

            print("---------- Round Key Guess ----------")
            for rd in range(r_in):
                rk_status = [round(key_V[rd][i].Xn) for i in range(8)]
                status_str = " ".join(map(str, rk_status))
                print(f"RK[{rd:2}] | {status_str}")
            for rd in range(r_out):
                rk_status = [round(state_W_l[rd][2 * i + 1].Xn) for i in range(8)]
                status_str = " ".join(map(str, rk_status))
                print(f"RK[{(rd+r_in+r_dist):2}] | {status_str}")

            print("---------- Key Basis ----------")
            for rd in range(r_in + r_dist + r_out):
                print("WK[{}]:".format(rd), end=" ")
                for i in range(20):
                    print(round(kb_state[0][calc_80(rd, i)].Xn), end=" ")
                print("\n")

    file.close()
    sys.stdout = sys.__stdout__


def Search_ds_mitm_attack(
    attack_type, r_dist, r_in, r_out, mx_obj, mx_off, mx_on, mx_data, ex_constr
):
    Reset_model()
    Build_distinguisher(r_dist, ex_constr)
    Build_key_recovery(attack_type, r_in, r_dist, r_out)
    Set_objective(mx_obj, mx_off, mx_on, mx_data)
    Start_solver()
    Print_result(r_in, r_dist, r_out)
    if TWINE.Status == 2:
        # print("Min_Obj: %g" % TWINE.ObjVal)
        return True
    return False


if __name__ == "__main__":
    # attack_type, r_dist, r_in, r_out, mx_obj, mx_off, mx_on, mx_data, [A, B]

    # [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    # A = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0]
    # B = [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0]
    A = [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    B = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0]
    Search_ds_mitm_attack(
        attack_type="TWINE-80",
        r_dist=11,
        r_in=4,
        r_out=6,
        mx_obj=0,
        mx_off=0,
        mx_on=0,
        mx_data=0,
        ex_constr=[A, B],
    )
