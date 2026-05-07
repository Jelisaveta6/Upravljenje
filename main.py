import os
from pulp import *
solver_path = "C:\\x64_win64\\cplex.exe" 
timeLimit = 300 
model = LpProblem("RailYardScheduling", LpMinimize) 

I_max = 5 
O_max = 5 
B_max = 5 
K_max = 5 
T_max = 120 
T = range(1, T_max+1) 
N_RK_max = 100 
Tmin = 15 
R_PRI = 4 
R_OTP = 3 
N_OV_max = 100 
N_OV_min = 50 
mi = 0.50 
H_RI = 5 
H_TP = 15 
H_FJ = 10 
M = 10000 
I = range(1, I_max+1) 
O = range(1, O_max+1) 
B = range(1, B_max+1) 
K = range(1, K_max+1) 
T = range(0, T_max+1) 

T_PRI = {1: 10, 2: 20, 3: 30, 4: 40, 5: 50} 
T_OTP = {1: 70, 2: 80, 3: 90, 4: 100, 5: 110} 
H_RJ = {1: 20, 2: 25, 3: 15, 4: 22, 5: 30} 

N_IT = {
    (1,1): 10, (2,1): 15, (3,1): 5, (4,1): 8, (5,1): 12,
    (1,2): 8, (2,2): 12, (3,2): 10, (4,2): 6, (5,2): 14,
    (1,3): 15, (2,3): 5, (3,3): 12, (4,3): 10, (5,3): 8,
    (1,4): 6, (2,4): 14, (3,4): 9, (4,4): 11, (5,4): 10,
    (1,5): 12, (2,5): 10, (3,5): 8, (4,5): 7, (5,5): 13
}

A_init = {1: 5, 2: 8, 3: 6, 4: 7, 5: 4} 

PHIkb = {
    (1,1): 1, (1,2): 0, (1,3): 0, (1,4): 0, (1,5): 0,
    (2,1): 0, (2,2): 1, (2,3): 0, (2,4): 0, (2,5): 0,
    (3,1): 0, (3,2): 0, (3,3): 1, (3,4): 0, (3,5): 0,
    (4,1): 0, (4,2): 0, (4,3): 0, (4,4): 1, (4,5): 0,
    (5,1): 0, (5,2): 0, (5,3): 0, (5,4): 0, (5,5): 1
} 

sigmabt = {
    (1,70): 1, (2,70): 0, (3,70): 0, (4,70): 0, (5,70): 0,
    (1,80): 0, (2,80): 1, (3,80): 0, (4,80): 0, (5,80): 0,
    (1,90): 0, (2,90): 0, (3,90): 1, (4,90): 0, (5,90): 0,
    (1,100): 0, (2,100): 0, (3,100): 0, (4,100): 1, (5,100): 0,
    (1,110): 0, (2,110): 0, (3,110): 0, (4,110): 0, (5,110): 1
} 

tau_s = LpVariable.dicts("tau_s", I, lowBound=1, upBound=T_max)
tau_e = LpVariable.dicts("tau_e", I, lowBound=1, upBound=T_max)
sii = LpVariable.dicts("s", (I, I), 0, 1, cat="Binary")
xit = LpVariable.dicts("x", (I, T), 0, 1, cat="Binary") 
Akt = LpVariable.dicts("Akt", (K, T), lowBound=0) 
At = LpVariable.dicts("At", T, lowBound=0) 
Dkt = LpVariable.dicts("D", (K, T), lowBound=0) 
Dt = LpVariable.dicts("Dt", T, lowBound=0) 
ykt = LpVariable.dicts("y", (K, T), 0, 1, cat="Binary") 
fu = LpVariable.dicts("fu", (B, T), lowBound=0)
fiz = LpVariable.dicts("fiz", (B, T), lowBound=0)
T_dep = list(T_OTP.values())
diz = LpVariable.dicts("diz", (B, T_dep), lowBound=0) 
nod = LpVariable.dicts("nod", O, lowBound=0, cat="Integer")
stbt = LpVariable.dicts("stbt", (B, T), lowBound=0, cat="Integer") 

for t in T: 
    model += At[t] == lpSum(N_IT[b, i] for i in I if T_PRI[i] <= t for b in B) 

for t in T: 
    model += Dt[t] == lpSum(diz[b][tprim] for b in B for tprim in T_dep if tprim <= t) 

model += lpSum(At[t] - Dt[t] for t in T) 

for i in I: 
    model += tau_s[i] >= T_PRI[i] + H_TP 

for i in I: 
    for j in I: 
        if i < j: 
            model += tau_s[j] >= tau_e[i] + H_RI - M*(1 - sii[i][j]) 
            model += tau_s[i] >= tau_e[j] + H_RI - M*(sii[i][j]) 

for k in K: 
    for t in T[1:]: 
        model += Akt[k][t] >= Akt[k][t-1] + lpSum(xit[i][t] * lpSum(N_IT.get((b,i),0) * PHIkb[(k,b)] for b in B) for i in I) 

for i in I: 
    model += tau_e[i] == lpSum(xit[i][t]*t for t in T) 

for i in I: 
    model += lpSum(xit[i][t] for t in T) <= 1 

for i in I: 
    for t in T: 
        if t < T_PRI[i]+H_TP: 
            model += xit[i][t] == 0 

for k in K: 
    for t in T: 
        model += Akt[k][t] - Dkt[k][t] <= N_RK_max 

for t in T: 
    if t + H_FJ - 1 <= T_max: 
        model += lpSum(ykt[k][tprim] for k in K for tprim in range(t, t+H_FJ)) <= 1 

for k in K: 
    for t in T: 
        if t >= 1: 
            model += Dkt[k][t]-Dkt[k][t-1] <= N_RK_max*ykt[k][t] 

for k in K: 
    for t in T: 
        if t >= 1: 
            model += Dkt[k][t]-Dkt[k][t-1] <= Akt[k][t-1]-Dkt[k][t-1] 

for b in B: 
    for t in T: 
        if t >= 1: 
            model += fu[b][t] == lpSum((Dkt[k][t]-Dkt[k][t-1])*PHIkb.get((k,b), 0) for k in K) 

for b in B: 
    for t in T[H_FJ:]: 
        model += fiz[b][t] <= fu[b][t - H_FJ] 

for b in B: 
    for t in T_dep: 
        if t-1 in T:
            model += fiz[b][t]+stbt[b][t-1] == stbt[b][t]+diz[b][t] 

for b in B: 
    for t in T[H_FJ:]: 
        model += fiz[b][t] <= N_RK_max * lpSum(ykt[k][t-H_FJ]*PHIkb.get((k,b), 0) for k in K) 

for b in B: 
    for t in T_dep: 
        model += diz[b][t] <= N_OV_max * sigmabt[(b,t)] 

for o in O: 
    t = T_OTP[o] 
    model += nod[o] == lpSum(diz[b][t] for b in B) 

for o in O: 
    model += nod[o] <= N_OV_max 

for o in O: 
    model += nod[o] >= N_OV_max*mi 

model.writeLP("modelMultiLevelClass.txt") 

if os.path.exists(solver_path): 
    if "cplex" in solver_path: 
        model.solve(apis.CPLEX_CMD(path=solver_path, timeLimit=timeLimit)) 
    elif "scip" in solver_path: 
        model.solve(apis.SCIP_CMD(path=solver_path, timeLimit=timeLimit)) 
        print("Solving with: cbc")
else: 
    model.solve(pulp.PULP_CBC_CMD(timeLimit=timeLimit))
    print("Solving with: cbc") 

for v in model.variables(): 
    if v.varValue is not None and abs(v.varValue) > 1e-6: 
        print(v.name, "=", v.varValue) 

print("=== MODEL INFO ===") 
print("Broj promenljivih:", len(model.variables())) 
print("Broj ograničenja:", len(model.constraints))