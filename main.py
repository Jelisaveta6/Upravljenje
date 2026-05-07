from pulp import *
import os
solver_path = "C:\\x64_win64\\cplex.exe"
timeLimit = 300
n = ["1111", "1112", "1113"]
m = ["1000", "1002", "1003"]
ai = {"1111": 10, "1112": 50, "1113": 150}
ti = {"1111": 30, "1112": 30, "1113": 30}
dj = {"1000": 100, "1002": 160, "1003": 250}
hi = {"1111": 5, "1112": 5, "1113": 5}
ui = {"1111": ["A", "B", "C"], "1112": ["D", "E"], "1113": ["F", "G"]}
uj = {"1000": ["C", "D"], "1002": ["A", "B"], "1003": ["E", "F", "G"]}
U = ["A", "B", "C", "D", "E", "F", "G"]
ju = {"1000": ["A","D","E","F","G"], "1002": ["A","B","D","E","F","G"], "1003": ["B","C","D","E","F","G"]}
F = ["r","w","l"]
efu = {("r", "A"): 10,("w","A"): 1250, ("l","A"): 150, ("r", "B"): 20, ("w","B"):2350, ("l","B"):300, ("r", "C"):15, ("w","C"):1750, ("l","C"):225, ("r", "D"): 10, ("w","D"):1250,("l","D"):150, ("r","E"):20,("w","E"):3350, ("l","E"):300, ("r","F"):20, ("w","F"):2250, ("l","F"):300, ("r","G"):20, ("w","G"):2550, ("l","G"):300}
t = {"A": [30], "B":[30], "C":[30], "D":[30], "E":[30], "F":[30], "G":[30]}
bvj = {"1000": 20, "1002": 20, "1003": 20}
gjf = {("1000","r"): 2000, ("1000","w"): 2000, ("1000","l"): 2000, ("1002", "r"): 2000, ("1002", "w"): 2000, ("1002", "l"): 2000, ("1003", "r"): 2000, ("1003", "w"):2000, ("1003", "l"):2000}
pjf = {("1000","r"): 20, ("1000","w"): 20, ("1000","l"): 20, ("1002", "r"): 20, ("1002", "w"): 20, ("1002", "l"): 20, ("1003", "r"): 20, ("1003", "w"):20, ("1003", "l"):20}
cju = {("1000", "A"): 100, ("1000", "B"): 100, ("1000", "C"): 100, ("1000", "D"): 100, ("1000", "E"): 200, ("1000", "F"): 250, ("1000", "G"): 290, ("1002", "A"): 300, ("1002", "B"): 300, ("1002", "C"): 300, ("1002", "D"): 300, ("1002", "E"): 300, ("1002", "F"): 300, ("1002", "G"): 300, ("1003", "A"): 300, ("1003", "B"): 300, ("1003", "C"): 300, ("1003", "D"): 300, ("1003", "E"): 300, ("1003", "F"): 300, ("1003", "G"): 300}
oj = {"1000": 10, "1002": 20, "1003": 30}
K = [1,2,3]
L = [1,2,3]
M = 1000
T = 25
Pik = LpVariable.dicts("Pik", [(i,k) for i in n for k in K], 0, 1, LpBinary)
Qjl = LpVariable.dicts("Qjl", [(j,l) for j in m for l in L], 0,1, LpBinary)
Xju = LpVariable.dicts("Xju", [(j,u) for j in m for u in U], 0, 1, LpBinary)
Yj = LpVariable.dicts("Yj", [(j) for j in m], 0,1, LpBinary)
Dlk = LpVariable.dicts("Dlk", [(l,k) for l in L for k in K], 0,1, LpBinary)
ROjf = LpVariable.dicts("ROjf", [(j,f) for j in m for f in F], 0,1,LpBinary)
Rk = LpVariable.dicts("Rk", K, lowBound=0)
Sl = LpVariable.dicts("Sl", L, lowBound=0)
model = LpProblem("Railcar_yards", LpMinimize)
model += (lpSum(cju[(j, u)]*Xju[(j, u)] for j in m for u in U) - lpSum(oj[(j)]*Yj[(j)] for j in m))
for k in K[1:]:
    model += (Rk[k] >= Rk[k-1] + lpSum(hi[i] * Pik[(i,k)] for i in n))
for k in K:
    model += (Rk[k] - lpSum((ai[(i)]+ti[(i)])*Pik[(i,k)] for i in n) >= 0)
for k in K:
    model += (lpSum(Pik[(i,k)] for i in n) == 1)
for i in n:
    model += (lpSum(Pik[(i,k)] for k in K) == 1)
for k in K:
    for l in L:
        model += (Sl[l] - Rk[k] - lpSum(hi[(i)]*Pik[(i,k)] for i in n) - T >= M*(Dlk[(l,k)]-1))
for k in K:
    for i in n:
        for j in m:
            for l in L:
                model += (M*(Dlk[(l,k)]+2 - Pik[(i,k)] - Qjl[(j,l)]) >= lpSum(Xju[(j,u)] for u in U if u in ui[i] and u in uj[j]))
for l in L[:-1]:
    model += (Sl[l+1]- Sl[l] - lpSum(bvj[j]*Qjl[(j,l)] for j in m)>= 0)
for l in L[:-1]:
    model += (Sl[l] <= lpSum((dj[(j)] - bvj[j] - bvj[j])*Qjl[(j,l)] for j in m) + M*(1-lpSum(Qjl[(j,l)] for j in m)))
for j in m:
    model += (lpSum(Qjl[(j,l)] for l in L) == Yj[j])
for l in L:
    model += (lpSum(Qjl[(j,l)] for j in m) == 1)
for u in U:
    model += (lpSum(Xju[(j,u)] for j in ju) == 1)
for j in m:
    for f in F:
        model += (lpSum(efu[(f,u)]*Xju[(j,u)] for u in U) >= pjf[(j,f)]*ROjf[(j,f)])
for j in m:
    model += (lpSum(ROjf[(j,f)] for f in F) >= Yj[(j)])
if os.path.exists(solver_path):
    if "cplex" in solver_path:
        model.solve(apis.CPLEX_CMD(path=solver_path, timeLimit=timeLimit))
    elif "scip" in solver_path:
        model.solve(apis.SCIP_CMD(path=solver_path, timeLimit=timeLimit))
else:
    model.solve(pulp.PULP_CBC_CMD(timeLimit=timeLimit))
print("Status: ", LpStatus[model.status])
model.writeLP('RailCarYards.txt')
for v in model.variables():
    if v.varValue is not None and abs(v.varValue) > 1e-6:
        print(v.name, "=", v.varValue)
print("===MODEL INFO===")
print("Broj promenljivih:", len(model.variables()))
print("Broj ograničenja:", len(model.constraints))