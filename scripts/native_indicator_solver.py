"""Gurobi native-indicator implementation for the bounded diagnostic only."""

from types import SimpleNamespace

import numpy as np


def solve_native(s, lo, hi, biomass_index, growth, witness, output):
    import gurobipy as gp

    n = len(lo)
    with gp.Env() as env, gp.Model("restricted_native_indicators", env=env) as m:
        m.Params.TimeLimit = 120
        m.Params.MIPGap = 1e-4
        m.Params.Threads = 4
        m.Params.FeasibilityTol = 1e-9
        m.Params.IntFeasTol = 1e-9
        m.Params.IntegralityFocus = 1
        m.Params.LogFile = str(output / "native.gurobi.log")
        v = m.addMVar(n, lb=np.minimum(lo, 0), ub=np.maximum(hi, 0), name="flux")
        y = m.addMVar(n, vtype=gp.GRB.BINARY, name="selected")
        m.addMConstr(s, v, "=", np.zeros(s.shape[0]))
        vv, yy = v.tolist(), y.tolist()
        for i in range(n):
            m.addGenConstrIndicator(yy[i], 0, vv[i] == 0)
            m.addGenConstrIndicator(yy[i], 1, vv[i] >= float(lo[i]))
            m.addGenConstrIndicator(yy[i], 1, vv[i] <= float(hi[i]))
        m.addConstr(vv[biomass_index] >= growth)
        m.setObjective(y.sum(), gp.GRB.MINIMIZE)
        v.Start = witness[:n]
        y.Start = witness[n:]
        m.optimize()
        has = m.SolCount > 0
        if has:
            m.write(str(output / "native.sol"))
        return SimpleNamespace(
            status=0 if m.Status == gp.GRB.OPTIMAL else 1,
            message=f"Gurobi status {m.Status}",
            x=np.r_[v.X, y.X] if has else None,
            fun=m.ObjVal if has else None,
            mip_gap=m.MIPGap if has else None,
            mip_dual_bound=m.ObjBound,
            mip_node_count=m.NodeCount,
        )
