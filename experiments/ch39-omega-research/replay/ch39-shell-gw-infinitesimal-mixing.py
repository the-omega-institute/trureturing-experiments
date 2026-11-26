"""Compute the sharp light-mass cost of an infinitesimal heavy-line rotation."""
from pathlib import Path
from runpy import run_path

import numpy as np
from flint import acb_mat, arb

rank = run_path(str(Path(__file__).with_name("ch39-shell-gw-rank-one-mixing.py")))
model, seven = rank["model"], rank["seven"]
adj, hs = model["adj"], model["hs"]


def tangent_gain(modes, x, v):
    u = x.conjugate()
    pu, pv = model["eye"](3) - u*adj(u), model["eye"](3) - v*adj(v)
    blocks = [pu*b*pv for b in modes]
    c = [(adj(u)*b*v)[0,0] for b in modes]
    metric = acb_mat([[hs(e,f)+z.conjugate()*w for f,w in zip(blocks,c)]
                      for e,z in zip(blocks,c)])
    assert seven["count"](seven["rows"](metric),arb(0)) == 0
    vectors = [pu*b*v for b in modes]
    W = acb_mat([[w[j,0] for w in vectors] for j in range(3)])
    gain = W*metric.solve(adj(W))
    proposal = float(np.linalg.eigvalsh(model["array"](gain))[-1])
    lower, upper = arb(proposal)-arb("2e-6"),arb(proposal)+arb("2e-6")
    assert seven["count"](seven["rows"](gain),lower)<3
    assert seven["count"](seven["rows"](gain),upper)==3
    return lower.sqrt().union(upper.sqrt())


def check(index):
    modes, _ = seven["spectral_band"](index)
    left, right = rank["mix"]["orbit_frame"](modes)
    modes = [adj(left)*b*right for b in modes[:5]]
    outside = model["complement"](modes)
    gains = []
    for p in range(3):
        for seed in rank["numeric_seeds"](outside,p):
            x,v = rank["certified_root"](outside,p,seed)
            gains.append(tangent_gain(modes,x,v))
    low = max(g.lower() for g in gains)
    high = max(g.upper() for g in gains)
    coefficient = 1/low.union(high)
    print("zero",index,"one-sided gain by root",[v.str(9) for v in gains],
          "limiting light-mass sum / mixing distance",coefficient.str(10))
    lo, hi = (".3252454", ".3252456") if index==43 else (".2853296", ".2853298")
    assert arb(lo) < coefficient < arb(hi)
    return coefficient


if __name__ == "__main__":
    for index in (43, 195):
        check(index)
