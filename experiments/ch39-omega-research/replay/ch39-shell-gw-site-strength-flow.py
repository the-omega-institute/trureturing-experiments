"""Balanced-site lower bound on the top singular value after first-threshold matching."""
from pathlib import Path
from runpy import run_path

from flint import acb_mat, arb

flow = run_path(str(Path(__file__).with_name("ch39-shell-gw-first-heavy-flow.py")))
fit = flow["fit"]
adj, rows, count = (fit[key] for key in ("adj", "rows", "count"))


def below(M, ceiling):
    return count(rows(M), ceiling) == M.nrows()


def trace(M):
    return sum(M[i,i].real for i in range(M.nrows()))


def site_gain(index,F,D):
    K = F*adj(D)*(D*adj(D)).inv()
    assert below(K*adj(K), arb(1))
    E = F-K*D
    A, B = acb_mat(3,3), acb_mat(3,3)
    pad = arb(".02" if index == 43 else ".035")
    for j in range(6):
        site = acb_mat([[E[r,9*j+k] for k in range(9)] for r in range(9)])
        H = site*adj(site)
        norm = trace(H).sqrt()
        S = acb_mat([[sum(H[3*r+k,3*s+k] for k in range(3))/norm for s in range(3)]
                     for r in range(3)])
        T = acb_mat([[sum(H[3*k+r,3*k+s] for k in range(3))/norm for s in range(3)]
                     for r in range(3)])
        S += pad*trace(S)/3*flow["eye"]
        T += pad*trace(T)/3*flow["eye"]
        tensor = acb_mat([[S[r//3,s//3]*T[r%3,s%3] for s in range(9)] for r in range(9)])
        assert below(H-tensor, arb(0))
        A += S
        B += T
    gain = arb(".0870" if index == 43 else ".0883")
    assert all(below(m,gain*arb(6).sqrt()) for m in (A,B))
    return gain


def check(index):
    _, up, down, _, F, D, _, _, _ = flow["at_threshold"](index)
    gain = site_gain(index,F,D)
    # e3 on each side of the running up matrix is a certified top singular-value lower bound.
    top_lo = abs(up[2,2])
    down_hi = fit["hs"](down,down).real.sqrt()
    required = (top_lo-down_hi)/gain
    frozen_limit = arb("5.64")
    assert not below(F*adj(F), arb(".5134")**2)
    assert (8*arb.pi()/3).sqrt()/arb(".5134") < frozen_limit
    assert required > arb("6.12" if index == 43 else "5.94")
    wave = 3*(required*arb(".5134"))**2/(16*arb.pi())
    print("zero",index,"balanced strength lower",required.str(8),"frozen full-spectrum wave lower",
          wave.str(8),"necessary peak-to-rms ratio when g <= 5.64",(required/frozen_limit).str(8))


if __name__ == "__main__":
    for index in (43,195):
        check(index)
