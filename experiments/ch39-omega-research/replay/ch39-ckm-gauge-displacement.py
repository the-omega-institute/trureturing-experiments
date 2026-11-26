"""Bound the change of the three gauge outputs caused by quark mixing."""
from pathlib import Path
from runpy import run_path

from flint import arb, arb_mat, ctx

volume = run_path(str(Path(__file__).with_name("ch39-volume-contenders.py")))
model, q, loop = volume["gauge"],volume["q"],volume["loop"]


def displacement(stop=33):
    low,high,logs = volume["input_bounds"]()
    boxes = [l.union(h) for l,h in zip(low,high)]
    delta = arb_mat(12,2) # light-light mixing, then unit third-family leakage
    def observe(step,g,y,traces):
        nonlocal delta
        jac = volume["flow_jacobian"](g,y)
        a = arb_mat([[abs(jac[i,j]).upper() for j in range(12)] for i in range(12)])
        norm = max(sum(a[i,j] for j in range(12)).upper() for i in range(12))
        source = arb_mat(12,2)
        for f in (0,1):
            opposite = y[3*(1-f):3*(2-f)]
            for i in range(3):
                row = 3+3*f+i
                if i < 2:
                    light = abs(opposite[1-i]-opposite[i])
                    heavy = abs(opposite[2]-opposite[i])
                    source[row,0] = (3*light/loop).upper()
                else:
                    heavy = max(abs(v-opposite[2]).upper() for v in opposite[:2])
                source[row,1] = (3*heavy/loop).upper()
        derivative = a*delta+source
        growth = step*(step*norm).exp()
        col = [max(derivative[i,j].upper() for i in range(12)) for j in range(2)]
        increment = arb_mat([[growth*c for c in col]]*12)
        # Integrate the nonnegative comparison system on this step. Iteration
        # preserves an upper bound for every component and both source columns.
        for _ in range(3):
            increment = step*(derivative+a*increment)
        end = delta+increment
        delta = arb_mat([[end[i,j].upper() for j in range(2)] for i in range(12)])
    model["evolve"](boxes,logs,q(stop),8192,diagonal_quarks=False,observe_flavors=observe)
    return [[delta[i,j]/r for j in range(2)] for i,r in enumerate(model["ratios"])]


def verify():
    ctx.prec = 128
    result = displacement()
    for i,(light,heavy) in enumerate(result):
        print(i+1,"normalized gauge deviation < light + (1-|V_tb|^2 bound)*heavy",
              light.upper().str(10),heavy.upper().str(10))
    print("max gauge deviation for |V_tb|^2 >= 0.997:",
          max((a+q(".003")*b).upper() for a,b in result).str(10))


if __name__ == "__main__":
    verify()
