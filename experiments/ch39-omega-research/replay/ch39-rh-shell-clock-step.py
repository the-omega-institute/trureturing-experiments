"""A pointwise obstruction at tick one and the exact first feasible zero prefix at tick two."""
from pathlib import Path
from runpy import run_path

from flint import arb, acb, arb_mat, ctx

readings = run_path(str(Path(__file__).with_name("ch39-finite-nbb-coherent-four-windows.py")))["readings"]
support = (1, 2, 9, 10, 21, 22, 23, 25)


def columns(g, step):
    log2, log3 = arb(2).log(), arb(3).log()
    return [g*g/(g*g+j*j) for j in range(4)]+[
        (step*s*g).cos() for s in (log2, log3, log3-log2, log3+log2)]


def joint_data():
    levels = readings(5)
    target = levels[:4]+[arb(0)]*4
    heights = [r.imag for r in acb.zeta_zeros(1, 25)]
    cols = [columns(g, 2) for g in heights]
    weights = arb_mat([cols[i-1] for i in support]).transpose().solve(
        arb_mat([[b] for b in target])).entries()
    return levels, target, heights, cols, weights


def main():
    ctx.prec, ctx.cap = 480, 3
    levels, target, heights, cols, weights = joint_data()
    first = heights[0]
    # D=(1-cos(a))(1-cos(b)). For g>=128, nonnegativity alone suffices.
    scale = arb(15)/4*first**2
    bound = arb(".057")
    assert scale/(1+128**2) < bound
    inspected = [v.imag for v in acb.zeta_zeros(1, 60)]
    assert inspected[-1]>128
    for g in inspected:
        deficit = (1-(g*arb(2).log()).cos())*(1-(g*arb(3).log()).cos())
        assert scale/(1+g*g)-deficit < bound
    excess = scale*(1-levels[1])-1-bound
    assert excess>arb(".1437")
    assert excess/3>arb(".0479") and excess/scale>arb(".000191")
    assert min(weights)>arb(".0178") and sum(weights).contains(1)
    # A separating affine functional on the first 24 columns; its value
    # is exactly zero on the seven defining columns and exactly one on target.
    defining = support[:-1]
    dual = arb_mat([cols[i-1] for i in defining]+[target]).solve(
        arb_mat([[0]]*7+[[1]])).entries()
    for i,col in enumerate(cols[:24],1):
        if i not in defining:
            assert sum(v*x for v,x in zip(dual,col))<0
    # The constant observable is one. Shifting down by a quarter gives a
    # strict separating gap even with dyadic rounding of the coefficients.
    rounded = [(v*2**24).floor().unique_fmpz() for v in dual]
    assert all(v is not None for v in rounded)
    sep = [arb(int(v))/2**24 for v in rounded]
    sep[0] -= arb(1)/4
    evaluate = lambda col: sum(v*x for v,x in zip(sep,col))
    assert all(evaluate(col)<arb("-.24") for col in cols[:24])
    assert evaluate(target)>arb(".74")
    # Fifth squared overlap is independent of the eight calibration equations.
    pred = sum(w*heights[i-1]**2/(16+heights[i-1]**2) for i,w in zip(support,weights))
    gap = levels[4]-pred
    assert not gap.contains(0)
    print("Tick one: at least one shell correlation exceeds .0479 in magnitude at the RH first reading.")
    print("Tick two: first feasible actual-zero prefix is 25; support", support)
    print("Probabilities", [w.str(12) for w in weights])
    print("Full RH minus joint model fifth squared overlap",gap.str(15))


if __name__ == "__main__":
    main()
