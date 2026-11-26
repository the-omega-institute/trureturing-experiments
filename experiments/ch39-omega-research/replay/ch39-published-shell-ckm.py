"""Compare the CKM table, count code edits, and bound the normalized family current."""
from itertools import product
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, ctx

root = Path(__file__).parent
model = run_path(str(root/"ch39-fixed-local-flavor.py"))
mixture = run_path(str(root/"ch39-clock-flavor-mixture.py"))
centers = tuple(map(arb, (".22704", ".04275", ".00377", "1.139")))
widths = tuple(map(arb, (".00082", ".00042", ".00008", ".023")))


def angles(response, codes):
    matrix = model["mixing"](response,[arb(code)/2048 for code in codes])
    s,t,u,phase = mixture["readouts"](mixture["direct_observables"](matrix))
    return s.asin(), t.asin(), u.asin(), phase


def passes(values):
    errors = [abs(v-c)/w for v,c,w in zip(values,centers,widths,strict=True)]
    if all(err < 1 for err in errors):
        return True
    assert any(err > 1 for err in errors)
    return False


def verify():
    ctx.prec = 192
    low, high = [[c+sign*w for c,w in zip(centers,widths,strict=True)] for sign in (-1,1)]
    box = [lo.union(hi) for lo,hi in zip(low,high,strict=True)]
    v = model["shell"]["standard_matrix"](*(a.sin() for a in box[:3]),box[3])
    assert all(abs(v[2][2]) > abs(v[i][j]) for i,j in product(range(3),repeat=2) if (i,j)!=(2,2))
    def floor(theta23, theta13):
        t, u = theta23.sin(), theta13.sin()
        return 3*(t*t+u*u-t*t*u*u).sqrt()/2
    assert floor(low[1],low[2]) > arb(".06371647")
    assert floor(high[1],high[2]) < arb(".06499150")
    original = model["codes"]
    response = model["native"]["plane_response"](acb.zeta_zero(11431).imag)
    assert (low[1].sin()-model["shell"]["centers"][1]-model["shell"]["sine_radius"]) > arb(".00043235")
    assert angles(response,original)[1] < low[1]
    witness = (*original[:3],original[3]-7)
    values = angles(response,witness)
    assert passes(values)
    print("one changed coefficient, four angles in radians",*(v.str(11) for v in values))
    mixing = model["mixing"](response,[arb(n)/2048 for n in witness])
    charge = acb_mat([[int(i==j)-1 for j in range(3)] for i in range(3)])
    down = mixing.conjugate().transpose()*charge*mixing
    magnitudes = [abs(down[i,j]) for i,j in ((0,1),(0,2),(1,2))]
    for value, rounded in zip(magnitudes, (".973921",".896081","1.099213"),strict=True):
        assert abs(value-arb(rounded)) < arb(".0000005")
    # The witness with one changed coefficient also gives a finite upper bound on the edit search.
    winners = []
    for cost in range(1,8):
        for change in product(range(-cost,cost+1), repeat=4):
            if sum(map(abs,change)) != cost:
                continue
            candidate = tuple(a+b for a,b in zip(original,change,strict=True))
            if passes(angles(response,candidate)):
                winners.append(change)
        if winners:
            break
    assert cost == 7 and winners == [(0,0,0,-7)]
    print("minimum total numerator change",cost,"all minimizers",winners)
    labels = []
    for n,z in enumerate(acb.zeta_zeros(10815,1278),10815):
        response = model["native"]["plane_response"](z.imag)
        if passes(angles(response,witness)):
            labels.append(n)
    assert labels == [11431]
    print("labels passing the four published marginal windows",labels)


if __name__ == "__main__":
    verify()
