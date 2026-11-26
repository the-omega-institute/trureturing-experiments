"""Sharp relative error for fitting three logarithmic polarizations with a gapped positive spectrum."""
from fractions import Fraction
from itertools import combinations
from pathlib import Path
from runpy import run_path

from flint import arb, arb_mat, ctx, fmpq, fmpq_poly

sharp = run_path(str(Path(__file__).with_name("ch39-sharp-hidden-spectral-gap.py")))
ball = sharp["ball"]


def system(xs):
    factors = [fmpq_poly([1,fmpq(a.numerator,a.denominator)]) for a in (x-1 for x in xs)]
    denominator = factors[0]*factors[1]*factors[2]
    rows = [denominator//f for f in factors]
    return arb_mat([[row[j] for j in range(3)] for row in rows])


def relative_boundary(xs, limit, divided_logs=None):
    basis = system(xs)
    truth = arb_mat([[(divided_logs[x] if divided_logs is not None else sharp["rational_log"](x)/ball(x-1))] for x in xs])
    covector = basis.transpose().solve(arb_mat([[0],[-limit],[1]]))
    contributions = [a*b for a,b in zip(covector.entries(),truth.entries(),strict=True)]
    assert all(not u.contains(0) for u in contributions)
    epsilon = sum(contributions)/sum(abs(v) for v in contributions)
    if epsilon < 0:
        # The unrestricted moments satisfy m2/m1 < limit. The positive
        # small-left-node Radau construction gives an exact finite witness.
        m0,m1,m2 = basis.solve(truth).entries()
        assert m2/m1 < limit and m1*m1 < m0*m2
        return arb(0),None,None,truth,basis
    assert 0 < epsilon < 1
    corner = arb_mat([[v*(1-epsilon if u>0 else 1+epsilon)]
                      for v,u in zip(truth.entries(),contributions,strict=True)])
    m0,m1,m2 = basis.solve(corner).entries()
    at_limit = m1/limit
    at_zero = m0-at_limit
    assert (m2-limit*m1).contains(0)
    if at_zero > 0 and at_limit > 0:
        return epsilon, at_zero, at_limit, truth, basis
    # The nearest admissible point is a single finite mass when even the
    # endpoint atom cannot match the ratio of the active high and low readings.
    ratios = [1/((1+ball(x-1)*limit)*v) for x,v in zip(xs,truth.entries())]
    hi,lo = max(range(3),key=lambda i:float(ratios[i])),min(range(3),key=lambda i:float(ratios[i]))
    assert hi > lo and all(ratios[hi]>r for i,r in enumerate(ratios) if i!=hi)
    assert all(ratios[lo]<r for i,r in enumerate(ratios) if i!=lo)
    epsilon = (ratios[hi]-ratios[lo])/(ratios[hi]+ratios[lo])
    raw_weight = 2/(ratios[hi]+ratios[lo])
    divided_weight = raw_weight
    for x in xs:
        divided_weight /= 1+ball(x-1)*limit
    return epsilon, arb(0), divided_weight, truth, basis


def finite_error(xs, limit, small, atom0, atom1, truth, basis):
    # Move the positive boundary weight from zero to a finite mass while
    # preserving the upper support limit. Convert divided moments back to data.
    model = basis*arb_mat([[atom0*small**i+atom1*limit**i] for i in range(3)])
    return max(abs(model[i,0]/truth[i,0]-1).upper() for i in range(3))


def verify():
    with ctx.workprec(320):
        for xs in ((Fraction(5,3),Fraction(7,3),Fraction(3)),
                   (Fraction(3,2),Fraction(2),Fraction(5,2)),
                   (Fraction(9,8),Fraction(5,4),Fraction(11,8))):
            for cut in (Fraction(1),Fraction(3,2)):
                limit = ball(1/(1+cut**2))
                eps,w0,w1,truth,basis = relative_boundary(xs,limit)
                errors = [finite_error(xs,limit,ball(Fraction(1,2**n)),w0,w1,truth,basis) for n in (12,20,30)]
                if w0 > 0:
                    assert errors[0]>errors[1]>errors[2]>eps
                else:
                    assert all(abs(v-eps) < arb("1e-30") for v in errors)
                assert errors[-1]-eps < arb("1e-8")
                print(tuple(map(str,xs)),"gap / M_N",str(cut),"finite endpoint attained",w0==0,"infimum relative error",eps.str(17),
                      "finite witness at 2^-30",errors[-1].str(16),flush=True)

        grid = [1+Fraction(i,8) for i in range(1,17)]
        divided_logs = {x:sharp["rational_log"](x)/ball(x-1) for x in grid}
        for cut in (Fraction(1),Fraction(3,2)):
            limit = ball(1/(1+cut**2))
            candidates = [(relative_boundary(xs,limit,divided_logs)[0],xs)
                          for xs in combinations(grid,3)]
            candidates.sort(key=lambda v:float(v[0].mid()), reverse=True)
            best,xs = candidates[0]
            expected = {(Fraction(9,8),Fraction(7,4),Fraction(3))} if cut==1 else {
                (grid[0],middle,grid[-1]) for middle in grid[1:-1]}
            winners = {ys for r,ys in candidates if r.overlaps(best)}
            assert winners == expected and len(candidates)==560
            assert all(r < best for r,ys in candidates if ys not in winners)
            if cut > 1:
                for ys in winners:
                    _,zero,_,truth,_ = relative_boundary(ys,limit,divided_logs)
                    ratios = [1/((1+ball(x-1)*limit)*v) for x,v in zip(ys,truth.entries())]
                    assert zero == 0 and ratios[0]<ratios[1]<ratios[2]
            print("grid gap",str(cut),"maximum",best.str(16),"first triple",tuple(map(str,xs)),
                  "number of maximizing triples",len(winners))



if __name__ == "__main__":
    verify()
