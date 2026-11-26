"""Recover the effective zero height from a finite tail of full one-loop moments."""
from collections import defaultdict
from fractions import Fraction as F
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx, fmpq, fmpq_mat, fmpq_poly
from sympy import Poly, Rational, factor_list, symbols

BOSON = run_path(str(Path(__file__).with_name("ch39-heavy-polarization-background.py")))
MODEL = BOSON["MODEL"]
q, x = symbols("q x")


def charged_scalar_factors(current=(5, 3, -8)):
    """Factor the exact signed characteristic product in units of lambda."""
    tables = [MODEL["mass_blocks"](r) for r in (0, 1, -1)]
    cache, weights = {}, defaultdict(F)
    for label, (base, metric) in tables[0].items():
        weight = sum(w*c for w, c in zip(MODEL["weight_indices"](label[1]), current))
        if not weight:
            continue
        pos, neg = [table[label][0] for table in tables[1:]]
        key = tuple(m.tobytes() for m in (base, pos, neg, metric))
        if key not in cache:
            size = len(base)
            inv = fmpq_mat(metric.tolist()).inv()
            matrices = [inv*fmpq_mat(m.tolist()) for m in (base, (pos-neg)//2, (pos+neg)//2-base)]
            def char(r):
                return (matrices[0]+r*matrices[1]+r*r*matrices[2]).charpoly()
            coeffs = [Poly(0, q, domain="QQ") for _ in range(size+1)]
            basis = Poly(1, q, domain="QQ")
            for r in range(size+1):
                poly = char(r)
                assert poly == char(-r)
                for j in range(size+1):
                    value = poly[j]
                    exact = Rational(int(value.p), int(value.q))
                    coeffs[j] += basis*((exact-coeffs[j].eval(r*r))/basis.eval(r*r))
                basis *= Poly(q-r*r, q, domain="QQ")
            characteristic = sum(p.as_expr()*x**j for j, p in enumerate(coeffs))
            cache[key] = factor_list(characteristic)[1]
        for factor, multiplicity in cache[key]:
            weights[factor] += weight*multiplicity
    # Discard massless modes, then remove the old weak triplet and color octet.
    weights.pop(x)
    weights[x-q] -= 2*current[1]
    weights[x-q-1] -= 3*current[2]
    return {p: w for p, w in weights.items() if w}


def flint_poly(p):
    return fmpq_poly([fmpq(int(v.p), int(v.q)) for v in reversed(p.all_coeffs())])


def inverse_sums(poly, count):
    """Newton power sums via the series of -p'/p at the origin."""
    deriv = -poly.derivative()
    out = []
    for j in range(count):
        val = deriv[j]-sum(poly[k]*out[j-k] for k in range(1, min(j, poly.degree())+1))
        out.append(val/poly[0])
    return out


def rational_example(factors, start, count):
    z, w, ell = fmpq(1), fmpq(64,25), fmpq(323,2)
    real, lam, sigma = fmpq(9), fmpq(1,12), fmpq(1,256)
    assert 2*real**2 == ell+fmpq(1,2)
    ratio = (3*z+w)/(fmpq(13,10)*z*ell)
    d2, h2 = fmpq(13,10)*ratio*z, fmpq(7,2)*ratio*z
    masses = [h2*(ell+1+2*real), h2*ell, d2*(ell+1-2*real), d2*ell]
    assert masses[-1] == 3*z+w
    last = start+count-1
    scalar = [fmpq(0)]*last
    for p, weight in factors.items():
        pol = Poly(p.subs(q, Rational(int(ratio.p), int(ratio.q))), x, domain="QQ")
        vals = inverse_sums(flint_poly(pol), last)
        for j, val in enumerate(vals, 1):
            scalar[j-1] += fmpq(weight.numerator, weight.denominator)*val/lam**j
    trace, det = (3+ratio)*z+w, (2+ratio)*z*w
    vector = inverse_sums(fmpq_poly([det, -trace, 1]), last)
    tail = []
    for j in range(start, last+1):
        scalar_sum = (scalar[j-1]+6/(lam*sigma*ratio)**j
                      -24/(lam*(1+sigma*ratio))**j-6)
        vector_sum = 6*(vector[j-1]-(3*z+w)**(-j))
        fermion = sum(wt/m**j for wt, m in zip((1,-1,1,-1), masses))
        # (j+1) times the original moment.
        tail.append(scalar_sum/24 + (-2*(j+1)-fmpq(5,8))*vector_sum+(j+1)*fermion)
    return tail, ell, masses


def recover_double_poles(sequence, start):
    capacity = len(sequence)//2
    hankel = fmpq_mat([[sequence[i+j] for j in range(capacity)] for i in range(capacity)])
    rank = hankel.rank()
    block = fmpq_mat([[sequence[i+j] for j in range(rank)] for i in range(rank)])
    coeffs = block.solve(fmpq_mat([[-sequence[i+rank]] for i in range(rank)])).entries()
    den = fmpq_poly([1, *reversed(coeffs)])
    assert all(sum(den[k]*sequence[j-k] for k in range(rank+1)) == 0
               for j in range(rank, len(sequence)))
    num = (den*fmpq_poly(sequence[:rank])).truncate(rank)
    result = []
    coordinate = fmpq_poly([0, 1])
    for factor, multiplicity in den.factor()[1]:
        assert multiplicity in (1, 2)
        if multiplicity == 2:
            # 2*y^(start-2)*num(y)/den''(y) is the coefficient of j+1.
            numerator = 2*coordinate**max(start-2, 0)*num
            derivative = den.derivative().derivative()*coordinate**max(2-start, 0)
            gcd, inv, _ = derivative.xgcd(factor)
            assert gcd == 1
            residue = (numerator*inv)%factor
            assert residue.degree() <= 0 and residue[0].q == 1
            slope = int(residue[0])
            if slope % 12:
                assert factor.degree() == 1 and slope % 12 in (1, 11)
                result.append((-factor[0]/factor[1], slope))
    return rank, sorted(result)


def verify():
    factors = charged_scalar_factors()
    assert sum(Poly(p, x).degree() for p in factors) == 9
    capacity = 9+3+2*(3+4)
    assert capacity == 26
    start = 5
    tail, ell, true_masses = rational_example(factors, start, 2*capacity)
    rank, labeled = recover_double_poles(tail, start)
    (E, ep), (D, dp), (P, pp), (W, wp) = labeled
    assert [ep%12, dp%12, pp%12, wp%12] == [1, 11, 11, 1]
    assert [W, P, E, D] == true_masses and dp == 11
    recovered = 2/(W/P+E/D-2)
    assert recovered == ell
    print("Rational collision case: Hankel rank", rank, "slopes", [ep, dp, pp, wp],
          "recovered ell", recovered)
    ctx.prec = 160
    for n in range(11485, 11491):
        gamma = acb.zeta_zero(n).imag
        ell = ((gamma/624)**2+arb(1)/4).sqrt()
        real = ((ell+arb(1)/2)/2).sqrt()
        imag = ((ell-arb(1)/2)/2).sqrt()
        D = acb(real, imag)
        amplitude = (arb(35)/13).sqrt()
        W, P, E = -amplitude*(D+1), -amplitude*D, D-1
        W2, P2, E2, D2 = [abs(m)**2 for m in (W, P, E, D)]
        assert E2 < D2 < P2 < W2
        recovered = 2/(W2/P2+E2/D2-2)
        height = 624*(recovered**2-arb(1)/4).sqrt()
        assert height.overlaps(gamma)
        print(n, "height from mass ratios", height.str(13))


if __name__ == "__main__":
    verify()
