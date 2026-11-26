"""Certified occupation and holonomy-change budgets for a uniform feasible register."""
from flint import acb, arb, ctx


def verify():
    ctx.prec = 144
    first, sure_first, sure_last, last = 10855, 10905, 12054, 12092
    target = 10875 - first
    size_floor = sure_last - sure_first + 1 + 2  # Additional witnesses 10871 and 10875.
    assert size_floor == 1152
    heights = [z.imag for z in acb.zeta_zeros(first, last-first+1)]
    assert all(a < b for a, b in zip(heights, heights[1:]))
    bare = [1/(4*g*g+1) for g in heights]
    oriented = [4*g*p for g, p in zip(heights, bare)]
    confidence = arb(99)/100
    for observable, bound in [(bare, "1.8243e-10"), (oriented, "3.9893e-6")]:
        # The largest possible uniform mean adds every low optional label to the
        # guaranteed block and omits every high optional label.
        stop = sure_last - first + 1
        mean = sum(observable[:stop])/stop
        assert observable[sure_first-first-1] > mean > observable[stop]
        output_floor = confidence*observable[target] + (1-confidence)*observable[-1]
        assert output_floor - mean > arb(bound)
        print("Minimum observable change", (output_floor-mean).lower().str(12))
    angle = confidence.sqrt().asin() - (1/arb(size_floor).sqrt()).asin()
    assert angle > arb("1.44116")
    assert 1/(confidence*size_floor) < arb(".000876824")
    print("Off-diagonal Hamiltonian action >=", angle.lower().str(10))


if __name__ == "__main__":
    verify()
