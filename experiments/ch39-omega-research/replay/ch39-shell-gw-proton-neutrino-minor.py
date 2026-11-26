"""Inverse tau-tau texture from the minimal renormalizable baryon-triality seesaw."""
from flint import arb, arb_poly, ctx

ctx.prec = 160
q = arb
solar = q("7.49e-5")


def analyse(name, angles, delta_deg, atm, brackets):
    s12, s23, s13 = map(q, angles)
    c12, c23, c13 = 1-s12, 1-s23, 1-s13
    re = (q(delta_deg)*arb.pi()/180).cos()
    first = s12*s23+c12*c23*s13-2*(s12*c12*s23*c23*s13).sqrt()*re
    third = c23*c13
    weights = [first, 1-first-third, third]
    atm = q(atm)
    offsets = [0, solar, atm] if atm > 0 else [-atm-solar, -atm, 0]
    factors = [arb_poly([d, 1]) for d in offsets]
    terms = [weight**2*factors[(i+1)%3]*factors[(i+2)%3] for i,weight in enumerate(weights)]
    poly = sum(2*term*terms[(i+1)%3]-term*term for i,term in enumerate(terms))
    signs = [1 if c > 0 else -1 if c < 0 else 0 for c in poly.coeffs()]
    assert 0 not in signs and signs[0] == -1
    assert sum(a!=b for a,b in zip(signs, signs[1:])) == len(brackets)
    for l,r in brackets:
        left, right = q(l), q(r)
        assert poly(left**2)*poly(right**2) < 0
        span = left.union(right)
        total = sum((span**2+f).sqrt() for f in offsets)
        print(name, "mass endpoint", l, r, "sum near", total.str(8))
    # D=d*I and d²=Q*mmax place all heavy states at or above Q.
    low = q(".010" if atm > 0 else ".070")
    assert poly(low**2) > 0
    masses = [(low**2+f).sqrt() for f in offsets]
    d = (q("4.352e13").union(q("4.473e13"))*masses[2 if atm>0 else 1]*q("1e-9")).sqrt()
    print("sample common Dirac mass GeV", d.str(7), "heavy mass span ratio",
          (max(masses)/min(masses)).str(7))


if __name__ == "__main__":
    # NuFIT 6.0 Table 1, first and second atmospheric data selections.
    analyse("NO first", (".307",".561",".02195"),177,".002534",
            [(".00639428",".00639429")])
    analyse("IO first", (".308",".562",".02224"),285,"-.002510",
            [(".05640080",".05640082")])
    analyse("NO second", (".308",".470",".02215"),212,".002513",
            [(".00582797",".00582799"),(".12513925",".12513927")])
    analyse("IO second", (".308",".550",".02231"),274,"-.002484",
            [(".06298374",".06298376")])
    # Envelope of both three-sigma IO boxes without assumed correlations.
    atm_min = q(".002421")
    d21_max = q(".0000805")
    u = (1-q(".597"))*(1-q(".02409"))
    floor = u*((atm_min-d21_max)/(1-2*u)).sqrt()
    total = floor+(floor**2+atm_min-d21_max).sqrt()+(floor**2+atm_min).sqrt()
    print("IO union of three-sigma boxes:", floor.lower().str(8), total.lower().str(8))
    assert floor>q(".041") and total>q(".168")
