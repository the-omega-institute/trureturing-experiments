"""Sharp bare mixing floor when each bridge family retains its shell ratio."""
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx, fmpq as Q, fmpq_mat


def check_family_order():
    eps, t = Q(1, 5), Q(7, 11)
    for h in ([Q(1), Q(5), Q(3)], [Q(1), Q(4), Q(3)]):
        products = []
        for cc in (True, False):
            block = fmpq_mat(3, 4)
            for a in range(3):
                block[(-a-int(cc)) % 3, a] = h[a] if cc else eps*h[a]
            block[2, 3] = t
            products.append((block*block.transpose()).det())
        assert products[0] == (h[1]*h[2])**2*(h[0]**2+t*t)
        assert products[1] == (eps**2*h[0]*h[2])**2*((eps*h[1])**2+t*t)
        assert (products[1] == eps**6*products[0]) == (h[0] == eps*h[1])


def verify():
    ctx.prec = 160
    check_family_order()
    ext = run_path(str(Path(__file__).with_name("ch39-selected-extended-flow.py")))
    v, state = ext["initial"]()
    original = state[2].sqrt()
    eps = (state[3]/state[2]).sqrt()
    cap = (4*arb.pi()).sqrt()
    assert 21698 < 1/eps < 21701
    for h in (cap, arb(3)):
        h0 = eps*h
        # Saturate h1=h2=h and match both heavy mass products.
        t = ((original**3/h**2)**2-h0**2).sqrt()
        assert t > eps*h0 and 0 < h0 < h
        assert (h**4*(h0*h0+t*t)-original**6).contains(0)
        assert ((eps**2*h0*h)**2*((eps*h)**2+t*t)-(eps*original)**6).contains(0)
        overlap = h0*h0/(h0*h0+t*t)
        assert (overlap-(eps*h**3/original**3)**2).contains(0)
        assert overlap < arb("2.126e-6")
        assert 2*arb.pi()*overlap/10 < arb("1.336e-6")
        gap = v*eps*h0/(4*(1+(arb(25)/24)**2))
        if h > 3:
            assert arb("6.6389e13") < t*v < arb("6.6390e13")
        else:
            assert arb("9.2697e13") < t*v < arb("9.2698e13")
            assert arb("4.5257e5") < gap < arb("4.5258e5")
        for m in (3, 9, 19):
            theta = 2*arb.pi()/(10*m)
            link = 1-overlap+overlap*acb(0, -theta).exp()
            assert 0 < -m*link.arg() < 2*arb.pi()*overlap/10
        print("Family cap",h.str(8),"mu / GeV",(t*v).str(11),
              "light target weight",overlap.str(10),"core gap / GeV",gap.str(10))


if __name__ == "__main__":
    verify()
