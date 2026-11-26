"""First zeta zero with a family-adjoint light doublet and the one-loop threshold window."""
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, arb_mat

fixed = run_path(str(Path(__file__).with_name("ch39-shell-gw-fixed-scale-density.py")))
full, flow, shell = (fixed[name] for name in ("full", "flow", "shell"))
count, rows = full["count"], full["rows"]
q = arb


def gauge_target():
    b = flow["beta"]
    inverse = [4*q.pi()/g for g in flow["start"]]
    coefficients = arb_mat([[b[0]-b[1], -q(1)/15], [b[1]-b[2], q(1)/6]])
    rhs = arb_mat([[inverse[0]-inverse[1]], [inverse[1]-inverse[2]]])
    t, required = coefficients.solve(2*q.pi()*rhs).entries()
    return 1000*t.exp(), required


def first_inversion(single, basis):
    angles = [(1-q(d)/26).acos() for d in (3,5,7)]
    assert count(rows(single), q(".01375939")) == 0
    assert count(rows(single), q(".01375942")) == 1
    for n, zero in enumerate(acb.zeta_zeros(1,114), 1):
        p,r = [zero.imag*(a-angles[0]) for a in angles[1:]]
        family = fixed["family"](single,basis,p,r)
        if n < 114:
            assert count(rows(family), q(".01375942")) == 0
    for cut,expected in [(".0124598",0),(".0124599",1),(".0126692",1),
                         (".0126694",2),(".019",2)]:
        assert count(rows(family), q(cut)) == expected
    print("First adjoint ground at positive zero", n, "height", zero.imag.str(18))
    return p,r,family


def chiral_window(p,r,width):
    selected = shell["yuk"]["model"]
    unit = shell["yuk"]["identity"](36)
    for sign in (1,-1,0):
        base,center = selected["action_kernel"](sign,p,r)
        tube = selected["action_kernel"](sign,p+q(0,width),r+q(0,width))[1]
        gap = q(".12")
        assert count(center,-gap) == count(center,gap) == 18-3*sign
        for block in selected["character_blocks"](base,center):
            assert count(block,gap) == 6-sign
        delta = sum(abs(z).abs_upper()**2 for z in (acb_mat(tube)-acb_mat(center)).entries()).sqrt()
        assert delta < gap
        scale = q(8)
        assert max(sum(abs(z) for z in row).upper() for row in center) < scale
        x = acb_mat(center)/scale
        b = gap/scale
        while 1-b > q("1e-14"):
            x = x*(3*unit-x*x)/2
            b = b*(3-b*b)/2
        start,orientation = (18,1) if sign >= 0 else (0,-1)
        comp = [[(unit[i,j]+orientation*x[i,j])/2 for j in range(start,start+18)]
                for i in range(start,start+18)]
        margin = q(".004")
        assert count(comp,margin) == 0
        # The vertical-line resolvent formula for the sign projection gives this bound.
        assert margin > (1-b)/2 + delta/(2*(gap-delta))


def threshold_window(single,family,width):
    Q,E = gauge_target()
    center = full["all_levels"](single,family)
    # Disjoint shortest edges, with unitary link partial derivatives of norm <= 1.
    # The adjoint change is <= 4w; both ends of a spectral gap can move.
    gaps = [d+q(0,(8*width).upper()) for d in center]
    assert min(gaps)>0
    def bounds(levels, lo, hi):
        def masses(M):
            return [(125**2+M**2*d).sqrt() for d in levels]
        lo,hi = q(lo),q(hi)
        assert sum((Q/m).log() for m in masses(lo)) > E
        assert sum((Q/m).log() for m in masses(hi)) < E
        first,last = masses(lo)[0],masses(hi)[-1]
        assert first > 1000 and last < Q
        print("M bracket",float(lo),float(hi),"first threshold above",float(first.lower()),
              "last below fraction of Q",float((last/Q).upper()))
        return first,last
    bounds(center, "1.4627e12", "1.4629e12")
    first,last = bounds(gaps,"1.43e12","1.50e12")
    assert first > q("1.62e10") and last < q(".075")*Q


def main():
    single,_,_,basis = shell["transported_spectra"]()
    p,r,family = first_inversion(single,basis)
    width = q("1e-5")
    assert q(".0126694")+4*width < q(".01375939")
    assert q(".01375942") < q(".019")-4*width
    assert q(".0126692")-q(".0124599") > 8*width
    chiral_window(p,r,width)
    threshold_window(single,family,width)
    assert (width/q.pi())**2 > q("1.0132e-11")
    print("Window contains two family levels below singlet and unpaired indices (+3,-3,0)")


if __name__ == "__main__":
    main()
