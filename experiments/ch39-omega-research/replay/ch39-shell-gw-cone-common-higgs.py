"""Heavy-subspace cost of common up/down Yukawa profiles in the eight-threshold cone."""
from pathlib import Path
from runpy import run_path

from flint import acb, arb

folder = Path(__file__).resolve().parent
rank = run_path(str(folder/"ch39-shell-gw-cone-flavor-rank.py"))
native = run_path(str(folder/"ch39-shell-gw-native-arc-reference.py"))
q = arb


def norm(a):
    return sum(abs(v)**2 for v in a).sqrt()


def apply(row, col):
    return sum(c*v for c,v in zip(row,col))


def fit():
    shell = rank["shell"]
    single,_,phases,_ = shell["transported_spectra"]()
    v = rank["normalized_trial"](single)
    l,err_l = rank["fermion"](1)
    r,err_r = rank["fermion"](-1)
    ref = native["p"]
    theta = native["relative_phase"](q(676))
    a = [ref[i]*phases[i]*l[i,0]*r[i,0] for i in range(6)]
    d = [c.conjugate()*acb(0,-(-1)**i*theta).exp() for i,c in enumerate(a)]
    y = [v[i,0]/ref[i] for i in range(6)]
    # Isolated-vector residuals in the rank certificate are below 1e-11 each.
    # Bound the exact row, ground overlap, and orthogonal row component.
    ua,da = apply(a,y),apply(d,y)
    assert norm(a)+q("2e-11") < q(".408126")
    assert norm(a)-q("2e-11") > q(".4081258")
    assert abs(ua)-q("3e-11") > q(".4081258")
    tail = [a[i]-ua*y[i].conjugate() for i in range(6)]
    assert norm(tail)+q("7e-11") < q("1e-9")
    assert abs(da)-q("3e-11") > q(".4005869")
    assert abs(da)+q("3e-11") < q(".4005872")
    assert err_l+err_r < q("4e-15")
    coupling = 4*q.pi()
    up = [q(x) for x in (".00000615",".00311",".8616")]
    down = [q(x) for x in (".0000135",".000268",".01401")]
    U,D = norm(up),norm(down)
    # |alpha|=|d y/a y|, while the norm on the complement is at most
    # sqrt(||d||²-|dy|²)+|alpha| ||a(I-P)||. Both are invariant under field frames.
    alpha_low = q(".4005869")/q(".408126")
    alpha_high = q(".408126")/q(".4081258")
    fixed_upper = (q(".408126")**2-q(".4005869")**2).sqrt()+alpha_high*q("1e-9")
    assert fixed_upper < q(".078084")
    uniform_upper = q(".408126")+alpha_high*q("1e-9")
    # Accommodate the fermion and band unitary identifications by 1e-12.
    # The light coefficient and full row here are exact at the center.
    eps = q("1e-12")
    assert (1+alpha_high)*(err_l+err_r+q("1e-16")) < eps
    heavy_fixed = ((alpha_low*U-D)/coupling-eps)/fixed_upper
    heavy_free = ((norm(down)**2-apply(up,down)**2/U**2).sqrt()/coupling-eps)/uniform_upper
    light = (U/coupling-q("1e-9")-eps)/q(".408126")
    assert heavy_fixed > q(".8475") and heavy_free > q("4.247e-5") and light > q(".1679")
    # Explicit two-row minimum norm right inverse at the origin.
    beta = apply(d,[v.conjugate() for v in a])/norm(a)**2
    transverse = [x-beta*y for x,y in zip(d,a)]
    # beta error < 5e-10 from row errors 2e-11 and denominator > .408^2.
    assert abs(beta)+q("5e-10") < q(".98154")
    assert norm(transverse)-q("1e-9") > q(".078081")
    witness = ((U/q(".4081258"))**2+((D+q(".98154")*U)/q(".078081"))**2).sqrt()
    rinv = (1/q(".4081258")**2+(1+q(".98154")**2)/q(".078081")**2).sqrt()
    assert rinv < 19
    assert witness/(1-rinv*eps) < q("11.212")
    # Weighted spectral variance remains after subtracting any common squared mass.
    split = 3*q("60410.25")**2/q("1e-18")-q("25618")**2
    for T,bound in [(heavy_fixed,q("2.57e22")),(heavy_free,q("1.28e18"))]:
        assert split*light*T/q("246.1")**2 > bound
    print("field norm in 45 high modes >",heavy_fixed.lower().str(7),"or",
          heavy_free.lower().str(7),"for arbitrary reference angle; light norm >",light.lower().str(7))
    print("all central-mass flavor targets admit a common full profile of norm <",witness.upper().str(8))


if __name__ == "__main__":
    fit()
