"""Minimum edge changes forced by the top/bottom hierarchy under sitewise adjoint conjugations."""
from pathlib import Path
from runpy import run_path

from flint import arb

orbit = run_path(str(Path(__file__).with_name("ch39-shell-gw-cone-adjoint-links.py")))
q,rank,shell = arb,orbit["rank"],orbit["shell"]


def certify():
    _,a,d,targets,_ = orbit["setup"]()
    common = orbit["common"]
    size = lambda M: shell["hs"](M,M).real.sqrt()
    U,D = [size(m) for m in targets]
    single,_,_,_ = shell["transported_spectra"]()
    spatial = rank["normalized_trial"](single)
    phases = [spatial[i,0]/orbit["native"]["p"][i] for i in range(6)]
    assert max(abs(v-phases[0]) for v in phases)+q("2e-11") < q("3e-11")
    # Center phases can be made common. Residual from the rows, equal-radius
    # replacement and rank-one approximation is absorbed in 1e-7 per output.
    error = q("1e-7")
    Rmax = 4*q.pi()/q(6).sqrt()
    assert 4*q.pi()*q("1.01e-9")+6*Rmax*(q("2e-10")+q("1e-14")) < error
    det = abs(a[0]*d[1]-a[1]*d[0])
    for i in (0,1):
        assert abs(d[i]) > q(".166") and abs(a[i]) < q(".167")
        assert (abs(d[i])*(U-error)-abs(a[i])*(D+error))/det > q("13.29")
        assert (abs(d[i])*(U+error)+abs(a[i])*(D+error))/det < q("13.75")
    assert abs(d[0]/d[1]-1)<q(".384")
    eps = (D+error)/abs(d[1])
    # ||X_e+X_o|| after normalization is at most |ratio-1|+2eps/L
    # up to the modulus mismatch; near-unit row ratio is checked separately.
    assert abs(abs(d[0]/d[1])-1)<q("1e-10") and eps<q(".086")
    chord = q(".384")+2*eps/q("13.29")+q("2e-10")
    angle = q.pi()-2*(chord/2).asin()
    for R,edge_min in [(Rmax,"1.126"),(12/q(6).sqrt(),"1.337")]:
        L = q("13.29")/R
        opening = ((L**2-3)/(2*L)).acos()
        assert 2*((angle-2*opening)/2).sin()-q("1e-15") > q(edge_min)
    internal = ((9-(q("13.75")/(12/q(6).sqrt()))**2)/6).sqrt()
    assert internal > q(".4325")
    # Frobenius cost of mixed links: 9-|tr W|² >= 1.5 ||Ad W - I||².
    assert (6*q("1.337")**2+4*q(".4325")**2)/72 > q(".159")
    print("all 6 cross-parity adjoint edge displacements above 1.126 at strength <=4pi;",
          "at strength 12 all six above 1.337 and >=4 long edges above .4325;",
          "mean normalized mixed trace loss > .159")


if __name__ == "__main__":
    certify()
