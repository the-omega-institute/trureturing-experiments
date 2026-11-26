"""Quartic auxiliary potential and balanced abelian Coulomb equation for shell matching."""
from pathlib import Path
from runpy import run_path

from flint import arb, arb_mat
import sympy as sp

flow = run_path(str(Path(__file__).with_name("ch39-shell-gw-first-heavy-flow.py")))
q = arb


def quartic_encoding():
    powers = (0,1,2,4,6,7)
    uv = sp.symbols("u0:6")
    coords = dict(zip(powers,(1,*uv[:5])))
    terms = {}
    for a in powers:
        for b in powers:
            terms[a+b] = coords[a]*coords[b]
    assert all(k in terms for k in range(14))
    x = sp.Symbol("x")
    coeff = sp.symbols("a0:54")
    blocks = [sum(coeff[i]*terms[i%14] for i in range(j*14,min(54,(j+1)*14)))
              for j in range(4)]
    # Coefficients stay symbolic, so one check covers all shell polynomials.
    s1,s2,s3 = sp.symbols("s1 s2 s3")
    t = uv[5]
    equations = [uv[1]-uv[0]**2, uv[2]-uv[1]**2, uv[3]-uv[1]*uv[2],
                 uv[4]-uv[0]*uv[3], t-uv[4]**2, s3-blocks[3],
                 s2-t*s3-blocks[2], s1-t*s2-blocks[1], t*s1+blocks[0]]
    variables = (*uv,s1,s2,s3)
    assert len(variables)==len(equations)==9
    assert all(sp.Poly(eq,*variables).total_degree()<=2 for eq in equations)
    mass_powers = {coords[i]:x**i for i in powers if i>0}
    replaced = [b.xreplace(mass_powers) for b in blocks]
    h = replaced[3]
    for b in reversed(replaced[:3]):
        h = x**14*h+b
    assert sp.expand(h-sum(c*x**i for i,c in enumerate(coeff)))==0


def energy_range():
    inv = [4*q.pi()/s for s in flow["start"]]
    b = flow["beta"]
    mat = arb_mat([[b[0]-b[1],-q(1)/15],[b[1]-b[2],q(1)/6]])
    t,E = mat.solve(2*q.pi()*arb_mat([[inv[0]-inv[1]],[inv[1]-inv[2]]])).entries()
    Q = 1000*t.exp()
    logC = 106*((q(1000)/125).log()+t)-2*E
    maxmass = q("2.435e18")
    forced = maxmass*(-logC/53).exp()
    copies = logC/(maxmass/Q).log()
    print("max light mass GeV with positive-charge masses below reduced Planck",forced.str(10),
          "balanced pair count required for every mass between Q and Planck exceeds",copies.str(10))
    assert forced < q(".011") and 230 < copies < 231


if __name__=="__main__":
    quartic_encoding()
    energy_range()
