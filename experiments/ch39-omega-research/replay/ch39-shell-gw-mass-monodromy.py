"""A full 53-cycle in the complex continuation of the one-loop shell matching mass."""
from math import gcd
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, ctx
import numpy as np

fixed = run_path(str(Path(__file__).with_name("ch39-shell-gw-fixed-scale-density.py")))
q, full, shell = (fixed[k] for k in ("q","full","shell"))


def spectral_symbol(index):
    single, _, _, G = shell["transported_spectra"]()
    phase = acb.zeta_zero(index).imag*((1-q(5)/26).acos()-(1-q(3)/26).acos())
    width = q("3e-8" if index == 43 else "1e-8")
    eye = shell["yuk"]["identity"](48)
    angle = q(2)*q.pi()/5
    central = [fixed["family"](single,G,phase,j*angle) for j in range(5)]
    tube = [fixed["family"](single,G,phase+arb(0,width),j*angle) for j in range(5)]
    # The singlet eigenvalue is independent of both angles. Certify a narrower interval.
    root = arb(float(np.linalg.eigvalsh(full["model"]["array"](single))[0]))
    lo, hi = root-q("1e-13"),root+q("1e-13")
    assert full["count"](full["rows"](single),lo)==0 and full["count"](full["rows"](single),hi)==1
    coeff = []
    err = q(0)
    for k in (0,1,2,-2,-1):
        omega = [(acb(0,-k*j*angle)).exp()/5 for j in range(5)]
        matrix = sum((c*v for c,v in zip(omega,central)), acb_mat(48,48))
        uncertain = sum((c*v for c,v in zip(omega,tube)), acb_mat(48,48))
        if k == 0:
            matrix -= root*eye
            uncertain -= lo.union(hi)*eye
        err += q(3)**k * sum(abs(z).abs_upper()**2 for z in (uncertain-matrix).entries()).sqrt()
        coeff.append(matrix*q(3)**k)
    return coeff, err.upper(), eye


def certify(index):
    ctx.prec = 160
    coeff, delta, eye = spectral_symbol(index)
    harmonics = (0,1,2,-2,-1)
    # 48-by-48 determinant with harmonics of degree at most 2 in each entry.
    N = 193
    roots = [acb(0,2*q.pi()*j/N).exp() for j in range(N)]
    samples = []
    movement = q(0)
    for j in range(N):
        matrix = sum((m*roots[(j*k)%N] for k,m in zip(harmonics,coeff)), acb_mat(48,48))
        center = matrix.det()
        inv = matrix.inv()
        inv_f = sum(abs(z).abs_upper()**2 for z in inv.entries()).sqrt()
        # D_p + eta*I about the reference node, for |eta| <= 1e-15.
        norm = inv_f*(delta+q(48).sqrt()*q("1e-15"))
        # Trace of each power in log det is bounded by norm**k via the
        # two Frobenius factors and the induced operator bound.
        assert norm < 1
        movement = arb.max(movement, (abs(center).abs_upper()*norm/(1-norm)).upper())
        samples.append(center)
    values = [sum(z*roots[(-j*k)%N] for j,z in enumerate(samples))/N for k in range(-96,97)]
    slope = sum(abs(k)*abs(c).abs_upper() for k,c in zip(range(-96,97),values))
    fine = 4096
    roots2 = [acb(0,2*q.pi()*j/fine).exp() for j in range(fine)]
    at = [sum(v*roots2[(j*k)%fine] for k,v in zip(range(-96,97),values)) for j in range(fine)]
    floor = min(abs(z).abs_lower() for z in at)-q.pi()*slope/fine
    span = sum(abs(v).abs_upper() for v in values)
    # The Lebesgue constant for odd trigonometric interpolation: directly
    # bound it by pairing the nearest grid point and the two distance tails.
    leb = 1+sum(q(1)/(k-(q(1)/2)) for k in range(1,97))
    print("window",index,"symbol perturbation",delta.str(6),"determinant min/max",floor.str(7),span.str(7),
          "interpolated relative determinant change",(movement*leb/floor).str(7),flush=True)
    assert floor>0
    assert slope*q.pi()*2/fine < floor
    winding = sum((at[(j+1)%fine]/z).arg() for j,z in enumerate(at))/(2*q.pi())
    assert winding.unique_fmpz()==10
    assert gcd(10,53)==1
    assert movement*leb < floor/10
    assert span < q("1e34")
    # All five shifted singlet gaps are between 3 and 7.
    # The earlier mass-law certificate gives log(C/h^53)>2484.
    radius = (q(2484).exp()/(7**5*span))**(q(1)/53)
    eps = q(".005")
    assert 1/((1-eps)*radius) < q("1e-15")
    singlet_shift = (1+q("1e-15")/3)**5-1
    total_error = (1+movement*leb/floor)*(1+singlet_shift)-1
    comparison = 2*53*eps-((1+eps)**53-1)
    assert (1+eps)**53*total_error < comparison
    assert eps < (q.pi()/53).sin()
    assert 53*q(7001).log() < 2484
    print("certified cycle winding",winding.str(6))


if __name__=="__main__":
    for n in (43,195):
        certify(n)
