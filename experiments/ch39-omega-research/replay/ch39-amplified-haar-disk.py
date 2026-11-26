"""Positive-definite Fourier and entropy bounds for the amplified full SU(5) three-face disk."""
from fractions import Fraction as F
from json import loads
from pathlib import Path
from runpy import run_path

ROOT = Path(__file__).resolve().parent
amplified = run_path(str(ROOT/'ch39-microscopic-center-amplification.py'))
model, q, power, t = [amplified[k] for k in ('model', 'q', 'POWER', 'T')]
arb = model['arb']


def verify():
    cell = loads((ROOT/'ch39-tick-separators.json').read_text())['prefix'][3199][0]
    lo, hi = (q(F(k, 1 << 400)) for k in (cell, cell+1))
    s = 1/(lo.union(hi)**2+q(F(1, 4)))
    # Unit Frobenius ball: an ambient 50-dimensional packing gives its Haar mass >= 6^-50.
    f_min = q(30+2*t-9*(1+t))
    log_partition = power*(1+s*f_min).log()-50*arb(6).log()
    f_mean = log_partition/(power*s)
    dim_holder = (2*arb(5).sqrt()+2*arb(10).sqrt())**2
    nonconstant = (f_mean/q(1+t))**3/dim_holder
    assert f_mean > 0 and nonconstant > 0
    # Minimum normalized Casimir among nontrivial SU(5) representations.
    fundamental_casimirs = [F(k*(5-k)*6, 10) for k in range(1, 5)]
    assert min(fundamental_casimirs)/24 == F(1, 10)
    extra = nonconstant/(10*(1+nonconstant))
    bare = 4*power*s/(1+q(30+2*t)*s)
    calibrated = 35*extra/(2*bare+extra)
    assert f_mean > q('11.884') and nonconstant > q('11.864')
    assert extra > q('.09222') and calibrated > q('.03891')
    print('Normalized expectation lower bound, nonconstant disk mass, extra stiffness, two-color residual:',
          f_mean.str(9), nonconstant.str(9), extra.str(9), calibrated.str(9))


if __name__ == '__main__':
    verify()
