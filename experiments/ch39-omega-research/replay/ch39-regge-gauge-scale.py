"""Test a single bare action prefactor and a matching-scale coordinate spacing."""
from fractions import Fraction as F
from json import loads
from pathlib import Path
from runpy import run_path
from flint import arb

ROOT = Path(__file__).resolve().parent
cross = run_path(str(ROOT/'ch39-monopole-free-matching.py'))
model, q = cross['model'], cross['q']
t = F(1, 15)
BOOK_PLANCK = q('1.2209e19')


def verify():
    last = loads((ROOT/'ch39-tick-separators.json').read_text())['prefix'][-1][0]
    height = q(F(last, 1 << 400)).union(q(F(last+1, 1 << 400)))
    bare = 4/(height**2+q(F(1, 4)+30+2*t))
    volume_factor = 3*arb(5).sqrt()/2
    g3 = q(F(13505, 4096))
    mu = 1000*q(F(429,16)).exp()
    ratio = (8*arb.pi()*g3/(bare*volume_factor)).sqrt()
    prediction = mu*ratio
    required_height = (4*volume_factor*(BOOK_PLANCK/mu)**2/(8*arb.pi()*g3)-q(F(1, 4)+30+2*t)).sqrt()
    assert q('4.0900e18') < prediction < q('4.0901e18')
    assert q('8.9106') < (BOOK_PLANCK/prediction)**2 < q('8.9108')
    assert q('11137.56') < required_height < q('11137.57')
    print('Witness: Planck/matching scale, predicted Planck, prefactor ratio, required height:',
          ratio.str(9), prediction.str(9), ((BOOK_PLANCK/prediction)**2).str(9), required_height.str(9))
    upper = []

    def observe(step, inverse, deriv):
        if step >= 23*64 and 0 in (inverse[1]-q('16/15')*inverse[2]):
            stop = q(F(step+1,64))
            upper.append(1000*stop.exp()*(8*arb.pi()*inverse[2]/(bare*volume_factor)).sqrt())

    centers = model['GAUGE']
    gauge = [q(v*F(999, 1000)).union(q(v*F(1001, 1000))) for v in centers]
    model['running']['integrate'](gauge, 33*64, observe)
    assert upper
    cap = max(x.upper() for x in upper)
    assert cap < q('4.430e18') and (BOOK_PLANCK/cap)**2 > q('7.5987')
    print('Full input box and all eligible matching scales, Planck upper:', cap.str(10),
          'minimum prefactor mismatch', ((BOOK_PLANCK/cap)**2).str(10))


if __name__ == '__main__':
    verify()
