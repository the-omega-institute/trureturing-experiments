"""Certified high-zero labels satisfying a common bare Newton and gauge prefactor."""
from fractions import Fraction as F
from pathlib import Path
from runpy import run_path
from flint import acb, arb

ROOT = Path(__file__).resolve().parent
geo = run_path(str(ROOT/'ch39-regge-gauge-scale.py'))
model, q = geo['model'], geo['q']
L = F(429, 16)
V = (F(24,25), F(16,15), F(1))
OFFSET = F(1, 4)+30+F(2, 15)


def verify():
    mu = 1000*q(L).exp()
    prefactor = (geo['BOOK_PLANCK']/mu)**2/(8*arb.pi())
    coeff = 4*prefactor*(3*arb(5).sqrt()/2)
    centers = model['GAUGE']
    radius = F(1, 1000)
    box = [q(c*(1-radius)).union(q(c*(1+radius))) for c in centers]
    def end(inputs):
        return model['running']['integrate'](inputs, int(64*L))[1]
    all_values = end(box)
    necessary_low = max((x/q(v)).lower() for x, v in zip(all_values,V))
    necessary_high = min((x/q(v)).upper() for x, v in zip(all_values,V))
    lower_face, upper_face = [], []
    for i, c in enumerate(centers):
        for side in (-1, 1):
            face = box.copy()
            face[i] = q(c*(1+side*radius))
            target_floor = end(face)[i]/q(V[i])
            if side == -1:
                lower_face.append(target_floor.lower())
            else:
                upper_face.append(target_floor.upper())
    sufficient_low, sufficient_high = max(upper_face), min(lower_face)
    assert necessary_low < sufficient_low < sufficient_high < necessary_high
    def height(g3):
        return (coeff/g3-q(OFFSET)).sqrt()
    outer = (height(necessary_high), height(necessary_low))
    inner = (height(sufficient_high), height(sufficient_low))
    def labels(interval):
        counts = [x.zeta_nzeros().unique_fmpz() for x in interval]
        assert all(v is not None for v in counts)
        return list(range(int(counts[0])+1, int(counts[1])+1))
    possible, guaranteed = labels(outer), labels(inner)
    assert possible == guaranteed == list(range(11485,11491))
    assert q('11134.68') < outer[0] < inner[0] < q('11134.92')
    assert q('11138.83') < inner[1] < outer[1] < q('11139.24')
    for j in guaranteed:
        z = acb.zeta_zero(j)
        assert z.real == q(F(1,2))
        g3 = coeff/(z.imag**2+q(OFFSET))
        assert sufficient_low < g3 < sufficient_high
    print('Necessary heights and labels:', [x.str(11) for x in outer], possible)
    print('Sufficient heights and labels:', [x.str(11) for x in inner], guaranteed)
    print('Common microscopic prefactor:', prefactor.str(11))


if __name__ == '__main__':
    verify()
