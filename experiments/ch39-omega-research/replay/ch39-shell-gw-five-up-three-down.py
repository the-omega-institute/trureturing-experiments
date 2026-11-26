"""Exclude five up modes with three down modes using the top and bottom left lines."""
from itertools import permutations
from pathlib import Path
from runpy import run_path

import numpy as np
from flint import acb_mat, arb

m = run_path(str(Path(__file__).with_name("ch39-shell-gw-asymmetric-bands.py")))
four = m["four"]
adj = m["adj"]
d = arb(".00042")
kappa = arb("3.06")


def certify_shapes():
    raw = m["images"](5)
    modes = m["orb"]["orthonormalize"](raw)
    y0,y1 = map(four["float_matrix"],raw[:2])
    left = four["unitary_trial"](np.linalg.eig(y1 @ np.linalg.inv(y0))[1])
    right = four["unitary_trial"](y0.conj().T @ left)
    L,R = four["exact_frame"](left),four["exact_frame"](right)
    modes = [adj(L)*b*R for b in modes]
    for i,j,_ in permutations(range(3)):
        assert sum(abs(b[i,j])**2 for b in modes[:3]).sqrt() < arb(".000415")
    for p in range(3):
        i,j = [a for a in range(3) if a != p]
        def entries(positions):
            return acb_mat([[b[r,c] for b in modes[3:]] for r,c in positions])
        other = entries([(i,j),(j,i)])
        for edge in (entries([(i,p),(j,p)]), entries([(p,i),(p,j)])):
            transfer = edge*other.inv()
            form = kappa**2*m["eye"](2) - adj(transfer)*transfer
            assert form[0,0].real > 0 and form.det().real > 0


if __name__ == "__main__":
    certify_shapes()
    ru,rd = arb(".00362"),arb(".01914")
    tol = arb("1e-6")
    assert arb(".00000615")/arb(".8616")*(1+tol)/(1-tol) < arb(".0000072")
    tail = (ru**2+arb(".0000072")**2).sqrt()
    root2 = arb(2).sqrt()
    assert arb(".000415")*(1+2*rd**2).sqrt() < d
    # Low-three down sector: its top line must be very close to one coordinate.
    db = root2*d/(1-rd-4*d)
    assert db < arb(".00061")
    x = arb(".044")
    mix = (1-((arb(".04275")+arb(0,tol)).cos()*(arb(".00377")+arb(0,tol)).cos())**2).sqrt()
    assert mix+db < x and mix-db > arb(".042")
    # x bounds the left-axis distance, y bounds the right-axis distance.
    fixed = kappa*(tail+root2*d)+root2*d
    y = (fixed+ru*x)/((1-x*x).sqrt()-kappa*x)
    assert y < arb(".016")
    tight = (kappa*x*y+fixed+ru*y)/(1-y*y).sqrt()
    assert tight < arb(".016") < mix-db
    print("required up left-line distance above",(mix-db).lower().str(7),
          "five-mode bound below",tight.upper().str(7))
