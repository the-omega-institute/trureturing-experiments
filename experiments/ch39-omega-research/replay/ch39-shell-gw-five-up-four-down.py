"""Complete the non-dominated band allocations by excluding five up and four down."""
from pathlib import Path
from runpy import run_path
from flint import arb

five = run_path(str(Path(__file__).with_name("ch39-shell-gw-five-up-three-down.py")))
four = five["four"]


if __name__ == "__main__":
    # Both programs use the polar frame of the first adjoint/singlet pencil.
    # Certify the down-sector entry bounds in the very same five-mode frame.
    import numpy as np
    raw = five["m"]["images"](5)
    y0,y1 = map(four["float_matrix"],raw[:2])
    L = four["unitary_trial"](np.linalg.eig(y1 @ np.linalg.inv(y0))[1])
    R = four["unitary_trial"](y0.conj().T @ L)
    L,R = four["exact_frame"](L),four["exact_frame"](R)
    smaller = five["m"]["orb"]["orthonormalize"](raw[:4])
    smaller = [five["adj"](L)*b*R for b in smaller]
    forward = [(0,1),(1,2),(2,0)]
    if abs(smaller[3][0,1]) < arb(".2"):
        forward = [(j,i) for i,j in forward]
    for i,j in forward:
        assert arb(".563") < abs(smaller[3][i,j]) < arb(".565")
        assert abs(smaller[3][j,i]) < arb(".125")
        for u,v in ((i,j),(j,i)):
            assert sum(abs(b[u,v])**2 for b in smaller[:3]).sqrt() < arb(".000415")
    five["certify_shapes"]()
    d,K = five["d"],five["kappa"]
    ru = arb(".00362")
    db = four["heavy_line_radius"](arb(".01914"),d)
    tol = arb("1e-6")
    mix = (1-((arb(".04275")+arb(0,tol)).cos()*(arb(".00377")+arb(0,tol)).cos())**2).sqrt()
    x = arb(".066")
    assert mix+db < x and mix-db > arb(".020")
    tail = (ru**2+arb(".0000072")**2).sqrt()
    fixed = K*(tail+arb(2).sqrt()*d)+arb(2).sqrt()*d
    y = (fixed+ru*x)/((1-x**2).sqrt()-K*x)
    tightened = (K*x*y+fixed+ru*y)/(1-y**2).sqrt()
    assert tightened < arb(".018") < mix-db
    print("required top-axis distance above",(mix-db).lower().str(8),
          "five-mode upper bound",tightened.upper().str(8))
