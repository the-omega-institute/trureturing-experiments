"""完整十九模在实际和平坦连接上的四维共享边积分，带傅里叶尾界。"""
from pathlib import Path
from runpy import run_path

import numpy as np
from flint import acb, acb_mat, arb, arb_poly, ctx


def angle_tables(poly, period, flat):
    derivative = poly.derivative()
    face, slope = np.empty(period), np.empty(period)
    for i in range(period):
        if i == 0:
            face[i], slope[i] = 1, 0
            continue
        if i == period//2 and flat:
            # The constant mode has b=1 in the flat graph.
            assert poly(1).contains(0)
            face[i], slope[i] = 0, 0
            continue
        if i == period//2:
            u, du = arb(1), arb(0)
        else:
            sine, cosine = (arb.pi()*i/period).sin_cos()
            u, du = sine*sine, sine*cosine
        value = poly(u)
        assert value > 0
        for arr, x in ((face, value**6), (slope, -derivative(u)*du/value)):
            arr[i] = float(x)
            assert x == 0 or abs(arb(arr[i])-x) < arb("1e-15")*abs(x)
    assert np.max(np.abs(slope)) < 1100
    return face, slope


def alias_bound(trace_bound):
    # Six copies of each eigenangle factor are the characteristic function
    # of independent signed integer steps. Bound their one-coordinate tails.
    t, u = arb(".183"), arb(".690")
    hyper = sum(d*((q*t).cosh()-1) for q,d in ((1,6),(4,3),(2,3),(3,2),(6,1)))
    small = 8*(u.cosh()-1)
    tail = 2*((3*trace_bound*hyper-840*t).exp()
              + 3*(3*trace_bound*small-206*u).exp())
    assert tail < arb("6e-32")
    return 256*tail


def verify():
    directory = Path(__file__).parent
    grid = run_path(str(directory/"ch39-three-mode-haar-grid.py"))
    shell = run_path(str(directory/"ch39-three-ray-chiral-boundary.py"))
    ctx.prec = 256
    height = acb.zeta_zero(11431).imag
    identity = acb_mat(19,19)
    for i in range(19):
        identity[i,i] = 1
    charges = (1,2,3,4,6)
    walks = {q:identity-shell["laplacian"](q,height)/36 for q in (0,*charges)}
    _, _, weyl = grid["tables"](840)
    bound = alias_bound(arb("11.501"))
    for steps in (1,2):
        data = {}
        for q, walk in walks.items():
            block = walk**steps
            assert sum(block[i,i].real for i in range(19)) < arb("11.501")
            raw = block.charpoly()
            assert all(c.imag.contains(0) for c in raw)
            poly = arb_poly([c.real for c in reversed(list(raw))])
            data[q] = angle_tables(poly, 840, q == 0)
        results = []
        for flat in (False, True):
            pairs = {q:data[0 if flat else q] for q in charges}
            z, norms = grid["grid_moments"](840,210,210,pairs,weyl)
            # Numerators before Cartan normalization. Each of the two
            # differentiated copies has multiplier bounded by K in modulus.
            vacuum = z+(-bound).union(bound)
            assert vacuum > 0
            moments = []
            for rms, k, terms in zip(norms,(684,152,152),(1,1,3)):
                err = bound*k*k*terms
                moments.append((z*rms*rms+(-err).union(err))/vacuum)
            y, w, c = moments
            ratios = [12*y/w, c/(4*w)]
            if steps == 2 and not flat:
                assert arb("60.29235") < ratios[0] < arb("60.29245")
                assert arb(".9987274") < ratios[1] < arb(".9987286")
            scale, strong = grid["one_loop"](ratios)
            assert arb(".07078") < strong < arb(".07104")
            results.append((ratios[0], ratios[1], strong))
            print("步数、平坦对照、超荷/弱、色/弱、匹配能标、强耦合",
                  steps, flat, *ratios, scale, strong, flush=True)
        changes = [x-y for x,y in zip(*results)]
        if steps == 2:
            assert arb(".0098") < changes[0] < arb(".0100")
            assert arb(".000081") < changes[1] < arb(".000083")
            assert arb("-.000037") < changes[2] < arb("-.000035")
        print("实际与平坦的三个差值", *changes, flush=True)


if __name__ == "__main__":
    verify()
