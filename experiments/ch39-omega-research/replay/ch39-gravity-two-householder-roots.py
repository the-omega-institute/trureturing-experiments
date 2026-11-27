"""右外壳的两条荷六源像射线以秩二厄米对合实化，并检验相容的三次根。"""
from itertools import product
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, ctx


def adj(matrix):
    return matrix.conjugate().transpose()


def norm(vector):
    return sum(abs(x)**2 for x in vector.entries()).sqrt()


def line(vector):
    return vector*adj(vector)*(1/(norm(vector)**2))


def gravity(root, *, angles=(), with_motion=False):
    joint = run_path(str(root / "ch39-one-generator-up-and-gravity.py"))
    base = run_path(str(root / "ch39-rank-two-generator-up-quark.py"))
    shell = run_path(str(root / "ch39-three-ray-chiral-boundary.py"))
    spectral = run_path(str(root / "ch39-adjacent-pvm-transport.py"))
    orbits = run_path(str(root / "ch39-up-quark-arbitrary-site-phase-floor.py"))["ORBITS"]
    _, parameters, gev, cutoffs = joint["verify"](root, return_parameters=True)
    a,b = parameters
    low, high, peak = cutoffs
    ctx.prec = 384
    n = base["N"]+1
    source = spectral["spectrum"](shell,6,n)[1]
    target = spectral["spectrum"](shell,6,n+1)[1]
    p = base["projector"](orbits)
    identity = base["eye"](19)
    e = ((gev/arb("1.22089e19")-low)/(high-low)).sqrt()
    basis = acb_mat([[1/arb(len(g)).sqrt() if i in g else 0 for g in orbits] for i in range(19)])
    fixed, offset, small, big, tangent, length = joint["mode_parameters"](p,basis,source,target)
    lo = offset + small*e + big*a
    near = offset + small*e.mid() + big*a.mid()
    radius = (length-norm(lo)**2).sqrt()
    y0 = fixed + near + tangent*(length-norm(near)**2).sqrt()
    y = y0
    delta = arb(0)
    variations = []
    for var, direction in ((e,small),(a,big)):
        derivative = (direction-tangent*((adj(lo)*direction)[0,0].real/radius))
        variations.append((derivative, abs(var-var.mid()).upper()))
        delta += norm(derivative).upper()*abs(var-var.mid()).upper()
        y += derivative*(var-var.mid())
    s,t = [joint["col"](source,j) for j in (17,18)]
    t3,t19 = [joint["align"](joint["col"](target,j),t) for j in (2,18)]
    z = t19*b + t3*(1-b*b).sqrt()
    original_y = y
    if angles:
        modes = (0,3,6,9,12,15,17)
        phases = dict(zip(modes, angles, strict=True))
        coefficients = adj(target)*y
        factors = [acb(0,phases.get(j,0)).exp() for j in range(19)]
        y = target*acb_mat([[coefficients[j,0]*factors[j]] for j in range(19)])
        rephased_s = target*acb_mat([[a*b.conjugate()] for a,b in zip((adj(target)*s).entries(),factors,strict=True)])
        central = adj(target)*y0
        y0 = target*acb_mat([[central[j,0]*factors[j]] for j in range(19)])
    images = (joint["align"](y,s), joint["align"](z,t))
    sources = (s,t)
    prior = [low + (high-low)*abs((adj(joint["col"](target,j))*original_y)[0,0])**2
             + (peak-low)*abs((adj(joint["col"](target,j))*z)[0,0])**2 for j in range(19)]
    branches = {}
    for signs in product((1,-1), repeat=2):
        diffs = [src-sgn*dst for src,dst,sgn in zip(sources,images,signs,strict=True)]
        assert all(norm(d) > 0 for d in diffs)
        assert (adj(diffs[0])*diffs[1])[0,0].contains(0)
        q = sum((line(d) for d in diffs),acb_mat(19,19))
        h = identity-2*q
        assert (adj(h)-h).contains(acb_mat(19,19)) and (h*h).contains(identity)
        assert q.trace().contains(2) and (q*q-q).contains(acb_mat(19,19))
        for src,dst,sgn in zip(sources,images,signs,strict=True):
            assert (h*src-dst*sgn).contains(acb_mat(19,1))
        proj = adj(target)*h
        actual = [low+(high-low)*abs(u)**2+(peak-low)*abs(v)**2
                  for u,v in zip((proj*s).entries(),(proj*t).entries(),strict=True)]
        assert all((a-b).contains(0) for a,b in zip(actual,prior,strict=True))
        branches[signs] = diffs
    print("两个源像内积的模", [abs((adj(src)*dst)[0,0]).str(10)
                                  for src,dst in zip(sources,images,strict=True)])
    print("实化符号分支数、每支一致的荷六目的层数", len(branches), len(prior))
    if with_motion:
        center = joint["align"](y0,s)
        overlap = abs((adj(y0)*s)[0,0])
        assert overlap > delta
        ray = rephased_s if angles else s
        overlap_delta = sum(abs((adj(ray)*derivative)[0,0]).upper()*dv for derivative,dv in variations)
        aligned_delta = delta + overlap_delta/(overlap-overlap_delta)
        print("荷六第一像及重相后的位移上界", delta.str(13), aligned_delta.str(13))
        return branches, s+center, aligned_delta
    return branches


if __name__ == "__main__":
    gravity(Path(__file__).resolve().parent)
