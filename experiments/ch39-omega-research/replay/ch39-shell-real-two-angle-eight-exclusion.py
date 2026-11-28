"""消去未观测角与 CP 相位，检验壳实协方差的含 theta23 双角八项闭窗。"""
from itertools import product
from pathlib import Path
from runpy import run_path
import sys

from flint import arb


def power(value, exponent, affine):
    interval = value.bound()
    assert interval > 0
    slope = exponent * value.c ** (exponent - 1)
    second = abs(exponent * (exponent - 1) * interval ** (exponent - 2))
    return affine(value.c ** exponent, value.g * slope,
                  abs(slope) * value.e + second * value.rad ** 2 / 2)


def recover(values, missing, targets, shared):
    """正逆两迹在缺失角的正弦平方及带号干涉项上是两个一次方程。"""
    affine = shared["Affine"]
    variables = shared["inputs"](values)
    masses = variables[:6]
    ratios = (masses[0]/masses[3]) * (masses[1]/masses[4]) * (masses[2]/masses[5])
    norm = power(ratios, arb(2)/3, affine)
    a, b, z = [x.sincos()[0] for x in variables[6:9]]
    a, b, z = a*a, b*b, z*z
    if missing == 6:
        a = affine(0)
    else:
        z = affine(0)
    constant = (
        ((1-a)*(1-z), a*(1-z), z),
        (a*(1-b)+(1-a)*b*z, (1-a)*(1-b)+a*b*z, b*(1-z)),
        (a*b+(1-a)*(1-b)*z, (1-a)*b+a*(1-b)*z, (1-b)*(1-z)),
    )
    lines = []
    for inverse, target in enumerate(targets):
        up = [x*x for x in masses[:3]]
        down = [(x*x).recip() for x in masses[3:]]
        scale = norm
        if inverse:
            up = [x.recip() for x in up]
            down = [x*x for x in masses[3:]]
            scale = norm.recip()
        interference = (up[1]-up[2])*(down[0]-down[1])
        if missing == 8:
            slope = (up[0]-b*up[1]-(1-b)*up[2])*(down[2]-(1-a)*down[0]-a*down[1])
        else:
            slope = (down[0]-down[1])*(-up[0]*(1-z)+up[1]*((1-b)-b*z)+up[2]*(b-(1-b)*z))
        offset = sum(up[i]*down[j]*constant[i][j] for i,j in product(range(3),repeat=2))
        lines.append((slope/interference, (scale*target-offset)/interference))
    (k, l), (inverse_k, inverse_l) = lines
    square = (l-inverse_l)/(k-inverse_k)
    cross = l-k*square
    if missing == 6:
        a = square
    else:
        z = square
    cosine = cross/((a*(1-a)*b*(1-b)*z).sqrt()*2)
    return square, cosine


def cover(shared, targets, missing, sign, key, bits):
    model = shared["model"]
    center = [0.]*9
    radius = [1.]*9
    for j, axis in enumerate((0,1,3)[:bits]):
        radius[axis] /= 2
        center[axis] += radius[axis] * (1 if (key >> j) & 1 else -1)
    stack = [(center, radius, 0)]
    leaves = deepest = 0
    while stack:
        center, radius, depth = stack.pop()
        values = [c+w*(arb(t)+arb(0,h)) for c,w,t,h in
                  zip(model["ref"][:9],model["tol"][:9],center,radius,strict=True)]
        square, cosine = recover(values, missing, targets, shared)
        x, y = square.bound(), cosine.bound()
        physical = not (x > 1 or x < 0 or y.lower() > 1 or y.upper() < -1)
        excluded = not physical
        if physical:
            angle = x.intersection(arb(0).union(arb(1))).sqrt().asin()
            phase = y.intersection(arb(-1).union(arb(1))).acos()
            values[missing] = angle
            if 0 < phase < arb.pi() and angle > 0:
                try:
                    excluded = shared["test"](values + [sign*phase])
                except (AssertionError, ValueError, ZeroDivisionError):
                    excluded = False
        if excluded:
            leaves += 1
            deepest = max(deepest,depth)
            continue
        assert depth < 45, (missing, sign, key, center, radius)
        priorities = (.002,.02,.02,.6,.3,.3,.06,.2,.2)
        weight = [f*r + 8*float(abs(cosine.g[0,i]).upper())
                  + (300 if missing==6 else 3e6)*float(abs(square.g[0,i]).upper())
                  for i,(f,r) in enumerate(zip(priorities,radius,strict=True))]
        weight[missing] = 0
        axis = max(range(9),key=weight.__getitem__)
        half = radius.copy()
        half[axis] /= 2
        for direction in (-1,1):
            child = center.copy()
            child[axis] += direction * half[axis]
            stack.append((child,half,depth+1))
    return leaves, deepest+bits


def verify(replay):
    shared = run_path(str(replay/"ch39-shell-real-nine-window-exclusion.py"),run_name="__main__")
    normalized = run_path(str(replay/"ch39-two-zero-sharp-quark-fit.py"))["normalized_traces"]
    targets = normalized(shared["up"], shared["down"])
    model = shared["model"]
    signs = (-1,1,1,-1,1,1,1,1,1)
    optimum = arb(".0574339450674740").union(arb(".0574339450674742"))
    corner = [v+s*w*optimum for v,w,s in zip(model["ref"][:9],model["tol"][:9],signs,strict=True)]
    actual = run_path(str(replay/"ch39-quadratic-gram-source.py"))["at_zero"]
    traces = normalized(actual(model,-4,2345666), actual(model,2,964890))
    known = run_path(str(replay/"ch39-nine-window-cp-cosine.py"))["phase"](corner,traces[0])
    for missing in (6,8):
        square, cosine = recover(corner,missing,traces,shared)
        assert square.bound().overlaps(corner[missing].sin()**2)
        assert cosine.bound().overlaps(known)
    full = [v+arb(0,w.upper()) for v,w in zip(model["ref"][:9],model["tol"][:9],strict=True)]
    for missing in (6,8):
        square, cosine = recover(full,missing,targets,shared)
        assert 0 < square.bound() < 1 and -1 < cosine.bound() < 1
    if len(sys.argv)>2:
        missing, sign, key = map(int,sys.argv[2:5])
        assert missing in (6,8) and sign in (-1,1)
        bits = 0 if (missing,sign)==(8,-1) else 3
        assert 0 <= key < 2**bits
        part = cover(shared,targets,missing,sign,key,bits)
        print("缺失角索引、CP符号、分区、叶数与总深度", missing,sign,key,*part)
        return
    for missing,sign in product((6,8),(-1,1)):
        bits = 0 if (missing,sign)==(8,-1) else 3
        parts = [cover(shared,targets,missing,sign,key,bits) for key in range(2**bits)]
        print("缺失角索引、CP符号、叶数与总深度", missing,sign,
              sum(n for n,_ in parts),max(depth for _,depth in parts),flush=True)


if __name__ == "__main__":
    verify(Path(sys.argv[1]).resolve() if len(sys.argv)>1 else Path(__file__).resolve().parent)
