"""原十九点壳实三角的完整频率、几何偏角、双钟反演与首例分辨率。"""
from collections import Counter
from itertools import combinations
from math import gcd
from pathlib import Path
from runpy import run_path
import sys

from flint import acb, arb, ctx


root = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).resolve().parent
graph = run_path(str(root / "ch39-arc-relative-determinant.py"))
R = graph["R"]
axis_x = ((0, 7, 16), (1, 6, 17))
axis_y = ((1, 8, 9), (0, 8, 10))
diagonal = (0, 7, 8)


def conjugate(pair):
    return pair[0], -pair[1]


def product(z, w):
    a, b = z
    c, d = w
    real, imag = a*c-R*b*d, a*d+b*c
    common = gcd(real, imag)
    return real//common, imag//common


def unit(pair):
    z = acb(pair[0], arb(R).sqrt()*pair[1])
    return z/abs(z)


def enumerate_faces():
    sites = graph["sites"]
    pole = (17, 17, 17)
    edges = {}
    for i, j, length, native in graph["edges"]:
        x, y = sites[i], sites[j]
        real = R*(R+graph["dot"](pole, x)+graph["dot"](x, y)+graph["dot"](y, pole))
        imag = -graph["det"](pole, x, y)
        z = (-real, -imag) if length == 200 else (real, imag)
        assert (native-unit(z)).contains(0)
        edges[i, j] = length, z
    faces = {}
    for i, j, k in combinations(range(19), 3):
        if not all(e in edges for e in ((i, j), (j, k), (i, k))):
            continue
        first, second, last = (edges[e] for e in ((i, j), (j, k), (i, k)))
        m = tuple(int(first[0] == clock)+int(second[0] == clock)-int(last[0] == clock)
                  for clock in (724, 804))
        z = product(product(first[1], second[1]), conjugate(last[1]))
        if m < (0, 0):
            m, z = tuple(-n for n in m), conjugate(z)
        faces[i, j, k] = m, z
    assert Counter(m for m, _ in faces.values()) == {(0, 0): 433, (1, 0): 42, (0, 1): 42, (1, -1): 12}
    return faces


def closest_critical(pair):
    # 余弦的最近驻点在 -arg(u) mod π。
    z = unit(pair)
    assert not z.real.contains(0)
    if z.real < 0:
        z = -z
    return -z.arg()


def reconstruction(faces):
    radius, aperture = arb.pi()/16, arb(".256255")
    for m, z in faces.values():
        if m != (0, 0):
            critical = closest_critical(z)
            assert abs(critical) < arb(".156") < radius
    scales = []
    for pair, direction in ((axis_x, (1, 0)), (axis_y, (0, 1))):
        plus, minus = (faces[t] for t in pair)
        assert plus[0] == minus[0] == direction and plus[1] == conjugate(minus[1])
        z = unit(plus[1])
        assert 0 < z.imag < z.real and 0 < z.arg() < z.arg()+radius < arb.pi()/2
        scales.append(1/z.imag)
    assert scales[0] < arb("6.4611") and scales[1] < arb("7.0278")
    assert (2**-12*scales[0]).asin() < arb(".001578")
    assert (2**-12*scales[1]).asin() < arb(".001716")
    triple = (axis_x[1], axis_y[0], diagonal)
    assert [faces[t][0] for t in triple] == [(1, 0), (0, 1), (1, -1)]
    u, v, w = [closest_critical(faces[t][1]) for t in triple]
    tau = arb(".02096")
    assert u-v > aperture+2*tau and u-w > radius+2*tau and v+w < -radius-2*tau
    assert max(radius+abs(u), radius+abs(v), aperture+abs(w)) < arb.pi()/2
    assert 4*(arb(2)**-16/tau.sin()).asin() < arb(".002913")
    print("四读数的两钟误差放大率", *(x.str(14) for x in scales))
    print("三读数对应的驻点 x、y、x-y", *(x.str(14) for x in (u, v, w)))


def actual_pairs(faces):
    examples = {}
    for n in (8236, 8359, 148341):
        clocks = graph["phases"](acb.zeta_zero(n).imag)
        examples[n] = {face: (unit(z)*clocks[724]**a*clocks[804]**b).real
                       for face, ((a, b), z) in faces.items()}
    best_face = (1, 6, 8)
    for i, j, low, high in ((8236, 8359, ".05796101246", ".05796101248"),
                            (8236, 148341, ".00468848189", ".00468848190")):
        delta = {face: abs(examples[i][face]-examples[j][face]) for face in faces}
        assert delta[best_face] > arb(low) and all(x < arb(high) for x in delta.values())
        print("零点、单面最大读数差的包围", i, j, low, high)
    assert 2*arb(2)**-6 < arb(".05796101246") < arb(".05796101248") < 2*arb(2)**-5


if __name__ == "__main__":
    ctx.prec = 192
    data = enumerate_faces()
    reconstruction(data)
    actual_pairs(data)
