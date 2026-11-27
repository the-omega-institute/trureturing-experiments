"""由第三枢轴的末端口和三个边界坐标，以整数点积构造站集。"""
from fractions import Fraction
from itertools import combinations, permutations
from math import isqrt
from pathlib import Path
from runpy import run_path


def dot(a, b):
    return sum(x*y for x, y in zip(a, b, strict=True))


def sub(a, b):
    return tuple(x-y for x, y in zip(a, b, strict=True))


def select(points, pole, boundary):
    centroid = tuple(Fraction(sum(column), 3) for column in zip(*boundary, strict=True))
    toward_pole = sub(pole, centroid)
    edges = [sub(a, b) for a, b in combinations(boundary, 2)]
    side2, = {dot(v, v) for v in edges}
    assert side2 > 0 and all(dot(toward_pole, v) == 0 for v in edges)
    # 三直角顶点距三角平面为边长除以 sqrt(6)，取极点对侧。
    ratio2 = side2/(6*dot(toward_pole, toward_pole))
    ratio = Fraction(isqrt(ratio2.numerator), isqrt(ratio2.denominator))
    assert ratio > 0 and ratio*ratio == ratio2
    center = tuple(c-ratio*v for c,v in zip(centroid, toward_pole, strict=True))
    vectors = [sub(b, center) for b in boundary]
    assert all(dot(a, b) == 0 for a, b in combinations(vectors, 2))
    normal = tuple(v-m for v, m in zip(boundary[-1], center, strict=True))
    length2 = sum(v*v for v in normal)
    assert length2 > 0
    scores = [sum(v*(x-m) for v, x, m in zip(normal, p, center, strict=True))**2 for p in points]
    support = tuple(i for i, score in enumerate(scores) if score >= length2**2)
    return support, center, normal, scores, length2


def verify(root):
    graph = run_path(str(root / "ch39-three-ray-chiral-boundary.py"))
    points = graph["geo"]["sites"]
    ports = (graph["soft"]["pole"], *graph["boundary_sites"]())
    assert ports == (14, 8, 7, 6)
    boundary = [points[i] for i in ports[1:]]
    pole = points[ports[0]]
    support, center, normal, scores, length2 = select(points, pole, boundary)
    assert center == (13, 13, 13) and normal == (0, 0, 10)
    assert support == (0, 2, 3, 5, 6, 9, 11, 13)
    assert min(scores[i] for i in support) == 100*length2
    assert max(scores[i] for i in range(len(points)) if i not in support) == 64*length2
    for order in permutations(range(3)):
        def permute(point):
            return tuple(point[j] for j in order)
        assert select(list(map(permute, points)), permute(pole), list(map(permute, boundary)))[0] == support
    backwards = select(points[::-1], pole, boundary)[0]
    assert tuple(len(points)-1-i for i in reversed(backwards)) == support

    # 独立从原始 Gram 区间核对末端口是第三源三个边界分量中唯一最轻者。
    from flint import ctx
    ctx.prec = 256
    gram = run_path(str(root / "ch39-native-rh-ray-pvm.py"))["gram_cholesky"](root)
    born = [x**2 for x in gram[3]]
    assert all(x.is_finite() for x in born) and born[0] > born[1] > born[2] > born[3] > 0
    print("末端口、中心、法向、所选站位", ports[-1], tuple(map(int, center)), tuple(map(int, normal)), support)
    print("归一化平方投影，支撑内最低、补集最高", 100, 64)
    print("原始第三源三条边界 Born 权", *(x.str(10) for x in born[1:]))
    return support


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)
