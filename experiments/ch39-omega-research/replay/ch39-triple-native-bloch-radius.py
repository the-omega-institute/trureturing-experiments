"""三个相邻实际编号上，上型平凡二阶共同正块的旋转不变径向界。"""
from fractions import Fraction
from itertools import combinations
from pathlib import Path
from runpy import run_path

from flint import acb_mat, arb


TRIPLE_RATIO_FLOORS = (".00711281",".00449434")


def dot(a,b):
    return sum(x*y for x,y in zip(a,b,strict=True))


def cross(a,b):
    return tuple(a[(i+1)%3]*b[(i+2)%3]-a[(i+2)%3]*b[(i+1)%3] for i in range(3))


def outer_planes(rows, model):
    rays = [model["bloch"](r) for r in rows if r[0,0].real+r[1,1].real>0]
    assert len(rays) == 7
    planes = []
    for normal,_ in model["facets"](rows):
        integer = tuple(int((x*2**40).mid().floor().unique_fmpz()) for x in normal)
        bound = int(max(dot(integer,r).upper() for r in rays).ceil().unique_fmpz())
        assert all(dot(integer,r) < bound for r in rays)
        planes.append((integer,bound))
    return planes


def radius_squared(planes):
    # 加入单位立方体，使有理外包有界；任一顶点有三条线性独立的相等约束。
    planes = list(planes)
    for j in range(3):
        for sign in (-1,1):
            n = tuple(sign if i==j else 0 for i in range(3))
            planes.append((n,1))
    radii = []
    for (a,h),(b,k),(c,l) in combinations(planes,3):
        ab,bc,ca = cross(a,b),cross(b,c),cross(c,a)
        det = dot(a,bc)
        if det == 0:
            continue
        num = tuple(h*x+k*y+l*z for x,y,z in zip(bc,ca,ab,strict=True))
        if det < 0:
            det = -det
            num = tuple(-x for x in num)
        if all(dot(n,num) <= bound*det for n,bound in planes):
            radii.append(Fraction(dot(num,num),det**2))
    assert radii
    best = max(radii)
    return (arb(best.numerator)/best.denominator).sqrt()


def verify(root):
    model = run_path(str(root/"ch39-free-four-site-native-positive-spectrum.py"))
    complex_model = run_path(str(root/"ch39-complex-four-site-unit-quark-amplitudes.py"))
    _,weight = complex_model["verify"](root)
    center = 2_306_888_100_085
    blocks = {i:model["native"](root,center+i) for i in (-1,0,1)}
    targets = []
    for q,ww in enumerate(weight):
        raw = sum((p*w for p,w in zip(blocks[0][q],ww,strict=True)),acb_mat(4,4))
        target = model["bloch"](raw)
        targets.append(dot(target,target).sqrt())
    print("中心目标上、下平凡块的 Bloch 半径",*[r.str(15) for r in targets])
    assert targets[0] > arb(".99605")

    for q in (0,1):
        planes = {i:outer_planes(rows[q],model) for i,rows in blocks.items()}
        limits = ((".989440",".992967",".996527",".985875") if q==0 else
                  (".996815",".998590",".999086",".991052"))
        for numbers,limit in zip(((-1,0),(-1,1),(0,1),(-1,0,1)),limits,strict=True):
            radius = radius_squared([plane for i in numbers for plane in planes[i]])
            assert radius < arb(limit)
            ratio = (1-radius)/(1+radius)
            print("部门、编号偏移、Bloch 半径上界、特征值比下界",q,numbers,
                  radius.upper().str(14),ratio.lower().str(12))
            if len(numbers)==3:
                assert ratio>arb(TRIPLE_RATIO_FLOORS[q])
            if q==0 and -1 in numbers:
                assert radius < targets[0]


if __name__=="__main__":
    verify(Path(__file__).resolve().parent)
