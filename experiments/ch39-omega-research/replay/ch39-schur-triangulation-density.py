"""Certified RH density for degree three stable under all diagonal choices, at zero energy."""
from functools import lru_cache
from itertools import product
from pathlib import Path
import runpy
from flint import acb, arb

model = runpy.run_path(str(Path(__file__).with_name("ch39-shell-schur-flux.py")))
coefficients = model["schur_data"](0)
bond = lambda i,j: model["bond"](coefficients,i,j)
x = bond(7,6)
y = bond(11,12)[0,0]
pi = arb.pi()


def loop(vertices):
    poly = {(0,0): acb(1)}
    for i,j in zip(vertices,vertices[1:]+vertices[:1]):
        new = {}
        for (a,b), v in poly.items():
            for (c,d),w in bond(i,j).items():
                if w == 0:
                    continue
                key = a+c,b+d
                new[key] = new.get(key,acb(0))+v*w
        poly = new
    return poly


def terms(poly, factor=1):
    constant = (poly[0,0]*factor).imag
    result = []
    for a,b in sorted(set(poly) | {(-a,-b) for a,b in poly}):
        if a < 0 or a == 0 and b <= 0:
            continue
        plus = factor*poly.get((a,b),acb(0))
        minus = factor*poly.get((-a,-b),acb(0))
        result.append((a,b,(plus+minus).imag,(plus-minus).real))
    return constant, result


def center_radius(poly):
    t, harmonics = terms(poly, acb(0,1))
    return t - sum((c**2+s**2).sqrt() for a,b,c,s in harmonics)


# The triangle along the outer ring has principal angle in (-pi/2,pi/2).
for f in [(6,12,11),(7,12,11)]:
    assert center_radius(loop(f)) > 0
split = [loop((6,7,12)),loop((6,7,11))]
tests = [terms(x),terms(x,-y)] + [terms(v) for v in split]


@lru_cache(maxsize=131072)
def trig(level,center,width):
    scale = 1 << level
    angle = pi*(center%(2*scale)+arb(0,width))/scale
    return angle.cos(),angle.sin()


def test_value(test,level,i,j):
    value,harmonics = test
    for a,b,c,s in harmonics:
        cosine,sine = trig(level,(2*i+1)*a+(2*j+1)*b,abs(a)+abs(b))
        value += c*cosine + s*sine
    return value


def density(depth):
    valid = excluded = 0
    todo = [(0,0,0)]
    while todo:
        level,i,j = todo.pop()
        good = True
        for test in tests:
            value=test_value(test,level,i,j)
            if value < 0:
                excluded += 4**(depth-level)
                break
            if not value > 0:
                good = False
        else:
            if good:
                valid += 4**(depth-level)
            elif level < depth:
                todo.extend((level+1,2*i+a,2*j+b) for a,b in product((0,1),repeat=2))
    return arb(valid)/4**depth, 1-arb(excluded)/4**depth


def actual_degrees(n):
    g = model["geometry"]
    h = acb.zeta_zero(n).imag
    p,q = (model["shell"][key]*h for key in ("a","b"))
    def face_angle(f):
        value=model["evaluate"](loop(f),p,q)
        assert not value.imag.contains(0) or value.real>0
        return value.arg()
    triangles=sum(face_angle(f) for f in g["big"] if len(f)==3)
    quad=[f for f in g["big"] if len(f)==4]
    m=int(((triangles+sum(face_angle(f) for f in quad))/(2*pi)).unique_fmpz())
    def split_angle(f,rotate):
        if rotate:
            f=f[1:]+f[:1]
        return face_angle((f[0],f[1],f[2]))+face_angle((f[0],f[2],f[3]))
    choices=[[split_angle(f,rotate) for rotate in (0,1)] for f in quad]
    deg={int(((triangles+sum(seq))/(2*pi)).unique_fmpz()) for seq in product(*choices)}
    return m,sorted(deg)


if __name__ == "__main__":
    low,high=density(14)
    assert 0 < low < high < arb(".005")
    print("density of degree three on polygons and all eight triangulations",
          low.str(15),high.str(15))

    for n, result in [(1,(0,[0])), (5,(3,[0,1,2,3])), (6,(3,[0,1,2,3])),
                      (35,(0,[0,1,2,3])), (68,(3,[3]))]:
        assert actual_degrees(n) == result
