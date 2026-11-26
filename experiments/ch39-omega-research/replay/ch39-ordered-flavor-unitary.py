"""同一规范均值与端口六点权重下，有限有序味演化的单一幺正矩阵。"""
from itertools import permutations
from pathlib import Path
from runpy import run_path

from flint import acb, arb

previous = run_path(str(Path(__file__).with_name("ch39-zero-port-phase-branch.py")))
flavor = previous["box"]["flavor"]
study = previous["study"]
ids = (11431, 11304, 11372, 11415, 11430, 11432)


def evolution(response, duration, coefficients):
    eta, h32, h152, alpha = coefficients
    ident = flavor["identity"]
    action = eta*response + h32*flavor["arcs"][32] + h152*flavor["arcs"][152] + alpha*ident
    mean = action.trace().real/2
    reduced = action-mean*ident
    gap = (reduced[0,0].real**2 + abs(reduced[0,1])**2).sqrt()
    assert gap > 0
    angle = duration*gap
    return acb(0,duration*mean).exp()*(angle.cos()*ident + acb(0,angle.sin()/gap)*reduced)


def readings(steps, order):
    wave = flavor["identity"]
    for n in order:
        wave = wave*steps[n]
    basis = flavor["basis"]
    matrix = flavor["radial"] + basis*wave*basis.transpose()
    a,b,u = [abs(matrix[i,j])**2 for i,j in ((0,1),(1,2),(0,2))]
    quartet = matrix[0,1]*matrix[1,2]*matrix[0,2].conjugate()*matrix[1,1].conjugate()
    phase = (quartet + a*b*u/(1-u)).arg()
    return [(a/(1-u)).sqrt(),(b/(1-u)).sqrt(),u.sqrt(),phase]


def verify():
    cost,strength,*_ = previous["setup"]()
    heights = [acb.zeta_zero(n).imag for n in ids]
    xs = [study["center"]/(h*h+study["offset"])-1 for h in heights]
    xs[0] = arb(0)
    r = arb(1)/3
    fixed = [arb(n)/10**9 for n in (203009,710202,6405509)]
    remainder = r-sum(fixed)
    left = (-remainder*xs[-1]-sum(x*w for x,w in zip(xs[1:4],fixed)))/(xs[-2]-xs[-1])
    weights = [1-r,*fixed,left,remainder-left]
    assert all(w>0 for w in weights)
    price = strength*sum(w*cost(x)[0] for w,x in zip(weights,xs))
    assert price < arb("1.542272e-29") < arb("1.5432819e-29")
    inputs = [study["native"]["plane_response"](h) for h in heights]
    centers = previous["box"]["centers"]
    cube = [c+arb(0,arb(2)**-28) for c in centers]
    thresholds = [arb(".000055")]*3 + [arb(".006")]
    example = (11372,11304,11432,11431,11430,11415)
    for coefficients in (centers,cube):
        steps = {n:evolution(t,w,coefficients) for n,t,w in zip(ids,inputs,weights)}
        accepted = 0
        for order in permutations(ids):
            quartet = [n for n in order if n not in (11304,11372)]
            predicted = quartet[0]==11432 or quartet[1:3]==[11432,11431]
            values = readings(steps,order)
            errors = [abs(x-y) for x,y in zip(values,flavor["shell"]["centers"])]
            if predicted:
                assert all(e<cut for e,cut in zip(errors,thresholds)),order
                accepted += 1
            else:
                assert any(e>cut for e,cut in zip(errors,thresholds)),order
        assert accepted == 240
        print("Example CKM",*[x.str(15) for x in readings(steps,example)])
    print("Accepted orders and port deviation",accepted,price.str(19))


if __name__ == "__main__":
    verify()
