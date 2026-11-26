"""五枚正零点排除只含实端口的舒尔选点。"""
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, arb_mat, ctx

port = run_path(str(Path(__file__).with_name("ch39-single-port-ground-selection.py")))
study = port["study"]
native = study["native"]
adjoint = port["adjoint"]
labels = (147,11432,11491,11559,11763)
order = 7


def sector(n):
    h = acb.zeta_zero(n).imag
    x = study["center"]/(h*h+study["offset"])-1
    phases = study["geo"]["phases"](h)
    shell,block = study["phases_to_block"](phases[724],phases[804])
    diagonal = acb_mat(2,2)
    for i,k in enumerate(native["portals"]):
        diagonal[i,i] = shell[k,native["origin"]]
    source = native["source"]*diagonal*native["columns"].transpose()
    return x,block,source


def taylor(sector,z):
    x,block,source = sector
    shifted = acb_mat(block)
    for j in range(shifted.nrows()):
        shifted[j,j] += z
    coeff = []
    column = source
    for k in range(order):
        column = shifted.solve(column)
        t = (-1)**(k+1)*adjoint(source)*column
        coeff.append([t[0,0].real,t[1,1].real,t[0,1].real,t[0,1].imag])
    return x,coeff


def weights_on_box(data,mid,radius):
    expansions = [taylor(s,mid) for s in data]
    _,ref = expansions[-1]
    scales = [expansions[0][0]]+[arb(1)]*4
    matrices = []
    for k in range(order):
        cols = []
        for (x,coeff),scale in zip(expansions,scales):
            head = [arb(1),x] if k==0 else [arb(0),arb(0)]
            cols.append([v/scale for v in head+[a-b for a,b in zip(coeff[k][:3],ref[k][:3])]])
        matrices.append(arb_mat(cols).transpose())
    inverse = matrices[0].mid().inv().mid()
    xi = arb(0,radius)
    powers = [arb(1)]
    for k in range(1,order):
        powers.append(powers[-1]*xi)
    conditioned = sum((p*(inverse*m) for p,m in zip(powers,matrices)),arb_mat(5,5))
    gap = arb("1.3")-4*arb.pi()/32+mid
    tail = 4*(radius/gap)**order/(gap-radius)
    for i in range(5):
        amp = sum(abs(inverse[i,k]) for k in (2,3,4))
        for j,scale in enumerate(scales):
            conditioned[i,j] += arb(0,(amp*tail/scale).upper())
    rhs = inverse*arb_mat([[1],[0],[0],[0],[0]])
    try:
        sol = conditioned.solve(rhs)
    except ZeroDivisionError:
        return None
    weights = [sol[j,0]/scales[j] for j in range(5)]
    imaginary = arb(0)
    for weight,(_,coeff) in zip(weights,expansions):
        entry = sum(p*(v[3]-u[3]) for p,v,u in zip(powers,coeff,ref))
        imaginary += weight*(entry+arb(0,tail.upper()))
    return weights,imaginary


def verify():
    ctx.prec = 192
    data = [sector(n) for n in (*labels,11431)]
    exact,_ = weights_on_box(data,arb(0),arb(0))
    assert all(w>0 for w in exact)
    stack = [(arb(1),arb(1))]
    count = 0
    floors = [arb(1)]*5
    while stack:
        mid,r = stack.pop()
        result = weights_on_box(data,mid,r)
        if result is None or not result[1] < -arb(".008") or any(not w>0 for w in result[0]):
            assert r > arb(2)**-14
            stack.extend([(mid-r/2,r/2),(mid+r/2,r/2)])
            continue
        ws,_ = result
        count += 1
        floors = [min(prev,w.lower()) for prev,w in zip(floors,ws)]
    print("Weights at zero energy",*[w.str(12) for w in exact])
    print("Closed interval [0,2], boxes and componentwise lower bounds",count,
          *[w.str(8) for w in floors])
    T = native["plane_response"](acb.zeta_zero(11431).imag)
    # Gershgorin gives a convenient sufficient range for every real unit port.
    radius = max(abs(T[j,j])+abs(T[j,1-j]) for j in range(2))
    assert 1/study["center"]+radius < 1


if __name__ == "__main__":
    verify()
