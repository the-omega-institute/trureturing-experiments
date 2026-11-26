"""在整个双弧相位环面上给中心的无迹球壳响应构造支撑平面。"""
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, arb_mat, ctx

study = run_path(str(Path(__file__).with_name("ch39-gauge-flavor-variance.py")))
native = study["native"]
size = 48
degree = 19
curvature = arb(1)/1000


def response_grid(shift=None):
    phases = study["geo"]["phases"](acb.zeta_zero(11431).imag)
    nodes = [acb(0,2*arb.pi()*k/size).exp() for k in range(size)]
    inside, portals, columns = (native[k] for k in ("inside","portals","columns"))
    origin = native["origin"]
    rows = []
    for a in nodes:
        row = []
        for b in nodes:
            m,block = study["phases_to_block"](a*phases[724],b*phases[804])
            if shift is not None:
                for k in range(block.nrows()):
                    block[k,k] += shift
            det = block.det().real
            assert det > 0
            solution = block.solve(native["source"])
            exchange = acb_mat([[-m[origin,i]*solution[inside.index(i),k]*m[j,origin]
                                 for k,j in enumerate(portals)] for i in portals])
            t = columns*exchange*columns.transpose()
            row.append((det,t))
        rows.append(row)
    return nodes,rows


def samples():
    nodes,rows = response_grid()
    def xyz(t):
        return ((t[0,0].real-t[1,1].real)/2, t[0,1].real, t[0,1].imag)
    det0,t0 = rows[0][0]
    xyz0 = xyz(t0)
    return [[[(v-w)*d/det0 for v,w in zip(xyz(t),xyz0)]
             +[d/det0*(4-2*nodes[i].real-2*nodes[j].real)]
             for j,(d,t) in enumerate(row)] for i,row in enumerate(rows)]


def fourier(series):
    return [v/len(series) for v in acb.dft(series)]


def derivative(series):
    coeff = fourier(series)
    return sum((acb(0,k)*coeff[k%size]).real for k in range(-degree,degree+1))


def polynomial():
    values = samples()
    # det(B) has degree at most 18. Each entry of adj(B) has degree
    # at most 17. The two external phase factors increase this to 19.
    gradients = [[derivative([row[j] for row in axis]) for j in range(3)]
                 for axis in (values[0],[row[0] for row in values])]
    system = arb_mat([row[:2] for row in gradients])
    a,b = system.solve(arb_mat([[row[2]] for row in gradients])).entries()
    print("Supporting coefficients",a.str(18),b.str(18),flush=True)
    scores = [[a*x+b*y-z-curvature*d for x,y,z,d in row] for row in values]
    return a,b,torus_coefficients(scores)


def torus_coefficients(scores):
    first = [fourier(row) for row in scores]
    full = [fourier([row[j%size] for row in first]) for j in range(-degree,degree+1)]
    coeff = {(i,j):full[j+degree][i%size] for i in range(-degree,degree+1) for j in range(-degree,degree+1)}
    return coeff


def certify(coeff, *, stationary_origin=True):
    # 递归方格上的 Taylor 下界；原点附近直接检验 Hessian。
    active = {ij:v for ij,v in coeff.items() if abs(v).upper()>arb("1e-13")}
    errors = {}
    for p,q in ((0,0),(1,0),(0,1),(2,0),(1,1),(0,2),(3,0),(2,1),(1,2),(0,3)):
        errors[p,q] = sum(abs(v).upper()*abs(i)**p*abs(j)**q for (i,j),v in coeff.items() if (i,j) not in active)
    def global_bound(p,q):
        return errors[p,q]+sum(abs(v).upper()*abs(i)**p*abs(j)**q for (i,j),v in active.items())
    m20,m11,m02 = [global_bound(*ij) for ij in ((2,0),(1,1),(0,2))]
    radius = arb(0)
    if stationary_origin:
        third = [global_bound(p,3-p) for p in range(4)]
        h0,h1,h2 = [sum((v*acb(0,i)**p*acb(0,j)**q).real for (i,j),v in coeff.items())
                    for p,q in ((2,0),(1,1),(0,2))]
        radius = arb(1)/8
        while True:
            left = h0-radius*(third[3]+third[2])
            cross = abs(h1)+radius*(third[2]+third[1])
            right = h2-radius*(third[1]+third[0])
            if left>0 and right>0 and left*right>cross**2:
                break
            radius /= 2
            assert radius>arb(2)**-18
        print("Local strictly convex square half-width",radius.str(9),"active modes",len(active),flush=True)
    stack = [(0,0,1)]
    visited = 0
    zero_approx = arb(0)
    while stack:
        i,j,level = stack.pop()
        # 方格中心为 (i,j)*pi/level，半边长为 pi/level。
        width = arb.pi()/level
        u,v = i*width,j*width
        if abs(u)+width < radius and abs(v)+width < radius:
            continue
        zu,zv = acb(0,u).exp(),acb(0,v).exp()
        powers_u = {k:zu**k for k in {m for m,n in active}}
        powers_v = {k:zv**k for k in {n for m,n in active}}
        value,du,dv = zero_approx,zero_approx,zero_approx
        for (m,n),c in active.items():
            term = c*powers_u[m]*powers_v[n]
            value += term.real
            du -= m*term.imag
            dv -= n*term.imag
        lower = (value-errors[0,0]-width*(abs(du)+abs(dv)+errors[1,0]+errors[0,1])
                 -width**2*(m20+2*m11+m02)/2)
        visited += 1
        if lower>0:
            continue
        assert level < 2**18,(i,j,level)
        stack.extend((2*i+p,2*j+q,2*level) for p in (-1,1) for q in (-1,1))
    print("Certified torus cells outside half-width",radius.str(9),"visited",visited)


if __name__ == "__main__":
    ctx.prec = 192
    study["uniform_strip"]()
    run_path(str(Path(__file__).with_name("ch39-zero-haar-clocks.py")))
    before, center, after = [z.imag for z in acb.zeta_zeros(11430, 3)]
    assert before < center < after
    a,b,coeff = polynomial()
    certify(coeff)
