"""只用三种原味基态的有向概率恢复平均壳响应并给有限时长界。"""
from fractions import Fraction
from math import factorial, prod
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, arb_mat, ctx

torus = run_path(str(Path(__file__).with_name("ch39-exposed-flavor-response.py")))
study = torus["study"]
plane = study["flavor"]["basis"]
edges = ((0,1),(1,2),(2,0))
radial = (1,5,29)
c = [arb(v)/arb(867).sqrt() for v in (29,1,5)]
centers = list(map(arb, (".2","-1.009","-.205","0")))
radii = list(map(arb, (".052",".007",".009",".031")))


def fourier_range(rows):
    first = [acb.dft(row) for row in rows]
    modes = [acb.dft(col) for col in zip(*first)]
    n = len(rows)
    base = modes[0][0].real/n**2
    tail = sum(abs(z) for i,col in enumerate(modes) for j,z in enumerate(col) if i or j)
    return arb(base, (tail/n**2).upper())


def check_torus():
    _, rows = torus["response_grid"]()
    det0 = rows[0][0][0]
    eta,h,k,alpha = study["coefficients"]
    P,Q = [study["flavor"]["arcs"][d] for d in (32,152)]
    static = h*P+k*Q+alpha*study["flavor"]["identity"]
    values = []
    for row in rows:
        data = []
        for d, T in row:
            block = eta*T + static
            matrix = plane*block*plane.transpose()
            coords = [matrix[i,j].real for i,j in edges]+[block[0,1].imag]
            data.append([d/det0]+[d/det0*(v-u) for v,u in zip(coords,centers)])
        values.append(data)
    ranges = [fourier_range([[entry[k] for entry in row] for row in values])
              for k in range(5)]
    assert ranges[0] > arb(".80")
    assert all(abs(v)/ranges[0].lower() < w for v,w in zip(ranges[1:],radii))
    assert alpha-eta*arb(9)/4 > 0 and alpha+h+k < 9
    assert 2*eta*arb(9)/4 < 4
    print("Torus coordinate radii",*[v.str(6) for v in radii])


def jacobian(v):
    *d, y = v
    rows = [[2*d[i]*int(i==j) for j in range(3)]+[2*c[i]**2*y] for i in range(3)]
    last = []
    for i in range(3):
        j,k = (i+1)%3, (i+2)%3
        last.append(y*(c[j]*d[k]+c[k]*d[j]))
    last.append(sum(c[i]*d[(i+1)%3]*d[(i+2)%3] for i in range(3))-3*prod(c)*y*y)
    return arb_mat(rows+[last])


def max_row_norm(m):
    return max(sum(abs(m[i,j]).upper() for j in range(m.ncols())) for i in range(m.nrows()))


def identification(A, B):
    box = [arb(v,r.upper()) for v,r in zip(centers,radii)]
    J = jacobian(box)
    preconditioner = J.mid().inv().mid()
    identity = arb_mat([[int(i==j) for j in range(4)] for i in range(4)])
    error = max_row_norm(identity-preconditioner*J)
    assert error < arb(".42")
    lip = max_row_norm(preconditioner)/(1-error)
    # H p = 0 supplies the real diagonal entries from the three cyclic edges.
    norm = arb(1)
    for idx,(i,j) in enumerate(edges):
        H = acb_mat(3,3)
        H[i,j] = H[j,i] = 1
        H[i,i] = -arb(radial[j])/radial[i]
        H[j,j] = -arb(radial[i])/radial[j]
        T = plane.transpose()*H*plane
        norm += abs(A*(T[0,0]-T[1,1]).real/2+B*T[0,1].real)
        assert abs(plane[i,0]*plane[j,1]-plane[i,1]*plane[j,0]-c[idx]) < arb("1e-32")
    total = norm*lip/study["coefficients"][0]
    assert total < 33
    print("Jacobian residual and response-to-probability coefficient",error.str(8),total.str(8))


def signal(times):
    def exact(q):
        return arb(q.numerator)/q.denominator
    ws = [abs(prod(v/(v-u) for v in times if v!=u)) for u in times]
    n = len(times)
    bounds = []
    for k in (2,3):
        tail = (4*exact(prod(times))*9**(n+k-1)*(9*exact(max(times))).exp()
                /factorial(n+k-1))
        amplification = exact(sum(w/t**k for w,t in zip(ws,times)))
        bounds.append((arb("1.59808e-8")/33-tail)/amplification)
    return min(v.lower() for v in bounds)


def integer_example():
    previous = run_path(str(Path(__file__).with_name("ch39-six-duration-probability.py")))
    evolution = previous["ordered"]["evolution"]
    flavor = study["flavor"]
    def response(n):
        return study["native"]["plane_response"](acb.zeta_zero(n).imag)
    time = arb(3)/16
    wave = flavor["identity"]
    for n,tick in zip(previous["labels"],previous["ticks"]):
        wave = evolution(response(n),time*tick/10**12,study["coefficients"])*wave
    pure = evolution(response(11431),time,study["coefficients"])
    def probability(gate):
        full = flavor["radial"]+plane*gate*plane.transpose()
        return abs(full[1,0])**2
    diff = probability(wave)-probability(pure)
    assert diff > arb("3.32e-6")
    print("Six-segment clock, first to second flavor, t=3/16",diff.str(12))


def verify():
    ctx.prec = 192
    check_torus()
    A,B,_ = torus["polynomial"]()
    identification(A,B)
    times = [Fraction(x,4096) for x in (38,53,96,164,252,351,455,555,642,710,753,768)]
    finite = signal(times)
    single = signal([Fraction(1,2**40)])
    assert finite > arb("1.49e-17") and single > arb("3.17e-47")
    print("Native transition probability gaps, twelve durations and one",finite.str(9),single.str(9))
    integer_example()


if __name__ == "__main__":
    verify()
