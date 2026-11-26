"""规范均值固定时，四项味窗对无中心零点混合的相对方差约束。"""
from math import factorial, prod
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, ctx

ctx.prec = 160
root = Path(__file__).parent
flavor = run_path(str(root / "ch39-fixed-local-flavor.py"))
gauge = run_path(str(root / "ch39-high-zero-planck-window.py"))
native = flavor["native"]
geo = native["geo"]
coefficients = [arb(k)/2048 + arb(0, arb(2)**-28) for k in flavor["codes"]]
offset = gauge["q"](gauge["OFFSET"])
center = acb.zeta_zero(11431).imag**2 + offset
tau = arb("1.141").tan()
small_u = arb(".003677")**2
large_u = arb(".003787")**2
first_lower = arb(".224955")**2
u0 = arb(".003732")**2
a0, b0 = [(1-u0)*arb(v)**2 for v in (".22501", ".04183")]
da, db, du = arb(".000026"), arb(".0000047"), arb(".000000414")


def data(response):
    v = flavor["mixing"](response, coefficients)
    j = (v[0,1]*v[1,2]*v[0,2].conjugate()*v[1,1].conjugate()).imag
    return [abs(v[i,k])**2 for i,k in ((0,1),(1,2),(0,2),(2,0))] + [j]


def dual(values):
    a, b, u, v, j = values
    lower = ((1-2*u0)*u + u0**2-du**2 - (1-small_u)*v
             -(a0+b0)*u-u0*(a+b)+(a0+b0)*u0-du*(da+db)
             +(1+2*small_u)*(b0*a+a0*b-a0*b0-da*db))
    return arb("53.798")*(a-first_lower*(1-u)) + arb("10303.7")*(j-tau*lower/2)


def phases_to_block(z, w):
    m = geo["matrix"](0, {724:z, 804:w})
    inside = native["inside"]
    return m, acb_mat([[m[i,j] for j in inside] for i in inside])


def grid_score(z, w):
    m, block = phases_to_block(z, w)
    portals, inside = native["portals"], native["inside"]
    x = block.solve(native["source"])
    small = acb_mat([[-m[0,i]*x[inside.index(i),k]*m[j,0] for k,j in enumerate(portals)]
                     for i in portals])
    columns = native["columns"]
    return dual(data(columns*small*columns.transpose()))


def uniform_strip():
    nodes = [acb(0, 2*arb.pi()*i/32).exp() for i in range(32)]
    for z in nodes:
        for w in nodes:
            _, block = phases_to_block(z, w)
            shifted = block - arb("1.3")*acb_mat([[int(i==j) for j in range(18)] for i in range(18)])
            for k in range(1,19):
                minor = acb_mat([[shifted[i,j] for j in range(k)] for i in range(k)])
                assert minor.det().real > 0
            columns = block.solve(native["source"])
            assert sum(abs(c)**2 for c in columns.entries()) < arb(".41")**2
    # 实角最近网点的内部块扰动。每种相位的边矩阵每行、列至多有两项。
    perturb = 4*arb.pi()/32
    gap = arb("1.3") - perturb
    projected = arb(".41")/(1-perturb/arb("1.3"))
    radius = arb(1)/4
    exp_r = radius.exp()
    delta = 2*(exp_r-1)
    assert gap > delta
    compressed_change = projected**2*delta / (1-delta/gap)
    # 两个门点的相位与其逆分别从左、右乘。两列单位法线的算子范数平方至多二。
    response_change = 2*((exp_r**2-1)*projected + exp_r**2*compressed_change)
    matrix_bound = (abs(coefficients[0]).upper()*response_change).exp()
    # 复角上的概率及其反射解析延拓至多 B²，四元积虚部至多 B⁴。
    q = matrix_bound**2
    # dual 对传入的五个读数仿射，基向量之差给出其精确线性系数。
    constant = dual([arb(0)]*5)
    linear = [dual([arb(int(i==k)) for i in range(5)])-constant for k in range(5)]
    bound = abs(constant) + q*sum(abs(x) for x in linear[:4]) + abs(linear[4])*q**2
    print("Strip upper bound", bound.upper().str(10), flush=True)
    return bound.upper()


def torus_bound():
    bound = uniform_strip()
    size, order = 384, 12
    nodes = [acb(0,2*arb.pi()*i/size).exp() for i in range(size)]
    scores = [[grid_score(z,w) for w in nodes] for z in nodes]
    highest = max(x.upper() for row in scores for x in row)
    second = max(abs(scores[(i+di)%size][(j+dj)%size]-2*scores[i][j]
                     + scores[(i-di)%size][(j-dj)%size]).upper()
                 for i in range(size) for j in range(size) for di,dj in ((0,1),(1,0)))
    positions = list(range(1-order, order+1))
    norm, curvature_weight = arb(0), arb(0)
    for i in range(128):
        x = (arb(i)/128).union(arb(i+1)/128)
        weights = [abs(prod((x-k)/(j-k) for k in positions if k!=j)) for j in positions]
        norm = max(norm, sum(weights).upper())
        curvature_weight = max(curvature_weight,
                               sum(abs(j*(j-1))*w/2 for j,w in zip(positions,weights)).upper())
    remainder = bound*factorial(order)**2*(8*arb.pi()/size)**(2*order)
    upper = highest + (1+norm)*(curvature_weight*second + remainder)
    print("Grid max, second difference, torus upper", highest.str(10), second.str(10),
          upper.upper().str(10), flush=True)
    assert upper < arb(".06")


def spectrum_bound():
    assert abs(small_u.union(large_u)-u0) < du
    for ref, radius, lo, hi in ((a0,da,".224955",".225065"),(b0,db,".041775",".041885")):
        interval = (1-small_u.union(large_u))*arb(lo).union(arb(hi))**2
        assert abs(interval-ref) < radius
    def gap(height):
        return center/(height**2+offset)-1
    for n in (8999,14001):
        x = gap(acb.zeta_zero(n).imag)
        assert (x > 0 if n==8999 else x < -arb(".01"))
        assert x*x + arb(".01018")*x > arb(".06681")
    central_score = dual(data(native["plane_response"](acb.zeta_zero(11431).imag)))
    assert central_score < arb(".006985")
    floor = None
    for start in range(9000,14001,256):
        count = min(256,14001-start)
        zeros = acb.zeta_zeros(start,count)
        assert len(zeros) == count
        for n,z in enumerate(zeros,start):
            if n == 11431:
                continue
            x = gap(z.imag)
            score = dual(data(native["plane_response"](z.imag)))
            lower = (x*x+arb(".01018")*x-score).lower()
            floor = lower if floor is None else min(floor,lower)
    print("Inner band, variance lower at zero center weight", floor.str(12), flush=True)
    return floor


def witness():
    labels = [10762, 11432, 11831, 12105]
    h = [acb.zeta_zero(n).imag for n in labels]
    x = [center/(t*t+offset)-1 for t in h]
    # 固定两个有理权重，再由总重与均值解出另两项。
    p, q = arb("0.31087432"), arb("0.05404096")
    left = -(p*x[1]+q*x[2]+(1-p-q)*x[3])/(x[0]-x[3])
    weights = [left,p,q,1-p-q-left]
    assert all(w>0 for w in weights)
    obs = [data(native["plane_response"](t)) for t in h]
    a,b,u,v,j = [sum(w*y for w,y in zip(weights,col)) for col in zip(*obs)]
    ss, tt = a/(1-u), b/(1-u)
    real = ((1-u)*(u-v)-u*(a+b)+(1+u)/(1-u)*a*b)/2
    outputs = [ss.sqrt(), tt.sqrt(), u.sqrt(), acb(real,j).arg()]
    radii = [flavor["shell"]["sine_radius"]]*3 + [flavor["shell"]["phase_radius"]]
    error = max((abs(v-c)/r).upper() for v,c,r in zip(outputs,flavor["shell"]["centers"],radii))
    variance = sum(w*g*g for w,g in zip(weights,x))
    print("Weights", *(w.str(12) for w in weights), "max flavor fraction",error.str(12),
          "relative variance", variance.str(12), flush=True)
    assert error < 1 and variance < arb(".006839892")


if __name__ == "__main__":
    witness()
    torus_bound()
    assert spectrum_bound() > arb(".0068")
