"""四个零点的确切取等分布与全谱方差的锐下界，固定中心味系数。"""
from pathlib import Path
from runpy import run_path

from flint import acb, arb, arb_mat, arb_poly

study = run_path(str(Path(__file__).with_name("ch39-gauge-flavor-variance.py")))
study["coefficients"][:] = [arb(k)/2048 for k in study["flavor"]["codes"]]
S, tau = study["first_lower"], study["tau"]
lo, hi = study["small_u"], study["large_u"]
support = (10762, 11432, 11831, 12105)


def row(zero):
    h = zero.imag
    x = study["center"]/(h*h+study["offset"])-1
    return x, study["data"](study["native"]["plane_response"](h))


def attaining_distribution(rows, constraints, targets):
    gaps, obs = zip(*rows)
    a,b,u,v,j = zip(*obs)
    size = len(rows)
    matrix = arb_mat([[1]*size, list(gaps), *constraints,
                      [ai-S*(1-ui) for ai,ui in zip(a,u)], list(u)])
    base = matrix.solve(arb_mat([[t] for t in (1,0,*targets,0,lo)]))
    speed = matrix.solve(arb_mat([[0]]*(size-1)+[[hi-lo]]))
    weights = [arb_poly([base[i,0],speed[i,0]]) for i in range(size)]
    U = arb_poly([lo,hi-lo])
    B,V,J = [sum(w*y for w,y in zip(weights,col)) for col in (b,v,j)]
    P = (1-S)*U*(1-U)-(1-U)*V+(S-(1-S)*U)*B
    phase = J-(tau/2)*P
    assert phase.degree() <= 2 and phase(0)*phase(1)<0
    assert 0 not in phase.derivative()(arb(0).union(1))
    left,right = arb(0),arb(1)
    for _ in range(100):
        mid = (left+right)/2
        value = phase(mid)
        assert 0 not in value
        if (value>0)==(phase(0)>0):
            left = mid
        else:
            right = mid
    t = left.union(right)
    probs = [w(t) for w in weights]
    assert all(w>0 for w in probs) and J(t)>0
    second = B(t)/(1-U(t))
    assert arb(".041775")**2 < second < arb(".041885")**2
    return probs, [p(t) for p in (U,B,V,J)]


def candidate():
    rows = [row(acb.zeta_zero(n)) for n in support]
    probs, reference = attaining_distribution(rows, [], [])
    minimum = sum(w*x**2 for w,(x,_) in zip(probs,rows))
    print("Unique candidate",*(w.str(18) for w in probs),"variance",minimum.str(20), flush=True)
    return rows, reference, minimum


def phase_tests(U, reference):
    uc,bc,vc,jc = reference
    cu = (1-S)*(1-U-uc)+vc-(1-S)*bc
    cb = S-(1-S)*U
    def tests(o):
        a,b,u,v,j = o
        phase = j-jc-tau/2*(cu*(u-uc)-(1-U)*(v-vc)+cb*(b-bc))
        return a-S*(1-u),phase
    return tests


def dual_at(U, reference, rows):
    tests = phase_tests(U, reference)
    m = arb_mat([[1,x,*tests(o)] for x,o in rows])
    c = m.solve(arb_mat([[x*x] for x,o in rows]))
    assert c[2,0]>0 and c[3,0]>0
    def linear(o):
        first, phase = tests(o)
        return c[2,0]*first + c[3,0]*phase
    return m.det(), c, linear


def compare_torus(phase_coefficient, linear):
    # 全环面的旧仿射函数小于 0.06。消去四重积系数后，用径向不动轴约束第二行至第三列的味跃迁概率。
    alpha = phase_coefficient/arb("10303.7")
    assert alpha > 0
    def difference(o):
        return linear(o)-alpha*study["dual"](o)
    constant = difference([arb(0)]*5)
    coefficients = [difference([arb(int(i==k)) for i in range(5)])-constant for k in range(5)]
    radial = study["flavor"]["native"]["p"]
    assert tuple(radial) == (1,5,29)
    # |V_cb| = |((V-I)e_3)_2| 至多为两倍的平面投影长度。
    bounds = [1,arb(104)/867,1,1]
    upper = alpha*arb(".06")+constant+abs(coefficients[4])
    upper += sum(max(arb(0),k.upper())*b for k,b in zip(coefficients,bounds))
    return upper


def verify():
    rows, reference, minimum = candidate()
    duals = [dual_at(U,reference,rows) for U in (lo,hi)]
    assert all(det > 0 for det,_,_ in duals)
    assert -hi + arb(".041775")**2*(1-hi) > 0
    margins = []
    for _,c,linear in duals:
        upper = compare_torus(c[3,0],linear)
        assert upper < arb(".095")
        slope = c[1,0]
        for n in (8999,15001):
            h = acb.zeta_zero(n).imag
            x = study["center"]/(h*h+study["offset"])-1
            assert (2*x > slope if n==8999 else 2*x < slope)
            assert x*x-slope*x > minimum+upper
    for start in range(9000,15001,256):
        count = min(256,15001-start)
        zeros = acb.zeta_zeros(start,count)
        assert len(zeros) == count
        for n,z in enumerate(zeros,start):
            if n==11431 or n in support:
                continue
            x,obs = row(z)
            for _,c,linear in duals:
                margins.append((x*x-minimum-c[1,0]*x-linear(obs)).lower())
    assert min(margins) > 0
    print("Endpoint duals, minimum strict finite slack", min(margins).str(12),
          "relative standard deviation", minimum.sqrt().str(17), flush=True)


if __name__ == "__main__":
    study["torus_bound"]()
    verify()
