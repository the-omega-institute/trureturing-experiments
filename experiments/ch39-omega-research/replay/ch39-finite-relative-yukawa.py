"""Relative edge-mode mass readout, exact parity obstruction and feasibility threshold."""
from itertools import combinations, permutations

from sympy import I, Matrix, Rational as Q, eye, ones, zeros


def rotor(t):
    return (1-t*t)/(1+t*t), 2*t/(1+t*t)


def benchmark():
    a,b = rotor(Q(57,500))
    c,d = rotor(Q(107,5000))
    e,f = rotor(Q(377,200000))
    re,im = rotor(Q(16,25))
    phase = re+I*im
    V = (Matrix([[1,0,0],[0,c,d],[0,-d,c]])
         *Matrix([[e,0,f*phase.conjugate()],[0,1,0],[-f*phase,0,e]])
         *Matrix([[a,b,0],[-b,a,0],[0,0,1]]))
    H = Matrix([[1,1,1],[1,-1,-1],[-1,1,-1],[-1,-1,1]])/2
    u = [Q(v) for v in (".00000615",".00311",".8616")]
    d = [Q(v) for v in (".0000135",".000268",".01401")]
    assert (V.H*V-eye(3)).applyfunc(lambda x:x.expand()) == zeros(3)
    assert H.T*H == eye(3) and H*H.T == eye(4)-ones(4)/4
    return H, H*Matrix.diag(*u)*H.T, (H*V*Matrix.diag(*d)*V.T*H.T).applyfunc(lambda x:x.expand())


def norm_squared(z):
    return (z*z.conjugate()).expand()


def main():
    axes = set(range(4))
    pairs = list(combinations(range(4),2))
    comp = lambda p: tuple(sorted(axes-set(p)))
    roots = {p:eye(4)[:,p[0]]-eye(4)[:,p[1]] for p in pairs}
    forms = {p:v*v.T for p,v in roots.items()}
    groups = {p:[] for p in pairs}
    for point in permutations((0,1,3,4)):
        key = tuple(sorted((point.index(0),point.index(1))))
        inversions = sum(point[i]>point[j] for i,j in pairs)
        groups[key].append((-1)**inversions)
        assert set(comp(key)) == {point.index(3),point.index(4)}
    assert all(sorted(signs)==[-1,-1,1,1] for signs in groups.values())
    H, up, down = benchmark()
    for p in pairs:
        assert (H.T*(forms[p]+forms[comp(p)])*H).is_diagonal()
    for Y in (up,down):
        assert sum((-Y[p]*forms[p] for p in pairs), zeros(4)) == Y
    # At g=5 each independent pair mean must lie inside its unit disk.
    radius2 = Q(25,16*24)
    disk_data = []
    for sign in (1,-1):
        k = {p:-up[p]-sign*down[p].conjugate() for p in pairs}
        for p in pairs:
            if p>comp(p):
                continue
            ev,od = (k[p]+k[comp(p)])/2, (k[p]-k[comp(p)])/2
            e2,o2 = norm_squared(ev),norm_squared(od)
            cross = abs((ev.conjugate()*od+od.conjugate()*ev).expand())
            disk_data.append((e2,cross,o2))
    def feasible(beta):
        return all((radius2-e2)*beta**2-cross*beta-o2>=0 for e2,cross,o2 in disk_data)
    assert any(o2>0 for _,_,o2 in disk_data) and feasible(1)
    denominator = 10**10
    lo,hi = 0,denominator
    while hi-lo>1:
        mid = (hi+lo)//2
        if feasible(Q(mid,denominator)):
            hi = mid
        else:
            lo = mid
    assert not feasible(Q(lo,denominator)) and feasible(Q(hi,denominator))
    print("Exact fixed-benchmark threshold in [",lo,",",hi,"] /",denominator)
    single_g2 = 16*24*max(e2+cross+o2 for e2,cross,o2 in disk_data)
    lo,hi = 0,5*denominator
    while hi-lo>1:
        mid = (hi+lo)//2
        if Q(mid,denominator)**2 >= single_g2:
            hi = mid
        else:
            lo = mid
    assert Q(lo,denominator)**2<single_g2<=Q(hi,denominator)**2
    print("One observed minimal edge selects one cover; g threshold in [",lo,",",hi,"] /",denominator)
    print("Independent opposite parity pairs per group realize the two disk means.")


if __name__ == "__main__":
    main()
