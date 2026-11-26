"""Native charge-five triangle flux and independent positive endpoint profiles."""
from itertools import combinations, permutations, product
from pathlib import Path
from runpy import run_path
from flint import acb, acb_mat, arb, ctx

ctx.prec = 160
radius2 = 867
pole = (17,17,17)
outer = sorted(set(permutations((1,5,29))))
inner = sorted(set(permutations((13,13,23))))
points = outer+inner+sorted(v for v in product(range(30),repeat=3)
                            if sum(t*t for t in v)==radius2 and v not in outer+inner)
n = len(points)
assert n==19
laplacian = acb_mat(n,n)
transfer = {}


def dot(x,y):
    return sum(a*b for a,b in zip(x,y))


for i,j in combinations(range(n),2):
    x,y = points[i],points[j]
    d = 2*(radius2-dot(x,y))
    if d>804:
        continue
    det = sum(pole[k]*(x[(k+1)%3]*y[(k+2)%3]-x[(k+2)%3]*y[(k+1)%3]) for k in range(3))
    t = acb(radius2*(radius2+dot(pole,x)+dot(x,y)+dot(y,pole)),-arb(radius2).sqrt()*det)
    link = (t/abs(t)*(-1 if d==200 else 1))**5
    transfer[i,j] = link
    transfer[j,i] = link.conjugate()
    laplacian[i,i] += 1
    laplacian[j,j] += 1
    laplacian[i,j] = -link
    laplacian[j,i] = -link.conjugate()
assert len(transfer)==270
green = laplacian.inv()


def local_choice(triangle, word):
    # With a point first source, the two unturned entries on opposite sides
    # of the target ray supply a positive exact solution for the second.
    for site in triangle:
        for turn in (word,word.conjugate()):
            for i in range(n):
                reference = -turn*green[i,site]
                rays = [green[i,j]*reference.conjugate() for j in range(n) if j!=site]
                if all(v.real>0 for v in rays) and any(v.imag>0 for v in rays) and any(v.imag<0 for v in rays):
                    return True
    # Exclude all mixed first-source profiles as well, for either orientation.
    miss = -word if (-word).imag>0 else -word.conjugate()
    for site in triangle:
        for other in range(n):
            for j,k in ((site,other),(other,site)):
                for i in range(n):
                    for l in range(i,n):
                        area = miss*(green[i,j]*green[l,k].conjugate()+green[l,j]*green[i,k].conjugate())
                        assert area.imag>0
    return False


def verify():
    ref = green[6,7].real
    assert arb('.02208284702124')<ref<arb('.02208284702125')
    for i,j in product(range(n),repeat=2):
        assert abs(green[i,j].imag)<green[i,j].real/5
        if i in (6,7,8) and j in (6,7,8) and i!=j:
            continue  # a single orbit up to complex conjugation
        assert green[i,j].real>ref
    # Coordinates 0 and 1 can be interchanged. This reflects 0 <-> 2 and
    # 10 <-> 16, conjugating the parallel-transport matrix.
    value = (green[0,10]+green[2,10])/2
    ratio = value.imag/value.real
    assert arb('.12776994016813')<ratio<arb('.12776994016814')
    # A non-unit multiple of exp(i * maximum spread) suffices for signs.
    phase = value**2
    exceptional = {}
    for p in permutations(range(3)):
        move = lambda i: points.index(tuple(points[i][k] for k in p))
        odd = sum(p[i]>p[j] for i,j in combinations(range(3),2))%2
        j,k = (move(10),move(16)) if odd else (move(16),move(10))
        exceptional[j,k] = tuple(sorted((move(0),move(2))))
    assert len(exceptional)==3
    # For each ordered pair of columns, bound the signed area of their
    # positive mixtures. The three exceptional 2x2 blocks are positive
    # multiples of (a_i-a_l)^2 by reflection and equality at the midpoint.
    for j,k in product(range(n),repeat=2):
        block = exceptional.get((j,k),())
        for i in range(n):
            for l in range(i,n):
                coeff = (phase*(green[i,j]*green[l,k].conjugate()+
                                green[l,j]*green[i,k].conjugate())).imag
                if i!=l and i in block and l in block:
                    assert coeff<0
                else:
                    assert coeff>0
    cycle = (0,6,7,8,6,0)
    word = acb(1)
    for i,j in zip(cycle,cycle[1:]): word *= transfer[i,j]
    miss = -word
    assert miss.real>0 and miss.imag>0
    half = miss.imag/(1+miss.real)
    assert arb('.13196065051392')<half<arb('.13196065051393')
    assert half>ratio and word.real<0
    split = 1/(2*(1-word.real))
    assert 0<split and 2*split<1
    # The pole is invariant under coordinate permutations. Its three Green
    # entries against the inner orbit are positive real and equal.
    pole_index = points.index(pole)
    assert all(green[pole_index,j].real>0 for j in (6,7,8))
    # Overlapping continuous source paths cover every angular spread down
    # to zero. All tested column ratios have positive real part.
    u = green[5,6]/green[5,7]
    v = green[5,10]/green[5,16]
    assert 0<v.imag/v.real<u.imag/u.real<phase.imag/phase.real
    allowed = {}
    remote_only = []
    total = 0
    for i,j,k in combinations(range(n),3):
        cycle_edges = [(i,j),(j,k),(k,i)]
        if any(edge not in transfer for edge in cycle_edges):
            continue
        total += 1
        t = transfer[i,j]*transfer[j,k]*transfer[k,i]
        # |tan(delta/2)| <= ratio, delta measured from the negative axis.
        # Use the cosine comparison, which also handles t = 1 exactly.
        criterion = (1+ratio**2)*t.real+1-ratio**2
        assert criterion>0 or criterion<0
        if criterion<0:
            lengths = tuple(sorted(2*(radius2-dot(points[a],points[b])) for a,b in cycle_edges))
            allowed[lengths] = allowed.get(lengths,0)+1
            if not local_choice((i,j,k),t): remote_only.append(lengths)
    assert total==529
    assert allowed=={(200,244,724):6,(200,244,804):6,(52,172,200):6,
                     (12,200,292):6,(52,200,372):6,(68,68,200):3,
                     (52,52,200):3,(72,200,392):6,(56,56,200):3,(136,136,200):3}
    assert len(remote_only)==6
    print('Local one-return cancellation on 42 triangles. Transport required on',sorted(set(remote_only)))
    print(sum(allowed.values()),'of',total,'unoriented elementary triangles permit one-return cancellation.')
    print('Unphased minimum',ref.str(15),'maximum phase half-angle tangent',ratio.str(15))
    print('Native single-word miss half-angle tangent',half.str(15),'two-word weights',split.str(13))


def physical():
    # Use the independent archived eighth-zero fit for the absolute scale.
    root = Path(__file__).resolve().parent
    if not (root/'ch39-quaternion-low-masses.txt').is_file():
        run_path(str(root/'ch39-quaternion-scalar-levels.py'),run_name='__main__')
    mass = run_path(str(root/'ch39-su5-wilson-mass.py'))
    flux = run_path(str(root/'ch39-su5-shell-flux.py'))
    _,_,tau,inverse = mass['antiparallel_crossing']()
    internal = 1000*tau.exp()/mass['r']
    coeff = 4*arb.pi()/inverse.upper()*flux['volume'].lower()*green[6,7].real.lower()/internal.upper()**2
    rate = arb('6.582119569e-25')/(arb('2.4e34')*31557600)
    mp,mpi = arb('.938272'),arb('.134977')
    bound = (32*arb.pi()*rate/(mp*(1-(mpi/mp)**2)**2)).sqrt()/coeff
    print('Unphased geometry coefficient floor GeV^-2',coeff.lower().str(12),
          'effective hadronic factor upper GeV^2',bound.upper().str(12))
    assert coeff>arb('9.71e-32') and bound<arb('.1015')


if __name__=='__main__':
    import sys
    verify()
    if '--physical' in sys.argv: physical()
