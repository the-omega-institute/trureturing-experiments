"""Optimal invariant site sets and one-loop thresholds of shell gauge links."""
from itertools import combinations, permutations, product
from pathlib import Path
from runpy import run_path

from flint import arb, ctx, fmpz_mat

root = Path(__file__).parent
geometry = run_path(str(root/"ch39-intrinsic-shell-comparisons.py"))
squared = geometry["squared"]


def networks(points, group, vertices, cutoff):
    n = len(vertices)
    edges = [(i,j) for i,j in combinations(vertices,2) if 0<squared(points[i],points[j])<=cutoff]
    orbits = sorted({tuple(sorted({tuple(sorted((p[i],p[j]))) for p in group})) for i,j in edges})
    for bits in product((False,True),repeat=len(orbits)):
        chosen = [edge for active,orbit in zip(bits,orbits) if active for edge in orbit]
        if len(chosen)<n-1:
            continue
        lap = fmpz_mat(n,n)
        for i,j in chosen:
            a,b = vertices.index(i),vertices.index(j)
            lap[a,a] += 1
            lap[b,b] += 1
            lap[a,b] -= 1
            lap[b,a] -= 1
        trees = int(fmpz_mat([[lap[i,j] for j in range(n-1)] for i in range(n-1)]).det())
        if trees:
            yield n,len(chosen),trees,max(int(lap[i,i]) for i in range(n)),bits


def connected(points, vertices, cutoff):
    reached = {vertices[0]}
    frontier = list(reached)
    while frontier:
        i = frontier.pop()
        nearby = {j for j in vertices if 0<squared(points[i],points[j])<=cutoff} - reached
        reached.update(nearby)
        frontier.extend(nearby)
    return len(reached)==len(vertices)


def site_choices():
    points = sorted(p for p in product(range(30),repeat=3) if sum(x*x for x in p)==867)
    pole = (17,17,17)
    points.append(pole)
    group = [[19 if i==19 else points.index(tuple(p[a] for a in perm)) for i,p in enumerate(points)]
             for perm in permutations(range(3))]
    orbits = sorted({tuple(sorted({p[i] for p in group})) for i in range(20)})
    assert [len(o) for o in orbits]==[6,6,3,3,1,1]
    best = {}
    for n,cutoff in ((15,152),(16,152),(14,68)):
        active = []
        for mask in product((False,True),repeat=6):
            vertices = sorted(i for enabled,orbit in zip(mask,orbits) if enabled for i in orbit)
            if len(vertices)!=n:
                continue
            assert not connected(points,vertices,cutoff-1)
            if connected(points,vertices,cutoff):
                active.append((vertices,mask))
        assert len(active)==(2 if n==16 else 1)
        vertices,mask=active[-1]
        assert mask[:4]=={15:(1,1,1,0),16:(1,1,1,0),14:(0,1,1,1)}[n]
        best[n]=list(networks(points,group,vertices,cutoff))
    return [best[n] for n in (15,16,14)]


def mass_solution(graphs):
    rows = []
    for k,(n,e,trees,_,_) in zip((4,2,2),graphs):
        rows.append([arb(k)/n*z for z in (n-1-arb(21)/2*arb(n*trees).log(), -21*(n-1), 2*e-n+1)])
    (a,b,c),(d,e,f) = [[x-y for x,y in zip(row,rows[-1])] for row in rows[:2]]
    det = b*f-c*e
    assert abs(det)>0
    x,y = (c*d-a*f)/det,(a*e-b*d)/det
    common = rows[0][0]+rows[0][1]*x+rows[0][2]*y
    return x,y,common


def verify():
    ctx.prec = 160
    choices = site_choices()
    assert [len(group) for group in choices]==[2,16,4]
    sparse = [min(group,key=lambda g:(g[1],g[3])) for group in choices]
    assert [g[1:4] for g in sparse]==[(15,9,4),(15,1,5),(15,12,3)]
    x,y,_ = mass_solution(sparse)
    assert arb("1051.55") < (y-x).exp() < arb("1051.56")
    winners=[]
    for graphs in product(*choices):
        x,y,common = mass_solution(graphs)
        distance = abs(y-x)
        if distance < arb("1.084221"):
            assert arb("1.0842205")<distance and y>x
            assert [g[1:3] for g in graphs]==[(18,243),(18,27),(15,12)]
            winners.append((graphs,x,y,common))
        else:
            assert distance > arb("1.084221")
    assert len(winners)==6
    graphs,x,y,common = min(winners,key=lambda w:tuple(g[3] for g in w[0]))
    assert [g[3] for g in graphs]==[4,5,3]
    assert arb("2.95713") < (y-x).exp() < arb("2.95714")
    data=run_path(str(root/"ch39-published-gauge-inputs.py"))
    shared=data["model"]["target"](11431,data["q"]("26.81725683017028"))
    g=(shared/15+common/(48*arb.pi()**2))**(-arb(1)/2)
    assert arb("2.06889") < g < arb("2.06890")
    assert 4*g*g/(16*arb.pi()**2) < arb(".10843")
    print("128 invariant graph triples; optimal M_S/(gv)",(y-x).exp().str(12),
          "gv/Q",x.exp().str(12),"M_S/Q",y.exp().str(12),"site coupling",g.str(12))


if __name__=="__main__":
    verify()
