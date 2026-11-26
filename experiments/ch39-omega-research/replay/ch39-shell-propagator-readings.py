"""Cyclic observability of charge-five endpoint currents at RH-selected momenta."""
from functools import lru_cache
from pathlib import Path
from runpy import run_path
import json
from flint import acb, acb_mat, arb, arb_mat, ctx

root = Path(__file__).resolve().parent
native = run_path(str(root/'ch39-shell-current-cyclic-protection.py'))
points = native['points']
groups = []
for i, v in enumerate(points):
    if any(i in orbit for orbit in groups):
        continue
    groups.append(sorted({points.index(v[k:]+v[:k]) for k in range(3)}))
pole = next(j for j,o in enumerate(groups) if len(o)==1)
source = acb_mat([[int(j==pole)] for j in range(7)])
ctx.prec = 768
prefix = json.loads((root/'ch39-tick-separators.json').read_text())['prefix']

@lru_cache(None)
def height(i):
    h = prefix[i-1][0]
    return arb(h).union(arb(h+1))/2**400

@lru_cache(None)
def row_edges():
    rows = []
    for j, orbit in enumerate(groups):
        i = orbit[0]
        for k, other in enumerate(groups):
            for m in other:
                if (i,m) not in native['links']:
                    continue
                a,b = native['link'](i,m)
                phase = acb(a,arb(3).sqrt()*b)/arb(a*a+3*b*b).sqrt()
                length = sum((u-v)**2 for u,v in zip(points[i],points[m]))
                rows.append((j,k,length,phase))
    return rows

@lru_cache(None)
def density(length):
    R = native['radius2']
    angle = (1-arb(length)/(2*R)).acos()
    top = (1-arb(804)/(2*R)).acos()
    return 104*arb.pi()**2*angle/(100*460**2*top)

def pencil(s,continuous):
    matrix = acb_mat([[s if j==k else 0 for k in range(7)] for j in range(7)])
    terms = {}
    for j,k,length,phase in row_edges():
        if length not in terms:
            if continuous and s!=0:
                h = (s*density(length)).sqrt()
                terms[length] = (h/h.tanh(),h/h.sinh())
            else:
                terms[length] = (arb(1),arb(1))
        diag,off = terms[length]
        matrix[j,j] += diag
        matrix[j,k] -= off*phase
    return matrix

def sample(s,continuous):
    solution = pencil(s,continuous).solve(source)
    return [solution[j,0].conjugate() for j in range(7)]

def real_row(row):
    # Conjugate pairs of C3-orbits are (0,1), (3,4).
    return [row[j].real for j in (0,2,3,5,6)]+[row[j].imag for j in (0,3)]

def check_samples(shifts,continuous):
    rows = [sample(s,continuous) for s in shifts]
    samples = acb_mat(rows)
    det = abs(samples.det())
    assert det>0
    real_rows = [real_row(row) for row in rows]
    assert abs(arb_mat([r[:5] for r in real_rows[:5]]).det())>0
    assert abs(arb_mat([r[5:] for r in real_rows[:2]]).det())>0
    static = acb_mat([sample(arb(0),continuous)]).transpose()
    coeff = samples.transpose().solve(static)
    factor = sum(abs(coeff[i,0]) for i in range(7))
    recovery = samples.inv()
    condition = max(sum(abs(recovery[i,j]) for j in range(7)) for i in range(7))
    return factor,condition

def main():
    assert groups==[[0,3,4],[1,2,5],[6,7,8],[9,15,16],[10,13,17],[11,12,18],[14]]
    static = pencil(arb(0),False)
    columns = []
    v = source
    for _ in groups:
        columns.append([v[j,0] for j in range(7)])
        v = static*v
    cyclic = abs(acb_mat(columns).det())
    assert arb(17638680)<cyclic<arb(17638681)
    for continuous in (False,True):
        for name,shifts in [('inverse',[1/height(i)**2 for i in range(1,8)]),
                            ('scaled',[(height(i)/height(1))**2 for i in range(1,8)])]:
            sf,cf=check_samples(shifts,continuous)
            lower,upper = {(False,'inverse'):('618.75632','618.75633'),
                           (True,'inverse'):('618.75642','618.75643'),
                           (False,'scaled'):('1497.74159','1497.74160'),
                           (True,'scaled'):('1498.25878','1498.25879')}[continuous,name]
            assert arb(lower)<sf<arb(upper)
            print(continuous,name,'static recovery',sf.str(15),'orbit recovery',cf.str(15))

if __name__=='__main__': main()
