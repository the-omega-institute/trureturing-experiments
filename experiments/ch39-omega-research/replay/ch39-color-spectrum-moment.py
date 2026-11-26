"""Third static spectral moment of a fundamental color probe on the same original edges."""
from itertools import combinations
from pathlib import Path
import runpy
from flint import acb_mat

root = Path(__file__).resolve().parent
color = runpy.run_path(str(root/'ch39-eighth-color-holonomy.py'))
m = runpy.run_path(str(root/'ch39-gaussian-momentum.py'))
q = m['real']
g,h = color['generators']
identity = color['identity'](3)
edges = {(i,j):(d,p) for i,j,d,p in m['shell']['edges']}

def third_moment(c):
    # Trace all unoriented triangles once. The two orientations contribute six rooted walks.
    z = m['x']*m['y'] + c*m['su']*m['sv']
    shift = q(0)
    for i,j,k in combinations(range(19),3):
        if any(e not in edges for e in ((i,j),(j,k),(i,k))):
            continue
        word = edges[i,j], edges[j,k], edges[i,k]
        lengths = [d for d,p in word]
        changed = [d for d in lengths if d in (724,804)]
        if not changed:
            continue
        assert len(changed)<=2
        signs = [sign if d in (724,804) else 0 for d,sign in zip(lengths,(1,1,-1))]
        if len(changed)==2:
            assert sum(signs)==0
        half = (m['x'] if changed[0]==724 else m['y']) if len(changed)==1 else (q(1) if changed[0]==changed[1] else z)
        matrices = [{724:g,804:h}.get(d,identity) for d in lengths]
        cycle = matrices[0]*matrices[1]*color['adjoint'](matrices[2])
        trace = sum(cycle[t,t] for t in range(3))
        phase = word[0][1]*word[1][1]*word[2][1].conjugate()
        shift -= 12*half*(phase*(trace-3)).real
    return shift

if __name__=='__main__':
    margins = [third_moment(q(c)) for c in ('-.1855','-.127')]
    print('colored minus tripled trace cube on the ground-state band',
          min(v.lower() for v in margins).str(12),max(v.upper() for v in margins).str(12))
    assert all(v>0 for v in margins) or all(v<0 for v in margins)
