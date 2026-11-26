"""Sharp bounded-error recovery of the static shell amplitude from RH heights."""
from pathlib import Path
from runpy import run_path
from flint import arb, arb_mat

root = Path(__file__).resolve().parent
model = run_path(str(root/'ch39-shell-propagator-readings.py'))
height, groups = model['height'],model['groups']

def row(s,continuous):
    return model['real_row'](model['sample'](s,continuous))

def dot(v,w):
    return sum(a*b for a,b in zip(v,w,strict=True))

def shift(i):
    return (height(i)/height(1))**2

def tail(s0,dual,continuous):
    # On the full mesh D is diagonal and O has absolute row sums at most 18.
    # The RHS is uniform within each orbit; conjugate orbit sums are paired.
    weights = [abs(dual[0])+abs(dual[5]),abs(dual[2])+abs(dual[6])]
    largest = max(weights[0]/6,weights[1]/6,abs(dual[1])/3,
                  abs(dual[3])/3,abs(dual[4]))
    cutoff=10**6
    assert abs(dual[4])/cutoff+18*largest/(cutoff*(cutoff-18))<1
    # Validated real intervals from the end of the height list to the norm bound.
    # Evaluate through the first edge of the diagonal Neumann expansion and
    # bound the paths using two or more edges by the absolute row sums.
    s0=int(s0.lower().floor().unique_fmpz())
    while s0<cutoff:
        right=min(cutoff,s0+s0//40)
        s=arb(s0).union(arb(right))
        h=model['pencil'](s,continuous)
        # All couplings from the pole are real, including the orbit row sums.
        # Quotient diagonal h[j,j] may have within-orbit off-diagonal terms,
        # so use the full physical diagonal for the one-edge term.
        diag=[s for _ in groups]
        for j,k,length,phase in model['row_edges']():
            if continuous:
                x=(s*model['density'](length)).sqrt()
                diag[j]=diag[j]+x/x.tanh()
            else:diag[j]=diag[j]+1
        p=model['pole']
        sums=[dual[0]/6,dual[0]/6,dual[1]/3,dual[2]/6,dual[2]/6,dual[3]/3,dual[4]]
        linear=dual[4]/diag[p]
        for j in range(7):
            if j!=p:linear-=h[p,j].real*sums[j]/(diag[j]*diag[p])
        rem=largest*18**2/(s0**2*(s0-18))
        assert abs(linear)+rem<1
        s0=right

def verify(continuous,ids,lo,hi):
    sample=arb_mat([row(shift(i),continuous) for i in ids])
    static=arb_mat([row(arb(0),continuous)]).transpose()
    coeff=sample.transpose().solve(static).entries()
    signs=[(-1)**j for j in range(7)]
    assert all(c*sign>0 for c,sign in zip(coeff,signs,strict=True))
    dual=sample.solve(arb_mat([[s] for s in signs])).entries()
    factor=sum(abs(c) for c in coeff)
    assert arb(lo)<factor<arb(hi)
    for i in range(1,len(model['prefix'])+1):
        if i not in ids:
            value=abs(dot(row(shift(i),continuous),dual))
            assert value<1
    tail(shift(len(model['prefix'])),dual,continuous)
    print('Continuous edges' if continuous else 'Vertex masses',ids,'amplification',factor.str(15))

if __name__=='__main__':
    verify(False,[1,2,4,8,16,36,120],'60.71351','60.71352')
    shared=[1,2,4,8,16,36,115]
    verify(True,shared,'61.17696','61.17697')
    vertex=arb_mat([row(shift(i),False) for i in shared])
    coeff=vertex.transpose().solve(arb_mat([row(arb(0),False)]).transpose())
    assert arb('60.72885')<sum(abs(c) for c in coeff.entries())<arb('60.72886')
