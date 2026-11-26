"""Positive mass splitting from the original inner triangle in the minimal complete SU(5) tower."""
from pathlib import Path
import runpy
from flint import arb

root = Path(__file__).resolve().parent
all_clocks = runpy.run_path(str(root/'ch39-su5-all-clocks.py'))
m = all_clocks['local']
g, q, D = m['g'], m['q'], m['D']
flow = m['flow']
vec = m['vectors']
phase = runpy.run_path(str(root/'ch39-ground-holonomy.py'))['phase']
penalties = {charge:1-(phase**(1+charge)).real for charge in m['indices']}
r=flow['r']
def projection(spectra,kappa,ratio=r):
    w = {charge:sum((ratio/(s+kappa*penalties[charge]).sqrt()).log().max(q(0)) for s in eigs)
         for charge,eigs in spectra.items()}
    target = [D*x/2-q(8)/3-q('7/3')*sum(q(row[j])*w[c] for c,row in m['indices'].items())
              for j,x in enumerate(g['G0'])]
    h = ratio.log()-vec['log_product']/2
    tau = m['dot'](m['time_projection'],target)
    residual = m['dot'](m['flux_projection'],target)-h
    neutral = 18*ratio.log()-(q(int(19*vec['trees'])).log()+(-18*vec['error']).union(q(0)))/2
    inverse = 4*arb.pi()*(g['G0'][0]-2*(m['running']['beta'][0]*tau+q(8)/3
                +q('7/3')*sum(q(row[0])*w[c] for c,row in m['indices'].items())-q('7/2')*5*(neutral+h))/D)
    return residual,tau,inverse

def antiparallel_crossing():
    xs = {c:read(m['u'],q('.00001'),True,m['v'],arb.pi()) for c,read in m['readers'].items()}
    low,high = q('2.738'),q('2.755')
    assert projection(xs,low)[0] > 0 and projection(xs,high)[0] < 0
    k = low.union(high)
    derivative = q(0)
    for c,row in m['indices'].items():
        delta = k*penalties[c]
        assert all(val+delta<r*r for val in xs[c])
        coeff = q('7/3')*m['dot'](m['flux_projection'],list(map(q,row)))
        derivative += coeff * sum(penalties[c]/(2*(s+delta)) for s in xs[c])
    residual,tau,inverse = projection(xs,k)
    # Compute a lower bound on all intermediate inverse couplings by dropping vector
    # antiscreening. The remaining cumulative one-loop loss is convex in time.
    for j in range(3):
        matter = sum(q('7/3')*q(m['indices'][c][j])*sum((r/(s+k*penalties[c]).sqrt()).log() for s in eigs)
                     for c,eigs in xs.items())
        loss = max(q(0),(m['running']['beta'][j]*tau+matter+q(8)/3).upper())
        assert 4*arb.pi()*(g['G0'][j]-2*loss/D)>1
    assert tau>5 and derivative < 0 and inverse>1
    print('crossing time, unified inverse, residual derivative',
          [(str(a.lower().str(12)),str(a.upper().str(12))) for a in (tau,inverse,derivative)])
    return xs,k,tau,inverse

if __name__=='__main__':
    antiparallel_crossing()
    for c,sign in ((q(0),-1),(q('.3'),1)):
        xs = {ch:read(m['u'],q('.00001'),True,m['v'],c.acos()) for ch,read in m['readers'].items()}
        margin = sign*projection(xs,q(3))[0]
        print('unit-axis phase face, positive margin',c.str(8),margin.lower().str(12))
        assert margin>0
