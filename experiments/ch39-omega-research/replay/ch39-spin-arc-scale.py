"""Test exact r = longest/shortest original arc at every perturbative crossing of the spin-balanced tower."""
from pathlib import Path
import runpy
from flint import arb

root = Path(__file__).resolve().parent
m = runpy.run_path(str(root/'ch39-spin-balanced-rge.py'))
q,D,g,B = m['q'],m['D'],m['g'],m['B']
base = m['m']['data']['g']
cap = 4*arb.pi()
theta = m['m']['data']['theta']
arc_ratio = theta(804)/theta(12)
import sys
broaden = __name__=='__main__' and '--finite-arc' in sys.argv[1:]
r = arc_ratio*q('.999').union(q('1.001')) if broaden else arc_ratio
if broaden:
    h=q('.01')
    mass_ratio_error=(1+arb.pi()**2*h*h/6)/(1-h*h)**2
    assert 1/mass_ratio_error > q('.999') and mass_ratio_error < q('1.001')
w = 18*r.log()-sum(m['shape'])
tF = m['m']['data']['m']['model']['geometry'](r)
assert all(t<1 for t in tF)
pt = m['m']['project'](m['d'],m['beta'])
dot = m['m']['dot']
proj = base['proj']
c0,cS,cF,cY = [proj(pt,A) for A in (g['B'],B,g['BF'],g['D'])]
const = D*dot(pt,g['G0'])/2
gs = [base['nonneg'](cap)]*3
ys = [base['nonneg'](3*cap)]*3
T = base['nonneg'](g['G0'][0]*D/(2*(g['beta'][0]-sum(g['D'][0])*3*cap/D)))
def narrow(T,gs,shapes=None,exposure=None):
    shapes=m['shape'] if shapes is None else shapes
    exposure=w if exposure is None else exposure
    for j in range(7):
        # Positive densities are known; clipping at the input can only reduce total exposure.
        wl = base['nonneg'](exposure.upper()) if not T.lower()>max((r.log()-s).upper() for s in shapes) else exposure
        wf = base['nonneg'](q('1/3')) if not T.lower()>1 else q('1/3')
        pred=const-((dot(c0,gs)-dot(cY,ys))*T+dot(cS,gs)*wl+dot(cF,gs)*wf)/D
        if not T.overlaps(pred):return None
        T=base['cut'](T,pred)
        updated=[]
        for i in range(3):
            sm = g['beta'][i]+(dot(g['B'][i],gs)-dot(g['D'][i],ys))/D
            l = m['d'][i]+dot(B[i],gs)/D
            f = 8+dot(g['BF'][i],gs)/D
            # Convex cumulative threshold exposure is below the chord from 0 to its known total.
            max_loss = max(q(0),(sm*T+l*wl+f*wf).upper())
            min_loss = min(q(0),(sm*T).lower())
            hiG = g['G0'][i]-2*min_loss/D
            loG = max(g['G0'][i]-2*max_loss/D,1/gs[i].upper())
            updated.append(base['cut'](gs[i],1/loG.lower().union(hiG.upper())))
        gs=updated
    return T
if __name__=='__main__':
    T=narrow(T,gs)
    print('before integrating projected flow',T.lower().str(12),T.upper().str(12))
    assert T>5 and T<46
    survivors=[]
    for j in range(5,46):
        trial=q(j).union(q(j+1))
        trial=narrow(trial,gs)
        if trial is None:continue
        pred_t,_,_,_=m['projected'](trial,w,4096)
        if pred_t.overlaps(trial):
            survivors.append(base['cut'](trial,pred_t))
    assert survivors
    T=min(v.lower() for v in survivors).union(max(v.upper() for v in survivors))
    for _ in range(3):
        pred_t,pred_w,_,_=m['projected'](T,w)
        T=base['cut'](T,pred_t)
    print('necessary time, exposure upper',T.lower().str(11),T.upper().str(11),pred_w.upper().str(11))
    print('geometric exposure lower',w.lower().str(13),'excess lower',(w-pred_w).lower().str(13))
    assert w>pred_w
