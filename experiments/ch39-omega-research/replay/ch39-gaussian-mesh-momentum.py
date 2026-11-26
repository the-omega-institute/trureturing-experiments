"""Transfer the four-dimensional determinant slope to arbitrary edge-preserving h<=0.01 meshes."""
from pathlib import Path
import runpy
from flint import arb

base=runpy.run_path(str(Path(__file__).with_name('ch39-gaussian-momentum.py')))
real=base['real'];shell=base['shell']
h=real('.01')
rho_max=shell['rho']*shell['theta'](804)
assert rho_max*1024<real('.05')
assert max(sum(i in (a,b) and d==804 for a,b,d,_ in shell['edges']) for i in range(19))==2
assert sum(d==804 for _,_,d,_ in shell['edges'])==6
def cap(v):return min(real(1),v.upper())
def mesh_slope_error(a,b):
    t=rho_max*b
    # s*coth(s/2) <= 2+s^2/6 by the partial fractions for coth.
    full=t*h*h*(2+t/6)/4
    off=full/2
    diag=18*t*h/2
    assert full<=t*h/2
    den=a+real('.08')
    # First derivative in c bounded from the second derivative in beta using both endpoints.
    before=24*((diag*cap(36/(a+4))+18*off)/den**2+off*cap(18/(a+4))/den)
    second=96*(diag/den**3+off/den**2)
    return arb.pi()*(before+second)/2
if __name__=='__main__':
    # Phase-independent static floor from the three inner sites and six unclocked 244 spokes.
    es={(i,j):(d,z) for i,j,d,z in shell['edges']}
    for i in range(9):
        internal=sum(i in pair and pair[1]<9 for pair in es)
        assert shell['V'][i]+sum(i in pair for pair in es)==internal
    assert all(shell['V'][i]>959990 for i in range(9,19))
    assert all(sum(ds==244 and i in pair and 6<=pair[1]<9 for pair,(ds,_) in es.items())==1 for i in range(6))
    assert all(sum(ds==244 and pair[0]<6 and pair[1]==i for pair,(ds,_) in es.items())==2 for i in range(6,9))
    flux=es[6,7][1]*es[7,8][1]*es[6,8][1].conjugate()
    assert abs(1-flux)**2/9>real('.44')
    assert real('.44')/5-real('.08')-18**2/(real(959990)-real('.08'))>0
    lo,hi=base['quadrature'](mesh_slope_error)
    assert lo > real('.11')
    pref=1/arb.pi()**2
    print('one-loop slope/M^4 lower and anti-to-parallel energy gap/M^4 lower',
          (pref*lo).lower().str(12),(2*pref*lo).lower().str(12))
