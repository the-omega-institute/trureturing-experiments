"""The prespecified triangle mass coefficient three on a phase rectangle."""
from pathlib import Path
import runpy
from fractions import Fraction
from flint import arb

root = Path(__file__).resolve().parent
vac = runpy.run_path(str(root/'ch39-su5-wilson-vacuum.py'))
strip = runpy.run_path(str(root/'ch39-unit-clock-strip.py'))
mass,m,q = vac['mass'],vac['m'],vac['q']
unit = vac['u8']
base_u,base_v = mass['m']['u'],mass['m']['v']
degree = 20
first = Fraction(1,1024)
cutoff = 8192
weights = [(q(0).union(q(first)),q(first**2/2))]
for power in range(-10,13):
    half = q(Fraction(2)**power)/2
    for j in range(degree):
        abscissa,w = arb.legendre_p_root(degree,j,weight=True)
        k = half*(3+abscissa)
        weights.append((k,k*half*w))
tables = {}
for p,dim in vac['net'].items():
    delta = 3*(1-(vac['phase']**p).real)
    tab = [(w,vac['sectors'][p](*unit['node_enclosure'](k+delta))) for k,w in weights]
    # Uniform positive lower bound on the real part of the endpoint form: use .009
    # from the arbitrary-clock nine-edge cycle bound (p=1 uses the inner triangle).
    # |k| <= 5 Re(k) implies |k|/L^2 <= 5/(4*.009), L=Re(k)+.009.
    err = 4*q(12000)*cutoff/2**(2*degree) + m['tail_bound'](cutoff)
    tab.append((arb(0,err.upper()),None))
    tables[p] = tab

def gradient(u,v,c):
    return sum(dim*strip['gradient'](tables[p],u,v,c,0) for p,dim in vac['net'].items())

if __name__=='__main__':
    u=base_u+arb(0,q('.00001'))
    worst=q(100)
    for j in range(16):
        v=base_v+q('.65')+q(Fraction(j,160)).union(q(Fraction(j+1,160)))
        for i in range(8):
            c=q(Fraction(i-4,4)).union(q(Fraction(i-3,4)))
            worst=min(worst,gradient(u,v,c).lower())
    print('phase strip continuous slope lower',worst.str(12),flush=True)
    assert worst>q(2)
    # For mesh transfer reinstate the .08 static floor by a direct uniform inertia check.
    flow = mass['flow']
    for ch in (1,2,6):
        matrix = mass['m']['readers'][ch].__globals__['low_matrix']
        for j in range(4):
            v=base_v+q('.65')+q(2*j+1)/80
            for i in range(64):
                beta=(2*i+1)*arb.pi()/128
                h=matrix(base_u,v,beta)
                assert flow['inertia'](h,q('.25'))==0
    assert q('.25')-4*(arb.pi()/128+q(1)/80+q('.00001'))-flow['high_error']>q('.08')
    finite = worst-sum(vac['net'].values())*(q('.07976')+2*m['tail_bound'](cutoff))
    print('phase strip finite slope lower',finite.lower().str(12),flush=True)
    assert finite>0
    times=[]
    for offset,sign in ((q('.65'),-1),(q('.75'),1)):
        xs={ch:reader(base_u,q('.00001'),True,base_v+offset,arb.pi()) for ch,reader in mass['m']['readers'].items()}
        residual,tau,inv=mass['projection'](xs,q(3))
        assert sign*residual>0
        print('phase face margin',(sign*residual).lower().str(12),flush=True)
    for j in range(8):
        offset=q('.65')+q(2*j+1)/160
        width=q(1)/160+q('.00001')
        xs={ch:reader(base_u,width,True,base_v+offset,arb.pi()) for ch,reader in mass['m']['readers'].items()}
        _,tau,inverse=mass['projection'](xs,q(3))
        assert tau>5 and inverse>1
        for i in range(3):
            matter=sum(q('7/3')*q(mass['m']['indices'][ch][i])*
                       sum((mass['r']/(s+3*mass['penalties'][ch]).sqrt()).log().max(q(0)) for s in eigs)
                       for ch,eigs in xs.items())
            loss=max(q(0),(mass['m']['running']['beta'][i]*tau+matter+q(8)/3).upper())
            assert 4*arb.pi()*(mass['g']['G0'][i]-2*loss/mass['D'])>1
        times.append(tau)
    print('crossing time enclosure',min(t.lower() for t in times).str(12),max(t.upper() for t in times).str(12))
