"""Enclose every perturbative high-scale gauge crossing for the two-flavor arc model."""
from pathlib import Path
from fractions import Fraction
import runpy
from flint import arb
root=Path(__file__).parent
m=runpy.run_path(str(root/'ch39-global-gauge-crossings.py'))
g=m['g'];q=g['q'];cap=m['cap']
T,W=m['intro']()
# Intro also holds for thresholds capped at 1 TeV: Z_eff <= 1/3. The complete
# fermion one-loop coefficient drops out of both projections. For current T>13
# the fermion exposure length 1 already implies Z_eff=1/3 exactly.
assert T.lower()>13 and W.upper()<71
# If the lightest scalar lay below 1 TeV, the first 18 would lie below 14 TeV
# by the d=256 comparison. Thus truncated logs alone would violate the upper bound.
theta=lambda d:(1-q(d)/1734).acos()
assert q('10.753')/q('1.5')<14**2
assert q(256)/q('1.5')*theta(804)/theta(724)<14**2
assert 18*(T.lower()-q(14).log())>W.upper()
gs=[m['nonneg'](cap)]*3
for _ in range(6):T,W,gs=m['enclose'](T,W,gs,[3*cap]*3)
assert W.upper()<q('50.85') and T.lower()>q('25.59') and T.upper()<q('27.89')
# The same extreme-shape integrator, now retaining the light/subtracted inverse
# couplings and applying an affine intersection operator to all possible roots.
env=(root/'ch39-exposure-envelope.py').read_text()
env=env[:env.index('def certify(')]
a=env.index('def go(');b=env.index('    assert W<52',a)
env=env[:a]+'''def terminal(tau,W,steps=2500):
'''+env[b:]
a=env.index('    R=4*arb.pi()',env.index('def terminal('))
env=env[:a]+'    return xr,4*arb.pi()*G[0]\n'
h={'__file__':str(root/'ch39-exposure-envelope.py')}
exec(compile(env,'terminal_exposures','exec'),h)
eq=g['eff'];beta=g['beta'];hf=g['hf'];D=m['D']
# Project the terminal equations directly onto scalar exposure and time.
# Annihilate the common inverse coupling; project tau off eff and W off beta.
def projection(v,target):
    w=[v[(i+2)%3]-v[(i+1)%3] for i in range(3)]
    denom=m['dot'](w,target)
    return [x/denom for x in w]
PT=projection(eq,beta);PW=projection(beta,eq)
for j in range(5):
    xr,alpha=h['terminal'](T,W)
    normalized=[D*x/2-f/3 for x,f in zip(xr,hf)]
    T=m['cut'](T,m['dot'](PT,normalized))
    W=m['cut'](W,m['dot'](PW,normalized))
    print(j,T.lower().str(10),T.upper().str(10),W.lower().str(10),W.upper().str(10),
          alpha.lower().str(10),alpha.upper().str(10))
assert T>q('26.837') and T<q('26.899') and W<q('49.222') and alpha<q('38.364')
