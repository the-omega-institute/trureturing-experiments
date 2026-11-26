"""Sharp exchange coefficient for two aligned, identically distributed endpoint currents."""
from pathlib import Path
import runpy
from flint import arb, arb_mat

root = Path(__file__).resolve().parent
flux = runpy.run_path(str(root/'ch39-su5-shell-flux.py'))
mass = runpy.run_path(str(root/'ch39-su5-wilson-mass.py'))
q = flux['q']
P = arb_mat([[flux['green'][i,j].real for j in range(19)] for i in range(19)])
y = P.solve(arb_mat([[1] for _ in range(19)]))
norm = sum(y[i,0] for i in range(19))
minimum = 1/norm
pole = flux['shell']['sites'].index((17,17,17))
# The original edge phase has zero oriented area on every edge incident to the pole.
for i,j,ds,z in flux['shell']['edges']:
    if pole in (i,j):
        assert ds!=200
        assert flux['shell']['det']((17,17,17),flux['shell']['sites'][i],flux['shell']['sites'][j])==0
        assert z.real>0
assert y[pole,0].contains(0)
assert all(y[i,0]>0 for i in range(19) if i!=pole)
assert minimum>q('.0351096')
# For one chirality and a massless positron in relativistic normalization:
# Gamma = mp*(1-mpi^2/mp^2)^2/(32*pi) * |C W|^2.
# The bound constrains the combined low-energy flavor, running and hadron factor W.
_,k,tau,inverse = mass['antiparallel_crossing']()
internal = 1000*tau.exp()/mass['r']
coeff = 4*arb.pi()/inverse.upper()*flux['volume'].lower()*minimum.lower()/internal.upper()**2
year = q('365.25')*86400
rate = q('6.582119569e-25')/(q('2.4e34')*year)
mp,mpi = q('.938272'),q('.134977')
amplitude = (32*arb.pi()*rate/(mp*(1-(mpi/mp)**2)**2)).sqrt()
effective = amplitude/coeff.lower()
assert coeff>0 and effective>0
print('minimum aligned Green factor',minimum.str(12))
print('tree coefficient floor in GeV^-2',coeff.lower().str(12))
print('allowed effective flavor-running-hadron factor in GeV^2',effective.upper().str(12))
