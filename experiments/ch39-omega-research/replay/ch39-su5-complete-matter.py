"""One-loop test with minimally completed matter representations on the same SU(5) Wilson background."""
from fractions import Fraction
from pathlib import Path
from types import FunctionType
import runpy

from flint import acb, arb

root = Path(__file__).resolve().parent
vectors = runpy.run_path(str(root/'ch39-su5-vector-threshold.py'))
flow = runpy.run_path(str(root/'ch39-spin-phase-closure.py'))
q,D = flow['q'],flow['D']
running = flow['flow']
g = running['g']
gamma = acb.zeta_zero(8).imag
u,v = (flow['shell'][s]*gamma for s in ('a','b'))

def twisted_spectrum(charge):
    # Reuse the original Hermitian inertia and min-max bounds with a private shell geometry.
    shell = dict(flow['shell'])
    shell['edges'] = [(i,j,d,z**(1+charge)) for i,j,d,z in shell['edges']]
    local = dict(flow['spectrum'].__globals__,shell=shell)
    def bind(f):
        return FunctionType(f.__code__,local,f.__name__,f.__defaults__,f.__closure__)
    local['low_matrix'] = bind(flow['low_matrix'])
    return bind(flow['spectrum'])

indices = {
    1: ('1/10','3/2','1'),
    -3: ('3/10','1/2','0'),
    -4: ('4/5','0','1/2'),
    2: ('1/5','0','1/2'),
    6: ('3/5','0','0')
}
assert all(sum(q(row[j]) for row in indices.values()).contains(2) for j in range(3))
readers = {charge:twisted_spectrum(charge) for charge in indices}
vec = [5,3,2]
# beta*tau - (7/2) T_X H = D*G0/2 - common + (7/3) sum T_s shape_s.
# Here shape_s = sum log(s_i)/2; full multiplets share the coefficient of log(r).
vec_slope = [-q('7/2')*i for i in vec]
proj,dot = running['m']['project'],running['m']['dot']
time_projection = proj(vec_slope,running['beta'])
flux_projection = proj(running['beta'],vec_slope)

def solve_axis(axis,width):
    spectra = {charge:reader(u,width,True,v,axis) for charge,reader in readers.items()}
    shape = {charge:sum(s.log()/2 for s in eigs) for charge,eigs in spectra.items()}
    target = [D*x/2-q(8)/3 + q('7/3')*sum(q(row[j])*shape[charge] for charge,row in indices.items())
              for j,x in enumerate(g['G0'])]
    flux = dot(flux_projection,target)
    time = dot(time_projection,target)
    logr = flux + vectors['log_product']/2
    return time,flux,logr

if __name__=='__main__':
    times,fluxes,scales = [],[],[]
    count = 64
    for j in range(count):
        half = arb.pi()/(2*count)
        # Weyl bounds use at most four moving clock edges incident to a vertex.
        time,flux,scale = solve_axis((2*j+1)*half,half+q('.00001'))
        times.append(time)
        fluxes.append(flux)
        scales.append(scale)
    for name,values in [('time',times),('required vector contrast',fluxes),('required log ratio',scales)]:
        print(name,min(v.lower() for v in values).str(12),max(v.upper() for v in values).str(12))
    assert min(v.lower() for v in fluxes)>100 and vectors['exposure']<q('2.12')
