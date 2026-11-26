"""Transfer the all-phase two-flavor obstruction to arbitrary nonuniform finite subdivisions."""
from pathlib import Path
from itertools import combinations,product
from flint import arb
root=Path(__file__).parent
h=arb(1)/100
d=arb(384)
gap=arb('.92')/(1+arb.pi()**2/d*h*h*(arb('.92')+10)/6)
assert gap>arb('.919')
assert 4/(1+4*arb.pi()**2*h*h/6)>arb('1.62')**2
pts=[v for v in product(range(30),repeat=3) if sum(i*i for i in v)==867]
theta=lambda z:(1-arb(z)/1734).acos()
def WD(s):
    total=arb(0)
    for a,b in combinations(pts,2):
        ds=sum((x-y)**2 for x,y in zip(a,b))
        if ds>804:continue
        z=s*(1-h*h)*(theta(ds)/theta(804)).sqrt()
        if z>1:total+=2*z.log()
        else:assert z<1
    return total
assert WD(arb('1.62'))>54 and WD(arb(2))>101
assert d/arb('.919')*theta(804)/theta(724)/(1-h*h)**2<22**2
print('finite spectral floor and Dirichlet budgets',gap.lower().str(12),WD(arb('1.62')).lower().str(12),WD(arb(2)).lower().str(12))
s=(root/'ch39-all-clock-angles.py').read_text()
# The established continuum floor is .92; the finite clock problem uses .919.
s=s.replace("gap=sector['gap_at'](dstar,xmin,1024)","gap=sector['gap_at'](dstar,xmin,1024)\nxmin=q('.919')")
s=s.replace("q('.92')","q('.919')").replace("xmin=q('.919')\ngap","xmin=q('.92')\ngap")
g={'__file__':str(root/'ch39-all-clock-angles.py'),'__name__':'finite_all_phases'}
exec(compile(s,'all_finite','exec'),g)
