"""No opposite-chirality overlap zero modes for either uniform magnetic factor, for every flat twist."""
from pathlib import Path
import runpy
from flint import arb, acb_mat

m=runpy.run_path(str(Path(__file__).with_name('ch39-four-dirac-product.py')))
n=8
gap=arb(4)/5
shift=4*(arb.pi()/(2*n*n)).sin()
# The perturbation is distributed along all forward and backward hops. Along the path
# from the center to the parameter use its arc length in operator norm for the projector.
path_length=2*arb.pi()/(n*n)
rotation=path_length/(2*(gap-shift))
factor=gap/3
for _ in range(6):
    factor=factor*(3-factor**2)/2
error=(1-factor).upper()
for flux in (1,3):
    h=m['magnetic_kernel'](n,flux,arb(0),arb(0))
    ident=m['eye'](2*n*n)
    app=h/3
    for _ in range(6):
        app=(app*(3*ident-app*app))/2
    restricted=acb_mat([[(app[2*i,2*j]+int(i==j))/2 for j in range(n*n)] for i in range(n*n)])
    threshold=arb(1)/50
    assert m['inertia'](restricted,threshold)==0
    overlap=(threshold-error/2).sqrt()-rotation
    print('flux and wrong-chirality positive-projector norm lower',flux,overlap.lower().str(12))
    assert overlap>0
