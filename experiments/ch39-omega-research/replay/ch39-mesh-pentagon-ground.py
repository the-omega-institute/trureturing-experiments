"""Interior curvature transfer and the five-edge ground state on every admissible finite mesh."""
from fractions import Fraction
from pathlib import Path
import runpy

from flint import arb

root = Path(__file__).resolve().parent
curv = runpy.run_path(str(root / 'ch39-gaussian-axis-curvature.py'))
m, mesh = curv['m'], curv['mesh']
q = m['real']
unit = runpy.run_path(str(root / 'ch39-unit-pentagon.py'))
center, p2 = unit['pentagon']()
K = 8192
error = q(0)
for a,b in [(0,1)]+[(2**j,2**(j+1)) for j in range(13)]:
    for j in range(256):
        l = q(Fraction(a)+Fraction((b-a)*j,256))
        r = q(Fraction(a)+Fraction((b-a)*(j+1),256))
        error += (r*r-l*l)*mesh['mesh_slope_error'](l,r)/2
s = abs(m['sv']).upper()
transfer = s*error + 2*m['tail_bound'](K)
def local_curvature_loss(c):
    return s*error*(2/arb.pi()+c)/(1-c*c) + 2*curv['curvature_tail'](K)

def fit_mesh_crossing():
    # At c=-.25, certify a negative phase residual wherever the time projection can vanish.
    f = runpy.run_path(str(root/'ch39-spin-phase-closure.py'))
    u,v = (m['shell'][key]*m['g'] for key in ('a','b'))
    xs = f['spectrum'](u,q('.00001'),True,v,q('-.25').acos())
    T = q('26.98').union(q('27.025'))
    for j in range(3):
        pred,res,_,_ = f['run'](T,xs,8192 if j<2 else 16384,return_target=True)
        assert T.overlaps(pred)
        T = T.intersection(pred)
    assert res < 0
    # To the right of c0 both contributions to the slope are positive.
    # Between -.65 and c0 use a direct second-derivative bound.
    sinwidth = q('.65')
    at_zero = q('-.1567')
    [slope] = unit['slopes_at']([at_zero])
    floor = slope - q('.19632')*(q('.09')-at_zero) - transfer
    required = floor.lower()/(center+q('.25'))
    hessian = required - q('.19632') - local_curvature_loss(sinwidth)
    assert hessian > 0
    # Left of -.65 compare the two slopes relative to the prescribed stationary point.
    # Phi'' >= -.19632 on the entire continuum and each slope error is <=transfer.
    margin = (required-q('.19632'))*(sinwidth-q('.25'))-2*transfer
    print('phase face gap, required quadratic curvature, local Hessian and left slope margin',
          (-res).lower().str(12),required.lower().str(12),hessian.lower().str(12),margin.lower().str(12),flush=True)
    assert hessian > 1 and margin > q('.3')
    return required, hessian

if __name__ == '__main__':
    combined = p2 - q('.19632')
    displacement = transfer/combined
    finite_lo = q('-.1567') - displacement
    finite_hi = q('-.1559') + displacement
    worst_c = max(abs(finite_lo).upper(),abs(finite_hi).upper())
    curvature = combined - local_curvature_loss(worst_c)
    print('unit mesh ground interval and its curvature',finite_lo.lower().str(12),finite_hi.upper().str(12),
          curvature.lower().str(12),flush=True)
    assert finite_lo>q('-.1855') and finite_hi<q('-.127') and curvature>q('2.44')
    fit_mesh_crossing()
