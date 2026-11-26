"""A common heat-kernel portal in the zero-35 window, with full-tower tree bounds."""
from pathlib import Path
import runpy

from flint import arb

root = Path(__file__).resolve().parent
model = runpy.run_path(str(root / "ch39-dense-six-mode-portal.py"))
q = model["q"]
data = model["data"]
Q,D,g2,lh,msq = (model[k] for k in ("Q","D","g2","lh","msq"))
ys,ls = model["ys"],model["ls"]
R2 = data.scale_lower.union(data.scale_upper)**2
es = [y*R2 for y in ys]
pure, fermion = model["pure"], model["fermion"]
target = 2*arb.pi()**2*(lh-model["ll"])/3


def heat_match(time):
    t = time if isinstance(time,arb) else q(time)
    p = [(-t*e).exp() for e in es]
    weights = [w*w for w in p]
    f = (target / (-sum(l*w for l,w in zip(ls,weights))/2)).sqrt()
    mixed = sum(y*q("144/5")*g2*(3-l)*f*w for y,l,w in zip(ys,ls,p))
    sunset = sum(24*y*((l-2)**2+1)*f**2*w for y,l,w in zip(ys,ls,weights))
    hard = pure+mixed+sunset+fermion
    response = sum(y*(1-l)*w for y,l,w in zip(ys,ls,p))
    b = (hard+msq*(D/Q)**2)/(8*D*response)
    # The omitted tower starts at E/d >= R^2; the mass-normalized finite
    # Dirichlet k-th level is at least k^2/3 for every boundary-zero edge.
    # The Dirichlet subspace has codimension 38, including both flavors.
    v = (-t/3).exp()
    tower = (-t*data.scale_lower**2).exp()*(38+270*v/(1-v))
    psq = sum(weights).upper()+tower.upper()
    pmax = max(p).upper()
    # Heavy interactions c*(sum_i p_i |phi_i|^2)^2 in each color sector.
    color = b-f
    neutral = b+3*f
    cc = max(q(0), (-color).upper())**2/(4*lh)+q(".001")
    cs = q(".005")
    ac = 2*cc*(6*psq+pmax**2)
    an = 2*cs*(2*psq+pmax**2)
    # Shifting the radial Hermitian block by +-cap bounds the full trace
    # space by its Schur complement, including the weighted diagonal.
    cap = q("25.12")
    feasible = max(ac,an) < cap
    if feasible:
        schur = cap-6*lh-12*color**2*psq/(cap-ac)-4*neutral**2*psq/(cap-an)
        feasible = schur>0
    assert max(4*cc*pmax**2, 4*cs*pmax**2, abs(color)*pmax, abs(neutral)*pmax) < cap
    return f,b,pmax,ac,an,feasible, tower, psq, hard*(Q/D)**2/msq+1


if __name__ == "__main__":
    theta = lambda chord: (1-q(chord)/1734).acos()
    t = theta(804)/theta(12)
    f,b,pmax,ac,an,ok,tail,psq,sens = heat_match(t)
    assert ok and q("25.12") < 8*arb.pi()
    for name,val in (("arc length ratio",t), ("f",f), ("b",b), ("trace tail",tail),
                     ("heat square trace upper",psq),("color block upper",ac),
                     ("mass sensitivity",sens), ("relative portal square",(b/f)**2/3)):
        print(name,val.lower().str(12),val.upper().str(12))
    # Cover a continuous range of heat times, recalibrating the two amplitudes
    # at each point. Cellwise interval evaluation includes both endpoints.
    for j in range(1600):
        t = (q(8)+q(j)/200).union(q(8)+q(j+1)/200)
        assert heat_match(t)[5]
