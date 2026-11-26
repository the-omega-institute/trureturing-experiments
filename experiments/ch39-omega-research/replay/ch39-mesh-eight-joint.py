"""A two-loop gauge and Higgs crossing on the noncommuting ground state of any h<=0.01 mesh."""
from pathlib import Path
import runpy
from flint import arb

root = Path(__file__).resolve().parent
flow = runpy.run_path(str(root/'ch39-spin-phase-closure.py'))
unit = runpy.run_path(str(root/'ch39-unit-pentagon.py'))
m,q = unit['model'],unit['q']
u,v = (m['shell'][s]*m['g'] for s in ('a','b'))
r = flow['r']
r0,r1 = r*q('.999'),r*q('1.001')
t0,t1 = q('26.994'),q('27.008')
T = t0.union(t1)
D = flow['D']
light = runpy.run_path(str(root/'ch39-five-singlet-floor.py'))

def higgs(xs, inverse):
    ratios = r0.union(r1)
    y = [x/ratios**2 for x in xs]
    lh, delta, ms = light['lh'], light['delta'], light['msq']
    Q = 1000*T.exp()
    coupling = 4*arb.pi()/inverse
    tower = flow['data']['m']['model']['geometry'](ratios)
    hard = q('36/5')*coupling**2*sum((1+12*t)*(-2*t).exp() for t in tower)
    response = q(0)
    for mass in y:
        log = mass.log()
        hard += mass*coupling**2*((q('78/25')+3*q('12/25'))*(3*log**2/2+2*log+q('7/2'))+q('156/25')*(1-6*log))
        response += mass*(1-log)
    pair = min(s.lower() for s in y[:2]).union(max(s.upper() for s in y[:2]))
    z = (pair.lower()/pair.upper()).union(q(1))
    log = pair.log()
    sunset = pair*((1+z)*((log-2)**2+1)/2-z*(2-log)*z.log()+(1-z)*q(0).union(arb.pi()**2/6))
    positive = (D*delta/(q('6/5')*(-log)/2)).sqrt()
    negative = positive/5
    hard += q('6/5') * positive**2 * sunset
    charged = (hard+ms*(D/Q)**2)/(29*D*response)
    radial = 12*(negative**2/(4*lh)+q('.0001'))
    limit = 8*arb.pi()
    for sign in (-1,1):
        minor = (limit-sign*radial)*(limit-sign*6*lh-(2*positive**2+1044*charged**2)/limit)-10*negative**2
        assert limit-sign*radial>0 and minor>0
    assert 0<charged and charged<limit and positive<limit and radial/6<limit
    rank_five = D*delta/(8*lh*(-pair.lower().log()/2))
    return positive, charged, rank_five, 1+hard*(Q/D)**2/ms

def refine_spectra(c0,c1):
    # Union fine spectral boxes level by level before enclosing the integral flow.
    enclosures = []
    shape_sum = None
    parts = 48
    for j in range(parts):
        left = (c0+(c1-c0)*j/parts).acos()
        right = (c0+(c1-c0)*(j+1)/parts).acos()
        xs = flow['spectrum'](u,(left-right)/2+q('.00001'),True,v,(left+right)/2)
        shape = sum(x.log()/2 for x in xs)
        shape_sum = shape if shape_sum is None else shape_sum.union(shape)
        if enclosures:
            enclosures = [old.union(new) for old,new in zip(enclosures,xs)]
        else:
            enclosures = xs
    return enclosures, shape_sum

if __name__ == '__main__':
    results = []
    lo,hi = q('-.1855'),q('-.127')
    for j in range(8):
        xs,shape_sum = refine_spectra(lo+(hi-lo)*j/8,lo+(hi-lo)*(j+1)/8)
        faces = ((t0,r0.union(r1),0,1),(t1,r0.union(r1),0,-1),(T,r0,1,1),(T,r1,1,-1))
        margins = []
        for t,scale,k,sign in faces:
            trial = t
            if k:
                for _ in range(2):
                    pred = flow['run'](trial,xs,8192,scale,return_target=True)[0]
                    assert trial.overlaps(pred)
                    trial = trial.intersection(pred)
            projections = flow['run'](trial,xs,32768 if k else 4096,scale)
            value = projections[k]
            if k:
                # The final Arb subtraction encloses A-W for the two operand intervals.
                # Recover bounds on A from R+sup(W) and R+inf(W) before subtracting W_fine.
                wide = projections[2]
                fine = 18*scale.log()-shape_sum
                upper_witness = (wide.mid()+wide.rad()).lower()
                lower_witness = (wide.mid()-wide.rad()).upper()
                target_lo = (value.lower()+upper_witness).lower()
                target_hi = (value.upper()+lower_witness).upper()
                assert target_lo<=target_hi
                value = target_lo.union(target_hi)-fine
            margins.append((sign*value).lower())
        assert min(margins)>0, (j,margins)
        inverse = flow['run'](T,xs,4096,r0.union(r1))[3]
        results.append((min(margins),inverse,*higgs(xs,inverse)))
    for i,name in enumerate(('gauge face margin','unified inverse','positive portal','charged portal','rank-five bound','sensitivity')):
        print(name,min(row[i].lower() for row in results).str(12), max(row[i].upper() for row in results).str(12))
