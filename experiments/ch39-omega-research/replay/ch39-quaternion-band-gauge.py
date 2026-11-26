"""Uniform two-loop crossing and hard mass obstruction while the Dirichlet scale runs from 30M up."""
from pathlib import Path
from fractions import Fraction
import runpy
from flint import arb
root=Path(__file__).parent
v=runpy.run_path(str(root/'ch39-quaternion-gauge.py'))
q=v['q']
raw=(root/'ch39-two-gauge-certificate.py').read_text().split('faces=[')[0]
left=raw.index('    u=(W+sum(logs[:15]))/15')
right=raw.index('    X=G0[:]',left)
raw=raw[:left]+'''    u0=(W+sum(logs))/18
    assert u0.exp()<32 and min(u0.lower(),vstar.lower())>logs[-1].upper()
    u1=(W+sum(logs)+44*vstar)/62
    u=min(u0.lower(),u1.lower()).union(u0.upper())
    exposures=geometry(u.exp())
'''+raw[right:]
raw=raw.replace('E,n=exposure(before+u,arb(0),logs[:15])','''E,n=exposure(before+u,arb(0),logs)
        extra=before+u-vstar
        E+=arb(0).union((44*positive(extra)).upper())
        n+=0 if extra<0 else arb(0).union(44)''')
raw=raw.replace("k1=before+q('.9999');k2=before+q('.0001')","k1=before+exposures[0];k2=before+exposures[1]")
model={'__name__':'band_certificate'}
exec(compile(raw,'band_gauge','exec'),model)
boxes=[list(map(int,line.split())) for line in (root/'ch39-quaternion-band-boxes.txt').read_text().splitlines()]
xs=[q(Fraction(a,10**6)).union(q(Fraction(b,10**6))) for a,b in boxes]
model.update(logs=[w.log()/2 for w in xs],vstar=q('29.9').log(),geometry=v['geometry'])
# Below 32M and d >= 900, only 6+6 edge Dirichlet first levels can fit, each doubled.
angle=lambda d:(1-q(d)/1734).acos()
assert q(30)*(angle(804)/angle(632)).sqrt()>32
# At energy 29.9^2 all edges are below their first pole. Removing the 18 low vertex
# components leaves a positive end matrix: V>=959990, max degree 18.
k=q('29.9')*arb.pi()/30
assert q(959990)-q('29.9')**2-36*k/k.sin()>0
def certificate(a,b,c,d):
    faces=[(a,a,c,d),(b,b,c,d),(a,b,c,c),(a,b,d,d),(a,b,c,d)]
    data=[model['run'](*face) for face in faces]
    for row in data:print([(float(t.lower()),float(t.upper())) for t in row])
    assert data[0][0]<0 and data[1][0]>0 and data[2][1]>0 and data[3][1]<0
    r,alpha=data[4][2:]
    assert alpha<q('38.3316')
    L=(r/xs[0].sqrt()).log()
    assert L<q('3.21')
    print('log width',L, 'first two fermion exposures',v['geometry'](r))
    w=runpy.run_path(str(root/'ch39-spectrum-width-obstruction.py'))
    _,_,gap=w['bound']('3.21')
    # g2, T and width hypotheses certified on this rectangle.
    Q=1000*q(c).union(q(d)).exp()
    mass=gap*(Q/(16*arb.pi()**2))**2
    print('mass floor incl quartic matching and tree contact',mass.lower().str(13))
if __name__=='__main__':
    import sys
    certificate(*(sys.argv[1:] or ['21.461','21.477','26.863','26.875']))
