"""Two-loop gauge crossing for two flavor components with geometrically determined Dirac levels."""
from pathlib import Path
from fractions import Fraction
from flint import arb,ctx
root=Path(__file__).parent
# Reuse the RG coefficients and the interval integrator. Replace only the number
# of scalar levels and the fermion spectrum entering the exposure coordinates.
src=(root/'ch39-two-gauge-certificate.py').read_text().split('faces=[')[0]
src=src.replace('sum(logs[:15]))/15','sum(logs[:18]))/18')
src=src.replace('logs[14]<u and u<logs[15]','logs[17]<u and u<logs[18]')
src=src.replace('arb(0),logs[:15]','arb(0),logs[:18]')
src=src.replace('X=G0[:]','exposures=geometry(u.exp())\n    X=G0[:]')
src=src.replace("k1=before+q('.9999');k2=before+q('.0001')","k1=before+exposures[0];k2=before+exposures[1]")
model={'__name__':'quaternion_rg'}
exec(compile(src,'gauge_interval','exec'),model)
q=model['q']
def mass(s):
    return s*(1+s*s)
def ratio(s):
    return arb(2).sqrt()*(1+2*s*s)/(1+s*s)
def dirac_s(r):
    def value(s):
        u=q(s)
        return mass(u)*mass(arb(2).sqrt()*u)-r*r/arb(1).exp()
    def endpoint(low):
        l,h=Fraction(1),Fraction(6)
        for _ in range(52):
            m=(l+h)/2
            y=value(m)
            if (y<0 if low else not y>0):l=m
            else:h=m
        return (l if low else h)
    lo,hi=endpoint(True),endpoint(False)
    assert value(lo)<0 and value(hi)>0
    return q(lo).union(q(hi))
def geometry(r):
    s=dirac_s(r)
    t2=((1-ratio(s.upper()).log())/2).union((1-ratio(s.lower()).log())/2)
    t1=1-t2
    assert 0<t2 and t2<t1 and t1<1
    assert arb(3).sqrt()/arb(2).sqrt()>t2.exp()
    for a,b in [(724,216),(804,292)]:
        theta=lambda d:(1-q(d)/1734).acos()
        beta=theta(a)/theta(b)
        assert beta*(1+beta**2*s.lower()**2)/(1+s.lower()**2)>t1.exp()
    return t1,t2
model['geometry']=geometry
rows=[list(map(int,l.split())) for l in (root/'ch39-quaternion-low-masses.txt').read_text().splitlines()]
shift=q('.0401')
assert 4*q('.01')*q('.03')/q('.03').sin()<shift
assert (q(16)*q('1.04')*arb.pi()**2/460**2).sqrt()<q('.03')
xs=[q(Fraction(a,10**6)).union(q(Fraction(b,10**6)))+arb(0,shift.upper()) for a,b in rows]
model['logs']=[x.log()/2 for x in xs]+[q(100).log()]
def run_rect(a,b,c,d):
    faces=[(a,a,c,d),(b,b,c,d),(a,b,c,c),(a,b,d,d),(a,b,c,d)]
    data=[model['run'](*f) for f in faces]
    for row in data:print([(float(x.lower()),float(x.upper())) for x in row],flush=True)
    assert data[0][0]<0 and data[1][0]>0 and data[2][1]>0 and data[3][1]<0
    r,alpha=data[4][2:]
    print('fermion exposures',[(z.lower().str(12),z.upper().str(12)) for z in geometry(r)])
    return r,alpha
if __name__=='__main__':
    import sys
    run_rect(*(sys.argv[1:] or ['21.466','21.473','26.865','26.873']))
