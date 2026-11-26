"""Rank-one source at the common clock vertex, at zeta zero 68."""
from itertools import product
from pathlib import Path
import runpy
from flint import arb,acb,acb_mat

geom=runpy.run_path(str(Path(__file__).with_name("ch39-shell-schur-flux.py")))
shell=geom["shell"]
W=geom["wall"]
J=geom["inside"]
block=geom["block"]
pi=arb.pi()
inner=(6,7,8)
outer=(11,12,18)
faces=geom["geometry"]["big"]
ref=acb.zeta_zero(68).imag


def data(p,q):
    m=shell["matrix"](0,{724:acb(0,p).exp(),804:acb(0,q).exp()})
    c=block(m,J,J).inv()
    b=block(m,W,J)
    column=b*block(c,range(len(J)),(J.index(0),))
    k=block(m,W,W)-b*c*b.conjugate().transpose()
    return k,column,c[J.index(0),J.index(0)].real


def all_faces():
    items=set(faces)
    for f in faces:
        if len(f)==4:
            for g in (f,f[1:]+f[:1]):
                items.update([(g[0],g[1],g[2]),(g[0],g[2],g[3])])
    return sorted(items)


def loop_poly(f,k,v):
    poly=[acb(1)]
    for a,b in zip(f,f[1:]+f[:1]):
        i,j=W.index(a),W.index(b)
        const=-k[i,j]
        linear=-v[i,0]*v[j,0].conjugate()
        new=[acb(0)]*(len(poly)+1)
        for m,z in enumerate(poly):
            new[m]+=const*z
            new[m+1]+=linear*z
        poly=new
    # The leading term is real by the closed rank-one loop factorization.
    assert poly[-1].imag.contains(0)
    poly[-1]=acb(poly[-1].real)
    return poly


def value(poly,t):
    z=poly[-1]
    for w in reversed(poly[:-1]):
        z=z*t+w
    return z


def degrees(polys,t):
    def angle(f):
        z=value(polys[f],t)
        assert z.real>0 or not z.imag.contains(0)
        return z.arg()
    fixed=sum(angle(f) for f in faces if len(f)==3)
    quad=[f for f in faces if len(f)==4]
    def split(f):
        return angle((f[0],f[1],f[2]))+angle((f[0],f[2],f[3]))
    base=int(((fixed+sum(angle(f) for f in quad))/(2*pi)).unique_fmpz())
    cuts={int(((fixed+sum(terms))/(2*pi)).unique_fmpz())
          for terms in product(*[(split(f),split(f[1:]+f[:1])) for f in quad])}
    return base,sorted(cuts)


def fourier_columns(k,v):
    omega=acb(-arb(1)/2,arb(3).sqrt()/2)
    root=arb(3).sqrt()
    # Family basis f_l(j)=omega^(jl)/sqrt(3).
    eig=[sum(-k[W.index(6),W.index(j)]*omega**(l*i) for i,j in enumerate(outer)) for l in range(3)]
    def transform(sites):
        return [sum(v[W.index(j),0]*omega**(-l*i) for i,j in enumerate(sites))/root for l in range(3)]
    return eig,transform(inner),transform(outer)


def cut_endpoint(polys,face,left,right):
    # Each endpoint here is a transverse negative-axis crossing of a single face.
    central=arb("-.0071").union(arb(".0062"))
    poly=polys[face]
    derivative=[j*w for j,w in enumerate(poly) if j]
    assert not value(derivative,central).imag.contains(0) and value(poly,central).real<0
    lo,hi=arb(left),arb(right)
    assert value(poly,lo).imag*value(poly,hi).imag<0
    for _ in range(48):
        mid=(lo+hi)/2
        if value(poly,lo).imag*value(poly,mid).imag<0:
            hi=mid
        else:
            assert value(poly,mid).imag*value(poly,hi).imag<0
            lo=mid
    return lo.union(hi)


def jarlskog(e,a,b,t):
    y=acb_mat([[e[i]*(i==j)-t*a[i]*b[j].conjugate() for j in range(3)] for i in range(3)])
    hh=y*y.conjugate().transpose()
    numerator=abs((hh[0,1]*hh[1,2]*hh[2,0]).imag)
    transfer=max(abs(z).upper() for z in e)
    correction=abs(t)*sum(abs(z)**2 for z in a).sqrt()*sum(abs(z)**2 for z in b).sqrt()
    shift=2*transfer*correction+correction**2
    lower=upper=arb(1)
    for i in range(3):
        for j in range(i):
            separation=abs(abs(e[i])**2-abs(e[j])**2)
            assert separation>2*shift
            lower*=separation-2*shift
            upper*=separation+2*shift
    return numerator/upper,numerator/lower


if __name__=="__main__":
    p,q=(shell[key]*ref for key in ("a","b"))
    k,v,g=data(p,q)
    assert g>0
    polys={f:loop_poly(f,k,v) for f in all_faces()}
    left_face,right_face=(7,14,8),(6,8,14)
    left=cut_endpoint(polys,left_face,"-.0071","-.0070")
    right=cut_endpoint(polys,right_face,".0061",".0062")
    interval=left.union(right)
    for face,poly in polys.items():
        if face not in (left_face,right_face):
            z=value(poly,interval)
            assert z.real>0 or not z.imag.contains(0)
    assert degrees(polys,arb(0))==(3,[3])
    assert degrees(polys,arb("-.0071"))==(2,[2]) and degrees(polys,arb(".0062"))==(2,[2])
    print("maximal connected robust component in effective strength",left.str(15),right.str(15))
    print("same component in original diagonal strength",(left/(1-g*left)).str(15),(right/(1-g*right)).str(15))
    radius=arb(".00001")
    k,v,_=data(p+arb(0,radius),q+arb(0,radius))
    for f in all_faces():
        z=value(loop_poly(f,k,v),arb(0,arb(".001")))
        assert z.real>0 or not z.imag.contains(0)
    e,a,b=fourier_columns(k,v)
    assert all(abs(x-y)>arb(".065") for i,x in enumerate(e) for y in e[:i])
    assert all(abs(x-y.conjugate())>arb(".01") for x in e for y in e)
    assert min(abs(z).lower() for z in a)>arb(".036") and min(abs(z).lower() for z in b)>arb(".192")
    lo,hi=jarlskog(e,a,b,arb(".001"))
    assert lo>arb("4.8e-15") and hi<arb("5.3e-15")
    print("orthonormal ring labels, Jarlskog magnitude between",lo.lower().str(12),hi.upper().str(12))
    assert (radius/pi)**2>arb("1.013e-11")
