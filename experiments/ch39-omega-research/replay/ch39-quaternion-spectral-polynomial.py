"""Recover the spectral determinant polynomial of SU(2) clock pairs on the 19-site graph."""
from pathlib import Path
import runpy
from flint import arb,acb,acb_mat,arb_mat,ctx
p=Path(__file__).parent
m=runpy.run_path(str(p/'ch39-three-determinant-certificate.py'))
src=m['src']
one=m['one'];x=m['x'];y=m['y'];z=m['z']
add=m['plus'];mul=m['mul'];times=m['scale']
def poly_add(a,b):
    c=a.copy()
    for k,v in b.items():
        c[k]=c.get(k,0)+v
    return c
# Avoid any branch on an interval containing zero during polynomial operations.
def poly_mul(a,b):
    c={}
    for i,u in a.items():
        for j,v in b.items():
            k=tuple(t+s for t,s in zip(i,j))
            c[k]=c.get(k,0)+u*v
    return c
def minus(a,b):
    return poly_add(a,times(b,-1))
xyz=[x,y,z]
dots=[[None]*3 for _ in range(3)]
for j,v in enumerate(xyz):dots[j][j]=minus(one,poly_mul(v,v))
dots[0][1]=dots[1][0]=minus(z,poly_mul(x,y))
dots[0][2]=dots[2][0]=minus(poly_mul(x,z),y)
dots[1][2]=dots[2][1]=minus(x,poly_mul(y,z))
def polynomial(u):
    po=src['sites']
    I=[6,7,8]; O=[i for i in range(19) if i not in I]
    H=src['matrix'](arb(u),{}) if u>=0 else m['real_energy_matrix'](arb(-u))
    def block(rows,cols):return acb_mat([[H[i,j] for j in cols] for i in rows])
    invA=block(O,O).inv()
    B=block(I,I)
    W=[acb_mat(16,3) for _ in range(3)]
    for ii,i in enumerate(O):
        for jj,j in enumerate(I):
            d=sum((a-b)**2 for a,b in zip(po[i],po[j]))
            ix={724:1,804:2}.get(d,0)
            W[ix][ii,jj]=H[i,j]
    blocks=[[wi.conjugate().transpose()*invA*wj for wj in W] for wi in W]
    orbit=[I[0]]
    for _ in range(2):
        pt=po[orbit[-1]]
        orbit.append(po.index(pt[1:]+pt[:1]))
    omega=acb(-arb(1)/2,arb(3).sqrt()/2)
    out=one.copy()
    for s in range(3):
        V=acb_mat(3,1)
        for i,j in enumerate(orbit):V[I.index(j),0]=omega**(s*i)/arb(3).sqrt()
        tV=V.conjugate().transpose()
        def el(a):return (tV*a*V)[0,0]
        c0=el(B)-sum(el(blocks[j][j]) for j in range(3))
        assert c0.imag.contains(0)
        c=[-el(blocks[i][j]) for i,j in ((0,1),(0,2),(1,2))]
        re=times(one,c0.real)
        for a,t in zip(c,xyz):
            re=poly_add(re,times(t,2*a.real))
        det=poly_mul(re,re)
        for i in range(3):
            for j in range(3):
                det=minus(det,times(dots[i][j],4*c[i].imag*c[j].imag))
        out=poly_mul(out,det)
    norm=sum(out.values())
    assert not norm.contains(0)
    out={idx:v/norm for idx,v in out.items()}
    return minus(out,one)
def evaluate(poly,t):
    return sum(v*t[0]**i*t[1]**j*t[2]**k for (i,j,k),v in poly.items())
def direct(u,U,V):
    H=src['matrix'](arb(u),{}) if u>=0 else m['real_energy_matrix'](arb(-u))
    M=acb_mat(38,38)
    for i in range(19):
        for j in range(19):
            for a in range(2):M[2*i+a,2*j+a]=H[i,j]
    for i,j,d,_ in src['edges']:
        if d not in (724,804):continue
        g=U if d==724 else V
        for a in range(2):
            for b in range(2):
                M[2*i+a,2*j+b]=H[i,j]*g[a,b]
                M[2*j+b,2*i+a]=M[2*i+a,2*j+b].conjugate()
    value=M.det()/H.det()**2-1
    assert value.imag.contains(0)
    return value.real
if __name__=='__main__':
    a,b=arb(3)/7,arb(4)/9
    U=acb_mat([[acb(0,a).exp(),0],[0,acb(0,-a).exp()]])
    V=acb_mat([[b.cos(),b.sin()],[-b.sin(),b.cos()]])
    for shift in (0,15000,-210000):
        f=polynomial(shift)
        val=evaluate(f,[a.cos(),b.cos(),a.cos()*b.cos()])
        expected=direct(shift,U,V)
        assert val.overlaps(expected) and val.rad()<arb('1e-33')
        print(shift,len(f),val.str(12))
