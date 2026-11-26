"""18 low levels for the arc-normal pair and a scalar spectral gap above them."""
from pathlib import Path
from fractions import Fraction
import runpy
from flint import arb,acb,acb_mat,ctx
root=Path(__file__).parent
m=runpy.run_path(str(root/'ch39-quaternion-spectral-polynomial.py'))
src=m['src']
q=lambda f:arb(f.numerator)/f.denominator
p1=src['a']*src['h0'];p2=src['b']*src['h0']
s,x=p1.sin_cos();t,y=p2.sin_cos()
c=(arb(8741547)/12563572).sqrt();v=(1-c*c).sqrt()
U=acb_mat([[acb(x,s),0],[0,acb(x,-s)]])
V=acb_mat([[acb(y,c*t),t*v],[-t*v,acb(y,-c*t)]])
# Work in cyclic irreducible sectors of the shell graph. Spin link has a cyclic-invariant pole.
orbits=[];seen=set()
for i,z in enumerate(src['sites']):
    if i in seen:continue
    orb=[i]
    for zz in (z[1:]+z[:1],z[2:]+z[:2]):
        j=src['sites'].index(zz)
        if j not in orb:orb.append(j)
    orbits.append(orb);seen.update(orb)
sectors=[]
omega=acb(-arb(1)/2,arb(3).sqrt()/2)
for r in range(3):
    os=[O for O in orbits if len(O)==3 or r==0]
    T=acb_mat(38,2*len(os))
    for j,O in enumerate(os):
        for h,i in enumerate(O):
            for a in range(2):T[2*i+a,2*j+a]=omega**(r*h)/arb(len(O)).sqrt()
    sectors.append(T)
def count(s):
    H=m['m']['real_energy_matrix'](s)
    mat=acb_mat(38,38)
    for j in range(19):
        for a in range(2):mat[2*j+a,2*j+a]=H[j,j]
    for i,j,d,z in src['edges']:
        G=U if d==724 else V if d==804 else acb_mat([[1,0],[0,1]])
        for a in range(2):
            for b in range(2):
                mat[2*i+a,2*j+b]=H[i,j]*G[a,b]
                mat[2*j+b,2*i+a]=mat[2*i+a,2*j+b].conjugate()
    assert s<arb(460)**2/arb('1.04')
    n=0
    for T in sectors:
        K=T.conjugate().transpose()*mat*T
        ids=list(range(K.nrows()))
        while ids:
            j=max(ids,key=lambda j:float(abs(K[j,j].real).mid()))
            ids.remove(j)
            pivot=K[j,j].real
            assert pivot>0 or pivot<0
            n+=int(pivot<0)
            for h,i in enumerate(ids):
                scale=K[i,j]/pivot
                for k in ids[h:]:
                    z=K[i,k]-scale*K[j,k]
                    K[i,k]=z.real if i==k else z
                    if i!=k:K[k,i]=z.conjugate()
    return n
if __name__=='__main__':
    assert count(q(Fraction(1,10)))==0 and count(q(Fraction(16)))==18
    # Cutting off the low-endpoint subspace leaves positive boundary diagonals
    # at a test energy well below the first Dirichlet pole.
    bound=q(Fraction(100))**2
    assert count(bound)==18
    rows=[]
    for j in range(1,19):
        lo,hi=Fraction(1,10),Fraction(16)
        for _ in range(20):
            mid=(lo+hi)/2
            if count(q(mid))<j:lo=mid
            else:hi=mid
        rows.append((int((10**6*lo).__floor__()),int((10**6*hi).__ceil__())))
    out=root/'ch39-quaternion-low-masses.txt'
    out.write_text(''.join(f'{lo} {hi}\n' for lo,hi in rows))
    print(rows)
