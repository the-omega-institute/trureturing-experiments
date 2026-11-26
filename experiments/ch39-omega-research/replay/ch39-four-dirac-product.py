"""A gapped product overlap kernel and its factor indices on two finite magnetic lattices."""
from pathlib import Path
import runpy
from flint import acb, acb_mat, arb, ctx

ctx.prec = 256
def eye(n):
    return acb_mat([[int(i==j) for j in range(n)] for i in range(n)])
def adjoint(a):
    return a.conjugate().transpose()
def magnetic_kernel(n, flux, u, v):
    phase = 2*arb.pi()*flux
    gamma = (acb_mat([[0,1],[1,0]]),acb_mat([[0,acb(0,-1)],[acb(0,1),0]]))
    ident = eye(2)
    size = 2*n*n
    d = eye(size)
    for x in range(n):
        for y in range(n):
            forward = ((n*((x+1)%n)+y,(-phase*y/n**2+u/n)),
                       (n*x+(y+1)%n,(phase*x/n if y==n-1 else 0)+v/n))
            a = 2*(n*x+y)
            for g,(to,angle) in zip(gamma,forward):
                b=2*to
                z=acb(0,angle).exp()
                for i in range(2):
                    for j in range(2):
                        d[a+i,b+j] -= z*(ident[i,j]-g[i,j])/2
                        d[b+i,a+j] -= z.conjugate()*(ident[i,j]+g[i,j])/2
    for i in range(1,size,2):
        for j in range(size):
            d[i,j]=-d[i,j]
    return d

def inertia(h,shift):
    mat=acb_mat(h)
    size=h.nrows()
    for i in range(size):
        mat[i,i]-=shift
    ids=list(range(size))
    negative=0
    while ids:
        j=max(ids,key=lambda t:float(abs(mat[t,t].real).mid()))
        ids.remove(j)
        pivot=mat[j,j].real
        assert pivot<0 or pivot>0
        negative+=int(pivot<0)
        for idx,i in enumerate(ids):
            multiplier=mat[i,j]/pivot
            for k in ids[idx:]:
                a=mat[i,k]-multiplier*mat[j,k]
                mat[i,k]=a.real if k==i else a
                if k!=i:
                    mat[k,i]=a.conjugate()
    return negative

if __name__=='__main__':
    indices=[]
    n=8
    gap=arb(4)/5
    perturb=4*(arb.pi()/(2*n)).sin()
    for flux in (1,3):
        h=magnetic_kernel(n,flux,arb(0),arb(0))
        low=inertia(h,-gap)
        high=inertia(h,gap)
        print('flux and inertia counts at gap faces',flux,low,high,flush=True)
        assert low==high and gap>perturb
        indices.append(low-n*n)
    print('all flat twists, index pair and product',indices,indices[0]*indices[1],
          'Wilson kernel gap lower',(gap-perturb).lower().str(12))
    assert indices[0]*indices[1]==3
