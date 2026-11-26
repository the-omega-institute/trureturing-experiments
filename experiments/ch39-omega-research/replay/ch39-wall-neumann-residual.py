"""Certified harmonic approximation for two scalar Maxwell moment potentials."""
from math import comb
from pathlib import Path
import runpy
import numpy as np
from flint import arb, acb, ctx

work = Path(__file__).resolve().parent
fe = runpy.run_path(str(work/"ch39-wall-maxwell-moments.py"))
unif = fe["uniform"]
geo = fe["geometry"]
alpha,beta,gamma,delta = (unif[k] for k in ("alpha","beta","gamma","delta"))
angles = [alpha,beta,gamma,delta]
ctx.prec = 128
pi,ii = arb.pi(),acb(0,1)
a = arb(18).sqrt()
b = a+2
C = [acb(0),acb(a,5*arb(2).sqrt()),acb(b,7*arb(2).sqrt()),acb(b+7*(arb(2)/3).sqrt())]
S = [C[(i+1)%4]-c for i,c in enumerate(C)]
unit = [s/abs(s) for s in S]
dirs = [u-unit[(i-1)%4] for i,u in enumerate(unit)]
dirs = [d/abs(d) for d in dirs]
center = acb(arb(11)/2,3)
scale = 8
centroid = acb(a+arb(19)/18)
powers = [(1/(3*alpha)-1,1/(3*alpha)),(1/beta-1,0),(1/gamma-1,0),(4,-1)]
fields = {}
for face in geo["geometry"]["big"]:
    tag = "middle" if len(face)==4 else "top" if 14 in face else "base"
    fields[tag] = geo["geometry"]["lifted"](face,geo["negative"])/geo["area"](face)
B = [fields[k] for k in ("top","middle","base")]
order, modes = 28,6


def fit():
    roots = np.array(list(map(complex,C)))
    directions = np.array(list(map(complex,dirs)))
    t = (1-np.cos(np.arange(121)*np.pi/120))/2
    system, rhs = [],[]
    for i in range(4):
        side=complex(S[i])
        n=1j*side/abs(side)
        pts=roots[i]+t*side
        z=(pts-complex(center))/scale
        pol=np.column_stack([n*z**j for j in range(order)])
        basis=[pol.real,-pol.imag]
        for v in range(4):
            for k in range(1,modes+1):
                if v in (i,(i+1)%4):
                    basis.append(np.zeros(len(t)))
                else:
                    g=((pts-roots[v])/(directions[v]*scale))**(k/float(angles[v])-1)
                    basis.append((n*(-1j)**k/directions[v]*g).real)
        system.extend(np.column_stack(basis))
        top=np.real((roots[i]-roots)*n.conjugate())
        top[[i,(i+1)%4]]=0
        deriv=-2*n.real*(float(B[0])*pts.real+float(B[1]-B[0])*np.maximum(0,pts.real-float(a))
                         +float(B[2]-B[1])*np.maximum(0,pts.real-float(b)))
        rhs.extend(-(top/np.maximum(abs(pts[:,None]-roots)**2,1e-60))
                   @np.array(powers,dtype=float)-np.column_stack([deriv,0*deriv]))
    coeff=np.linalg.lstsq(np.array(system),np.array(rhs),rcond=1e-13)[0]
    return [[acb(float(coeff[j,l]),float(coeff[j+order,l])) for j in range(order)]
            for l in range(2)], [[arb(float(v)) for v in row] for row in coeff[2*order:]]


polys, extra = fit()


def positive(x):
    return x.max(0)


def particular(x):
    return -B[0]*x**2-(B[1]-B[0])*positive(x-a)**2-(B[2]-B[1])*positive(x-b)**2


def pprime(x):
    return -2*(B[0]*x+(B[1]-B[0])*positive(x-a)+(B[2]-B[1])*positive(x-b))


def falling(e,m):
    value=arb(1)
    for j in range(m):
        value *= (e-j)/(j+1)
    return value


def edge_geometry(i,root):
    side = S[i]
    off = C[i]-root
    t = (-(off/side).real).max(0).min(1)
    return abs(off+t*side),abs(off).max(abs(off+side))


def normal_error(i,t):
    z=C[i]+S[i]*t
    n=ii*unit[i]
    h=(z-center)/scale
    results=[]
    for coefficients in polys:
        value=acb(0)
        for coef in reversed(coefficients):
            value=value*h+coef
        results.append((n*value).real)
    results[0]+=n.real*pprime(z.real)
    for v in range(4):
        if v in (i,(i+1)%4):
            continue
        along=(z-C[v])/(dirs[v]*scale)
        basic=along**(1/angles[v])
        # Jump-free branch on the convex wedge, and exactly zero normal on its two sides.
        part=1/along
        for k in range(1,modes+1):
            part *= -ii*basic
            term=(n*part/dirs[v]).real
            for l in range(2):
                results[l]+=extra[v*modes+k-1][l]*term
        term=(n/(z-C[v])).real
        for l in range(2):
            results[l]+=powers[v][l]*term
    return results


def nth_bound(i,m):
    side=abs(S[i])
    _,hi=edge_geometry(i,center)
    bounds=[sum(abs(co[j])*comb(j,m)*(side/scale)**m*(hi/scale)**(j-m)
                for j in range(m,order)) for co in polys]
    for v in range(4):
        if v in (i,(i+1)%4):
            continue
        lo,hi=edge_geometry(i,C[v])
        assert lo>0
        for k in range(1,modes+1):
            e=k/angles[v]-1
            bound=abs(falling(e,m))*(side/scale)**m*((lo/scale)**(e-m)).max((hi/scale)**(e-m))
            for l in range(2):
                bounds[l]+=abs(extra[v*modes+k-1][l])*bound
        for l in range(2):
            bounds[l]+=abs(powers[v][l])*side**m/lo**(m+1)
    return bounds


def residual(cells=64,m=8):
    # Chebyshev interpolation at m roots has each cardinal bounded by (2m-1)/m.
    # The product of node distances on one cell is bounded by 2^(1-2m)/cells^m.
    nodes=[(1+((2*j+1)*pi/(2*m)).cos())/2 for j in range(m)]
    maxima=[arb(0),arb(0)]
    for i in range(4):
        sample=[arb(0),arb(0)]
        for j in range(cells):
            for node in nodes:
                err=normal_error(i,(j+node)/cells)
                sample=[max(x,abs(y).upper()) for x,y in zip(sample,err)]
        bounds=nth_bound(i,m)
        for l in range(2):
            bound=(2*m-1)*sample[l]+bounds[l]/(2**(2*m-1)*cells**m)
            maxima[l]=max(maxima[l],bound.upper())
    return maxima


def green_bound():
    t=unif["t"]
    roots=[-1,-ii,(t-ii)/(t+ii),1]
    factor=2*unif["scale"]*abs(1+ii)**(beta-1)*abs(t+ii)**(gamma-1)
    total=arb(0)
    for i,angle in enumerate(angles):
        local=factor*(2/pi)**(angle-1)
        for j in range(4):
            if j!=i:
                local *= (abs(roots[j]-roots[i])/2)**(angles[j]-1)
        total+=2*local*pi**angle*(1/angle**2+arb(2).log()/angle)/pi
    return total


if __name__=="__main__":
    errors=residual()
    operator=green_bound()
    print("supremum normal residuals",*[v.str(10) for v in errors])
    print("Neumann Green bound and errors for normalized potentials",operator.str(10),
          (operator*errors[0]).str(10),(2*operator*errors[1]).str(10))
    assert operator<293 and operator*errors[0]<arb(".002") and 2*operator*errors[1]<arb(".0002")
