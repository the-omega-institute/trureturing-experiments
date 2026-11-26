"""Continuous-metric integral enclosures for the seven-site degree-three Maxwell model."""
from pathlib import Path
import runpy
from flint import arb,acb

h=runpy.run_path(str(Path(__file__).with_name("ch39-wall-neumann-residual.py")))
C,S,angles,dirs = (h[k] for k in ("C","S","angles","dirs"))
modes,order,ii = (h[k] for k in ("modes","order","ii"))
scale,center,centroid = (h[k] for k in ("scale","center","centroid"))
powers,extra,polys = (h[k] for k in ("powers","extra","polys"))
B=h["B"]


def positive_power(lo,hi,e):
    if e==0:
        return arb(1)
    assert e>0
    left = lo.lower()
    lower = left**e if left>0 else arb(0)
    return lower.union(hi.upper()**e)


def gradient_limit(cells=4096):
    # Bound separately on each side and use the maximum modulus principle.
    maxima = [arb(0),arb(0)]
    for i in range(4):
        side=S[i]
        for j in range(cells):
            t=arb(j).union(j+1)/cells
            z=C[i]+side*t
            u=(z-center)/scale
            values=[]
            for coeff in polys:
                value=acb(0)
                for c in reversed(coeff):
                    value=value*u+c
                values.append(value)
            # At incident corners use the exact radial argument to avoid evaluating
            # a complex fractional power on an interval containing zero.
            for v in range(4):
                local=(z-C[v])/(dirs[v]*scale)
                incident=v in (i,(i+1)%4)
                if incident:
                    direction=side if i==v else -side
                    fraction=t if i==v else arb(1)-t
                    magnitude=abs(side)*fraction/scale
                    unit=direction/(abs(direction)*dirs[v])
                else:
                    base=local**(1/angles[v])
                    term=1/local
                for k in range(1,modes+1):
                    if incident:
                        e=k/angles[v]-1
                        part=(-ii)**k * unit**e * positive_power(magnitude,magnitude,e)
                    else:
                        term *= -ii*base
                        part=term
                    for l in range(2):
                        values[l]+=extra[v*modes+k-1][l]*part/dirs[v]
            for l in range(2):
                maxima[l]=max(maxima[l],abs(values[l]).upper())
    return maxima


def primitive(z):
    u=(z-center)/scale
    values=[]
    for coeff in polys:
        value=acb(0)
        for j in range(order-1,-1,-1):
            value=value*u+coeff[j]/(j+1)
        values.append((value*u*scale).real)
    for v in range(4):
        w=(z-C[v])/(dirs[v]*scale)
        base=-ii*w**(1/angles[v])
        term=acb(1)
        for k in range(1,modes+1):
            term *= base
            real=(term*scale*angles[v]/k).real
            for l in range(2):
                values[l]+=extra[v*modes+k-1][l]*real
    return values


def moment_bounds(n):
    m,harmonic_u = gradient_limit()
    print("harmonic derivative bounds",m.str(8),harmonic_u.str(8),flush=True)
    anchor=primitive(centroid)[1] + sum(d[1]*abs(centroid-root).log() for d,root in zip(powers,C))
    # A common harmless shift sets the regular part of D to zero at the marked center.
    gauge=primitive(centroid)[0]+h["particular"](centroid.real)
    a,b=h["a"],h["b"]
    M,N = acb(a),acb(b)
    pieces=[(C[0],M,C[1]),(M,N,C[1]),(C[1],N,C[2]),(N,C[3],C[2])]
    dp=max(abs(h["pprime"](v)).upper() for v in (arb(0),a,b,C[3].real))
    second=2*max(map(abs,B))
    exps=[[arb(d[0])+2*k*d[1] for d in powers] for k in range(3)]
    sums=[arb(0),arb(0),arb(0)]
    for a0,b0,c0 in pieces:
        a1,b1=(b0-a0)/n,(c0-a0)/n
        off=(a1+b1)/3
        r=max(abs(v).upper() for v in [off,a1-off,b1-off])
        mean_sq=sum(abs(v)**2 for v in (off,a1-off,b1-off))/12
        area=abs((a1.conjugate()*b1).imag)/2
        for i in range(n):
            for j in range(n-i):
                node=a0+i*a1+j*b1
                centers=[node+off]
                if i+j+1<n:
                    centers.append(node+2*off)
                for z in centers:
                    reg,ur=primitive(z)
                    reg+=h["particular"](z.real)-gauge
                    ur-=anchor
                    distances=[abs(z-root) for root in C]
                    nearest=[d-r for d in distances]
                    boundary=min((-(z-c)/s).imag.lower()*abs(s).lower() for c,s in zip(C,S))
                    interior=boundary-r
                    for k,e in enumerate(exps):
                        analytic=m+2*k*harmonic_u
                        flat_grad=dp+analytic
                        mid_reg=reg+2*k*ur
                        if interior>0 and min(nearest)>0:
                            mid=(mid_reg+sum(p*d.log() for p,d in zip(e,distances))).exp()
                            grad=flat_grad+sum(p/d for p,d in zip(e,nearest))
                            curvature=second+analytic/interior+sum(p/d**2 for p,d in zip(e,nearest))
                            relative=(grad*r).exp()*(grad**2+curvature)*mean_sq/2
                            if relative < arb(".5"):
                                sums[k]+=area*mid*arb(1,relative.abs_upper())
                                continue
                        span=flat_grad*r
                        enclosure=(mid_reg+arb(0,span.abs_upper())).exp()
                        for d,p in zip(distances,e):
                            enclosure*=positive_power(d-r,d+r,p)
                        sums[k]+=area*enclosure
    # D is fixed up to an arbitrary constant. U is pinned at the physical mark.
    err_d,err_u=h["residual"]()
    operator=h["green_bound"]()
    return [v*(arb(0,(operator*(err_d+4*k*err_u)).abs_upper())).exp() for k,v in enumerate(sums)]


def coherence_upper(moments):
    # Compare endpoint triples. A convex combination of positive triples yields a
    # convex combination of their normalized triples, so the norm is bounded at a vertex.
    from itertools import product
    omega=acb(-arb(1)/2,arb(3).sqrt()/2)
    return max((abs(v[0]+omega*v[1]+omega.conjugate()*v[2])/sum(v)).upper()
               for v in product(*[(s.lower(),s.upper()) for s in moments]))


def reference_boxes():
    return [arb(a).union(b) for a,b in ((7700000,7930000),(4730000,4880000),(6080000,6330000))]


if __name__=="__main__":
    values=moment_bounds(112)
    print("certified integrals in a common gauge",*[(v.lower().str(8),v.upper().str(8)) for v in values],flush=True)
    assert all(box.contains(value) for box,value in zip(reference_boxes(),values))
    bound=coherence_upper(values)
    print("scalar moment contrast at most",bound.str(10))
    assert min(values)>0 and bound<arb(".24")
