"""Recover the four clock correlations from the normalized nine-point shell table."""
from itertools import combinations, product

from flint import arb, arb_mat, ctx


def main():
    ctx.prec=224
    points=[p for p in product(range(5),repeat=3) if sum(x*x for x in p)==17]
    vectors=[[1/arb(2).sqrt()]+[v/arb(34).sqrt() for v in p] for p in points]
    pairs=(((0,3),(1,2)),((0,2),(1,3)),((2,3),),((0,1),))
    # Phase labels ordered as 6,1,2,3; each sine cancels by the signed-pair law.
    def linear(u,v):
        return [sum(u[i]*v[j]+u[j]*v[i] for i,j in slots) for slots in pairs]
    def dot(u,v):
        return sum(x*y for x,y in zip(u,v))
    rows=[]
    obs=[]
    for i,j in combinations(range(len(vectors)),2):
        u,v=vectors[i],vectors[j]
        h=dot(u,v)
        a,b,c=linear(u,v),linear(u,u),linear(v,v)
        row=[2*h*s-h*h*(r+t) for s,r,t in zip(a,b,c)]
        for k in range(4):
            row.append(a[k]**2-h*h*b[k]*c[k])
            for l in range(k+1,4):
                row.append(2*a[k]*a[l]-h*h*(b[k]*c[l]+b[l]*c[k]))
        rows.append(row)
        obs.append((i,j,h,a,b,c))
    M=arb_mat(rows)
    inverse=(M.transpose()*M).solve(M.transpose())
    K=max(sum(abs(inverse[i,j]) for j in range(36)) for i in range(4))
    assert K<74
    # If all shell vectors normalize, |measured^2-target^2|<=2*sigma.
    # Raw squared norms <= 1+3t from the absolute row sum of the phase Gram.
    # t/(1+3t)^2 is increasing to 1/3 and then decreasing to 1/16 at t=1.
    t_lower=arb(".00020563448")
    assert min(t_lower/(1+3*t_lower)**2,arb(1)/16)/(2*74)>arb(".000001387")
    # At the optimal phase endpoint of the arithmetic tradeoff all four
    # correlations have the same positive sign; compute its actual arc table.
    t=arb(".0002056344835116").union(arb(".0002056344835118"))
    errs={}
    for i,j,h,a,b,c in obs:
        measured=(h+t*sum(a))/((1+t*sum(b))*(1+t*sum(c))).sqrt()
        errs[points[i],points[j]]=measured-h
    maximum=max(abs(v) for v in errs.values())
    assert maximum<arb(".000236818") and maximum>arb(".000236817")
    marked=errs[(0,1,4),(2,3,2)]
    shortest=errs[(2,2,3),(2,3,2)]
    assert marked>arb(".0000983107") and marked<arb(".0000983108")
    assert shortest>arb(".0000219965") and shortest<arb(".0000219966")
    print("Nine normalized shell states force the four unsquared phase overlaps, inverse row norm",K.str(10))
    print("For exact offsets 0..41 some shell overlap error exceeds 1.387e-6.")
    print("At the equal-phase-error optimum: marked arc shift",marked.str(11),"one shortest arc shift",shortest.str(11))


if __name__=="__main__":
    main()
