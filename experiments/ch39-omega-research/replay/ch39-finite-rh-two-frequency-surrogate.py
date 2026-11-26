"""Minimal two-pair Stieltjes quadrature for four RH nodes and its fifth-node discrepancy."""
from pathlib import Path
from runpy import run_path

from flint import acb,arb,arb_mat,arb_series,ctx

record = run_path(str(Path(__file__).with_name("ch39-finite-rh-pivot-mass.py")))


def main():
    ctx.prec,ctx.cap = 224,3
    normalized = [arb(1)]
    for j in (1,2,3):
        num = record["bins"][j][0]
        normalized.append((arb(num).union(arb(num+1))/record["den"])**2)
    W = arb_mat([[1,14,49,36],[1,13,36,0],[1,10,9,0],[1,5,4,0]])
    m = W.solve(arb_mat([[x] for x in normalized])).entries()
    m0,m1,m2,m3=m
    assert m0>0 and m0*m2>m1*m1
    a,b = arb_mat([[m0,m1],[m1,m2]]).solve(arb_mat([[-m2],[-m3]])).entries()
    z = [(-b+s*(b*b-4*a).sqrt())/2 for s in (-1,1)]
    probs=[]
    for x,y in (z,z[::-1]):
        probs.append((m1-m0*y)/(x-y)*(1+x)*(1+4*x)*(1+9*x))
    assert all(q>0 for q in probs) and sum(probs).contains(1)
    gamma=[1/x.sqrt() for x in z]
    assert gamma[0]>arb("48.94271811") and gamma[0]<arb("48.94271813")
    assert gamma[1]>arb("14.62763635") and gamma[1]<arb("14.62763637")
    ordinates = [r.imag for r in acb.zeta_zeros(1,15)]
    for x in gamma:
        assert any(left<x and x<right for left,right in zip(ordinates,ordinates[1:]))
    predicted = sum(p/(1+16*x) for p,x in zip(probs,z))
    def xi(j):
        s=arb_series([arb(1)/2+j,1])
        return s*(s-1)/2*((-s/2)*arb.pi().log()).exp()*(s/2).gamma()*s.zeta()
    anchor,extra=xi(0),xi(4)
    actual = (extra[1]/(8*extra[0]))/(anchor[2]/anchor[0])
    err = actual-predicted
    assert err>arb("7.87824e-8") and err<arb("7.87825e-8")
    # Also evaluate the physical five zero indices of the shell calibration
    # with unit weights in the spectral sum, independently of the fitted shell frame.
    heights = [ordinates[i-1] for i in (4,5,8,12,15)]
    L=[sum(1/(g*g+j*j) for g in heights) for j in range(3)]
    H=arb_mat([[arb(1) if i==j else (i*L[i]+j*L[j])/((i+j)*(L[i]*L[j]).sqrt())
                for j in range(3)] for i in range(3)])
    r=((1-H[0,2]**2)/H.det()).sqrt()
    assert r>arb("3291.4675") and r<arb("3291.4677")
    print("Quadrature frequencies",*[x.str(12) for x in gamma],"normalized probabilities",
          *[x.str(12) for x in probs],"fifth-node squared-overlap excess",err.str(11))
    print("Five unit-weight shell-calibration zeros give middle-node ratio",r.str(12))


if __name__=="__main__":
    main()
