"""Projective Riccati flow of the light Higgs quartic below the first heavy doublet."""
from pathlib import Path
from runpy import run_path
from flint import acb_mat, arb, arb_mat

folder = Path(__file__).resolve().parent
low = run_path(str(folder/"ch39-shell-gw-coupled-light-flow.py"))
flow,q,loop = low["flow"],low["q"],low["loop"]


def linear_coefficients(ys,gs,p):
    gram = []
    out = []
    for k in range(p):
        H = [sum((ys[j][i]*flow["adj"](ys[k-j][i]) for j in range(k+1)),acb_mat(3,3))
             for i in range(3)]
        gram.append(H)
        traces = sum(w*sum(m[i,i].real for i in range(3)) for m,w in zip(H,(3,3,1)))
        trace4 = sum(w*sum((gram[r][i]*gram[k-r][i])[j,j].real
                     for r in range(k+1) for j in range(3)) for i,w in enumerate((3,3,1)))
        b = 4*traces-9*gs[k][1]-q(9,5)*gs[k][0]
        c = -2*trace4+q(3,8)*sum(2*gs[r][1]*gs[k-r][1]+
            (gs[r][1]+q(3,5)*gs[r][0])*(gs[k-r][1]+q(3,5)*gs[k-r][0]) for r in range(k+1))
        out.append(arb_mat([[b,c],[-24 if k==0 else 0,0]]))
    return out


def coeffs(m,g,V,p):
    if p==1:
        A = linear_coefficients([m],[g],1)
    else:
        A = linear_coefficients(*low["jets"](m,g,p-1),p)
    co = [V]
    for k in range(p):
        co.append(sum((A[j]*co[k-j] for j in range(k+1)),arb_mat(2,2))/(k+1))
    return co,A[0][0,1]


def step(m,g,V,h):
    m1,g1,tube,gtube = low["advance"](m,g,h)
    pad = arb(0,".15")
    box = arb_mat([[pad,pad],[pad,pad]])
    (v,speed),forcing = coeffs(tube,gtube,V+box,1)
    assert all(abs(h*z)<q(".149") for z in speed.entries())
    assert forcing<0
    inside = V+arb(0).union(h)*speed
    co,_ = coeffs(m,g,V,3)
    rem = coeffs(tube,gtube,inside,4)[0][-1]
    for a in reversed(co):
        rem = a+h*rem
    return m1,g1,rem,inside


def value(V,lam):
    denom = V[1,0]*lam+V[1,1]
    assert denom>0
    return (V[0,0]*lam+V[0,1])/denom


def march(end,initial=None):
    m,g = flow["inputs"](),flow["start"]
    V = arb_mat([[1,0],[0,1]])
    h = q(1,1024) if end>0 else -q(1,4096)
    N = int(float((end/h).lower()))
    bracket = {}
    for i in range(N+1):
        stride = h if i<N else end-i*h
        m,g,V,tube = step(m,g,V,stride)
        assert tube[0,0]>0
        assert tube[1,1]+tube[1,0]*q(".08").union(q(".16"))>0
        if initial is not None:
            output = value(V,initial)
            scale = 1000*(loop*(i*h+stride)).exp()
            if output>0:
                bracket["positive"] = scale
            if output<0 and "negative" not in bracket:
                bracket["negative"] = scale
    return V,bracket


def check():
    # Eq. (4) of Antusch, Hinze and Saad quotes the quartic with a factor 1/4
    # in the potential at M_Z; the gauge and flavor boundary is their Table 2.
    mz = q("91.1876")
    V,_ = march((mz/1000).log()/loop)
    ref = q(".55696").union(q(".56010"))/4
    input_quartic = (V[0,1]-ref*V[1,1])/(ref*V[1,0]-V[0,0])
    print("1 TeV quartic from MZ",float(input_quartic.lower()),float(input_quartic.upper()))
    assert q(".0950") < input_quartic < q(".0960")
    for n,width in ((43,"3e-8"),(195,"1e-8")):
        _,m = flow["masses"](n)
        gap = (m[0]**2-125**2)/q("1.2e12")**2
        eps = 4*q(width)/gap
        b = ((1+eps)/(1-eps)).sqrt()
        orig = m[0]*(1/b).union(b)
        first = (125**2+(orig**2-125**2)*(2*q(".311").union(q(".338"))).exp()).sqrt()
        M,edge = march((first/1000).log()/loop,input_quartic)
        needed = -M[0,1]/M[0,0]
        result = value(M,input_quartic)
        lo,hi = ((".10926",".10936") if n==43 else (".11050",".11059"))
        assert q(lo)<needed<q(hi) and q("-.0225")<result<q("-.0190")
        assert q("9.03e7")<edge["positive"] and edge["negative"]<q("1.436e8")<first
        print("zero",n,"critical input",float(needed.lower()),float(needed.upper()),
              "quartic at first heavy",float(result.lower()),float(result.upper()),
              "zero between GeV",float(edge["positive"].lower()),float(edge["negative"].upper()))


if __name__=="__main__":
    check()
