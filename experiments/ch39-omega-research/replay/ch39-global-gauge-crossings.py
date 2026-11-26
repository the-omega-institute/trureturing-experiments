"""A priori bounds on *any* perturbative crossing, projected off the unified inverse coupling."""
from fractions import Fraction
from pathlib import Path
from flint import arb,ctx
root=Path(__file__).parent
g={}
exec(compile((root/'ch39-two-gauge-certificate.py').read_text().split('def run(al,')[0],'gauge','exec'),g)
q=g['q']
D=16*arb.pi()**2
cap=4*arb.pi()
P=[[q(Fraction(n,67)) for n in (5,3,-8)],
   [q(Fraction(n,134)) for n in (-115,333,-218)]]
beta,ds,B,bd,bf,dy=[g[k] for k in ('beta','d','B','BD','BF','D')]
def box(a,b):return q(a).union(q(b))
def nonneg(h):return arb(0).union(h.upper())
def dot(a,b):return sum(x*y for x,y in zip(a,b))
def proj(row,mat):return [sum(a*v[j] for a,v in zip(row,mat)) for j in range(3)]
rc=[[proj(p,mat) for mat in (B,bd,bf,dy)] for p in P]
base=[D*dot(row,g['G0'])/2 for row in P]
def intercept():
    floor=beta[0]-sum(dy[0])*3*cap/D
    assert floor>3
    T=(g['G0'][0]*D/(2*floor)).upper()
    # All gauge contributions in the hypercharge row are nonnegative.
    W=(g['G0'][0]*D/(2*ds[0])).upper()
    return nonneg(T),nonneg(W)
def cut(z,v):
    a=max(z.lower(),v.lower())
    b=min(z.upper(),v.upper())
    assert a<b
    return a.union(b)
def induced(T,W,Ig,Iy,Is,If):
    return [b-(dot(k0,Ig)-dot(k3,Iy)+dot(k1,Is)+dot(k2,If))/D
            for b,(k0,k1,k2,k3) in zip(base,rc)]
def intro():
    T,W=intercept()
    for _ in range(6):
        ig=[nonneg(cap*T)]*3
        iy=[nonneg(3*cap*T)]*3
        ins=[nonneg(cap*W)]*3
        inf=[nonneg(cap/3)]*3
        t,w=induced(T,W,ig,iy,ins,inf)
        T,W=cut(T,t),cut(W,w)
    return T,W
def enclose(T,W,gs,ymax,segments=256):
    # The minimum scalar energy at each inertia follows from the d=256 threshold
    # split. Derive a uniform (loose) width for this W before enclosing the flow.
    highlogs=sum((q(Fraction(int(l.split()[1]),10**6))).log()/2 for l in (root/'ch39-quaternion-band-boxes.txt').read_text().splitlines())/18
    assert W.upper()<q(100)
    # Dirichlet exposure at scaled ratio 2 exceeds 101, uniformly in d.
    width=max((W.upper()/18+highlogs-q('1.5').log()/2).upper(),(2*(q(256)/q('1.5')).sqrt()).log().upper())
    assert T.lower()>width+1
    ig=[arb(0)]*3; active=[None]*3;ferm=[None]*3
    yi=[nonneg(y) for y in ymax]
    for j in range(segments):
        frac=box(Fraction(j,segments),Fraction(j+1,segments))
        time=T*frac
        b=T*(frac-1)
        E=nonneg(W*g['positive'](1+b/width))
        Z=nonneg(g['positive'](1+b)/3)
        zz=[]
        for i in range(3):
            # Exact one-loop terms and cumulative bounds on *all earlier* corrections.
            rg=dot(B[i],gs)-dot(dy[i],yi)
            v=g['G0'][i]-2*(beta[i]*time+ds[i]*E+8*Z)/D-2*(rg*time+dot(bd[i],gs)*E+dot(bf[i],gs)*Z)/D**2
            if v>0:
                z=cut(gs[i],1/v)
            else:
                z=gs[i]
            zz.append(z)
            ig[i]+=T/segments*z
            if not b+width<0:active[i]=z if active[i] is None else active[i].union(z)
            if not b+1<0:ferm[i]=z if ferm[i] is None else ferm[i].union(z)
        if j==0:g_global=zz[:]
        else:g_global=[v.union(w) for v,w in zip(g_global,zz)]
    assert all(x is not None for x in active+ferm)
    QG=[x*W for x in active]
    QF=[x/3 for x in ferm]
    t,w=induced(T,W,ig,[x*T for x in yi],QG,QF)
    return cut(T,t),cut(W,w),g_global
def output(T,W,gs):
    for v in [T,W]+gs:print(v.lower().str(9),v.upper().str(9),end=' | ')
    print()
if __name__=='__main__':
    T,W=intro()
    gs=[nonneg(cap)]*3
    output(T,W,gs)
    for j in range(6):
        T,W,gs=enclose(T,W,gs,[3*cap]*3)
        output(T,W,gs)
