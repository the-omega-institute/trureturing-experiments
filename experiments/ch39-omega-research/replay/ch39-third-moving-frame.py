"""三条独立非负边权下的全域两列轻态截面。"""
from functools import cache
from pathlib import Path
from runpy import run_path
import sys
import numpy as np
from flint import acb,acb_mat,arb,ctx

FIXED=((0,2),(5,14))
INNER=(6,7,8)
OUTER=tuple(i for i in range(19) if i not in INNER)
ZERO=(0,0)
def adj(A):
    return A.conjugate().transpose()
def sub(A,rows,cols):
    return acb_mat([[A[i,j] for j in cols] for i in rows])
def num(A):
    return np.array(A.tolist(),complex)
def star(pol):
    return {(-a,-b):adj(M) for (a,b),M in pol.items()}
def prod(p,q):
    result={}
    for (a,b),M in p.items():
        for (c,d),N in q.items():
            k=(a+c,b+d)
            v=M*N
            result[k]=v+result[k] if k in result else v
    return {k:v for k,v in result.items() if any(z!=0 for z in v.entries())}
def add(p,q):
    return {k:p[k]+q[k] if k in p and k in q else (p[k] if k in p else q[k]) for k in p.keys()|q.keys()}
def eye(n):
    return acb_mat([[int(i==j) for j in range(n)] for i in range(n)])
def frame(base,links,edges):
    # 三条边构成森林，每个连通分量给一条沿边平行运输的单位列。
    groups=[]
    for i in range(19):
        if any(i in g for g in groups):continue
        queue=[i];wave={i:(0,0,acb(1))}
        for j in queue:
            for f,t in edges:
                if f!=j and t!=j:continue
                l=t if j==f else f
                if l in wave:continue
                atoms=[(power,-M[l,j]) for power,M in ((ZERO,base),((1,0),links[0]),((0,1),links[1]),((-1,0),adj(links[0])),((0,-1),adj(links[1]))) if M[l,j]!=0]
                assert len(atoms)==1
                (a,b),c=atoms[0]
                s,t0,z=wave[j];wave[l]=(s+a,t0+b,c*z);queue.append(l)
        groups.append(wave)
    assert len(groups)==16
    coeff={}
    for j,g in enumerate(groups):
        for i,(a,b,z) in g.items():
            M=coeff.setdefault((a,b),acb_mat(19,16))
            M[i,j]=z/arb(len(g)).sqrt()
    return coeff


@cache
def phase(n,k):
    return acb(0,arb.pi()*k/n).exp()
@cache
def variation(n):
    return 4*arb(2).sqrt()*(arb.pi()/(2*n)).sin()+6*arb("1e-8")


def gate_deleted(base,links,edges):
    M,X,Y=[acb_mat(c) for c in (base,*links)]
    for i,j in edges:
        M[i,i]-=1;M[j,j]-=1
        for A in (M,X,Y):A[i,j]=A[j,i]=0
    D=sub(M,OUTER,OUTER);C=sub(M,INNER,INNER)
    eig=[z.real for z in D.eig(multiple=True)]
    hs=[sub(A,OUTER,INNER) for A in (M,X,Y)]
    @cache
    def schur(c,n):
        t=arb(c)/100+variation(n)
        assert all(z<t or z>t for z in eig)
        need=16-sum(z>t for z in eig)
        R=(D-t*eye(16)).inv()
        center=C-t*eye(3)
        for H in hs:center-=adj(H)*R*H
        return need,center,[-adj(hs[i])*R*hs[j] for i,j in ((0,1),(0,2),(1,2))]
    def check(c,n,i,j):
        need,mid,(x,y,z)=schur(c,n)
        if need<=0:return True
        if need>3:return False
        u=phase(n,2*i+1);v=phase(n,2*j+1)
        F=u*x+v*y+u.conjugate()*v*z
        S=mid+F+adj(F)
        for a,b in ((0,1),(0,2),(1,2)):
            if S[a,a].real>0:
                if need==1:return True
                if S[a,a].real*S[b,b].real>abs(S[a,b])**2:
                    if need==2 or S.det().real>0:return True
        return sum(z.real>0 for z in S.eig(multiple=True))>=need
    return check,[num(A) for A in (M,X,Y)]


def certify(base,links,e,rounds):
    edges=(*FIXED,e)
    full={ZERO:base,(1,0):links[0],(0,1):links[1],(-1,0):adj(links[0]),(0,-1):adj(links[1])}
    P=frame(base,links,edges)
    M=prod(prod(star(P),full),P)
    Q=acb_mat(np.linalg.eigh(num(M[ZERO]))[1][:,:2].tolist())
    X={ZERO:Q}
    step=add({ZERO:28*eye(16)},{k:-v for k,v in M.items()})
    for _ in range(rounds):
        X=prod(step,X)
    gram=prod(star(X),X)
    energy=prod(prod(star(X),M),X)
    # 预留千分之一的全域 Rayleigh 余量。滤波矩阵与能量均按同一 Laurent 多项式精确计算。
    @cache
    def pencil(c):
        return add({k:(arb(c)/100-arb(".001"))*v for k,v in gram.items()},{k:-v for k,v in energy.items()})
    check,deleted=gate_deleted(base,links,edges)
    MD,MX,MY=deleted
    ng={k:num(v) for k,v in gram.items()};ne={k:num(v) for k,v in energy.items()}
    def evaluate(pol,angle):
        return sum(z*np.exp(1j*(a*angle[0]+b*angle[1])) for (a,b),z in pol.items())
    def frame_positive(c,n,i,j):
        poly=pencil(c)
        F=acb_mat(poly[ZERO]);error=[arb(0)]*4
        for (a,b),z in poly.items():
            if not (a>0 or (a==0 and b>0)):continue
            mid=phase(n,(2*i+1)*a+(2*j+1)*b)
            v=z*mid
            F+=v+adj(v)
            rad=arb.pi()*(abs(a)+abs(b))/n
            for ii in range(2):
                for jj in range(2):error[2*ii+jj]+=rad*(abs(z[ii,jj])+abs(z[jj,ii]))
        d0=F[0,0].real-error[0];d1=F[1,1].real-error[3]
        return d0>0 and d1>0 and d0*d1>(abs(F[0,1])+error[1])**2
    todo=[(i,j,8) for i in range(8) for j in range(8)];leaves=0;finest=8
    while todo:
        i,j,n=todo.pop()
        a=np.pi*(2*i+1)/n;b=np.pi*(2*j+1)/n
        off=np.exp(1j*a)*MX+np.exp(1j*b)*MY
        top=np.linalg.eigvalsh(MD+off+off.conj().T)[3]
        bottom=max(np.linalg.eigvals(np.linalg.solve(evaluate(ng,(a,b)),evaluate(ne,(a,b)))).real)
        c=int(round(100*(.4*top+.6*bottom)))
        if frame_positive(c,n,i,j) and check(c,n,i,j):
            leaves+=1;finest=max(finest,n)
        else:
            if n>=1024:raise ArithmeticError((e,i,j,n,top,bottom))
            todo.extend((2*i+x,2*j+y,2*n) for x in (0,1) for y in (0,1))
    return leaves,finest


def obstruction(base,links):
    # 已观察到的单个反例点，只需 12 次单位根。
    u=acb(0,-1)
    v=acb(arb(3).sqrt(),1)/2
    es=(*FIXED,(0,13))
    P=acb_mat(19,16)
    for (a,b),z in frame(base,links,es).items():
        P+=u**a*v**b*z
    off=u*links[0]+v*links[1]
    full=base+off+adj(off)
    deleted=acb_mat(full)
    for i,j in es:
        deleted[i,i]-=1
        deleted[j,j]-=1
        deleted[i,j]=deleted[j,i]=0
    def spectrum(A):
        eig=A.eig(multiple=True)
        assert all(z.imag.contains(0) for z in eig)
        levels=sorted((x.real for x in eig),key=float)
        assert all(a<b for a,b in zip(levels,levels[1:]))
        return levels
    low=spectrum(adj(P)*full*P)[1]
    high=spectrum(deleted)[3]
    assert low-high>arb(".01526")
    print("第三边 (0,13) 荷一，单位辐边反例的压缩第二值、删边第四值及次序反转量",low.str(15),high.str(15),(low-high).str(14),flush=True)


def verify(root):
    shell=run_path(str(root/"ch39-three-ray-chiral-boundary.py"))
    light=run_path(str(root/"ch39-native-light-projector-pump.py"))
    ctx.prec=112
    tasks={1:[((1,9),1),((3,9),1),((3,8),1),((1,7),1)],
           2:[((3,8),1),((1,7),1)],-3:[((1,7),1)]}
    for q,entries in tasks.items():
        base,links=light["clock_matrices"](shell,q)
        if q==1:
            obstruction(base,links)
        for e,rounds in entries:
            counts=certify(base,links,e,rounds)
            print("第三边、补证的荷、滤波次数、叶格数和最细分格",e,q,rounds,*counts,flush=True)

if __name__=="__main__":
    root = Path(sys.argv[1]).resolve() if len(sys.argv)>1 else Path(__file__).resolve().parent
    del sys.argv[1:]
    verify(root)
