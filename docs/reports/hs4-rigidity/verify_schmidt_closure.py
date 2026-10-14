#!/usr/bin/env python3
"""HS4 direct Schmidt-block classification: finite checks, not a Lean proof.

SymPy checks the fixed coordinate and polynomial identities. NumPy tests
unitary rounding and the stated stability inequalities on seeded samples.
The universal mathematical proof is in HIGUCHI_SUDBERY_RIGIDITY.md, §§8-9.
"""
from __future__ import annotations
import itertools
import json
import math
import platform
import numpy as np
import sympy as sp


def require(test: bool, label: str) -> None:
    if not test:
        raise RuntimeError(label)


def main() -> None:
    exact = numerical = negative = 0
    def eq(a, b, label):
        nonlocal exact
        d = a-b
        good = (all(sp.simplify(x) == 0 for x in d)
                if isinstance(d, sp.MatrixBase) else sp.simplify(d) == 0)
        require(good, label)
        exact += 1
    I, rt = sp.I, sp.sqrt
    eye = sp.eye(2)
    X, Y, Z = (sp.Matrix([[0,1],[1,0]]), sp.Matrix([[0,-I],[I,0]]), sp.diag(1,-1))
    sing = sp.Matrix([0,1,-1,0])/rt(2)
    Q = sp.Matrix.hstack(I*sing, *[sp.kronecker_product(p,eye)*sing for p in (X,Y,Z)])
    YY = sp.kronecker_product(Y,Y)
    eq(Q.H*Q,sp.eye(4),'magic basis orthonormal')
    eq(YY*sp.conjugate(Q),Q,'magic basis spin-flip fixed')
    eq(Q.T*YY*Q,sp.eye(4),'spin-flip bilinear identity')
    w=(-1+I*rt(3))/2
    hp=sp.zeros(16,1)
    for i in (3,12): hp[i]=1/rt(6)
    for i in (5,10): hp[i]=w/rt(6)
    for i in (6,9): hp[i]=w**2/rt(6)
    hm=sp.conjugate(hp)
    def tensor_flat(v): return sp.Matrix(4,4,list(v))
    for sign,h in ((1,hp),(-1,hm)):
        canonical=sp.diag(sign*I/rt(2),1/rt(6),1/rt(6),1/rt(6))
        eq(Q*canonical*Q.T,-tensor_flat(h),f'canonical endpoint {sign}')
    aa,bb,cc,dd=sp.symbols('a b c d',real=True)
    vv=aa*eye-I*(bb*X+cc*Y+dd*Z)
    M=(Q.H*sp.kronecker_product(vv,vv)*Q).applyfunc(sp.expand)
    n=aa**2+bb**2+cc**2+dd**2
    R=sp.Matrix([[aa**2+bb**2-cc**2-dd**2,2*(bb*cc-aa*dd),2*(bb*dd+aa*cc)],
                 [2*(bb*cc+aa*dd),aa**2-bb**2+cc**2-dd**2,2*(cc*dd-aa*bb)],
                 [2*(bb*dd-aa*cc),2*(cc*dd+aa*bb),aa**2-bb**2-cc**2+dd**2]])
    eq(M,sp.diag(n,R),'quaternion collective-rotation polynomial')
    eq(R*R.T,n**2*sp.eye(3),'Rodrigues orthogonality polynomial')
    eq(R.det(),n**3,'Rodrigues orientation polynomial')
    # Generic (real, imaginary) triplet coordinates, not sampled values.
    xs=sp.symbols('x0:9',real=True); ys=sp.symbols('y0:9',real=True)
    us=[x+I*y for x,y in zip(xs,ys)]
    q=-sp.Rational(1,2)+sum(u*u for u in us)/6
    mass=sum(x*x+y*y for x,y in zip(xs,ys))
    eq(sp.re(q)+sum(y*y for y in ys)/3,(mass-3)/6,'realness defect polynomial')
    # Exact global-state spin polynomial for every complex coefficient matrix.
    A=sp.Matrix(4,4,sp.symbols('z0:16'))
    eq(sp.trace((Q*A*Q.T).T*YY*(Q*A*Q.T)*YY),sum(x*x for x in A),'full spin polynomial')
    def chi(v):
        ab=tensor_flat(v)
        ac=sp.Matrix(4,4,[v[k] for k in [0,1,4,5,2,3,6,7,8,9,12,13,10,11,14,15]])
        return sp.simplify(sp.im(ab.det()*sp.conjugate(ac.det())))
    cp,cm=chi(hp),chi(hm)
    eq(cp,rt(3)/864,'positive determinant chirality')
    eq(cm,-rt(3)/864,'negative determinant chirality')
    qn=np.array(Q.evalf(),complex)
    hyn=np.array(sp.kronecker_product(YY,YY),complex)
    hpn=np.array(hp.evalf(),complex).ravel()
    hmn=np.array(hm.evalf(),complex).ravel()
    rng=np.random.default_rng(1406416382)
    def norm(v): return float(np.linalg.norm(v))
    def unitary(k):
        z=rng.normal(size=(k,k))+1j*rng.normal(size=(k,k))
        u,r=np.linalg.qr(z)
        diag=np.diag(r); return u@(np.diag(diag/abs(diag)))
    def marginal(v,keep):
        axes=list(keep)+[i for i in range(4) if i not in keep]
        f=v.reshape((2,)*4).transpose(axes).reshape(2**len(keep),-1)
        return f@f.conj().T
    def metrics(v):
        evs=[np.linalg.eigvalsh(marginal(v,(0,j))) for j in (1,2,3)]
        pur=sum(float(e@e) for e in evs)
        f=sum(float(-np.sum(e[e>0]*np.log2(e[e>0]))) for e in evs)/3
        return pur,f,abs(v@hyn@v)**2,evs
    def check(condition,label):
        nonlocal numerical
        require(bool(condition),label); numerical+=1
    fstar=1+.5*math.log2(3)
    rootmaxerr=0.0
    # All O(3) components, transported by arbitrary product local unitaries.
    for sign in (1,-1):
        for _ in range(80):
            o,_=np.linalg.qr(rng.normal(size=(3,3)))
            o[:,0]*=sign*np.linalg.det(o)
            z=np.zeros((4,4),complex);z[0,0]=1j/math.sqrt(2);z[1:,1:]=o/math.sqrt(6)
            v=(qn@z@qn.T).ravel()
            u4=np.array([[1]],complex)
            for j in range(4): u4=np.kron(u4,unitary(2))
            v=u4@v
            p,f,tau,evs=metrics(v)
            check(abs(p-1)<2e-12,'classified purity')
            check(abs(f-fstar)<2e-12,'classified entropy')
            check(tau<1e-24,'classified spin zero')
            check(all(norm(e-np.array([1/6,1/6,1/6,1/2]))<2e-12 for e in evs),'classified three spectra')
    # Unitary -> real orthogonal rounding for every sampled triplet block.
    for _ in range(200):
        u=unitary(3)
        a,_,bh=np.linalg.svd(u.real); o=a@bh
        qval=-.5+np.sum(u*u)/6
        d2=norm(u-o)**2/6
        check(d2<=abs(qval)+2e-12,'block rounding bound')
        check(abs(np.linalg.det(o)**2-1)<2e-12,'rounding orthogonal target')
    B=1.5*(5*math.log(3)-6)
    scoef=math.sqrt(162*math.log(2));mcoef=math.sqrt(3*math.log(2)/(-B))
    ecoef=6*math.sqrt(2)*mcoef+(12+math.sqrt(6))*scoef
    kcoef=math.sqrt(2)*mcoef+2*ecoef
    delta0=1/(5832*math.log(2))
    def gap(x):
        if x==0:return (4-3*math.log(3))/8
        return ((4-3*math.log(3))/8+(-11+2*math.log(2)+12*math.log(3))/2*x+
                (36-39*math.log(3))/2*x*x+18*(math.log(3)-1)*x**3+x*math.log(x))
    for x in np.linspace(0,1,2001):
        check(gap(x)+2e-14>=min(abs(x-1/6),abs(x-.5))**2/54,'sharpened scalar gap')
    # Global vector rounding. Nonlocal comparison rotations are used only
    # to construct a nearby reference; they are not a state-preparation claim.
    def bell_gauge(v):
        u,_,vh=np.linalg.svd(v.reshape(2,2))
        bell=(u@vh)/math.sqrt(2)
        q0=qn[:,0].reshape(2,2)
        local=2*q0@bell.conj().T
        return np.kron(local,np.eye(2))
    def align(v):
        q0=qn[:,0]; c=np.vdot(q0,v)
        require(abs(c.imag)<1e-10,'positive overlap in gauge')
        c=float(c.real); rem=v-c*q0;s=norm(rem)
        if s<1e-13:return np.eye(4,dtype=complex)
        e=rem/s
        return (np.eye(4)+(c-1)*(np.outer(q0,q0.conj())+np.outer(e,e.conj()))+
                s*(np.outer(q0,e.conj())-np.outer(e,q0.conj())))
    accepted=0
    for scale in (1e-5,1e-4,1e-3,1e-2,3e-2,1e-1):
        for _ in range(35):
            noise=rng.normal(size=16)+1j*rng.normal(size=16)
            noise-=hpn*np.vdot(hpn,noise);noise/=norm(noise)
            v=(hpn+scale*noise)/math.sqrt(1+scale*scale)
            u4=np.array([[1]],complex)
            for j in range(4):u4=np.kron(u4,unitary(2))
            v=u4@v
            p,f,tau,evs=metrics(v);delta=max(0.,fstar-f)
            s=scoef*math.sqrt(delta);m=mcoef*math.sqrt(delta);e=ecoef*math.sqrt(delta)
            if delta>=delta0:continue
            accepted+=1
            check(norm(evs[0]-np.array([1/6,1/6,1/6,1/2]))<=s+3e-11,'nearest-node count and spectrum')
            F=v.reshape(4,4);l,vals,rh=np.linalg.svd(F)
            lv,rv=l[:,0],rh[0,:]
            gl,gr=bell_gauge(lv),bell_gauge(rv)
            G=gl@F@gr.T
            wl,wr=align(gl@lv),align(gr@rv)
            F1=wl@G@wr.T
            ul,_,ur=np.linalg.svd(F1)
            F2=ul@np.diag(np.sqrt([.5,1/6,1/6,1/6]))@ur
            check(norm(G-F2)<=e+3e-11,'global rounding vector distance')
            C=qn.conj().T@F2@qn.conj()
            phase=(1j/math.sqrt(2))/C[0,0]
            C*=phase; G*=phase; F2*=phase
            U=math.sqrt(6)*C[1:,1:]
            a,_,b=np.linalg.svd(U.real);O=a@b
            target=np.zeros((4,4),complex);target[0,0]=1j/math.sqrt(2);target[1:,1:]=O/math.sqrt(6)
            T=qn@target@qn.T
            bound=e+math.sqrt(math.sqrt(2)*m+2*e)
            check(norm(G-T)<=bound+3e-10,'global fourth-root stability')
            q2=F2.ravel()@hyn@F2.ravel()
            check(norm(F2-T)**2<=abs(q2)+3e-10,'post-rounding block bound')
            rootmaxerr=max(rootmaxerr,norm(C[0,1:])+norm(C[1:,0]))
    # Deliberate invalid statements must be detected.
    z=np.diag([1/math.sqrt(2),1/math.sqrt(6),1/math.sqrt(6),1/math.sqrt(6)])
    test=(qn@z@qn.T).ravel();_,f,tau,_=metrics(test)
    require(abs(f-fstar)>1e-4 and tau>.99,'missing spin condition negative control');negative+=1
    cluster=np.zeros(16,complex);cluster[[0,3,12,15]]=[.5,.5,.5,-.5]
    p,f,tau,_=metrics(cluster)
    require(abs(p-1)<1e-12 and tau<1e-24 and abs(f-fstar)>1e-4,'missing spectrum negative control');negative+=1
    require(sp.simplify(cp-cm)!=0,'LU-orbit collapse negative control');negative+=1
    print(json.dumps({'status':'passed','exact_checks':exact,'numerical_checks':numerical,
                      'classified_orbit_samples':160,'unitary_rounding_samples':200,
                      'small_deficit_samples':accepted,'negative_controls':negative,
                      'max_block_off_diagonal_residual':rootmaxerr,
                      'spectral_constant':scoef,'one_body_constant':mcoef,
                      'rounding_constant':ecoef,'spin_error_constant':kcoef,
                      'small_deficit_threshold':delta0,
                      'global_holder_constant':max(math.sqrt(2)/delta0**.25,ecoef*delta0**.25+math.sqrt(kcoef)),
                      'python':platform.python_version(),'numpy':np.__version__,'sympy':sp.__version__,
                      'scope':'symbolic identities and bounded numerical tests; no Lean compilation or universal-proof certification'},indent=2))

if __name__=='__main__':main()
