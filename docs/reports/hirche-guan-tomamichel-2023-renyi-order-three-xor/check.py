#!/usr/bin/env python3
"""Hirche–Guan–Tomamichel arXiv:2305.02589, alpha = 3: numerical check that
H3(X1+X2|B1B2)_tau = -1/2 log[(4 K1 K2 - K1 - K2 + 1)/3], K_i = exp(-2 H3(X_i|B_i)),
and that this equals both branches of Conjecture V.5 (BSC-PSC bound) at alpha = 3,
for random noncommuting density matrices, including rank-deficient supports.
H3(A|B)_rho = -1/2 log Tr[(rho_B^{-1/3} rho_AB rho_B^{-1/3})^3], inverse powers on supp rho_B."""
import sys, numpy as np
rng=np.random.default_rng(int(sys.argv[1]) if len(sys.argv)>1 else 7)
def rand_density(d, rank):
    G=rng.normal(size=(d,rank))+1j*rng.normal(size=(d,rank)); R=G@G.conj().T; return R/np.trace(R).real
def mpow(R,s,tol=1e-12):
    w,U=np.linalg.eigh((R+R.conj().T)/2); f=np.array([x**s if x>tol else 0.0 for x in w]); return (U*f)@U.conj().T
def ptrace_X(rho,dx,db):  # trace out the first (classical) factor
    return rho.reshape(dx,db,dx,db).trace(axis1=0,axis2=2)
def H3(rho,dx,db):
    M=np.kron(np.eye(dx),mpow(ptrace_X(rho,dx,db),-1/3)); T=M@rho@M
    return -0.5*np.log(np.trace(T@T@T).real)
def cq(s0,s1):
    e0=np.diag([1,0]); e1=np.diag([0,1]); return 0.5*np.kron(e0,s0)+0.5*np.kron(e1,s1)
def xor_state(s10,s11,s20,s21):
    d1,d2=s10.shape[0],s20.shape[0]; out=np.zeros((2*d1*d2,)*2,dtype=complex)
    S=[[s10,s11],[s20,s21]]
    for z in (0,1):
        blk=sum(0.25*np.kron(S[0][z^x2],S[1][x2]) for x2 in (0,1))
        e=np.zeros((2,2)); e[z,z]=1; out+=np.kron(e,blk)
    return out,d1*d2
def h3(p): return -0.5*np.log(p**3+(1-p)**3)
def h3inv(H):
    lo,hi=0.0,0.5
    for _ in range(200):
        m=(lo+hi)/2
        if h3(m)<H: lo=m
        else: hi=m
    return (lo+hi)/2
def conv(p,q): return p*(1-q)+(1-p)*q
worst=0.0; n=0
for trial in range(400):
    d1,d2=rng.integers(1,5),rng.integers(1,5)
    ss=[rand_density(d, rng.integers(1,d+1)) for d in (d1,d1,d2,d2)]
    H1=H3(cq(ss[0],ss[1]),2,d1); H2=H3(cq(ss[2],ss[3]),2,d2)
    tau,db=xor_state(*ss); Hout=H3(tau,2,db)
    K1,K2=np.exp(-2*H1),np.exp(-2*H2)
    formula=-0.5*np.log((4*K1*K2-K1-K2+1)/3)
    if H1+H2<=np.log(2): conj=h3(conv(h3inv(H1),h3inv(H2)))
    else: conj=H1+H2-np.log(2)+h3(conv(h3inv(np.log(2)-H1),h3inv(np.log(2)-H2)))
    assert -1e-9<=H1<=np.log(2)+1e-9 and -1e-9<=H2<=np.log(2)+1e-9
    worst=max(worst,abs(Hout-formula),abs(Hout-conj)); n+=1
print(f"trials={n} max_abs_err={worst:.3e}")
sys.exit(0 if worst<1e-8 else 1)
