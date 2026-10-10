#!/usr/bin/env python3
"""Bunth–Pitrik–Titkos–Virosztek arXiv:2402.13150v4, §6.1: "whether the Wasserstein complexity (eq:qw-compl-def)
is convex or not is an open question". C_W(Φ) = max_ρ d_A(ρ, Φ(ρ)) with d_A^2 = D_A^2(ρ,ω) - (D_A^2(ρ,ρ)+D_A^2(ω,ω))/2,
D_A^2(ρ,ω) = min over couplings Π ≥ 0 on H⊗H* with tr_{H*}Π = ω, tr_H Π = ρ^T of tr(C Π), C = Σ_j (A_j⊗I - I⊗A_j^T)^2.
Qubit, A = (X, Y, Z). Channels Φ_± = Ad(U_±), U_± = sqrt(1-t) I ± i sqrt(t) Z, midpoint Ψ = (1-t) id + t Ad(Z), t = 1/4.
Checks (cvxpy SDP, SCS/Clarabel):
 (1) the self-cost formula D^2(σ,σ) = 8 - 4 (tr sqrt σ)^2 on random σ;
 (2) d^2(ρ, Φ_±(ρ)) <= 8t = 2 on 300 random ρ (upper bound of Lemma 2, which is proved for all ρ);
 (3) d^2(|+><+|, Ψ(|+><+|)) = 1 + sqrt 3 > 2, so C_W(Ψ) > (C_W(Φ_+) + C_W(Φ_-))/2."""
import sys, numpy as np, cvxpy as cp
rng=np.random.default_rng(3)
I2=np.eye(2); X=np.array([[0,1],[1,0]],complex); Y=np.array([[0,-1j],[1j,0]]); Z=np.diag([1,-1]).astype(complex)
A=[X,Y,Z]
C=sum((np.kron(a,I2)-np.kron(I2,a.T))@(np.kron(a,I2)-np.kron(I2,a.T)) for a in A)
def ptrace2(P):  # trace out second factor (H*): 4x4 -> 2x2
    return [[P[0,0]+P[1,1],P[0,2]+P[1,3]],[P[2,0]+P[3,1],P[2,2]+P[3,3]]]
def D2(rho,om):
    P=cp.Variable((4,4),hermitian=True)
    cons=[P>>0]
    t2=ptrace2(P); t1=[[P[0,0]+P[2,2],P[0,1]+P[2,3]],[P[1,0]+P[3,2],P[1,1]+P[3,3]]]
    rT=rho.T
    for i in range(2):
        for j in range(2):
            cons+= [t2[i][j]==om[i,j], t1[i][j]==rT[i,j]]
    prob=cp.Problem(cp.Minimize(cp.real(cp.trace(C@P))),cons); prob.solve(solver=cp.CLARABEL if 'CLARABEL' in cp.installed_solvers() else cp.SCS)
    return prob.value
def sqrtm(r):
    w,U=np.linalg.eigh((r+r.conj().T)/2); return (U*np.sqrt(np.clip(w,0,None)))@U.conj().T
def d2(rho,om): return D2(rho,om)-(D2(rho,rho)+D2(om,om))/2
def rand_state():
    G=rng.normal(size=(2,2))+1j*rng.normal(size=(2,2)); R=G@G.conj().T; return R/np.trace(R).real
bad=0
for _ in range(20):
    s=rand_state(); v=D2(s,s); f=8-4*np.trace(sqrtm(s)).real**2
    if abs(v-f)>1e-5: bad+=1
print('self-cost formula ok:', bad==0)
t=0.25
U=lambda sgn: np.sqrt(1-t)*I2+sgn*1j*np.sqrt(t)*Z
worst=0
for _ in range(300):
    r=rand_state()
    for sgn in (1,-1):
        u=U(sgn); worst=max(worst,d2(r,u@r@u.conj().T))
for r in [np.array([[.5,.5],[.5,.5]],complex), np.array([[.5,-.5j],[.5j,.5]]), np.eye(2)/2]:
    for sgn in (1,-1):
        u=U(sgn); worst=max(worst,d2(r,u@r@u.conj().T))
print('max d^2 over sampled rho for Phi_+-:', round(worst,6), '(bound 8t = 2)'); bad+= worst>2+1e-5
plus=np.array([[.5,.5],[.5,.5]],complex); psi=(1-t)*plus+t*Z@plus@Z
val=d2(plus,psi); print('d^2(+, Psi(+)) =',round(val,6),' 1+sqrt3 =',round(1+np.sqrt(3),6)); bad+= abs(val-(1+np.sqrt(3)))>1e-4
print('bad=',bad); sys.exit(1 if bad else 0)
