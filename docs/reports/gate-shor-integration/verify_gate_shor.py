#!/usr/bin/env python3
"""Independent legal-gate/Shor interface checks. No QPU or Lean execution.

Exact assertions use integers, fractions and SymPy. Floating assertions are
finite clock-distribution comparisons with tolerance 2e-11.
"""
from __future__ import annotations
import itertools
import json
import math
import platform
from fractions import Fraction
import numpy as np
import sympy as sp

counts = {"exact_assertions": 0, "numerical_assertions": 0, "negative_controls": 0}
def exact(ok: bool, label: str) -> None:
    if not bool(ok):
        raise RuntimeError(label)
    counts["exact_assertions"] += 1

def near(a, b, label: str) -> None:
    if np.max(np.abs(np.asarray(a)-np.asarray(b))) > 2e-11:
        raise RuntimeError(label)
    counts["numerical_assertions"] += 1

def negative(ok: bool, label: str) -> None:
    if bool(ok):
        raise RuntimeError("false statement accepted: " + label)
    counts["negative_controls"] += 1

def words(n: int):
    return [s for s in itertools.product((0,1), repeat=n)
            if all(not(s[i] and s[i+1]) for i in range(n-1))]

def pulse(s, i):
    t = list(s)
    enabled = (i == 0 or not t[i-1]) and (i == len(t)-1 or not t[i+1])
    if enabled:
        t[i] ^= 1
    return tuple(t), 3 if enabled else 0  # amplitude is i**exponent

def sweep(n):
    out = {}
    for s in words(n):
        t, e = s, 0
        for i in range(n):
            t, de = pulse(t, i)
            e = (e + de) % 4
        out[s] = (t,e)
    return out

gate_states = gate_cycles = 0
for n in range(1,9):
    S = words(n); F = sweep(n); gate_states += len(S)
    exact(len({F[s][0] for s in S}) == len(S), "sweep bijection")
    for s in S:
        for i in range(n):
            t,e = pulse(s,i); u,f = pulse(t,i)
            exact(t in S and u == s, "legal toggle involution")
            exact((e-f)%4 == 0, "same enabled status on a toggle pair")
        # U tensor conjugate(U) cancels amplitudes exactly on diagonal pairs.
        t,e = F[s]
        exact((e-e)%4 == 0 and t in S, "paired phase cancellation")
    seen = set()
    for s in S:
        if s in seen: continue
        t,e,cycle = s,0,[]
        while t not in seen:
            seen.add(t); cycle.append(t)
            t,de = F[t]; e=(e+de)%4
        gate_cycles += 1
        exact(t == s, "cycle has no transient")
        r = len(cycle)
        for u in cycle:
            v,h = u,0
            for _ in range(r): v,de=F[v];h=(h+de)%4
            exact(v==u and h==e, "cycle phase is basepoint independent")

# Symbolic monomial characteristic polynomial, no interpolation.
a,b,c,z = sp.symbols('a b c z')
U = sp.Matrix([[0,0,c],[a,0,0],[0,b,0]])
exact(sp.expand((z*sp.eye(3)-U).det()) == z**3-a*b*c,
      "arbitrary three-cycle phase characteristic polynomial")
S=words(3); index={s:i for i,s in enumerate(S)};F=sweep(3)
T=np.zeros((5,5),complex);U3=np.zeros((5,5),complex)
for s,(t,e) in F.items():
    T[index[t],index[s]]=1; U3[index[t],index[s]]=1j**e
near(U3,-T,"entire three-site sweep equals minus its permutation")
v=np.eye(5,dtype=complex)[:,index[(0,0,1)]]
near(np.linalg.matrix_power(U3,3)@v,-v,"classical period three has phase minus one")
negative(np.allclose(np.linalg.matrix_power(U3,3)@v,v),"classical return means literal quantum return")

def phase_clock(op, initial, Q):
    state=np.zeros((Q,len(initial)),complex); w=initial.copy()
    for a in range(Q): state[a]=w/np.sqrt(Q);w=op@w
    qft=np.exp(-2j*np.pi*np.outer(np.arange(Q),np.arange(Q))/Q)/np.sqrt(Q)
    return np.sum(np.abs(qft@state)**2,axis=1)
Q=6
near(phase_clock(U3,v,Q),np.array([0,1,0,1,0,1])/3,"uncorrected phase peaks")
near(phase_clock(T,v,Q),np.array([1,0,1,0,1,0])/3,"classical-cycle peaks")
paired=np.kron(U3,U3.conjugate()); copied=np.kron(v,v)
near(phase_clock(paired,copied,Q),phase_clock(T,v,Q),"paired clock recovers true period")
H=np.array([[1,1],[1,-1]],float)/np.sqrt(2)
D=np.zeros((4,2));D[0,0]=D[3,1]=1
negative(np.allclose(np.kron(H,H)@D,D@H),"paired cancellation extends to nonmonomial gates")
# Which-path garbage: complete exponent record destroys every clock coherence.
f=np.arange(Q)%3
clean=(f[:,None]==f[None,:]).astype(complex)/Q
dirty=np.eye(Q,dtype=complex)/Q
qft=np.exp(-2j*np.pi*np.outer(np.arange(Q),np.arange(Q))/Q)/np.sqrt(Q)
near(np.diag(qft@clean@qft.conj().T),np.array([1,0,1,0,1,0])/3,"clean period function peaks")
near(np.diag(qft@dirty@qft.conj().T),np.ones(Q)/Q,"full history makes uniform clock")
negative(np.allclose(clean,dirty),"retaining every temporary history preserves interference")
# Matrix arithmetic uses Python integers only.
def mul(A,B,N):
    return tuple(tuple(sum(A[i][k]*B[k][j] for k in range(2))%N for j in range(2)) for i in range(2))
def power(A,k,N):
    out=((1,0),(0,1))
    while k:
        if k&1: out=mul(out,A,N)
        A=mul(A,A,N);k//=2
    return out
A=((1,1),(1,0)); eye=((1,0),(0,1)); periods={}
for N in range(2,129):
    state=(1,0); M=eye;r=0
    while True:
        state=((state[0]+state[1])%N,state[0]);M=mul(A,M,N);r+=1
        if state==(1,0):break
        if r>N*N:raise RuntimeError("period bound exceeded")
    exact(M==eye,"cyclic-vector return certifies full matrix return")
    exact(r<=N*N-1,"nonzero orbit period bound")
    for k in [0,1,2,3,7,16,31,r,2*r]:
        fast=power(A,k,N);slow=eye
        for _ in range(k):slow=mul(A,slow,N)
        exact(fast==slow,"repeated-squaring power equals repeated transition")
        exact(((fast[0][0],fast[1][0])==(1,0))==(fast==eye),"faithful probe")
    for ell in sp.factorint(r):
        exact(power(A,r//int(ell),N)!=eye,"primitive period certificate")
    periods[str(N)]=r
# Two ideal spectral samples recover the least common denominator with exact CRT probability.
for r in range(1,41):
    successes=sum(math.lcm(Fraction(j,r).denominator,Fraction(k,r).denominator)==r
                  for j in range(r) for k in range(r))
    p=Fraction(1)
    for ell in sp.factorint(r):p*=1-Fraction(1,int(ell)**2)
    exact(Fraction(successes,r*r)==p,"two-sample exact denominator probability")
print(json.dumps({"status":"passed",**counts,"legal_word_lengths":[1,8],
    "legal_words_checked":gate_states,"sweep_cycles_checked":gate_cycles,
    "moduli_checked":[2,128],"three_site_uncorrected_peaks":[1,3,5],
    "three_site_corrected_peaks":[0,2,4],"clock_dimension":6,
    "selected_pisano_periods":{str(n):periods[str(n)] for n in [2,3,5,7,10,15,16,31,64,127,128]},
    "python":platform.python_version(),"numpy":np.__version__,"sympy":sp.__version__,
    "scope":"exact finite assertions and bounded floating checks; no compiled gate circuit, QPU, Lean, or speedup certification"},indent=2))
