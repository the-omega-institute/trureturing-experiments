#!/usr/bin/env python3
"""Ferenczi-Zamboni, "Clustering and Arnoux-Rauzy words", arXiv:2305.18986v1; Adv. Appl. Math. 153 (2024) 102621.
Conjecture 1: "In a given AR language Lambda satisfying (Rabc), the word in Proposition [long] is the longest clustering
word, or, if this fails, the longest clustering word conjugate to a standard AR word."
Proposition [long]: v = S_{mu_y}^{mu_y - mu + 1} M_{mu_y} (S, M = shortest, middle standard blocks at stage mu_y).

Family: directive prefix a^m b c a c b a, continued by abcabc (only these stages are used; any continuation using each rule infinitely often has the same stages). U_m = A^2 D A D with A = phi(a), C = phi(c),
D = C A, phi = sigma_a^m o sigma_b o sigma_c. Checks for m = 1..4:
 (1) the source's stages lambda_2, lambda_a, lambda_b, mu_a, mu_b, x, y, mu and the word v of Proposition [long]
     (computed from the AR rules, not from the closed form), |v| = 18m + 22, v perfectly clustering;
 (2) U_m is a factor of a standard block of Lambda; |U_m| = 19m + 24;
 (3) BWT(U_m) for a < c < b equals b^12 c^7 a^(19m+5) (perfect clustering), U_m primitive;
 (4) U_m is conjugate to a standard AR word (of the AR language with prefix a^m b c a a c c);
 (5) |U_m| > |v|."""
import sys
bad=0
def rep(l,ok):
    global bad; print(l,ok); bad+=0 if ok else 1
def blocks(directive, K):
    A,B,C='a','b','c'; out=[(A,B,C)]
    for k in range(K):
        r=directive[k]
        if r=='a': A,B,C=A,B+A,C+A
        elif r=='b': A,B,C=A+B,B,C+B
        else: A,B,C=A+C,B+C,C
        out.append((A,B,C))
    return out
def bwt(w,order):
    rk={ch:i for i,ch in enumerate(order)}
    rots=sorted((w[i:]+w[:i] for i in range(len(w))), key=lambda s:[rk[c] for c in s])
    return ''.join(r[-1] for r in rots)
def perfectly_clustering(w):
    for order in ('abc','acb','bac','bca','cab','cba'):
        B=bwt(w,order); rev=order[::-1]
        target=''.join(ch*w.count(ch) for ch in rev)
        if B==target and len(set(w))>1: return order
    return None
def primitive(w):
    n=len(w); return all(w!=w[:d]*(n//d) for d in range(1,n) if n%d==0)
for m in range(1,5):
    d=('a'*m+'bcacba'+'abc'*2)
    st=blocks(d,len(d))
    lam2=d.index('c'); laa=max(i for i in range(lam2) if d[i]=='a'); lab=max(i for i in range(lam2) if d[i]=='b')
    mua=next(i for i in range(lam2+1,len(d)) if d[i]=='a'); mub=next(i for i in range(lam2+1,len(d)) if d[i]=='b')
    x,y=('a','b') if mua<mub else ('b','a'); muy=max(mua,mub)
    # mu: first stage of the run of rules (z), z in {x, c}, just before mu_y
    j=muy-1; z=d[j]; assert z in (x,'c')
    while j-1>=0 and d[j-1]==z: j-=1
    mu=j
    S,M,L=sorted(st[muy],key=len)
    v=S*(muy-mu+1)+M
    ok=len(v)==18*m+22 and perfectly_clustering(v) is not None and primitive(v)
    rep(f'm={m} (1) stages lam2={lam2} mu_a={mua} mu_b={mub} x={x} y={y} mu={mu}; |v|={len(v)}',ok)
    A='a'+'b'+'a'*m+'c'+'a'*m+'b'+'a'*m; C='c'+'a'*m+'b'+'a'*m; D=C+A
    U=A+A+D+A+D
    big=st[len(d)-1]
    ok=len(U)==19*m+24 and any(U in w for w in big)
    rep(f'm={m} (2) U factor of a standard block, |U|={len(U)}',ok)
    ok=bwt(U,'acb')=='b'*12+'c'*7+'a'*(19*m+5) and primitive(U)
    rep(f'm={m} (3) BWT(U) for a<c<b = b^12 c^7 a^{19*m+5}, primitive',ok)
    d2='a'*m+'bcaacc'+'abc'*2; st2=blocks(d2,len(d2))
    conj={U[i:]+U[:i] for i in range(len(U))}
    ok=any(w in conj for k in range(len(d2)) for w in st2[k])
    rep(f'm={m} (4) U conjugate to a standard AR word',ok)
    rep(f'm={m} (5) |U|={len(U)} > |v|={len(v)}',len(U)>len(v))
print('bad=',bad); sys.exit(1 if bad else 0)
