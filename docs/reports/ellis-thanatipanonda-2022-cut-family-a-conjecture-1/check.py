"""Ellis-Thanatipanonda, "The Arithmetic-Periodicity of cut for C={1,2c}", arXiv:2203.02457v1 (Discrete Appl. Math. 322 (2022)
391-403), Section 7, Family A, Conjecture 1: for C = {1,3} u X or {1,3} u X u Y (X a nonempty set of even numbers >= 4 with least
element 2c, Y a nonempty set of odd numbers >= 5), the nim-sequence of cut is (0,1)^c(+2), i.e. G(n) = 2 floor((n-1)/(2c)) + ((n-1) mod 2).

Check (exact, brute force over all options = partitions of n into d+1 positive parts, d in C): for c = 2..5, Y empty and Y
nonempty, random X and Y (fixed seed), n <= 26: the computed nim-sequence equals the formula.
"""
import random
def brute(C,N):
    from functools import lru_cache
    G=[0]*(N+1)
    def parts(n,p,mx):
        if p==0:
            if n==0: yield ()
            return
        for h in range(min(n-(p-1),mx),0,-1):
            for rest in parts(n-h,p-1,h): yield (h,)+rest
    for n in range(2,N+1):
        opts=set()
        for d in C:
            for pt in parts(n,d+1,n):
                x=0
                for h in pt: x^=G[h]
                opts.add(x)
        g=0
        while g in opts: g+=1
        G[n]=g
    return G
def F(n,c):
    q=(n-1)//(2*c); return 2*q+((n-1)%2)
random.seed(220302457)
ok=True; count=0
for c in (2,3,4,5):
    for withY in (False, True):
        for trial in range(3):
            X={2*c}|{x for x in range(2*c+2,30,2) if random.random()<0.3}
            Y={y for y in range(5,30,2) if random.random()<0.3} if withY else set()
            if withY and not Y: Y={5}
            C=sorted({1,3}|X|Y)
            G=brute(C,26)
            good=all(G[n]==F(n,c) for n in range(1,27))
            ok&=good; count+=1
            print(f"c={c} C={C}: formula holds for n<=26: {good}")
print("families checked:",count)
print("RESULT", "PASS" if ok else "FAIL")
raise SystemExit(0 if ok else 1)
