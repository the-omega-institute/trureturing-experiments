"""Checks for Gallardo-Rahavandrainy, arXiv:2308.01181v2, Conjecture 3.11: for A = (M^2+M+1)^n, n >= 2, the odd Collatz sequence
has length n+1. Part 1: the source process in F2[x] (bit-vector polynomials) for n = 2..80, with the degree sequences of the
source table for n <= 12. Part 2: in the y-process (y = x(x+1)), h(q f^2) = 2h(f)+1 and h(f^2) = 2h(f) for every odd f of
degree <= 10, and h(q^n) = n for n = 0..199."""
def mul(a,b):
    r=0
    while b:
        if b&1: r^=a
        a<<=1; b>>=1
    return r
def divmod2(a,b):
    q=0; db=b.bit_length()
    while a and a.bit_length()>=db:
        s=a.bit_length()-db; q|=1<<s; a^=b<<s
    return q,a
X=0b10; X1=0b11; M=0b111
def oddpart(p):
    while p and not p&1: p>>=1          # divide by x
    while True:
        q,r=divmod2(p,X1)
        if r: break
        p=q
    return p
def length(A):
    A=oddpart(A); seq=[A]
    while A!=1:
        A=oddpart(1^mul(M,A)); seq.append(A)
        if len(seq)>100000: return None
    return len(seq), [s.bit_length()-1 for s in seq]
Q=1^mul(M,M)^M      # M^2+M+1
ok=True
for n in range(2,81):
    A=1
    for _ in range(n): A=mul(A,Q)
    L,deg=length(A)
    if n<=12: print(n,deg,L)
    ok&= (L==n+1)
print('Conjecture 3.11 holds for n=2..80:',ok)
ok_source=ok
m=0b11; q=0b111
def T(f):
    g=1^mul(m,f)
    while g and not g&1: g>>=1
    return g
def h(f):
    k=0
    while f!=1:
        f=T(f); k+=1
        if k>10**6: return None
    return k
ok1=ok2=True; cnt=0
for f in range(1,1<<11):
    if not f&1: continue      # odd: f(0)=1
    cnt+=1
    hf=h(f); f2=mul(f,f)
    ok1&= h(mul(q,f2))==2*hf+1
    ok2&= h(f2)==2*hf
print('odd f of degree <= 10:',cnt,'h(q f^2)=2h(f)+1:',ok1,' h(f^2)=2h(f):',ok2)
# and the y-process lengths of q^n
f=1; good=True
for n in range(0,200):
    good&= h(f)==n; f=mul(f,q)
print('h(q^n)=n for n=0..199:',good)
print("RESULT", "PASS" if (ok_source and ok1 and ok2 and good) else "FAIL")
