"""Phunphayap-Noppakaew-Pongsriiam (arXiv:2303.13017v1) Question 3. Out(s_b, n) = {s_b(N) : N >= 1, n | N}.
Claim checked: with d = gcd(b-1, n) and s = s_b(n), every element of Out is a multiple of d; every multiple of d at or
above the explicit threshold (b-1)[(s/d)(L+1) - 1] (L least with b^L > n) lies in Out; and Out has elements not
divisible by s beyond the threshold exactly when s > d. Exact digit-automaton reachability over (N mod n, s_b(N)),
for b in {2,3,4,5,10} and n in 1..30, 99, 121, 343 (165 pairs). Cofiniteness in dN gives natural density 1/d.
"""
from math import gcd
def digsum(n,b):
    s=0
    while n: s+=n%b; n//=b
    return s
def out_set(b,n,Smax):
    # reachable (residue, digitsum) over numbers N>=1 built most-significant-digit first
    reach=set()
    frontier=set()
    for dgt in range(1,b):
        if dgt<=Smax: frontier.add((dgt%n,dgt))
    reach|=frontier
    while frontier:
        new=set()
        for (r,s) in frontier:
            for dgt in range(b):
                t=(r,s)
                ns=s+dgt
                if ns>Smax: continue
                st=((r*b+dgt)%n, ns)
                if st not in reach: reach.add(st); new.add(st)
        frontier=new
    return {s for (r,s) in reach if r==0}
ok=True; rows=[]
for b in (2,3,4,5,10):
    for n in list(range(1,31))+[99,121,343]:
        d=gcd(b-1,n); s=digsum(n,b)
        L=1
        while b**L<=n: L+=1
        thr=(b-1)*((s//d)*(L+1)-1)
        Smax=thr+6*s+20
        O=out_set(b,n,Smax)
        all_mult_d = all(u%d==0 for u in O)
        cof = all(u in O for u in range(max(d,thr),Smax+1,d) if u>=thr)
        missing=[u for u in range(d,Smax+1,d) if u not in O]
        nonmult_s = [u for u in O if u%s!=0]
        inf_pred = s>d
        # A ∩ Out nonempty beyond threshold iff s>d
        tail_nonmult = any(u%s!=0 for u in O if u>=thr)
        good = all_mult_d and cof and (tail_nonmult==inf_pred)
        ok &= good
        if not good: print('FAIL',b,n,d,s,thr,missing[:10])
        rows.append((b,n,d,s,thr,max(missing) if missing else 0))
print("checked", len(rows), "pairs (b,n); all Out elements are multiples of d, all multiples of d >= threshold attained, tail non-multiples of s exist iff s>d:", ok)
print("largest missing multiple of d vs proved threshold (sample):", [r for r in rows if r[1] in (7,12,99,121,343)][:12])
print("RESULT", "PASS" if ok else "FAIL")
raise SystemExit(0 if ok else 1)
