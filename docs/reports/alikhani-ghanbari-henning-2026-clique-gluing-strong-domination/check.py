"""Alikhani-Ghanbari-Henning, Quaest. Math. 49 (2026) 195-211, Problem 1, Inequality (4):
gamma_st(G1 u_Kr G2) >= gamma_st(G1) + gamma_st(G2) - 1 - Sigma_{X,Y}, with X_i = N_{G1}(u_i) minus A.
Checks (exact, by exhaustive minimum strong dominating sets):
  1. exhaustive: r = 2, 3, connected G1, G2 on r+1..5 vertices containing the clique on 0..r-1, identity gluing:
     no violation of (4), and the proof's construction D_i = (D minus A, restricted to side i) + {max-degree clique
     vertex of G_i} + (all external clique neighbours when D meets A) is strong dominating in G_i;
  2. random: r = 2, 3, 4, graphs on r+1..r+6 vertices with random clique matchings: no violation of (4);
  3. the source's tightness family (t = 2, r = 2, 3): equality in (4);
  4. arXiv v1 (2306.01608v1): the family H_r (r = 2, 3, 4) violates v1 Theorem 4 (r = 2) and v1 Question 1 (6).
"""
import itertools, random, sys
def gst(n, adj):
    deg=[bin(a).count('1') for a in adj]
    # cover[v] = bitmask of dominators allowed for v (neighbors w with deg[v]<=deg[w]) plus v itself
    cov=[(1<<v)|sum(1<<w for w in range(n) if adj[v]>>w&1 and deg[v]<=deg[w]) for v in range(n)]
    best=[n]
    order=sorted(range(n), key=lambda v: bin(cov[v]).count('1'))
    def rec(chosen, k):
        if k>=best[0]: return
        for v in order:
            if not cov[v]&chosen: break
        else:
            best[0]=k; return
        c=cov[v]
        while c:
            w=c&-c; c^=w
            rec(chosen|w, k+1)
    rec(0,0); return best[0]
def connected(n, adj):
    seen=1; fr=1
    while fr:
        nx=0
        for v in range(n):
            if fr>>v&1: nx|=adj[v]
        fr=nx&~seen; seen|=nx
    return seen==(1<<n)-1
def randgraph(n, r, p):
    adj=[0]*n
    for a in range(r):
        for b in range(a+1,r): adj[a]|=1<<b; adj[b]|=1<<a
    for a in range(n):
        for b in range(a+1,n):
            if b>=r and random.random()<p: adj[a]|=1<<b; adj[b]|=1<<a
    return adj
def glue(n1,a1,n2,a2,r):
    n=n1+n2-r; m=lambda v: v if v<r else v+n1-r
    adj=[0]*n
    for v in range(n1): adj[v]|=a1[v]
    for v in range(n2):
        for w in range(n2):
            if a2[v]>>w&1: adj[m(v)]|=1<<m(w)
    return n, adj
def sigma(n1,a1,r):
    return sum(bin(a1[i]>>r).count('1') for i in range(r))
def graphs(n, r):
    free=[(a,b) for a in range(n) for b in range(a+1,n) if b>=r]
    for mask in range(1<<len(free)):
        adj=[0]*n
        for a in range(r):
            for b in range(a+1,r): adj[a]|=1<<b; adj[b]|=1<<a
        for k,(a,b) in enumerate(free):
            if mask>>k&1: adj[a]|=1<<b; adj[b]|=1<<a
        if connected(n,adj): yield adj
def mindom(n, adj):
    # return one minimum strong dominating set (as bitmask) by brute force over sizes
    deg=[bin(a).count('1') for a in adj]
    cov=[(1<<v)|sum(1<<w for w in range(n) if adj[v]>>w&1 and deg[v]<=deg[w]) for v in range(n)]
    for k in range(n+1):
        for D in itertools.combinations(range(n),k):
            m=sum(1<<d for d in D)
            if all(cov[v]&m for v in range(n)): return m
def isdom(n, adj, m):
    deg=[bin(a).count('1') for a in adj]
    return all(m>>v&1 or any(adj[v]>>w&1 and m>>w&1 and deg[v]<=deg[w] for w in range(n)) for v in range(n))
def iden(n2, a2): return a2
ok = True
# 1. exhaustive
for r, N in ((2, 5), (3, 5)):
    viol = cnt = fail = 0
    for n1 in range(r + 1, N + 1):
        for n2 in range(r + 1, N + 1):
            L1 = list(graphs(n1, r)); L2 = list(graphs(n2, r))
            g1s = [gst(n1, a) for a in L1]; g2s = [gst(n2, a) for a in L2]
            for a1, g1 in zip(L1, g1s):
                for a2, g2 in zip(L2, g2s):
                    n, g = glue(n1, a1, n2, a2, r); cnt += 1
                    D = mindom(n, g); G = bin(D).count('1')
                    if G < g1 + g2 - 1 - sigma(n1, a1, r) - sigma(n2, a2, r): viol += 1
                    off = n1 - r; hitA = any(D >> v & 1 for v in range(r))
                    for nn, aa, emb in ((n1, a1, lambda v: v), (n2, a2, lambda v: v if v < r else v + off)):
                        dd = [bin(x).count('1') for x in aa]
                        m = 1 << max(range(r), key=lambda v: dd[v])
                        for v in range(r, nn):
                            if D >> emb(v) & 1: m |= 1 << v
                            if hitA and any(aa[j] >> v & 1 for j in range(r)): m |= 1 << v
                        if not isdom(nn, aa, m): fail += 1
    ok &= viol == 0 and fail == 0
    print(f"exhaustive r={r} n_i<={N}: pairs={cnt} violations_of_(4)={viol} construction_failures={fail}")
# 2. random
random.seed(20261020)
for r in (2, 3, 4):
    trials = viol = 0; slack = None
    while trials < 1500:
        n1 = random.randint(r + 1, r + 6); n2 = random.randint(r + 1, r + 6)
        a1 = randgraph(n1, r, random.choice([.2, .3, .45])); a2 = randgraph(n2, r, random.choice([.2, .3, .45]))
        if not (connected(n1, a1) and connected(n2, a2)): continue
        perm = list(range(r)); random.shuffle(perm); mp = perm + list(range(r, n2)); b2 = [0] * n2
        for v in range(n2):
            for w in range(n2):
                if a2[v] >> w & 1: b2[mp[v]] |= 1 << mp[w]
        n, g = glue(n1, a1, n2, b2, r); trials += 1
        d = gst(n, g) - (gst(n1, a1) + gst(n2, b2) - 1 - sigma(n1, a1, r) - sigma(n2, b2, r))
        slack = d if slack is None else min(slack, d); viol += d < 0
    ok &= viol == 0
    print(f"random r={r}: trials={trials} violations_of_(4)={viol} min_slack={slack}")
# 3. tightness family
def tightG1(t, r):
    E = [(a, b) for a in range(r) for b in range(a + 1, r)]; n = r
    for i in range(t):
        w = n; n += 1; E.append((0, w))
        for arm in range(t + r - 1):
            a, b, c = n, n + 1, n + 2; n += 3; E += [(w, a), (a, b), (b, c)]
    adj = [0] * n
    for a, b in E: adj[a] |= 1 << b; adj[b] |= 1 << a
    return n, adj
for r in (2, 3):
    n1, a1 = tightG1(2, r); n, g = glue(n1, a1, n1, a1, r); g1 = gst(n1, a1); gg = gst(n, g); s = 2 * sigma(n1, a1, r)
    ok &= gg == 2 * g1 - 1 - s
    print(f"tightness t=2 r={r}: gamma_st(G1)={g1} gamma_st(G)={gg} Sigma={s} equality={gg == 2 * g1 - 1 - s}")
# 4. arXiv v1 family H_r
def H(r):
    E = [(a, b) for a in range(r) for b in range(a + 1, r)] + [(0, r)]
    for j in range(r): E += [(r, r + 1 + j), (r + 1 + j, 2 * r + 1 + j)]
    adj = [0] * (3 * r + 1)
    for a, b in E: adj[a] |= 1 << b; adj[b] |= 1 << a
    return 3 * r + 1, adj
for r in (2, 3, 4):
    n1, a1 = H(r); n, g = glue(n1, a1, n1, a1, r); g1 = gst(n1, a1); gg = gst(n, g)
    deg = [bin(x).count('1') for x in a1]; psi = sum(deg[v] for v in range(r))
    v1rhs = 2 * g1 + 2 * r * r - 2 * r + 1 - 2 * psi; p1rhs = 2 * g1 - 1 - 2 * sigma(n1, a1, r)
    ok &= gg < v1rhs and gg >= p1rhs
    print(f"v1 family H_{r}: gamma_st(H)={g1} gamma_st(G)={gg} Psi={psi} v1_bound={v1rhs} violated={gg < v1rhs} (4)_bound={p1rhs}")
print("RESULT", "PASS" if ok else "FAIL")
raise SystemExit(0 if ok else 1)
