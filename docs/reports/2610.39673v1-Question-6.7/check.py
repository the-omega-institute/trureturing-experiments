#!/usr/bin/env python3
"""Exact finite checks and seeded samples for BGS Question 6.7 (no Lean proof).
All reported values are integers, rational numbers, or certified rational intervals.
The congruence closure derives only original-presentation equalities; it adds no zero
or truncation relations. A finite answer is accepted only after constructing and
checking a complete associative multiplication table satisfying the presentation.
"""
from collections import Counter, deque
from fractions import Fraction as Q
from itertools import combinations, product
from math import comb, factorial
from random import Random


def words(n, h):
    return [w for k in range(1, h + 1) for w in product(range(n), repeat=k)]


def factors(w):
    return {w[i:j] for i in range(len(w)) for j in range(i + 1, len(w) + 1)}


def choose(a, b):
    return comb(a, b) if 0 <= b <= a else 0


def avoid(N, m, k):
    return Q(choose(choose(N-k, 2), m), choose(choose(N, 2), m))


def pattern(N, m, absent, present):
    return sum((-1)**j * comb(present, j) * avoid(N, m, absent+j)
               for j in range(present+1))


def exp_negative_bounds(x, terms=80):
    # For x>=0 and terms+2>x, exp(x) has a positive-term geometric tail bound.
    partial = sum((x**j / factorial(j) for j in range(terms+1)), Q(0))
    nxt = x**(terms+1) / factorial(terms+1)
    tail = nxt / (1-x/Q(terms+2))
    return 1/(partial+tail), 1/partial


def rational_interval_text(lo, hi, digits=9):
    scale = 10**digits
    return f'[{lo.numerator*scale//lo.denominator}/{scale}, ' \
           f'{-((-hi.numerator*scale)//hi.denominator)}/{scale}]'


def finite_checks():
    count = 0
    # Source Section 2: exact ball and unordered-relation counts.
    for n in range(2, 6):
        for l in range(1, 5):
            B = words(n, l)
            assert len(B) == (n**(l+1)-n)//(n-1)
            assert len(list(combinations(B, 2))) == comb(len(B), 2)
            count += 1
    print('FIDELITY source ball and unordered-pair identities:', count, 'cases')
    # Source Proposition 6.6's hypergeometric avoidance formula, and the proposed
    # joint endpoint law, independently checked by enumerating every edge subset.
    cases = 0
    for N in range(3, 7):
        E = list(combinations(range(N), 2))
        for m in range(1, min(4, len(E))+1):
            hist = Counter()
            for R in combinations(E, m):
                endpoints = {u for edge in R for u in edge}
                hist[tuple(u in endpoints for u in range(min(N, 3)))] += 1
            total = comb(len(E), m)
            for mask in product([False, True], repeat=min(N, 3)):
                assert Q(hist[mask], total) == pattern(N, m, mask.count(False), mask.count(True))
                cases += 1
    print('FIDELITY source Proposition 6.6 avoidance and joint endpoint patterns:', cases, 'exact cases')
    # Finite factorial-language atom identity as an exact polynomial in q:
    # enumerate every mark set on all words of length <=3 over two letters.
    U = words(2, 3)
    hist = Counter()
    for bits in range(1 << len(U)):
        marks = {U[i] for i in range(len(U)) if bits >> i & 1}
        F = frozenset(w for w in U if factors(w).isdisjoint(marks))
        hist[F, len(marks)] += 1
    checked = 0
    for F in {F for F, _ in hist}:
        if any(len(w) == 3 for w in F):
            continue  # Only a genuinely closed finite language, not a truncated one.
        boundary = {w for w in U if w not in F and (factors(w)-{w}) <= F}
        for q in [Q(1, 4), Q(1, 2), Q(2, 3)]:
            actual = sum(count * (1-q)**k * q**(len(U)-k)
                         for (G, k), count in hist.items() if G == F)
            assert actual == q**len(F) * (1-q)**len(boundary)
            checked += 1
    print('SETTLEMENT finite factorial-language mass formula:', checked, 'exact evaluations')
    # Exact finite-pattern limits are numerically fenced with rational exp bounds.
    qlo, qhi = exp_negative_bounds(Q(4, 3))
    plo, phi = (1-qhi)**3, (1-qlo)**3
    print('LIMIT n=3 q=exp(-4/3):', rational_interval_text(qlo, qhi))
    print('LIMIT n=3 triviality=(1-exp(-4/3))^3:', rational_interval_text(plo, phi))
    print('ENDPOINT exact all-three-generators-marked probabilities:')
    for l in range(2, 7):
        N = sum(3**j for j in range(1, l+1))
        v = pattern(N, 3**l, 0, 3)
        print(' n=3 l=', l, 'probability=', v, sep='')
    return plo, phi


class Closure:
    def __init__(self, n, h, rels):
        self.n, self.h = n, h
        self.ws = words(n, h)
        self.index = {w:i for i,w in enumerate(self.ws)}
        N = len(self.ws)
        self.parent = list(range(N))
        self.size = [1]*N
        self.transitions = []
        for w in self.ws:
            self.transitions.append([self.index[(a,)+w] for a in range(n)] +
                                    [self.index[w+(a,)] for a in range(n)]
                                    if len(w) < h else [-1]*(2*n))
        queue = deque((self.index[u], self.index[v]) for u,v in rels)
        while queue:
            a,b = queue.popleft(); a,b = self.root(a),self.root(b)
            if a == b: continue
            if self.size[a] < self.size[b]: a,b = b,a
            self.parent[b] = a; self.size[a] += self.size[b]
            for k,(u,v) in enumerate(zip(self.transitions[a], self.transitions[b])):
                if u == -1: self.transitions[a][k] = v
                elif v != -1: queue.append((u,v))
        self.rels = rels

    def root(self, x):
        while self.parent[x] != x:
            self.parent[x] = self.parent[self.parent[x]]; x = self.parent[x]
        return x

    def finite_table(self):
        reps = {}
        for i,w in enumerate(self.ws): reps.setdefault(self.root(i), w)
        roots = list(reps)
        # No assertion of infinitude is made when this finite certificate fails.
        if any(-1 in self.transitions[r] for r in roots) or len(roots)>100:
            return None
        pos = {r:i for i,r in enumerate(roots)}
        def append(r,w):
            for a in w: r = self.root(self.transitions[r][self.n+a])
            return r
        T = [[pos[append(a,reps[b])] for b in roots] for a in roots]
        K = len(T)
        assert all(T[T[a][b]][c] == T[a][T[b][c]]
                   for a in range(K) for b in range(K) for c in range(K))
        gens = [pos[self.root(self.index[(a,)])] for a in range(self.n)]
        def evaluate(w):
            ans = gens[w[0]]
            for a in w[1:]: ans = T[ans][gens[a]]
            return ans
        assert all(evaluate(w) == pos[self.root(i)] for i,w in enumerate(self.ws))
        assert all(evaluate(u) == evaluate(v) for u,v in self.rels)
        return T, reps, roots, pos


def simulation(plo, phi):
    rng = Random(261039673)
    trials = 160
    reductions = 0
    for l in range(2, 7):
        B = words(3,l); N = len(B); m = 3**l
        # Sampling distinct unordered edges without constructing the complete graph.
        pairs = [(a,b) for a in range(N) for b in range(a+1,N)] if l <= 3 else None
        hist = Counter(); unresolved = 0; trivial = 0; marked_all = 0; mismatches = 0
        for _ in range(trials):
            if pairs is not None:
                E = rng.sample(pairs,m)
            else:
                E = set()
                while len(E)<m:
                    a,b = rng.sample(range(N),2); E.add((min(a,b),max(a,b)))
            R = [(B[a],B[b]) for a,b in E]
            endpoints = {w for edge in R for w in edge}
            marked = all((a,) in endpoints for a in range(3))
            marked_all += marked
            C = Closure(3,l+1,R)
            cert = C.finite_table()
            if cert is None:
                unresolved += 1; continue
            T,reps,roots,pos = cert
            K=len(T); hist[K]+=1; trivial += K == 1
            mismatches += ((K == 1) != marked)
            zeros = [z for z in range(K) if all(T[z][a] == z == T[a][z] for a in range(K))]
            if len(zeros) != 1: continue
            z=zeros[0]
            # Independently check the deterministic reduction on every sampled case
            # satisfying its hypotheses, without treating failed hypotheses as failure.
            r = next((r for r in range(1,l+1) if all(pos[C.root(C.index[w])] == z
                      for w in C.ws if len(w)==r)), None)
            if r is None or any(len(u)<r and len(v)<r for u,v in R): continue
            allowed = [w for w in C.ws if len(w)<r and factors(w).isdisjoint(endpoints)]
            assert len(allowed)+1 == K
            assert len({pos[C.root(C.index[w])] for w in allowed}) == len(allowed)
            assert all(pos[C.root(C.index[w])] != z for w in allowed)
            for u in allowed:
                for v in allowed:
                    want = pos[C.root(C.index[u+v])] if u+v in allowed else z
                    assert T[pos[C.root(C.index[u])]][pos[C.root(C.index[v])]] == want
            reductions += 1
        low,high = Q(trivial,trials),Q(trivial+unresolved,trials)
        print('SIM n=3 l=',l,' trials=',trials,' finite_size_hist=',dict(sorted(hist.items())),
              ' unresolved=',unresolved,' certified_trivial_probability_interval=',(str(low),str(high)),
              ' all_generators_marked=',str(Q(marked_all,trials)),
              ' finite_trivial_vs_marked_mismatches=',mismatches,sep='')
        # These are observations, not a statistical convergence claim.
    print('SETTLEMENT deterministic monomial reduction checked on',reductions,'certified samples')


if __name__ == '__main__':
    plo,phi = finite_checks()
    simulation(plo,phi)
    print('ALL_CHECKS_OK')
