# Derived from the closed unit24259 readiness program.
# Only image, oracle and Encoding definitions retained; no campaign entry point.
# SPDX-License-Identifier: Apache-2.0
"""Exact bounded binary automaton synthesis by counterexample-guided refinement.

Original target: Wang--Dong, arXiv:2610.39833v1, Conjecture 6.
The fixed n,m,s,r,S,A and witness length limits restrict candidate search only.
No finite UNSAT result establishes the conjecture. Solver answers are untrusted.
"""
import argparse, collections, hashlib, json, math, resource, time
from pathlib import Path
import z3


def image(table, mask, letter):
    out = 0
    while mask:
        bit = mask & -mask
        out |= 1 << table[letter][bit.bit_length()-1]
        mask -= bit
    return out


def apply(table, mask, word):
    for a in word:
        mask = image(table, mask, a)
    return mask


def strong(table):
    n = len(table[0])
    for reverse in (False, True):
        seen, todo = {0}, [0]
        while todo:
            q = todo.pop()
            nxt = ([x for x in range(n) if any(table[a][x] == q for a in (0, 1))]
                   if reverse else [table[a][q] for a in (0, 1)])
            for x in nxt:
                if x not in seen:
                    seen.add(x); todo.append(x)
        if len(seen) != n:
            return False
    return True


def oracle(table, S, A, B, certificate=False):
    """Complete image BFS, including epsilon. Merging equal subsets is exact.

    A violation is (A not subset image) OR (image size < |S|).
    Determinism makes future images depend only on the current subset.
    """
    n, s = len(table[0]), S.bit_count()
    assert len(table) == 2 and all(len(row) == n for row in table)
    assert all(0 <= t < n for row in table for t in row)
    assert B >= 0 and A and A != S and A & S == A
    queue = collections.deque([S]); records = {S: (0, ())}; levels = [0]*(B+1)
    edges = []
    while queue:
        mask = queue.popleft(); depth, word = records[mask]; levels[depth] += 1
        if mask & A != A or mask.bit_count() < s:
            return {'survives': False, 'word': list(word), 'mask': mask,
                    'visited': len(records), 'processed': sum(levels), 'levels': levels}
        if depth == B:
            continue
        for a in (0, 1):
            nxt = image(table, mask, a)
            if certificate: edges.append([mask, a, nxt])
            if nxt not in records:
                records[nxt] = (depth+1, word+(a,)); queue.append(nxt)
    result = {'survives': True, 'word': None, 'visited': len(records),
              'processed': sum(levels), 'levels': levels}
    if certificate:
        result['nodes'] = [[mask, d] for mask, (d, _) in sorted(records.items())]
        result['edges'] = edges
    return result


class Encoding:
    def __init__(self, n, S, A, B, reach_bound, reset_bound, timeout_ms=20000):
        self.n, self.S, self.A, self.B = n, S, A, B
        self.reach_bound, self.reset_bound = reach_bound, reset_bound
        self.bits = max(1, (n-1).bit_length()); self.serial = 0
        self.solver = z3.Solver(); self.solver.set(timeout=timeout_ms)
        self.declared = collections.Counter(); self.cuts = []
        self.table = [[self.bv(f't_{a}_{q}', self.bits) for q in range(n)] for a in (0, 1)]
        for row in self.table:
            for t in row:
                if n < (1 << self.bits):
                    self.solver.add(z3.ULT(t, z3.BitVecVal(n, self.bits)))
        self.onehots = [[z3.BitVecVal(1,n) << z3.ZeroExt(n-self.bits,t)
                         for t in row] for row in self.table]
        # Exact root reachability and reverse reachability in <= n-1 edges.
        for reverse in (False, True):
            prev = z3.BitVecVal(1,n)
            for step in range(n-1):
                nxt = self.bv(f'sc_{int(reverse)}_{step}',n)
                imgs = []
                for a in (0, 1):
                    if reverse:
                        terms = [z3.If(z3.Extract(0,0,z3.LShR(prev,z3.ZeroExt(n-self.bits,self.table[a][q]))) == 1,
                                       z3.BitVecVal(1 << q,n),z3.BitVecVal(0,n)) for q in range(n)]
                        imgs.append(self.union(terms))
                    else:
                        imgs.append(self.img(prev,a))
                self.solver.add(nxt == (prev | imgs[0] | imgs[1])); prev = nxt
            self.solver.add(prev == z3.BitVecVal((1 << n)-1,n))
        self.reach_images, self.reach_letters = self.word_path('reach',reach_bound)
        self.solver.add(z3.Or([v == z3.BitVecVal(S,n) for v in self.reach_images]))
        self.reset_images, self.reset_letters = self.word_path('reset',reset_bound)
        last = self.reset_images[-1]
        # Deterministic total images are nonempty. Bit trick means singleton.
        self.solver.add(last != 0, last & (last-1) == 0)
        self.add_cut(())

    def bv(self,name,width):
        self.declared[f'BV{width}'] += 1
        return z3.BitVec(name,width)

    def boolean(self,name):
        self.declared['Bool'] += 1
        return z3.Bool(name)

    @staticmethod
    def union(terms):
        out = terms[0]
        for term in terms[1:]: out = out | term
        return out

    def img(self, mask, letter):
        return self.union([z3.If(z3.Extract(q,q,mask) == 1,self.onehots[letter][q],z3.BitVecVal(0,self.n))
                           for q in range(self.n)])

    def word_path(self,prefix,bound):
        prev = z3.BitVecVal((1 << self.n)-1,self.n)
        images, letters = [prev], []
        for i in range(bound):
            a = self.boolean(f'{prefix}_a_{i}'); nxt = self.bv(f'{prefix}_s_{i}',self.n)
            self.solver.add(nxt == z3.If(a,self.img(prev,1),self.img(prev,0)))
            letters.append(a); images.append(nxt); prev = nxt
        return images,letters

    def add_cut(self,word):
        mask = z3.BitVecVal(self.S,self.n)
        for a in word: mask = self.img(mask,a)
        # Cardinality == s, not merely preservation of A. Includes epsilon.
        self.solver.add(mask & z3.BitVecVal(self.A,self.n) == self.A)
        self.solver.add(z3.Sum([z3.If(z3.Extract(q,q,mask) == 1,1,0) for q in range(self.n)]) == self.S.bit_count())
        self.cuts.append(list(word))

    def decode(self,model):
        table = [[model.eval(t,model_completion=True).as_long() for t in row] for row in self.table]
        rw = [int(z3.is_true(model.eval(a,model_completion=True))) for a in self.reach_letters]
        zw = [int(z3.is_true(model.eval(a,model_completion=True))) for a in self.reset_letters]
        full = (1 << self.n)-1
        reaches = [i for i in range(len(rw)+1) if apply(table,full,rw[:i]) == self.S]
        resets = [i for i in range(len(zw)+1) if apply(table,full,zw[:i]).bit_count() == 1]
        assert strong(table) and reaches and resets
        return {'table':table,'reaching_word':rw[:reaches[0]],'reset_word':zw[:resets[0]]}

    def counts(self):
        return {'declared_variables':dict(self.declared),'total_declared_variables':sum(self.declared.values()),
                'top_level_assertions':len(self.solver.assertions()),'no_success_cuts':len(self.cuts),
                'cut_word_lengths':[len(w) for w in self.cuts]}

