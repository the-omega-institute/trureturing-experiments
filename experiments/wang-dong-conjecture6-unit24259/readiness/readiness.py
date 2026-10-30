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


def selfcheck():
    # Complete cross-check against literal word enumeration on every binary
    # two-state table, all nonempty proper A in S=Q, and B=0..3.
    import itertools
    checks = 0
    for entries in itertools.product(range(2),repeat=4):
        table = [list(entries[:2]),list(entries[2:])]
        for A in (1,2):
            for B in range(4):
                brute = all((lambda v: (v & A == A and v.bit_count() == 2))(apply(table,3,w))
                            for length in range(B+1) for w in itertools.product((0,1),repeat=length))
                assert oracle(table,3,A,B)['survives'] == brute
                checks += 1
    # Explicit attainable boundary: a is a 4-cycle; b maps 3 to 0.
    # S={0,1,2}, A={0}: shortest success is one letter (avoidance), not >1.
    table = [[1,2,3,0],[0,1,2,0]]
    o = oracle(table,7,1,1)
    assert strong(table) and o['word'] == [0] and not o['survives']
    # Valid small exact synthesis with every bounded word imposed up front.
    enc = Encoding(3,3,1,0,3,5,5000)
    result = enc.solver.check(); assert result == z3.sat,result
    model = enc.decode(enc.solver.model()); assert oracle(model['table'],3,1,0)['survives']
    # Fix all two-state tables and compare symbolic image and strong-connectivity
    # encoding to direct computations; no witness claims for non-SC tables.
    symbolic_checks = 0
    for entries in itertools.product(range(2),repeat=4):
        enc = Encoding(2,3,1,0,0,3,2000)
        table = [list(entries[:2]),list(entries[2:])]
        for a in (0,1):
            for q in range(2): enc.solver.add(enc.table[a][q] == table[a][q])
        expected = strong(table) and any(apply(table,3,w).bit_count() == 1
                         for length in range(4) for w in itertools.product((0,1),repeat=length))
        got = enc.solver.check(); assert got != z3.unknown
        assert (got == z3.sat) == expected
        for mask in range(1,4):
            for a in (0,1):
                expr = enc.img(z3.BitVecVal(mask,2),a)
                substitutions = [(enc.table[x][q],z3.BitVecVal(table[x][q],enc.bits)) for x in (0,1) for q in range(2)]
                assert z3.simplify(z3.substitute(expr,*substitutions)).as_long() == image(table,mask,a)
                symbolic_checks += 1
    return {'oracle_literal_comparisons':checks,'symbolic_image_comparisons':symbolic_checks,
            'fixed_table_encoding_comparisons':16,'small_synthesis':'sat, decoded and replayed',
            'avoidance_boundary_control':'ac=1, not strict ac>1'}


def main():
    p = argparse.ArgumentParser(); p.add_argument('--out',required=True); p.add_argument('--seconds',type=float,default=80)
    p.add_argument('--smt',required=True); args = p.parse_args()
    z3.set_param('memory_max_size',1536)
    start = time.monotonic(); result = {'solver':'Z3 '+z3.get_version_string(),'selfcheck':selfcheck()}
    t = time.monotonic(); enc = Encoding(19,(1 << 13)-1,1,17,126,220,20000)
    result['encoding_build_seconds'] = time.monotonic()-t
    result['initial_counts'] = enc.counts()
    t = time.monotonic(); text = enc.solver.sexpr(); Path(args.smt).write_text(text)
    result['smtlib_bytes'] = len(text.encode()); result['smtlib_sha256'] = hashlib.sha256(text.encode()).hexdigest()
    result['serialization_seconds'] = time.monotonic()-t
    del text
    result['checks'] = []; result['candidate'] = None
    while time.monotonic()-start < args.seconds:
        t = time.monotonic(); status = enc.solver.check()
        check = {'solver_status':str(status),'solver_seconds':time.monotonic()-t,
                 'statistics':{str(k):v for k,v in enc.solver.statistics()},'counts':enc.counts()}
        result['checks'].append(check)
        if status == z3.unknown:
            check['reason'] = enc.solver.reason_unknown(); break
        if status == z3.unsat: break
        model = enc.decode(enc.solver.model()); t = time.monotonic()
        checked = oracle(model['table'],enc.S,enc.A,enc.B,True)
        check['oracle_seconds'] = time.monotonic()-t
        check['oracle'] = {k:v for k,v in checked.items() if k not in ('edges','nodes')}
        if checked['survives']:
            result['candidate'] = dict(model,n=19,m=6,s=13,r=1,S=list(range(13)),A=[0],B=17,certificate=checked)
            break
        enc.add_cut(checked['word'])
    result['final_counts'] = enc.counts(); result['elapsed_seconds'] = time.monotonic()-start
    result['peak_rss_bytes_macos'] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    result['complete_constraint_system'] = 'bounded admissibility plus lazily separated exact all-word no-success constraints; readiness run may end before closure'
    Path(args.out).write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'out':args.out,'elapsed_seconds':result['elapsed_seconds'],'candidate':result['candidate'] is not None,'checks':len(result['checks'])}),flush=True)

if __name__ == '__main__': main()
