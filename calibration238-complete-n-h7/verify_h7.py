#!/usr/bin/env python3
"""Offline exact verification of the normalized seven-coordinate certificate.

Python standard library only. No solver, generator, network, source repository,
or stored geometry trace is used. See README.md for the precise proof scope.
"""
import argparse
from fractions import Fraction as F
from functools import reduce
from itertools import combinations, product
from math import gcd, lcm
from pathlib import Path
import hashlib
import json


def require(ok, message):
    if not ok:
        raise ValueError(message)


def primitive(v):
    d = reduce(gcd, v, 0)
    require(d != 0, 'zero homogeneous vector')
    return tuple(x // abs(d) for x in v)


def homogeneous(x):
    q = lcm(*(t.denominator for t in x))
    return primitive(tuple(int(t*q) for t in x) + (q,))


def rank(rows):
    """Exact integer row elimination, with primitive reduction at each step."""
    pivots = {}
    for source in rows:
        r = tuple(source)
        for j, p in sorted(pivots.items()):
            if r[j]:
                a, b = p[j], r[j]
                r = tuple(a*x-b*y for x, y in zip(r, p))
                if any(r):
                    r = primitive(r)
        if any(r):
            j = next(j for j, t in enumerate(r) if t)
            pivots[j] = r
    return len(pivots)


def solve(a, b):
    m = [list(map(F, r)) + [F(t)] for r, t in zip(a, b)]
    n = len(b)
    for j in range(n):
        k = next((k for k in range(j, n) if m[k][j]), None)
        if k is None:
            return None
        m[j], m[k] = m[k], m[j]
        c = m[j][j]
        m[j] = [v/c for v in m[j]]
        for k in range(n):
            if k != j and m[k][j]:
                c = m[k][j]
                m[k] = [v-c*w for v, w in zip(m[k], m[j])]
    return homogeneous([r[-1] for r in m])


def slack(a, b, v):
    return b*v[-1] - sum(x*y for x, y in zip(a, v))


def rows_for(L):
    words = list(product((0, 1), repeat=L))
    rows = []
    for wi, w in enumerate(words):
        if 0 in w:
            n = w.index(0)+1
            rows.extend([(wi, 'A', n, {0}), (wi, 'C', n, {n}),
                         (wi, 'P', n, {0, n})])
        else:
            rows.append((wi, 'T', 0, {0}))
    return words, rows


def score_check(weights, words, rows, pairs):
    """Reconstruct the full event reward and every independent N cell maximum."""
    reward = [[F(0) for _ in range(8)] for _ in words]
    beta, K, const = [F(0)]*7, F(0), F(0)
    seen = set()
    for ri, s in weights:
        require(type(ri) is int and 0 <= ri < len(rows) and ri not in seen,
                'invalid or repeated score row')
        seen.add(ri)
        t = F(s)
        require(t >= 0, 'negative score')
        wi, kind, n, event = rows[ri]
        K += t
        if kind == 'C':
            beta[n-1] += t
        elif kind == 'A':
            beta[n-1] -= t
        if kind in ('A', 'P', 'T'):
            const += t
        for z in event:
            reward[wi][z] += t
    maxima = []
    for i, j in pairs:
        cm = []
        for a, b in product((0, 1), repeat=2):
            cm.append(max(sum(reward[k][z] for k, w in enumerate(words)
                              if (w[i-1], w[j-1]) == (a, b))
                          for z in range(8)))
        maxima.append(cm)
    BN = max(map(sum, maxima))
    return beta, K, const, maxima, BN


def geometry(A, B, supplied):
    require(len(A) == len(B) == 481 and all(len(a) == 7 for a in A),
            'halfspace dimensions')
    require(all(type(v) is int for a in A for v in a) and
            all(type(v) is int for v in B), 'noninteger halfspace')
    domainA = [[1,0,0,0,0,0,0],[-1,1,0,0,0,0,0],
               [0,-1,1,0,0,0,0],[0,0,-1,1,0,0,0],
               [0,0,0,-1,1,0,0],[0,0,0,0,-1,1,0],
               [0,0,0,0,0,-1,1],[0,0,0,0,0,0,-1],
               [0,-1,-1,-1,-1,-1,8],[1,1,1,1,1,1,-8]]
    require(A[:10] == domainA and B[:10] == [1,0,0,0,0,0,0,0,0,1],
            'H7 domain differs from statement')
    vertices = {}
    bases = 0
    for I in combinations(range(10), 7):
        bases += 1
        v = solve([A[i] for i in I], [B[i] for i in I])
        if v is not None and all(slack(a, b, v) >= 0
                                 for a, b in zip(A[:10], B[:10])):
            mask = sum(1 << i for i in range(10) if slack(A[i], B[i], v) == 0)
            vertices[v] = mask
    initial = set(vertices)
    require(len(initial) == 12, 'initial domain vertices')
    expected = {homogeneous([F(0)]*7), homogeneous([F(1)]+[F(0)]*6)}
    for k in range(3, 8):
        b = [F(1,k-1)]*k + [F(0)]*(7-k)
        expected.add(homogeneous([b[0],b[1],*b[3:],F(0)]))
        e = F(k-2,k+1)
        b = [F(1)]*k + [e]*(7-k)
        expected.add(homogeneous([b[0],b[1],*b[3:],e]))
    require(initial == expected, 'initial domain vertex formulas')
    rank_cache = {}
    def edge(mask):
        if bin(mask).count('1') < 6:
            return False
        if mask not in rank_cache:
            indices = [i for i in range(len(A)) if (mask >> i) & 1]
            rank_cache[mask] = rank([A[i] for i in indices]) == 6
        return rank_cache[mask]
    crossings, peak = 0, len(vertices)
    for h in range(10, len(A)):
        inside, outside, on = [], [], []
        for v, m in vertices.items():
            s = slack(A[h], B[h], v)
            (inside if s > 0 else outside if s < 0 else on).append((v, m, s))
        new = {v: m for v, m, s in inside}
        new.update({v: m | (1 << h) for v, m, s in on})
        for v, vm, sv in outside:
            for u, um, su in inside:
                common = vm & um
                if edge(common):
                    crossings += 1
                    p = primitive(tuple(su*x-sv*y for x, y in zip(v, u)))
                    require(p[-1] > 0, 'intersection denominator')
                    pm = common | (1 << h)
                    if p in new:
                        require(new[p] == pm, 'inconsistent active set')
                    new[p] = pm
        vertices = new
        require(vertices, 'empty clipping stage')
        peak = max(peak, len(vertices))
    supplied_map = {}
    for rec in supplied:
        v = homogeneous(list(map(F, rec['x'])))
        require(v not in supplied_map, 'duplicate supplied vertex')
        active = [i for i, (a,b) in enumerate(zip(A,B)) if slack(a,b,v) == 0]
        require(all(slack(a,b,v) >= 0 for a,b in zip(A,B)), 'vertex containment')
        require(rec['active'] == active, 'supplied full active set')
        require(rank([A[i] for i in active]) == 7, 'nonvertex active rank')
        supplied_map[v] = sum(1 << i for i in active)
    require(vertices == supplied_map, 'incomplete or incorrect full vertex list')
    origin = next(iter(vertices))
    require(rank([tuple(v[i]*origin[-1]-origin[i]*v[-1] for i in range(7))
                  for v in vertices]) == 7, 'final affine dimension')
    final_edges = sum(edge(vertices[u] & vertices[v]) for u,v in combinations(vertices,2))
    return {'domain_bases': bases, 'domain_vertices':len(initial),
            'cuts_clipped':len(A)-10, 'crossing_edges':crossings,
            'peak_vertices':peak, 'vertices':len(vertices), 'affine_dimension':7,
            'final_edges':final_edges,
            'vertices_with_exactly_seven_active_rows':sum(bin(m).count('1')==7 for m in vertices.values())}


def lifts_check(c, words):
    require(len(c['lifts']) == len(c['vertices']), 'lift cardinality')
    counts = {'lifts':0, 'schedule_rows':0, 'cell_vectors':0,
              'word_laws':0, 'independent_event_rows':0}
    for vi, lift in enumerate(c['lifts']):
        require(lift['vertex_index'] == vi and lift['x'] == c['vertices'][vi]['x'],
                'lift vertex binding')
        x = list(map(F, lift['x']))
        b, e = [x[0],x[1],x[1],*x[2:6]], x[6]
        schedules = []
        for rec in lift['schedules']:
            i,j = rec['pair']
            require(type(i) is int and type(j) is int and 1 <= i < j <= 7,
                    'schedule addresses')
            mass = F(rec['mass'])
            require(mass >= 0 and len(rec['cells']) == 3, 'schedule mass/cells')
            cells = []
            for sparse in rec['cells']:
                v, seen = [F(0)]*8, set()
                for z,s in sparse:
                    require(type(z) is int and 0 <= z <= 7 and z not in seen,
                            'cell label')
                    seen.add(z); v[z] = F(s)
                require(all(t >= 0 for t in v) and sum(v) == mass,
                        'counterfactual cell budget')
                cells.append(v)
            schedules.append((i,j,mass,cells))
        require(sum(m for i,j,m,cs in schedules) == 1, 'global mass')
        for w in words:
            p = [F(0)]*8
            for i,j,m,cs in schedules:
                v = cs[0 if w[i-1] == 0 else 1 if w[j-1] == 0 else 2]
                p = [a+t for a,t in zip(p,v)]
            require(sum(p) == 1 and min(p) >= 0, 'word probability')
            if 0 in w:
                n = w.index(0)+1
                require(p[0] >= 1-b[n-1]-e, 'base row')
                require(p[n] >= b[n-1]-e, 'correct row')
                require(p[0]+p[n] >= 1-e, 'independent union row')
                counts['independent_event_rows'] += 3
            else:
                require(p[0] >= 1-e, 'positive row')
                counts['independent_event_rows'] += 1
            counts['word_laws'] += 1
        counts['lifts'] += 1
        counts['schedule_rows'] += len(schedules)
        counts['cell_vectors'] += 3*len(schedules)
    return counts


def example_check(ex, words, rows, pairs):
    low = ex['lower']
    weights = []
    for typ,word,t in low['typed_nonzero_rows']:
        wi = words.index(tuple(map(int,word)))
        ri = next(k for k,r in enumerate(rows) if r[0]==wi and r[1]==typ)
        weights.append((ri,str(t)))
    beta,K,const,cm,BN = score_check(weights,words,rows,pairs)
    require(beta == low['beta'] and K == low['K'] and const == low['constant_demand']
            and cm == low['N_four_cell_maxima'] and BN == low['B_N'], 'example full N score')
    amp = ex['amplitudes_numerators_over_64']
    require(sum(a*t for a,t in zip(amp,beta)) == 694 and BN-const == 1
            and F(694-64,64*42)==F(ex['optimal_error']), 'example lower equality')
    correct, outside, base = [set() for _ in range(7)], [0]*7, [-1000]*7
    schedules = ex['upper']['schedules']; d = ex['upper']['denominator']
    require(d == 64 and len(schedules)==14 and sum(r[2] for r in schedules)==d,
            'example schedule law')
    for i,j,m,*cs in schedules:
        require(1<=i<j<=7 and m>=0 and len(cs)==3 and
                all(sum(v.values())==m and min(v.values())>=0 for v in cs), 'example cells')
    laws = {}
    for w in words:
        p=[0]*8
        for i,j,m,*cs in schedules:
            for z,t in cs[0 if w[i-1]==0 else 1 if w[j-1]==0 else 2].items():
                p[int(z)]+=t
        laws[w] = p
        require(sum(p)==64,'example word mass')
        if 0 in w:
            n=w.index(0)
            correct[n].add(p[n+1]);outside[n]=max(outside[n],64-p[0]-p[n+1])
            base[n]=max(base[n],64-amp[n]-p[0])
            require(p[0]>=64-amp[n]-15 and p[n+1]>=amp[n]-15
                    and p[0]+p[n+1]>=49, 'example event row')
        else:
            require(p==ex['upper']['positive_law_counts'],'example positive law')
    zero = (0,)*7
    c0 = laws[zero]
    ell, Q = [], {}
    for i in range(7):
        w = tuple(int(j==i) for j in range(7))
        ell.append([a-b for a,b in zip(laws[w],c0)])
    for i,j in combinations(range(7),2):
        w = tuple(int(k in (i,j)) for k in range(7))
        Q[i,j] = [laws[w][z]-c0[z]-ell[i][z]-ell[j][z] for z in range(8)]
    H = [0]*8
    for i in range(7):
        V = [sum(max(0,-Q[i,j][z]) for j in range(i+1,7)) for z in range(8)]
        H = [a+max(0,v-l) for a,v,l in zip(H,V,ell[i])]
    require(c0 == ex['upper']['criterion109']['c'] and
            H == ex['upper']['criterion109']['H'] == c0, 'example coefficient criterion')
    require([sorted(v) for v in correct]==ex['upper']['correct_counts'] and
            outside==ex['upper']['max_outside_pair_counts'] and
            base==ex['upper']['max_base_deficits'], 'example fibre bounds')
    return {'schedules':14,'counterfactual_cells':42,'words':128,'event_rows':382,
            'score_rows':22,'N_pairs':21,'cell_label_maximizations':21*4*8,
            'optimal_error':'15/64','weighted_amplitude_numerator_over_64':694}


def scalar_obstruction():
    words=list(product((0,1),repeat=3));f=lambda w:(1-w[0])*w[1]+w[0]*w[2]
    total=avoid=best=0
    schedules=[()]+[(i,) for i in range(3)]+list(combinations(range(3),2))
    for S in schedules:
        for table in product(range(4),repeat=2**len(S)):
            total+=1
            def out(w):
                a=0
                for i in S:a=2*a+w[i]
                return table[a]
            if all(out(w)!=0 for w in words if f(w)):
                avoid+=1;best=max(best,sum(out(w)==0 for w in [(0,0,1),(1,1,0)]))
    require((total,avoid,best)==(820,354,1),'scalar obstruction enumeration')
    return {'deterministic_descriptions':total,'avoid_forbidden_zero_labels':avoid,
            'maximum_two_word_zero_score':best,'target_score':'3/2'}


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--certificate',required=True,type=Path)
    p.add_argument('--example',required=True,type=Path)
    p.add_argument('--output',required=True,type=Path)
    args=p.parse_args()
    raw=args.certificate.read_bytes();c=json.loads(raw)
    require(c['schema']=='research161-H7-complete-normalized-common-allocation-v1'
            and c['L']==7 and c['cell_order_N']==['00','01','10','11']
            and c['cell_order_O']==['0','10','11'], 'certificate scope')
    require(c['reduced_coordinates']==['b1','b2=b3','b4','b5','b6','b7','e'] and
            c['row_order']=='lexicographic seven-bit words, A,C,P at each zero word; T at 1111111',
            'coordinate/row convention')
    words,rows=rows_for(7);pairs=list(combinations(range(1,8),2))
    require(len(c['cuts'])==471 and len(c['vertices'])==1606,'certificate sizes')
    for h,cut in enumerate(c['cuts'],10):
        beta,K,const,cm,BN=score_check(cut['weights'],words,rows,pairs)
        require(beta==cut['beta'] and K==cut['K'] and const==cut['const']
                and cm==cut['cell_maxima'] and BN==cut['BN'], 'cut full-catalog support')
        unscaled=[beta[0],beta[1]+beta[2],*beta[3:],-K,BN-const]
        scale=lcm(*(t.denominator for t in unscaled))
        v=primitive(tuple(int(t*scale) for t in unscaled))
        require(list(v[:-1])==cut['normal']==c['A'][h] and v[-1]==cut['rhs']==c['B'][h],
                'cut positive primitive normalization')
    result={'certificate_sha256':hashlib.sha256(raw).hexdigest(),
            'certificate_bytes':len(raw),'arithmetic':'exact integers and fractions',
            'scores':{'cuts':471,'N_pairs_per_cut':21,'cells_per_pair':4,
                      'labels_per_cell':8,'event_rows':len(rows),'words':len(words)},
            'geometry':geometry(c['A'],c['B'],c['vertices']),
            'lifts':lifts_check(c,words),
            'example':example_check(json.loads(args.example.read_text()),words,rows,pairs),
            'scalar_budget_obstruction':scalar_obstruction()}
    args.output.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(json.dumps(result,sort_keys=True))

if __name__=='__main__':
    main()
