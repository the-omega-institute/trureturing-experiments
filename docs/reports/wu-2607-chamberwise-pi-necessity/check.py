#!/usr/bin/env python3
"""Exact checks for Wu's global/chamber-wise PI distinction and an MID obstruction.

Only rational matrices, polynomial coefficients and finite incidence data are
compared. Equal symbolic Gaussian kernels certify density formulas; this program
is not a proof of equality of pushforward measures or mixture identifiability.
"""
from fractions import Fraction as F
from collections import Counter, defaultdict
from itertools import product
import sympy as sp

I = sp.eye(2)
N = -I
R = sp.diag(1, -1)
S = sp.diag(-1, 1)
rules = {'I': I, 'N': N, 'R': R, 'S': S}
means = (sp.Matrix([1, 0]), sp.Matrix([0, 1]))
covs = (sp.diag(1, 2), sp.diag(3, 4))
weights = (F(1, 3), F(2, 3))
cut = (-4, -3, -1, 0, 1, 3, 4)
intervals = tuple(zip((None,) + cut, cut + (None,)))
samples = tuple(F(b - 1) if a is None else F(a + 1) if b is None
                else F(a + b, 2) for a, b in intervals)
cells = tuple(product(range(8), repeat=2))

def theta(A, k):
    return tuple(A * means[k]), tuple(A * covs[k] * A.T)

def cell_of(z):
    return tuple(sum(int(bool(z[c] > a)) for a in cut) for c in range(2))

def rule_at(which, x, y):
    if -1 < x < 1 and 0 < y < 1:
        return 'I' if which == 'f' else 'S'
    if -1 < x < 1 and -1 < y < 0:
        return 'N' if which == 'f' else 'R'
    if 3 < x < 4 and -1 < y < 1:
        return 'R'
    if -1 < x < 1 and 3 < y < 4:
        return 'S'
    if (3 < x < 4 and 3 < y < 4) or (-4 < x < -3 and -4 < y < -3):
        return 'N'
    return 'I'

# Fidelity: Appendix E, Proposition 18, d=1, sigma^2=1,
# mu_1=1, mu_2=-9, rules z and z+10 on the two half-lines.
assert 1 == -9 + 10
source_images = ((None, F(0)), (F(10), None))
assert source_images[0][1] < source_images[1][0]
source_parameters = ((1, 1), (-9, 1), (11, 1), (1, 1))
assert len(set(source_parameters)) == 3
assert len(set(source_parameters[:2])) == len(set(source_parameters[2:])) == 2
print('FIDELITY: Appendix E Proposition 18: injective; chamber-wise PI; global collision (1,1)=(2,2)')

active = {w: {c: set() for c in cells} for w in ('f', 'g')}
preimages = {w: {c: [] for c in cells} for w in ('f', 'g')}
assignment = {w: {} for w in ('f', 'g')}
for w in ('f', 'g'):
    for c in cells:
        z = sp.Matrix([samples[c[0]], samples[c[1]]])
        r = rule_at(w, *z)
        assignment[w][c] = r
        A = rules[r]
        assert abs(A.det()) == 1
        out = cell_of(A * z)
        # A is diagonal with signs, so its image of each entire open cell
        # is exactly the following cell, including both unbounded intervals.
        expected = tuple(c[j] if A[j, j] == 1 else 7 - c[j] for j in range(2))
        assert out == expected
        assert r not in active[w][out], 'same rule counted twice at a regular output'
        active[w][out].add(r)
        preimages[w][out].append((c, r))
    counts = Counter(assignment[w].values())
    assert set(counts) == set(rules)
    print('REDUCED', w, '64 input cells; four distinct nonempty rules;', dict(sorted(counts.items())))

for k, D in enumerate(covs):
    assert D == D.T and D[0, 0] > 0 and D.det() > 0
assert theta(I, 0) != theta(I, 1)
assert sum(weights) == 1 and all(a > 0 for a in weights)
print('GMM: positive weights; sum=1; distinct components; covariance principal minors positive')

x, y = sp.symbols('x y', real=True)
z = sp.Matrix([x, y])
def exponent(k, A):
    v = A*z-means[k]
    return sp.expand(-(v.T*covs[k].inv()*v)[0]/2)
for k, A, B in [(0,R,I),(0,S,N),(1,S,I),(1,R,N)]:
    assert sp.expand(exponent(k,A)-exponent(k,B)) == 0
print('DENSITY: four component-specific reflection identities, as exact quadratic polynomials')

pi_failures = {}
for w in ('f', 'g'):
    failures=[]
    for r, k in product(rules, range(2)):
        for s, ell in product(rules, range(2)):
            if (r,k) < (s,ell) and theta(rules[r],k)==theta(rules[s],ell):
                failures.append((r,k+1,s,ell+1))
    assert failures
    pi_failures[w]=failures
    for c in cells:
        pars=[theta(rules[r],k) for r in active[w][c] for k in range(2)]
        assert len(pars)==len(set(pars)), (w,c,pars)
    print('PI',w,': chamber-wise PI on all 64 cells; global collisions', failures)

for c in cells:
    signatures=[]
    for w in ('f','g'):
        sig=defaultdict(F)
        for r in active[w][c]:
            for k in range(2):sig[theta(rules[r],k)]+=weights[k]
        signatures.append(dict(sig))
    assert signatures[0]==signatures[1], (c,signatures)
    n=len(preimages['f'][c])
    assert n==len(preimages['g'][c])
    print('CELL',c,'sample',str(samples[c[0]]),str(samples[c[1]]),
          'f',','.join(sorted(active['f'][c])) or 'empty',
          'g',','.join(sorted(active['g'][c])) or 'empty',
          'preimage_count',n,'density_kernels_equal')
qplus=[c for c in cells if -1 < samples[c[0]] < 1 and 0 < samples[c[1]] < 1]
assert len(qplus)==2
assert all(len(preimages[w][c])==2 for w in ('f','g') for c in qplus)
print('NONINJECTIVE: every regular point of Q+ has two distinct interior preimages, for both maps')

# A second, measure-theoretic route: localize R for component 1 and S for
# component 2 to the invariant square Q. Check g = f composed with the
# respective localized symmetry as matrix identities on every input cell.
for k, reflection in enumerate((R,S)):
    for c0 in cells:
        z0=sp.Matrix([samples[c0[0]],samples[c0[1]]])
        in_square=all(-1 < z0[j] < 1 for j in range(2))
        T=reflection if in_square else I
        At=rules[rule_at('f',*(T*z0))]
        Ag=rules[rule_at('g',*z0)]
        assert At*T==Ag
    print('LOCALIZED SYMMETRY component',k+1,': g=f composed with local reflection on all 64 input cells')


# Merge adjacent grid cells with the same active set. These exact incidence
# components give the chamber interiors modulo null grid lines.
for w in ('f','g'):
    todo=set(cells); chambers=[]
    while todo:
        start=min(todo); todo.remove(start); comp={start}; stack=[start]
        while stack:
            c=stack.pop()
            for d in cells:
                if d in todo and sum(abs(c[j]-d[j]) for j in range(2))==1 and active[w][c]==active[w][d]:
                    todo.remove(d);comp.add(d);stack.append(d)
        chambers.append(comp)
    print('CHAMBERS',w,[(tuple(sorted(active[w][min(c)])),len(c)) for c in chambers])
    witnesses={}
    for c in cells:
        for j in range(2):
            d=tuple(c[t]+(1 if j==t else 0) for t in range(2))
            if d not in active[w]:continue
            for left,right in [(c,d),(d,c)]:
                leaving=active[w][left]-active[w][right]
                if len(leaving)==1:
                    r=next(iter(leaving));witnesses.setdefault(r,(left,right))
    assert set(witnesses)==set(rules)
    # The three prescribed marker faces are explicit codimension-one witnesses.
    for face,inside,outside,branch in [
        ('x=4,y=1/2',(F(7,2),F(1,2)),(F(5),F(1,2)),'R'),
        ('x=1/2,y=4',(F(1,2),F(7,2)),(F(1,2),F(5)),'S'),
        ('x=4,y=7/2',(F(7,2),F(7,2)),(F(5),F(7,2)),'N')]:
        assert active[w][cell_of(inside)]=={branch}
        assert active[w][cell_of(outside)]=={'I'}
        assert sum(abs(cell_of(inside)[j]-cell_of(outside)[j]) for j in range(2))==1
        print('USB',w,face,branch,'uniquely leaves; I uniquely enters')

# Conditional algebraic step of MST: Gaussian-mixture identifiability must
# separately justify fixing both components from unequal positive weights.
a,b,c,d=sp.symbols('a b c d',real=True)
A=sp.Matrix([[a,b],[c,d]])
eqs=[]
for D in covs:eqs += list(A*D*A.T-D)
solutions=sp.solve(eqs,(a,b,c,d),dict=True)
assert len(solutions)==4
accepted=[]
for sol in solutions:
    B=A.subs(sol)
    if B*(means[0]-means[1])==means[0]-means[1]:accepted.append(B)
assert accepted==[I]
assert means[0]-I*means[0]==sp.zeros(2,1)
print('MST ALGEBRA: four common covariance stabilizers; only I fixes the mean difference; translation=0')
print('MST OBLIGATION: finite Gaussian-mixture identifiability, not assumed or established by this finite check')
assert sp.Matrix([F(1,2),F(1,2)]) != S*sp.Matrix([F(1,2),F(1,2)])
print('MID OBSTRUCTION: f != g on (0,1)^2; conditional on MST no common affine link')
print('LID DISTINCTION: identical latent laws admit alpha=Id; this is not a failure of LID')
print('FAMILY: all kernel equalities are coefficientwise per component, hence hold for every lambda in (0,1)')
print('FAMILY: component-fixing from weights uses lambda != 1/2; the same MST algebra then applies')
print('ALL_CHECKS_OK')
