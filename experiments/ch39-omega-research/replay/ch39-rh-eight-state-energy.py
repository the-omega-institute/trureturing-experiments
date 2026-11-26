"""Low positive energies in the finite eight-quaternion four-cycle action."""
from collections import Counter
from functools import lru_cache
from itertools import product
from math import prod
from pathlib import Path
from runpy import run_path

zero=run_path(str(Path(__file__).with_name("ch39-rh-heap-zero-fields.py")))
states,klein,four,cosets,cid,even,odd=zero["quotient"]()
compose=zero["compose"]
CUTOFF=20


def squares():
    index={p:i for i,p in enumerate(states)}
    edges={}
    for p,g in product(states,four):
        q=compose(p,g)
        r=compose(q,g)
        pair=tuple(sorted((index[p],index[r])))
        edges.setdefault(pair,set()).add(index[q])
    groups=[[] for _ in cosets]
    for (a,b),middle in edges.items():
        assert len(middle)==2
        i=cid[states[a]]
        groups[i].append((a,b,*sorted(middle)))
    assert len(edges)==36 and all(len(x)==6 for x in groups)
    return groups


@lru_cache(None)
def sign_polynomial(constraints):
    histogram=Counter()
    for signs in product((0,1),repeat=4):
        errors=sum((signs[a]^signs[b])!=v for a,b,v in constraints)
        if 8*errors<=CUTOFF:
            histogram[8*errors]+=1
    return histogram


def low_spectrum(labels,groups):
    additive=0
    polys=[]
    for coset,group in zip(cosets,groups):
        conditions=[]
        for a,b,c,d in group:
            if labels[a]!=labels[b] or (labels[a]==labels[c])!=(labels[a]==labels[d]):
                additive+=4
                if additive>CUTOFF:
                    return {}
            else:
                conditions.append((coset.index(states[a]),coset.index(states[b]),labels[a]!=labels[c]))
        polys.append(sign_polynomial(tuple(conditions)))
    counts={additive:1}
    for poly in polys:
        result=Counter()
        for a,n in counts.items():
            for b,m in poly.items():
                if a+b<=CUTOFF:
                    result[a+b]+=m*n
        counts=result
    return counts


def canonical(seq):
    labels={}
    return tuple(labels.setdefault(v,len(labels)) for v in seq)


def main():
    groups=squares()
    spectra=Counter()
    patterns={canonical(v) for v in product(range(4),repeat=4)}
    for pattern in sorted(patterns):
        mismatches=sum(pattern[i]!=pattern[j] for i in range(4) for j in range(i))
        if 4*mismatches>CUTOFF:
            continue
        # Left translations permute the six cosets transitively. Two impure
        # cosets already cost >=24. With no impurity, use just one location.
        arity=len(set(pattern))
        orbit=(6 if mismatches else 1)*prod(range(4-arity+1,5))
        for rest in product(range(4),repeat=5):
            labels=[pattern[cosets[0].index(p)] if cid[p]==0 else rest[cid[p]-1] for p in states]
            for e,n in low_spectrum(labels,groups).items():
                spectra[e]+=orbit*n
    assert spectra[0]==256+9216 and set(spectra)=={0,16,20}
    assert spectra[16]==442368
    print("Exact low-energy multiplicities among 8^24 fields:",sorted(spectra.items()))
    # A rational finite Boltzmann rule, with one binary factor per energy unit.
    # The complement starts at energy 24. Compare with all states as an upper bound.
    from fractions import Fraction as Q
    for bits in (2,4,8):
        r=Q(1,2**bits)
        known=sum(n*r**e for e,n in spectra.items())
        tail=(8**24-sum(spectra.values()))*r**24
        print("r=2^(-%d), lower probability of a three-line zero field:"%bits,float(9216/(known+tail)))
    assert 9216/(known+tail)>Q(9729,10000)


if __name__=="__main__":
    main()
