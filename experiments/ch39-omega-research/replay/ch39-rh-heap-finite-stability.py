"""Rational finite classification tests for four-cycle fields with measurement error."""
from fractions import Fraction as Q
from itertools import permutations,product
from pathlib import Path
from runpy import run_path

zero=run_path(str(Path(__file__).with_name("ch39-rh-heap-zero-fields.py")))
mul,conj=zero["mul"],zero["conj"]
compose=zero["compose"]
states,klein,four,cosets,cid,even,odd=zero["quotient"]()
trails=[(p,compose(p,g),compose(p,compose(g,g))) for p,g in product(states,four)]
edges=list(product(even,odd))
basis=((1,0,0,0),(0,1,0,0),(0,0,1,0))


def dot(a,b):
    return sum(x*y for x,y in zip(a,b))


def classify(field):
    refs=[field[c[0]] for c in cosets]
    near={edge for edge in edges if 2*abs(dot(refs[edge[0]],refs[edge[1]]))>1}
    if len(near)==9:
        return (), tuple(dot(refs[0],v)>0 for v in refs[1:])
    assert len(near)==3 and all(sum(j in e for e in near)==1 for j in range(6))
    return tuple(sorted(near)), tuple(dot(refs[i],refs[j])>0 for i,j in sorted(near))


def representatives():
    for match in [None,*permutations(odd)]:
        labels=[0]*6
        if match is not None:
            for i,(v,w) in enumerate(zip(even,match)):
                labels[v]=labels[w]=i
        signs={}
        for i,c in enumerate(cosets):
            p=c[0]
            signs[p]=1
            for g in four:
                q=compose(p,g)
                signs[compose(q,g)]=1 if labels[i]==labels[cid[q]] else -1
        for flips in product((-1,1),repeat=6):
            if match is None and flips[0]<0 or match is not None and any(flips[i]<0 for i in even):
                continue
            yield {p:tuple(flips[cid[p]]*signs[p]*x for x in basis[labels[cid[p]]]) for p in states}


def perturbations():
    # Different rational tangent rotations at all sites, including distinct axes.
    fields={}
    for j,p in enumerate(states):
        v=[Q((j+2*k)%5-2,2**25) for k in (0,1,2)]
        size=dot(v,v)
        fields[p]=tuple([(1-size)/(1+size)]+[2*x/(1+size) for x in v])
    return fields


def main():
    eps=Q(1,1000)
    assert 12*eps<2 and 1-(50*eps)**2/2>10*eps
    assert 40*eps<1 and (Q(40,21)-1)**2<1-40*eps
    assert 2*(35*eps)**2<1 and 35*Q(1,2**20)<Q(1,2**14)
    # Four chord errors around the cycle imply |d^4-1| <= 4 eps. The
    # scalar threshold 1/2 on the six projective representatives has slack.
    assert 1-(10*eps)**2/2>Q(1,2)>10*eps
    changes=perturbations()
    keys=set()
    for field in representatives():
        label=classify(field)
        assert label not in keys
        keys.add(label)
        sample={p:mul(field[p],changes[p]) for p in states}
        assert all(dot(v,v)==1 for v in sample.values()) and classify(sample)==label
        for p,q,r in trails:
            pred=mul(mul(sample[q],conj(sample[p])),sample[q])
            difference=tuple(x-y for x,y in zip(pred,sample[r]))
            assert dot(difference,difference)<Q(1,2**40)
    assert len(keys)==80
    print("All 80 inequivalent classes decoded after rational site perturbations, 144 tests per field.")
    print("Worst-case theorem: local chord residual <= epsilon <= 1/1000 gives distance <= 35 epsilon.")


if __name__=="__main__":
    main()
