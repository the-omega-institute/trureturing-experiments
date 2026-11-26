"""Exhaust simple unitary embeddings with canonical matching and an adjoint loop parameter below one."""
from functools import lru_cache
from itertools import product
from collections import Counter
from math import prod
from flint import fmpz_mat,arb

# Modes are SU(3)_X, SU(4), SU(2)_L, SU(2)_R, SU(2)_S. Family is omitted,
# so any true product-group embedding restricts to a case included here.
ceiling=int((4*arb.pi()*arb("40.9906")).floor().unique_fmpz())
assert ceiling==515
levelsets=[(x,c,c,r,s) for c in range(1,ceiling//11+1) for x in range(1,(3*c)//16+1)
           for r in range(1,c+1) for s in range(1,9*c//4+1) if 9*c==48*x+9*r+4*s]
max_dim=max(ceiling//z[1] for z in levelsets)
@lru_cache(None)
def weights(top):
    # GT restriction of U(n); weights are the counts of each letter in the tableau.
    if len(top)==1:
        return [(top[0],)]
    return [lower+(sum(top)-sum(row),)
            for row in product(*(range(lo,hi+1) for hi,lo in zip(top,top[1:])))
            for lower in weights(row)]
def simple(n):
    reps=[]
    def traverse(prefix):
        if len(prefix)<n-1:
            i=0
            while True:
                t=prefix+(i,)+(0,)*(n-2-len(prefix))
                d=prod(sum(t[k:l])+l-k for k in range(n) for l in range(k+1,n))//prod(
                    l-k for k in range(n) for l in range(k+1,n))
                if d>max_dim:
                    break
                traverse(prefix+(i,))
                i+=1
        else:
            top=tuple(sum(prefix[j:]) for j in range(n))
            w=weights(top)
            index=sum((v[0]-v[1])**2 for v in w)//2
            reps.append((len(w),index,prefix))
    traverse(())
    return reps
reps=(simple(3),simple(4),simple(2),simple(2),simple(2))
samples=((1,2,0,0),(1,2,1,0),(1,2,0,1),(1,2,1,1),(1,0,2,1),(0,0,1,1),(0,0,0,1))
def light(t):
    a,b,w,h=t
    col=(a,b,-a-b)
    return 3*sum(v**5 for v in ([x+k*w+h for x in col for k in (-1,1)]
                               +[-x+k*h for x in col for k in (-4,2)]+[w-3*h,-w-3*h,6*h,0]))
def canonical(rep):
    return min(rep,rep[::-1])
def blocks(target,nmax):
    found=[]
    def visit(current,d):
        i=len(current)
        if i==5:
            cost=tuple(t[1]*(d//t[0]) for t in current)
            if any(cost):
                found.append((cost,d,tuple(t[2] for t in current)))
            return
        for t in reps[i]:
            dd=d*t[0]
            if dd<=nmax and t[2]==canonical(t[2]) and t[1]*d<=target[i] and all(
                    z[1]*(dd//z[0])<=target[k] for k,z in enumerate(current)):
                visit(current+(t,),dd)
    visit((),1)
    return found
def suffix_solver(blocks,target,nmax):
    # Generate multisets in pivot order by exhausting one factor at a time.
    @lru_cache(None)
    def walk(rem,room,start):
        if not any(rem):
            return 1
        if room<=0:
            return 0
        j=next(i for i,r in enumerate(rem) if r)
        total=0
        for k in range(start,len(blocks)):
            c,d,_=blocks[k]
            if not c[j] or d>room or any(a>b for a,b in zip(c,rem)):
                continue
            nr=tuple(b-a for a,b in zip(c,rem))
            total+=walk(nr,room-d,k if nr[j] else 0)
        return total
    def gen(rem,room,start,prefix):
        if not any(rem):
            yield prefix
            return
        j=next(i for i,r in enumerate(rem) if r)
        for k in range(start,len(blocks)):
            c,d,_=blocks[k]
            if not c[j] or d>room or any(a>b for a,b in zip(c,rem)):
                continue
            nr=tuple(b-a for a,b in zip(c,rem))
            ns=k if nr[j] else 0
            if walk(nr,room-d,ns):
                yield from gen(nr,room-d,ns,prefix+(k,))
    return walk,gen
@lru_cache(None)
def moment(block):
    top=[tuple(sum(z[k:]) for k in range(len(z)+1)) for z in block]
    out=[]
    for a,b,w,h in samples:
        coords=((4*h,4*h,-8*h),(a+h,b+h,-a-b+h,-3*h),(w,-w),(3*h,-3*h),(2*h,-2*h))
        ev=[[sum(i*j for i,j in zip(weight,row)) for weight in weights(t)] for t,row in zip(top,coords)]
        e=[sum(z) for z in product(*ev)]
        out.append((sum(x**3 for x in e),sum(x**5 for x in e)))
    return out
def orientations(block):
    return [tuple(v) for v in product(*[sorted(set((p,p[::-1]))) for p in block])]
def examine():
    checked=Counter()
    viable=Counter()
    for target in levelsets:
        nmax=ceiling//target[1]
        choices=blocks(target,nmax)
        walk,gen=suffix_solver(choices,target,nmax)
        if walk(target,nmax,0):
            p2=[target[1]*(a*a+b*b+(a+b)**2+2*w*w+30*h*h) for a,b,w,h in samples]
            for ids in gen(target,nmax,0,()):
                choices_by_block=[orientations(choices[k][2]) for k in ids]
                for combination in product(*choices_by_block):
                    rows=[moment(blk) for blk in combination]
                    M=fmpz_mat([[sum(row[i][1] for row in rows),
                                 p2[i]*sum(row[i][0] for row in rows),light(t)]
                                for i,t in enumerate(samples)])
                    checked[target]+=1
                    if M.rank()<3:
                        viable[target]+=1
        walk.cache_clear()
    print("levels with unitary spectra inside adjoint budget",sorted(checked.items()))
    print("passing fifth moment",sorted(viable.items()),"tested with orientations",sum(checked.values()))
    return checked,viable
if __name__=="__main__":
    checked,viable=examine()
    assert sum(checked.values())==22964 and not viable
