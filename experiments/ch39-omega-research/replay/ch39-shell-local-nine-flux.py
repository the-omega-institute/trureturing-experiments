"""Minimum degree of an indexed chord law when every unmarked edge is also prescribed."""
from itertools import combinations,product
from math import comb,isqrt
from fractions import Fraction as F

def classes(r,cap):
    sites=[p for p in product(range(isqrt(r)+1),repeat=3) if sum(x*x for x in p)==r]
    ds={sum((x-y)**2 for x,y in zip(v,w)) for v,w in combinations(sites,2)}
    return sorted(d//2 for d in ds if d<=cap)
sets=(classes(17,14),classes(867,804))
def solve9(rows,values):
    mat=[[x%9 for x in row]+[y%9] for row,y in zip(rows,values)]
    m,n=len(mat),len(rows[0])
    order=list(range(n))
    used=0
    def pivot(mod,k):
        for i in range(k,m):
            for j in range(k,n):
                if mat[i][j]%3:
                    mat[i],mat[k]=mat[k],mat[i]
                    order[j],order[k]=order[k],order[j]
                    for row in mat:
                        row[j],row[k]=row[k],row[j]
                    inv=pow(mat[k][k],-1,mod)
                    mat[k]=[x*inv%mod for x in mat[k]]
                    for a in range(k+1,m):
                        b=mat[a][k]
                        mat[a]=[(x-b*y)%mod for x,y in zip(mat[a],mat[k])]
                    return True
        return False
    while pivot(9,used):
        used+=1
    free=used
    for row in mat[used:]:
        if row[-1]%3:
            return None
        row[:]=[z//3 for z in row]
    while pivot(3,free):
        free+=1
    if any(row[-1]%3 for row in mat[free:]):
        return None
    vals=[0]*n
    for i in reversed(range(free)):
        vals[i]=(mat[i][-1]-sum(a*b for a,b in zip(mat[i][i+1:n],vals[i+1:])))%(9 if i<used else 3)
    sol=[0]*n
    for i,j in enumerate(order):
        sol[j]=vals[i]
    assert all((sum(a*b for a,b in zip(row,sol))-t)%9==0 for row,t in zip(rows,values))
    return sol
pairs=[(d,ref,i,j) for i,(ds,ref,tag) in enumerate(zip(sets,(3,122),((6,7),(362,402))))
       for d in ds if d!=ref for j in [tag.index(d) if d in tag else -1]]
def build(n):
    return [[comb(d,k)-comb(ref,k) for k in range(1,n+1)] for d,ref,_,_ in pairs]
def solve(n):
    rows=build(n)
    out=[]
    for a,b in ((1,1),(1,-1)):
        signs=(a,b)
        vals=[signs[i]*((2,1),(4,6))[i][j] if j>=0 else 0 for _,_,i,j in pairs]
        coeff=solve9(rows,vals)
        if coeff is not None:
            out.append((signs,[(k+1,z if z<=4 else z-9) for k,z in enumerate(coeff) if z]))
    return out
lo,hi=0,402
while lo+1<hi:
    mid=(lo+hi)//2
    if solve(mid):
        hi=mid
    else:
        lo=mid
assert hi==54 and not solve(hi-1)
print("class counts",*[len(s) for s in sets],"minimum local binomial degree",hi)
print("first witnesses by shell orientation",solve(hi))
