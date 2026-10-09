# Literal check: S312(r,k,d) = #{ r x k top-left corners of 312-avoiding square ASMs that are ASRs with d nonempty rows }.
import itertools, sys
sys.setrecursionlimit(10000)
from check import S, Sch
def rows(m, colsum):
    # all rows: entries in {-1,0,1}, nonzero alternate starting and ending with 1 (row sum 1), new column sums in {0,1}
    out=[]
    def rec(j, cur, last, acc):
        if j==m:
            if last==1: out.append(tuple(acc))
            return
        for v in (0,1,-1):
            if v==1 and last==1: continue
            if v==-1 and last!=1: continue
            ns=colsum[j]+v
            if ns<0 or ns>1: continue
            acc.append(v); rec(j+1, cur, v if v!=0 else last, acc); acc.pop()
    rec(0,None,-1,[])
    return out
def asms(m):
    res=[]
    def rec(i, colsum, M):
        if i==m:
            if all(c==1 for c in colsum): res.append(tuple(M))
            return
        for row in rows(m, colsum):
            ns=[colsum[j]+row[j] for j in range(m)]
            rec(i+1, ns, M+[row])
    rec(0,[0]*m,[])
    return res
def avoids312(M):
    m=len(M); ones=[(i,j) for i in range(m) for j in range(m) if M[i][j]==1]
    for a in ones:
        for b in ones:
            if b[0]<=a[0] or b[1]>=a[1]: continue
            for c in ones:
                if c[0]>b[0] and b[1]<c[1]<a[1]: return False
    return True
def is_asr(R):
    r=len(R); k=len(R[0]) if r else 0
    for row in R:
        nz=[v for v in row if v]
        if nz and (nz[0]!=1 or nz[-1]!=1 or any(nz[t]==nz[t+1] for t in range(len(nz)-1))): return False
    for j in range(k):
        nz=[R[i][j] for i in range(r) if R[i][j]]
        if nz and (nz[0]!=1 or any(nz[t]==nz[t+1] for t in range(len(nz)-1))): return False
    return True
M=int(sys.argv[1]) if len(sys.argv)>1 else 6
corners={}
for m in range(1,M+1):
    A=[a for a in asms(m) if avoids312(a)]
    print("m",m,"312-avoiding ASMs",len(A),"Sch(m-1)",Sch(m-1)); sys.stdout.flush()
    for a in A:
        for r in range(0,m+1):
            for k in range(0,m+1):
                C=tuple(tuple(a[i][:k]) for i in range(r))
                corners.setdefault((r,k),set()).add(C)
ok=0;bad=[]
for (r,k),cs in sorted(corners.items()):
    if r+k>M: continue
    for d in range(0,min(r,k)+1):
        cnt=sum(1 for C in cs if is_asr(C) and sum(1 for row in C if any(row))==d)
        if cnt!=S(r,k,d): bad.append((r,k,d,cnt,S(r,k,d)))
        else: ok+=1
print("checked (r,k,d) with r+k <=",M,": ok",ok,"bad",bad[:10])
