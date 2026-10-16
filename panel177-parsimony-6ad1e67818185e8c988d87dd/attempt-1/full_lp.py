import numpy as np, scipy, json, itertools, pathlib
from scipy.optimize import linprog
from scipy.sparse import lil_matrix, csr_matrix
ROOT=pathlib.Path(__file__).parent

def model(L,kind):
    words=np.array(list(itertools.product([0,1], repeat=L)))
    n=np.array([next((i+1 for i,t in enumerate(w) if t==0),0) for w in words])
    schedules=[()]+[(i,) for i in range(L)]+list(itertools.combinations(range(L),2))
    groups=[]; nv=0
    for s in schedules:
        nc=1 if not s else 2 if len(s)==1 else 3 if kind=='O' else 4
        groups.append((s,nv,nc));nv+=nc*(L+1)
    eq=lil_matrix((sum(nc-1 for _,_,nc in groups)+1,nv+1));rhs=[];ii=0
    for s,start,nc in groups:
        eq[-1,start:start+L+1]=1
        for c in range(1,nc):
            eq[ii,start:start+L+1]=-1;eq[ii,start+c*(L+1):start+(c+1)*(L+1)]=1;rhs.append(0);ii+=1
    rhs.append(1)
    nr=3*(2**L-1)+1;A=lil_matrix((nr,nv+1));base=np.zeros(nr);slope=np.zeros((nr,L));rows=[];rr=0
    for wi,w in enumerate(words):
        events=[(0,),(n[wi],),(0,n[wi])] if n[wi] else [(0,)]
        for typ,E in enumerate(events):
            rows.append((wi,typ))
            base[rr]=1 if typ in (0,2) else 0
            if not n[wi]:base[rr]=1
            else:slope[rr,n[wi]-1]=-1 if typ==0 else 1 if typ==1 else 0
            A[rr,-1]=-1
            for s,start,nc in groups:
                cell=0 if not s else w[s[0]] if len(s)==1 else (0 if w[s[0]]==0 else 1+w[s[1]]) if kind=='O' else 2*w[s[0]]+w[s[1]]
                for z in E:A[rr,start+cell*(L+1)+z]=-1
            rr+=1
    c=np.zeros(nv+1);c[-1]=1
    return dict(L=L,kind=kind,words=words,n=n,groups=groups,eq=eq.tocsr(),rhs=np.array(rhs),A=A.tocsr(),base=base,slope=slope,c=c,rows=rows)

def solve(m,b):
    return linprog(m['c'],A_ub=m['A'],b_ub=-m['base']-m['slope']@b,A_eq=m['eq'],b_eq=m['rhs'],bounds=(0,None),method='highs')

if __name__=='__main__':
    rng=np.random.default_rng(177);L=8;ms={k:model(L,k) for k in ['N','O']};out=[]
    for i in range(24):
        v=np.sort(rng.integers(16,65,size=L))[::-1];v[2]=v[1]
        rs={k:solve(m,v/64) for k,m in ms.items()}
        row={'i':i,'prefix_numerators':v.tolist(),'N':rs['N'].fun,'O':rs['O'].fun,'gap':rs['O'].fun-rs['N'].fun,'status':{k:r.status for k,r in rs.items()}}
        out.append(row);print(json.dumps(row),flush=True)
        for k,r in rs.items():
            np.savez_compressed(ROOT/f'lp-{i:02d}-{k}.npz',x=r.x,dual=r.ineqlin.marginals,equality_dual=r.eqlin.marginals)
    (ROOT/'lp-summary.json').write_text(json.dumps({'numpy':np.__version__,'scipy':scipy.__version__,'results':out},indent=2)+'\n')
