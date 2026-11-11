"""Bounded complete-catalog exploration, all actual finite rows, common mass.
Predeclared: strict class-gap >1e-8 is only a candidate and needs exact
certification and original transport. Absence of a gap proves no containment.
64 deterministic random/structured profiles at each L=8,9,10; seed177031.
No fixed-law replacement or normalized affine mass substitution is used.
"""
from pathlib import Path
import json, itertools, sys, platform
import numpy as np
import scipy
from scipy.sparse import coo_matrix, hstack
from scipy.optimize import linprog
ROOT=Path(__file__).resolve().parent

def model(L,C):
    words=list(itertools.product((0,1),repeat=L)); ns=[next((i+1 for i,b in enumerate(w) if b==0),0) for w in words]
    schedules=[()]+[(i,) for i in range(L)]+list(itertools.combinations(range(L),2))
    eqr=[];eqc=[];eqv=[];er=[];ec=[];ev=[];catalog=[];col=0;eqrow=1
    for S in schedules:
        nc=1 if not S else 2 if len(S)==1 else 3 if C=='O' else 4
        t=col;col+=1;eqr.append(0);eqc.append(t);eqv.append(1)
        cellcols=[]
        for cell in range(nc):
            inds=list(range(col,col+L+1));col+=L+1;cellcols.append(inds)
            eqr.extend([eqrow]*(L+2));eqc.extend(inds+[t]);eqv.extend([1]*(L+1)+[-1]);eqrow+=1
        for wi,w in enumerate(words):
            cell=0 if not S else w[S[0]] if len(S)==1 else (0 if w[S[0]]==0 else 1+w[S[1]]) if C=='O' else 2*w[S[0]]+w[S[1]]
            inds=cellcols[cell]
            if ns[wi]:
                row=3*wi; er.extend([row,row+1,row+2,row+2]);ec.extend([inds[0],inds[ns[wi]],inds[0],inds[ns[wi]]]);ev.extend([-1.]*4)
            else:
                er.append(3*wi);ec.append(inds[0]);ev.append(-1.)
        catalog.append({'schedule':[i+1 for i in S],'mass_col':t,'cell_cols':cellcols})
    neq=eqrow;nrows=3*(len(words)-1)+1
    A=coo_matrix((ev,(er,ec)),shape=(nrows,col)).tocsr();E=coo_matrix((eqv,(eqr,eqc)),shape=(neq,col)).tocsr()
    objective=np.zeros(col+1);objective[-1]=1.
    A=hstack([A,-np.ones((nrows,1))],format='csr');E=hstack([E,np.zeros((neq,1))],format='csr')
    eqrhs=np.zeros(neq);eqrhs[0]=1.
    return {'A':A,'E':E,'eqrhs':eqrhs,'objective':objective,'ns':ns,'catalog':catalog,'L':L,'C':C}

def solve(m,a):
    d=[]
    for n in m['ns']:
        d.extend([1-a[n-1],a[n-1],1.] if n else [1.])
    res=linprog(m['objective'],A_ub=m['A'],b_ub=-np.array(d),A_eq=m['E'],b_eq=m['eqrhs'],bounds=(0,None),method='highs')
    if not res.success:raise RuntimeError(res.message)
    return res

def main():
    rng=np.random.default_rng(177031);results=[];best=None
    for L in (8,9,10):
        models={C:model(L,C) for C in ('N','O')}
        for k in range(64):
            if k<16:
                a=np.sort(rng.random(L))[::-1]
            elif k<32:
                a=np.sort(rng.beta(2.,2.,L))[::-1];a[1]=a[2]
            elif k<48:
                z=rng.uniform(.15,.85);a=np.clip(z+np.sort(rng.normal(0,.08,L))[::-1],.01,.99)
            else:
                a=np.sort(rng.integers(1,65,L)/64.)[::-1]
            rn=solve(models['N'],a);ro=solve(models['O'],a);gap=ro.fun-rn.fun
            row={'L':L,'index':k,'a':a.tolist(),'N_error':float(rn.fun),'O_error':float(ro.fun),'gap':float(gap)};results.append(row)
            if best is None or gap>best['gap']:best=row
            if gap>1e-8:
                np.savez_compressed(ROOT/('candidate-L%d-k%d.npz'%(L,k)),a=a,N=rn.x,O=ro.x,O_dual=ro.ineqlin.marginals)
                (ROOT/('catalog-L%d.json'%L)).write_text(json.dumps({C:models[C]['catalog'] for C in models}))
            if k%16==15:print(json.dumps({'completed_L':L,'through':k,'maximum_gap':best['gap']}),flush=True)
    output={'predeclared_count':192,'completed_count':len(results),'seed':177031,'models':'complete nonnegative common-cell mass, all labels/all words/A-C-P-T rows, fixed profile minimize uniform TV','python':sys.version,'scipy':scipy.__version__,'numpy':np.__version__,'platform':platform.platform(),'best':best,'results':results,'claim':'floating-point exploration only; no universal theorem or exact original separator'}
    (ROOT/'search_complete_results.json').write_text(json.dumps(output,indent=2)+'\n');print(json.dumps({'best':best,'count':len(results)}),flush=True)
if __name__=='__main__':main()
