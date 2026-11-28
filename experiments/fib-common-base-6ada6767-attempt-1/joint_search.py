# Publication repair (2026-10-23): resolve original suppliers within this checkout.
"""Joint four-endpoint discovery with p and v free. All accepted rows remain uncertified."""
from discover import *
import warnings,sys
warnings.filterwarnings('ignore')
LIMIT=int(sys.argv[1]) if len(sys.argv)>1 else 1000
rng=np.random.default_rng(6767)
def worlds(z):
    R,B,k,h=z[:4];S,C,l,t=z[4:]
    return [(R-h,B-k*h),(R+h,B+k*h),(S-t,C-l*t),(S+t,C+l*t)]
def equations(z):
    w=[source(*x) for x in worlds(z)]
    h=z[3];t=z[7]
    pj=(w[0][0]+w[1][0])/2;pm=(w[2][0]+w[3][0])/2
    # Saturation by p-24 is legitimate only for off-24 candidates.
    return np.r_[(w[1]-w[0])[[0,1,2]]/(2*h)/(pj-24),
                 (w[3]-w[2])[[0,1,3]]/(2*t)/(pm-24),
                 (w[0]-w[2])[:2]/(pj-24)]/np.array([1,2,5,1,2,5,1,2])
def valid(z):
    return abs(z[3])>1e-5 and abs(z[7])>1e-5 and all(0<r<1 and b>0 for r,b in worlds(z))
def mkpair(ww):
    (r,b),(s,c)=ww
    return [(r+s)/2,(b+c)/2,(c-b)/(s-r),(s-r)/2]
if __name__=='__main__':
    data=json.loads((HERE/'parent_curves.json').read_text())
    jseeds=[mkpair([q[:2] for q in d['endpoints']]) for d in data['J'][::20]]
    mseeds=[mkpair([q[:2] for q in d['endpoints']]) for d in data['M'][::20]]
    # Established positive-J arc centre, followed by its other endpoint search.
    import ast
    code=ast.parse((Path(__file__).resolve().parents[1]/'fib-self-calibration/quantitative-arc/quantitative_certificate.py').read_text())
    blocks=next(ast.literal_eval(q.value) for q in code.body if isinstance(q,ast.Assign) and any(isinstance(n,ast.Name) and n.id=='blocks' for n in q.targets))
    for lo,hi,tc,bc in blocks:
        for pp in np.linspace(lo,hi,4):
            tt=2*(pp-lo)/(hi-lo)-1;bb=np.polynomial.polynomial.polyval(tt,np.array(bc)*1e-20);rr=bb*np.polynomial.polynomial.polyval(tt,np.array(tc)*1e-20)
            w=source(rr,bb)
            for guess in [[1.125-rr,.125],[.85,.12],[.2,.13]]:
                f=lambda rb:(source(*rb)-w)[[0,1]]/np.array([100,200])
                sol=root(f,guess,jac=lambda rb:jac(f,rb),tol=1e-10)
                r1,b1=sol.x
                if 0<r1<1 and b1>0 and np.linalg.norm(sol.x-[rr,bb])>.01 and np.max(np.abs(f(sol.x)))<1e-8 and abs(source(r1,b1)[2]-w[2])<1e-3:
                    jseeds.append(mkpair([(rr,bb),(r1,b1)]));break
    out={'status':'numerical only','roots':[],'best':[],'starts':0,'seed_counts':[len(jseeds),len(mseeds)]};t0=time.monotonic();best=1e20
    for i in range(LIMIT):
        if i<len(jseeds)*len(mseeds):
            a=jseeds[i%len(jseeds)];b=mseeds[(i//len(jseeds))%len(mseeds)];z=np.r_[a,b]
            if i>=len(jseeds):z=z*(1+rng.normal(0,.08,8))
        else:
            def randpair():
                R=rng.uniform(.1,.9);h=rng.uniform(.001,min(R,1-R)*.95);B=10**rng.uniform(-4,1);k=B/h*rng.uniform(-.95,.95);return [R,B,k,h]
            z=np.r_[randpair(),randpair()]
        try:
            sol=root(equations,z,jac=lambda zz:jac(equations,zz),tol=1e-10,options={'maxfev':500})
            e=float(np.max(np.abs(equations(sol.x))));good=valid(sol.x)
            if good and e<best:
                best=e;row={'start':i,'z':sol.x.tolist(),'error_scaled':e,'worlds':[list(w) for w in worlds(sol.x)],'responses':[source(*w).tolist() for w in worlds(sol.x)]};out['best'].append(row);print('BEST',i,e,row,flush=True)
            if good and e<1e-7 and abs(source(*worlds(sol.x)[0])[0]-24)>1e-4:
                ww=np.array(worlds(sol.x));
                if not any(np.linalg.norm(ww-np.array(q['worlds']))<1e-3 for q in out['roots']):out['roots'].append(row);print('ROOT',row,flush=True)
        except (ValueError,FloatingPointError,ZeroDivisionError,OverflowError):pass
        out['starts']=i+1
        if (i+1)%100==0:
            (HERE/'joint_candidates.json').write_text(json.dumps(out,indent=2)+'\n');print('COUNT',i+1,'seconds',time.monotonic()-t0,'roots',len(out['roots']),flush=True)
        if len(out['roots'])>=8:break
    (HERE/'joint_candidates.json').write_text(json.dumps(out,indent=2)+'\n')
