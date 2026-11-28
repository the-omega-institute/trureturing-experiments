"""Bounded free-base discovery of collision branches at finitely fixed separations."""
from discover import *
import warnings
warnings.filterwarnings('ignore')
rng=np.random.default_rng(10236767)
out={'status':'numerical only','J':[],'M':[],'starts':0}
hvals=[.01,.05,.1,.2,.3,.4]
for term,key in [(2,'J'),(3,'M')]:
    for h in hvals:
        for i in range(120):
            R=rng.uniform(h+.00001,1-h-.00001);B=10**rng.uniform(-4,1);k=B/h*rng.uniform(-.95,.95);z=[R,B,k]
            def f(z):
                R,B,k=z;w0=source(R-h,B-k*h);w1=source(R+h,B+k*h);p0=(w0[0]+w1[0])/2
                return (w1-w0)[[0,1,term]]/(2*h)/(p0-24)/np.array([1,2,5])
            try:
                sol=root(f,z,jac=lambda z:jac(f,z),tol=1e-10,options={'maxfev':500});err=float(np.max(np.abs(f(sol.x))));R,B,k=sol.x
                if valid(sol.x,h) and err<1e-7:
                    wm=source(R-h,B-k*h);wp=source(R+h,B+k*h)
                    row={'h':h,'z':sol.x.tolist(),'base':wm[:2].tolist(),'endpoints':[[float(R-h),float(B-k*h)],[float(R+h),float(B+k*h)]],'responses':[wm.tolist(),wp.tolist()],'error_scaled':err}
                    if abs(wm[0]-24)>1e-4 and not any(abs(h-q['h'])<1e-8 and np.linalg.norm(sol.x-q['z'])<1e-3 for q in out[key]):
                        out[key].append(row);print(key,row,flush=True)
            except (ValueError,ZeroDivisionError,OverflowError):pass
            out['starts']+=1
    (HERE/'pair_catalogue.json').write_text(json.dumps(out,indent=2)+'\n')
print('TOTAL',out['starts'],len(out['J']),len(out['M']))
