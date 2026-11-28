"""The final finite discovery batch tests the allowed two-world overlap stratum."""
from discover import *
import warnings
warnings.filterwarnings('ignore')
rng=np.random.default_rng(6767005)
def eq(z):
    a=np.asarray(response(z[0],z[1])[:4]);b=np.asarray(response(z[2],z[3])[:4]);sep=z[2]*z[3]-z[0]*z[1]
    return (b-a)/sep/((a[0]+b[0])/2-24)/np.array([1,2,5,5])
out={'status':'numerical only','starts':0,'candidates':[]}
for i in range(1060):
    z=np.ravel([[rng.uniform(.001,.999)/(b:=10**rng.uniform(-5,2)),b] for j in range(2)])
    try:
        sol=root(eq,z,jac=lambda z:jac(eq,z),tol=1e-10,options={'maxfev':500});e=float(np.max(np.abs(eq(sol.x))));z=sol.x
        rr=[z[0]*z[1],z[2]*z[3]]
        if 0<min(rr)<max(rr)<1 and min(z[1],z[3])>0 and abs(rr[0]-rr[1])>1e-5 and e<1e-8:
            ws=[np.asarray(response(z[2*j],z[2*j+1])[:4]) for j in range(2)]
            err=float(np.max(np.abs(ws[0]-ws[1])))
            if err<1e-5:out['candidates'].append({'z':z.tolist(),'absolute_residual':err,'response':[w.tolist() for w in ws]});print(out['candidates'][-1],flush=True)
    except (ValueError,ZeroDivisionError,OverflowError):pass
    out['starts']=i+1
(HERE/'same_pair_candidates.json').write_text(json.dumps(out,indent=2)+'\n')
print('FINAL_DISCOVERY_BATCH',out['starts'],'candidate_count',len(out['candidates']),flush=True)
