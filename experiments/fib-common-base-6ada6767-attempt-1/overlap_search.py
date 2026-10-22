"""Three-world overlap discovery in the already published t=r/b regular chart."""
from discover import *
import warnings
warnings.filterwarnings('ignore')
rng=np.random.default_rng(6767003)
def resp(t,b):return np.asarray(response(t,b)[:4])
def equations(z):
    w=[resp(z[2*i],z[2*i+1]) for i in range(3)];p=(w[0][0]+w[1][0]+w[2][0])/3
    return np.r_[(w[0]-w[1])[[0,1,2]],(w[0]-w[2])[[0,1,3]]]/(p-24)/np.array([1,2,5,1,2,5])
def sourceworlds(z):return [[z[2*i]*z[2*i+1],z[2*i+1]] for i in range(3)]
def distinct(z):
    w=sourceworlds(z);return np.linalg.norm(np.array(w[0])-w[1])>1e-4 and np.linalg.norm(np.array(w[0])-w[2])>1e-4
def valid(z):return all(0<r<1 and b>0 for r,b in sourceworlds(z)) and distinct(z)
seed=np.array([2.0783,.00055,.6847/.32344,.32344,.9723/.02448,.02448])
out={'status':'numerical only','candidates':[],'real_roots_outside_domain':[],'starts':0}
for i in range(1000):
    if i<500:z=seed*(1+rng.normal(0,.1,6))
    else:z=np.ravel([[rng.uniform(.001,.999)/(b:=10**rng.uniform(-4,1)),b] for k in range(3)])
    try:
        sol=root(equations,z,jac=lambda zz:jac(equations,zz),tol=1e-11,options={'maxfev':500});e=float(np.max(np.abs(equations(sol.x))))
        if e<1e-8 and distinct(sol.x):
            row={'z':sol.x.tolist(),'error_scaled':e,'worlds':sourceworlds(sol.x),'responses':[resp(*sol.x[2*k:2*k+2]).tolist() for k in range(3)]}
            key='candidates' if valid(sol.x) else 'real_roots_outside_domain'
            if not any(np.linalg.norm(np.array(q['z'])-sol.x)<.01 for q in out[key]):out[key].append(row);print(key,row,flush=True)
    except (ValueError,ZeroDivisionError,OverflowError):pass
    out['starts']=i+1
(HERE/'overlap_candidates.json').write_text(json.dumps(out,indent=2)+'\n')
print('counts',len(out['candidates']),len(out['real_roots_outside_domain']),flush=True)
