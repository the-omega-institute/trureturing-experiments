"""Free-base four-endpoint joint search in the reused t,b source chart."""
from discover import *
import warnings
warnings.filterwarnings('ignore')
rng=np.random.default_rng(6767004)
def resp(t,b): return np.asarray(response(t,b)[:4])
def rb(z):return [[z[2*i]*z[2*i+1],z[2*i+1]] for i in range(4)]
def equations(z):
    w=[resp(*z[2*i:2*i+2]) for i in range(4)];rs=[r for r,b in rb(z)];pj=(w[0][0]+w[1][0])/2;pm=(w[2][0]+w[3][0])/2
    return np.r_[(w[1]-w[0])[[0,1,2]]/(rs[1]-rs[0])/(pj-24),(w[3]-w[2])[[0,1,3]]/(rs[3]-rs[2])/(pm-24),(w[0]-w[2])[:2]/(pj-24)]/np.array([1,2,5,1,2,5,1,2])
def distinct(z):
    w=np.array(rb(z));return np.linalg.norm(w[0]-w[1])>1e-4 and np.linalg.norm(w[2]-w[3])>1e-4
def good(z):return all(0<r<1 and b>0 for r,b in rb(z)) and distinct(z)
def tb(ws):return np.array([[r/b,b] for r,b in ws]).ravel()
cat=json.loads((HERE/'pair_catalogue.json').read_text());cur=json.loads((HERE/'parent_curves.json').read_text());jseeds=[];mseeds=[]
for key,term,coll in [('J',2,jseeds),('M',3,mseeds)]:
    coll.extend(tb([q[:2] for q in row['endpoints']]) for row in cur[key][::20])
    for row in cat[key]:
        if max(abs(row['responses'][0][i]-row['responses'][1][i]) for i in [0,1,term])<1e-5 and min(q[1] for q in row['endpoints'])>1e-5:coll.append(tb(row['endpoints']))
out={'status':'numerical only','starts':0,'candidates':[],'inadmissible_intersections':[]}
for i in range(4000):
    if i<2000:
        js=jseeds[-1] if i%2 else jseeds[rng.integers(len(jseeds))];ms=mseeds[rng.integers(len(mseeds))];z=np.r_[js,ms];z=z*(1+rng.normal(0,.25,8))
    else:z=np.ravel([[rng.uniform(.001,.999)/(b:=10**rng.uniform(-5,2)),b] for k in range(4)])
    try:
        sol=root(equations,z,jac=lambda zz:jac(equations,zz),tol=1e-10,options={'maxfev':500});e=float(np.max(np.abs(equations(sol.x))))
        if e<1e-8 and distinct(sol.x):
            w=[resp(*sol.x[2*k:2*k+2]) for k in range(4)];absres=max(abs(w[0][0]-w[k][0])+abs(w[0][1]-w[k][1]) for k in range(1,4));absres=max(absres,abs(w[0][2]-w[1][2]),abs(w[2][3]-w[3][3]))
            if absres<1e-4:
                key='candidates' if good(sol.x) else 'inadmissible_intersections'
                if not any(np.linalg.norm(sol.x-np.array(q['z']))<.05 for q in out[key]):
                    row={'start':i,'z':sol.x.tolist(),'worlds':rb(sol.x),'error_scaled':e,'absolute_residual_max':float(absres),'responses':[q.tolist() for q in w]};out[key].append(row);print(key,row,flush=True)
    except (ValueError,ZeroDivisionError,OverflowError):pass
    out['starts']=i+1
    if (i+1)%500==0:print('COUNT',i+1,flush=True);(HERE/'joint_tb_candidates.json').write_text(json.dumps(out,indent=2)+'\n')
    if len(out['candidates'])>=8:break
(HERE/'joint_tb_candidates.json').write_text(json.dumps(out,indent=2)+'\n')
