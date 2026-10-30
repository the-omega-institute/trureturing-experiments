import itertools,json,time

def literal(a,b):
 if any(x>y for x,y in zip(a,b)):return False
 n=len(a)
 return all(a[i]-a[j]<=b[i]-b[j] for i in range(n) for j in range(i+1,n) if all(a[k]>a[j] for k in range(i,j)))
def compute(bounds):
 t=time.monotonic(); vs=list(itertools.product(*(range(k+1) for k in bounds)));N=len(vs)
 le=[[literal(a,b) for b in vs] for a in vs]
 h=[[-1]*N for _ in vs]
 for p in range(N):
  h[p][p]=0
  for q in range(p+1,N):
   if le[p][q]:h[p][q]=1+max(h[p][r] for r in range(p,q) if le[p][r] and le[r][q])
 fs=[]
 for b in range(N):
  c=[0]*(sum(vs[b])+1)
  for p in range(N):
   for q in range(N):
    if le[p][q] and le[q][b]:c[h[p][q]]+=1
  fs.append(c)
 import subprocess
 d=json.loads(subprocess.check_output(['./probe','control',*map(str,bounds)],text=True))
 assert d['dp_height_controls']==sum(map(sum,le))
 assert d['last_coefficients']==fs[-1]
 for b,c in zip(vs,fs):
  z=json.loads(subprocess.check_output(['./probe','single',*map(str,b)],text=True))
  assert z['last_coefficients']==c,(b,c,z)
 return {'bounds':bounds,'vertices':N,'literal_comparabilities':sum(map(sum,le)),'complete_F_comparisons':N,'elapsed_seconds':time.monotonic()-t}
res=[]
for b in [[3,3,3],[2,2,2,2],[1,1,1,1,1,1],[2,1],[1,2,0,2]]:
 r=compute(b);res.append(r);print(json.dumps(r),flush=True)
 with open('independent-controls.json','w') as f:json.dump(res,f)
