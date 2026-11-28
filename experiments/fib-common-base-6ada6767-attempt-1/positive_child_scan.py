"""Continue the positive J pair and test independently reconstructed companions.
Only candidate data; numerical roots do not establish an entire-child inventory.
"""
from discover import *
from full_parent import world,secant,pair_worlds,valid
from companions import inventory
cat=json.loads((HERE/'pair_catalogue.json').read_text())
seed=next(q for q in cat['J'] if q['base'][0]>0 and max(abs(q['responses'][0][i]-q['responses'][1][i]) for i in [0,1,2])<1e-5)
R,B,k=seed['z'];h=seed['h']
def drift(r,b):
 j=b+r*(1-r)/b;n=j+r;T=1+b;U=1+2*b;return 24-6*(8*b-1)*(5*b+2)/(5*T*U)/(1+6/T+1/U-4/j-4/n)
ds=[drift(*w) for w in seed['endpoints']];z0=np.r_[R,B,k,(ds[0]+ds[1])/2,(ds[1]-ds[0])/(2*h)]
out=[];steps=0
for hs in [np.linspace(.3,.0005,180),np.linspace(.3,.345,180)]:
 z=z0.copy()
 for idx,h in enumerate(hs):
    f=lambda zz:secant(np.r_[zz,h],2)/np.array([30,30,100,200,500])
    sol=root(f,z,jac=lambda zz:jac(f,zz),tol=1e-11,options={'maxfev':500});err=float(np.max(np.abs(f(sol.x))));steps+=1
    if not valid(np.r_[sol.x,h]) or err>1e-7:print('STOP',h,err,sol.x,flush=True);break
    z=sol.x;ww=pair_worlds(np.r_[z,h]);p,v=world(*ww[0])[:2]
    if idx%4==0 or idx==len(hs)-1:
      inv=inventory(float(p),float(v),seeds=True);diffs=[]
      for i in range(len(inv)):
       for j in range(i):diffs.append({'pair':[j,i],'J':inv[i]['response'][2]-inv[j]['response'][2],'M':inv[i]['response'][3]-inv[j]['response'][3]})
      out.append({'h':float(h),'z':z.tolist(),'base':[float(p),float(v)],'J_worlds':[list(w) for w in ww],'companions':inv,'differences':diffs})
      print('CHILD',h,p,v,'worlds',len(inv),'M_diffs',[d['M'] for d in diffs],flush=True)
(HERE/'positive_child_scan.json').write_text(json.dumps({'steps':steps,'data':out},indent=2)+'\n')
