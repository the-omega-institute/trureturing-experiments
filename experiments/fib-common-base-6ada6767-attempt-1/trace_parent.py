from discover import *
from full_parent import world,secant,pair_worlds,valid

def trace(terminal,z,hs):
    data=[]
    for h in hs:
        f=lambda zz:secant(np.r_[zz,h],terminal)/np.array([30,30,100,200,500])
        sol=root(f,z,jac=lambda zz:jac(f,zz),tol=1e-10)
        err=float(np.max(np.abs(f(sol.x))))
        if not valid(np.r_[sol.x,h]) or err>1e-7:
            print('STOP',terminal,h,sol.message,err,sol.x,flush=True);break
        z=sol.x
        ww=pair_worlds(np.r_[z,h]);ans=[world(*w) for w in ww]
        data.append({'h':float(h),'z':z.tolist(),'base':ans[0][:2].tolist(),'endpoints':[w.tolist() if hasattr(w,'tolist') else w for w in ww],'responses':[w.tolist() for w in ans],'error_scaled':err})
    return data
if __name__=='__main__':
    raw=json.loads((HERE/'curves.json').read_text());out={}
    for name,terminal,hs in [('J',2,np.linspace(.001,.099,500)),('M',3,np.linspace(.11875169891073695,.492,500)),('M_back',3,np.linspace(.11875169891073695,.0001,200))]:
        rr=raw['J' if terminal==2 else 'M'][0];R,B,k=rr['z'];h=rr['h'];ww=rr['endpoints'];ds=[-response(r/b,b)[-1]*0 for r,b in []]
        def drift(r,b):
            jj=b+r*(1-r)/b;nn=jj+r;T=1+b;U=1+2*b;Lam=1+6/T+1/U-4/jj-4/nn
            return 24-6*(8*b-1)*(5*b+2)/(5*T*U)/Lam
        ds=[drift(*w) for w in ww];z=np.r_[R,B,k,(ds[0]+ds[1])/2,(ds[1]-ds[0])/(2*h)]
        out[name]=trace(terminal,z,hs)
        d=out[name];print(name,'count',len(d),'last',d[-1] if d else None,flush=True)
    (HERE/'parent_curves.json').write_text(json.dumps(out,indent=2)+'\n')
