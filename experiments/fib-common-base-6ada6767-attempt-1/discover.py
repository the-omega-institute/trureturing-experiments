# Publication repair (2026-10-23): resolve original suppliers within this checkout.
"""Free-base collision discovery, using XI 6.2 and the existing source implementation.
Numerical data are candidates only; no root is certified by this program.
"""
from pathlib import Path
import runpy,json,time
import numpy as np
from scipy.optimize import root
HERE=Path(__file__).resolve().parent
REUSED=(Path(__file__).resolve().parents[1]/'fib-self-calibration/fixed-b-and-negative-child/negative-j-seed.py')
old=runpy.run_path(str(REUSED)); response=old['response']
def source(r,b): return np.array(response(r/b,b)[:4])
def jac(fun,x):
    x=np.asarray(x,dtype=float); ret=[]
    for k in range(len(x)):
        z=x.astype(complex);z[k]+=1e-25j;ret.append(np.imag(fun(z))/1e-25)
    return np.array(ret).T

def pair(z,h,terminal):
    R,B,k=z; w0=source(R-h,B-k*h);w1=source(R+h,B+k*h)
    return (w1-w0)[[0,1,terminal]]/(2*h)
def valid(z,h):
    R,B,k=z
    return 0<R-h<R+h<1 and min(B-k*h,B+k*h)>0 and np.all(np.isfinite(z))
def trace(terminal,z,h0,hs):
    data=[]
    for h in hs:
        f=lambda zz:pair(zz,h,terminal)/np.array([100,200,500])
        sol=root(f,z,jac=lambda zz:jac(f,zz),tol=1e-10)
        err=float(np.max(np.abs(f(sol.x))))
        if not valid(sol.x,h) or err>1e-7: break
        z=sol.x
        R,B,k=z;ww=[source(R-h,B-k*h),source(R+h,B+k*h)]
        data.append({'h':float(h),'z':z.tolist(),'base':ww[0][:2].tolist(),'endpoints':[[float(R-h),float(B-k*h)],[float(R+h),float(B+k*h)]],'responses':[w.tolist() for w in ww],'error_scaled':err})
    return data
if __name__=='__main__':
    jm=np.array([.8995309197957799,.12856648200849686,-1.218941298820346]);
    mm=np.array([(.4+.6375033978214739)/2,(.040782959048828384+.0409169842106061)/2,(.0409169842106061-.040782959048828384)/(.6375033978214739-.4)])
    out={'reuse':str(REUSED),'J':trace(2,jm,.001,np.linspace(.001,.1,300)), 'M':trace(3,mm,.11875169891073695,np.linspace(.11875169891073695,.49,500)),'M_back':trace(3,mm,.11875169891073695,np.linspace(.11875169891073695,.0001,200))}
    (HERE/'curves.json').write_text(json.dumps(out,indent=2)+'\n')
    for name in ['J','M','M_back']:
        d=out[name];print(name,'count',len(d),'first',d[0] if d else None,'last',d[-1] if d else None,flush=True)
