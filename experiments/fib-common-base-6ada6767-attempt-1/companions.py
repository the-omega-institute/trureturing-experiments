"""Numerical companion enumeration by the existing IX/X degree-18 carrier.
The carrier is necessary; only residual-checked strict reconstructed sources are retained.
This program makes no completeness/exclusion claim from floating-point roots.
"""
from discover import *
from scipy.optimize import root
co=json.loads((HERE/'inverse_coefficients.json').read_text())
def coeffs(P,p,v):
    out=np.zeros(max(e[0] for e,a in P)+1)
    for (i,j,k),a in P:out[i]+=a*p**j*v**k
    return out
def peval(P,x,p,v):return np.polynomial.polynomial.polyval(x,coeffs(P,p,v))
def inverse(x,d,p,v):
    h3=(24-d)/x;h4=18-d;h5=(p-d)/(x+1);h6=(v-d)/(x+2)
    H0=20*h3-45*h4+36*h5-10*h6;H1=10*h3-20*h4+15*h5-4*h6;H2=4*h3-6*h4+4*h5-h6
    B0=294/5-2*d-H2;U=(-d*H0+B0*H2)/(d*d+B0*B0);b=(U-1)/2;j=U*(24-d)/H1;r=x*U-j
    return r,b

def inventory(p,v,seeds=True):
    cs=coeffs(co['R'],p,v);xs=np.roots(cs[::-1]);worlds=[];guesses=[]
    for x in xs:
        if abs(x.imag)<1e-4 and x.real>0:
            x=x.real;a=peval(co['alpha'],x,p,v);z=peval(co['zeta'],x,p,v)
            if a:
                r,b=inverse(x,-z/a,p,v)
                if 0<r<1 and b>0:guesses.append([r,b])
    if seeds:guesses.extend([[r,b] for r in np.linspace(.001,.999,16) for b in np.geomspace(1e-5,100,20)])
    for rb in guesses:
        f=lambda z:(source(*z)[:2]-[p,v])/np.array([100,200])
        try:
            sol=root(f,rb,jac=lambda z:jac(f,z),tol=1e-11,options={'maxfev':500});r,b=sol.x
            if 0<r<1 and b>0 and np.max(np.abs(f(sol.x)))<1e-8 and not any(np.linalg.norm(sol.x-np.array(w['source']))<1e-5 for w in worlds):
                worlds.append({'source':[float(r),float(b)],'response':source(r,b).tolist(),'error':float(np.max(np.abs(f(sol.x))))})
        except (ValueError,ZeroDivisionError,OverflowError):pass
    return worlds
if __name__=='__main__':
    bases=[]
    cur=json.loads((HERE/'parent_curves.json').read_text())
    for key in ['J','M']:
        for row in cur[key][::max(1,len(cur[key])//12)]:bases.append((key,row['h'],row['base']))
    cat=json.loads((HERE/'pair_catalogue.json').read_text())
    for key,t in [('J',2),('M',3)]:
        for row in cat[key]:
            if max(abs(row['responses'][0][i]-row['responses'][1][i]) for i in [0,1,t])<1e-6 and min(w[1] for w in row['endpoints'])>1e-4:bases.append((key,row['h'],row['base']))
    out=[]
    for key,h,(p,v) in bases:
        ww=inventory(p,v)
        diffs=[]
        for i in range(len(ww)):
            for j in range(i):diffs.append({'pair':[j,i],'J':ww[i]['response'][2]-ww[j]['response'][2],'M':ww[i]['response'][3]-ww[j]['response'][3]})
        row={'branch':key,'h':h,'base':[p,v],'worlds':ww,'terminal_differences':diffs};out.append(row)
        print(key,h,p,v,'worlds',len(ww),'diffs',diffs,flush=True)
    (HERE/'companion_data.json').write_text(json.dumps(out,indent=2)+'\n')
