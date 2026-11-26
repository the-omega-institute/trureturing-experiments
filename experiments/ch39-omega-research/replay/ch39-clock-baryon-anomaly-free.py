"""Anomaly cancellation for the opposite embedding of baryon triality in the clock torus."""
from collections import defaultdict
from fractions import Fraction
from pathlib import Path
import runpy
from flint import fmpz_mat

root=Path(__file__).resolve().parent
m=runpy.run_path(str(root/"ch39-clock-baryon-anomalies.py"))
sp=runpy.run_path(str(root/"ch39-clock-baryon-spectators.py"))

def shorten(vector,kernel,cost):
    # Babai nearest plane on the dimension-weighted exact kernel. The final
    # character and anomaly sums are verified independently below.
    weighted=fmpz_mat([[v*c for v,c in zip(row,cost)] for row in kernel]).lll().tolist()
    weighted=[list(map(int,row)) for row in weighted]
    ortho=[]
    den=[]
    for row in weighted:
        z=list(map(Fraction,row))
        for u,s in zip(ortho,den):
            dot=sum(x*y for x,y in zip(row,u))/s
            z=[x-dot*y for x,y in zip(z,u)]
        ortho.append(z)
        den.append(sum(x*x for x in z))
    v=[int(x)*c for x,c in zip(vector,cost)]
    for row,u,s in reversed(list(zip(weighted,ortho,den))):
        k=round(sum(x*y for x,y in zip(v,u))/s)
        v=[x-k*int(y) for x,y in zip(v,row)]
    assert all(x%c==0 for x,c in zip(v,cost))
    return [x//c for x,c in zip(v,cost)]

def spectator(offset,desired,mult,degree=6,step=36):
    rows=[]
    labels=[]
    for p in range(degree+1):
        for q in range(degree+1-p):
            if (p-q)%3:
                continue
            w=sp["weights"](p,q)
            for a in range(offset-2*step,offset+2*step+1,step):
                spec={a+4*h:v for h,v in w.items()}
                assert all(y%6==0 for y in spec)
                rows.append((spec,sum(a*h*h*v for h,v in w.items()),sum(h**3*v for h,v in w.items())))
                labels.append((p,q,a))
    ys=sorted(set().union(*(v[0] for v in rows)))
    n=len(ys)
    mat=fmpz_mat([[s.get(y,0) for y in ys]+[mix,cubic] for s,mix,cubic in rows])
    h,u=mat.hnf(transform=True)
    assert u*mat==h
    match=[i for i in range(h.nrows()) if all(h[i,j]==0 for j in range(n)) and
           (int(h[i,n]),int(h[i,n+1]))==desired]
    assert len(match)==1
    kernel=[list(u.tolist()[i]) for i,row in enumerate(h.tolist()) if not any(row)]
    cost=[(p+1)*(q+1)*(p+q+2)//2 for p,q,a in labels]
    vector=shorten([mult*int(u[match[0],i]) for i in range(len(rows))],kernel,cost)
    return {label:v for label,v in zip(labels,vector) if v}

def check():
    fields=[("Q",6,2,3,1,1),("uc",3,1,0,-4,-1),("dc",3,1,0,2,-1),
            ("L",2,0,1,-3,0),("ec",1,0,0,6,0),("nc",1,0,0,0,0)]
    anomaly=[0]*6
    labels={}
    for name,dim,tc,tw,y,b in fields:
        h=(-b-8*y)%9
        labels[name]=h
        src={0:{(0,0):3},3:m["A"],6:m["dual"](m["A"])}[h]
        p=m["shift"](src,y)
        rank,a,_,a3,_,_,h3,ah2=m["moments"](p)
        assert rank==3
        anomaly=[v+u for v,u in zip(anomaly,[tc*a,tw*a,dim*a,dim*a3,dim*ah2,dim*h3])]
    assert anomaly[-2:]==[-24948,486]
    answer=defaultdict(int)
    for block in (spectator(6,(324,486),77),spectator(12,(0,1458),-26)):
        for key,val in block.items():
            answer[key]+=val
    post=defaultdict(int)
    total=0
    for (p,q,a),v in answer.items():
        for h,n in sp["weights"](p,q).items():
            post[a+4*h,h%9]+=v*n
            total+=abs(v)*n
            vals=[0,0,a,a**3,a*h*h,h**3]
            anomaly=[r+v*n*z for r,z in zip(anomaly,vals)]
    assert not any(post.values()) and not any(anomaly)
    assert (3*labels["Q"]+labels["L"])%9!=0
    print("residual weights",labels)
    print("spectator dimension with multiplicity",total,"nonzero signed irreps",sum(v!=0 for v in answer.values()))
    print("continuous anomaly totals",anomaly)
    return answer
if __name__=="__main__":
    check()
