"""Exact spectator anomaly lattice for the displayed finite SU(3) and external charge window."""
from collections import Counter
from flint import fmpz_mat

def weights(p,q):
    data=Counter()
    for i in range(q,p+q+1):
        for j in range(q+1):
            data[3*(i+j)-2*(p+2*q)]+=i-j+1
    return data

def anomaly_lattice(triality,offset):
    rows=[]
    for p in range(7):
        for q in range(7-p):
            if (p-q)%3!=triality:
                continue
            mult=weights(p,q)
            for a in range(offset-72,offset+73,36):
                spectrum={a+4*h:v for h,v in mult.items()}
                mix=sum(a*h*h*v for h,v in mult.items())
                cubic=sum(h**3*v for h,v in mult.items())
                rows.append((spectrum,mix,cubic))
    ys=sorted(set().union(*(s for s,_,_ in rows)))
    n=len(ys)
    mat=fmpz_mat([[s.get(y,0) for y in ys]+[mix,cubic] for s,mix,cubic in rows]).hnf()
    # Row Hermite form generates the entire integer row lattice. The rows
    # whose leading entry is beyond the spectrum span its exact zero fiber.
    return [[int(mat[i,n]),int(mat[i,n+1])] for i in range(mat.nrows())
            if all(mat[i,j]==0 for j in range(n)) and (mat[i,n] or mat[i,n+1])]

def check():
    rows=[]
    for tr in range(3):
        for a in range(36):
            yclass=(a+4*tr)%6
            dim=(1,6,3,2,3,6)[yclass]
            # The SM center quotient imposes these dimension divisibilities.
            rows.extend([dim*m,dim*c] for m,c in anomaly_lattice(tr,a))
    h=fmpz_mat(rows).hnf()
    basis=[row for row in h.tolist() if any(row)]
    assert basis==[[324,486],[0,1458]]
    target=(-96552,1944)
    assert target[0]%324==0 and (target[1]-486*(target[0]//324))%1458!=0
    print("enclosing lattice (A H^2,H^3) in finite spectator window",basis)
    print("required correction outside this lattice",target)

if __name__=="__main__":
    check()
