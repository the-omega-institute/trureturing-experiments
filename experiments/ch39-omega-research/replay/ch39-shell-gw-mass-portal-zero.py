"""A single signed bulk portal can cancel the complete transported light quadratic coefficient."""
from pathlib import Path
from runpy import run_path
from flint import acb, arb, arb_mat

mass=run_path(str(Path(__file__).with_name("ch39-shell-gw-full-mass-transport.py")))
base=mass["base"]
q,loop=mass["q"],mass["loop"]
unit=mass["unit"]


def flow_matrix(classes,quart,gs,ym,gradient=None):
    n=len(classes)
    G=arb_mat(n,n)
    damp=q(3,2)*(3*gs[1]+q(3,5)*gs[0])
    for i,num in classes.items():
        G[i,i]=12*quart["l",i]
        if gradient is None:
            G[i,i]+=-damp+(q(0).union(6*ym) if i==0 else 0)
        if num>1:
            G[i,i]+=2*(num-1)*mass["portal"](quart,i,i)
        for j in classes:
            if i!=j:
                G[i,j]=2*num*mass["portal"](quart,i,j)
    return G


def advance_matrix(G,dt):
    n=G.nrows()
    ident=arb_mat([[int(i==j) for j in range(n)] for i in range(n)])
    P,term=ident,ident
    norm=max(sum(abs(G[i,j]).upper() for j in range(n)) for i in range(n))
    bound=dt.upper()*norm
    for k in range(1,5):
        term=dt*G*term/k
        P+=term
    err=arb(0,(bound**5*bound.exp()/120).upper())
    return P+arb_mat(n,n,[err]*(n*n))


def calculate(zeta):
    events=[]
    base["check"](bulk_portal=zeta,mass_path=events)
    values=[q("-1e6").union(q("1e6"))/unit]
    derivative=[q(0)]
    response=[q(1)]
    R={}
    for op,*args in events:
        if op=="threshold":
            k,classes,quart,M=args
            g=1 if k<=2 else 2
            if len(values)<len(classes):
                values.append(q(0))
                response.append(q(0))
                derivative.append(q(0))
            for key in quart:
                if key not in R:
                    R[key]=q(int(key==("a",(0,2))))
            for i,n in classes.items():
                num=n-int(i==g)
                if num:
                    values[i]+=num*mass["portal"](quart,i,g)*M**2/(unit*loop)
                    derivative[i]+=num*mass["portal"](R,i,g)*M**2/(unit*loop)
            values[g]+=M**2/unit
            continue
        classes,quart,gs,ym,dt=args
        # The quartic beta function is at most quadratic in quart, so its
        # imaginary part at quart+i*R is the exact directional derivative.
        pad=q(".018")
        probe={j:acb(v,R[j]+arb(0,pad)) for j,v in quart.items()}
        dR={j:v.imag for j,v in base["rates"](classes,probe,gs,ym).items()}
        assert all(abs(dt*v)<pad for v in dR.values())
        Rt={j:v+q(0).union(dt)*dR[j] for j,v in R.items()}
        J=flow_matrix(classes,quart,gs,ym)
        H=flow_matrix(classes,Rt,gs,ym,gradient=True)
        n=len(classes)
        block=arb_mat([[J[i%n,j%n] if i//n==j//n else H[i%n,j%n] if i>=n else 0
                        for j in range(2*n)] for i in range(2*n)])
        advance=advance_matrix(block,dt)
        both=(advance*arb_mat([[v] for v in values+derivative])).entries()
        values,derivative=both[:n],both[n:]
        response=(advance_matrix(J,dt)*arb_mat([[v] for v in response])).entries()
        R={j:v+dt*dR[j] for j,v in R.items()}
    return values[0]*unit,derivative[0]*unit,response[0]


if __name__=="__main__":
    left,right=q("-.009"),q("-.001")
    lower,_,_=calculate(left)
    upper,_,_=calculate(right)
    print("signed mass faces",lower.str(6),upper.str(6))
    assert lower<q("-1.2e23") and upper>q("1.4e23")
    _,slope,response=calculate(left.union(right))
    print("mass derivative per unit portal",float(slope.lower()),float(slope.upper()),
          "sensitivity to low mass",float(response.lower()),float(response.upper()))
    assert slope>q("1e25") and response>0
