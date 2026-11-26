"""A zero UV mass eigenvalue with an explicitly nonzero finite polynomial shell vertex."""
from pathlib import Path
from runpy import run_path
from flint import arb, arb_mat

folder=Path(__file__).resolve().parent
signed=run_path(str(folder/"ch39-shell-gw-mass-portal-zero.py"))
base=signed["base"]
mass=signed["mass"]
q,loop,unit=signed["q"],signed["loop"],signed["unit"]


def reference_step(classes,quart,gs,ym,dt,values,minimum,maximum):
    n=len(classes)
    H=signed["flow_matrix"](classes,quart,gs,ym)
    big=arb_mat(3*n,3*n)
    for i,num in classes.items():
        for j in classes:
            big[i,j]=H[i,j]
        same=mass["portal"](quart,i,i) if num>1 else q(0)
        for r in (n,2*n):
            big[r+i,r+i]=H[i,i]-2*num*same
            for j in classes:
                big[r+i,j]=2*(same if j==i else mass["portal"](quart,i,j))
    state=values+minimum+maximum
    scale=max(abs(v).upper() for v in state)
    norm=max(sum(abs(big[i,j]).upper() for j in range(3*n)) for i in range(3*n))
    move=scale*((dt.upper()*norm).exp()-1)
    frob=2*sum(num*(max(abs(minimum[i]).upper(),abs(maximum[i]).upper())+move)**2
               for i,num in classes.items()).sqrt()
    out=(signed["advance_matrix"](big,dt)*arb_mat([[v] for v in state])).entries()
    values[:],minimum[:],maximum[:]=out[:n],out[n:2*n],out[2*n:]
    return frob


def face(zeta,power):
    gain=q(81,977)**power
    events=[]
    T,ZF,A,v0,raw,time=base["check"](gain,mass_path=events,bulk_portal=zeta)
    Y=(12*T).sqrt()
    YH=arb(12).sqrt()*raw.sqrt()*gain*(6*T*time).exp()
    ceiling=q(".020")
    Zrate=12*A+3*ceiling+4*Y**2
    Zforce=(8*Y*ZF+96*Y**3)*YH
    Zmax=Zforce*((Zrate*time).exp()-1)/Zrate
    assert Zmax<ceiling and Zmax<6*v0
    Aq=2*A+Zmax+2*Y**2
    t,E=q(0),q(0)
    values=[q("-1e6").union(q("1e6"))/unit]
    minimum,maximum=values.copy(),values.copy()
    for op,*args in events:
        if op=="threshold":
            k,groups,quart,M=args
            g=1 if k<3 else 2
            m=M**2/unit
            if len(values)<len(groups):
                values.append(q(0))
                minimum.append(m)
                maximum.append(m)
            for i,num in groups.items():
                old=num-int(i==g)
                if old:
                    inc=mass["portal"](quart,i,g)*m/loop
                    values[i]+=old*inc
                    minimum[i]+=inc
                    maximum[i]+=inc
            values[g]+=m
            minimum[g]=arb.min(minimum[g],m)
            maximum[g]=arb.max(maximum[g],m)
            continue
        groups,quart,gs,ym,dt=args
        size=reference_step(*args,values,minimum,maximum)
        if len(groups)==3:
            t+=dt.upper()
            errZ=Zforce*((Zrate*arb.min(t,time.upper())).exp()-1)/Zrate
            # Euclidean mass perturbations are measured in 10^24 GeV².
            forcing=(errZ+4*Y*YH)*size
            step=(Aq*dt.upper()).exp()
            E=(E*step+forcing*(step-1)/Aq).upper()
    print("portal",zeta.str(4),"power",power,"quartic error",Zmax.upper().str(7),
          "mass error in 1e24",E.upper().str(7),"auxiliary light mass",values[0].str(6),
          "min heavy",min(v.lower() for v in minimum[1:]).str(7))
    return values[0],min(v.lower() for v in minimum[1:]),E,max(v.upper() for v in maximum)


if __name__=="__main__":
    lo,hi=q("-.009"),q("-.001")
    left,hleft,err_l,_=face(lo,5)
    right,hright,err_r,_=face(hi,5)
    assert left+err_l<q("-.095") and min(right.lower(),hright)-err_r>q(".11")
    _,bulk,error,maximum=face(lo.union(hi),5)
    assert bulk-error>q(".49") and error/bulk<q(".055")
    flow=base["flow"]
    shell=flow["full"]["model"]
    sing,family,_,_=shell["m"]["transported_spectra"]()
    ds=flow["full"]["all_levels"](sing,family)
    width=q("12e-8")
    positive_ratio=(ds[0]+width)/(ds[-1]-width)
    reverse_ratio=(ds[-1]-ds[-2]+2*width)/(ds[-1]-width)
    high_ratio=(bulk-error)/(maximum+error)
    assert positive_ratio<q(".000307") and reverse_ratio<q(".002")
    for ratio in (positive_ratio,reverse_ratio):
        absolute=(bulk-error-ratio*(maximum+error))/(2*(1+ratio))
        fractional=(high_ratio-ratio)/(2*(1+ratio))
        print("native end gap ratio above bounded by",ratio.upper().str(7),"affine approximation errors exceed",
              absolute.lower().str(8),"in units of 1e24; relative to spectral span",fractional.lower().str(8))
    assert high_ratio>q(".0273")
    print("largest mass squared in 1e24 below",(maximum+error).upper().str(8))
    print("heavy compression lower in units of 1e24",(bulk-error).lower().str(8),
          "upper orthogonal amplitude of a tuned zero",(error/bulk).upper().str(8))
