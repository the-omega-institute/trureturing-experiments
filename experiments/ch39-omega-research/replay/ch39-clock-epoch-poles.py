"""Which uniform prime-2 epoch lengths can still open the third one-loop pole?"""
from fractions import Fraction as Q
from itertools import combinations
from math import isqrt
from pathlib import Path
from runpy import run_path

root=Path(__file__).resolve().parent
source=run_path(str(root/"ch39-consecutive-dephasing.py"))
model,B,root_cell,base,word=[source[k] for k in ("model","B","root_cell","base","word")]
unit=word["DEN"]
grid=1 << 52


def shell_interval(v):
    return v.lo*unit//1,-((-v.hi*unit)//1)


def kernels(lengths,ids,p):
    ns=run_path(str(root/"ch39-native-clock-dephasing.py"))["exponents"](ids)
    phases=list(map(word["power"],ns))
    halves=[]
    for a,b in phases:
        c=B(Q(a[0],unit),Q(a[1],unit))
        x=root_cell((B(1)+c)/B(2))
        y=root_cell((B(1)-c)/B(34))
        assert b[0]*b[1]>0
        if b[1]<0:y=B(-y.hi,-y.lo)
        halves.append((shell_interval(x),shell_interval(y)))
    n=len(ids)
    mats={k:[[0]*n for _ in ids] for k in lengths}
    for i,v in enumerate(p):
        h=v.lo*grid//1
        assert v.hi*grid<h+2
        for mat in mats.values():mat[i][i]=h
    high=max(lengths)
    for i,j in combinations(range(n),2):
        h=word["multiply"](halves[i],word["conjugate"](halves[j]))
        rot=word["multiply"](phases[i],word["conjugate"](phases[j]))
        l,r=h[1]
        f=root_cell(p[i]*p[j])/B(Q(l,unit),Q(r,unit))
        a,b=shell_interval(f)
        even=word["ONE"]
        for k in range(1,high+1):
            if k%2==0:even=word["multiply"](even,rot)
            if k not in mats:continue
            v=word["multiply"](even,h)[1] if k%2 else even[1]
            nums=(v[0]*a,v[1]*a,v[0]*b,v[1]*b)
            den=k*unit**2
            floor=min(nums)*grid//den
            assert max(nums)*grid < (floor+2)*den
            mats[k][i][j]=mats[k][j][i]=floor
    return mats


def main():
    ids,z,p,*_=model["calibration"]()
    n=len(ids)
    m=B(1)/sum((v/x for v,x in zip(p,z)),B(0))
    rec=base["data"]("ch39-clock-epoch-poles.json")
    opens={int(k):v for k,v in rec["open"].items()}
    closed={int(k):v for k,v in rec["closed"].items()}
    mats=kernels(set(closed)|set(opens),ids,p)
    critical=B(Q(".00082657012848"),Q(".00082657012851"))
    parent=run_path(str(root/"ch39-phase-support-parent.py"))["parent_endpoint"]
    weights=model["base_weights"]()
    q2=base["physical_inputs"]()[2]
    mean=sum(z,B(0))/B(n)
    for c,side,verdict in ((critical.lo,-1,(-1,144)),(critical.hi,1,(1,119))):
        mass2=q2*z[2]/B(c)
        mass=isqrt((mass2.hi if side>0 else mass2.lo)//1)+int(side>0)
        assert mass**2>mass2.hi if side>0 else mass**2<mass2.lo
        assert parent(mass,ids,[weights],kinetic_scale=B(c)/mean)==verdict
    for k,rows in opens.items():
        den=1 << (2*rec["shift"])
        assert len(rows)==n and len({len(r) for r in rows})==1
        assert all(sum(x*x for x in r)<den for r in rows)
        mat=mats[k]
        total=den*sum(mat[i][i] for i in range(n))
        total+=2*sum(mat[i][j]*sum(x*y for x,y in zip(rows[i],rows[j])) for i,j in combinations(range(n),2))
        upper=Q(total+2*n*n*den,den*grid)*m.hi
        assert upper<critical.lo,(k,float(upper))
    for k,y in closed.items():
        assert len(y)==n and Q(sum(y),2**32)*m.lo>critical.hi
        mat=[row.copy() for row in mats[k]]
        for i,yi in enumerate(y):
            mat[i][i]-=yi*(grid//2**32)+2*n
        source["positive"](mat)
    good=set(closed)
    tail_start=max(opens)+1
    period=min(good)
    for k in range(1,tail_start+period):
        if any(j<k and k-j in good for j in good):good.add(k)
        assert (k in opens) != (k in good)
    assert all(k in good for k in range(tail_start,tail_start+period))
    # Adding the first closed length proves the entire tail by induction.
    print("First universally closed length:",min(good),"last length that can open:",max(opens))
    print("Earlier isolated closed lengths:",sorted(k for k in good if k<tail_start))


if __name__=="__main__":
    main()
