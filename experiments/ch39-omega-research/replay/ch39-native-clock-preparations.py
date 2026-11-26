"""Conditional effective scale from coherent preparation by actual five-tick shell words."""
from fractions import Fraction as Q
from pathlib import Path
from runpy import run_path

root=Path(__file__).resolve().parent
model=run_path(str(root/"ch39-clock-schur-mediators.py"))
B,root_cell,base=[model[k] for k in ("B","root_cell","base")]
rounded=model["clock"]["rounded"]


def compute():
    ids,z,p,rest,k,gamma,theta=model["calibration"]()
    m=B(1)/sum((v/x for x,v in zip(z,p)),B(0))
    second=sum((v*x for v,x in zip(p,z)),B(0))
    d=B(1)-(second-m)/(B(2)*m*theta*root_cell(p[0])*(root_cell(p[1])+B(k)*root_cell(rest)))
    im=root_cell((B(1)-d*d)/B(17))
    h=(B(Q(1,3)),B(Q(10,27)))
    U=[(B(3),B(0)),(B(Q(-1,3)),B(3))]
    for v in p[2:]:
        b=B(-k)*root_cell(v/rest)
        U.append(tuple(b*x for x in h))
    diagonal=sum((v*(a*a+b*b) for v,(a,b) in zip(p,U)),B(0))
    prefix=base["data"]("ch39-tick-separators.json")["prefix"]
    word=base["word"]
    powers=[[-43235114615559591407865398948523*(((prefix[i-1][0]>>112)*r+2**311)//2**312)
             for r in base["rates"]] for i in ids]
    phases=[(a,b,b-a,b+a) for a,b in powers]
    def amplitude(n):
        return tuple(B(Q(l,word["DEN"]),Q(h,word["DEN"])) for l,h in word["power"](n))
    for label,column in zip(("prime 2","prime 3","relative","product"),zip(*phases)):
        native=[amplitude(n) for n in column]
        results=[]
        for alternate in (False,True):
            phases_at=list(native)
            if alternate:
                a,b=native[0]
                phases_at[0]=(a*d+B(17)*b*im,b*d-a*im)
            norm=B(0)
            for axis in zip(*U):
                coords=[rounded(root_cell(v)*u) for v,u in zip(p,axis)]
                real=sum((s*a for s,(a,_) in zip(coords,phases_at)),B(0))
                imag=sum((s*b for s,(_,b) in zip(coords,phases_at)),B(0))
                norm+=real*real+B(17)*imag*imag
            results.append(rounded(m*(B(1)+theta*(diagonal-norm))))
        first = [root_cell(p[1])*y + B(k)*sum((v*x for v,x in zip(p[2:],coordinate[2:])),B(0))/root_cell(rest)
                 for y,coordinate in ((native[1][0],[v[0] for v in native]),(native[1][1],[v[1] for v in native]))]
        a,b = native[0]
        factor = B(2)*m*theta*root_cell(p[0])
        s_re = rounded(factor*(a*first[0]+B(17)*b*first[1]))
        s_im = rounded(factor*(a*first[1]-b*first[0]))
        center = results[0]-s_re
        ideal = center+d*s_re-B(17)*im*s_im
        assert (ideal-results[1]).lo <= 0 <= (ideal-results[1]).hi
        yield label,results,(center,s_re,s_im)


def main():
    samples = {
        "prime 2": (41658990000000,42098110000000,119),
        "prime 3": (49511270000000,64376240000000,119),
        "relative": (54414820000000,77866230000000,144),
        "product": (52433700000000,43618000000000,119),
    }
    parent = run_path(str(root/"ch39-phase-support-parent.py"))["parent_endpoint"]
    weights = model["base_weights"]()
    ids=sorted(weights)
    prefix=base["data"]("ch39-tick-separators.json")["prefix"]
    mean=sum((model["clock"]["source"]["amb"]["upper"]["height_cell"](prefix[i-1][0]) for i in ids),B(0))/B(len(ids))
    clock_word=base["word"]
    def check(scale,center,count):
        bracket = (center-10000000, center+10000000)
        for mass,sign in zip(bracket,(1,-1)):
            assert parent(mass,ids,[weights],kinetic_scale=scale/mean)==(sign,count)
        return bracket
    for label,results,(center,a,b) in compute():
        m0,m1,n1=samples[label]
        bracket0=check(results[0],m0,119)
        bracket1=check(results[1],m1,n1)
        print(label, "original",bracket0,119,"rephased",bracket1,n1)
        if label=="relative":
            for n in range(-11,12):
                (l,h),(u,v)=clock_word["power"](n)
                c=center+a*B(Q(l,clock_word["DEN"]),Q(h,clock_word["DEN"]))-B(17)*b*B(Q(u,clock_word["DEN"]),Q(v,clock_word["DEN"]))
                assert c.lo>Q(".00082658")
            (l,h),(u,v)=clock_word["power"](12)
            c=center+a*B(Q(l,clock_word["DEN"]),Q(h,clock_word["DEN"]))-B(17)*b*B(Q(u,clock_word["DEN"]),Q(v,clock_word["DEN"]))
            limits=check(c,69725810000000,144)
            assert Q(".0008224905")<c.lo<=c.hi<Q(".0008224906")
            print("Relative clock and first-height word C^12:",limits,144)


if __name__=="__main__":
    main()
