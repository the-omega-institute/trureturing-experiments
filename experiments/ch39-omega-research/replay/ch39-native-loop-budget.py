"""Twelve shortest-triangle return words are necessary and sufficient to open the third pure pole."""
from fractions import Fraction as Q
from math import isqrt
from pathlib import Path
from runpy import run_path

root=Path(__file__).resolve().parent
native=run_path(str(root/"ch39-native-clock-preparations.py"))
base,B,root_cell,rounded=[native[k] for k in ("base","B","root_cell","rounded")]
model=native["model"]


def add(v,w):
    return tuple(rounded(x+y) for x,y in zip(v,w))


def scale(c,v):
    return tuple(rounded(c*x) for x in v)


def product(v,w):
    a,b=v
    c,d=w
    return rounded(a*c-B(17)*b*d),rounded(a*d+b*c)


def dot(v,w):
    return rounded(v[0]*w[0]+B(17)*v[1]*w[1])


def main():
    ids,z,p,r,k,gamma,theta=model["calibration"]()
    weights=model["base_weights"]()
    m=B(1)/sum((v/x for x,v in zip(z,p)),B(0))
    mean=sum(z,B(0))/B(len(z))
    threshold=B(Q(".00082658"))
    pivot,_,q2,_=base["physical_inputs"]()
    mass2=q2*z[2]/threshold
    mass=isqrt(mass2.hi//1)+1
    assert mass**2>mass2.hi
    parent=run_path(str(root/"ch39-phase-support-parent.py"))["parent_endpoint"]
    assert parent(mass,ids,[weights],kinetic_scale=threshold/mean)==(1,119)
    H=B(125**2)/(z[2]*q2)
    response=sum((weights[i]/(B(1)-H*zi) for i,zi in zip(ids,z)),B(0))
    assert H.hi<1/z[0].hi and ((B(1)-pivot)*response).hi<1
    word=base["word"]
    def turn(n):
        return tuple(B(Q(a,word["DEN"]),Q(b,word["DEN"])) for a,b in word["power"](n))
    moves={n:turn(-n) for n in range(-12,13)}
    step=root_cell(B(2)*(B(1)-turn(1)[0]))
    prefix=base["data"]("ch39-tick-separators.json")["prefix"]
    powers=[[-43235114615559591407865398948523*(((prefix[i-1][0]>>112)*rate+2**311)//2**312)
             for rate in base["rates"]] for i in ids]
    coef=B(k)/root_cell(r)
    g=B(gamma)*coef*coef
    diag=g*sum((v*v for v in p[2:]),B(0))
    pmax=max(v.hi for v in p[2:])
    two_registers = {(11,-1):{}, (12,0):{}}
    for name,ns in zip(("prime 2","prime 3","relative","product"),zip(*[(a,b,b-a,b+a) for a,b in powers])):
        q=[turn(n) for n in ns]
        tail=(B(0),B(0))
        for v,phase in zip(p[2:],q[2:]):
            tail=add(tail,scale(v,phase))
        first={n:scale(root_cell(p[0]),product(q[0],v)) for n,v in moves.items()}
        second={n:scale(root_cell(p[1]),product(q[1],v)) for n,v in moves.items()}
        def cost(i,j,left):
            a,b=first[i],second[j]
            v=add(a,b)
            rest=rounded(B(2)*dot(a,b)+B(2)*coef*dot(v,tail)-g*dot(tail,tail)+diag)
            slope=add(scale(coef,v),scale(B(0)-g,tail))
            radius=step*B(pmax*left)
            pad=B(2)*root_cell(dot(slope,slope))*radius+g*radius*radius
            return rounded(m*(B(1)+theta*(rest-pad)))
        for (i,j),table in two_registers.items():
            table[name]=cost(i,j,0)
        low=min(cost(i,j,11-abs(i)-abs(j)).lo for i in range(-11,12)
                for j in range(-11+abs(i),12-abs(i)))
        assert low>Q(".00086")
        if name=="relative":
            winners=[]
            for i in range(-12,13):
                for j in range(-12+abs(i),13-abs(i)):
                    left=12-abs(i)-abs(j)
                    c=cost(i,j,left)
                    if c.lo>threshold.hi:continue
                    assert left==0
                    above=isqrt((q2*z[2]/c).hi//1)+1
                    upper_side=parent(above,ids,[weights],kinetic_scale=c/mean)
                    print("Threshold test",i,j,float(c.lo),upper_side,flush=True)
                    if upper_side==(1,119):
                        continue
                    winners.append((i,j))
                    # For each surviving register certify a physical crossing above pole 3.
                    right=q2*z[2]/c
                    below=isqrt(right.lo//1)
                    assert below**2<right.lo
                    assert parent(below,ids,[weights],kinetic_scale=c/mean)==(-1,144)
                    central = 69621690000000 if i==11 else 69725810000000
                    for endpoint,sign in ((central-10000000,1),(central+10000000,-1)):
                        assert parent(endpoint,ids,[weights],kinetic_scale=c/mean)==(sign,144)
            assert winners==[(11,-1),(12,0)]
            print("Optimal relative-clock height registers:",winners)
        else:
            assert min(cost(i,j,12-abs(i)-abs(j)).lo for i in range(-12,13)
                       for j in range(-12+abs(i),13-abs(i))) > threshold.hi
    critical=B(Q(".00082657012848"),Q(".00082657012851"))
    for c,above_pole,expected in ((critical.lo,False,(-1,144)),(critical.hi,True,(1,119))):
        mass_cell=q2*z[2]/B(c)
        boundary=isqrt((mass_cell.hi if above_pole else mass_cell.lo)//1)+int(above_pole)
        assert boundary**2 > mass_cell.hi if above_pole else boundary**2 < mass_cell.lo
        assert parent(boundary,ids,[weights],kinetic_scale=B(c)/mean)==expected
    for key,table in two_registers.items():
        relative=table["relative"]
        allowed={label:(critical-relative)/(v-relative) for label,v in table.items() if label!="relative"}
        worst=allowed["prime 2"]
        assert all(worst.hi<v.lo for label,v in allowed.items() if label!="prime 2")
        if key==(11,-1):assert Q(".00074272")<worst.lo<=worst.hi<Q(".00074273")
        else:assert Q(".00200482")<worst.lo<=worst.hi<Q(".00200484")
        print("Allowed contaminant fractions",key,[(label,float(v.lo),float(v.hi)) for label,v in allowed.items()])
    w1,w2=weights[ids[0]],weights[ids[1]]
    rem=B(1)-w1-w2
    cos1,cos11,cos12=[turn(n)[0] for n in (1,11,12)]
    loss_a=B(2)*w1*(B(1)-w1)*(B(1)-cos12)
    loss_b=B(2)*(w1*w2*(B(1)-cos12)+w1*rem*(B(1)-cos11)+w2*rem*(B(1)-cos1))
    assert Q(".05592030")<loss_b.lo<=loss_b.hi<Q(".05592031")
    assert Q(".06450148")<loss_a.lo<=loss_a.hi<Q(".06450149")
    print("Projective losses (12,0) and (11,-1):",float(loss_a.lo),float(loss_b.lo))
    print("All four native preparations and every 50-height integer register of cost <=11 keep pole 3 closed.")
    print("Least total native return count: 12. There are exactly two optimal registers in the four preparations.")


if __name__=="__main__":
    main()
