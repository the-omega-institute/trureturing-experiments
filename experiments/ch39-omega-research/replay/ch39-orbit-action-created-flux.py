"""A self-selected, cyclically symmetric stiffness orbit creates three flux units."""
from pathlib import Path
import runpy
from flint import acb,arb,acb_mat

work=Path(__file__).resolve().parent
fit=runpy.run_path(str(work/"ch39-hierarchy-in-zero-flux.py"))
flux=runpy.run_path(str(work/"ch39-source-created-flux.py"))
s=fit["s"]


def verify_zero(n,amount,selected):
    height=acb.zeta_zero(n).imag
    radius=arb(0,arb(".000001"))
    model=fit["model"](s["a"]*height+radius,s["b"]*height+radius)
    m,orbits=model
    g=m.inv()
    powers=[]
    for orbit in orbits:
        coeff=acb_mat([[v[i,0] for v in orbit] for i in range(19)])
        response=amount*coeff.conjugate().transpose()*g*coeff
        for i in range(3):
            response[i,i]+=1
        powers.append(response.det().real)
    gap=min((val-powers[selected]).lower() for j,val in enumerate(powers) if j!=selected)
    assert gap>0
    _,_,old=fit["evaluate"](model,[arb(0)]*4,(0,0))
    _,_,new=fit["evaluate"](model,[amount*int(j==selected) for j in range(4)],(0,0))
    assert flux["degrees"](old)[0]==(0,[0])
    result,margin=flux["degrees"](new)
    assert result==(3,[3]) and margin>0
    values=sorted((abs(z) for z in fit["symbols"](-fit["block"](new,fit["ii"],fit["oo"]))),
                  key=lambda x:float(x.mid()))
    assert values[0]<values[1]<values[2]
    print("zero, selected orbit, determinant ratio lead, angular margin, middle/heavy",
          n,selected,gap.str(12),margin.str(12),(values[1]/values[2]).str(9))


if __name__=="__main__":
    verify_zero(10,arb(".05"),3)
    verify_zero(527,arb(1),2)
