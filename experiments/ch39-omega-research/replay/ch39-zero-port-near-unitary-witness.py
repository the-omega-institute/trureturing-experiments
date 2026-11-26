"""保持三分之二中心权重与精确规范均值时，连续减小第五项味读数的缺口。"""
from pathlib import Path
from runpy import run_path

from flint import acb, arb, arb_poly

family = run_path(str(Path(__file__).with_name("ch39-zero-port-support-switches.py")))
box, study, sharp, prev, edges = (family[k] for k in ("box","study","sharp","prev","edges"))


def verify():
    cost,strength,_,*unused = prev["setup"]()
    ids = (11431,11432,11415,11430,11372,11304,50,60495)
    corner = [c+sgn*arb(2)**-28 for c,sgn in zip(box["centers"],(1,-1,1,-1))]
    coords,prices,rows = [],[],[]
    for n in ids:
        height = acb.zeta_zero(n).imag
        x = arb(0) if n==11431 else study["center"]/(height**2+study["offset"])-1
        coords.append(x)
        prices.append(cost(x)[0] if n!=11431 else arb(0))
        rows.append(box["flavor_data"](study["native"]["plane_response"](height),corner)[0])
    r = arb(1)/3
    first,_ = edges["candidate"](coords[:6],rows[:6],r)
    low = -r*coords[-1]/(coords[-2]-coords[-1])
    assert 0<low<r
    last = [1-r,0,0,0,0,0,low,r-low]
    first = [*first,0,0]
    means = []
    for column in zip(*rows):
        start,finish = (family["dot"](weights,column) for weights in (first,last))
        means.append(arb_poly([start,finish-start]))
    a,b,u,d,j = means
    one = 1-u
    # 实、虚部使用同一正分母。半径平方以同一多项式形式计算。
    re = one**2*(u-d)-u*one*(a+b)+(1+u)*a*b
    im = 2*one*j
    lo_phase = im-arb("1.141").tan()*re
    hi_phase = im-arb("1.153").tan()*re
    band = arb(0).union(1)
    assert lo_phase.derivative()(band)>0 and hi_phase(band)<0
    assert re(band)>0 and im(band)>0
    for col,small,large in ((a,".224955",".225065"),(b,".041775",".041885")):
        assert (col-arb(small)**2*one)(1)>0 and (col-arb(large)**2*one)(band)<0
    assert sharp["lo"]<u(band)<sharp["hi"]
    # 乘以 4(1-U)^4 后，R² 的分母恰为 4ABU(1-U-A)(1-U-B)。
    numerator = re**2+im**2
    denominator = 4*a*b*u*(one-a)*(one-b)
    assert denominator(band)>0
    derivative = numerator.derivative()*denominator-numerator*denominator.derivative()
    assert derivative(band)>0,derivative(band).str(9)
    ratio = numerator(1)/denominator(1)
    price = strength*family["dot"](last,prices)
    print("50, 60495 absolute weights",low.str(18),(r-low).str(18))
    print("Port deficit",price.str(19),"radius squared",ratio.str(19),"unit deficit",(1-ratio).str(15))
    print("Endpoint sines",*((v(1)/one(1)).sqrt().str(13) for v in (a,b)),u(1).sqrt().str(13),
          "CP phase",acb(re(1),im(1)).arg().str(13))
    assert 0<1-ratio<arb("1.7e-7") and price>strength*family["dot"](first,prices)


if __name__=="__main__":
    verify()
