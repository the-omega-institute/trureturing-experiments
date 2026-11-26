"""Finite inverse temperatures sufficient to recover the nearest zeta-zero label."""
from pathlib import Path
from runpy import run_path
from flint import arb, acb, ctx

root = Path(__file__).parent
thermal = run_path(str(root/"ch39-fermion-addition-temperature.py"))
gauge = run_path(str(root/"ch39-floating-scale-winner.py"))
s = thermal["finite"]
weight, hardy = s["weight"], s["hardy"]


def below(lo,hi,beta,coherent,ceiling):
    stack = [(lo,hi)]
    epsilon = arb("3.65e-22")
    while stack:
        a,b = stack.pop()
        g = a.union(b)
        w = weight(g)
        if w<ceiling:
            continue
        if beta:
            t = max(arb(0),abs(hardy(g)).lower()-epsilon)
            h = thermal["free_energy"](t,beta)
            if w/(1+thermal["load"](h,coherent))<ceiling:
                continue
        assert b-a>arb("1e-9")
        mid=(a+b)/2
        stack.extend([(a,mid),(mid,b)])


def check(index,low,high,coherent,first_beta,witness):
    left,zero,right = [z.imag for z in acb.zeta_zeros(index-1,3)]
    left,right = (left+zero)/2,(right+zero)/2
    epsilon = arb("3.65e-22")
    ref = weight(zero)/(1+thermal["load"](epsilon,coherent))
    below(arb(low),left.upper(),first_beta,coherent,ref)
    if right.lower()<high:
        below(right.lower(),arb(high),first_beta,coherent,ref)
    assert s["source"]["modulus"](-arb(59)/108)/(arb(high)**2+arb(1)/4) < ref
    if witness is not None:
        t = abs(hardy(witness))+epsilon
        comparisons = [(first_beta/2,0)]
        if not coherent:
            comparisons = [(arb(16),arb(16)),(arb(8),0)]
        elif index==10854:
            comparisons = [(arb(1)/d,arb(1)/d) for d in (2,4,8,16,32,64)]
            comparisons.append((arb(1)/128,0))
        for b,max_inside_beta in comparisons:
            rival = weight(witness)/(1+thermal["load"](thermal["free_energy"](t,b),coherent))
            below(left.lower(),right.upper(),max_inside_beta,coherent,rival)
    print(index,"coherent" if coherent else "trace","nearest-zero cell guaranteed from beta",first_beta)


def feasible_lower_competitor():
    q = gauge["q"]
    base = q("26.861182485")
    old = acb.zeta_zero(10854).imag
    height = q("10593.748")
    offset = q(gauge["fixed"]["OFFSET"])
    drift = q("-.0000095")
    time = base+((old**2+offset)/(height**2+offset)).log()/2+drift
    common = gauge["target"](10854,base)*(-2*drift).exp()
    logs = list(map(q,("-25.4446915838 -12.9928374813 -1.5528331663 -23.8189018241 "
                       "-17.8423094899 -10.1208976879 -25.6505068564 -14.9499484326 "
                       "-9.2850641599").split()))
    inv,y,_ = gauge["evolve"]([common*r for r in gauge["ratios"]],logs,-time,131072,diagonal_quarks=True)
    assert 23<time<33
    assert all(abs(1/(v.sqrt()*c)-1)<q(".001020") for v,c in zip(inv,gauge["centers"]))
    assert all(q(".99")<(v/2).exp()/c<q("1.01") for v,c in zip(y,gauge["yukawa_centers"]))
    return height


def verify():
    ctx.prec,ctx.cap = 128,2
    first_rival = arb("10617.645332791202380").union(arb("10617.645332791202395"))
    check(10875,10603,10622,False,arb(32),first_rival)
    check(10875,10603,10622,True,arb(1),first_rival)
    check(10854,10567,10613,True,arb(1),feasible_lower_competitor())


if __name__ == "__main__":
    verify()
