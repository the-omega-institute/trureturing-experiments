"""在三壳联合见证的 N+4 荷六截止上，复算含低能初值的两项中性规范残差。"""
from pathlib import Path
from runpy import run_path

from flint import arb, ctx


def evaluate(root,cutoffs,weights):
    def lib(name):
        return run_path(str(root/("ch39-"+name+".py")))
    gauge=lib("all-field-signed-crossing")
    neutral=lib("full-neutral-free-cutoff")
    number=2_306_888_100_085+4
    # joint_system 的五个旧试验位置不参与读数，这里只读取它的物种、初值与端点。
    data=gauge["joint_system"](root,"1.176397e14",order=4,last=8,number=number)
    signs,masses=lib("rotation-portal-top-energy")["rotation_signs"](root,number)
    ctx.prec=192
    rho=(arb("246.22")/arb("1.22089e19"))**2/2
    species,_=neutral["response_rows"](data,masses,signs,rho,arb(".42995401872"))
    assert len(species)==119 and len(cutoffs)==len(weights)==19 and sum(weights).contains(1)
    p0,p1=data["p0"],data["p1"]
    log=(p1/p0).log()
    initial1,initial2,initial3=data["initial"]
    eight_thirds=arb(8)/3
    start=(5*initial1/3+initial2-eight_thirds*initial3,initial2-initial3)
    total=list(start)
    rule=gauge["angular_rule"](12,32)
    coeffs=[(kind,mass,(qq-eight_thirds*colour,mixed-colour)) for kind,mass,(qq,mixed,colour) in species]
    for x,c in zip(cutoffs,weights,strict=True):
        assert x>0
        node=[arb(0),arb(0)]
        for kind,mass,beta in coeffs:
            expo=gauge["angular_exposure"](kind,mass,x,p0,p1,log,rule)/(4*arb.pi())
            for j,b in enumerate(beta):
                node[j]-=b*expo
        for j,v in enumerate(node):
            total[j]+=c*v
    print("两项低能初值",*[v.str(15) for v in start],flush=True)
    print("第四壳 A_QQ−8A_c/3、A_QT3−A_c",*[v.str(15) for v in total],flush=True)
    assert all(not v.contains(0) for v in total)
    return total


def verify(root):
    full=run_path(str(root/"ch39-fourth-full-ckm-gravity.py"))
    cutoffs,weights=full["verify"](root,return_cutoffs=True)
    return evaluate(root,cutoffs,weights)


if __name__=="__main__":
    verify(Path(__file__).resolve().parent)
