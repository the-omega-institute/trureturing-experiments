"""保持第七、八壳四式，两组重新校准的四轮门分别精确匹配第九壳的一项引力矩。"""
from itertools import product
from pathlib import Path
from runpy import run_path

from flint import arb, arb_mat

ANGLES = ("2.5625",".96875","-1.86713344996825",".203125",
          ".734375","2.45423078856012","-1.75","2.4375",
          ".88786584120882",".51926804188806","1.59505036668143","-.0625")
ACTIVE = (2,5,8,9,10)


VACUUM_ANGLES = ("-2.83184964678738","-1.65625","-.484375","-.21875",
                 "-2.29296746286498","2.546875",".984375","2.203125",
                 "-2.05288937607464","-.16553068565108","-1.65625","-1.78616970185026")
VACUUM_ACTIVE = (0,4,8,9,11)


def experiment(root, fifth, columns, levels):
    def lib(name):
        return run_path(str(root/("ch39-"+name+".py")))
    geometry = []
    for n in (2_306_888_100_092,2_306_888_100_093,2_306_888_100_094):
        _,_,_,basis,_,weights,heat,field = lib("fifth-seventh-shortest-site-gates")["geometry"](
            root,fifth,columns,levels,number=n)
        geometry.append((basis.conjugate().transpose(),weights,heat,field))
    transport = lib("two-shell-ordered-bond-gates")["transport"]
    gravity = lib("shortest-common-complex-gate")["gravity_system"]
    def system(point, values, parameters, active):
        mats = transport(fifth,point,parameters=parameters,active=active)
        residual, jacobian = [], []
        for left,weights,heat,field in geometry:
            v,*directions = [(left*m)*columns for m in mats]
            f,j=gravity(v,directions,levels,weights,heat,field,values)
            residual.extend(f.entries())
            jacobian.extend(j.tolist())
        return residual,jacobian
    def certify(parameters,active,selected,remaining,radius,n,maximum_ratio):
        center=[arb(parameters[i]) for i in active]
        f,j = system(center,True,parameters,active)
        pre=arb_mat([j[i] for i in selected]).mid().inv()
        assert not pre.det().contains(0)
        eye=arb_mat([[int(i==j) for j in range(5)] for i in range(5)])
        offset=max(v.abs_upper() for v in (pre*arb_mat([[f[i]] for i in selected])).entries())
        box=[v+arb(0,radius) for v in center]
        _,wide=system(box,False,parameters,active)
        value=(arb(1) if remaining==4 else arb("1e-122"))+f[remaining]
        value+=sum(w*(b-c) for w,b,c in zip(wide[remaining],box,center,strict=True))
        maximum=arb(0)
        for item in product(range(n),repeat=5):
            sub=[c+radius*((2*j+1)-n)/n+arb(0,(radius/n).upper())
                 for c,j in zip(center,item,strict=True)]
            _,local=system(sub,False,parameters,active)
            err=eye-pre*arb_mat([local[i] for i in selected])
            defect=max(sum(abs(a).upper() for a in row) for row in err.tolist())
            maximum=max(maximum,defect)
        ratio=offset/radius+maximum
        print("余项、各轴格数、中心偏移、导数误差、自映比与余项读数",
              remaining,n,*[v.str(12) for v in (offset,maximum,ratio,value)], flush=True)
        assert ratio<arb(maximum_ratio)
        return box,value
    first=certify(ANGLES,ACTIVE,(0,1,2,3,4),5,arb(".00006"),3,".873")
    second=certify(VACUUM_ANGLES,VACUUM_ACTIVE,(0,1,2,3,5),4,arb(".00005"),2,".824")
    assert arb("-.01014") < first[1] < arb("-.01006")
    assert arb(".7960") < second[1] < arb(".7980")
    return first,second


def verify(root):
    source=run_path(str(root/"ch39-fifth-full-positive-joint.py"))
    return experiment(root,*source["verify"](root,return_operator="spectral"))


if __name__=="__main__":
    verify(Path(__file__).resolve().parent)
