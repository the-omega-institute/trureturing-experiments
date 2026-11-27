"""每条边施一次按长度共用的两站酉门，认证第七与第八壳的四个引力条件。"""
from itertools import combinations
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, arb_mat


LAYERS = ((0,((6,11),(7,12),(8,18))),
          (1,((0,2),(1,3),(4,5))),
          (2,((6,9),(7,15),(8,16))),
          (2,((6,13),(7,10),(8,17))))
ANGLES = ("-1.75","1","-2.62014134706321","-.5",
          "-.60479820190850","-.5","1.25","-2",
          "-1.63282812895125","1.5",".75","-2.62411886066799")
ACTIVE = (2,4,8,11)


def block(angles):
    theta,alpha,beta,gamma = angles
    s,c = theta.sin_cos()
    imaginary = acb(0,1)
    phases = [(imaginary*t).exp() for t in (gamma+alpha,gamma-beta,gamma+beta,gamma-alpha)]
    values = [u*v for u,v in zip(phases,(c,-s,s,c),strict=True)]
    derivatives = [[u*v for u,v in zip(phases,(-s,-c,c,-s),strict=True)]]
    for row in ((1,0,0,-1),(0,-1,1,0),(1,1,1,1)):
        derivatives.append([imaginary*v*k for v,k in zip(values,row,strict=True)])
    return values,derivatives


def transport(fifth, point, *, parameters=ANGLES, active=ACTIVE):
    angles = list(map(arb,parameters))
    for key,val in zip(active, point,strict=True):
        angles[key] = val
    gates = [block(angles[4*j:4*j+4]) for j in range(3)]
    # 在第五壳谱标架的十九个列向量上同时传播 U 及其四个导数。
    mats = [acb_mat(fifth)]+[acb_mat(19,19) for _ in active]
    for index,pairs in LAYERS:
        gate,jac = gates[index]
        for a,b in pairs:
            for j in range(19):
                src = [(m[a,j],m[b,j]) for m in mats]
                for k,(m,(u,v)) in enumerate(zip(mats,src,strict=True)):
                    d1,d2 = gate[0]*u+gate[1]*v,gate[2]*u+gate[3]*v
                    if k and active[k-1]//4==index:
                        dg = jac[active[k-1]%4]
                        l,r = src[0]
                        d1 += dg[0]*l+dg[1]*r
                        d2 += dg[2]*l+dg[3]*r
                    m[a,j],m[b,j] = d1,d2
    return mats


def experiment(root, fifth, columns, levels, *, parameters=ANGLES, active=ACTIVE,
               half_width=".00005", maximum_ratio=".656"):
    def lib(name):
        return run_path(str(root/("ch39-"+name+".py")))
    assert len(parameters)==12 and len(active)==len(set(active))==4
    sites = lib("three-ray-chiral-boundary")["geo"]["sites"]
    expected = set()
    for i,j in combinations(range(19),2):
        d = sum((a-b)**2 for a,b in zip(sites[i],sites[j],strict=True))
        if d<=52:
            expected.add((d,i,j))
    used = [( (12,32,52)[k],a,b) for k,edges in LAYERS for a,b in edges]
    assert len(used)==len(set(used))==12 and set(used)==expected
    geometry = []
    for n in (2_306_888_100_092,2_306_888_100_093):
        _,_,_,basis,_,weights,heat,field = lib("fifth-seventh-shortest-site-gates")["geometry"](
            root,fifth,columns,levels,number=n)
        geometry.append((basis.conjugate().transpose(),weights,heat,field))
    gravity = lib("shortest-common-complex-gate")["gravity_system"]
    def system(point,values):
        mats = transport(fifth,point,parameters=parameters,active=active)
        fs,js = [],[]
        for left,w,heat,field in geometry:
            p,*deriv = [(left*x)*columns for x in mats]
            f,j=gravity(p,deriv,levels,w,heat,field,values)
            fs.extend(f.tolist())
            js.extend(j.tolist())
        return arb_mat(fs),arb_mat(js)
    center=[arb(parameters[i]) for i in active]
    f,j=system(center,True)
    pre=j.mid().inv()
    assert not pre.det().contains(0)
    eye=arb_mat([[int(i==j) for j in range(4)] for i in range(4)])
    move=max(v.abs_upper() for v in (pre*f).entries())
    radius=arb(half_width)
    box=[v+arb(0,radius) for v in center]
    _,jac=system(box,False)
    norm=max(sum(abs(v).upper() for v in row) for row in (eye-pre*jac).tolist())
    ratio=move/radius+norm
    print("每边一次的联合盒半径、预条件中心位移、导数偏差、自映比",
          *[v.str(14) for v in (radius,move,norm,ratio)],flush=True)
    assert ratio < arb(maximum_ratio)
    return box,ratio


def verify(root):
    source=run_path(str(root/"ch39-fifth-full-positive-joint.py"))
    return experiment(root,*source["verify"](root,return_operator="spectral"))


if __name__=="__main__":
    verify(Path(__file__).resolve().parent)
