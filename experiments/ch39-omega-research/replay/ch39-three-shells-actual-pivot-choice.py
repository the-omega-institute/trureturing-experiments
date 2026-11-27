"""在已经校准的前两步上逐壳读取三个原生荷六枢轴的牛顿项与真空项。"""
from pathlib import Path
from runpy import run_path

from flint import acb_mat, arb, ctx


def verify(root):
    def lib(name):
        return run_path(str(root/("ch39-"+name+".py")))

    n=2_306_888_100_085
    x,y,_,_,(m,h,k),_weights,_thermal=lib("fourth-shell-two-gravity-rules")["shell_data"](root)
    native=lib("native-rh-ray-pvm")
    moments,heat,fields=lib("three-index-free-atoms")["broken_fields"](root,[n+2,n+3])
    spec=lib("adjacent-pvm-transport")["spectrum"]
    shell=lib("three-ray-chiral-boundary")
    _,source=spec(shell,6,n+1)
    source=acb_mat([[source[i,17],source[i,18]] for i in range(19)])
    ctx.prec=256
    x3=x**3
    for number,transport,field in zip((n+2,n+3),(x3,y**3*x3),fields,strict=True):
        levels,frame=spec(shell,6,number)
        ov=frame.conjugate().transpose()*transport*source
        sites=[(m+(h-m)*abs(ov[i,0])**2+(k-m)*abs(ov[i,1])**2)**2 for i in range(19)]
        columns=[moments(field,heat,v) for v in sites]
        _,ordered,pivots,_=next(row for row in native["shell_weights"](root,number) if row[0]==6)
        assert all(a.overlaps(b) for a,b in zip(levels,ordered,strict=True))
        for j,pivot in enumerate(pivots,1):
            c=[v.real for v in pivot]
            assert sum(c).contains(1)
            K,W=[sum(w*v[b] for w,v in zip(c,columns,strict=True)) for b in (0,1)]
            print("零点序号、实际枢轴、K 与 W",number,j,K.str(13),W.str(13),flush=True)
            if j==1:
                assert K.contains(1) and W.contains(arb("1e-122"))
            else:
                assert not K.contains(1) and not W.contains(arb("1e-122"))


if __name__=="__main__":
    verify(Path(__file__).resolve().parent)
