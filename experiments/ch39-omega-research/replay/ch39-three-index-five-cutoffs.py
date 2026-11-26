"""五个可自由放置的正截止节点精确满足三个相邻实际零点的七项引力矩。"""
from pathlib import Path
from runpy import run_path
from flint import arb, arb_mat


def column(values):
    return arb_mat([[v] for v in values])


def verify(root):
    previous=run_path(str(root/"ch39-three-index-free-atoms.py"))
    moment,fields=previous["spectra"](root)
    slope=previous["column_slope"]
    target=column([1,1,arb("1e-122"),0,0,0,0])
    x=list(map(arb,("1.2206072485718666e-10","1.6223370759615275e-9",
                    "1.67474638503468617727543046e-8","9.01879755551460057162233542e-8",
                    ".5234083107326523")))
    w=list(map(arb,("6741408.15073752042486947927","-17893364.2379894603848872284",
                    "11224564.6228782044526679114","-72607.5356262519173580395030",
                    "-1.25752921227472164202739073e-8")))
    variable=(2,3)
    scale=[abs(v) for v in w]
    def jac(values,derivs,weights):
        cols=[[s*a for a in v] for s,v in zip(scale,values)]
        cols.extend([weights[j]*x[j]*a for a in derivs[j]] for j in variable)
        return arb_mat(cols).transpose()
    v=[moment(a) for a in x]
    dv={j:slope(fields,x[j]) for j in variable}
    J=jac(v,dv,w)
    center=arb_mat([[a.mid() for a in row] for row in J.tolist()])
    pre=center.inv()
    assert not pre.det().contains(0)
    res=pre*(arb_mat(v).transpose()*column(w)-target)
    radius=arb("1e-23")
    unit=arb(0,radius.upper())
    x_box=x.copy()
    for j in variable:
        x_box[j]=x[j]*(1+unit)
    w_box=[a+s*unit for a,s in zip(w,scale)]
    dv_box={j:slope(fields,x_box[j]) for j in variable}
    v_box=v.copy()
    for j in variable:
        v_box[j]=[a+x[j]*unit*b for a,b in zip(v[j],dv_box[j])]
    box_jac=jac(v_box,dv_box,w_box)
    identity=arb_mat([[int(i==j) for j in range(7)] for i in range(7)])
    defect=identity-pre*box_jac
    rate=max(sum(abs(a).upper() for a in row) for row in defect.tolist())
    residual=max(abs(a).upper() for a in res.entries())
    print("预乘中心残差、全盒压缩率、坐标半宽",residual.str(16),rate.str(16),radius.str(3))
    assert rate<1 and residual+rate*radius<radius
    assert x_box[0]>0 and x_box[-1]<1 and all(a<b for a,b in zip(x_box,x_box[1:]))
    assert all(s*a>0 for s,a in zip((1,-1,1,-1,-1),w_box))
    cost=-w_box[1]-w_box[3]-w_box[4]
    print("五个实际截止节点",*[a.str(22) for a in x_box])
    print("相应的五个非零权重",*[a.str(22) for a in w_box])
    print("负质量",cost.str(22))
    assert arb("17965971.7736157")<cost<arb("17965971.7736158")


if __name__=="__main__":
    verify(Path(__file__).resolve().parent)
