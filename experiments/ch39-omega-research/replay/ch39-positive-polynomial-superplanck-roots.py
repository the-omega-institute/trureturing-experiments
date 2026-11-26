"""两种实际首枢轴的相邻幂三条件根在普朗克上限之外。"""
from pathlib import Path
from runpy import run_path

from flint import arb, arb_mat, ctx


def certify(root):
    gravity=run_path(str(root/"ch39-three-index-free-atoms.py"))
    moments,heat,(rows,)=gravity["broken_fields"](root,[2_306_888_100_085])
    native=run_path(str(root/"ch39-native-rh-ray-pvm.py"))
    direct=run_path(str(root/"ch39-neutral-mz-cutoff.py"))["low_weak_column"]
    high=run_path(str(root/"ch39-full-neutral-free-cutoff.py"))["weak_column"]
    ctx.prec=224
    mz,mp=arb("91.1876"),arb("1.22089e19")
    p0,p1=(mz/mp)**2,(arb("1.176397e14")/mp)**2
    rho=(arb("246.22")/mz)**2/2
    w=arb(".63811")**2*rho/2
    z=w+3*arb(".467453")**2*rho/10
    top=arb(".8616")**2*rho
    for q,levels,direct_weights,_ in native["shell_weights"](root):
        if q not in (-4,6):
            continue
        weights=[v.real for v in direct_weights[0]]
        r=[x/levels[-1] for x in levels[:-1]]+[arb(1)]
        degree,guess=(108,("1.19485",".74533")) if q==-4 else (62,("1.05608",".20685"))
        powers=[x**degree for x in r]
        def evaluate(kt,s,jacobian):
            k,t=kt
            shape=[y*(1-t*(1-x)) for x,y in zip(r,powers)]
            nodes=[p0*s+k*y for y in shape]
            if not jacobian:
                cols=[moments(rows,heat,x) for x in nodes]
                return arb_mat([[sum(c*col[i] for c,col in zip(weights,cols))-goal]
                                for i,goal in enumerate((1,arb("1e-122")))])
            matrix=arb_mat(2,2)
            for c,v,a,radial,y in zip(weights,nodes,shape,r,powers):
                dk=sum(u*(-m/v).exp() for _,u,m in rows)/(2*arb.pi())
                dw=-v*sum(g*(-m/v).exp() for g,_,m in rows)/(4*arb.pi())
                for i,f in enumerate((dk,dw)):
                    matrix[i,0]+=c*f*a
                    matrix[i,1]+=c*f*k*y*(radial-1)
            return matrix

        s_lo,s_hi=map(arb,(".3244972",".3244975") if q==-4 else (".3272101",".3272104"))
        center=list(map(arb,guess))
        for _ in range(3):
            change=evaluate(center,arb(".5"),True).solve(evaluate(center,arb(".5"),False))
            center=[(a-v.mid()).mid() for a,v in zip(center,change.entries())]
        box_size=arb("2e-9")
        box=[x+arb(0,box_size.upper()) for x in center]
        assert 1<box[0] and 0<box[1]<1
        pre=evaluate(center,arb(".5"),True).inv()
        whole=s_lo.union(s_hi)
        defect=arb_mat([[1,0],[0,1]])-pre*evaluate(box,whole,True)
        residual=pre*evaluate(center,whole,False)
        for i in range(2):
            assert abs(residual[i,0])+box_size*sum(abs(defect[i,j]) for j in range(2))<box_size
        def weak(s):
            k,t=box
            nodes=[p0*s+k*y*(1-t*(1-x)) for x,y in zip(r,powers)]
            total=arb(0)
            for c,x in zip(weights,nodes):
                if x/p0 < 1:
                    # 小截止处用已证单调性，以两个点值包住整区间。
                    y=x/p0
                    points=[direct(w,z,top,p0,p1,x_over_p0=v) for v in (y.lower(),y.upper())]
                    unit=points[0].union(points[1])
                else:
                    unit=high(x,p0,p1)
                total+=c*unit
            return total,nodes
        left,low_nodes=weak(s_lo)
        right,high_nodes=weak(s_hi)
        qmin = mp*low_nodes[0].union(high_nodes[0]).sqrt()
        qmax = low_nodes[-1].union(high_nodes[-1]).sqrt()
        assert sum(v < p0 for v in high_nodes) == (3 if q==-4 else 1)
        assert sum(v < p0 for v in low_nodes) == (3 if q==-4 else 1)
        print(q,"κ、t、h/p0 左右端、D 左右端",box[0].str(13),box[1].str(13),s_lo,s_hi,left.str(12),right.str(12),flush=True)
        assert left>0>right
        print("最低截止 GeV、最高截止与普朗克之比",qmin.str(12),qmax.str(12),flush=True)


if __name__ == "__main__":
    certify(Path(__file__).resolve().parent)
