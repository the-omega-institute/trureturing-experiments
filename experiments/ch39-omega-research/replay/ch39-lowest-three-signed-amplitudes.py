"""仅传播最低三层、先压缩实振幅再平方的质量逆谱与顶底混合。"""
from itertools import product
from pathlib import Path
from runpy import run_path

import numpy as np
from flint import arb, acb_mat, arb_mat, ctx


def adj(a):
    return a.conjugate().transpose()


def column(values):
    return arb_mat([[x] for x in values])


def inverse_branches(rays, low, middle):
    v = acb_mat([[u[i,0] for u in rays] for i in range(3)])
    gram = adj(v)*v
    k,l,d = [1-abs(gram[i,j])**2 for i,j in ((0,1),(0,2),(1,2))]
    g=gram.det().real
    assert g>arb(".21") and min(k,l,d)>0
    assert 0<low<middle<1
    p=[u*adj(u) for u in rays]
    def vals(a,b,c):
        return column((a+b+c,k*a*b+l*a*c+d*b*c,g*a*b*c))
    def jac(a,b,c):
        return arb_mat([[1,1,1],[k*b+l*c,k*a+d*c,l*a+d*b],
                        [g*b*c,g*a*c,g*a*b]])
    models=[]
    for light_sign,middle_sign in product((1,-1),repeat=2):
        x,y=light_sign*low,middle_sign*middle
        x0,y0=x.mid(),y.mid()
        target=column((1+x0+y0,x0+y0+x0*y0,x0*y0))
        a1,a2,a3=[float(z) for z in target.entries()]
        triples=[]
        # 近似对称的 Gram 数据仅用于生成牛顿初值，后续方程用实际三个二阶主子式。
        kk=float((k+l)/2)
        for a in sorted(np.roots([kk,-kk*a1,a2,-float(d/g)*a3])):
            assert a.imag==0
            for b in sorted(np.roots([1,a-a1,a3/(float(g)*a)])):
                assert b.imag==0
                triples.append(list(map(arb,(float(a),float(b),float(a1-a-b)))))
        boxes=[]
        for point in triples:
            for _ in range(3):
                step=jac(*point).solve(vals(*point)-target)
                point=[(v-e.mid()).mid() for v,e in zip(point,step.entries(),strict=True)]
            pre=jac(*point).inv()
            residual=pre*(vals(*point)-target)
            dx,dy=x-x0,y-y0
            # 将三个相关的质量多项式先乘逆雅可比，再独立加入两种质量的偏差。
            residual-=pre*column((1,1+y0,y0))*dx+pre*column((1,1+x0,x0))*dy
            residual-=pre*column((0,1,1))*(dx*dy)
            radii=[(arb("1.6")*abs(z).upper()+arb("1e-18")).upper() for z in residual.entries()]
            box=[v+arb(0,r.upper()) for v,r in zip(point,radii,strict=True)]
            defect=arb_mat([[int(i==j) for j in range(3)] for i in range(3)])-pre*jac(*box)
            for j in range(3):
                assert abs(residual[j,0])+sum(abs(defect[j,i])*radii[i] for i in range(3))<radii[j],(
                    light_sign,middle_sign,box,j,residual[j,0],radii)
            assert all(any(not s.overlaps(t) for s,t in zip(box,b,strict=True)) for b in boxes)
            boxes.append(box)
            ref=sum((w*P for w,P in zip(point,p,strict=True)),acb_mat(3,3))
            levels,vecs=ref.eig(right=True,algorithm="rump")
            idx=max(range(3),key=lambda i:levels[i].real.mid())
            eig=levels[idx]
            assert eig.imag.contains(0) and eig.real>arb(".99")
            u=acb_mat([[vecs[j,idx]] for j in range(3)])
            u=u*(1/sum(abs(z)**2 for z in u.entries()).sqrt())
            delta=abs(1-eig)+sum(r*abs((adj(z)*u)[0,0]) for r,z in zip(radii,rays,strict=True))
            radius=(delta/(1-abs(y).upper())).asin().upper()
            models.append((u,radius))
    return models


def verify(root):
    model=run_path(str(root/"ch39-three-source-positive-ckm.py"))
    up,down=model["source_rays"](root)
    windows=[((".00000601",".00000629"),(".00306",".00316"),(".8573",".8659")),
             ((".0000133",".0000137"),(".000265",".000271"),(".01392",".01410"))]
    branches=[]
    for rays,rows in zip((up,down),windows,strict=True):
        small,mid,large=[arb(a).union(arb(b)) for a,b in rows]
        models=inverse_branches(rays[:3],small/large,mid/large)
        assert len(models)==24
        print("最低三层，每部门分支数及重方向最大跟踪角误差",
              len(models),max(r for _,r in models).str(12),flush=True)
        branches.append(models)
    tiny,large=[],[]
    for (u,ru),(v,rv) in product(*branches):
        angle=model["distance"](u,v)
        if angle+ru+rv < arb(".0419483"):
            tiny.append((angle+ru+rv).upper())
        elif angle-ru-rv > arb(".0434"):
            large.append((angle-ru-rv).lower())
        else:
            raise AssertionError((angle,ru,rv))
    assert (len(tiny),len(large))==(192,384)
    assert max(tiny).cos()>arb(".9992027") and min(large).cos()<arb(".5902170")
    print("两类组合数、近角上界、远角下界、两端顶底模",len(tiny),len(large),
          max(tiny).str(12),min(large).str(12),max(tiny).cos().str(12),min(large).cos().str(12))


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)
