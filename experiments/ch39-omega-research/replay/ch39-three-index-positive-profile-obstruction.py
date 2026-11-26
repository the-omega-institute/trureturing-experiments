"""三枚相邻实际零点的破缺引力差矩对任何共同非负累计截止轮廓的分离证书。"""
from pathlib import Path
from runpy import run_path

from flint import arb


def coefficients(fields,d):
    # 四个坐标依次为 2π(ΔK′,ΔW′/c)，两个相邻高度分别减中心高度。
    c=arb("1e-10")
    pairs=((-d[0]-d[2],-d[1]-d[3]),(d[0],d[1]),(d[2],d[3]))
    return [(height,i,z,k*a,-w*g/(2*c)) for height,(field,(k,w)) in enumerate(zip(fields,pairs))
            for i,(g,a,z) in enumerate(field[4:])]


def verify(root):
    api=run_path(str(root/"ch39-common-cutoff-negative-mass.py"))
    _,fields=run_path(str(root/"ch39-three-index-free-atoms.py"))["spectra"](root)
    d=list(map(arb,("1","1",".47",".495")))
    radius=arb(".06")
    directions=[d]
    for j in range(4):
        for sign in (-1,1):
            directions.append([v+(sign*radius if i==j else 0) for i,v in enumerate(d)])
    # 最轻的两个行是同一编号中仅相差一个门户移质的标量与费米伴侣。
    plus=fields[2][4:]
    pair=(82,83)
    light,partner=pair if plus[82][2]<plus[83][2] else tuple(reversed(pair))
    assert plus[light][2]<plus[partner][2]
    zstar=plus[light][2]
    edge=arb("1e-14")
    cut=arb("1e-5")
    families=[]
    tails=[]
    for direction in directions:
        original=coefficients(fields,direction)
        tail=arb(0)
        good=arb(0)
        rows=[]
        for background,i,z,k,w in original:
            dz=arb(0) if (background,i)==(2,light) else z-zstar
            rows.append((dz,k,w))
            if background==2 and i in pair:
                assert k+arb(0).union(edge)*w>0
                if i==light:good=(k+arb(0).union(edge)*w).lower()
            else:
                assert dz>arb("1.75e-13")
                tail+=(abs(k)+edge*abs(w))*(-dz/edge).exp()
        tails.append((good-tail).lower())
        families.append((rows,api["long_series"](rows,27)))
    assert min(tails)>arb(".27333328")
    # 共同正因子 exp(−z*/x) 已移除；余下函数在极小 x 端也具有可检验的尺度。
    counts=[]
    for label,(rows,large) in enumerate(families):
        count=0
        for left,right,fn in ((edge,cut,lambda box:api["small_derivative"](box,rows)),(cut,arb(1),large)):
            stack=[(left,right)]
            while stack:
                lo,hi=stack.pop()
                value=fn(lo.union(hi))
                if value>0:
                    count+=1
                else:
                    assert hi-lo>arb("1e-25"),(label,lo,hi,value)
                    mid=((lo+hi)/2).mid()
                    stack.extend(((lo,mid),(mid,hi)))
        counts.append(count)
    print("最轻编号加一的质量平方比、归一分离量尾部统一下界",zstar.str(20),min(tails).str(16))
    print("中心与八个扰动方向的正区间数、四行归一化孔径",counts,radius.str(4))


if __name__=="__main__":
    verify(Path(__file__).resolve().parent)
