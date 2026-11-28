"""第 120 枚两端和第 163 枚右端的单站临界电流及局部动量导数。"""
from pathlib import Path
from runpy import run_path
import sys

from flint import acb,arb,ctx

SITES=(0,2,3,5,6,9,11,13)
# 编号、零移位端点、两活动站、内侧第三站、内侧方向、集中度与七行的分离泛函。
CASES=(
    (120,"-.086956806487102",0,5,3,1,".88","-.0098483475630",
     ((-87,206),(-263,625),(-66,1365),(-400,6270),(-116,4136),(-38,361),(-155,-15))),
    (120,".064357764508482",2,3,5,-1,".88",".0167062395219",
     ((-285,-708),(-148,-351),(-289,-5885),(-103,106),(-153,4),(-17,-114),(-28,-3721))),
    (163,".070074949093800",2,3,5,-1,".89",".0082888591060",
     ((-251,-596),(-70,-166),(-20339,-359680),(-56,-12),(-154,20),(-49,-513),(-123,-4059))),
)


def area(z,w):
    return (z.conjugate()*w).imag


def verify(root):
    shell=run_path(str(root/"ch39-three-ray-chiral-boundary.py"))
    green_slope=run_path(str(root/"ch39-proton-endpoint-speed.py"))["green_and_slope"]
    ctx.prec=192
    for index,center,r,k,third,sign,minimum,expected,normals in CASES:
        h=acb.zeta_zero(index).imag
        t=arb(center)
        lo,hi=t-arb("1e-14"),t+arb("1e-14")
        for offset,direction in ((lo,-1),(hi,1)):
            g=shell["laplacian"](5,h+offset).inv()
            assert direction*area(g[6,r],g[6,k])>0
        g,dh=green_slope(shell,h+lo.union(hi))
        z,w=g[6,r],g[6,k]
        assert (z.conjugate()*w).real<0
        def area_slope(derivative):
            return area(derivative[6,r],w)+area(z,derivative[6,k])
        assert area_slope(dh)>0
        rate=-area_slope(-g*g)/area_slope(dh)
        assert abs(rate-arb(expected))<arb("2e-12")

        outside=t-sign*arb("1e-7")
        start,stop=(outside,arb(0)) if sign==1 else (arb(0),outside)
        for i in range(2880):
            cell=(start+(stop-start)*i/2880).union(start+(stop-start)*(i+1)/2880)
            g,dh=green_slope(shell,h+cell)
            z,w,v=g[6,r],g[6,k],g[6,third]
            assert area(dh[6,r],w)+area(z,dh[6,k])>0
            assert sign*area(w,v)>0 and sign*area(v,z)>0

        g=shell["laplacian"](5,h+t+arb(0,"1e-7")).inv()
        f=arb(minimum)
        good=min(g[i,j].real.lower() for i in SITES for j in SITES if {i,j} not in ({6,3},{6,5}))
        bad=min(g[i,j].real.lower() for i in SITES for j in SITES if {i,j} in ({6,3},{6,5}))
        c=good/(good-bad)
        assert c>arb(1)/2 and (c+(2*c-1).sqrt())/2>f
        ref=g[6,r]
        assert abs(ref)>0
        normal=-sign*ref.conjugate()/abs(ref)
        assert all(sign*area(g[6,j],ref)>0 for j in SITES if j not in (r,k))
        for i,(x,y) in zip((i for i in SITES if i!=6),normals,strict=True):
            for m in SITES:
                for j in SITES:
                    v=f*g[6,j]+(1-f)*g[m,j]
                    assert (normal*g[i,j]).imag-(x*v.real+y*v.imag)/1000>arb(".000005")
        print(index,"左端" if sign==1 else "右端",lo.union(hi).str(16),"动量导数",rate.str(15),
              "发射站 6，接收站",r,k)


if __name__=="__main__":
    verify(Path(sys.argv[1]).resolve() if len(sys.argv)>1 else Path(__file__).resolve().parent)
