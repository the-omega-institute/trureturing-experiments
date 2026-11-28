"""荷五第 163 枚零点的左相消边界由双站配双站的判别式折点给出。"""
from itertools import combinations, combinations_with_replacement
from pathlib import Path
from runpy import run_path
import sys

from flint import acb,arb,ctx


def area(z,w):
    return (z.conjugate()*w).imag


def endpoint_entries(matrix):
    z,w = matrix[6,0],matrix[6,5]
    return z,w,matrix[2,0]-z,matrix[2,5]-w


def coefficients(matrix):
    z,w,u,v = endpoint_entries(matrix)
    return area(z,w),area(u,w)+area(z,v),area(u,v)


def coefficient_slopes(matrix,slope):
    z,w,u,v = endpoint_entries(matrix)
    dz,dw,du,dv = endpoint_entries(slope)
    return (area(dz,w)+area(z,dw), area(du,w)+area(u,dw)+area(dz,v)+area(z,dv),
            area(du,v)+area(u,dv))


def discriminant(matrix):
    a,b,c = coefficients(matrix)
    return b*b-4*a*c


def discriminant_slope(matrix,slope):
    a,b,c=coefficients(matrix)
    da,db,dc=coefficient_slopes(matrix,slope)
    return 2*b*db-4*(da*c+a*dc)


def intervals(left,right,count):
    for i in range(count):
        yield (left+(right-left)*i/count).union(left+(right-left)*(i+1)/count)


def verify(root):
    shell=run_path(str(root/"ch39-three-ray-chiral-boundary.py"))
    slope=run_path(str(root/"ch39-proton-endpoint-speed.py"))["green_and_slope"]
    supports=run_path(str(root/"ch39-proton-next-79-heights.py"))
    ctx.prec=192
    sites=(0,2,3,5,6,9,11,13)
    height=acb.zeta_zero(163).imag

    def at(offset):
        return shell["laplacian"](5,height+offset).inv()

    left,right=arb("-.0798856832281540"),arb("-.0798856832281535")
    assert discriminant(at(left))<0<discriminant(at(right))
    fold=left.union(right)
    green,dh=slope(shell,height+fold)
    a,b,c=coefficients(green)
    weight=-b/(2*c)
    assert arb(".0224648659957")<weight<arb(".0224648659959")
    dh_disc=discriminant_slope(green,dh)
    dp_disc=discriminant_slope(green,-green*green)
    assert dh_disc>arb(".000000025")
    rate=-dp_disc/dh_disc
    assert arb("-.018691929096")<rate<arb("-.018691929092")

    # 排除折点邻域内的另 405 个 2×2 支撑及全部单行、单列。
    neighborhood=at(fold+arb(0,"1e-8"))
    for i in sites:
        assert supports["separated"]([neighborhood[i,j] for j in sites])
    excluded=0
    for (i,k),(j,m) in combinations_with_replacement(combinations(sites,2),2):
        if {(i,k),(j,m)}=={(2,6),(0,5)}:
            continue
        supports["exclude_square"](neighborhood[i,j],neighborhood[i,m],neighborhood[k,j],neighborhood[k,m])
        excluded+=1
    assert excluded==405

    # 自折点至单站三角形出现处，判别式与固定站 6 的有向面积皆递增。
    join_start,join_end=arb("-.079886"),arb("-.079848")
    for offset in intervals(join_start,join_end,64):
        g,h=slope(shell,height+offset)
        a,b,c=coefficients(g)
        da,_,_=coefficient_slopes(g,h)
        assert da>0 and discriminant_slope(g,h)>0
        assert c<0 and 0 < -b/(2*c) < arb(".03")
        for x in intervals(arb(0),arb(".03"),16):
            z,w,u,v=endpoint_entries(g)
            assert ((z+x*u).conjugate()*(w+x*v)).real<0

    # 对所有接收站核实单站发射的最早可行时刻；随后同一三角形通至原零点。
    row_left,row_right=arb("-.07984834698844"),arb("-.07984834698841")
    assert coefficients(at(row_left))[0]<0<coefficients(at(row_right))[0]
    for offset in intervals(left,row_right,64):
        g=at(offset)
        base=g[6,0]
        assert all(area(g[6,j],base)>0 for j in sites if j not in (0,5))
    for offset in intervals(row_left,arb(0),2600):
        g,h=slope(shell,height+offset)
        z,w,v=[g[6,j] for j in (0,5,3)]
        assert area(h[6,0],w)+area(z,h[6,5])>0
        assert area(w,v)>0 and area(v,z)>0
    print("第 163 枚左边界、站 2 的临界发射权、对 p 的导数",fold.str(16),weight.str(15),rate.str(14))
    print("站 6 单独发射重新可行的边界、两端调权产生的外移",(row_left.union(row_right)).str(15),
          (row_left.union(row_right)-fold).str(12))


if __name__=="__main__":
    verify(Path(sys.argv[1]).resolve() if len(sys.argv)>1 else Path(__file__).resolve().parent)
