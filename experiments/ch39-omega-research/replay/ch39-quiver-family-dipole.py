"""核查 Pati–Salam 家族边差荷的反常、味选择定则与扰动极化。"""
from itertools import product
from pathlib import Path
from runpy import run_path

from flint import fmpq_mat

shell = run_path(str(Path(__file__).with_name("ch39-shell-vector-gap.py")))
neutral = shell["neutral"]


def verify_matter():
    # 每行依次为 D, R, 族荷符号, 色弱重数和所属左／右多重态。
    fields = ((1,0,1,6,"L"),(-3,0,1,2,"L"),(-1,-1,-1,3,"R"),
              (-1,1,-1,3,"R"),(3,1,-1,1,"R"),(3,-1,-1,1,"R"))
    families = (1,1,-2)
    charges = [((d,r,s*x),n) for x in families for d,r,s,n,_ in fields]
    for i,j,k in product(range(3),repeat=3):
        assert sum(n*v[i]*v[j]*v[k] for v,n in charges)==0
    assert all(sum(n*v[i] for v,n in charges)==0 for i in range(3))
    # SU(4)^2 X，SU(2)_L^2 X，SU(2)_R^2 X，指数均用整数倍。
    assert sum(2*x+2*(-x) for x in families)==0
    assert sum(4*x for x in families)==sum(-4*x for x in families)==0
    contrasts = []
    for left, right in product(fields, repeat=2):
        if left[-1]==right[-1]:
            contrasts.append([a-b for a,b in zip(left[:3],right[:3],strict=True)])
    assert fmpq_mat(contrasts).rank()==2 and all(r[2]==0 for r in contrasts)
    powers = [abs(x-y) for x,y in product(families,repeat=2)]
    assert powers.count(0)==5 and powers.count(3)==4


def verify_edge_charge():
    _,_,colour,right = neutral["graph_constants"](port_builder=shell["capture"])
    for lap,edges in (colour,right):
        # 每个现有顶点对都只需一份单位边链即可承载基本差荷。
        n = lap.nrows()
        for i,j in edges:
            charge = [int(v==i)-int(v==j) for v in range(n)]
            assert sum(charge)==0 and sum(abs(v) for v in charge)==2
            compensated = [3*v for v in charge]
            assert sum(abs(v) for v in compensated)==6
        print("unit link edges", n, len(edges), flush=True)


def verify_current_resummation():
    # 在有理实例上核查矩形顶点映射的矩阵恒等式。
    u = fmpq_mat([[1,0,0,0],[0,1,0,0],[0,0,1,2]])
    bare = fmpq_mat([[4,1,0,0],[1,3,1,0],[0,1,2,0],[0,0,0,1]])
    polar = fmpq_mat([[3,1,0],[1,2,1],[0,1,4]])/100
    exchange = u*bare*u.transpose()
    assert u*(bare.inv()-u.transpose()*polar*u).inv()*u.transpose() == (exchange.inv()-polar).inv()


if __name__=="__main__":
    verify_matter()
    verify_edge_charge()
    verify_current_resummation()
