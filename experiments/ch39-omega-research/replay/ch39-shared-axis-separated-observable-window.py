"""三壳逐点可行尺度处，共同物理转角对独立观测量拟合的分离证书。"""
from pathlib import Path
from runpy import run_path
from flint import arb,ctx

def verify(root,lib):
    interval=lib["interval"]
    slopes=(interval(".0181999",".0182001"),interval(".1129995",".1130005"))
    profiles={}
    def collect(vals,polar):
        tu,td,r,p,l=vals
        # 归一化正二阶块的 Bloch 半径至多为一。
        if r>1:
            return True
        if not (r>arb(".9958") and polar<arb("-.9745")):return False
        if not (tu.rad()<arb(".000003") and td.rad()<arb(".0000005") and r.rad()<arb(".0004") and p.rad()<arb(".0004") and l.rad()<arb(".0004")):return False
        vals=(tu,td,r.intersection(interval(".9958","1")),p,l)
        key=int(float(tu.mid())*100000)
        if key in profiles:vals=[a.union(b) for a,b in zip(profiles[key],vals,strict=True)]
        profiles[key]=vals
        return True
    lib["certify_invariants"](root,slopes,collect)
    fixed=run_path(str(root/"ch39-free-four-site-native-positive-spectrum.py"))
    center=2_306_888_100_085
    up,down=fixed["native"](root,center-1)
    left=lib["native_rays"](up);support=lib["forms"](down)
    total=0;least=arb(100)
    for key,v in sorted(profiles.items()):
        for j,visited,margin in lib["exclude_rays"](v,left,support,allowed=(16,18)):
            total+=visited
            if margin is not None:least=min(least,margin)
    dist=2*((2*(1-arb(".9958"))).sqrt()+(2*(1-arb(".9745"))).sqrt())
    print("不变量包围组数、方位角盒、严格边际、共同轴所需最大弦长",len(profiles),total,least.str(12),dist.str(12))
    for node in (center,center+1):
        rays=lib["native_rays"](fixed["native"](root,node)[0])
        shortest=min(sum((a-b)**2 for a,b in zip(x,y,strict=True)).sqrt().lower() for j,_,x in left if j in (16,18) for _,_,y in rays)
        assert shortest>dist
        print("实际编号、左邻保留射线到本壳任一射线的最小弦长",node,shortest.str(12))

if __name__=="__main__":
    root=Path(__file__).resolve().parent
    lib=run_path(str(root/"ch39-left-invariant-character-window.py"))
    verify(root,lib)
    witness=run_path(str(root/"ch39-common-embedding-central-right-ckm.py"))
    n=witness["CENTRAL"]
    geom=(".0182",".113",(".15674",".06997",".97496","-.14136"))
    for node,phase,down in ((n,"0",(1,6,10,16)),(n+1,".029266",(6,7,10,12))):
        witness["verify_nodes"](root,((node,phase),),".74",geom,((6,7,10,13),down))
    left=(".0182",".113",("-.21668","-.53516",".71837","-.38807"))
    witness["verify_nodes"](root,((n-1,"0"),),".74",left,((1,10,16,18),(6,7,10,18)))
