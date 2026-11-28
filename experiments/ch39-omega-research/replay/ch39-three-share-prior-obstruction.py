"""同一重谱批中三种高能正例需要的公共伯努利参考概率。"""
from pathlib import Path
from runpy import run_path

from flint import arb

ROOT = Path(__file__).resolve().parent
TARGETS = (2346088,2346090,2346094)


def verify():
    def load(name):
        return run_path(str(ROOT/f"ch39-{name}.py"))
    model = load("three-share-energy-interval")
    split = load("one-of-nineteen-forward")["search"]
    cover = load("three-share-global-energy")["cover"]
    shared = load("two-share-global-forward")["context"](TARGETS,20)
    *_,baseline,increments,loss = shared
    assert (baseline+sum(increments[i] for i in model["LINES"]))/20+loss < 1
    sm,sp = arb(".95"),arb("1.05")
    low,high = arb(10)**13,arb(10)**17
    info,_,_,_,_ = load("chiral-minimal-neutral-gravity")["chiral"]["source"]["mean"]["inputs"]()
    beta = [v/(2*arb.pi()) for v in info["beta"]]
    expected = ([4,5,6],[5,6],[7,8])
    for heavy,success in zip(TARGETS,expected,strict=True):
        point,mz = model["setup"](heavy,batch_context=shared)
        def no_physics(s):
            return point(s,physical_only=True) is model["EXCLUDED"]
        outside = cover(((arb(1)/20,sm),(sp,arb(1000)/7)),no_physics)
        allowed = []
        for k in range(1,19):
            z = (arb(k)/(19-k)).log()
            roots,_ = split(point,z)
            for root in roots:
                energy,common = point(root,common_readout=True)
                ell = (energy/mz).log()
                couplings = [common,*[common+b*ell for b in beta]]
                checks = [energy-low,high-energy,*couplings]
                assert all(not v.contains(0) for v in checks)
                if all(v>0 for v in checks):
                    allowed.append(k)
        assert allowed == success
        print("编号、标尺窗外排除盒数、十九态允许计数",heavy,outside,allowed,flush=True)

        if heavy==2346090:
            continue
        def prior_separated(s):
            res = point(s,physical_only=True)
            if res is model["EXCLUDED"]:
                return True
            if res is None:
                return False
            energy,shares,_ = res
            if energy<low or energy>high:
                return True
            if not all(0<p<1 for p in shares):
                return False
            _,_,_,slopes,_ = point(s,derivatives=True)
            total = sum(slopes)
            if not total>0:
                return False
            z = sum(d*(p/(1-p)).log() for p,d in zip(shares,slopes,strict=True))/total
            q = 1/(1+(-z).exp())
            if heavy==2346088:
                return q<arb(".323")
            return q>arb(".330")
        count = cover(((sm,sp),),prior_separated)
        print("全部高能物理点的参考概率严格界",heavy,count,"<0.323" if heavy==2346088 else ">0.330",flush=True)


if __name__=="__main__":
    verify()
