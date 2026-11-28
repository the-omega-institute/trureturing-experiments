"""由目标零点的荷六折叠谱计数，在预登记的十二处预测高能交汇。"""
from itertools import combinations
from pathlib import Path
from runpy import run_path

from flint import acb,arb,ctx

ROOT = Path(__file__).resolve().parent
TRAIN = (2345712,2346088,2346090,2346094)
ALLOWED = ((1,3),(4,6),(5,6),(7,8))
TARGETS = tuple(range(2346124,2346147,2))
BANDS = (0,1,13,16)


def load(name):
    return run_path(str(ROOT/f"ch39-{name}.py"))


def levels(shell,heavy,charge):
    height = acb.zeta_zero(heavy).imag
    delta = arb(23)/25*arb(2)**-20/(height**2+shell["soft"]["point"]["study"]["offset"])
    roots = shell["laplacian"](charge,height).eig(multiple=True)
    assert len(roots)==19 and all(v.imag.contains(0) for v in roots)
    vals = [v.real for v in roots]
    folded = [v*abs(v-15+delta/3) for v in vals]
    norm = []
    for seq in (vals,folded):
        seq.sort(key=lambda v:float(v.mid()))
        assert all(a<b for a,b in zip(seq,seq[1:])) and seq[0]>0
        norm.append([v/seq[-1] for v in seq[:-1]]+[arb(1)])
    return norm


def prefixes(seq,d):
    result = []
    for i in range(d+1):
        cutoff = arb(i)/d
        assert all(v<=cutoff or v>cutoff for v in seq)
        result.append(sum(v<=cutoff for v in seq))
    return result


def count(prefix,bands):
    a,b,c,d = bands
    return prefix[b]-prefix[a]+prefix[d]-prefix[c]


def verify_selection(shell):
    data = {(charge,fold):[] for charge in (1,2,3,4,6) for fold in (0,1)}
    for charge in (1,2,3,4,6):
        for heavy in TRAIN:
            for fold,seq in enumerate(levels(shell,heavy,charge)):
                data[charge,fold].append(seq)
    for d in (4,8,16,32):
        valid = []
        for key,seqs in data.items():
            grids = [prefixes(v,d) for v in seqs]
            for bands in combinations(range(d+1),4):
                if all(lo<=count(pr,bands)<=hi for pr,(lo,hi) in zip(grids,ALLOWED,strict=True)):
                    valid.append((*key,bands))
        print("训练分母与满足全部计数的配方数",d,len(valid))
        if d<32:
            assert not valid
        else:
            assert len(valid)==33 and min(valid)==(6,1,BANDS)
            assert [count(prefixes(seq,d),BANDS) for seq in data[6,1]]==[3,4,5,7]
            return valid


def verify():
    shell = load("three-ray-chiral-boundary")
    ctx.prec = 320
    verify_selection(shell)
    ks = {n:count(prefixes(levels(shell,n,6)[1],32),BANDS) for n in TARGETS}
    print("目标谱给出的预测计数",ks,flush=True)
    geometry = load("three-share-energy-interval")
    entropy = load("one-of-nineteen-forward")
    locate = load("three-share-forward-feasibility")["locate"]
    shared = load("two-share-global-forward")["context"](TARGETS,20)
    *_,baseline,increments,loss = shared
    small = (baseline+sum(increments[i] for i in geometry["LINES"]))/20+loss
    assert small<1
    model,_,_,_,_ = load("chiral-minimal-neutral-gravity")["chiral"]["source"]["mean"]["inputs"]()
    beta = [v/(2*arb.pi()) for v in model["beta"]]
    low,high = arb(10)**13,arb(10)**17
    selected,available = [],[]
    for heavy,k in ks.items():
        point,mz = geometry["setup"](heavy,batch_context=shared)
        roots,boxes = entropy["search"](point,(arb(k)/(19-k)).log()) if 0<k<19 else ([],0)
        print("目标、计数、内点数、覆盖盒",heavy,k,len(roots),boxes,flush=True)
        for s in roots:
            energy,common = point(s,common_readout=True)
            ell = (energy/mz).log()
            couplings = [common,*[common+b*ell for b in beta]]
            checks = [ell,(arb("1.22089e19")/energy).log(),*couplings]
            assert all(not v.contains(0) for v in [*checks,energy-low,high-energy])
            physical = all(v>0 for v in checks)
            in_window = physical and low<energy<high
            if in_window:
                selected.append(heavy)
            print("驻点能标、物理、高能",heavy,energy.str(10),physical,in_window,flush=True)
        witness,boxes = locate(point,geometry["EXCLUDED"],low,high)
        if witness is not None:
            available.append(heavy)
        print("放开计数高能域、覆盖或见证前盒数",heavy,witness is not None,boxes,flush=True)
    assert list(ks.values()) == [7,2,7,4,7,5,6,6,7,5,7,2]
    assert not selected and available == [2346136,2346144]
    print("选中、全部高能、漏选",selected,available,sorted(set(available)-set(selected)),flush=True)


if __name__=="__main__":
    verify()
