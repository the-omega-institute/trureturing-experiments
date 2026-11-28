"""五背景预登记的几何排序及加入第六实际背景时的交汇检验。"""
from pathlib import Path
from runpy import run_path

from flint import arb

ROOT = Path(__file__).resolve().parent
PARTNERS = (323,285,18,15)


def load(name):
    return run_path(str(ROOT/f"ch39-{name}.py"))


def ranking(scores,fifth):
    chosen = (*PARTNERS,fifth)
    assert all(scores[i]>scores[j] for i,j in zip(chosen,chosen[1:]))
    assert all(scores[fifth]>value for i,value in enumerate(scores) if i not in chosen)
    gap = scores[PARTNERS[-1]]-scores[fifth]
    print("前五名、分数、第四至第五名的分差",chosen,[scores[i].str(12) for i in chosen],gap.str(12))
    return gap


def verify():
    paired = load("shared-twelve-node-gauge-gap")["paired_lines"]
    shell = load("native-rh-ray-pvm")["shell_weights"]
    edges = [paired(ROOT,shell,2306888100085+j) for j in range(6)]
    assert all(len(row)==361 for row in edges)
    scores = [[abs(w)*(b-a) for a,b,w in row] for row in edges]
    training = [sum(row)/5 for row in zip(*scores[:5],strict=True)]
    assert ranking(training,17)>arb("49.90")
    assert ranking(scores[5],300)>arb("77.93")

    _,lib,(background,) = load("two-zero-cutoff-gravity")["gravity_fields"](
        ROOT,(2345676,),tight_neutral=True)
    model,_,_,initial,_ = lib["chiral"]["source"]["mean"]["inputs"]()
    beta = [b/(2*arb.pi()) for b in model["beta"]]
    fields = [(0,0,0),*load("signed-gravity-unification")["gauge_weights"](model)]
    kernel = load("continuous-split-roots")
    sectors = [(edges[5],kernel["series"](background,fields),background)]
    direction = (-77073971362,-100000000000,34478833972,-100000000000,0)
    projection = load("three-background-calibration")["projection"]
    bound,_ = projection(sectors,fields,beta,initial,model,kernel,[direction],None,set(),
                         scale_millionths=[6000])
    norm = sum(abs(arb(n)/10**11) for n in direction[2:4])
    print("第六组达到引力双目标后，所有比较能标的最大规范差绝对值下界",(bound/norm).str(15))
    assert bound/norm > arb("0.33504")


if __name__=="__main__":
    verify()
