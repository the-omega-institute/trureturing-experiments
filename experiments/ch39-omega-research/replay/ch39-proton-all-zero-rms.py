"""固定十四项真空权时，96 个活跃面任意加权的自然计数均方最小值。"""

from pathlib import Path
from runpy import run_path
import sys

from flint import arb, arb_mat, ctx


FREE_SUPPORT = (0,1,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27)
CANCEL_SUPPORT = (0,1,6,7,8,10,11,12,13,14,15,16,17,18,19,20,22,24,25,26,27,31)


def minimum(gram, constraints, support, forced):
    n = len(gram)
    k = len(constraints[0])
    zero = arb(0)
    rhs = [arb(i==0) for i in range(k)]
    system = [([gram[i][j] for j in support]+[-z for z in constraints[i]]) for i in support]
    system += [[constraints[j][h] for j in support]+[zero]*k for h in range(k)]
    sol = arb_mat(system).solve(arb_mat([[z] for z in [zero]*len(support)+rhs])).entries()
    masses = [zero]*n
    for i,z in zip(support,sol):
        assert z > arb(".006")
        masses[i] = z
    multipliers = sol[len(support):]
    exclusions = []
    for i in range(n):
        dual = sum(a*b for a,b in zip(gram[i],masses))-sum(a*b for a,b in zip(constraints[i],multipliers))
        if i in support:
            assert dual.contains(0)
        elif i not in forced:
            exclusions.append(dual)
    margin = min(r.lower() for r in exclusions)
    assert margin > arb("3.4e-8" if not forced else "5.8e-6")
    for h,target in enumerate(rhs):
        assert (sum(p*c[h] for p,c in zip(masses,constraints))-target).contains(0)
    value = sum(x*sum(a*y for a,y in zip(row,masses)) for x,row in zip(masses,gram))
    assert (value-multipliers[0]).contains(0)
    print("强制零组、支撑数、非活跃组对偶余量下界、均方、均方根",forced,len(support),
          margin.str(12),value.str(22),value.sqrt().str(22))
    for i,p in enumerate(masses):
        print(i,p.str(16),end="  " if (i+1)%4 else "\n")
    return value


def verify(root):
    load = lambda name: run_path(str(root/("ch39-"+name+".py")))
    types,rows,delta = load("fixed-word-proton-faces")["signals"](root)
    ctx.prec = 256
    shell = load("three-ray-chiral-boundary")
    model = load("mincut-global-winding")
    faces = load("mincut-real-triangle-tomography")
    local = load("mincut-reference-ckm-action")
    flavor = load("chiral-neutral-flavor-phase")
    x = arb("-.01046")
    y = x-local["reference_angles"](flavor)[0]
    weights = load("height-independent-flavor-vacuum")["minimum_order_weights"](model,model["terms"](shell),x,y)
    assert len(types)==32 and len(weights)==14 and (sum(weights.values())-1).contains(0)
    labels = [faces["unit"](pair) for _,pair,_ in types]
    gram = [[arb(0) for _ in types] for _ in types]
    for i,(clock,_,_) in enumerate(types):
        for j in range(i+1):
            if clock == types[j][0]:
                phase = labels[i]*labels[j].conjugate()
                entry = sum(w*w*(phase**q).real for q,w in weights.items())/2
                gram[i][j] = gram[j][i] = entry
    # 严格正定同时保证最优组概率唯一。用 Sylvester 判据逐个认证顺序主子式。
    spectral_floor = arb("1e-8")
    shifted = [[v-spectral_floor if i==j else v for j,v in enumerate(row)]
               for i,row in enumerate(gram)]
    assert all(arb_mat([row[:k] for row in shifted[:k]]).det()>0 for k in range(1,33))
    plain = minimum(gram,[(arb(1),) for _ in types],FREE_SUPPORT,())
    forced = tuple(i for i,(key,row) in enumerate(zip(types,rows)) if key[2])
    assert forced==(2,21) and all(rows[i][1]<0 for i in forced)
    assert delta>0 and all(rows[i][1].is_zero() for i in range(32) if i not in forced)
    constraints = [(arb(1),ref,diff) for ref,_,diff in rows]
    cancel = minimum(gram,constraints,CANCEL_SUPPORT,forced)
    assert arb(".00013605795350898")<plain<arb(".00013605795350901")
    assert arb(".00018518829899377")<cancel<arb(".00018518829899382")
    print("相消约束引起的均方与均方根相对增幅", (cancel/plain-1).str(16),
          ((cancel/plain).sqrt()-1).str(16))
    return arb_mat(gram)


if __name__=="__main__":
    root = Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).resolve().parent
    verify(root)
