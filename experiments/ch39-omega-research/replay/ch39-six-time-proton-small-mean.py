"""七枚实际三窗零点的自由概率保持六项时间，配固定真空作用得到很小的十二面均方。"""
from pathlib import Path
from runpy import run_path
import sys

from flint import acb, arb, arb_mat, ctx


SUPPORT=(15447127,142269606336,142300003154,295280111868,296004209875,13217728165506,143921420592856)


def six_time_support(root, support):
    def read(name):
        return run_path(str(root / ("ch39-" + name + ".py")))

    clock=read('arithmetic-height-bandwidth')
    bands=read('six-time-height-bands')
    dictionary=read('seventh-time-zero-dictionary')
    physical=read('first-three-mode-observables')
    ctx.prec=176
    heights={n:acb.zeta_zero(n).imag for n in support}
    ctx.prec=512
    coordinates={n:1/(1+(h/bands['scale'])**2) for n,h in heights.items()}
    times=clock['times']()[:7]
    times[6]+=arb('1.360900010e-8')
    denominator,moments=bands['radial']['compact_moments'](times)
    _,prob=dictionary['certified_rule'](support,coordinates,times,moments,denominator)
    print('七枚实际零点概率',*(p.str(15) for p in prob),flush=True)
    ctx.prec=224
    model,initial=physical['model_inputs']()
    chiral=physical['joint']['chiral']
    for h in heights.values():
        rows = chiral['source']['spectral_supplier'](model)(h)
        layers = {q: chiral['ordered']([v for v, charge in rows if charge == q])
                  for q in (0, *model['charges'])}
        _, alpha, curve = physical['layer_observables'](model, initial, layers)
        assert arb('.117')<alpha<arb('.119') and arb('.998')<curve<arb('1.002')
    return heights, prob


def verify(root):
    def read(name):
        return run_path(str(root / ("ch39-" + name + ".py")))

    shell=read('three-ray-chiral-boundary')
    lib=read('mincut-global-winding')
    fixed=read('height-independent-flavor-vacuum')
    angles=read('mincut-reference-ckm-action')
    flavor=read('chiral-neutral-flavor-phase')
    global_model=read('reference-flavor-global-winding')
    heights, prob = six_time_support(root, SUPPORT)
    x=arb('-.01046')
    y=x-angles['reference_angles'](flavor)[0]
    tri=lib['terms'](shell)
    cap=arb(23)/25*arb(2)**-20/14**2
    delta=arb(0).union(cap)
    phase={n:shell['geo']['phases'](h) for n,h in heights.items()}
    gaps={n:(p[724]/p[804]).arg() for n,p in phase.items()}
    original=read('proton-factor-zero')
    frozen,solved=original['FIXED'],original['SOLVED']
    bare={q:lib['coeff'](tri,q,arb(0))[2]/12 for q in (*frozen,*solved)}
    def factor(q,t):
        return (bare[q]*acb(0,q*t).exp()).real
    cols={q:[arb(1),*fixed['gradient_columns'](lib,tri,q,x,y),factor(q,x-y)] for q in bare}
    weights={q:arb(v) for q,v in frozen.items()}
    matrix=arb_mat([cols[q] for q in solved]).transpose()
    assert not matrix.det().contains(0)
    rhs=arb_mat([[int(j==0)-sum(w*cols[q][j] for q,w in weights.items())] for j in range(len(solved))])
    weights.update(zip(solved,matrix.solve(rhs).entries(),strict=True))
    assert min(weights.values())>0
    modes=[(q,[w*c for c in lib['coeff'](tri,q,delta)]) for q,w in weights.items()]
    global_model['conjugation_symmetry'](root)
    drift=global_model['global_cover'](lib,tri,weights,modes,x,y)
    assert drift<0
    factors=[sum(w*factor(q,gaps[n]) for q,w in weights.items()) for n in SUPPORT]
    mean=sum(p*z**2 for p,z in zip(prob,factors,strict=True))
    assert arb('3.3452621807725e-6') < mean < arb('3.3452621807727e-6')
    print('七个十二面因子、六时间实际均方',*[v.str(14) for v in (*factors,mean)],flush=True)


if __name__=='__main__':
    verify(Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).resolve().parent)
