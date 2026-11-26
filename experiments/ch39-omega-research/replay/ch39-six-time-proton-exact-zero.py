"""七枚实际三窗正零点精确给出六项时间，共用作用权在指定编号上消去十二面因子。"""
from pathlib import Path
from runpy import run_path
import sys

from flint import acb, arb, arb_mat, ctx


SUPPORT=(28681727, 142245205585, 142285199398, 295860731185, 295891313832, 13219305981909, 143921323287490)
FIXED={1: '0.022924098228705', 3: '0.021081065426690', 35: '0.277801422813467'}
SOLVED=(104, 105, 83, 30, 67, 33, 69, 40, 5, 32, 81, 38)


def verify(root):
    def read(name):
        return run_path(str(root / ("ch39-" + name + ".py")))

    shell=read('three-ray-chiral-boundary')
    lib=read('mincut-global-winding')
    fixed=read('height-independent-flavor-vacuum')
    angles=read('mincut-reference-ckm-action')
    flavor=read('chiral-neutral-flavor-phase')
    global_model=read('reference-flavor-global-winding')
    heights, _ = read('six-time-proton-small-mean')['six_time_support'](root, SUPPORT)
    x=arb('-.01046')
    y=x-angles['reference_angles'](flavor)[0]
    tri=lib['terms'](shell)
    cap=arb(23)/25*arb(2)**-20/14**2
    delta=arb(0).union(cap)
    phase={n:shell['geo']['phases'](h) for n,h in heights.items()}
    gaps={n:(p[724]/p[804]).arg() for n,p in phase.items()}
    bare={q:lib['coeff'](tri,q,arb(0))[2]/12 for q in (*FIXED,*SOLVED)}
    def factor(q,t):
        return (bare[q]*acb(0,q*t).exp()).real
    cols={q:[arb(1),*fixed['gradient_columns'](lib,tri,q,x,y),*[factor(q,gaps[n]) for n in SUPPORT]] for q in bare}
    weights={q:arb(v) for q,v in FIXED.items()}
    matrix=arb_mat([cols[q] for q in SOLVED]).transpose()
    assert not matrix.det().contains(0)
    rhs=arb_mat([[int(j==0)-sum(w*cols[q][j] for q,w in weights.items())] for j in range(len(SOLVED))])
    weights.update(zip(SOLVED,matrix.solve(rhs).entries(),strict=True))
    assert min(weights.values())>0
    assert all(abs(gaps[n]-gaps[k])>arb('1e-5') and abs(gaps[n]+gaps[k])>arb('1e-5')
               for j,n in enumerate(SUPPORT) for k in SUPPORT[j+1:])
    print('七个实际角差与最小绕数权',*[gaps[n].str(13) for n in SUPPORT],min(weights.values()).str(13))
    modes=[(q,[w*c for c in lib['coeff'](tri,q,delta)]) for q,w in weights.items()]
    global_model['conjugation_symmetry'](root)
    drift=global_model['global_cover'](lib,tri,weights,modes,x,y)
    assert drift<0
    for n in SUPPORT:
        assert sum(w*factor(q,gaps[n]) for q,w in weights.items()).contains(0)
    print('七枚支撑上的十二面因子全部严格为零，六时间均方等于零。')


if __name__=='__main__':
    verify(Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).resolve().parent)
