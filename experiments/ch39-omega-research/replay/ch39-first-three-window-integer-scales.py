"""第 8359 枚正零点的最小割、正整数质量标尺和参考味误差。"""
from pathlib import Path
from runpy import run_path
import sys

from flint import acb, arb, ctx


def verify(root):
    physical = run_path(str(root/'ch39-first-three-mode-observables.py'))
    scales = run_path(str(root/'ch39-joint-integer-scale-pairs.py'))
    shell = run_path(str(root/'ch39-three-ray-chiral-boundary.py'))
    cp = run_path(str(root/'ch39-three-window-two-probability-cp.py'))['cp']
    model, initial = physical['model_inputs']()
    ctx.prec, ctx.cap = 224, 3
    joint = physical['joint']
    n = 8359
    height = acb.zeta_zero(n).imag
    eig = model['shell']['shell_spectra']((0,*model['charges']),zero_index=n)
    center = sum(eig[0])/18
    layers = {q:joint['chiral']['ordered']([v*abs(v-center) for v in vals]) for q,vals in eig.items()}
    delta=arb(23)/25*arb(2)**-20/(height**2+shell['soft']['point']['study']['offset'])
    cut=9+(1-delta)**2
    third=max(layers[q][2].upper() for q in model['charges'])
    fourth=min(layers[q][3].lower() for q in model['charges'])
    print('高度、最高第三质量、最低第四质量',*(v.str(18) for v in (height,third,fourth)),flush=True)
    assert third<cut<fourth
    energy,alpha,curve=physical['layer_observables'](model,initial,layers)
    assert arb('.117')<alpha<arb('.119') and arb('.998')<curve<arb('1.002')
    unit=arb('2.435e18')/724**2
    a=(arb('246.22')**2/(arb('1.22089e19')*32*unit*layers[2][15]))**2
    c=(9*arb(3).log()-12*arb(2).log())/2
    main=c*a*a/(arb.pi()*arb('1e-122'))
    lo=32*(main*(1-arb('1e-57'))/arb('1.002')).root(4)
    hi=32*(main/arb('.998')).root(4)
    kmin=(lo*1000000).lower().floor()/1000000
    kmax=(hi*1000000).upper().ceil()/1000000
    print('必要有效倍率带、名义交点及两窗',*(z.str(15) for z in (kmin,kmax,energy,alpha,curve)),flush=True)
    flavor=cp['flavor']
    coordinate=flavor['singlet'](n)[1]
    reference=flavor['singlet'](11661)[1]
    data=cp['two']['spectral_data']()
    polys=cp['cp_data'](data,flavor['reference_matrix']())
    def read(v):
        return cp['two']['probabilities'](v,data)[0][0],cp['cp_value'](v,polys)
    relative=[u/v-1 for u,v in zip(read(coordinate),read(reference),strict=True)]
    print('六面坐标差与 U、J 有符号相对误差',*(v.str(20) for v in (coordinate-reference,*relative)),flush=True)
    assert abs(relative[0])>arb('2e-6') and abs(relative[1])>arb('1.2e-6')
    times={q:[(unit*v/model['mz']).log() for v in layers[q][3:]] for q in (1,3)}
    print('首个奇阈值',min(v.lower() for ts in times.values() for v in ts).str(15),flush=True)
    scales['verify_scale_layers'](model,initial,joint,layers,kmin,kmax)


if __name__=='__main__':
    verify(Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).resolve().parent)
