"""中央三特征块的高零点反例，沿固定双钟角盒比较软极点退调。"""
from pathlib import Path
from runpy import run_path
import sys
from flint import acb,arb,ctx

def verify(root):
    ctx.prec=256
    shell=run_path(str(root/'ch39-three-ray-chiral-boundary.py'))
    pump=run_path(str(root/'ch39-native-light-projector-pump.py'))
    bases=list(pump['rotation_bases'](shell))
    adj=pump['adjoint']
    def lows(q,height,phases=None):
        lap=shell['laplacian'](q,height,phases=phases)
        out=[]
        for block in bases:
            vs=(adj(block)*lap*block).eig(algorithm="rump")
            assert all(v.is_finite() and v.imag.contains(0) for v in vs)
            vals=sorted((v.real for v in vs),key=lambda x:float(x.mid()))
            assert 0<vals[0]<vals[1]
            out.append(vals[0])
        return out
    for q,n in ((-4,2345666),(2,964890)):
        h=acb.zeta_zero(n).imag
        ctr=lows(q,h)
        unit=shell['geo']['phases'](h)
        near={d:z*acb(0,arb(0,arb('1e-27').upper())).exp() for d,z in unit.items()}
        lifts=[x-y for x,y in zip(lows(q,2*h,near),ctr,strict=True)]
        print(n,q,[v.str(13) for v in lifts],flush=True)
        assert all(v>arb('4.67e-22' if q==-4 else '5.46e-22') for v in lifts)
        offset = -39 if q==-4 else 14
        other = lows(q,acb.zeta_zero(n+offset).imag)
        gaps = [a-b for a,b in zip(ctr,other,strict=True)]
        assert all(g>arb('.025' if q==-4 else '.014') for g in gaps)
        print('下降见证编号',n+offset,[g.str(11) for g in gaps])
if __name__=='__main__':
    verify(Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).resolve().parent)
