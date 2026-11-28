"""六质量与三个混合角对壳实源 CP 相位的区间检验。"""
from itertools import product
from pathlib import Path
from runpy import run_path
import sys

import numpy as np
from flint import arb


def cosine_bound(values,inverse,spectrum):
    up,down=values[:3],values[3:6]
    ss,tt,ww=[v.sin() for v in values[6:9]]
    a,b,z=ss**2,tt**2,ww**2
    # CKM 元素模方在 cos(delta)=0 处的常数项。
    constant = [[values[6].cos()**2*values[8].cos()**2,a*values[8].cos()**2,z],
                [a*(1-b)+(1-a)*b*z,(1-a)*(1-b)+a*b*z,b*(1-z)],
                [a*b+(1-a)*(1-b)*z,(1-a)*b+a*(1-b)*z,(1-b)*(1-z)]]
    amplitude=2*ss*tt*ww*values[6].cos()*values[7].cos()
    exponent=arb(2)/3
    if inverse:
        exponent=-exponent
    u=[(v**2/(up[(i+1)%3]*up[(i+2)%3]))**exponent for i,v in enumerate(up)]
    d=[(v**2/(down[(i+1)%3]*down[(i+2)%3]))**(-exponent) for i,v in enumerate(down)]
    value=sum(u[i]*d[j]*constant[i][j] for i,j in product(range(3),repeat=2))
    slope=(u[1]-u[2])*(d[0]-d[1])*amplitude
    assert not slope.contains(0)
    return (spectrum-value)/slope


def verify(replay, signs):
    library=run_path(str(replay/"ch39-shell-real-ten-observation-window.py"))
    orbit=run_path(str(replay/"ch39-two-zero-orbit-condition.py"))
    model,up,down=orbit["calibration"]()
    shell=run_path(str(replay/"ch39-shell-pencil-six-mass-window.py"))
    chart,_=shell["shell_frame"](model)
    adj=model["adj"]
    relative=up*down.inv()
    eigen,basis=model["eig"](adj(relative)*relative)
    eigen[0]=abs(relative.det())**2/(eigen[1]*eigen[2])
    assert 0 < eigen[0] < eigen[1] < eigen[2]
    normalized=relative/abs(relative.det())**(arb(1)/3)
    spectra=[sum(abs(v)**2 for v in matrix.entries()) for matrix in (normalized,normalized.inv())]
    normals=[library["axes"](adj(chart)*matrix*basis) for matrix in (down.inv(),adj(down))]
    standard=model["load"]("shell-clock-ckm-blocks")["standard_matrix"]
    source_values=[v+arb(0,t.upper()) for v,t in zip(model["ref"][:9],model["tol"][:9],strict=True)]
    cosine=cosine_bound(source_values,False,spectra[0])
    assert arb(".23937") < cosine < arb(".63364")
    assert arb(".85") < cosine.upper().acos() and cosine.lower().acos() < arb("1.35")
    phase_count=128
    table=library["phase_table"](phase_count)
    priority=[1,10,5,1100,550,560,260,430,990,1010]
    for sign in signs:
        center=[*model["ref"][:9],sign*arb("1.10")]
        width=[*model["tol"][:9],arb(".25")]
        stack=[([0.]*10,[1.]*10,0,np.arange(phase_count**2))]
        closed=trace_closed=deepest=0
        while stack:
            point,radius,depth,cells=stack.pop()
            values=[v+arb(c)*w+arb(0,(arb(r)*w).upper())
                    for c,r,v,w in zip(point,radius,center,width,strict=True)]
            observed=values[9].cos()
            bounds=[cosine_bound(values,inv,s) for inv,s in enumerate(spectra)]
            if not bounds[0].overlaps(bounds[1]) or any(not x.overlaps(observed) for x in bounds):
                trace_closed+=1
                continue
            cosine=bounds[0].intersection(bounds[1]).intersection(observed)
            values[9]=(sign*cosine.acos()).intersection(values[9])
            rows=library["coefficients"](values,eigen,standard,normals,adj)
            cells=library["surviving_cells"](rows,table,cells)
            if 0<len(cells)<=32 and depth>=6:
                cells=np.array([i for i in cells if not library["certify_phase_box"](
                    rows,7,int(i)//phase_count,int(i)%phase_count)],dtype=np.int64)
            if not len(cells):
                closed+=1
                deepest=max(deepest,depth)
                continue
            assert depth<36,(point, radius)
            def split_priority(i):
                return priority[i]*radius[i]*(.25/.023 if i==9 else 1)
            j=max(range(10),key=split_priority)
            smaller=radius.copy()
            smaller[j]/=2
            for direction in (-1,1):
                shifted=point.copy()
                shifted[j]+=direction*smaller[j]
                stack.append((shifted,smaller,depth+1,cells))
        print("CP 分支、壳实性排除盒、谱迹排除盒、实性排除最大深度",sign,closed,trace_closed,deepest)


if __name__=="__main__":
    verify(Path(sys.argv[1]).resolve() if len(sys.argv)>1 else Path(__file__).resolve().parent,(-1,))
