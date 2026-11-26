"""Calibrate the finite positive-weight witnesses at 96 and 272 binary places."""
from pathlib import Path
from runpy import run_path
import json
from flint import acb, arb, arb_mat, ctx


def integer_floor(value):
    k=value.floor().unique_fmpz()
    assert k is not None
    return int(k)


def data(config):
    root=Path(__file__).resolve().parent
    words=run_path(str(root/"ch39-shortest-triangle-clock-words.py"))
    rates=run_path(str(root/"ch39-two-prime-height-law.py"))["rate_numerators"]
    ctx.prec,ctx.cap=1792,3
    levels=run_path(str(root/"ch39-rh-shell-first-conflict.py"))["read"](54)
    den=1<<288
    result={"depth":288,"levels":[integer_floor(x*den) for x in levels],"cases":{}}
    for tick,src in config.items():
        tick=int(tick)
        ids=src["support"]; n=len(ids)
        offsets=src["indices"]
        assert len(offsets)+4==n
        heights=[acb.zeta_zero(i).imag for i in ids]
        bins=[integer_floor(v*den) for v in heights]
        c=rates(tick)
        phases=[[] for _ in c+c]
        for h in bins:
            u,v=[-43235114615559591407865398948523*((h*s+2**311)//2**312) for s in c]
            for row,k in zip(phases,(u,v,v-u,v+u)):
                lo,hi=words["power"](k)[0]
                row.append((arb(lo)/words["DEN"]).union(arb(hi)/words["DEN"]))
        rows=[[1/(1+(j/g)**2) for g in heights] for j in offsets]+phases
        mid=arb_mat([[t.mid() for t in row] for row in rows])
        inv=mid.inv(); error=inv*arb_mat(rows)-arb_mat([[int(i==j) for i in range(n)] for j in range(n)])
        norm=max(sum(abs(error[i,j]).upper() for j in range(n)).upper() for i in range(n))
        assert norm<arb(2)**-300
        trial=(inv*arb_mat([[levels[j]] for j in offsets]+[[0]]*4)).entries()
        radius=norm*max(abs(v).upper() for v in trial)/(1-norm)
        # The integer votes themselves are subsequently checked. Choosing
        # midpoints here does not require proving the exact weights round uniquely.
        nums=[integer_floor(v.mid()*den) for v in trial[:-1]]
        nums.append(den-sum(nums))
        assert min(nums)>0 and min(x.lower() for x in trial)-radius>0
        result["cases"][str(tick)]={"indices":offsets,"rows":[list(row) for row in zip(ids,bins,nums)]}
    return result


if __name__=="__main__":
    root=Path(__file__).resolve().parent
    original=json.loads((root/"ch39-minimum-tick-finite-inputs.json").read_text())
    config={k:{"support":[r[0] for r in v["rows"]],"indices":v["indices"]} for k,v in original["cases"].items()}
    assert data(config)==original
    print("Zeros, 54 target cells, two positive-weight witnesses and tick calibrations certified.")

