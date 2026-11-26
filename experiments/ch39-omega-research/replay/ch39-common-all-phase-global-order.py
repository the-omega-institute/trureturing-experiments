"""同一相位与同一作用在所有高度给出合窗参考味全局真空的全域最小阶。"""
import json
from decimal import Decimal
from pathlib import Path
from runpy import run_path
import sys

from flint import acb, arb, ctx


def rotation(angle):
    return acb(0,angle).exp()


def separation(modes, sign, distance, record):
    _, left, right, directions, rivals = record
    gx, rx, gy, ry = map(arb,directions)
    points = [(arb(x),arb(y),arb(w)) for x,y,w in rivals]
    assert all(w>0 for _,_,w in points)
    total = sum(w for _,_,w in points)
    constants = []
    oscillatory = []
    for q,(bare,response) in enumerate(modes,1):
        # S(θ)-S(z) 的常数基线相消。四个有理乘子乘未缩放的原始梯度。
        reference = sum(w*sum((v*rotation(q*t)).real for v,t in zip(bare,(x,y,x-y),strict=True))
                        for x,y,w in points)
        united = [(-total-acb(0,q*g))*a-acb(0,q*r)*b
                  for g,r,a,b in zip((gx,gy,gx-gy),(rx,ry,rx-ry),bare,response,strict=True)]
        offset = rotation(q*sign*distance)
        constants.append(reference+(united[2]*offset).real)
        oscillatory.append(united[0]+united[1]*offset.conjugate())
    # x=s t, y=s(t-D_s)。四分量分离余量沿等味线仅余单个谐波。
    pending = [(arb(left),arb(right))]
    count = 0
    worst = arb(100)
    while pending:
        l,r = pending.pop()
        t = l.union(r)
        margins = [v+(c*rotation(q*sign*t)).real for q,(v,c) in enumerate(zip(constants,oscillatory),1)]
        low = min(z.lower() for z in margins)
        if low>arb(".1"):
            count += 1
            worst = worst.min(low)
        else:
            assert r-l>arb("1e-10"), (sign,t,low,margins.index(min(margins,key=float)))
            mid=(l+r)/2
            pending.extend(((l,mid),(mid,r)))
    return count,worst


def verify(root):
    cutoff=105
    shell = run_path(str(root/"ch39-three-ray-chiral-boundary.py"))
    lib = run_path(str(root/"ch39-mincut-global-winding.py"))
    previous = run_path(str(root/"ch39-monotone-winding-action.py"))
    corners = run_path(str(root/"ch39-curvature-window-density-cap.py"))
    ctx.prec=192
    # 荷二折叠给两条窄带；两钟同时加 π 的提升由荷一至少四个轻模排除。
    corners["restrict_corners"]()
    bands = [(arb(l),arb(r)) for l,r in previous["BANDS"]]
    distances = dict(zip((1,-1),previous["verify_phases"](root,shell,bands),strict=True))
    tri=lib["terms"](shell)
    modes=[]
    keys=((1,0),(0,1),(1,-1))
    for q in range(1,cutoff+1):
        bare=lib["coeff"](tri,q,arb(0))
        response=[sum((z**q for z,spoke in tri[k] if spoke),acb(0)) for k in keys]
        modes.append((bare,response))
    file=Path(__file__).with_name("ch39-common-all-phase-global-certs.json")
    rows=json.loads(file.read_text(),parse_float=str)
    used=0
    for sign in (1,-1):
        for l,r in previous["BANDS"]:
            # 用精确十进有理数核对所用闭段无缝覆盖整个已证外区间。
            parts=[row for row in rows if row[0]==sign and Decimal(l)<=Decimal(row[1])<Decimal(r)]
            assert Decimal(parts[0][1])==Decimal(l) and Decimal(parts[-1][2])==Decimal(r)
            assert all(Decimal(x[2])==Decimal(y[1]) for x,y in zip(parts,parts[1:]))
            count=0
            worst=arb(100)
            for row in parts:
                assert Decimal(row[1])<Decimal(row[2])
                limit=104 if sign==-1 and Decimal(l)<Decimal(".075") else 105
                n,z=separation(modes[:limit],sign,distances[sign],row)
                count+=n
                worst=worst.min(z)
            used+=len(parts)
            print("等味支、窄带、排除至绕数、证书与末段数、严格正余量",sign,l,r,limit,len(parts),count,worst.str(12))
    assert used==len(rows)


if __name__ == "__main__":
    root=Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).resolve().parent
    verify(root)
