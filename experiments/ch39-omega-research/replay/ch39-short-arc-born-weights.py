"""两步开路和三步闭路的相干圆弧输运对逆平方谱的最大调制。"""
from collections import defaultdict
from pathlib import Path
from runpy import run_path

from flint import acb_mat, acb, arb, ctx

root = Path(__file__).parent
physical = run_path(str(root/"ch39-first-three-mode-observables.py"))
geometry = run_path(str(root/"ch39-arc-relative-determinant.py"))
indices = (8359, 12336, 19872, 22027, 28325, 48067, 58214, 2529893, 2562426,
           4035321, 4193973, 6736068, 6964994, 25215595, 25408663)
charges = (1, 2, 3, 4, 6)


def words(max_steps, closed):
    adj = defaultdict(list)
    for u, v, length, _ in geometry["edges"]:
        a, b = int(length == 724), int(length == 804)
        adj[u].append((v, a, b))
        adj[v].append((u, -a, -b))
    result = set()
    for start in range(19):
        layer = {(start, 0, 0)}
        for _ in range(max_steps):
            layer |= {(v, a+r, b+s) for u, a, b in layer for v, r, s in adj[u]}
        result.update((a,b) for v,a,b in layer if not closed or v == start)
    return sorted(result)


def records():
    data = []
    for index in indices:
        _, coupling, curve = physical["observables"](index)
        assert arb(".117") < coupling < arb(".119") and arb(".998") < curve < arb("1.002")
        ctx.prec = 256
        height = acb.zeta_zero(index).imag
        clocks = geometry["phases"](height)
        assert height < 10**8
        data.append((height, clocks[724], clocks[804]))
    return data


def gram_bounds(data, modes, charge):
    first = data[0][0]
    matrix = acb_mat([[(z**(charge*a))*(w**(charge*b)) for a,b in modes] for h,z,w in data])
    for j in range(1, len(data)):
        weight = first/(arb(".06")*data[j][0])
        for k in range(len(modes)):
            matrix[j,k] *= weight
    inverse = (matrix.conjugate().transpose()*matrix).inv()
    return 2*sum(abs(inverse[j,k]) for j in range(len(modes)) for k in range(len(modes)))


def verify():
    modes = words(2, False)
    assert modes == words(3, True) and len(modes) == 7
    assert (arb("2899267716.58")/10**8)**2 > 100
    # 从非负五项矩核扣去首枚后，对 x>=100 的其余支撑作总量估计。
    floor = run_path(str(root/"ch39-neutral-five-scale-floor.py"))
    ctx.prec, ctx.cap = 768, 3
    first, second = [acb.zeta_zero(i).imag for i in (8359, 12336)]
    bound = arb(117387924448)
    cap = ((2903616618/second)**2).upper()
    _, _, _, contact, den, *_ = floor["envelope"](cap)
    times = floor["model"]["shared"]["clock"]["clock"]["times"]()
    moments = floor["neutral"]["contrast"]["moments_five"](times)
    expected = sum(a*v for a, v in zip(contact.coeffs(), moments, strict=True))
    p = arb("4.65271e-11")
    remaining = (expected-p*contact(bound)/den(bound))*den(arb(100))/contact(arb(100))
    assert remaining < arb(".06")**2*p

    data = records()
    ctx.prec = 384
    caps = (370416000, 6202100, 594637, 116420, 12124)
    for charge, cap in zip(charges, caps, strict=True):
        limit = gram_bounds(data, modes, charge)
        assert limit < cap
        print("荷、全相位最大效率与首枚效率之比上界", charge, limit.str(13))


if __name__ == "__main__":
    verify()
