"""原边拉普拉斯和边界三角给出的有限有理味作用与重态泄漏界。"""
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, ctx


local = run_path(str(Path(__file__).with_name("ch39-local-light-generator.py")))
mix, adjoint, shell = local["mix"], local["adjoint"], local["mix"]["shell"]
controls = ((1,2,0,"2982.27174342341187158521"),
            (0,2,0,"83.01907414532314321401"),
            (0,2,1,"184.67882174109678886202"),
            (0,2,0,"83.01907414532314321401"),
            (0,1,0,"22756.86729024121849952320"))


def boundary_sources():
    sites = shell["geo"]["sites"]
    r = [sites.index((y,z,x)) for x,y,z in sites]
    a = shell["boundary_sites"]()[0]
    ring = (a, r[a], r[r[a]])
    omega = acb(-1,arb(3).sqrt())/2
    return [[acb_mat([[omega**(k*((ring.index(i)-s) % 3))/3 if i in ring else 0]
                      for i in range(19)]) for k in range(3)] for s in range(2)]


def filtered_word(lap, center, sources):
    identity = local["identity"](19)
    t = arb(3)/4 * lap * (lap - center*identity)
    response = (identity + t*t).inv()
    # 三次翻倍将二次响应升至十六次，分母谱均落在 [1/2,1]。
    for _ in range(3):
        a = response*response
        rest = identity - response
        response = (a + rest*rest).solve(a)
    columns = [[response*v for v in rows] for rows in sources]
    gate = identity
    for i,j,site,t in controls:
        x,y = columns[site][i],columns[site][j]
        block = arb(t) * (x*adjoint(y) - y*adjoint(x))
        gate = gate*(identity-block).solve(identity+block)
    return gate


def frobenius(matrix):
    return sum(abs(z)**2 for z in matrix.entries()).sqrt()


def verify():
    ctx.prec = 448
    sources = boundary_sources()
    calibrated, _ = mix["calibration"](mix["boundary_frames"](mix["light_modes"](28325)))
    assert all((i,j,s)==(a,b,r) and abs(arb(t)-v)<arb("1e-20")
               for (i,j,s,t),(a,b,r,v) in zip(controls,calibrated,strict=True))
    error = max((arb(".161")*3/4)**16, (arb("10.95")*3/4)**-16)
    total = 8*error*sum(abs(arb(t)) for _,_,_,t in controls) + arb("1e-18")
    assert total < arb("4.88e-10")
    print("完整十九维酉乘积相对理想谱投影门的统一算子误差上界", total.str(11))
    for n in (11661,11934,12207,22027,28325,48067,58214):
        h = acb.zeta_zero(n).imag
        lap = shell["laplacian"](1,h)
        neutral = shell["laplacian"](0,h)
        center = sum(neutral[k,k].real for k in range(19))/18
        bands = sorted([abs(z.real*(z.real-center)) for z in lap.eig(multiple=True)],key=float)
        assert bands[2] < arb(".161") and bands[3] > arb("10.95")
        gate = filtered_word(lap,center,sources)
        frame = mix["light_modes"](n)
        for j in range(3):
            v = frame[shell["boundary_sites"]()[0],j]
            factor = v.conjugate()/abs(v)
            for i in range(19):
                frame[i,j] *= factor
        output = gate*frame
        light = adjoint(frame)*output
        heavy = output-frame*light
        frames = mix["boundary_frames"](mix["light_modes"](n))
        expected = mix["read"](frames,calibrated,1)
        assert frobenius(light-expected) < arb("1e-11")
        probability = frobenius(heavy)**2/3
        print(n, "均匀三轻态入重概率",probability.str(7))
        assert 0 < probability < arb("1.70e-30")


if __name__ == "__main__":
    verify()
