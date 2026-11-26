"""两个实际壳谱在有理手征顶角下的费米子一圈有效势，减去同谱轻模型。"""
from pathlib import Path
from runpy import run_path

from flint import acb_mat, arb, ctx


resolver = run_path(str(Path(__file__).with_name("ch39-chiral-mass-resolvent-yukawa.py")))
filtered = resolver["filtered"]
model = filtered["model"]
adjoint = model["adjoint"]


def light_directions(mass):
    values, right = (adjoint(mass)*mass).eig(right=True,algorithm="rump")
    assert all(v.real>0 and not any(v.overlaps(w) for w in values[:j]) for j,v in enumerate(values))
    left = acb_mat(3,3)
    masses = [v.real.sqrt() for v in values]
    for j,m in enumerate(masses):
        length = sum(abs(right[i,j])**2 for i in range(3)).sqrt()
        for i in range(3):
            right[i,j] /= length
        for i in range(3):
            left[i,j] = sum(mass[i,k]*right[k,j] for k in range(3))/m
    return left, right, masses


def trace_and_log(parts, matrices, background, scale, power=1, response=None):
    moment, log_moment, uniform = arb(0), arb(0), arb(0)
    for (left,right), y, colors in zip(filtered["sectors"],matrices,(3,3,1)):
        ml, mr = parts[left][1], parts[right][1]
        s = background**2
        fl, fr = [[(s/(s+m*m))**power for m in vals] for vals in (ml,mr)]
        factor = background/arb(2).sqrt()
        mass = acb_mat([[factor*y[i,j]*(fl[i]*fr[j] if response is None else response(ml[i],mr[j]))
                         for j in range(19)] for i in range(19)])
        light = acb_mat([[mass[i,j] for j in range(3)] for i in range(3)])
        u,v,ms = light_directions(light)
        lh = adjoint(u)*acb_mat([[mass[i,j] for j in range(3,19)] for i in range(3)])
        hl = acb_mat([[mass[i,j] for j in range(3)] for i in range(3,19)])*v
        hh = acb_mat([[mass[i,j] for j in range(3,19)] for i in range(3,19)])
        ranges = [(ms,mr[3:],lh),(ml[3:],ms,hl),(ml[3:],mr[3:],hh)]
        contribution, perturb2 = arb(0), arb(0)
        for rows,cols,block in ranges:
            for i,a in enumerate(rows):
                for j,b in enumerate(cols):
                    e2 = abs(block[i,j])**2
                    x,z = a*a,b*b
                    assert not x.overlaps(z)
                    kernel = 2*(x*x*((x/scale**2).log()-1)-z*z*((z/scale**2).log()-1))/(x-z)
                    contribution += e2*kernel
                    perturb2 += e2
                    moment += 2*colors*(x+z)*e2
        # 轻与含重块分开平方，避免百万亿 GeV 阈值或共同轻块的数值消去。
        ll = acb_mat([[ms[i] if i==j else 0 for j in range(3)] for i in range(3)])
        norm2 = resolver["norm_squared"]
        excess4 = (2*norm2(ll*adjoint(hl)) + norm2(adjoint(hl)*hl)
                   + 2*norm2(adjoint(ll)*lh+adjoint(hl)*hh)
                   + norm2(adjoint(lh)*lh+adjoint(hh)*hh))
        moment += colors*excess4
        radius = max(m.upper() for m in (*ml,*mr,*ms)) + perturb2.sqrt()
        assert radius>scale
        third_derivative = 24*(scale+radius*(radius**2/scale**2).log())+16*radius
        remainder = third_derivative*(2*perturb2)**arb("1.5")/12
        assert remainder < arb("1e-36")
        log_moment += colors*(contribution+arb(0,remainder.upper()))
        abs_log = max(abs((radius/mu).log()*2).upper() for mu in (arb("91.1876"),arb("1.22089e19")))
        second_derivative = radius**2*(12*(1+abs_log)+4)
        uniform += colors*second_derivative*perturb2/2
    return moment,log_moment,uniform


def verify():
    ctx.prec = 768
    products = filtered["products"]
    q = arb(1000)
    v = arb("246.22")
    loop = 16*arb.pi()**2
    roots_and_moments = []
    for n,w,z in ((22027,model["WITNESS"],model["HIGGS"]),
                  (28325,products["SECOND_WITNESS"],products["SECOND_HIGGS"])):
        parts, matrices = filtered["source"](n,model["column"](w),model["column"](z))
        for power in (1,2):
            moment,logs,uniform = trace_and_log(parts, matrices, v, q, power)
            zero_scale = q*(logs/(2*moment)).exp()
            planck = arb("1.22089e19")
            reference_density = arb("1e-122")*planck**4/(8*arb.pi())
            band = [(-logs + 2*moment*(mu/q).log())/loop for mu in (arb("91.1876"),planck)]
            assert moment > 0
            if power == 1:
                assert -logs/loop < arb("-1.93e-10")
                roots_and_moments.append((zero_scale,moment))
            else:
                assert all(abs(b) < arb("2.49e-9")*reference_density for b in band)
                assert uniform/loop < arb(".000234")*reference_density
                print("全电弱背景及减法尺度上的势差绝对上界/目标密度",
                      (uniform/loop/reference_density).str(10))
            print(n,power,"四次迹差、1 TeV 有效势差 GeV^4、势差过零尺度 GeV",
                  moment.str(12),(-logs/loop).str(12),zero_scale.str(12))
            print("从 Z 到普朗克尺度的增量与原目标密度之比",*( (b/reference_density).str(9) for b in band))

    (q1,t1),(q2,t2) = roots_and_moments
    best = 2/loop*t1*t2/(t1+t2)*abs((q2/q1).log())
    print("一次响应在任意共同减法尺度上的最小最大残差 GeV^4",best.str(12))
    assert best > arb("1e-14")


if __name__ == "__main__":
    verify()
