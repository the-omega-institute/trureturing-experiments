"""任意非负三辐边刚度、双时钟相位和激活尺度的矢量阈值增量界。"""
from pathlib import Path
from runpy import run_path

from flint import acb_mat, arb, ctx

star = run_path(str(Path(__file__).with_name("ch39-neutral-clock-star-threshold.py")))
geo, pole, orbit = star["geo"], star["pole"], star["orbit"]


def positive(matrix):
    for n in range(1, matrix.nrows()+1):
        minor = acb_mat([[matrix[i, j] for j in range(n)] for i in range(n)])
        assert minor.det().real > 0


def verify():
    ctx.prec = 224
    edges = geo["edges"]
    U = acb_mat([[int(i == x)-int(i == pole) for x in orbit] for i in range(19)])
    neutral = acb_mat([[arb(1)/19]*19 for _ in range(19)])
    lower = acb_mat(19, 19)
    clock_bound = acb_mat(19, 19)
    for i, j, d, z in edges:
        if pole in (i, j) and (i in orbit or j in orbit):
            continue
        neutral[i, i] += 1
        neutral[j, j] += 1
        neutral[i, j] -= 1
        neutral[j, i] -= 1
        if d in (724, 804):
            clock_bound[i, i] += 2
            clock_bound[j, j] += 2
        else:
            lower[i, i] += 1
            lower[j, j] += 1
            lower[i, j] = -z**5
            lower[j, i] = lower[i, j].conjugate()
    positive(neutral)
    positive(lower)

    def gram(matrix):
        return U.transpose()*matrix.solve(U)

    k0, kmax, kmin = gram(neutral), gram(lower), gram(lower+clock_bound)
    positive(kmin-arb(".805")*k0)
    positive(arb("1.383")*k0-kmax)
    coeff = 7/(8*arb.pi())
    eye = acb_mat([[int(i == j) for j in range(3)] for i in range(3)])
    delete = coeff*max(((eye+k).det().real.log().upper() for k in (kmax, k0)))
    assert delete < arb(".16260")
    gram_det = (U.transpose()*U).det().real
    low = coeff*max(((36**3*k.det().real/gram_det).log().upper() for k in (kmax, k0)))
    high = 3*coeff*(max(-arb(".805").log().lower(), arb("1.383").log().upper())+arb(2).log())
    print("Deletion, low-cutoff and high-cutoff total upper bounds", delete.str(12),
          (delete+low).str(12), (delete+high).str(12))
    assert delete+low < arb("1.435") and delete+high < arb("1.014")

    radius = (geo["theta"](804)/geo["theta"](12))*arb(".999").union(arb("1.001"))
    assert radius > 6
    # With no spokes all levels are active here. Neutral's artificial unit
    # eigenvalue gives the positive pseudodeterminant, without adding exposure.
    bare_upper = coeff*(2*radius.log()-(lower.det().real/neutral.det().real).log())
    bare_lower = coeff*(2*radius.log()-((lower+clock_bound).det().real/neutral.det().real).log())
    tail = 3*coeff*(1+36/radius**2).log()
    total_upper = bare_upper-3*coeff*arb(".805").log()+tail
    total_lower = bare_lower-3*coeff*arb("1.383").log()-tail
    print("Absolute weak minus color vector threshold at geometric cutoff",
          total_lower.lower().str(12), total_upper.upper().str(12))
    assert total_lower > arb("-.15") and total_upper < arb("2.31")


if __name__ == "__main__":
    verify()
