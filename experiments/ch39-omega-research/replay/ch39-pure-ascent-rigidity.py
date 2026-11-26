"""Interpolation witnesses for exact pure-state dynamics and the radius-58 identity obstruction."""
from flint import acb, acb_mat, arb, ctx

BASE, TOP, TARGET = 10855, 12092, 10875


def certification():
    ctx.prec = 192
    probes, scores = {}, {}
    for n,z in enumerate(acb.zeta_zeros(BASE,TOP-BASE+1),BASE):
        g = z.imag
        w = acb(g,arb(1)/2)/acb(g,-arb(1)/2)
        d = arb(2).sqrt()*acb(0,g*arb(2).log()).exp()
        probes[n] = w,d
        scores[n] = abs(1-w)**2*abs((d-1)*(d-2)*(9+4*d))
    wt, dt = probes[TARGET]
    f = {n: abs((1+wt.conjugate()*w)/2)**4 * abs((1+dt.conjugate()*d)/3)**6
         for n,(w,d) in probes.items()}
    assert all(n == TARGET or (abs(w-wt) > 0 and abs(d-dt) > 0) for n,(w,d) in probes.items())

    def evaluation(labels, a, b):
        rows = []
        for n in labels:
            w,d = probes[n]
            shifted = 100000*(w-acb(1,arb(".00009")))
            unit = d/arb(2).sqrt()
            rows.append([shifted**j*unit**k for j in range(a+1) for k in range(b+1)])
        return abs(acb_mat(rows).det())

    labels = [10905+33*i for i in range(35)]
    assert max(labels) <= 12054
    determinant = evaluation(labels, 4, 6)
    assert determinant > arb("8e-53")

    forced = []
    for n in range(10905,12055):
        if all(j == n or scores[j] < scores[n] or f[j] < f[n]
               for j in range(max(BASE,n-58),min(TOP,n+58)+1)):
            forced.append(n)
    assert len(forced) == 17
    independent = [TARGET]+forced[:11]
    assert all(f[n] > 0 for n in independent)
    span = evaluation(independent, 2, 3)
    assert span > 0
    print("Purity test on 35 guaranteed states, determinant modulus", determinant.str(12))
    print("17 additional forced fixed points, 12 spanning vectors, determinant modulus", span.str(12))


if __name__ == "__main__":
    certification()
