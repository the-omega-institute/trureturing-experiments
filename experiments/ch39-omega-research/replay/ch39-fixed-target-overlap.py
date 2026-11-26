"""Neighbor radius required by target-overlap monotonicity on the nonorthogonal probes."""
from flint import acb, arb, ctx

FIRST, LAST, TARGET = 10855, 12092, 10875
GUARANTEED = set(range(10905,12055)) | {10871,TARGET}


def verify():
    ctx.prec = 160
    phases = {}
    for n,z in enumerate(acb.zeta_zeros(FIRST,LAST-FIRST+1),FIRST):
        g = z.imag
        phases[n] = (acb(g,arb(1)/2)/acb(g,-arb(1)/2), arb(2).sqrt()*acb(0,g*arb(2).log()).exp())
    wt,dt = phases[TARGET]
    fidelity, score = {}, {}
    for n,(w,d) in phases.items():
        fidelity[n] = abs((1+wt.conjugate()*w)/2)**4 * abs((1+dt.conjugate()*d)/3)**6
        score[n] = abs(1-w)**2*abs((d-1)*(d-2)*(9+4*d))/13
    trapped = 11391
    rivals = [fidelity[n] for n in score if n not in (TARGET,trapped) and score[n]>score[trapped]]
    assert all(score[n]<score[trapped] or score[n]>score[trapped] for n in score if n!=trapped)
    assert fidelity[trapped] > arb(".99999837") and fidelity[trapped] > max(rivals)
    assert trapped-TARGET == 516
    for n in score:
        if n == TARGET:
            continue
        assert any(score[j] > score[n] and fidelity[j] > fidelity[n]
                   for j in GUARANTEED if abs(j-n)<=516)
    print("Source-target fidelity", fidelity[trapped].str(15), "best other uphill destination", max(rivals).str(15))
    print("Smallest radius allowing a simultaneous increase of both readings from every other state: 516.")


if __name__ == "__main__":
    verify()
