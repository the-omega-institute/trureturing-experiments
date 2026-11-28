"""Bound the derivative of the published fixed-clock two-loop crossing map."""
from pathlib import Path
from runpy import run_path

from flint import arb, arb_mat, ctx

ROOT = Path(__file__).parent
fixed = run_path(str(ROOT / "ch39-chiral-clock-unit-readout.py"))
ladder = run_path(str(ROOT / "ch39-portal-residual-resolution-ladder.py"))
two, sm = fixed["two"], fixed["sm"]


def bound(model, initial, read, index, *, wide):
    ctx.prec = 192
    levels = ladder["levels_of"](model, read, index)
    start = (arb("2.435e18")/724**2/model["mz"]).log()
    table = ladder["table_of"](model, levels, start)
    pi = arb.pi()
    limit = arb("29.040") + arb(0, ".014")
    strong = 1/arb(".1184") + arb(0, ".04")
    if wide:
        limit = arb("29.040") + arb(0, ".100")
        strong = (1/arb(".119")).union(1/arb(".117"))
    even = [*levels[0], *(h for q in (2,4,6) for h in levels[q][3:])]
    unit = arb("2.435e18")/724**2
    assert all(32*h*unit > model["mz"]*limit.exp() for h in even)
    definite = [r for r in table if r[0] < limit]
    possible = [r for r in table if not (r[0] > limit)]
    segments = fixed["gauge_segments"](limit.upper(), possible, 1024)
    baseline = [sm["BETA"][i]+sum(r[1][i] for r in definite) for i in range(3)]
    matrix = arb_mat([[baseline[i], 2*pi, -2*pi*int(i == 2)] for i in range(3)])
    transform = matrix.inv()
    schedule = ((b, B) for _, b, B, steps in segments for _ in range(steps))
    sensitivity = [arb(i == 2) for i in range(6)]
    yukawa = [[9,3,2],[3,9,2],[6,6,5]]
    lip = arb(0)

    def observe(t, step, gauge, y):
        nonlocal sensitivity, lip
        b, B = next(schedule)
        # X_i are inverse gauge couplings and h_i are logarithms of Yukawas.
        # Differentiate the six-dimensional flow with respect to X_3(0).
        J = arb_mat([
            [B[i][j]*gauge[j]**2/(8*pi**2) for j in range(3)]
            + [sm["D"][i][j]*y[j]/(16*pi**3) for j in range(3)]
            for i in range(3)
        ] + [
            [sm["C"][i][j]*gauge[j]**2/(4*pi) for j in range(3)]
            + [yukawa[i][j]*y[j]/(16*pi**2) for j in range(3)]
            for i in range(3)
        ])
        pad = step.upper()/4
        enclosure = [v + arb(0, pad) for v in sensitivity]
        ds = (J * arb_mat([[v] for v in enclosure])).entries()
        assert all(abs(step*v) < pad for v in ds)
        if not (t+step < limit.lower() or t > limit.upper()):
            flow = [-(sum(B[i][j]*gauge[j] for j in range(3)))/(8*pi**2)
                    + sum(sm["D"][i][j]*y[j] for j in range(3))/(32*pi**3)
                    for i in range(3)]
            # At a threshold use both adjacent smooth pieces. The one-loop
            # slope above the baseline is exactly the derivative of the ramps.
            slopes = arb_mat([[2*pi*flow[i]-(b[i]-baseline[i]),
                               2*pi*(enclosure[i]-int(i==2))] for i in range(3)])
            derivative = transform * slopes
            factor = arb(".04") / arb(".014")
            rows = [abs(derivative[0,0]) + factor*abs(derivative[0,1]),
                    abs(derivative[2,0])/factor + abs(derivative[2,1])]
            lip = lip.max(rows[0]).max(rows[1])
        sensitivity = [v+step*dv for v,dv in zip(sensitivity, ds)]

    two["flow"]([*initial[:2],strong],segments,sm["D"],sm["C"],record=observe)
    target = arb(".42") if wide else arb(".152")
    assert next(schedule,None) is None and 0 < lip < target, lip
    print(index, "thresholds", len(definite), len(possible), "Lipschitz bound", lip.upper().str(8), flush=True)
    return lip


if __name__ == "__main__":
    import sys
    model, _, _, initial, _ = fixed["chiral"]["source"]["mean"]["inputs"]()
    ctx.prec = 192
    read = fixed["chiral"]["source"]["spectral_supplier"](model)
    for index in (list(map(int,sys.argv[1:])) or [*ladder["FIT"], *ladder["hold"]["HOLD"]]):
        for wide in (False, True):
            bound(model, initial, read, index, wide=wide)
