"""Size of the singlet dimension-five correction and the single-Higgs running of the minor."""
from pathlib import Path
from runpy import run_path

from flint import arb, ctx

smdr = run_path(str(Path(__file__).with_name("ch39-shell-gw-smdr-quartic.py")))
q = arb


def planck_term():
    d1 = q(".002421")-q(".0000805")
    d2 = q(".002421")
    u = (1-q(".597"))*(1-q(".02409"))
    bound = q(".120")
    lo, hi = q(0), q(".030")
    for _ in range(40):
        mid = (lo+hi)/2
        total = mid+(mid**2+d1).sqrt()+(mid**2+d2).sqrt()
        if total < bound:
            lo = mid
        else:
            assert total > bound
            hi = mid
    residual = u/hi-(1-u)/d1.sqrt()
    minimum_c = residual*d2.sqrt()*q("2.435e18")/q("4.473e13")
    print("mass-sum benchmark, IO highest possible m3 and minimum inverse residual",
          hi.str(8), residual.lower().str(8), "dimension-five coefficient exceeds", minimum_c.lower().str(8))
    assert minimum_c > 20000


def running():
    ctx.cap = 5
    beta = smdr["evaluator"](5)
    y = smdr["init"](".967", ".55853")
    h = q(1)/1536
    loop = smdr["loop"]
    first = q("7.41e10").union(q("7.66e10"))
    min_t, max_t = ((a/q("91.1876")).log()/loop for a in (first.lower(),first.upper()))
    n = int(float((max_t/h).upper()))+1
    total_min = total_max = q(0)
    for step in range(n):
        duration = arb.max(q(0), arb.min(h, max_t-step*h))
        guaranteed = arb.max(q(0), arb.min(h, min_t-step*h))
        # Validate one initial broad tube; two nested interval substitutions
        # reduce its width for the two small charged lepton Yukawa couplings.
        tube = [v+arb(0,".022") for v in y]
        for refinement in range(3):
            speed = beta(tube)
            assert all(abs(h*r)<q(".022") for r in speed)
            tube = [v+q(0).union(h)*r for v,r in zip(y,speed)]
        gap = tube[11]**2-tube[12]**2
        assert gap>0
        total_min += 2*guaranteed*gap.lower()**2/loop
        total_max += 2*duration*gap.upper()**2/loop
        y = smdr["advance"](beta,y,h)
    print("single-light two-loop log drift", total_min.lower().str(9), total_max.upper().str(9))
    assert total_min > q("1e-16") and total_max < q("5e-16")


if __name__ == "__main__":
    planck_term()
    running()
