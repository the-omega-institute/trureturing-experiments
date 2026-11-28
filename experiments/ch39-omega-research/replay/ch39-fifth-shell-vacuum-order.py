"""固定壳分布且在源重谱满足牛顿目标时，认证目标重谱真空矩的变化。"""
from itertools import permutations
from math import factorial
from pathlib import Path
from runpy import run_path

from flint import arb, arb_mat

ROOT = Path(__file__).resolve().parent
HOLDOUT = tuple(range(2_345_700, 2_345_723, 2))


def vacuum_intervals(load, numbers, *, source_only=False, selected_pairs=None):
    context = load("shared-large-scale-gap")["load"]()
    fields, kernel = context[1], context[-2]
    shell = load("fourth-fifth-shell-exchange")["FIFTH_SHELL"]
    edges = load("shared-twelve-node-gauge-gap")["paired_lines"](
        ROOT, load("native-rh-ray-pvm")["shell_weights"], shell)
    scale = arb(7)/1000
    moments = load("fifth-shell-transfer")["shell_moments"](edges, kernel, scale)
    _, _, spectra = load("two-zero-cutoff-gravity")["gravity_fields"](ROOT, numbers, tight_neutral=True)
    floor = scale*min(a.lower() for a, _, _ in edges)
    assert all(max(z.upper() for _, _, z in rows) < min(arb(1), floor) for rows in spectra)
    inventory = [sum(g*z for g, _, z in rows) for rows in spectra]
    coupling = 2*sum(alpha for _, alpha, _ in spectra[0])
    assert coupling < 0
    series = [kernel["series"](rows, fields) for rows in spectra]
    remainder = 2*sum(abs(w) for _, _, w in edges)/(4*arb.pi()*factorial(kernel["ORDER"]+1))
    vacuum_tail = remainder*sum(abs(g) for g, _, _ in spectra[0])
    newton_tail = remainder*sum(abs(2*a) for _, a, _ in spectra[0])
    vacuum = arb_mat([[0,1,0,0,0]])
    newton = arb_mat([[1,0,0,0,0]])
    actual_moments = None
    if source_only:
        frozen = load("fifth-shell-frozen-witnesses")
        layout = frozen["data"]("ch39-fifth-shell-quadratic-layout.json")
        assert layout["heavy"] == numbers[0]
        actual_moments = frozen["witness_moments"](load, shell, edges, kernel, layout)
    pairs = ((0,j) for j in range(1,len(numbers))) if source_only else permutations(range(len(numbers)),2)
    if selected_pairs is not None:
        pairs = [(numbers.index(a), numbers.index(b)) for a,b in selected_pairs]
    for i,j in pairs:
        lead = (inventory[j]-inventory[i])/coupling
        coeff = vacuum*(series[j]-series[i])-lead*(newton*series[i])
        error = 2*vacuum_tail+abs(lead)*newton_tail
        change = lead+(coeff*moments)[0,0]+arb(0,error.upper())
        actual = None
        if actual_moments is not None:
            difference = vacuum*(series[j]-series[i])
            actual = (difference*actual_moments)[0,0]+arb(0,(2*vacuum_tail).upper())
            assert change.contains(actual)
        yield numbers[i], numbers[j], lead, change, actual


def verify():
    load = run_path(str(ROOT/"ch39-four-background-crossing.py"))["load"]
    numbers = load("fifth-shell-frozen-witnesses")["NUMBERS"]
    for a,b,lead,change,_ in vacuum_intervals(load, numbers):
        print(a,b,"主项、真空矩差区间",lead.str(11),change.lower().str(10),change.upper().str(10),flush=True)
        assert change*lead > 0 and abs(change) > arb("2.4e-9")
        assert abs(change-lead) < abs(lead)/8

    worst, witnessed = arb(0), arb(0)
    for a,b,lead,change,actual in vacuum_intervals(load, (numbers[-1], *HOLDOUT), source_only=True):
        relative = abs(change-lead)/abs(lead)
        assert abs(lead) > 0
        assert relative < arb(1)/8 and lead*change > 0
        worst = worst.max(relative.upper())
        witnessed = witnessed.max((abs(actual-lead)/abs(lead)).upper())
        print("前向", b, lead.str(10), change.lower().str(10), change.upper().str(10),
              actual.str(10), relative.upper().str(8), flush=True)

    assert worst < arb(".02721") and witnessed < arb(".00012045")
    print("十二处全域、冻结见证最大相对误差上界", worst.str(10), witnessed.str(10))


if __name__=="__main__":
    verify()
