"""五荷共用三态与逐荷独立三态的单步平均能量差，见同目录说明。"""
from functools import cache
from pathlib import Path
from runpy import run_path

from flint import acb, arb

SHELL = run_path(str(Path(__file__).with_name("ch39-torus-cross-charge-light-leakage.py")))
CHARGES = SHELL["CHARGES"]
DAGGER = SHELL["dagger"]
DENOMINATOR = 4096


def clock_blocks():
    blocks = {}
    for q in CHARGES:
        fixed, first, second, soft = SHELL["parts"](q)
        # For either clock, the phase-dependent block joins disjoint vertex sets.
        # Its spectral norm is bounded by sqrt(max row sum * max column sum).
        for clock in (first, second):
            row_sums = [sum(abs(clock[i, j]) for j in range(19)) for i in range(19)]
            col_sums = [sum(abs(clock[j, i]) for j in range(19)) for i in range(19)]
            assert all(row == 0 or col == 0 for row, col in zip(row_sums, col_sums))
            assert all(row < arb("1.0001") for row in row_sums)
            assert all(col < arb("2.0001") for col in col_sums)
        assert all(sum(abs(soft[i,j]) for j in range(19)) <= 6 for i in range(19))
        blocks[q] = [(DAGGER(b)*fixed*b, DAGGER(b)*first*b, DAGGER(b)*second*b)
                     for b, _ in SHELL["BASES"]]
    return blocks


def light_energy(matrices):
    levels = [z.real for matrix in matrices for z in matrix.eig(multiple=True)]
    assert len(levels) == 19 and all(z.is_finite() for z in levels)
    lower = sum(sorted(z.lower() for z in levels)[:3])/3
    upper = sum(sorted(z.upper() for z in levels)[:3])/3
    return lower.union(upper)


def verify():
    blocks = clock_blocks()
    tau = 2*arb.pi()
    # Norm per differentiated clock adjacency < 1.415. The 3--4 gap of each
    # charge exceeds 0.915 throughout the torus at zero softening.
    norm = arb("1.415")
    assert (arb("1.0001")*arb("2.0001")).sqrt() < norm
    curvature = sum(q*q for q in CHARGES)*(2*norm + 2*norm**2/arb(".915"))/len(CHARGES)
    soft_error = arb("12e-8")

    @cache
    def phase(index):
        return acb(0, tau*index/DENOMINATOR).exp()

    def deficit(u, v):
        charges = []
        for q in CHARGES:
            a, b = u**q, v**q
            matrices = []
            for fixed, first, second in blocks[q]:
                clock = first*a + second*b
                matrices.append(fixed+clock+DAGGER(clock))
            charges.append(matrices)
        shared = [sum(sector[1:], sector[0])/len(CHARGES) for sector in zip(*charges)]
        return light_energy(shared)-sum(light_energy(m) for m in charges)/len(CHARGES)

    @cache
    def vertex(i, j):
        return deficit(phase(i), phase(j))

    witness_phase = acb(0, tau/48).exp()
    witness = deficit(witness_phase, witness_phase.conjugate())
    assert witness+soft_error < arb(".09236")
    print("upper witness", witness, flush=True)

    def check_cell(i, j, width):
        values = [vertex((i+di)%DENOMINATOR, (j+dj)%DENOMINATOR)
                  for di in (0, width) for dj in (0, width)]
        error = curvature*(tau*width/DENOMINATOR)**2/4 + soft_error
        if all(value-error > arb(".09") for value in values):
            return 1
        assert width > 1, (i, j, values)
        half = width//2
        return sum(check_cell(i+di, j+dj, half) for di in (0, half) for dj in (0, half))

    width = DENOMINATOR//32
    leaves = sum(check_cell(i, j, width)
                 for i in range(0, DENOMINATOR, width)
                 for j in range(0, DENOMINATOR, width))
    print("certified cells", leaves, "vertices", vertex.cache_info().currsize)


if __name__ == "__main__":
    verify()
