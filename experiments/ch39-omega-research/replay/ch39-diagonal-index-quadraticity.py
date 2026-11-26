"""Stabilizer character and the parallelogram identity for diagonal link fields."""
from itertools import combinations, product
from pathlib import Path
from runpy import run_path


def character(weights, index, modulus):
    return [-index[-x % modulus]+sum(index[(x+y) % modulus]
            for j, y in enumerate(weights) if i != j) for i, x in enumerate(weights)]


def verify():
    root = Path(__file__).parent
    flux = run_path(str(root/"ch39-local-triality-mixed-flux.py"))
    index = (0, sum(flux["EXPECTED"]), sum(flux["EXPECTED"]))
    assert index == (0, 16, 16)
    weights = (0, 1, -1, 1, -1)
    generator = (1, -1, 0, 0, 0)
    k = character(weights, index, 3)
    assert k == [64, 16, 16, 16, 16]
    charge = sum(u*v for u, v in zip(k, generator))
    assert charge == 48 and 3*charge == 144
    # Every heavy 10 and conjugate 10 pair contributes opposite continuous charge.
    pairs = list(combinations(range(5), 2))
    heavy = [(weights[i]+weights[j], generator[i]+generator[j]) for i, j in pairs]
    assert sum(q*index[s % 3] - q*index[-s % 3] for s, q in heavy) == 0
    for m, table, bad in [(3, index, 40), (5, (0, 50, 25, 25, 50), 460)]:
        witnesses = [w for w in product(range(m), repeat=5) if sum(w) % m == 0]
        assert len(witnesses) == m**4
        assert sum(len(set(character(w,table,m))) != 1 for w in witnesses) == bad
        for x,y in product(range(m), repeat=2):
            k = character((0,x,-x,y,-y),table,m)
            assert k[0]-k[1] == 2*table[x]+2*table[y]-table[(x+y)%m]-table[(x-y)%m]
    # Polarize over the rationals: a universal quadratic index cancels every weight.
    # Test coefficients separately on a generic determinant-one diagonal 5-tuple.
    from sympy import expand, symbols
    vals = list(symbols("x0:4"))
    vals.append(-sum(vals))
    assert all(expand(-x*x+sum((x+y)**2 for j,y in enumerate(vals) if j!=i)
                      -sum(y*y for y in vals)) == 0 for i,x in enumerate(vals))
    assert 4**3*16 <= 1250 < 4**4*16
    print("Continuous stabilizer charge 48 per generation, 144 for three; no discrete links.")
    print("Obstructed diagonal flux assignments: 40/81 at order 3; 460/625 at order 5.")


if __name__ == "__main__":
    verify()
