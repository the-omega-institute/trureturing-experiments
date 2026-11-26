"""Finite-clock memory ladder and full visible poles of the rational relay block."""
from pathlib import Path
from runpy import run_path

from flint import fmpq_mat
from sympy import cancel, symbols

mem = run_path(str(Path(__file__).with_name("ch39-three-clock-minimal-memory.py")))


def verify():
    d, columns, *_ = mem["coefficients"]()
    b = fmpq_mat([row[:19] for row in columns.tolist()])
    power, basis = b, b
    ranks = []
    for _ in range(4):
        ranks.append(basis.rank())
        power = d*power
        basis = fmpq_mat([a+b for a,b in zip(basis.tolist(),power.tolist(),strict=True)])
    assert ranks == [19,38,53,57]
    char = d.charpoly()
    repeat = char.gcd(char.derivative())
    assert char.degree() == 57 and repeat.degree() == 19
    assert repeat.gcd(repeat.derivative()).degree() == 0 and char[0] != 0
    assert (char//repeat).degree() == 38
    x, z = symbols("x z")
    ladder = []
    for m in range(2,11):
        basis = sum(z**(m-1-j)*x**j for j in range(m))
        assert cancel((z-x)*basis-(z**m-x**m)) == 0
        fixed = ranks[min(3,(m-1)//2)]
        free = ranks[min(3,(m-2)//2)]
        ladder.append((fixed,free))
    assert ladder == [(19,19),(38,19),(38,38),(53,38),(53,53),
                      (57,53),(57,57),(57,57),(57,57)]
    print("Block Krylov dimensions", ranks, "clock grid minima, fixed/free", ladder)
    print("Visible relay poles: 19 simple and 19 double, none at zero.")


if __name__ == "__main__":
    verify()
