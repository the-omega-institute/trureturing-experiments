"""Sharp charged singlet multiplicity for the Gaussian completion at the crossing."""
from pathlib import Path
from runpy import run_path

from flint import arb, fmpq, fmpq_mat

finite = run_path(str(Path(__file__).with_name("ch39-shell-gw-finite-mass-eigenvalue.py")))
q, base = finite["q"], finite["base"]


def disk_moment_floor():
    # bidegrees (i,j) with 0 <= i,j <= 8; angular integration groups i-j.
    least = fmpq(1)
    for difference in range(9):
        n = 9-difference
        block = fmpq_mat([[fmpq(1, 1+difference+i+j) for j in range(n)] for i in range(n)])
        for k in range(1, n+1):
            assert fmpq_mat([[block[i,j] for j in range(k)] for i in range(k)]).det() > 0
        inverse = block.inv()
        least = min(least, 1/max(sum(abs(inverse[i,j]) for j in range(n)) for i in range(n)))
    assert least > fmpq(1, 10**12)
    return least


def certify():
    events = []
    z = q("-.00675").union(q("-.0035"))
    T, B, A, margin, raw, time = base["check"](q(81, 977)**5, mass_path=events, bulk_portal=z)
    Y = (12*T).sqrt()
    heavy = q(12).sqrt()*raw.sqrt()*q(81, 977)**5*(6*T*time).exp()
    rate = 12*A+3*q(".020")+4*Y**2
    E = (8*Y*B+96*Y**3)*heavy*((rate*time).exp()-1)/rate
    assert E < min(q(".020"), 6*margin)
    smallest = q(1)
    for op, *args in events:
        if op == "threshold":
            k, _, st, _ = args
        else:
            _, st, _, _, _ = args
        if k != 53:
            continue
        aa, bb, ab = (st["a", pair] for pair in ((1,1), (2,2), (1,2)))
        assert (aa+bb)/2-ab > 0
        mu = min((st["l", 1]-aa/2).lower(), (st["l", 2]-bb/2).lower(),
                 *(val.lower()/2 for (kind, _), val in st.items() if kind == "b"))
        smallest = arb.min(smallest, (2*mu-E*q(2).sqrt()/6).lower())
    floor = disk_moment_floor()
    print("opposite-fiber witness on rank-two indefinite Gram matrices above", smallest.str(8))
    print("monomial disk moment floor above", floor)
    assert smallest > q(".001")


if __name__ == "__main__":
    certify()
