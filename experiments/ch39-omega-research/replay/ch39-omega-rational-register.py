"""Check a finite projective register using only the Python standard library."""
from fractions import Fraction as F
from pathlib import Path

DEN = 2**48


def fidelity(k, chart, t):
    # w=C(k,DEN); D=(1+i)C(a,b). C is the signed Cayley chart from volume I.
    a, b = (t, DEN) if chart == 0 else (DEN, t)
    p = F(k*k, k*k + DEN*DEN)
    re_d = F(b*b - a*a - 2*a*b, a*a + b*b)
    return p**2 * ((3-2*re_d)/6) * ((6-4*re_d)/15) * ((113+72*re_d)/291)


def verify():
    table = Path(__file__).with_name("ch39-omega-rational-register.tsv")
    scores = {}
    for row in table.read_text().splitlines()[1:]:
        n, k, chart, t = map(int, row.split())
        assert chart in (0, 1) and abs(t) <= DEN
        scores[n] = fidelity(k, chart, t)
    assert list(scores) == list(range(10855, 12093))
    assert all(scores[10875] > x for n,x in scores.items() if n != 10875)
    peaks = sum(scores[i-1] < scores[i] > scores[i+1] for i in range(10906,12054))
    assert peaks == 213
    guaranteed = set(range(10905,12055)) | {10871,10875}
    for n,x in scores.items():
        if n != 10875:
            assert any(j in guaranteed and scores[j] > x for j in range(n-58,n+59))
    assert all(scores[11499] > scores[j] for j in range(11442,11557) if j != 11499)
    print("Rational register: 1238 rows, unique maximum 10875, 213 interior peaks, radius 58.")


if __name__ == "__main__":
    verify()
