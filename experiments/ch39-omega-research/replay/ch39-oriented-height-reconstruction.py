"""Recover one independent rational orientation coefficient at each finite RH height."""
import json
from pathlib import Path


def phases(chart, denominator, prime):
    sign, n = chart
    inv = pow((denominator**2+n*n) % prime, -1, prime)
    return (sign*(denominator**2-n*n)*inv % prime, 2*sign*n*denominator*inv % prime)


def determinant(rows, prime):
    rows = [r.copy() for r in rows]
    n = len(rows)
    assert all(len(r) == n for r in rows)
    det = 1
    for j in range(n):
        pivot = next((i for i in range(j,n) if rows[i][j] % prime), None)
        if pivot is None:
            return 0
        if pivot != j:
            rows[j], rows[pivot] = rows[pivot], rows[j]
            det = -det
        d = rows[j][j] % prime
        det = det*d % prime
        inv = pow(d, -1, prime)
        for i in range(j+1,n):
            coeff = rows[i][j]*inv % prime
            for k in range(j+1,n):
                rows[i][k] = (rows[i][k]-coeff*rows[j][k]) % prime
    return det


def main():
    table = json.loads(Path(__file__).with_name("ch39-five-tick-record.json").read_text())
    den = 2**table["bits"]
    prime = 65537
    for channels in (1, 2, 4):
        rows = []
        for j in range((46 + channels-1)//channels):
            block = [[] for _ in range(channels)]
            for _, h, w, p, q in table["rows"]:
                ca, sa = phases(p, den, prime)
                cb, sb = phases(q, den, prime)
                # Imaginary parts of a, b, conj(a)b, ab.
                sines = [sa, sb, ca*sb-sa*cb, ca*sb+sa*cb]
                weight = w*h*h*pow(den*(h*h+(j*den)**2) % prime, -1, prime) % prime
                for row, sine in zip(block, sines):
                    row.append(weight*sine % prime)
            rows.extend(block)
        det = determinant(rows[:46], prime)
        assert det != 0
        print(channels, len(rows)//channels, det)


if __name__ == "__main__":
    main()
