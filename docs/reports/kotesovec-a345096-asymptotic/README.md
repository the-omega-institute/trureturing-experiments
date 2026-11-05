# Exact checks for OEIS A345096

The source is [OEIS A345096](https://oeis.org/A345096), revision 10, whose NAME defines a floor/ceiling double sum and whose FORMULA conjectures a(n) ~ n²/2. Wesley Ivan Hurt contributed the entry; Václav Kotěšovec stated the conjecture.

Run with Python 3, with no third-party packages:

```sh
python3 check.py
curl -L https://oeis.org/A345096/b345096.txt -o /tmp/b345096.txt
python3 check.py --bfile /tmp/b345096.txt
```

The stored output uses the second invocation. The source table's SHA-256 is 30a4fef1f856bf78aa26601ec7f05f9e77cfa8955d55856f87fccb54d3526b2f. The table is available under the [OEIS license](https://oeis.org/LICENSE); the finite source values and formula used here are attributed to that entry.

`check.py` uses integers and `Fraction` throughout. It reproduces a(4)=15, the first ten terms, the prime-input formula at every prime through 1000, and all 1000 linked table entries when the b-file is supplied. The literal NAME formula and the divisibility reduction agree for n=1…80. For n=1…1000 it checks the grouped excess bound, the divisor-pairing bound, and the squared relative-error bound, without approximating square roots.

These are finite consistency and source-fidelity checks. They do not prove the conjectured infinite limit or certify a Lean proof.
