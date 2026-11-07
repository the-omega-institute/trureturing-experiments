# Hardin's bounded mixed-difference array recurrences

`check.py` counts the original 0..3 arrays in OEIS A250538–A250544 using exact integer arithmetic. Run `python3 check.py` with Python 3 and NumPy. NumPy object arrays hold arbitrary-precision Python integers; no floating-point arithmetic is used.

`source-data.json` contains the attributed OEIS entry text and b-files, with source URLs and b-file SHA-256 values. OEIS material is subject to the [OEIS End-User License Agreement](https://oeis.org/LICENSE). `output.txt` is the program's complete deterministic output.

The checks reproduce all 957 b-file terms of A250538–A250543: the first four files through n=210, A250542 through n=87, and A250543 through n=30. The displayed source terms and all signed recurrence coefficients are parsed and matched. Each sequence counts (n+1)-by-(k+1) arrays; k=2,...,7 means actual widths 3,...,8.

For A250544, the check matches 140 of the 144 antidiagonal table entries using these counts and transposition. Four entries with both dimensions exceeding 8 remain untested; all available entries belonging to columns k=1,...,7 are tested. The distinct order-7 formula for k=1 is checked through n=32, and the common order-16 formula is checked for every sequence's available range, including the natural n=0 extension at n=16.

A separate actual-row transfer calculation retaining all offsets agrees for widths 2, 3 and 4 through n=24. Exact Horner matrix evaluation verifies the proposed annihilating polynomials for the 7-state width-2 and 37-state width-3 transitions. These finite computations corroborate the preregistration outline; they are not a Lean proof and do not establish the unbounded assertions.
