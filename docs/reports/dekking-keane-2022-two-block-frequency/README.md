# Dekking–Keane arXiv:2202.13548v1 Conjecture 4: frequency of 1 in the Thue–Morse two-block fixed point

`check.py N` (plain Python, integer arithmetic) generates x^(00) from the source's relations x_{3n} = x_{2n}, x_{3n+1} = x_{2n+1}, x_{3n+2} = 1 − x_{2n+1} with x_0 = x_1 = 0, and checks the relations for every n < N/3. It then:
- prints the source's prefix 001110101101110010110001101110001;
- computes the coefficient functional d_k of T^k, where T(u, v) = (u, v, −v) on ±1 pairs, for k = 0..10;
- checks E_k = Σ d_k² = 3^(k−1) and the cyclic adjacent correlation R_k = 0 for every k ≥ 1;
- reports the frequency of 1 and max |S_N|/N^α with α = log √6 / log 3.

`check.log` holds `check.py 1594323` (N = 3¹³), with frequency 0.500266 and exit 0.
