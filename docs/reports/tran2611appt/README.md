# APPT spectral-volume asymptotics checks

Source: Anh T. Tran, *Purity maxima, absolute separability, and spectral volumes of absolutely PPT states*, arXiv:2611.01492v1, Section 7 and the open determination in Section 8. The comparison with the qubit precursor uses arXiv:2609.18568v1, Theorem C and Section 4.

Run with Python 3.12 or later:

```sh
python3 docs/reports/tran2611appt/check.py
python3 docs/reports/tran2611appt/monte_carlo.py
```

`check.py` uses only integers and rational arithmetic. It checks the source's simplex-vertex coefficients, CAS factorial/product identities, outer-volume section formula and dominant summand cancellation for 63 pairs (2 ≤ m ≤ 8, m ≤ n ≤ m+8). It reproduces the precursor's exact outer formula for m=2 and 2 ≤ n ≤ 8, and checks beta_2=9/2. It checks every extreme-slot permutation for eight rational spectra at m=2,3 (17,304 matrices), by all principal minors. These spectra satisfy the proposed local PSD certificate. It checks finite exponential-moment identities used in the asymptotic argument. Finite checks do not establish the universal volume asymptotic.

The independent tilted variables require the support restriction B+Z ≥ 0. At m=n=2 their total normalization T_m is 1, whereas V_U=1/2; the unconditioned tilted law is therefore not the exact conditional outer law. The exact integral identity is V_K=T_m Q(B+Z≥0 and x in K), for K contained in the outer region. The proposed argument must control the excluded support with concentration. The finite mean of B need not be positive; its positive limit after division by d is an asymptotic assertion.

`monte_carlo.py` uses floating point exponential sampling under this tilted law with the support restriction. It is a separate diagnostic, with seed 261101492 and 200,000 samples at each n=2,...,8. CAS and U estimates are compared with their exact formulas using a prespecified six-standard-error tolerance. APPT is tested using the exact qubit criterion (evaluated in floating point). This is not exact evidence, a rigorous probability interval, or a proof of the general asymptotic.

Outputs are stored beside the programs. Both programs use the Python standard library only. No third-party implementation is used.
