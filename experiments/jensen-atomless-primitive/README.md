# An atomless primitive source with a signed first-defect node

The program constructs the same normalized positive polynomial density on
[-1,1], extracts its entire Fourier square-variable source, and checks its
finite Jensen rows with exact rational algebra and Sturm root counts.

Run from this directory:

    uv run --python 3.12 --with sympy==1.14.0 python exact_source.py --output results.json

The source has twelve initial rows with positive simple roots. Its degree
thirteen row has eleven positive simple roots and one nonreal conjugate
pair. At the unique actual degree-twelve node in (131/2,6551/100), the
unchanged derivative-source response has ratio strictly between -3/2
and -7/5.

The data also verify the density's Bernstein certificate and fifteen
coefficients of its derivative identity. The infinite positive-shift
statement uses the classical Bessel product and Laguerre--Polya/Jensen
suppliers. Arbitrarily long normal prefixes require the analytical argument in the owning Romik
source-interface note at Library/Analytic/romik2021orthogonal.md of
https://github.com/the-omega-institute/trureturing.
The numeric program does not certify those infinite statements.

The model is not the actual theta density, actual Xi exterior, full-growth
kernel configuration or required moving-band witness. It supplies no RH
or Robin counterexample, Lean certification or mathematical-priority claim.
Classical beta, Bessel and Laguerre results are reused.
