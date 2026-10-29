# Bala fractional factorial checks

Run `python3 check.py`. The recorded stdout is `output.txt`.

The source is Peter Bala, *Some integer ratios of factorials*, August 2016,
page 5, https://oeis.org/A276098/a276098.pdf, and OEIS A276100, A276101,
A276102. Fractional factorials mean Gamma(x+1).

For every index from 0 through 400 in each sequence, the checker cancels
formal Gamma symbols using the functional equation, compares the resulting
exact rational number with an independently formed progression product,
checks positivity and denominator one, and compares every ordinary-index
subsequence value with the classical Chebyshev factorial ratio. It reproduces
the first six OEIS terms of each sequence. Its digests encode all 401 values
as decimal integers, one per line, with a final newline.

The shared floor-count identity is tested by explicit counting for indices
0 through 80 and moduli 2 through 100 coprime to the sampling denominator.
The six displayed initial values of Radcliffe's distinct A364175 sequence
are also reproduced. Only Python integers and `fractions.Fraction` are used.

These are finite checks of transcription, identities, and proposed settlement
behavior. They prove neither universal integrality nor literature priority.
