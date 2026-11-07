# Alpha unpredictability and divergent Cesaro averages

Run `python3 check.py` with Python 3 and SymPy. All evidence uses integers,
Fractions, and exact symbolic algebra. Successful output ends in `ALL_CHECKS_OK`.

The source is Akhmet and Tilessov, *Integral of an Alpha Unpredictable Function*,
Mathematics 14(18), 3391, DOI https://doi.org/10.3390/math14183391, Section 3,
Remark 1. The journal question is absent from arXiv:2608.09809v2.

The fidelity check reproduces the Section 2 exponential-kernel mass 1/2 and the
associated linear ODE, which give the stated bound for inputs in [0,1].

The proposed counterexample uses a binary word built by concatenating every word
of length n and then an alternating constant block of length n^2(L_n+1). The
program checks exact frequency estimates through stage 40. It independently
expands stages 1 through 5 to cross-check its compressed representation.

The reflected step function is convolved with a normalized even kernel supported
in [-1/4,1/4]. Five return shifts are checked by complete finite central-block
comparisons, supplemented by exact convolution samples. The same shifts give
separation intervals: support-containment certificates establish the full interval
plateaux, supplemented by exact rational samples. The two mean conventions agree
for the even extension. Three affine parameter choices check the family scaling.

The exactly integrated test kernel is (15/4)(1-16u^2)^2 on its support. It is not a
C-infinity kernel. The support certificates also apply to a normalized smooth
bump. Smoothness, all infinite limits, and the complete cluster interval require
the proposed mathematical argument; these finite checks do not prove them. At
integer endpoints the test-kernel smoothing correction is exactly
(5/128)(a_T-a_(T-1)).

`output.txt` is the program output. No source-paper PDF or model conversation is
included.
