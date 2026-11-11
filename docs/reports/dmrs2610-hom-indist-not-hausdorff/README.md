# Exact checks for the homomorphism distinguishing topology

Run `python3 check.py` with Python 3.12 or newer. The script uses only the standard library and integer or rational arithmetic; `output.txt` is its output. Successful completion prints `ALL_CHECKS_OK`.

The fidelity check reproduces the source's connected-guest identity `hom(Z,G+H)=hom(Z,G)+hom(Z,H)` on all 44 connected labelled simple guests of orders one through four and all 144 ordered pairs of simple targets of orders zero through three. The source is van Dobben de Bruyn, Marquès, Roberson, Seppelt, Spitzer and Zeman, *Topological and Geometric Perspectives on Homomorphism Indistinguishability*, arXiv:2610.07927v1, Section 3, the equation labelled `eq:connected` in its TeX.

The settlement checks compare direct homomorphism enumeration with integer matrix-power traces for every labelled simple target of order at most four and cycle lengths three through eight. Two explicit pairs have exactly equal characteristic polynomials and agree on cycle lengths four, six, eight and ten: the star K1,4 versus C4 plus one isolated vertex, and K1,6 versus K2,3 plus two isolated vertices. A different-order pair checks that zero eigenvalue multiplicities do not affect positive-length traces. The script also checks sparse moment signatures on the 76 small targets and exact nonsingularity of 16 positive-base moment matrices.

These finite checks support the proposed general argument. They do not prove rigidity from infinitely many moments, the cofinite induced topology, or non-Hausdorffness. Those require the general proof proposed in the preregistration.
