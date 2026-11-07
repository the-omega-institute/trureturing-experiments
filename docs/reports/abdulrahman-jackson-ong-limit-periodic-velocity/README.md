# Limit-periodic harmonic velocity: exact checks

Source: Houssam Abdul-Rahman, Thomas A. Jackson, Darren C. Ong,
[Bottleneck Effects and Harmonic-Type Velocity Bounds for Periodic Quantum Walks](https://arxiv.org/abs/2606.27740v1), Section 2.3, Theorem 2.5 and remarks (D)–(E).

Run `python3 check.py` with Python 3 and SymPy. The checked output is in `output.txt`.
Arithmetic uses integers, rational numbers and exact symbolic algebra throughout.
A star in the output denotes an adjoint; the coin is
`[[a,b],[-conjugate(b),conjugate(a)]]` and the walk is coin followed by shift.

The fidelity check reproduces remark (D): the transmitting left half-line sends
an initially left-moving state at zero to minus time, whereas the periodic
harmonic bound for the displayed approximants is 3/4 when delta is 3/5.

The proposed settlement uses an edge coordinate with increments
`h_j = min(R, 1/abs(a_j))`, extended by R at zero. The checker verifies its
increment identity symbolically modulo the unit-sphere equation, on 480 local
coin/cap cases, and on 3,600 full-line basis-vector cases. The latter cover every
transmission word of lengths 1 through 4 over {0, 3/5, 1}, three caps, and a fixed
varying complex phase pattern. It also checks exact two-sided period increments
in all 360 profile/cap cases.

These checks support an infinite-dimensional argument; finite enumeration does
not establish the limit-periodic theorem. Uniform-approximation means, operator
domains, telescoping, and the all-state limsup estimate are mathematical proof
obligations for the later Lean delivery. No Lean proof is included here.

The report targets only the ordinary reciprocal harmonic-mean bound. It does
not settle the neighboring-coin refinement in equation (2.21).
