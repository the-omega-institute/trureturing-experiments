# Relational quantum automaton certificate

This experiment is the executable finite slice of the relational quantum
automaton proposed in the paper. It instantiates the smallest state space
needed to make the closure chain explicit:

1. one normalized unit state;
2. a binary amplitude split (a,b) in C^2, with |a|^2 + |b|^2 = 1;
3. three real relational coordinates
   (x,y,z) = (2 Re(conj(a)b), 2 Im(conj(a)b), |a|^2 - |b|^2),
   which lie on the unit sphere.

The transition alphabet is H, S, X, Z. The program enumerates the projective
orbit of |0>, checks the gate relations H^2 = X^2 = Z^2 = I and S^4 = I,
and verifies normalization, the Bloch-sphere identity, and the Born
readout formula p0 = (1+z)/2, p1 = (1-z)/2. Global phase is quotiented
when counting states, so the orbit contains six physical readout states.

Run from this directory:

    python3 verify.py
    python3 -m unittest -v test_verify.py

The default command writes results.json. It uses only the Python standard
library and double precision with explicit tolerances. It certifies these
finite identities for the declared gate set. It does not claim to derive all
of quantum mechanics, prove the infinite network limit, or establish the
C128 and FIB theorems. Those are separate contracts in the main
[trureturing](https://github.com/the-omega-institute/trureturing) repository.


## Finite closure trace

The certificate also emits a closure_trace. Each round reports the number of
visible projective records, the relation-edge count, the PASS/WAIT/FAIL status
of R0 through R5, and whether the finite family has reached a fixed point.

For the declared H/S/X/Z generator set, the trace is:

| round | records | edges | passed rules | fixed point |
| ---: | ---: | ---: | --- | :---: |
| 0 | 1 | 0 | R0, R4 | no |
| 1 | 3 | 1 | R0, R1, R4 | no |
| 2 | 5 | 3 | R0, R1, R2, R4 | no |
| 3 | 6 | 6 | R0, R1, R2, R3, R4, R5 | yes |

This is a finite fixed-point certificate for the six-state projective orbit.
R3 and R5 are still network proxies for universal successor-kernel and
independent-composition theorems. The output does not claim the full infinite
closure or a kernel proof of quantum mechanics.
