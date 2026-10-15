# Exact common-unbiased root cover on an enlarged strict-X rectangle

This experiment certifies one explicitly specified two-real-parameter rectangle. It does not classify the whole strict-X family and does not establish the global dimension-six four-MUB conjecture.

## Exact statement

Write alpha=a+i b and set

    |a+1/5| <= 2^-24,   |b| <= 2^-24.

For p_alpha(t)=(1+a)t^3+b t^2+(a-3)t+b, order its three real roots increasingly and take them in the order (smallest, largest, middle). Apply the same ordering to p_minus_alpha. Their Cayley images c(t)=(1+i t)/(1-i t) are r0,r1,r2 and s0,s1,s2. Define

    A=circ(r0*r1,r1,1), B=circ(s0*s1,s1,1),
    H=[[A,B],[B*,-A*]], H0=L H,

where L reorders rows as (2,0,1,4,5,3), with zero-based indices. At alpha=-1/5 this gives the exact regular seed with b0=(-3+4i)/5 and e0=(-2+i sqrt(21))/5.

Every common-unbiased ray is contained in one of 60 validated uniqueness boxes. The complete orthogonality graph is contained in the supplied 114-edge graph: one canonical six-vertex block, disjoint from a bipartite graph on the other 54 vertices. The canonical six branches are eigenrays of the simultaneous three-cycle S.

Consequently every actual common-unbiased orthogonal triple has total mode collision exactly 3. Every actual six-vector completion has symmetry affinity exactly 2 and uses the same six canonical rays, so two such completions cannot be mutually unbiased. These are uniform statements over the displayed closed rectangle, not sampled parameter tests.

The previous preserved certificate covered half-width 2^-32. The present half-width is 256 times larger; the absolute parameter rectangle is still small. The complete older proof is reused and attributed in provenance.json. No world-novelty claim is made.

## Reproduce

Requirements: Python 3 standard library and g++ supporting C++17 and __int128. Run from the repository root:

    python3 experiments/mub-strict-x-rectangle/verify_enlarged_patch.py \
      experiments/mub-strict-x-rectangle/input.json \
      --output ./mub-strict-x-replay --jobs 4

The output directory must not already exist. The input's parameter_radius records the historical 2^-32 center certificate. Its historical claim-status fields are preserved metadata, not acceptance evidence; centers, charts and allowed edges are untrusted proposals that the verifier must check. The new driver's explicit target is 2^-24; it first reconstructs and revalidates the old local root boxes for symmetry anchoring, then constructs the enlarged matrix enclosure and executes the new cover. It never accepts a stored PASS report as an input.

The published replay in results/delta-2neg24/ has:

- all 32 closed signed-Cayley charts covered;
- 5,992,724 visited nodes, zero pending and zero unresolved boxes;
- 60 strictly validated uniqueness roots, radius 2^-12;
- 1,770 distinct-ray pair checks and 1,656 forbidden-pair checks;
- every forbidden squared overlap greater than 1/500000.

Node counts are diagnostics, not mathematical novelty scores. The cap 1 centre preflight intentionally reports INCOMPLETE after validating all 60 local root boxes. It is a negative control, not an incomplete enlarged global traversal. The final 32 enlarged charts must all report COVERED.

The independent audit program in audit/ uses a separately implemented arbitrary-precision integer interval arithmetic. Its dyadic preconditioners are untrusted proposals, not assumed inverses. It verifies the enlarged root inclusions, graph bounds, report identities, fresh center fixation, and center-box nesting. From the experiment directory, run:

    python3 audit/verify_independent.py \
      --result results/delta-2neg24 --input input.json \
      --preconditioners audit/preconditioners.txt \
      --output independent-verification.json

The audit checks the fresh global chart reports and reviewed cover source; it does not independently replay the global search trace. Run without Python optimization (`-O`).

## Proof interpretation

Every unit coordinate is represented in one of the two compact signed-Cayley charts epsilon*(1+i t)/(1-i t), t in [-1,1]. Their 32 products cover the dephased phase torus including boundaries. Sign changes use t -> -1/t where defined.

The first five components of |H0* u|^2-6 are solved; the sixth follows from the exact Hadamard Gram identity. A residual exclusion or disjoint Krawczyk image rejects a box. Intersection with the image retains every zero. Closed bisection covers the parent. A box may otherwise terminate only inside one of the validated uniqueness neighborhoods, including a valid chart transition.

For each root box X, center m, and proposed rational C, the uniform checks are

    K=m-C f(m)+(I-C J(X))(X-m) strictly inside X,
    ||I-C J(X)||_infinity < 1.

The contraction theorem gives one root in each box for every exact parameter. All arithmetic used for acceptance is outward-rounded integer or Fraction arithmetic; floating-point inversion only proposes C. The computed uniform bounds include

    max contraction <= 385012914919/2^40 < 351/1000,
    max displacement <= 206508646/2^40 < (77/100)*2^-12.

The graph stage checks every pair at the same uniform parameter enclosure. Allowed edges need not persist: removing edges cannot create a triangle outside the canonical block.

Canonical shift fixation propagates analytically. The old center roots are freshly certified shift-fixed. Their radius 2^-16 boxes, using the same dyadic centers, lie strictly inside the newly validated radius 2^-12 boxes, identifying the actual branches at the center. A continuous simple-root branch's shift-fixed locus is closed; local uniqueness makes it open. The parameter rectangle is connected, so fixation holds throughout. This replaces an unnecessarily conservative shifted-interval-image test without assuming global modewise orthogonality.

Each shift eigenray has mode probabilities a permutation of (1,0,0). The graph therefore gives collision 3 for every orthogonal triple and affinity (6-2)/2=2 for every actual completion.

## Provenance and verification limits

The historical source is pinned to trureturing commit 63a005da4c695c3136da0bac01e3c138e7338840; its global cover was introduced at 3b0a12fea60aa281645b786cb5e1e87db58afdd0. The four historical source files and input.json retain exact original bytes. provenance.json records their Git blob and SHA-256 identities.

The only C++ arithmetic-core change is the uniqueness-box radius. The new driver checks that byte-level relationship before execution. It adds the enlarged parameter enclosure, fresh symmetry anchoring, complete traversal, and stronger graph margin. A separate independent arithmetic implementation checks the root and graph stages.

This is an external exact computational interval proof together with the stated analytic transfer. No Lean elaboration, kernel admission, full strict-X cover, or source-wide originality claim is asserted. Repository licensing applies; no additional license terms are introduced here.
