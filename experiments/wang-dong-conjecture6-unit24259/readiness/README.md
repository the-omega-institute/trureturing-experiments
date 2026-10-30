# Wang–Dong Conjecture 6: bounded synthesis readiness

This private experiment concerns the full literal Wang–Dong Conjecture 6,
arXiv:2610.39833v1, with the source scope stated below.
The source event asserts both binary and arbitrary-alphabet maxima equal
`min(n-r, 3*m-1)` for every admissible `m>=4,n,r`. The source is cited only; no paper text is included. See ../SOURCE-LICENSE.md. Original code here is
Apache-2.0 under the repository LICENSE. The isolated, unvendored dependency
`z3-solver==4.15.4.0` reports MIT License in its installed package metadata;
upstream is https://github.com/Z3Prover/z3.

The exact one-sided search is for binary total tables on 19 states, m=6,
S={0,...,12}, A={0}, s=13, r=1 and B=17. A candidate must be strongly connected,
have a reset word of length at most 220 and an actual reaching word Qv=S of
length at most 126, and preserve both A inclusion and cardinality 13 for every
original-letter word of length at most 17, including the empty word. Witness
limits restrict the search, not the source event. SAT admissibility alone is
not a candidate, and finite UNSAT would not settle the conjecture.

## Exact encoding and separation

All 38 transition entries are independent five-bit variables constrained only
to the state domain; no restricted-table mutation or transition template is
used. Subsets are 19-bit masks. Symbolic images union exactly the one-hot
successors of present states. Forward and reverse closures from state zero
through 18 edges encode strong connectivity. A symbolic binary word path from
Q contains S at some index 0 through 126, encoding actual image reachability.
Another path has singleton image at length 220. Total deterministic letters
preserve singleton cardinality, so an exactly-220 reset witness exists iff one
exists of length at most 220.

The initial SMT formula contains bounded admissibility and the empty-word
no-success constraint. Other word constraints are separated lazily. Each SAT
model is concretely replayed for totality, strong connectivity, reaching and
reset witnesses, then checked by a complete breadth-first search of subset
images through B. A violation means A is missing OR cardinality is below s;
the next exact cut requires BOTH inclusion and cardinality equal to s for that
word. Merging identical subsets at minimum depth is exact because all future
images depend only on the subset, and the shortest arrival leaves at least as
much remaining word budget. A complete surviving traversal exports its nodes
and both outgoing edges below depth B. A final survivor alone could become a
whole-source counterexample candidate. No such survivor was obtained.

## Measured result and limitations

`measurements.json` contains the measured finite result, not a transcript.
The initial formula has 766 declared variables (38 BV5, 382 BV19, 346 Boolean),
427 top-level assertions and 2,722,502 uncompressed SMT-LIB bytes. The compressed
initial formula is retained as `initial-admissibility.smt2.gz`; it is an
admissibility relaxation of the all-word search, with the exact oracle supplied
in the program. Formula construction took 0.546 seconds and serialization
1.438 seconds. Z3's first SAT check took 9.164 seconds. The decoded model passed
concrete admissibility checks but failed on word [0], whose image mask 110577
has size 12 while still containing A. Adding the exact cut exercised the
solver/oracle/refinement integration. The second check returned unknown
(reason `canceled`) at its 20-second timeout. There was no UNSAT conclusion.
The measured solver subprocess took 31.466 seconds; including supervision,
32.592 seconds. Peak process RSS was 1,272,627,200 bytes on macOS;
Z3 reported max memory 993.04 MB. These figures describe this one process and
its two checks, not the entire worker or any projected full search.

Finite verification compared the complete oracle with literal enumeration on
128 two-state cases, compared 96 symbolic images with direct images, and
compared the fixed-table admissibility encoding on all 16 binary two-state
tables. A separate small SAT model was decoded and replayed. The explicit
four-state cycle/merge control has ac=1, so it fails strict ac>1. These checks
establish tested behavior only. The separate closed fixed-table audit is
retained in ../audit.py and ../audit-results.json; it was not rerun during retention. The rejected 19-state table was not retained, and no exact
admissible full counterexample is claimed. Progress beyond the first cut,
complete all-word closure on a surviving model, independent replay, kernel
closure, admission, publication and delivery remain undone.

To reproduce, create an isolated Python 3.13 environment, install
`requirements.txt`, and invoke `readiness.py --out RESULT.json --smt FORMULA.smt2
--seconds 75`. Each solver check has a 20-second timeout and Z3 has a 1536 MB
memory setting; the measured run also had an external 125-second process
limit. The program's seconds option is checked between solver calls, so a
caller must supply an external deadline as well. Do not place caches or
virtual environments in the main mathematical repository. These artifacts
make no claim of conjecture resolution, KPI credit, general solver feasibility
or a complete search bound.
