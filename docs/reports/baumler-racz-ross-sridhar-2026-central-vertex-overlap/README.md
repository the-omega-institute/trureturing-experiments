# Central vertex overlap in uniform attachment trees

Source: Johannes Bäumler, Miklós Z. Rácz, Nathan Ross and Anirudh Sridhar,
*Correlated uniform attachment trees*, arXiv:2606.02472v1, Definition 1.4,
Definition 2.2, Lemma 2.6 (proof), and Section 1.4, third bullet.

Run `python3 check.py`; Python 3.12 or newer and its standard library suffice.
`check.out` is the captured output. All evidence uses exact integers and rational
numbers; no simulation or floating point is used.

The source fidelity check compares every fringe-size law with the explicit
(1,k) Pólya-urn recursion stated in the source. It also computes centrality
by literal edge deletion and graph traversal for all parent maps through n=6,
independently of the fringe-size formula used by the other checks.

For every parent map through n=8 and every threshold, the program checks the
candidate finite tail identities, total threshold mean, covariance bound,
variance bound and expected centrality/fringe discrepancy at most 2.
For n≤6, it checks the rank squeeze for every early-label cutoff, budget and
threshold, using the exact minimum and maximum across all possible cutoff ties.
The three larger examples use n=K² and exact alternating-series enclosures of
exp(-1/2) to compare normalized finite means with 2(1-exp(-1/2)).

These finite checks support a proposed general proof; they neither prove the
conjecture nor certify a uniform asymptotic error rate. The proposed positive
settlement uses a general urn/concentration argument, not finite enumeration.

Producer: codex-cli research seat, flight f912fe4df2ed2a3e6edb00a7,
prior disclosure repo-prior-exposed; search candidate supplied by
nyxid-oracle gptpro-search-s942. Original check code; Apache-2.0 under the
repository licence. No third-party program or Lean source was copied.
