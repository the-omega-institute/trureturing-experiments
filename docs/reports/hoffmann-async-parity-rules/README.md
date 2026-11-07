# Hoffmann asynchronous parity candidate rules

Source: Rolf Hoffmann, *Towards New Parity Rules with Asynchronous Updating Inspired by Fatès's Rule*, arXiv:2610.14026v1. The source defines a family of 24 rules in Section 3.3, using the active neighborhoods of Table 1 and simultaneous complementation of a selected adjacent pair.

Run with Python 3 (standard library only):

```sh
python3 docs/reports/hoffmann-async-parity-rules/check.py
```

`check.out` contains the exact output. All evidence uses integer arithmetic.

The fidelity check reproduces the four-by-four HC transition-count matrix for n=5 in Appendix 2 (divide each entry by 5 for its transition probability). It also checks all five block-count examples of Section 3.1 and the equivalence of HC's formula to its four active neighborhoods.

The graph check enumerates all 24 rules at n=1,3,5,7,9,11. Every selected-site transition preserves parity, both constants absorb, and reverse breadth-first search from the constants reaches every vertex. This covers 65,520 rule/configuration pairs.

A separate constructive check erases an even-length cyclic run two cells at a time using HC, and checks the prescribed sites against the four-neighborhood update. It runs on every binary word of odd size through 15 and verifies the reflected schedule under mirrored HC. The path length is checked against n². This covers 43,690 initial words in both orientations.

These are finite checks, not Lean kernel verification or an all-size correctness proof. No timing, asymptotic scaling, expected-time estimate, or simulation-derived probability is asserted by this program.
