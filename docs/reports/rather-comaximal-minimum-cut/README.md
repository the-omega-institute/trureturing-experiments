# Rather's comaximal minimum-cut classification

Source: Bilal Ahmad Rather, *Vertex connectivity of the nonzero nonunit core of the comaximal graph of Z_n*, [arXiv:2605.01565v1](https://arxiv.org/abs/2605.01565v1), Problem 8.1.

Run `python3 check.py`. The deterministic output is in `check-output.txt`; success ends with `ALL_CHECKS_OK`. Only Python's standard library and exact integer arithmetic are used.

The program constructs the nonzero nonunit residue graph directly from gcd, independently of the support-model proof. Fidelity checks reproduce Proposition 3.2's fiber sizes, Theorem 3.3's adjacency criterion, Lemma 5.1's disconnecting class, and Table 3's minimum-degree values for n = 6, 30, 42, 70, 210, 2310. The table's connectivity values are checked by exhaustive deletion only for the entries at most 70.

For every squarefree composite n at most 70, the program enumerates every deletion with size at most the proposed connectivity and finds exactly one disconnecting set, the largest-prime singleton zero-set class. It independently reproduces the search seat's count of 149844 deletions of size delta on n = 6, 10, 15, 30, 35, 42, 70.

These results check definitions and finite instances. The proposed positive settlement requires the general surviving-transversal argument; no universal conclusion follows from this enumeration.
