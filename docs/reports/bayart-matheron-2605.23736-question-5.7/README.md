# Bayart–Matheron Question 5.7: exact preregistration checks

Source: Frédéric Bayart and Étienne Matheron, *On the dynamics of composition operators: supercyclicity, odometers and translations*, arXiv:2605.23736v1, Question 5.7, Fact 1.2 and the basic-cylinder preimage formula preceding equation (1.2). Source text: https://arxiv.org/html/2605.23736v1 .

Run `python3 check.py`; Python 3's standard library is sufficient. `check.out` is the complete output. All assertions use integers and `fractions.Fraction`; no floating point enters the evidence.

Fidelity: 50,518 mixed-radix iterate checks reproduce Fact 1.2 for radices 2, 3 and 4, lengths 1–3, including wraparound, and check the source's cylinder-preimage formula.

Candidate: use binary probabilities μᵢ(0)=1/2+1/(4(i+1)), μᵢ(1)=1/2−1/(4(i+1)), i≥1. Exhaustion of all prefix atoms and cyclic powers for N=1–8 checks 87,380 probability ratios, the one-step maximum 5/3, the all-power maximum (2N+3)/3 and the exact inverse-density second moment. Telescoping and moment estimates are checked through N=128.

The finite checks support an analytic counterexample route on complex L²: uniformly bounded inverse-density moments allow Haar averaging on cylinder functions to extend continuously, giving a nonzero invariant functional and obstructing hypercyclicity. These checks do not prove infinite-product limits, density of cylinder functions, extension of the composition operator or non-hypercyclicity. Those obligations belong to the subsequent proof; no Lean proof or probe is included.
