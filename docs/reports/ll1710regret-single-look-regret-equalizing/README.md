# Single-look regret: failure of product weighting

Source: Thomas Lidbetter and Kyle Lin, *Searching for Multiple Objects in Multiple Locations*, arXiv:1710.05332v1, Section 4.2 (unnumbered conjecture), with the equalizing-property definition in Section 5. DOI of the journal version: 10.1016/j.ejor.2019.05.002. The journal full text has not been read: ScienceDirect returned `CPE00001` and the Elsevier full-text API returned `AUTHENTICATION_ERROR`. This report concerns the arXiv statement.

Run `python3 check.py` with Python 3 and SymPy. All evidence uses integers, rational numbers, and exact symbolic algebra; no floating point or numerical LP solver is used. The checked output is in `output.txt`. No Lean proof is supplied.

## Model and fidelity

The Hider chooses a k-element subset of n boxes. The Searcher chooses one of the n! fixed permutations and stops when every hidden object is found. Regret is the sum of costs of opened empty boxes, equivalently total search cost minus the sum of costs of occupied boxes. The program constructs the matrix both by simulating the search and by comparing positions in the permutation, and checks equality entry by entry.

The program reproduces the source's equalizing-strategy lemma (TeX label `le:Wnk`): product weights on every placement give expected regret k S_(k+1)/S_k against every permutation. It checks costs (2,3,5) with k=1,2 and costs (19,13,13,4,2) with k=2.

## Exact counterexample

For n=5, k=2, costs (19,13,13,4,2), placement order (12,13,14,15,23,24,25,34,35,45), the Hider distribution is (0,0,20,99,0,20,65,20,65,0)/289 and the regret value is 4975/289. The output contains all 120 lower bounds and a strictly positive eight-permutation Searcher mixture giving matching upper bounds on all ten placements. Its strict upper slacks force every optimum to vanish on 12,13,23,45. Normalization and five positive-mixture row equalities form a six-variable matrix with determinant 561340. Exact solution gives the displayed Hider distribution, proving uniqueness in ordinary finite-game arithmetic.

Both p14 and p24 are positive and equal, while their product weights have ratio 19/13. Therefore no support with normalized product weights can be optimal. The enumeration independently checks every one of the 1023 nonempty supports: the best product strategy has support {15,25,35,45}, guaranteed regret 843/49, and gap 148/14161 from the optimal value.

## A uniform family

Let c0=(19,13,13,4,2) and r=1/95726048. For any real c with |c_i-c0_i|≤r, all costs are positive. Every regret entry changes by at most 3r, since at most three empty boxes are opened. The fixed certified Hider distribution therefore guarantees at least 4975/289−3r.

For any nonempty placement support S, put w_H=∏_(i∈H)c0_i and w'_H=∏_(i∈H)c_i. With r≤1, each weight changes by at most 38r+r²≤39r. Thus Δ=Σ_(H∈S)|w'_H−w_H|≤390r, whereas W=Σ_(H∈S)w_H≥8. Both normalizers are positive. Direct subtraction of w'/W' and w/W gives an l1 distance at most 2Δ/W≤195r/2. Choose a permutation minimizing regret for the product distribution at c0. At c0 its entries lie in [0,45]. At c, its expected regret is at most 843/49+(3+45·195/2)r. Consequently the fixed certified Hider strategy beats every product-supported strategy by at least

148/14161 − (8787/2)r = 97614109/9381152704 > 0.

This argument covers every real vector in the closed cube, rather than merely sampled rational points. The program verifies the polynomial identity and all rational constants used in these uniform inequalities. Positive rescaling multiplies every payoff by the scale and every unnormalized product weight by its square, leaving normalized distributions unchanged. Thus every positive rescaling of the cube also consists of counterexamples. No uniqueness assertion is made away from c0.

## Scope and source discrepancy

This is the source's explicit fixed-permutation regret game. Adaptive-search variants are outside the claim. The compact candidate described its payoff as total search cost; that description is incorrect. Against the certified Hider distribution, the expected occupied-box cost is 5169/289, so the minimum expected total search cost is 10144/289, not 4975/289. The program checks and prints this discrepancy.
