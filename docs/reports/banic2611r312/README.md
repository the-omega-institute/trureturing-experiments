# Exact checks for Banić–Elezović cyclic dispersion

The target is arXiv:2611.14735v1, Remarks 3.12–3.13. The statistic is the sum of absolute successive multinomial count differences, including the wrap edge; its dispersion squared is Var(D)/N.

Run `python3 docs/reports/banic2611r312/check.py` with Python 3.12 or newer and SymPy 1.13 or newer. The adjacent `check.out` contains the complete output. Every evidentiary calculation uses integers, rational numbers or exact symbolic identities.

The program independently compares labelled histories with multinomial occupancy weights, reproduces the paper's raw edge moments, mean-transfer identity and cyclic variance decomposition, checks the N=4 variance profile through 24 bins, and checks its symbolic moment and derivative identities. It also counts the complete four-bin five-observation distribution. These checks give the strict reversals v(4,4)−v(3,4)=15199/186624 and v(4,5)−v(4,4)=−117/4096. The positive coefficients in the shifted derivative polynomial support the proposed proof that four bins uniquely maximize the N=4 profile.

Finite grid checks do not prove the all-bin multinomial moment formula. No Lean proof or mathematical settlement is claimed by this experiment.
