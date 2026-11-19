# Zenger Question 2

Run `python3 check.py` with Python 3.12 or newer. The evaluator uses exact
finite sets, Boolean truth values, and integer indices.

The fidelity case reproduces Lukas Zenger’s Lemma 3.10 and Figure 3.1, where
a total functional model falsifies distribution of box over disjunction.
The settlement checks evaluate the even-step consequences on every dynamic
model on one or two labeled worlds for two atoms. They also check the
odd/even defect chains and stabilization of every one-atom formula of
size at most seven on the tested range. The preregistration gives the
induction for arbitrary formulas and path lengths.

[Source thesis, arXiv:2607.13528v1](https://arxiv.org/abs/2607.13528v1)
