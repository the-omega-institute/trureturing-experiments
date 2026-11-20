# Salazar Problem 8.2

Run `python3 check.py` and compare stdout with `check.out`. The script uses Python's standard-library rational arithmetic.

The fidelity check compares the GRWS samples with the explicit ordinary Bernstein interpolant and its atomic Lévy representation, using a finite sum and an exact geometric tail. The settlement check tests the formula and geometric bound for the increments, then checks that the lower bound from positive density on a fixed interval eventually exceeds that upper bound. The interval checks use rational coordinates for the exponential of the negative endpoints and a rational lower bound for the interval length. The test grid includes (D = 0) and the constant (N = D) boundary. The preregistration gives the argument for arbitrary admissible real parameters.

[Salazar, arXiv:2607.09868v1, Definition 2.1, Theorem 7.1 and Problem 8.2](https://arxiv.org/abs/2607.09868v1)
