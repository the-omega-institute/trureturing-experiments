# Subcubic partial-Petrial counterexample

`check.py` computes the rank enumerator from the paper's diagonal-perturbation definition by exact Gaussian elimination over GF(2) on 100 small trees. It independently evaluates the 47-vertex witness and the 46-vertex leaf deletion with the rooted two-state recurrence, checks the published coefficient lists and total mass, and verifies the Jacobsthal tail identity and strict coefficient failure for the pendant-path family through 100 extensions. All arithmetic is integer arithmetic.

Run:

```sh
python3 check.py
```

The expected final line is `ALL_CHECKS_OK`.
