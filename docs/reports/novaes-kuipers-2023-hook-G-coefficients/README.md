# Hook G-coefficients of Novaes and Kuipers

Source: Marcel Novaes and Jack Kuipers, *Effect of a tunnel barrier on time delay statistics*, arXiv:2302.09030v1, Section II.6 (conjecture before Eq. (46)), coefficient definition Eq. (37). https://arxiv.org/abs/2302.09030v1

Run with Python 3, without third-party dependencies:

    python3 check.py

`check.py` counts Littlewood–Richardson tableaux directly, using weak rows, strict columns, prescribed content, and the lattice-word prefix condition. It sums Eq. (37) with integer coefficients, treating t = sqrt(gamma) as a formal indeterminate. Homogeneity bounds the partition sizes, and the empty partition is included.

`check.out` reproduces the source example for nu = (3,1,1,1): precisely rho = (3,1,1,1), (2,1,1), (3,1), (2), with G = 1, t^2, t^2, t^4, respectively. It compares the proposed closed form with direct LR computation for all 1296 ordered pairs of 36 nonempty hooks of sizes 1 through 8. All comparisons pass: 104 are nonzero and 1192 are zero. The cache contains 3504 LR triples. The last output line is ALL_CHECKS_OK.

These are finite checks of the stated identity, not a proof for unbounded partitions.

SHA256 of check.py: f2a554b93f9aeaecfe02fc63c7072f5bd99b40bf4d8173c19cc2640c60e842f1
