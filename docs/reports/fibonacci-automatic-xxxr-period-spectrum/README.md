# Primitive periods of binary words avoiding xx x^R

Run `python3 check.py`. It uses Python 3 integers and strings, no floating point and no dependencies. Successful output ends with `ALL_CHECKS_OK`; `output.txt` is the stored output.

The target is only the first assertion of Conjecture 11 in Du, Mousavi, Rowland, Schaeffer and Shallit, *Decision algorithms for Fibonacci-automatic words, II: Related sequences and avoidability*, TCS 657 (2017), 146–162, DOI [10.1016/j.tcs.2016.10.005](https://doi.org/10.1016/j.tcs.2016.10.005). It is Conjecture 47 in [arXiv:1406.0670v4](https://arxiv.org/abs/1406.0670v4).

The program reproduces all 48 entries with n ≤ 48 in Table 1 of Currie and Rampersad, *Growth rate of binary words avoiding xxx^R*, [arXiv:1502.07014v1](https://arxiv.org/abs/1502.07014v1), TCS 609 (2016), 456–468, DOI [10.1016/j.tcs.2015.11.004](https://doi.org/10.1016/j.tcs.2015.11.004). It also reproduces the periodic avoidance examples of DMRSS Theorem 12 for k = 1, 2, 3, 4, at periods 10, 42, 178, 754. The Rote-Fibonacci prefixes are generated using the source's transducer on the Fibonacci substitution word.

The complete search of primitive binary periods of length at most 48 gives exactly 2, 10, 42. Finite words are enumerated by appending one letter and rejecting newly created forbidden suffixes; every forbidden factor has a final letter, so this exhausts finite avoiders. Every periodic avoider must itself be a finite avoider. Each candidate's primitive root is computed directly.

The periodic decision has an exact cutoff. Let p be the least period. Starts can be reduced modulo p. If xx occurs with |x| ≥ p, comparing p letters shows that |x| is a period modulo p, hence p divides |x|. If the following x^R also occurs, the length-p block at that start is palindromic, and it already gives a forbidden pattern with root length p. Thus checking 1 ≤ |x| ≤ p suffices. An independent modular-index implementation agrees with the string-slicing implementation on every nonempty binary word of length at most 9.

Further bounded checks verify the suitable-pair equations for three substitution levels, the expected Fibonacci weights, examples of both terminal families, and the cyclic parsing implication from Currie–Rampersad's forbidden-factor list on all 131070 auxiliary words of lengths 1 through 16. The source's complement bar in the third suitable-pair equation is essential.

These computations support a structural proof plan; they do not establish the unbounded period classification. The extra prefix-conjugacy assertion at n = 0 is outside the target and fails because the Rote-Fibonacci prefix is 00 whereas the alternating primitive period is 01.

The check program is independently authored; the papers supply the definitions and comparison data.
