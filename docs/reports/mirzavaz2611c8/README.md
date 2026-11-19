# Möbius sign for multiplicative partitions at exponent pattern (7, 7)

Run the exact check with Python 3.12 or later.

```sh
python3 check.py
```

The program enumerates sorted multisets of nonzero exponent pairs. It obtains all upper covers by merging two entries and uses reachability to compute the ordinary integer Möbius recurrence from each endpoint. It reproduces the three diagonal entries `(2,2)`, `(3,3)` and `(4,4)` in the source's table. The `(7,7)` interval has 2,998 elements and bottom-to-top value 344, opposite the sign in Conjecture 8. The program also compares the cover graphs and computed values after relabelling to three distinct pairs of primes. For any distinct primes, unique prime factorization gives the same order isomorphism by replacing each vector with the corresponding prime-power product.

The expected output is in `check.out`.

Source. [Mirzavaziri and Yaqubi, arXiv:2611.18194v1, §4](https://arxiv.org/abs/2611.18194v1).
