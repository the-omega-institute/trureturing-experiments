# Kurkov A225620 partition-code recurrence

`python3 docs/reports/kurkov-a225620-partition-code-recurrence/check.py` runs offline with Python 3.12 or later and the standard library. All mathematical comparisons use exact integers. Successful execution ends with `ALL_CHECKS_OK`; `output.txt` contains the result.

The DATA fixture is the complete 66-term DATA field of [OEIS A225620](https://oeis.org/A225620), revision 98. Omar E. Pol's definition replaces each positive partition part j by j−1 zeros followed by one 1, concatenates the blocks, reverses the string, and reads it in base 2. The program implements these instructions independently of the cumulative-power formula. It reproduces the example C(2,2,1,1)=58 and all 66 published DATA terms exactly.

The second construction implements Mikhail Kurkov's COMMENTS conjecture: b(2^k)=2^k; otherwise b(n)=2^(binary_length(n)−popcount(n))·(2b(n−2^floor(log2(n)))+1), for positive n. It checks the identity b(n)=C(s(n)), injectivity, both decoding steps, the zero-count identity, and highest-bit deletion for every 1≤n<32768. Independently generated partitions of every total from 1 through 15 give 683 cases, with the standard partition counts 1,2,3,5,7,11,15,22,30,42,56,77,101,135,176. The tested b-values below 4096 are exactly the 271 partition codes below 4096; the same comparison also passes using only n<8192.

For decreasing set-bit exponents e₁>⋯>eₜ, the staircase transform is sᵢ=eᵢ−t+i+1; its inverse is n=Σᵢ2^(sᵢ+t−i−1). This is the known [A125106](https://oeis.org/A125106) encoding. Its bijectivity is discussed in Tikhonov et al., [LANTERN](https://arxiv.org/abs/2609.32264v1), Appendix F, item 1. The check program is original, using these mathematical definitions; it does not copy a third-party implementation.

These are finite numerical checks supporting a proposed unbounded proof. They are not a Lean kernel certificate. The separate A225620 FORMULA conjecture involving A228354 is outside the check.
