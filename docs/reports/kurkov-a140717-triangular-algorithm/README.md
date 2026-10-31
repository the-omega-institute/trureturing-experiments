# Kurkov's A140717 triangular algorithm

Run `python3 check.py`. The program uses exact integers and coefficient tuples for polynomials. It exits zero and prints `ALL_CHECKS_OK` after checking:

- OEIS A140717's printed rows for semilengths 0 through 6, including T(4,2)=5.
- The increasing peak-pair encoding and its explicit inverse for all 23,714 Dyck paths of semilength at most 10, preserving the sum of peak heights minus the number of peaks.
- The peak-pair rectangle recurrence by independent enumeration for 0 <= i < j <= 10.
- Every intermediate triangular-update state for vector lengths 1 through 30, its agreement with that recurrence, and stability of all saved output prefixes.
- The conjectured output identity by direct Dyck enumeration for all n <= min(m,10), m <= 30; row sums and maximum degrees through 30.
- Agreement with the inversion distribution on all 321-avoiding permutations through size 8.

These finite checks support a general proof outline; they do not prove the unbounded conjecture. The target is the first FORMULA conjecture in <https://oeis.org/A140717/internal>, by Mikhail Kurkov, dated August 30, 2024. The direct-update version also appears in UPD1 of <https://mathoverflow.net/questions/477465/elegant-algorithm-for-a140717>.

For the recurrence, let F(i,j) count sequences of pairs (a,b) with strictly increasing a and b, 1 <= b <= a < j, b <= i, weighted by the sum of a-b+1. Removing the last pair with b=i gives F(i,j)=F(i-1,j)+sum_{a=i}^{j-1} z^(a-i+1)F(i-1,a). The vector updates implement this recurrence. A nontrivial peak with u up-steps and d down-steps before its down-step gives (a,b)=(u-1,d+1); height-one peaks have weight zero and are reconstructed at the baseline.

Two errors in the source examples do not affect the target assertion: 2341 has three inversions, so the listed inversion-two example should be 2314; the printed Dyck word UDUUUDD has seven steps, whereas UDUUUDDD has semilength four and statistic two. Both corrections are checked explicitly.

Source definitions and numerical rows are attributed to OEIS A140717 (Emeric Deutsch, with Kurkov's conjecture), available under the OEIS Foundation's published license terms. The check program is independently written and uses only Python's standard library. Its captured stdout is `output.txt`.
