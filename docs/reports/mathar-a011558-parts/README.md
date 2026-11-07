# Mathar A011558 distinct-part partitions

Run `python3 check.py` with Python 3. It uses only the standard library and exact
integer arithmetic. `output.txt` is its complete output; success ends with
`ALL_CHECKS_OK` and exit status 0.

Sources: R. J. Mathar’s COMMENTS conjecture in
[OEIS A011558](https://oeis.org/A011558), revision 76, and Paul Barry’s
[OEIS A084215](https://oeis.org/A084215), revision 51. The short coefficient
lists embedded in the program are the entries’ DATA at offset zero. Source text
and data are attributed to their OEIS contributors; OEIS content is available
under [CC BY-SA 4.0](https://oeis.org/wiki/The_OEIS_End-User_License_Agreement).

The fidelity checks reproduce all 106 and 34 displayed source coefficients,
respectively, from the two original rational generating functions. They also
reproduce A084215’s stated initial values and doubling recurrence. Allowed parts
are obtained from that recurrence. Subset-sum DP and an independent enumeration
of all 8192 subsets of the 13 parts at most 10000 agree with A011558’s original
modulo-two coefficient recurrence for every `0 <= n <= 10000`, including the
empty partition at zero. Counts are zero for 2000 inputs and one for 8001 inputs.

These are finite checks, not a proof of the assertion for all natural numbers.
No Lean proof or kernel-verified resolution is contained here.
