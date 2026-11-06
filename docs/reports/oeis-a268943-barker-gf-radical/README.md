# Barker's OEIS A268943 formulas

Python 3, standard library only. Run:

```sh
python3 docs/reports/oeis-a268943-barker-gf-radical/check.py
```

The program independently counts words by their actual last letter and last adjacent-equality event label. It compares these counts with the proposed rational generating function for alphabets of sizes 2 through 10 through length 210. Separate exhaustive enumeration checks 430,982 words. Exact arithmetic in Q[s]/(s²−17) checks Barker's radical expression through index 210. The output ends with `ALL_CHECKS_OK`. These finite checks support a proposed proof; they are not a proof for all lengths or a Lean verification.

Source fidelity uses all 19 displayed terms, all 16 example columns, and all 210 terms in `oeis-b268943.txt`, downloaded from <https://oeis.org/A268943/b268943.txt>. The sequence is by R. H. Hardin; the two formulas are conjectures of Colin Barker in [OEIS A268943](https://oeis.org/A268943), revision 8 (2019-01-17 04:20:36 EST). Source data attribution: OEIS Foundation Inc., OEIS content under [CC BY-SA 4.0](https://oeis.org/LICENSE). The source fixture is unchanged. The language records overlapping adjacent equal pairs and requires consecutive event labels to advance by one modulo the alphabet size.

The checker also distinguishes Fried's unequal-event-label condition from the target cyclic-successor condition: at alphabet size 8 and length 4 the counts are 3976 and 3928, respectively. Source: Sela Fried, *Proofs of several conjectures from the OEIS*, JIS 28 (2025), 25.4.3, subsection “Restricted repeating letters”.
