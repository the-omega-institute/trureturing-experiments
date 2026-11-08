# Hardin–Plouffe OEIS A300344

Source: [OEIS A300344](https://oeis.org/A300344), revision 12, and its [history](https://oeis.org/history?seq=A300344). R. H. Hardin defines the binary array count; Simon Plouffe supplies the empirical generating function. The first 29 published terms and ten example boards are used for fidelity checks. The parent [A300350](https://oeis.org/A300350) labels its width-two recurrence empirical.

Run `python3 docs/reports/hardin-a300344-ladder/check.py` from the experiments repository root. Python 3's standard library suffices. The program uses only exact integers, prints `ALL_CHECKS_OK` on success, and exits zero. Its output is in [check.out](check.out).

Literal cell-neighbor enumeration for n=1 through 9 gives 1, 2, 6, 15, 37, 90, 223, 550, 1355. Every board in that enumeration is also checked against an independently written row predicate. The 16-state transfer model and the eight-state simultaneous-complement quotient agree through n=210. The program derives and checks the quotient matrix, initial and terminal vectors, six continuation vectors and the order-five annihilator. Both models reproduce all 29 displayed OEIS terms. The recurrence is checked for n=6 through 210 and the cleared generating-function identity through coefficient 210.

The two diagonal directions have the same counts, including all nine literal counts above. Column reflection followed by complementation when needed to restore the upper-left zero gives a bijection; the program checks its action on every valid enumerated board. The main-diagonal rule accepts only two of the ten published example boards, whereas the antidiagonal rule accepts all ten. The examples and the word “antidiagonally” fix the intended rule; the counts do not distinguish the directions.

These checks support a proposed all-length argument via the native counting bridge and the continuation identity. They are finite computational evidence and do not constitute a Lean proof or a claim that literature discovery is complete.
