# Layman's A082630 insertion recurrence

Exact finite audits of the simultaneous insertion rule at threshold 7/11.
Run with Python 3 and its standard library:

```sh
python3 docs/reports/layman-a082630/check.py
```

To reproduce the optional 1,000-term table comparison in `output.txt`:

```sh
curl -L https://oeis.org/A082630/b082630.txt -o /tmp/b082630.txt
python3 docs/reports/layman-a082630/check.py --bfile /tmp/b082630.txt
```

The script reproduces OEIS revision 68's complete rows 0–4, its stated
prefixes in rows 5–6, and all 35 displayed DATA entries. It compares the
simultaneous rows 0–19 on their certified stable prefixes and compares four
1,000-term constructions: the source generator, ordered terminal-gap traversal,
the floor-successor rule, and the proposed linear recurrence. It checks 998
alternating coefficients, 996 recurrence equations, exact rational interval
endpoints, and the 2-by-2 matrix identity. All threshold and mathematical
comparisons use integers or `fractions.Fraction`.

These are finite checks, not a proof of infinite stabilization or of the
all-index recurrence. `output.txt` is the program's successful output with the
optional b-file comparison; the table's SHA-256 is included in that output.

Sources: [OEIS A082630](https://oeis.org/A082630), revision 68 (NAME, DATA,
EXAMPLE, and Zhuorui He's generator);
[Layman's corrected note](https://oeis.org/A085376/a085376.txt);
and the general successor-leaf/floor argument in
[NiuVector/A085260-proof](https://github.com/NiuVector/A085260-proof/blob/bed0c63a10c1250e7cbc86e7c7bfdb24be361d80/proof_note_en.md).
The displayed source terms/examples and generator are attributed above;
OEIS content is supplied under the [OEIS end-user license](https://oeis.org/LICENSE).
No priority is claimed for these constructions.
