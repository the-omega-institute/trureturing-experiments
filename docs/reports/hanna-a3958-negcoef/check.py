#!/usr/bin/env python3
"""Exact finite checks for Hanna's A395835--A395838 infinitude family.

Source data: Paul D. Hanna, https://oeis.org/A395835 through A395838,
b-files retrieved 2026-10-25 (OEIS data CC BY-SA 4.0).
Finite checks are not a proof of the universal conjecture.
"""
from fractions import Fraction
from math import comb

# Frozen n=2..30 b-file prefixes, independently compared with coefficient recursion.
SOURCE = {395835: {2: 3, 3: 5, 4: 8, 5: 11, 6: 14, 7: 17, 8: 21, 9: 25, 10: 30, 11: 35, 12: 40, 13: 45, 14: 51, 15: 56, 16: 62, 17: 69, 18: 76, 19: 82, 20: 90, 21: 97, 22: 105, 23: 113, 24: 121, 25: 129, 26: 138, 27: 147, 28: 156, 29: 165, 30: 175}, 395836: {2: 1, 3: 3, 4: 382, 5: 15596, 6: 428751, 7: 1390598, 8: 1036726902, 9: 59595779523, 10: 374973153024631, 11: 390927908309636208, 12: 287315311734503241860, 13: 144479928997809087241170, 14: 783078130263953707885969214, 15: 86783826953710133179485700584, 16: 28300905706048462234112872097712, 17: 16491215013305141268915597516911007570, 18: 339716226362949998409129401707477318094330, 19: 62337937889673576811173913688605693804569880, 20: 55762111768900628413806628750323494167905435045410, 21: 389789723004835651368435416849276791305771097162509204, 22: 36539019129605977234569260343426550741354968338266780767480, 23: 2201554016572939288518687999681296945824547204096199754737261740, 24: 96372165427279184919473080907761110213781934985652723155293567209608, 25: 1927431814591928575381630176776381258010784556686191936687823009961848986, 26: 1040581724997178027006373754709215396525852527687530228012146938746011669843750, 27: 272185504764008739464581244248566310613249098487536745968938677178713934014520758636, 28: 47438177214334132379328683775784511693903725392589310758186943016446713907514886440089592, 29: 2367672712826989601118111188386705595846338640071564164715115525730619898475853971952053004770, 30: 9524960795164817939147713960124781395748183665347464515052971213578350807798273752989715103947357918}, 395837: {2: 1, 3: 6, 4: 65, 5: 1302, 6: 85425, 7: 7316310, 8: 978911960, 9: 196901769732, 10: 55193621225743, 11: 11189954978322600, 12: 10087627980546480826, 13: 26791549795221415615580, 14: 1343026219489938132656481, 15: 134516049370090208389365126747, 16: 496520746161771547927510217314512, 17: 2488082740706962234924700258811769650, 18: 2769667610250267949025753741296670821845, 19: 275563776057409155724789072615391313452353880, 20: 1881275716339907064792569917172704469468570102025, 21: 117229022718015377572671635027973233055126684937184940, 22: 2679684544652781081111275377114746142896690672676550693512, 23: 103174798733167100003310601245637165197613800708588723927006672, 24: 12227333941089181055419680520178562235150004866471783897562399069507, 25: 1544874406068066662277740626830576821161241275452601101225992261234306310, 26: 197654417135483168417851414415276624803364250992655595289491259586075591258420, 27: 35700231795010493292822030647347467625022059629761058658372041031650811193134507992, 28: 11721691064005751146114152651738619221665373015335763338116859277119958417635130058956378, 29: 4782744227213109695025645122916903774968480906182357887740328154866388255586206099125623528880, 30: 2569512867837604602750843140242817081163977832347558827050035777751083256272586890492526350628728299}, 395838: {2: 3, 3: 2, 4: 3, 5: 3, 6: 3, 7: 3, 8: 4, 9: 4, 10: 5, 11: 5, 12: 5, 13: 5, 14: 6, 15: 5, 16: 6, 17: 7, 18: 7, 19: 6, 20: 8, 21: 7, 22: 8, 23: 8, 24: 8, 25: 8, 26: 9, 27: 9, 28: 9, 29: 9, 30: 10}}
EXAMPLES = {
    2: [1,1,1,-1,-10,-33,-55,65,842],
    3: [1,1,2,4,6,-3,-82,-481,-2098],
    4: [1,1,2,5,13,32,66,65,-382,-3553,-20254],
    5: [1,1,2,5,14,41,120,338,870,1802,1302,-15596,-135056],
    6: [1,1,2,5,14,42,131,415,1310,4046,11951,32452,73428,85425,-428751,-4557495],
}

def coefficients_by_equation(n, degree):
    """c_m=[x^(m-1)]F^2-[x^(m-n)]F^(2n); both use earlier c's."""
    assert n >= 2
    powers = [[0]*(degree+1) for _ in range(2*n+1)]
    for row in powers:
        row[0] = 1
    for m in range(1, degree+1):
        powers[1][m] = powers[2][m-1] - (powers[2*n][m-n] if m >= n else 0)
        for k in range(2, 2*n+1):
            powers[k][m] = sum(powers[1][j]*powers[k-1][m-j] for j in range(m+1))
    return powers[1]

def coefficient_by_lagrange(n, m):
    # c_m = [t^m](1+t-t^n)^(2m+1)/(2m+1).
    numerator = sum((-1)**j*comb(2*m+1,j)*comb(2*m+1-j,m-n*j)
                    for j in range(m//n+1))
    value = Fraction(numerator, 2*m+1)
    assert value.denominator == 1, (n,m,value)
    return value.numerator

def main():
    for n, expected in EXAMPLES.items():
        assert coefficients_by_equation(n,len(expected)-1) == expected
    print('FIDELITY: all coefficients displayed in A395835 EXAMPLE, n=2..6, agree exactly')
    positions = {}
    checked_coefficients = 0
    for n in range(2,31):
        cutoff = SOURCE[395835][n]  # A test range, not a proposed uniform bound.
        c = coefficients_by_equation(n,cutoff)
        for m, cm in enumerate(c):
            assert cm == coefficient_by_lagrange(n,m), (n,m)
            checked_coefficients += 1
        first = next(m for m,cm in enumerate(c) if cm < 0)
        assert first == cutoff
        assert -c[first] == SOURCE[395836][n]
        assert c[first-1] == SOURCE[395837][n]
        positions[n] = first
        difference = first if n == 2 else first-positions[n-1]
        assert difference == SOURCE[395838][n]
        print(f'n={n} M={first} c_M={c[first]} c_previous={c[first-1]} difference_table={difference}')
    print(f'EXACT: recursion and Lagrange agree at {checked_coefficients} coefficients; all four b-files agree for n=2..30')
    # The quotient rule gives q'=(phi-2t phi')/phi^3.
    # Exact polynomial identity (no sampled identity or floating point).
    for n in range(2,31):
        phi = [1,1]+[0]*(n-2)+[-1]
        numerator = [(1-2*k)*a for k,a in enumerate(phi)]
        expected = [1,-1]+[0]*(n-2)+[2*n-1]
        assert numerator == expected
        for j in range(1,33):
            t = Fraction(j,32)
            v = 1+t-t**n
            assert v >= 1
            charge = 1-t+(2*n-1)*t**n
            assert charge > 0
            assert charge/v**3 > 0
    print('EXACT: q derivative numerator identity for n=2..30; positive on the rational grid t=j/32, j=1..32')
    print('SOURCE INDEXING: A395838 table seeds a(2)=M(2)=3, then a(n)=M(n)-M(n-1); A395835 FORMULA (3) instead shifts forward')
    print('SCOPE: finite source fidelity and route checks only; no uniform existence proof or asymptotic bound is certified')
    print('ALL_CHECKS_OK')

if __name__ == '__main__':
    main()
