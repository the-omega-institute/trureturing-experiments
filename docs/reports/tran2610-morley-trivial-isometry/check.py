#!/usr/bin/env python3
"""Exact independent checks for Tran, arXiv:2610.01216v1, section 12(4).

Only integers, Fraction intervals, and exact symbolic polynomial identities are
used. No numerical floating-point root or eigenvalue is evidence. The explicit
Latin square below was generated independently; no ancillary source was used.
"""
from fractions import Fraction as Q
from itertools import combinations
import hashlib
import json
import sympy as sp

LATIN = [
    [22, 15, 24, 39, 38, 21, 17, 13, 43, 44, 48, 31, 0, 2, 3, 12, 28, 26, 29, 4, 1, 46, 14, 7, 35, 45, 9, 30, 19, 40, 33, 49, 27, 42, 37, 8, 11, 25, 20, 5, 34, 16, 47, 23, 36, 32, 6, 41, 10, 18],
    [38, 49, 43, 13, 0, 23, 9, 15, 20, 34, 47, 18, 37, 25, 14, 16, 19, 39, 7, 40, 10, 27, 1, 41, 46, 4, 21, 42, 24, 28, 48, 5, 6, 22, 29, 3, 30, 12, 17, 32, 8, 11, 35, 26, 33, 45, 44, 31, 2, 36],
    [1, 13, 11, 2, 35, 7, 18, 21, 33, 23, 36, 12, 34, 22, 4, 15, 24, 8, 26, 5, 46, 9, 28, 3, 42, 29, 0, 47, 37, 20, 19, 6, 25, 27, 39, 43, 45, 32, 41, 30, 17, 10, 16, 31, 38, 44, 49, 40, 48, 14],
    [24, 41, 26, 29, 20, 8, 33, 22, 0, 16, 4, 7, 47, 3, 6, 45, 44, 43, 11, 46, 14, 15, 18, 9, 39, 28, 48, 2, 21, 31, 1, 12, 32, 13, 40, 42, 23, 37, 34, 35, 30, 19, 38, 36, 10, 49, 5, 27, 25, 17],
    [15, 26, 49, 43, 2, 18, 27, 44, 31, 17, 11, 23, 5, 30, 47, 1, 29, 33, 8, 9, 38, 7, 45, 37, 48, 13, 32, 20, 16, 14, 46, 35, 21, 25, 41, 12, 34, 42, 40, 28, 39, 3, 0, 10, 24, 19, 22, 4, 36, 6],
    [37, 25, 28, 38, 9, 22, 10, 30, 24, 20, 19, 36, 35, 7, 26, 42, 41, 13, 31, 44, 4, 21, 39, 11, 45, 43, 5, 12, 48, 8, 17, 2, 15, 3, 47, 23, 49, 6, 16, 18, 14, 1, 40, 33, 27, 46, 0, 29, 32, 34],
    [32, 34, 17, 20, 15, 42, 2, 37, 27, 1, 24, 25, 31, 36, 16, 9, 45, 48, 21, 26, 35, 12, 6, 40, 3, 8, 23, 5, 39, 33, 43, 4, 44, 14, 0, 10, 28, 11, 46, 38, 7, 49, 41, 22, 29, 30, 13, 19, 18, 47],
    [9, 35, 42, 5, 14, 49, 3, 43, 38, 36, 12, 10, 6, 39, 31, 23, 47, 40, 4, 32, 44, 1, 20, 29, 21, 16, 17, 45, 28, 41, 26, 11, 7, 33, 13, 37, 46, 2, 27, 0, 48, 30, 18, 25, 34, 22, 24, 15, 19, 8],
    [18, 6, 23, 48, 39, 16, 42, 32, 26, 22, 37, 41, 46, 49, 35, 13, 34, 20, 10, 12, 33, 4, 29, 5, 15, 19, 43, 11, 30, 9, 31, 21, 36, 40, 2, 14, 3, 47, 0, 27, 1, 8, 45, 24, 25, 38, 7, 17, 28, 44],
    [5, 7, 4, 49, 29, 2, 47, 6, 15, 38, 1, 20, 39, 11, 43, 40, 31, 28, 19, 27, 32, 42, 30, 45, 8, 22, 3, 26, 0, 48, 12, 36, 9, 35, 34, 16, 14, 44, 10, 25, 13, 41, 24, 18, 37, 23, 33, 46, 17, 21],
    [11, 32, 10, 4, 31, 6, 37, 45, 48, 7, 41, 30, 16, 23, 22, 28, 0, 12, 2, 19, 34, 24, 25, 18, 29, 42, 35, 33, 47, 3, 9, 39, 46, 20, 1, 49, 27, 17, 21, 14, 40, 44, 36, 38, 13, 5, 26, 8, 15, 43],
    [16, 39, 13, 41, 46, 10, 32, 34, 17, 5, 25, 22, 29, 0, 11, 35, 23, 45, 33, 43, 21, 30, 36, 1, 40, 14, 15, 9, 20, 38, 18, 3, 28, 26, 8, 4, 37, 24, 49, 19, 42, 2, 12, 44, 48, 47, 27, 6, 31, 7],
    [31, 11, 1, 26, 41, 14, 5, 4, 42, 39, 30, 35, 45, 44, 37, 21, 36, 7, 43, 18, 28, 40, 9, 24, 34, 49, 33, 48, 46, 32, 23, 0, 12, 8, 15, 47, 38, 3, 22, 16, 27, 13, 19, 6, 2, 17, 10, 20, 29, 25],
    [20, 12, 32, 21, 47, 43, 36, 19, 5, 3, 46, 2, 25, 34, 13, 49, 39, 0, 23, 30, 15, 22, 37, 14, 10, 38, 8, 27, 41, 18, 6, 31, 1, 17, 33, 35, 48, 26, 4, 42, 45, 7, 9, 40, 44, 29, 16, 11, 24, 28],
    [27, 31, 14, 24, 32, 4, 23, 7, 44, 33, 39, 26, 21, 35, 38, 0, 13, 10, 12, 45, 18, 3, 8, 46, 20, 6, 49, 40, 15, 5, 22, 42, 16, 34, 19, 17, 41, 28, 48, 11, 25, 47, 37, 43, 30, 9, 36, 2, 1, 29],
    [14, 4, 46, 3, 8, 48, 40, 18, 29, 12, 44, 1, 42, 9, 23, 10, 35, 15, 0, 24, 17, 20, 41, 39, 19, 32, 6, 36, 49, 37, 11, 22, 31, 30, 27, 38, 7, 16, 25, 47, 26, 21, 2, 34, 43, 33, 28, 45, 5, 13],
    [29, 38, 2, 35, 37, 30, 16, 5, 19, 21, 42, 49, 18, 20, 41, 17, 15, 24, 3, 13, 45, 0, 4, 22, 47, 1, 12, 46, 6, 43, 44, 32, 11, 39, 25, 31, 10, 40, 36, 48, 9, 28, 14, 27, 7, 8, 34, 33, 26, 23],
    [41, 29, 39, 11, 7, 15, 46, 27, 21, 31, 38, 3, 40, 17, 32, 33, 48, 19, 30, 49, 47, 25, 10, 16, 2, 9, 37, 13, 1, 23, 35, 8, 45, 18, 22, 44, 26, 0, 5, 20, 24, 6, 34, 28, 12, 14, 4, 36, 43, 42],
    [44, 9, 20, 19, 40, 27, 49, 1, 13, 26, 8, 17, 7, 48, 34, 31, 18, 36, 16, 25, 30, 32, 23, 43, 5, 21, 38, 37, 11, 45, 3, 41, 47, 28, 42, 6, 4, 10, 24, 29, 12, 46, 22, 2, 15, 0, 14, 39, 35, 33],
    [21, 24, 3, 16, 25, 40, 28, 10, 47, 0, 32, 39, 33, 45, 2, 14, 46, 41, 27, 42, 26, 35, 31, 17, 43, 34, 4, 19, 8, 13, 38, 48, 29, 6, 11, 5, 22, 30, 23, 37, 18, 36, 7, 49, 1, 15, 20, 9, 44, 12],
    [35, 23, 37, 6, 30, 20, 1, 31, 18, 46, 49, 33, 4, 24, 27, 38, 8, 42, 25, 34, 22, 5, 26, 12, 11, 2, 19, 21, 13, 0, 7, 10, 14, 29, 43, 48, 40, 36, 47, 3, 16, 15, 44, 32, 45, 39, 17, 28, 9, 41],
    [39, 14, 22, 18, 27, 24, 34, 11, 46, 30, 10, 45, 41, 29, 0, 8, 6, 32, 1, 47, 48, 28, 17, 19, 9, 12, 25, 7, 40, 42, 16, 37, 5, 49, 44, 26, 36, 33, 43, 2, 38, 4, 23, 15, 35, 3, 31, 21, 13, 20],
    [30, 48, 33, 47, 21, 29, 7, 8, 6, 49, 5, 24, 3, 38, 25, 44, 32, 1, 42, 14, 37, 26, 22, 15, 36, 20, 2, 35, 43, 39, 10, 13, 0, 45, 4, 41, 18, 34, 11, 46, 28, 31, 17, 9, 19, 12, 23, 16, 27, 40],
    [33, 16, 6, 40, 28, 45, 13, 2, 8, 32, 17, 27, 20, 14, 46, 4, 11, 29, 34, 36, 23, 44, 12, 38, 30, 15, 26, 10, 3, 21, 39, 43, 41, 24, 9, 22, 19, 18, 7, 49, 35, 25, 5, 48, 31, 37, 42, 1, 47, 0],
    [49, 40, 48, 15, 36, 19, 22, 47, 7, 2, 29, 4, 27, 42, 45, 11, 43, 37, 9, 3, 5, 16, 34, 31, 24, 30, 18, 1, 10, 26, 25, 14, 33, 41, 12, 20, 32, 39, 35, 6, 23, 17, 28, 0, 8, 13, 21, 44, 46, 38],
    [2, 37, 19, 8, 4, 1, 20, 48, 14, 41, 43, 28, 22, 12, 36, 46, 10, 47, 24, 35, 31, 23, 42, 34, 26, 17, 7, 6, 45, 25, 29, 33, 38, 44, 3, 18, 21, 49, 39, 13, 15, 5, 32, 30, 9, 11, 40, 0, 16, 27],
    [34, 19, 5, 10, 23, 9, 30, 16, 41, 43, 33, 21, 13, 26, 18, 47, 3, 31, 37, 0, 6, 2, 38, 28, 25, 7, 14, 8, 22, 15, 36, 44, 49, 4, 17, 39, 35, 45, 1, 40, 46, 24, 20, 42, 32, 27, 29, 48, 12, 11],
    [23, 5, 25, 1, 10, 35, 8, 12, 2, 48, 27, 46, 49, 40, 7, 43, 17, 21, 28, 31, 11, 29, 47, 6, 4, 39, 44, 38, 14, 22, 42, 45, 30, 9, 36, 34, 16, 41, 13, 26, 19, 0, 3, 37, 18, 20, 32, 24, 33, 15],
    [46, 42, 21, 9, 49, 12, 48, 17, 3, 15, 13, 29, 10, 28, 39, 22, 33, 14, 40, 6, 43, 45, 27, 47, 7, 23, 41, 16, 36, 30, 24, 19, 20, 0, 38, 25, 44, 35, 32, 8, 2, 18, 26, 1, 5, 31, 11, 37, 34, 4],
    [28, 2, 41, 33, 19, 36, 4, 20, 22, 27, 0, 42, 8, 32, 17, 29, 40, 11, 35, 1, 12, 31, 21, 23, 44, 5, 24, 25, 7, 34, 13, 18, 26, 10, 14, 45, 39, 46, 38, 9, 37, 48, 15, 47, 16, 43, 30, 3, 6, 49],
    [10, 27, 40, 32, 44, 38, 45, 36, 35, 42, 7, 13, 1, 47, 20, 26, 21, 22, 6, 11, 8, 19, 24, 2, 12, 0, 34, 28, 17, 46, 49, 15, 3, 37, 31, 33, 9, 43, 30, 23, 5, 29, 4, 39, 14, 48, 25, 18, 41, 16],
    [45, 43, 7, 36, 11, 37, 35, 38, 4, 14, 40, 32, 9, 31, 24, 39, 30, 16, 20, 41, 25, 18, 48, 49, 27, 33, 46, 44, 12, 10, 34, 29, 42, 2, 26, 19, 47, 8, 15, 17, 6, 23, 21, 5, 3, 28, 1, 13, 0, 22],
    [48, 44, 35, 7, 12, 46, 39, 49, 9, 47, 28, 37, 11, 13, 19, 6, 14, 18, 15, 20, 29, 17, 40, 36, 16, 41, 10, 23, 26, 2, 0, 27, 8, 32, 5, 1, 24, 22, 33, 34, 31, 42, 30, 45, 4, 25, 43, 38, 21, 3],
    [0, 46, 9, 31, 17, 32, 41, 23, 28, 11, 35, 48, 26, 21, 44, 2, 38, 6, 47, 8, 27, 39, 13, 25, 22, 40, 16, 24, 18, 4, 30, 7, 34, 1, 45, 15, 5, 29, 37, 36, 10, 43, 33, 20, 49, 42, 12, 14, 3, 19],
    [13, 28, 34, 25, 42, 33, 31, 39, 45, 19, 3, 11, 44, 8, 10, 48, 9, 17, 49, 16, 2, 38, 15, 30, 37, 27, 22, 29, 35, 1, 20, 23, 24, 46, 6, 32, 12, 5, 18, 7, 0, 14, 43, 4, 21, 36, 41, 47, 40, 26],
    [25, 8, 16, 34, 33, 31, 26, 35, 40, 24, 45, 15, 23, 18, 30, 19, 12, 3, 36, 37, 9, 10, 2, 48, 28, 47, 27, 0, 29, 44, 32, 17, 39, 38, 49, 11, 43, 7, 6, 4, 21, 22, 13, 14, 41, 1, 46, 42, 20, 5],
    [26, 21, 38, 30, 3, 47, 6, 24, 11, 8, 18, 19, 43, 27, 12, 25, 16, 34, 48, 17, 42, 41, 5, 20, 0, 31, 36, 32, 33, 49, 14, 1, 2, 23, 35, 46, 13, 4, 29, 15, 44, 37, 10, 7, 28, 40, 45, 22, 39, 9],
    [42, 36, 12, 17, 26, 44, 43, 14, 37, 45, 23, 0, 28, 41, 1, 32, 5, 30, 13, 39, 24, 33, 35, 8, 31, 10, 47, 22, 25, 19, 15, 46, 40, 48, 16, 7, 29, 20, 9, 21, 3, 38, 27, 11, 6, 4, 18, 34, 49, 2],
    [12, 30, 31, 23, 16, 41, 19, 9, 36, 13, 15, 40, 48, 43, 29, 34, 4, 46, 44, 21, 3, 37, 33, 0, 18, 25, 39, 49, 27, 6, 5, 28, 35, 47, 24, 2, 8, 14, 45, 22, 32, 26, 11, 17, 20, 10, 38, 7, 42, 1],
    [4, 10, 36, 44, 6, 34, 25, 26, 23, 28, 14, 9, 17, 19, 42, 37, 27, 35, 5, 7, 49, 47, 3, 21, 33, 24, 1, 39, 2, 11, 41, 30, 22, 12, 18, 13, 0, 38, 8, 43, 20, 40, 31, 29, 46, 16, 15, 32, 45, 48],
    [36, 22, 27, 14, 43, 25, 12, 33, 39, 10, 6, 34, 30, 16, 49, 41, 1, 5, 45, 38, 40, 8, 7, 13, 17, 48, 31, 15, 32, 47, 4, 26, 18, 19, 46, 28, 20, 9, 44, 24, 29, 35, 42, 3, 0, 21, 2, 23, 11, 37],
    [7, 20, 15, 0, 34, 39, 11, 42, 12, 4, 9, 8, 19, 1, 33, 3, 25, 27, 38, 2, 16, 43, 44, 26, 49, 18, 29, 41, 23, 36, 37, 40, 17, 5, 21, 24, 31, 13, 28, 10, 47, 45, 48, 46, 22, 6, 35, 30, 14, 32],
    [3, 17, 30, 42, 48, 13, 14, 46, 25, 18, 26, 43, 2, 37, 5, 27, 22, 44, 39, 10, 7, 6, 0, 33, 41, 36, 20, 34, 38, 24, 40, 16, 23, 11, 28, 21, 1, 15, 19, 45, 4, 32, 29, 12, 47, 35, 9, 49, 8, 31],
    [43, 33, 29, 22, 13, 5, 21, 40, 32, 9, 2, 44, 15, 6, 8, 36, 37, 25, 18, 23, 0, 49, 19, 27, 14, 11, 45, 3, 4, 35, 28, 24, 10, 7, 20, 30, 42, 48, 31, 39, 41, 12, 1, 16, 17, 34, 47, 26, 38, 46],
    [19, 45, 0, 28, 1, 3, 24, 25, 34, 35, 16, 5, 14, 4, 9, 18, 7, 49, 17, 15, 36, 48, 11, 32, 23, 37, 30, 31, 44, 29, 47, 38, 13, 43, 10, 40, 6, 27, 2, 41, 33, 20, 46, 21, 42, 26, 8, 12, 22, 39],
    [47, 1, 44, 45, 22, 0, 38, 28, 10, 29, 31, 16, 24, 46, 21, 20, 49, 9, 41, 48, 13, 14, 32, 42, 6, 3, 40, 18, 5, 7, 27, 34, 4, 15, 30, 36, 17, 23, 26, 12, 11, 33, 25, 8, 39, 2, 19, 43, 37, 35],
    [6, 3, 47, 37, 18, 17, 0, 29, 30, 25, 22, 14, 36, 10, 15, 5, 2, 23, 32, 28, 19, 11, 46, 44, 38, 26, 42, 4, 9, 16, 45, 20, 43, 21, 48, 27, 33, 1, 12, 31, 49, 34, 8, 13, 40, 41, 39, 35, 7, 24],
    [8, 47, 18, 12, 45, 28, 15, 0, 49, 40, 34, 6, 32, 5, 48, 7, 42, 38, 46, 33, 20, 13, 16, 4, 1, 35, 11, 43, 31, 17, 21, 9, 37, 36, 23, 29, 2, 19, 14, 44, 22, 27, 39, 41, 26, 24, 3, 25, 30, 10],
    [40, 0, 8, 46, 24, 26, 29, 3, 1, 37, 21, 38, 12, 15, 28, 30, 20, 4, 14, 22, 41, 36, 49, 10, 32, 44, 13, 17, 34, 27, 2, 47, 19, 16, 7, 9, 25, 31, 42, 33, 43, 39, 6, 35, 11, 18, 48, 5, 23, 45],
    [17, 18, 45, 27, 5, 11, 44, 41, 16, 6, 20, 47, 38, 33, 40, 24, 26, 2, 22, 29, 39, 34, 43, 35, 13, 46, 28, 14, 42, 12, 8, 25, 48, 31, 32, 0, 15, 21, 3, 1, 36, 9, 49, 19, 23, 7, 37, 10, 4, 30],
]


class I:
    """Closed rational interval, with exact outward enclosure operations."""
    def __init__(self, lo, hi=None):
        self.lo = Q(lo)
        self.hi = Q(lo if hi is None else hi)
        assert self.lo <= self.hi

    def __add__(self, other):
        other = interval(other)
        return I(self.lo + other.lo, self.hi + other.hi)

    __radd__ = __add__

    def __neg__(self):
        return I(-self.hi, -self.lo)

    def __sub__(self, other):
        return self + -interval(other)

    def __rsub__(self, other):
        return interval(other) + -self

    def __mul__(self, other):
        other = interval(other)
        vals = [a*b for a in (self.lo, self.hi) for b in (other.lo, other.hi)]
        return I(min(vals), max(vals))

    __rmul__ = __mul__

    def __truediv__(self, other):
        other = interval(other)
        assert other.lo > 0 or other.hi < 0
        return self * I(1/other.hi, 1/other.lo)

    def __rtruediv__(self, other):
        return interval(other) / self

    def __pow__(self, power):
        assert isinstance(power, int) and power >= 0
        if power == 0:
            return I(1)
        if power % 2:
            return I(self.lo**power, self.hi**power)
        return I(0 if self.lo <= 0 <= self.hi else min(self.lo**power, self.hi**power),
                 max(self.lo**power, self.hi**power))

    def contains_zero(self):
        return self.lo <= 0 <= self.hi

    def width(self):
        return self.hi - self.lo

    def __str__(self):
        # Decimal endpoints are terminating RATIONALS, rounded outwards.
        scale = 10**36
        lo = (self.lo*scale).__floor__()
        hi = (self.hi*scale).__ceil__()
        def decimal(x):
            sign = '-' if x < 0 else ''
            a, b = divmod(abs(x), scale)
            return f'{sign}{a}.{b:036d}'
        return '[' + decimal(lo) + ', ' + decimal(hi) + ']'


def interval(x):
    return x if isinstance(x, I) else I(x)


def phi(x):
    return (3*x - x**3)/2


def inverse_phi(b, steps=185):
    """Unique root in (1,2), since phi strictly decreases there from 1 to -1."""
    assert -1 < b < 1
    lo, hi = Q(1), Q(2)
    for _ in range(steps):
        mid = (lo + hi)/2
        if phi(mid) > b:
            lo = mid
        else:
            hi = mid
    assert phi(lo) >= b >= phi(hi)
    return I(lo, hi)


def inverse_interval(b):
    b = interval(b)
    return I(inverse_phi(b.hi).lo, inverse_phi(b.lo).hi)


def cycle_type(p):
    assert sorted(p) == list(range(len(p)))
    seen = set()
    parts = []
    for i in range(len(p)):
        if i in seen:
            continue
        j, length = i, 0
        while j not in seen:
            seen.add(j)
            j = p[j]
            length += 1
        parts.append(length)
    return tuple(sorted(parts))


def graph_certificate():
    q = len(LATIN)
    assert q == 50 and q > 4
    assert all(sorted(row) == list(range(q)) for row in LATIN)
    assert all(sorted(LATIN[r][c] for r in range(q)) == list(range(q)) for c in range(q))
    N = q*q
    rows = [[r*q+c for c in range(q)] for r in range(q)]
    cols = [[r*q+c for r in range(q)] for c in range(q)]
    syms = [[r*q+c for r in range(q) for c in range(q) if LATIN[r][c] == s]
            for s in range(q)]
    lines = rows + cols + syms
    masks = [sum(1 << v for v in line) for line in lines]
    closed = [0]*N
    for line, mask in zip(lines, masks):
        for v in line:
            closed[v] |= mask
    adj = [mask ^ (1 << i) for i, mask in enumerate(closed)]
    degree, lam, mu = 3*(q-1), q, 6
    assert all(mask.bit_count() == degree for mask in adj)
    pairs = 0
    for i, ai in enumerate(adj):
        for j in range(i):
            assert (ai & adj[j]).bit_count() == (lam if (ai >> j) & 1 else mu)
            pairs += 1
    print(f'Latin square q={q}, N={N}, SRG parameters=({N},{degree},{lam},{mu})')
    print(f'Exhaustive SRG common-neighbour pairs: {pairs}')
    # Any 5-clique has, at a fixed vertex, two of its four neighbours
    # sharing the same one of the three line types. These three vertices
    # lie on one line. The following exhaustive check proves every such
    # triple has precisely that line as its common CLOSED neighbourhood.
    # Consequently every clique of size >=5 lies on a line, and the
    # maximum q-cliques are precisely the 3q displayed lines.
    triples = 0
    for line, mask in zip(lines, masks):
        for a, b, c in combinations(line, 3):
            assert closed[a] & closed[b] & closed[c] == mask
            triples += 1
    # Disjoint maximum cliques are exactly the pairs in one parallel class.
    for i in range(3*q):
        for j in range(i):
            assert (masks[i] & masks[j]).bit_count() == (0 if i//q == j//q else 1)
    classes = [LATIN,
               [[LATIN[r][c] for r in range(q)] for c in range(q)],
               [[next(c for c in range(q) if LATIN[r][c] == s) for r in range(q)]
                for s in range(q)]]
    signatures = []
    for perms in classes:
        for i, f in enumerate(perms):
            signature = []
            for j, g in enumerate(perms):
                if i == j:
                    continue
                inv = [0]*q
                for x, y in enumerate(g):
                    inv[y] = x
                signature.append(cycle_type([inv[f[x]] for x in range(q)]))
            signatures.append(tuple(sorted(signature)))
    # Each partition gives half the lengths of the cycles in the induced
    # bipartite graph on two disjoint lines. It is intrinsic to the graph,
    # including when an automorphism interchanges the three line classes.
    assert len(set(signatures)) == 3*q
    # Exhaust all possible images of every maximum clique. Each is forced.
    candidates = [[j for j, t in enumerate(signatures) if t == s]
                  for s in signatures]
    assert candidates == [[i] for i in range(3*q)]
    # Every vertex is the intersection of its row and its column, so all
    # vertices are fixed. This certifies the FULL group, not an autotopy
    # subgroup preserving the named three coordinate classes.
    print(f'Exhaustive collinear triple tests: {triples}')
    print(f'Intrinsic maximum-clique signatures: {len(set(signatures))}/{3*q} distinct')
    print('All possible maximum-clique images exhausted; full graph automorphism group order = 1')
    print('Latin square sha256 =', hashlib.sha256(json.dumps(LATIN, separators=(',', ':')).encode()).hexdigest())
    return q


def latin_values(q, p):
    p = interval(p)
    e = Q(1, q)
    a = phi(p)
    b = -(e*e + 3*e*(1-e)*a)/((1-e)*(1-2*e))
    r = inverse_interval(b)
    U, B, C = p-r, 1+2*b-3*a, 1+2*r-3*p
    F = B*U*(U+2*e*C)-e*(a-b)*C**2
    hp, hm = q*(a-b)+B, B
    kp, km = q*U+C, C
    k0 = 1+3*(q-1)*p+(q-1)*(q-2)*r
    return locals()


def construction_certificate(q):
    # Exact symbolic identities for ALL parameters, before specializing.
    e, a, b, p, r = sp.symbols('e a b p r')
    U, B, C = p-r, 1+2*b-3*a, 1+2*r-3*p
    F = B*U*(U+2*e*C)-e*(a-b)*C**2
    assert sp.expand(F-e**2*(B*(U/e+C)**2-((a-b)/e+B)*C**2)) == 0
    bb = -(e**2+3*e*(1-e)*a)/((1-e)*(1-2*e))
    assert sp.cancel((1+3*(1/e-1)*a+(1/e-1)*(1/e-2)*b).subs(b, bb)) == 0
    qsym = sp.Symbol('q', nonzero=True)
    ee = 1/qsym
    bq = -(ee**2+3*ee*(1-ee)*a)/((1-ee)*(1-2*ee))
    bB, bC, bU = 1+2*bq-3*a, 1+2*r-3*p, p-r
    FF = bB*bU*(bU+2*ee*bC)-ee*(a-bq)*bC**2
    degree, d = 3*(qsym-1), 1-r
    # H = (1-b)I+(a-b)A+bJ.  Squaring KP in the exact SRG
    # algebra gives these three coefficients.  The residual identity
    # certifies every matrix entry simultaneously, for arbitrary q.
    residual_coefficients = [
        (1-bq)*bC**2-bB*(d**2+(degree-6)*bU**2),
        (a-bq)*bC**2-bB*(2*d*bU+(qsym-6)*bU**2),
        bq*bC**2-bB*(6*bU**2-(d+degree*bU)**2/qsym**2)]
    for coefficient, multiplier in zip(residual_coefficients, [-3*qsym,-qsym,3]):
        assert sp.cancel(coefficient-multiplier*FF) == 0
    print('Symbolic identities: H*1=0 and F=e^2*(h_minus*k_plus^2-h_plus*k_minus^2): exact')
    print('Every-entry identity: C^2*H-B*K*P*K = F*(-3*q*I-q*A+3*J): exact')
    lo, hi = Q(8,5), Q(13,8)
    fl, fh = latin_values(q,lo)['F'], latin_values(q,hi)['F']
    assert fl.hi < 0 < fh.lo
    print('Initial IVT bracket p = [8/5,13/8]')
    print('F(left) =', fl)
    print('F(right) =', fh)
    for _ in range(130):
        mid = (lo+hi)/2
        fm = latin_values(q,mid)['F']
        if fm.hi < 0:
            lo = mid
        elif fm.lo > 0:
            hi = mid
        else:
            raise AssertionError('insufficient inner root precision')
    z = latin_values(q,I(lo,hi))
    p, r, a, b = (z[k] for k in ('p','r','a','b'))
    assert 1 < p.lo < p.hi < r.lo < r.hi < 2
    assert z['C'].hi < Q(-1,5) and z['U'].hi < 0
    assert z['hp'].lo > 0 and z['hm'].lo > 0
    assert z['kp'].hi < 0 and z['km'].hi < 0 and z['k0'].lo > 0
    assert a.lo > b.hi
    tau = z['hm']/z['km']**2
    assert tau.lo > 0
    # At an actual IVT zero the exact symbolic identity above makes the
    # second nonconstant eigenspace equation zero, not merely small.
    residuals = [z['hm']-tau*z['km']**2,
                 z['hp']-tau*z['kp']**2,
                 phi(r)-b,
                 1+3*(q-1)*a+(q-1)*(q-2)*b]
    assert all(x.contains_zero() and x.width() < Q(1,10**32) for x in residuals)
    print('Certified p =', p)
    print('Certified r =', r)
    print('Normal Gram diagonal = 1 exactly; adjacent entry a =', a)
    print('Normal Gram nonadjacent entry b =', b)
    print('H nonconstant eigenvalues h_plus =', z['hp'], '; h_minus =', z['hm'])
    print('K nonconstant eigenvalues k_plus =', z['kp'], '; k_minus =', z['km'])
    print('tau =', tau)
    print('All four Gram/compatibility interval residuals contain 0, width < 1/10^32')
    t = q*q/z['k0']
    edge_squared = 2*t**2/tau
    assert edge_squared.lo > 0
    # All edges of M_i=(t/sqrt(tau))*P*e_i have this same positive
    # squared length. Signed-distance equations reduce exactly to
    # 1+t*(K_ji-k0/N)=t*K_ji, since t*k0/N=1.
    assert sp.cancel(1+sp.Symbol('t')*(sp.Symbol('x')-sp.Symbol('k0')/sp.Symbol('N'))
                     -sp.Symbol('t')*sp.Symbol('x')).subs(
                         sp.Symbol('t'),sp.Symbol('N')/sp.Symbol('k0')) == 0
    print('Every squared Morley edge = 2*t^2/tau in', edge_squared)
    print(f'Exact IVT certificate: a nondegenerate {q*q-1}-simplex; Morley regular; Isom = 1')


def johnson_fidelity():
    # Independent reproduction of Theorem 10.2, the certified k=9 row
    # of the table: p in [1.6676,1.6680], with a positive-to-negative
    # spectral compatibility residual. Paper formulas only; no ancillary code.
    k = 9
    degree, nondegree = 2*(k-2), (k-2)*(k-3)//2
    def values(p):
        p = interval(p)
        a = phi(p)
        b = -(1+degree*a)/nondegree
        r = inverse_interval(b)
        h1, h2 = 1+(k-4)*a-(k-3)*b, 1-2*a+b
        k1, k2 = 1+(k-4)*p-(k-3)*r, 1-2*p+r
        return p,r,h1,h2,k1,k2,h1*k2**2-h2*k1**2
    lo,hi = Q(4169,2500),Q(417,250)
    left,right = values(lo)[-1],values(hi)[-1]
    z = values(I(lo,hi))
    assert left.lo > 0 > right.hi
    assert z[0].lo > 1 and z[0].hi < z[1].lo and z[1].hi < 2
    assert z[2].lo > 0 and z[3].lo > 0 and z[4].hi < 0 and z[5].hi < 0
    print('Fidelity: source Theorem 10.2, k=9, dimension 35, p in [4169/2500,417/250]')
    print('Johnson compatibility at left =', left)
    print('Johnson compatibility at right =', right)
    print('Johnson h1,h2 positive; kappa1,kappa2 negative on the full interval: certified')


def main():
    johnson_fidelity()
    q = graph_certificate()
    construction_certificate(q)
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
