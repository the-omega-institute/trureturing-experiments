#!/usr/bin/env python3
"""Certify the local shape of an actual future-record pressure profile.

Global future-record status is an external published input. This program
checks the strict inventories and the two derivative signs, without RH.
"""

import argparse
import importlib.util
import json
from fractions import Fraction as Q
from pathlib import Path


LOG_SOURCE = Path(__file__).resolve().parents[1] / "robin-density-transport-3e47/certificate.py"
spec = importlib.util.spec_from_file_location("robin_density_certificate", LOG_SOURCE)
log_module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(log_module)
log_enclosure = log_module.log_enclosure
SCALE = 10**50


def ceiling(value):
    return -((-value.numerator) // value.denominator)


class Interval:
    """Outward integer fixed-point arithmetic; endpoints are multiples of 1/SCALE."""

    def __init__(self, value):
        value = Q(value) * SCALE
        self.lo = value.numerator // value.denominator
        self.hi = ceiling(value)

    @classmethod
    def raw(cls, lo, hi):
        if lo > hi:
            raise ArithmeticError("reversed interval")
        result = object.__new__(cls)
        result.lo, result.hi = lo, hi
        return result

    def __add__(self, other):
        other = interval(other)
        return self.raw(self.lo + other.lo, self.hi + other.hi)

    __radd__ = __add__

    def __neg__(self):
        return self.raw(-self.hi, -self.lo)

    def __sub__(self, other):
        return self + -interval(other)

    def __rsub__(self, other):
        return interval(other) + -self

    def __mul__(self, other):
        other = interval(other)
        products = [x*y for x in (self.lo, self.hi) for y in (other.lo, other.hi)]
        return self.raw(min(products)//SCALE, -((-max(products))//SCALE))

    __rmul__ = __mul__

    def inverse(self):
        if self.lo <= 0 <= self.hi:
            raise ArithmeticError("division across zero")
        return self.raw(SCALE*SCALE//self.hi, -((-SCALE*SCALE)//self.lo))

    def __truediv__(self, other):
        return self * interval(other).inverse()

    def __rtruediv__(self, other):
        return interval(other) * self.inverse()

    def log(self):
        # Reuse the existing rational atanh series with its explicit tail.
        low = log_enclosure(Q(self.lo, SCALE))[0] * SCALE
        high = log_enclosure(Q(self.hi, SCALE))[1] * SCALE
        return self.raw(low.numerator//low.denominator, ceiling(high))

    def bounds(self, digits=24):
        unit = 10**(50-digits)
        low = self.lo//unit
        high = -((-self.hi)//unit)

        def render(value):
            sign = "-" if value < 0 else ""
            integer, fractional = divmod(abs(value), 10**digits)
            return f"{sign}{integer}.{fractional:0{digits}d}"

        return [render(low), render(high)]


def interval(value):
    return value if isinstance(value, Interval) else Interval(value)


def price(t):
    t = interval(t)
    return 1/(t*t.log())


def layer(p, j):
    return interval(1 + Q(1, sum(p**k for k in range(1, j+1)))).log()/interval(p).log()


def positive(value, name):
    if value.lo <= 0:
        raise ArithmeticError(name)
    return value.bounds()


def inventory(clock, exponents):
    lam = price(clock)
    margins = []
    product = 1
    for p, a in zip((2, 3, 5, 7, 11), exponents):
        product *= p**a
        if a:
            margins.append(positive(layer(p, a)-lam, f"active {p}^{a}"))
        margins.append(positive(lam-layer(p, a+1), f"inactive {p}^{a+1}"))
    # r_(p,1) < 1/(p log p) <= 1/(13 log 13) for every prime p >= 13.
    tail_margin = positive(lam-price(13), "all remaining primes")
    return {"integer": product, "price": lam.bounds(), "local_margins": margins,
            "remaining_prime_margin": tail_margin}


def certificate():
    a = interval(2520).log()
    b = interval(5040).log()
    la = a.log()
    log2 = interval(2).log()
    reward = interval(Q(31, 30)).log()
    source = inventory(a, (3, 2, 1, 1, 0))
    if source["integer"] != 2520:
        raise ArithmeticError("wrong source inventory")
    # This exact comparison uses no Euler-constant approximation.
    record_gap = (b.log()/la).log()-reward
    positive(record_gap, "G(2520) > G(5040)")
    rows = []
    for value, sign in ((Q(11), -1), (Q(113, 10), 1)):
        t = interval(value)
        stock = inventory(t, (4, 2, 1, 1, 0))
        if stock["integer"] != 5040:
            raise ArithmeticError("wrong later inventory")
        lt = t.log()
        cap = (lt/la).log()+price(t)*(a-t)
        defect = reward-price(t)*log2
        q = cap-defect
        k = (1+lt)/(t*t*lt*lt)
        derivative = k*((t-a)*defect-log2*cap)/(cap*cap)
        positive(q, "positive future pressure rise")
        positive(defect, "positive supporting-line defect")
        positive(sign*derivative, "derivative sign")
        tangent = (lt/b.log()).log()+price(t)*(b-t)
        q_minus_tangent = q-tangent
        positive(q_minus_tangent, "full FR2 local margin")
        rows.append({"clock": str(value), "inventory": stock,
                     "cap": cap.bounds(), "pressure_rise": q.bounds(),
                     "ratio": (q/cap).bounds(), "ratio_derivative": derivative.bounds(),
                     "q_minus_H": q_minus_tangent.bounds()})
    return {
        "scope": "Local rational interval certificate, not certification of external analytic inputs or Lean.",
        "source": source, "log_source": a.bounds(),
        "log_G_source_minus_log_G_5040": record_gap.bounds(),
        "rows": rows,
        "external_global_input": {
            "url": "https://arxiv.org/pdf/1112.6010v2",
            "pdf_sha256": "7bf39ed12be2f708dc6c4dca747d7d9316470c11d8f09f071cdbc79e6f1148ae",
            "locator": "Theorem 5(i) proof, printed pp.14-15, equation (24), set R and Case 2",
            "statement": "G(m) < G(2520) for every integer m > 2520, by the displayed published inputs.",
            "status": "Reused published analytic and finite-verification premises; not rerun here."
        },
        "boundaries": [
            "The counterexample concerns monotonicity with no lower source-size restriction.",
            "It does not refute a version restricted to N > 5040 or sufficiently large sources.",
            "It proves neither a Robin counterexample nor RH."
        ]
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = certificate()
    args.output.write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps({"source": result["source"]["integer"],
                      "derivatives": [row["ratio_derivative"] for row in result["rows"]]}))
