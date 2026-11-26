"""Replay the six crossing signs and both critical signs with integer arithmetic."""
from fractions import Fraction
import json
from pathlib import Path


class Bounds:
    def __init__(self, left, right=None):
        self.lo = Fraction(left)
        self.hi = self.lo if right is None else Fraction(right)
        assert self.lo <= self.hi

    def __add__(self, other):
        return Bounds(self.lo+other.lo, self.hi+other.hi)

    def __sub__(self, other):
        return Bounds(self.lo-other.hi, self.hi-other.lo)

    def __mul__(self, other):
        values = [x*y for x in (self.lo, self.hi) for y in (other.lo, other.hi)]
        return Bounds(min(values), max(values))

    def __truediv__(self, other):
        assert other.lo*other.hi > 0
        return self*Bounds(1/other.hi, 1/other.lo)

    def __pow__(self, n):
        assert self.lo > 0 and isinstance(n, int) and n > 0
        return Bounds(self.lo**n, self.hi**n)

    def positive_scale(self):
        """Round a positive interval outward to bound sizes of products and sums."""
        assert self.lo > 0
        step = Fraction(2)**(self.lo.numerator.bit_length()-self.lo.denominator.bit_length()-150)
        return Bounds((self.lo//step)*step, -((-self.hi)//step)*step)

    def clipped_below_one(self):
        return Bounds(max(1, self.lo), max(1, self.hi))


def read_cell(triple):
    lo, hi, shift = triple
    return Bounds(lo, hi)*Bounds(Fraction(2)**shift)


def enclose_root(cell, power, target):
    assert cell.lo > 0 and cell.lo**power < target.lo <= target.hi < cell.hi**power
    return cell


def sign_for_sample(row, inputs):
    unit = row.get("unit", 2_100_000_000_000)
    poles = [c*Bounds(unit**2) for c in inputs["reference"]]
    rho = inputs["rho"]
    pivot, = inputs["pivot"]
    q02, = inputs["q02"]
    target, = inputs["target"]
    assert len(rho) == len(poles) == 21 and len(row["roots"]) == len(inputs["lambda"]) == 6
    assert all(x.lo > 0 for x in rho) and all(x.hi < y.lo for x, y in zip(poles, poles[1:]))
    if "unit" in row:
        q2 = poles[0]
        s = enclose_root(read_cell(row["ratio2"]), 4, (q02/q2)**11*Bounds(8**3*24))
        assert s.lo > 24
        pure_count = 0
    else:
        s = Bounds(row["ratio2"])
        t = (s/Bounds(8)).clipped_below_one()
        v = (s/Bounds(24)).clipped_below_one()
        q2 = enclose_root(read_cell(row["q2"]), 11, q02**11/(t**3*v))
        assert all(c.hi < q2.lo or c.lo > q2.hi for c in poles)
        pure_count = sum(c.hi < q2.lo for c in poles)
    mf2 = q2/s
    assert (mf2*Bounds(8)).lo > 1000**2 and poles[0].lo > 1000**2
    product = Bounds(1)
    active = 25*pure_count
    for c in poles[:pure_count]:
        product = (product*(q2/c)**25).positive_scale()
    for lam, n, block in zip(inputs["lambda"], (7, 6, 3, 3, 2, 2), row["roots"]):
        assert len(block) == pure_count+1
        load = lam*mf2
        def secular(x):
            point = Bounds(x)
            resolvent = sum((r/(point-c) for r, c in zip(rho, poles)), Bounds(0))
            return point-Bounds(125**2)-load*(pivot+point*resolvent)
        for index, entry in enumerate(block):
            root = read_cell(entry)
            assert root.lo > (poles[index-1].hi if index else 1000**2)
            assert root.hi < poles[index].lo
            assert secular(root.lo).hi < 0 and secular(root.hi).lo > 0
            if root.hi < q2.lo:
                active += n
                product = (product*(q2/root)**n).positive_scale()
            else:
                assert root.lo > q2.hi
    t = (s/Bounds(8)).clipped_below_one()
    v = (s/Bounds(24)).clipped_below_one()
    right = target*t**201*v**23
    left = product**11
    if left.hi < right.lo:
        return -1, active
    assert left.lo > right.hi
    return 1, active


def main():
    record = json.loads(Path(__file__).with_name("ch39-rh-clock-integer-data.json").read_text())
    inputs = {k: list(map(read_cell, v)) for k, v in record["inputs"].items()}
    rows = record["samples"]
    assert [row["ratio2"] for row in rows[:6]] == [
        "1.87758", "1.87768", "26.8343", "26.8364", "13116750", "13119360"]
    assert all("unit" not in row for row in rows[:6])
    assert [row["unit"] for row in rows[6:]] == [275_782_218_980, 275_782_218_981]
    result = [sign_for_sample(row, inputs) for row in rows]
    assert result == [(-1,119), (1,119), (1,71), (-1,71), (-1,23), (1,23), (1,23), (-1,23)]
    print("Eight exact rational comparisons verified:", result)


if __name__ == "__main__":
    main()
