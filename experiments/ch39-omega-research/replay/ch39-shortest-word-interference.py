"""Exact eight-history interference readout for the shortest matching quaternion orbit."""
from fractions import Fraction as Q
from itertools import product
from pathlib import Path
from runpy import run_path

root = Path(__file__).parent
words = run_path(str(root/"ch39-hurwitz-gauge-words.py"))
# Coordinate order is the fixed (hypercharge, weak, color) input order.
endpoint = (-13174679,6162390,8401872,13953440)


def verify():
    length = 21
    states = {}
    for bits in product((0,1),repeat=3):
        quat = (endpoint[0],*(v*(-1)**bit for bit,v in zip(bits,endpoint[1:])))
        sequence = words["factor"](quat,length)
        state = (1,0,0,0)
        for j in sequence:
            state = words["multiply"](state,words["GENERATORS"][j])
        assert state == quat
        states[bits] = state

    denominator = 8**2 * 5**length
    weights = {}
    for mask in states:
        amplitude = [0]*4
        for bits,quat in states.items():
            sign = (-1)**sum(a*b for a,b in zip(bits,mask))
            amplitude = [a+sign*b for a,b in zip(amplitude,quat)]
        weights[mask] = Q(sum(v*v for v in amplitude),denominator)
    scalar = weights[(0,0,0)]
    for j in range(3):
        mask = tuple(int(i == j) for i in range(3))
        assert weights[mask]/scalar == Q(endpoint[j+1],endpoint[0])**2
    assert all(not w for mask,w in weights.items() if sum(mask)>1)
    assert sum(weights.values()) == 1
    distance = sum(abs(w-Q(1,8)) for w in weights.values())/2
    assert Q(".54536") < distance < Q(".54537")
    print("Walsh ports (scalar, i, j, k) numerators:",*[v*v for v in endpoint],"denominator",5**length)
    print("Coherent versus dephased sign-readout total variation",distance)


if __name__ == "__main__":
    verify()
