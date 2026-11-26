"""Minimum squared-radius excursion for origin-free antipodes on the integer shell."""
from itertools import product


def round_of_moves(states,allowed):
    sums={tuple(a+b for a,b in zip(p,q)) for p,q in product(states,repeat=2)}
    return states|{p for p in allowed-states if any(tuple(a+b for a,b in zip(p,q)) in sums for q in states)}


def main():
    cube=set(product(range(-4,5),repeat=3))
    radius=lambda p:sum(x*x for x in p)
    shell={p for p in cube if radius(p)==17}
    starts={p for p in shell if min(p)>=0}
    antipodes={tuple(-i for i in p) for p in starts}
    states=starts
    for count in (15,18,24):
        states=round_of_moves(states,shell)
        assert len(states)==count
    assert states=={p for p in shell if next(x for x in p if x%2)>0}
    assert round_of_moves(states,shell)==states and not (states&antipodes)
    assert all(round_of_moves(states|{v},shell)==shell for v in shell-states)
    allowed={p for p in cube if radius(p) in (17,19)}
    states=starts
    for count in (18,42,72):
        states=round_of_moves(states,allowed)
        assert len(states)==count
    assert states==allowed and antipodes<=states
    print("Radius 17 closure is the odd-coordinate-positive half of 48 points.")
    print("All nine antipodes are reached at excursion two by three rounds on radii squared 17 and 19.")


if __name__=="__main__":
    main()
