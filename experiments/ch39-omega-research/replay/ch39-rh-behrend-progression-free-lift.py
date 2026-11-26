"""A three-term-progression-free register for finite origin-free angular comparisons."""
from itertools import permutations,product
from pathlib import Path
from runpy import run_path

grow=run_path(str(Path(__file__).with_name("ch39-rh-behrend-radial-budget.py")))["round_of_moves"]
radius=lambda p:sum(v*v for v in p)
universe=set(product(range(-5,6),repeat=3))
shell={p for p in universe if radius(p)==17}
seeds={p for p in shell if min(p)>=0}


def free(room):
    return all(tuple(2*b-a for a,b in zip(p,q)) not in room for p in room for q in room if p!=q)


def closure(room):
    current=seeds
    counts=[len(current)]
    for _ in range(len(room)):
        after=grow(current,room)
        if current==after:
            return current,counts
        current=after
        counts.append(len(current))
    raise AssertionError("finite increasing chain exceeded its cardinality")


def main():
    off={p for p in universe if radius(p)==19}
    assert not free(shell|off)
    successful=[]
    for v in off:
        room=shell|{v}
        assert free(room)
        if closure(room)[0]==room:
            successful.append(v)
    assert len(successful)==18
    coordinated={tuple(sorted(p)) for p in universe if sum(p)%2
                 and abs(radius(p)-17)<10 and radius(p)!=17}
    permitted=[set(permutations(shape)) for shape in coordinated
               if free(shell|set(permutations(shape)))]
    assert len(permitted)==4 and not shell<=closure(shell.union(*permitted))[0]
    permuted=set(permutations((-3,-3,3)))
    assert len(permuted)==3 and free(shell|permuted) and closure(shell|permuted)[0]==shell|permuted
    # Full signed permutation symmetry forces whole absolute-coordinate orbits.
    shapes={tuple(sorted(map(abs,p))) for p in universe if sum(p)%2 and abs(radius(p)-17)<10 and radius(p)!=17}
    candidates=[]
    for shape in shapes:
        orbit={p for p in universe if tuple(sorted(map(abs,p)))==shape}
        if free(shell|orbit):
            candidates.append(orbit)
    room=shell.union(*candidates)
    assert len(candidates)==1 and {tuple(sorted(map(abs,p))) for p in room-shell}=={(0,0,5)}
    assert not shell<=closure(room)[0]
    eight=set(product((-3,3),repeat=3))
    room=shell|eight
    reached,rounds=closure(room)
    assert reached==room and free(room) and len(room)==56
    # Integer encoding after adding 4. Base 17 exceeds twice the largest digit.
    encode=lambda p:sum((x+4)*17**j for j,x in enumerate(p))
    numbers={encode(p) for p in room}
    assert len(numbers)==56 and all(2*b-a not in numbers for a in numbers for b in numbers if a!=b)
    print("18 single radius-19 choices succeed. Coordinate permutations require excursion 10, with three bridge points sufficient.")
    print("Signed permutation symmetry and progression freedom: optimal squared radius excursion 10;",
          "eight radius-27 bridge points, counts",rounds)


if __name__=="__main__":
    main()
