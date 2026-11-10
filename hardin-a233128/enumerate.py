"""Exact complementary-pair row transfer and the registered recurrence certificate."""
import itertools
from fractions import Fraction
import json
from pathlib import Path

rows = [p for p in itertools.product(range(3), repeat=7)
        if all(x != y for x, y in zip(p, p[1:]))]
def bits(p):
    return tuple(int(p[j + 2] != p[j]) for j in range(5))
classes = sorted({min(bits(p), bits(p)[::-1]) for p in rows})
def class_of(p):
    return classes.index(min(bits(p), bits(p)[::-1])) + 1
representatives = [min(p for p in rows if class_of(p) == i + 1) for i in range(20)]
def scan(s, p):
    a, b = s
    weight = 1
    for x in p:
        weight *= (2 if x == 0 else (2 if a else 1) if x == 1
                   else (2 if b else 1) if a else 0)
        a |= x == 1
        b |= x == 2
    return (a, b), weight
transfer = []
for i, p in enumerate([representatives[0]] + representatives):
    row = [0] * 21
    for q in rows:
        if all(x != y for x, y in zip(p, q)):
            state, weight = scan((True, i != 0), q)
            row[class_of(q) if state[1] else 0] += weight
    transfer.append(row)
initial = [0] * 21
for q in rows:
    if q[0] == 0:
        state, weight = scan((False, False), q[1:])
        initial[class_of(q) if state[1] else 0] += weight
coefficients = [7168, -19070976, 26086473728, -21081846972416,
    10809917287956480, -3649692508960587776, 824336623994081902592,
    -124561639357723747287040, 12414326720213201465638912,
    -792885332473813431379558400, 31062140009178282014900486144,
    -704387882353381376417601683456, 8398185226512019784915658735616,
    -40564819207303340847894502572032]
distribution = [initial]
for _ in range(14):
    distribution.append([sum(distribution[-1][i] * transfer[i][k] for i in range(21))
                         for k in range(21)])
assert all(distribution[14][k] == sum(coefficients[j] * distribution[13-j][k]
           for j in range(14)) for k in range(21))
assert len(rows) == 192 and len(classes) == 20 and transfer[0][0] == 128
assert [sum(v) for v in distribution[:7]] == [528, 1495040, 4653056000,
    14800557965312, 47441169369530368, 152534173898504142848,
    491056502853059022422016]
matrix = [[Fraction(x) for x in row] for row in distribution[:14]]
rank = 0
for column in range(21):
    pivot = next((r for r in range(rank, len(matrix)) if matrix[r][column]), None)
    if pivot is None:
        continue
    matrix[rank], matrix[pivot] = matrix[pivot], matrix[rank]
    scale = matrix[rank][column]
    matrix[rank] = [x / scale for x in matrix[rank]]
    for r in range(rank + 1, len(matrix)):
        scale = matrix[r][column]
        matrix[r] = [x - scale * y for x, y in zip(matrix[r], matrix[rank])]
    rank += 1
assert rank == 14
assert 128**14 == sum(coefficients[j] * 128**(13-j) for j in range(14))
result = dict(krylov_rank=rank, recurrence_root_128=True,
              bits=classes, reps=representatives,
              sizes=[sum(class_of(p) == i + 1 for p in rows) for i in range(20)],
              first=initial, T=transfer, D=distribution, C=coefficients)
Path(__file__).with_name('data.json').write_text(json.dumps(result, indent=2) + '\n')
print('192 rows; 20 completed classes; 21 states; exact order-14 vector identity verified')
