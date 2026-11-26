"""Finite synthesis of both clocks from one shortest-triangle return on the shell.

Each word is (-1)^sign times an integer power of C. Integer intervals verify
the recorded outputs directly against the 46 RH rows. The generator C and
its inverse are rational functions of the same ordered shell loop.
"""
from functools import lru_cache
from math import isqrt
from fractions import Fraction as Q
import json
from pathlib import Path


DEPTH = 640
DEN = 1 << DEPTH
ONE = ((DEN,DEN), (0,0))


def interval(x):
    x = Q(x)*DEN
    return (x//1, -((-x)//1))


def add(x, y):
    return x[0]+y[0], x[1]+y[1]


def neg(x):
    return -x[1], -x[0]


def times(x, y):
    corners = [a*b for a in x for b in y]
    return min(corners)//DEN, -(min(-v for v in corners)//DEN)


def multiply(x, y):
    a,b = x
    c,d = y
    term = times(b,d)
    return (add(times(a,c), (-17*term[1],-17*term[0])),
            add(times(a,d), times(b,c)))


def conjugate(x):
    return x[0], neg(x[1])


def minus(x):
    return neg(x[0]), neg(x[1])


GENERATOR = interval(Q(35888,35937)), interval(Q(455,35937))


@lru_cache(maxsize=1024)
def power(n):
    if n == 0:
        return ONE
    if n < 0:
        return conjugate(power(-n))
    half = power(n//2)
    square = multiply(half,half)
    return multiply(square,GENERATOR) if n%2 else square


def rotation(word):
    sign,n = word
    return minus(power(n)) if sign else power(n)


def difference(x, y):
    return x[0]^y[0], x[1]-y[1]


def small_step():
    u = (0,1)
    v = (1,-1)
    while rotation(v)[0][1] <= 0:
        v = difference(v,u)
    assert rotation(v)[0][0] > 0 and rotation(v)[1][0] > 0
    for _ in range(12000):
        a,b = rotation(u)[1], rotation(v)[1]
        if a[1] < b[0]:
            u,v = v,u
            a,b = b,a
        assert a[0] > b[1] > 0
        if b[1] < DEN >> 107:
            assert b[0] > DEN >> 128
            return v
        u = difference(u,v)
    raise ArithmeticError("requested small turn was not reached")


def encode(chart, chart_den, step):
    sign,tcode = chart
    # Finite Pell calibration as an aid to finding the word. Replay checks
    # the final integer power against the source chart without this choice.
    n,m = 1,0
    for _ in range(20):
        n,m = 4*n+17*m,n+4*m
    t = Q(abs(tcode),chart_den)
    x = Q(m,n)*t
    target = interval((1-17*x*x)/(1+17*x*x)), interval(2*x/(1+17*x*x))
    unit = rotation(step)
    moves = [(1,unit)]
    while True:
        k = moves[-1][0]*2
        v = rotation((0,step[1]*k))
        if v[0][1] < 0:
            break
        assert v[0][0] > 0 and v[1][0] > 0
        moves.append((k,v))
    count,current = 0,ONE
    for k,v in reversed(moves):
        candidate = multiply(current,v)
        cross = add(times(candidate[0],target[1]), neg(times(candidate[1],target[0])))
        assert cross[0] > 0 or cross[1] < 0
        if cross[0] > 0:
            count += k
            current = candidate
    return [int(sign < 0) ^ ((count*step[0])%2), step[1]*count*(1 if tcode>=0 else -1)]


def standard(chart, den):
    sign,n = chart
    return sign*Q(den*den-n*n,den*den+n*n), sign*Q(2*n*den,den*den+n*n)


def check_output(rows, table):
    den = 1<<table["bits"]
    assert len(rows)==len(table["rows"])==46
    worst = signed = Q(0)
    root = isqrt(17*DEN*DEN)
    for (idx,wa,wb),(code,_,_,a,b) in zip(rows,table["rows"]):
        assert idx==code
        x,y = rotation(wa),rotation(wb)
        ca,sa = standard(a,den)
        cb,sb = standard(b,den)
        targets = [ca,cb,ca*cb+sa*sb,ca*cb-sa*sb]
        output = [x,y,multiply(conjugate(x),y),multiply(x,y)]
        worst = max(worst, max(abs(Q(bound,DEN)-t) for t,z in zip(targets,output) for bound in z[0]))
        sine_targets = [sa,sb,ca*sb-sa*cb,ca*sb+sa*cb]
        signed = max(signed, max(abs(Q(bound,DEN)-t) for t,z in zip(sine_targets,output)
                                 for bound in times(z[1],(root,root+1))))
    assert max(worst,signed) < Q(1,2**96)
    high = max(abs(w[1]) for r in rows for w in r[1:])
    print("Largest exponent bit length",high.bit_length(),"largest 184-entry error / 2^-96",float(worst*2**96), "sine error ratio",float(signed*2**96))


def main():
    root=Path(__file__).resolve().parent
    table=json.loads((root/"ch39-five-tick-record.json").read_text())
    step=small_step()
    rows=[[idx,encode(a,1<<table["bits"],step),encode(b,1<<table["bits"],step)]
          for idx,_,_,a,b in table["rows"]]
    check_output(rows,table)
    print("Small positive step (-1 sign, signed exponent)",step)
    return rows


if __name__ == "__main__":
    main()

