"""Finite rational evaluator for the 32-bit gauge inputs and 45-bit clock."""
from fractions import Fraction as Q
from pathlib import Path

SCALE = 2**128
INVERSE_LOOP = Q(2154863261815360613622218288436129902, SCALE)
STEPS, DEGREE = 256, 13


def matrix(lines):
    return [[Q(x) for x in line.split()] for line in lines]


COEFFICIENTS = {
    "B": matrix(["199/50 27/10 44/5", "9/10 35/6 12", "11/10 9/2 -26"]),
    "D": matrix(["17/10 1/2 3/2", "3/2 3/2 1/2", "2 2 0"]),
    "C": matrix(["17/20 9/4 8", "1/4 9/4 8", "9/4 9/4 0"]),
    "BETA": list(map(Q,("41/10","-19/6","-7"))),
}
RATIOS = [Q(24,25),Q(16,15),Q(1)]
GAUGE = list(map(Q,(".467453",".63811","1.0583")))
YUKAWA = list(map(Q,".00000615 .00311 .8616 .0000135 .000268 .01401 .0000028227 .00059465 .010101".split()))


def taylor(initial, order, running, inverse_loop, *, opposite_trace=None):
    # A polynomial vector field in u_i=g_i^2 and y_a=(Yukawa_a)^2.
    coeff = [initial]
    gauge_squares = []
    rates = []
    for k in range(order):
        term = coeff[k]
        u, y = term[:3], term[3:]
        traces = [sum(y[i:i+3]) for i in (0,3,6)]
        total = 3*(traces[0]+traces[1])+traces[2]
        gauge_rates, matter_rates = [], []
        for i in range(3):
            corr = sum(b*x-d*t for b,d,x,t in zip(running["B"][i],running["D"][i],u,traces))
            gauge_rates.append(2*((running["BETA"][i] if k == 0 else 0)+corr*inverse_loop)*inverse_loop)
            gy = sum(c*x for c,x in zip(running["C"][i],u))
            if k == 0 and "NEUTRAL_YUKAWA" in running:
                gy += 3*running["NEUTRAL_YUKAWA"][i]
            for a in range(3):
                mix = y[3*(1-i)+a] if i < 2 else 0
                if i < 2 and opposite_trace is not None:
                    mix = traces[1-i] if opposite_trace else 0
                matter_rates.append((3*(y[3*i+a]-mix)+2*(total-gy))*inverse_loop)
        rates.append(gauge_rates+matter_rates)
        gauge_squares.append([sum(coeff[j][i]*coeff[k-j][i] for j in range(k+1)) for i in range(3)])
        next_g = [sum(gauge_squares[j][i]*rates[k-j][i] for j in range(k+1))/(k+1) for i in range(3)]
        next_y = [sum(coeff[j][i]*rates[k-j][i] for j in range(k+1))/(k+1) for i in range(3,12)]
        coeff.append(next_g+next_y)
    return coeff


def truncate(x):
    return Q((x.numerator*SCALE)//x.denominator, SCALE)


def gauge_readout(codes, clock, yukawa):
    assert len(codes) == 3 and all(type(k) is int and 0 <= k < 2**32 for k in codes)
    assert type(clock) is int and 0 <= clock < 2**45
    assert len(yukawa) == 9 and all(type(y) in (int,Q) for y in yukawa)
    assert all(Q(99,100)*c <= y <= Q(101,100)*c for c,y in zip(YUKAWA,yukawa))
    time = 23 + Q(10*(2*clock+1),2**46)
    step = time/STEPS
    state = []
    for k,center in zip(codes,GAUGE):
        lo,hi = [1/(center*factor)**2 for factor in (Q(1001,1000),Q(999,1000))]
        state.append(truncate(1/(lo+(hi-lo)*Q(2*k+1,2**33))))
    state.extend(truncate(y*y) for y in yukawa)
    for _ in range(STEPS):
        terms = taylor(state,DEGREE,COEFFICIENTS,INVERSE_LOOP)
        value = terms[-1]
        for term in reversed(terms[:-1]):
            value = [a+step*b for a,b in zip(term,value)]
        state = [truncate(x) for x in value]
    return time,[1/(r*u) for r,u in zip(RATIOS,state)]


def target(time, numerator):
    # The register gives floor(2^96 G_n(23)); reduce the clock exponential.
    s = -2*(time-23)/128
    poly = Q(1)
    for j in range(20,0,-1):
        poly = 1+s*poly/j
    for _ in range(7):
        poly = truncate(poly*poly)
    return truncate(Q(numerator,2**96)*poly)


def register():
    path = Path(__file__).with_name("ch39-finite-targets.tsv")
    return {int(n):int(k) for n,k in (line.split() for line in path.read_text().splitlines()[1:])}


def residual_codes(time, normalized_inverse, numerator):
    common = target(time,numerator)
    return [round((x-common)*2**42) for x in normalized_inverse]


if __name__ == "__main__":
    time,end = gauge_readout([382269474,3503200819,2112544438],13349175109351,YUKAWA)
    values = register()
    expected = {11736:[-2509323613,-2509323623,-2509323624],
                11737:[10,-1,-2], 11738:[2973221074,2973221064,2973221063]}
    for label, signed_codes in expected.items():
        read = residual_codes(time,end,values[label])
        assert read == signed_codes
        print(label, "signed residual register in units of 2^-42", read)
