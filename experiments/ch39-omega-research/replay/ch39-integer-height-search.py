"""Find the best 64-bit height over the full rational outer band, with no zero list."""
from functools import cache
from heapq import heappop, heappush
from math import factorial
from pathlib import Path
from runpy import run_path

lib = run_path(str(Path(__file__).with_name("ch39-integer-hardy-readout.py")))
Ball, UNIT, GRID = (lib[k] for k in ("Ball","UNIT","GRID"))
ONE, ZERO, integer = (lib[k] for k in ("ONE","ZERO","integer"))
times, log_integer = lib["times"], lib["log_integer"]
ORDER = 56




@cache
def taylor(middle):
    g = integer(middle)+ONE.scale(1,2)
    terms = lib["dirichlet_terms"](g)
    real = [ZERO for _ in range(ORDER+1)]
    imag = [ZERO for _ in range(ORDER+1)]
    for n,(re,im) in enumerate(terms[1:4096],1):
        log_n = log_integer(n)
        for j in range(ORDER+1):
            real[j] += re
            imag[j] += im
            if j < ORDER:
                re,im = (im*log_n).scale(1,j+1), -(re*log_n).scale(1,j+1)
    # Uniform integral remainder on |delta|<=1/2, sum of absolute amplitudes < 128.
    log_half = log_integer(4096).scale(1,2)
    bound = integer(128)
    for j in range(1,ORDER+2):
        bound = (bound*log_half).scale(1,j)
    assert bound.hi*(10**38)<UNIT
    return g,real,imag,terms[2],bound.hi


def polynomial(coefficients,delta,err):
    value = coefficients[-1]
    for a in reversed(coefficients[:-1]):
        value = value*delta+a
    return value+Ball(-err,err)


def short_phase(angle):
    # Horner form for exp(i angle), |angle| <= log(2)/2.
    assert max(abs(angle.lo),abs(angle.hi)) < UNIT//2
    real,imag = ONE,ZERO
    for j in range(48,0,-1):
        factor = angle.scale(1,j)
        real,imag = ONE-factor*imag, factor*real
    error = -(-UNIT//(2**49*factorial(49)))
    return real+Ball(-error,error),imag+Ball(-error,error)


def modulus_squared(u):
    return ((integer(3)-u.scale(2)).square()*(integer(226)+u.scale(144))).scale(1,169)


CRITICAL = ONE.scale(-59,108)
ENVELOPE = modulus_squared(CRITICAL).sqrt()
ROOT_TWO = integer(2).sqrt().hi


def score_bounds(left,right,base):
    g = Ball(left*(UNIT//GRID),right*(UNIT//GRID))
    center,real,imag,two,err = taylor(base)
    delta = g-center
    phase2 = times(two,short_phase(-(delta*lib["LOG_TWO"])))
    u = phase2[0].scale(2)
    u = Ball(max(-ROOT_TWO,u.lo),min(ROOT_TWO,u.hi))
    # Clamp the cubic's critical point to the current real-part interval.
    def clip(x):
        return min(u.hi,max(u.lo,x))
    top = modulus_squared(Ball(clip(CRITICAL.lo),clip(CRITICAL.hi)))
    bottom = min(modulus_squared(Ball(x,x)).lo for x in (u.lo,u.hi))
    w = Ball(max(0,bottom),top.hi).sqrt()*(g.square()+ONE.scale(1,4)).reciprocal()

    phase_n = (ONE,ZERO)
    for _ in range(12):
        phase_n = times(phase_n,phase2)
    tail = times(phase_n,lib["tail_coefficient"](g))
    a = polynomial(real,delta,err)+tail[0]
    b = polynomial(imag,delta,err)+tail[1]
    amp = (a.square()+b.square()).sqrt()
    return w*(ONE+amp.scale(1,256)).reciprocal()


def verify():
    heap = []
    best, leader = -1, None
    evaluated = 0
    def insert(lo,hi,base):
        assert lo<hi
        box = score_bounds(lo,hi,base)
        if box.hi<best:
            return
        heappush(heap,(-box.hi,lo,hi,base))
    for k in range(10603,10622):
        insert(k*GRID,(k+1)*GRID-1,k)
    # Track singletons separately so close score intervals cannot conceal a tie.
    candidates = {}
    while heap:
        negative,lo,hi,base = heappop(heap)
        if -negative < best:
            break
        middle = (lo+hi)//2
        score = score_bounds(middle,middle,base)
        evaluated += 1
        if score.hi >= best:
            candidates[middle] = score
        if score.lo > best:
            best,leader = score.lo,middle
        for a,b in ((lo,middle-1),(middle+1,hi)):
            if a < b:
                insert(a,b,base)
            elif a==b:
                value = score_bounds(a,b,base)
                if value.hi >= best:
                    candidates[a] = value
                if value.lo > best:
                    best,leader = value.lo,a
    assert leader is not None and all(x==leader or v.hi<best for x,v in candidates.items())
    tail = ENVELOPE*(integer(10622).square()+ONE.scale(1,4)).reciprocal()
    assert tail.hi < best
    print("Global winning numerator",leader,"midpoints",evaluated,"retained scores",len(candidates))
    return leader


if __name__=="__main__":
    verify()
