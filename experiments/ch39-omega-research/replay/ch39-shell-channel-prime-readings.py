"""Test the two shell channel populations through finite prime-power sums."""

from math import comb, isqrt
from pathlib import Path
from runpy import run_path

from flint import acb_mat, arb, arb_mat, arb_series, ctx

root = Path(__file__).parent
windows = run_path(str(root / 'ch39-dyadic-prime-window-pivot.py'))


def log_derivative(s):
    series = arb_series([s, 1], prec=2).zeta()
    return series[1] / series[0]


def prime_prefixes(degree):
    limit = 2**degree
    sieve = bytearray([1]) * (limit+1)
    for p in range(2, isqrt(limit)+1):
        if sieve[p]:
            sieve[p*p::p] = bytes((limit-p*p)//p+1)
    buckets = [[arb(0) for _ in range(6)] for _ in range(degree+1)]
    for p in range(2, limit+1):
        if not sieve[p]:
            continue
        log_p = arb(p).log()
        power = p
        while power <= limit:
            m = arb(power)
            t = m.sqrt()
            values = (arb(1), 1/m, t, 1/(m*t), m*t, 1/(m*m*t))
            row = buckets[(power-1).bit_length()]
            for i, weight in enumerate(values):
                row[i] += weight*log_p
            power *= p
    for j in range(1, degree+1):
        buckets[j] = [a+b for a, b in zip(buckets[j], buckets[j-1], strict=True)]
    return buckets


def harmonic_measures(degree):
    half = arb(1)/2
    log_pi = arb.pi().log()
    prefixes = prime_prefixes(degree)
    cosine = []
    for i, a in enumerate((half, arb(1), arb(2))):
        b, c = half-a, half+a
        if i == 0:
            total = 1 + arb.const_euler()/2 - (4*arb.pi()).log()/2
        else:
            d_b, d_c = -log_derivative(b), -log_derivative(c)
            total = (1/c+1/(c-1)+((c/2).digamma()-log_pi)/2-d_c)/(2*a)
        column = [total]
        for n in range(1, degree+1):
            x = arb(2)**n
            if i == 0:
                poles = x*(1+arb.const_euler()-x.log()) - (2*arb.pi()).log()
            else:
                poles = x/((1-b)*(c-1)) + (x**b*d_b-x**c*d_c)/(2*a)
            prime_sum = (x**c*prefixes[n][2*i+1]-x**b*prefixes[n][2*i])/(2*a)
            trivial = sum(x**(-2*j)/((2*j+b)*(2*j+c)) for j in range(1, 161))
            # Both denominator factors have positive product greater than one for j >= 161.
            tail_bound = x**(-322)/(1-x**(-2))
            value = (poles+prime_sum+trivial+arb(0).union(tail_bound))/(2*x.sqrt())
            assert abs(value) < total
            column.append(value)
        cosine.append(column)
    # Partial fractions of the three positive kernels 1/(g^2+1/4),
    # 1/((g^2+1/4)(g^2+1)), and 4/((g^2+1/4)(g^2+1)(g^2+4)).
    measures = []
    for bare, at_one, at_two in zip(*cosine, strict=True):
        first = 4*(bare-at_one)/3
        second = 4*(first-(at_one-at_two)/3)/(arb(4)-arb(1)/4)
        measures.append((bare, first, second))
    return measures


def quadrature(moments, lobatto):
    # Lobatto interior nodes integrate y(1-y) times the original measure.
    modified = [a-b for a, b in zip(moments[1:-1], moments[2:], strict=True)] if lobatto else moments
    size = len(modified)//2
    h = arb_mat([[modified[i+j] for j in range(size)] for i in range(size)])
    g = arb_mat([[modified[i+j+1] for j in range(size)] for i in range(size)])
    roots = acb_mat(h.solve(g)).eig()
    nodes = sorted((v.real for v in roots), key=lambda x: float(x.mid()))
    assert 0 < nodes[0] and nodes[-1] < 1 and all(a < b for a,b in zip(nodes,nodes[1:]))
    if lobatto:
        nodes = [arb(0), *nodes, arb(1)]
    matrix = arb_mat([[y**j for y in nodes] for j in range(len(nodes))])
    weights = matrix.solve(arb_mat([[v] for v in moments[:len(nodes)]]))
    assert all(weights[j,0] > 0 for j in range(len(nodes)))
    return [(y, weights[j,0]) for j, y in enumerate(nodes)]


def rules(measures):
    degree = len(measures)-2
    root_two = arb(2).sqrt()
    result = []
    for col in zip(*measures, strict=True):
        filtered = [3*root_two*col[n]-2*(col[abs(n-1)]+col[n+1]) for n in range(degree+1)]
        moments = [(comb(2*j,j)*filtered[0]+sum(2*(-1)**n*comb(2*j,j-n)*filtered[n]
                    for n in range(1,j+1)))/4**j for j in range(degree+1)]
        result.append((quadrature(moments, True), quadrature(moments, False)))
    return result


def window_moments(k, ruleset):
    t = arb(2).sqrt()*k
    def evaluate(rule):
        return sum(w*((1-t)**2+4*t*(1-y)).nonnegative_part().sqrt() for y,w in rule)
    return [evaluate(low).union(evaluate(high)) for low, high in ruleset]


def readings(k, ruleset, targets):
    a, b, c = window_moments(k, ruleset)
    q, d = b/a, c/b
    return (1-q-targets[0], 1-4*q+3*q*d-targets[1], windows['schur_readout'](q, d))


def show(name, values):
    print(name, [(v.lower().str(10), v.upper().str(10)) for v in values])


def calibration(ruleset, q):
    def evaluate(v, derivative=False):
        values = []
        for pair in ruleset[:2]:
            ends = []
            for rule in pair:
                if derivative:
                    ends.append(sum(w*y/(1-v*y).sqrt() for y,w in rule))
                else:
                    ends.append(sum(w*(1-v*y).nonnegative_part().sqrt() for y,w in rule))
            values.append(ends[0].union(ends[1]))
        a, b = values
        return (b-q*a)/2 if derivative else q*a-b

    for i in range(1024):
        v = (arb(i)/1024).union(arb(i+1)/1024)
        if 964 <= i <= 968:
            assert evaluate(v, derivative=True) > 0
        else:
            gap = evaluate(v)
            assert gap < 0 if i < 964 else gap > 0
    small = arb('0.43638')
    large = arb('0.43642')
    for k, sign in ((small, -1), (large, 1)):
        t = arb(2).sqrt()*k
        assert sign*evaluate(4*t/(1+t)**2) > 0
    print('Two exact calibration solutions: small root in', small.union(large),
          'and large root in', (1/(2*large)).union(1/(2*small)))


def verify():
    ctx.prec, ctx.cap = 384, 3
    run_path(str(root/'ch39-intrinsic-shell-comparisons.py'))['verify']()
    ruleset = rules(harmonic_measures(16))
    targets = windows['arithmetic_times']()
    calibration(ruleset, 1-targets[0])
    values = []
    threshold = arb('266.6444299')
    for n in range(625):
        # A pure shift has exactly the same normalized spectrum as the unshifted window.
        k = arb(n)/(624-n) if n<624 else arb(0)
        delta1, delta2, pivot = readings(k, ruleset, targets)
        assert delta1 < 0 or delta1 > 0
        assert pivot > threshold if 163 <= n <= 365 else pivot < threshold
        values.append((delta1, delta2, pivot))

    best1 = abs(values[333][0])
    best_joint_upper = abs(values[333][1])
    assert best1 < best_joint_upper
    fiber_sizes = {16*i+36*j for i in range(13) for j in range(13)}
    assert 332 in fiber_sizes and 333 not in fiber_sizes and 333 == 111*3
    for n, (d1, d2, _) in enumerate(values):
        if n != 333:
            assert abs(d1) > best1
            assert abs(d1) > best_joint_upper or abs(d2) > best_joint_upper
        if n % 6 == 0 and n != 336:
            assert abs(d1) > abs(values[336][0])
        if n in fiber_sizes and n != 332:
            assert abs(d1) > abs(values[332][0])
    for n, bounds in {
        192: [('8.762e-6', '8.782e-6'), ('3.119e-5', '3.127e-5'), ('277.636', '277.641')],
        432: [('-2.267e-4', '-2.266e-4'), ('-8.94e-4', '-8.93e-4'), ('255.170', '255.172')],
        336: [('-8.61e-6', '-8.59e-6'), ('-3.71e-5', '-3.69e-5'), ('275.57', '275.59')],
        332: [('3.66e-6', '3.68e-6'), ('1.11e-5', '1.13e-5'), ('277.02', '277.03')],
        333: [('5.79e-7', '5.95e-7'), ('-9.52e-7', '-8.94e-7'), ('276.656', '276.660')],
    }.items():
        for value, (low, high) in zip(values[n], bounds, strict=True):
            assert arb(low) < value < arb(high)
        show(n, values[n])


if __name__ == '__main__':
    verify()
