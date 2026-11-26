"""有限公平位抽样和整数接受窗怎样选出十一点熵极小解的 CP 读数。"""
from pathlib import Path
from runpy import run_path

from flint import arb, arb_mat, ctx, fmpq, fmpq_mat

entropy = run_path(str(Path(__file__).with_name("ch39-neutral-entropy-weights.py")))
count = 2**120
tickets = 2**160
unit = 2**112
window = 2**12


def integer_part(value):
    integer = value.floor().unique_fmpz()
    assert integer is not None
    return int(integer)


def finite_tables(rows, times, reference):
    xs = [r[1] for r in rows]
    prior = [v/sum(xs) for v in xs]
    masses = [integer_part(tickets*v) for v in prior[1:]]
    masses.insert(0, tickets-sum(masses))
    assert min(masses) > 0
    tilt_error = max(abs((arb(n)/(tickets*v)).log()) for n, v in zip(masses, prior, strict=True))
    assert tilt_error < arb(2)**-90
    scores = [[integer_part(unit*v) for v in row[0][1:]] for row in rows]
    targets = [integer_part(unit*v) for v in [*times[1:6], arb(0), arb(0)]]
    flavor = [(integer_part(unit*r[2]), integer_part(unit*r[2]*(1-r[4]/reference))) for r in rows]
    return masses, scores, targets, flavor


def accepts(counts, scores, targets):
    assert len(counts) == len(scores) and min(counts) >= 0 and sum(counts) == count
    totals = [sum(k*v for k, v in zip(counts, values, strict=True))
              for values in zip(*scores, strict=True)]
    return all(abs(v-count*t) <= count*window for v, t in zip(totals, targets, strict=True))


def accepted_example(scores, targets):
    # 在整数读数表上作有理消元。末三项取此前求得的熵极小点邻域中的有理数。
    tail = [fmpq(v, 10**17) for v in (6458424719409188, 16361116312614066, 4697193110629481)]
    vectors = [[unit, *score] for score in scores]
    matrix = fmpq_mat(list(zip(*vectors[:8], strict=True)))
    residual = [t-sum(w*v[j] for w,v in zip(tail, vectors[8:], strict=True))
                for j,t in enumerate([unit, *targets])]
    first = matrix.solve(fmpq_mat([[v] for v in residual])).entries()
    weights = [*first, *tail]
    assert min(weights) > 0 and sum(weights) == 1
    counts = [int((count*v).floor()) for v in weights[1:]]
    counts.insert(0, count-sum(counts))
    assert accepts(counts, scores, targets)
    return counts


def verify():
    rows, times, _, _, _, reference = entropy["joint"]["columns_at"](2_903_600_000, entropy["ranges"]["indices"])
    ctx.prec = 640
    p = entropy["projection_weights"](rows, times)
    masses, scores, targets, flavor = finite_tables(rows, times, reference)
    assert sum(masses) == tickets and all(pair[0] > 0 for pair in flavor)
    assert not accepts([count]+[0]*10, scores, targets)
    sample = accepted_example(scores, targets)
    total_eta, total_cp = [sum(k*v for k, v in zip(sample, values, strict=True))
                           for values in zip(*flavor, strict=True)]
    assert fmpq(6838118,10**8) < fmpq(total_cp,total_eta) < fmpq(6838120,10**8)
    # 用八个独立约束列重建熵对偶梯度的常数项和七个斜率。
    basis = arb_mat([r[0] for r in rows[:8]])
    dual = basis.solve(entropy["column"]([(v/r[1]).log() for v, r in zip(p[:8], rows[:8], strict=True)]))
    norm = sum(abs(v) for v in dual.entries()[1:])
    assert norm < 2**48 and min(p) > arb("4e-11")
    assert all(abs(v) <= 1 for r in rows for v in r[0][1:])
    rounding = arb(22)/count
    quantization = arb(2)/unit
    tolerance = arb(window)/unit
    assert rounding+quantization < tolerance
    log_prior_error = arb(2)**-90
    band_error = norm*(tolerance+quantization)+log_prior_error
    witness_error = rounding**2/min(p)+norm*rounding+log_prior_error
    types_log = 11*arb(count+1).log()
    norm_cutoff = arb(2)**-24
    exponent = count*(norm_cutoff**2/2-band_error-witness_error)-2*types_log
    assert exponent/arb(2).log() > 2**70
    xs = [r[1] for r in rows]
    divergence = sum(v*(v*sum(xs)/x).log() for v,x in zip(p,xs,strict=True))
    bits_low = (divergence-band_error-types_log/count)/arb(2).log()
    bits_high = (divergence+witness_error+types_log/count)/arb(2).log()
    assert arb("50.08071989") < bits_low and bits_high < arb("50.08071991")
    loss = entropy["ranges"]["cp_loss"](rows, p, reference)
    sensitivity = max(abs(r[2]*(1-r[4]/reference-loss)) for r in rows)
    assert all(r[2] > 1 for r in rows)
    # 味分母至少一，向下量化分子分母只再引入有限表误差。
    assert abs(loss) < 1 and all(abs(1-r[4]/reference)<1 for r in rows)
    error = sensitivity*norm_cutoff+arb(2)/(unit-1)
    assert error < arb(".000000257")
    print("一组通过有限窗的整数次数", *sample)
    print("接受事件的每轮负对数概率（位）", bits_low.str(15), bits_high.str(15))
    print("条件 CP 相对亏额的中心与误差上界", loss.str(15), error.str(12))
    print("条件失败概率上界为 2 的负 2^70 次幂")


if __name__ == "__main__":
    verify()
