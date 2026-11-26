"""只激活十五维的首个弱三重态时，求正的三支一圈交汇与最少份数。"""
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, ctx

matter = run_path(str(Path(__file__).with_name("ch39-charge-matter-determinant-order.py")))
soft = matter["soft"]


def verify():
    ctx.prec = 224
    height = acb.zeta_zero(11431).imag
    star, neutral, charged = soft["vector_pencil"](height)
    geo = soft["geo"]
    radius = geo["theta"](804)/geo["theta"](12)
    assert radius**2 > 36
    lift = acb_mat([[arb(1)/19]*19 for _ in range(19)])
    neutral_exposure = 18*radius.log()-(star+neutral+lift).det().real.log()/2
    charged_exposure = 19*radius.log()-(star+charged).det().real.log()/2

    mz, em, weak = arb("91.1876"), arb("127.930"), arb(".23129")
    bare = (3*(1-weak)*em/5, weak*em, 1/arb(".1180"))
    beta = (arb(41)/10, -arb(19)/6, -7)
    base = [a+7*(j*neutral_exposure+k*charged_exposure)/(4*arb.pi())
            for a,j,k in zip(bare, (0,2,3), (5,3,2))]
    # The remaining active triplet has Dynkin indices (9/5,2,0).
    def projection(row):
        return 10*row[0]-9*row[1]-row[2]
    log_scale = 2*arb.pi()*projection(base)/projection(beta)
    values = [a-b*log_scale/(2*arb.pi()) for a,b in zip(base, beta)]
    mu, common = mz*log_scale.exp(), values[2]
    total = 3*arb.pi()*(values[1]-common)/4
    assert arb("1.03725489070050e14") < mu < arb("1.03725489070052e14")
    assert arb("82.05344460100") < common < arb("82.05344460101")

    matrices = {q:star+matter["nonspokes"](q, height) for q in (1, 4, 6)}
    floors = {q:soft["lowest"](m, arb(0), arb(4)) for q,m in matrices.items()}
    assert floors[6] < floors[4] < floors[1]
    # Enclose the next triplet eigenvalue to locate the first unwanted component.
    eigen = sorted((z.real for z in matrices[6].eig(multiple=True)), key=float)
    assert len(eigen) == 19 and eigen[0].overlaps(floors[6]) and eigen[1] > floors[4]
    maximum = -(1-(floors[4]-floors[6])/radius**2).log()/2
    count = int((total/maximum).ceil().unique_fmpz())
    assert count == 1770 and (count-1)*maximum < total < count*maximum
    mass2 = radius**2*(-2*total/count).exp()-floors[6]
    light_ratio = (mass2+floors[6]).sqrt()/radius
    next_ratio = (mass2+floors[4]).sqrt()/radius
    assert 0 < light_ratio < 1 < next_ratio
    exposure = -count*light_ratio.log()
    fitted = [a-2*t*exposure/(3*arb.pi()) for a,t in zip(values,(arb(9)/5,2,0))]
    assert all(abs(v-common) < arb("1e-13") for v in fitted)
    # Conditional slope above the match if the same broken-group running continues.
    weak_beta = beta[1]-arb(7)/2*(19*3+18*2)+arb(4)/3*count*2
    pole_ratio = (2*arb.pi()*common/weak_beta).exp()
    assert pole_ratio < arb("1.13")
    print("Matching scale GeV, common inverse, required copies", mu.str(16), common.str(15), count)
    print("Bare mass squared, light and next mass ratios", mass2.str(15), light_ratio.str(13),
          next_ratio.str(13), "one-loop pole scale ratio <=", pole_ratio.str(12))


if __name__ == "__main__":
    verify()
