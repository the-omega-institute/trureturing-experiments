"""冻结四个精确二次源的系数，检查相邻零点和共同对角输入处的正定性。"""
from pathlib import Path
from runpy import run_path

from flint import acb_mat, arb, arb_mat


def load(name):
    return run_path(str(Path(__file__).with_name(f"ch39-{name}.py")))


def coordinates(matrix):
    diagonal = [matrix[i,i].real for i in range(3)]
    off_diagonal = [matrix[i,j] for i,j in ((0,1),(0,2),(1,2))]
    return diagonal + [part for z in off_diagonal for part in (z.real,z.imag)]


def refine(center, outer_radius, equations, trig):
    x0,y0,_ = center
    x,y = x0,y0
    for _ in range(3):
        value,jacobian = trig(x,y)
        step = (equations*jacobian).solve(equations*value)
        x,y = (x-step[0,0]).mid(),(y-step[1,0]).mid()
    radius = arb("1e-12")
    pad = arb(0,radius.upper())
    assert abs(x-x0)+radius < outer_radius and abs(y-y0)+radius < outer_radius
    value,center_jacobian = trig(x,y)
    _,box_jacobian = trig(x+pad,y+pad)
    preconditioner = (equations*center_jacobian).inv()
    correction = preconditioner*equations*value
    error = arb_mat([[1,0],[0,1]]) - preconditioner*equations*box_jacobian
    assert all(abs(correction[i,0])/radius + sum(abs(error[i,j]) for j in range(2)) < 1
               for i in range(2))
    return x+pad,y+pad


def words(up,down,quadratic):
    identity = acb_mat([[int(i==j) for j in range(3)] for i in range(3)])
    return (identity,up,up*up,*quadratic["mixed_words"](up,down))


def coefficients_at_roots(model,up_site,down_site):
    quadratic = load("quadratic-gram-source")
    orbit = load("two-zero-orbit-condition")
    four = load("two-zero-quadratic-source")
    adj,trace = model["adj"],four["trace"]
    up,down = adj(up_site)*up_site,adj(down_site)*down_site
    # 依照同一顺序，把七个厄米基元限制到一个可逆的实七阶子式。
    rows = (0,2,3,5,6,7,8)
    columns = [coordinates(word) for word in words(up,down,quadratic)]
    minor = arb_mat([[col[i] for col in columns] for i in rows])
    assert minor.det() > arb("3.9e-20")
    frame,commutator_length = four["normal_frame"](up/trace(up).real,down/trace(down).real)
    parts,_ = orbit["orbit"](model,up_site,down_site)
    harmonics = four["orbit_coefficients"](parts,adj,frame)
    centers,radius = four["root_boxes"](harmonics,commutator_length)
    four["cover"](harmonics,centers,radius)
    equations = arb_mat([[harmonics[i,j] for j in range(7)] for i in (1,2)])
    for center in centers:
        x,y = refine(center,radius,equations,four["trigonometric_columns"])
        change = orbit["phase_sum"](parts,orbit["phases_at"](x,y))
        values = coordinates(change*adj(change))
        coeff = minor.solve(arb_mat([[values[i]] for i in rows])).entries()
        yield center,coeff


def verify():
    quadratic = load("quadratic-gram-source")
    model,_,_ = load("two-zero-orbit-condition")["calibration"]()
    adj = model["adj"]
    trace = load("two-zero-quadratic-source")["trace"]
    sites = [[quadratic["at_zero"](model,charge,n+i) for i in (-1,0,1)]
             for charge,n in ((-4,2345666),(2,964890))]
    grams = [[adj(a)*a for a in sector] for sector in sites]
    up,down = grams[0][1],grams[1][1]
    su,sd = arb(1)/1000,arb(1)/20
    assert trace(up).real > su and trace(down).real > sd
    roots = coefficients_at_roots(model,sites[0][1],sites[1][1])
    for root,(center,coeff) in enumerate(roots):
        scalar_words = (1,su,su**2,sd,2*su*sd,0,sd**2)
        scalar = sum(c*z for c,z in zip(coeff,scalar_words,strict=True))
        lower = (-393,-298,-212,-298)[root]
        assert lower < scalar < lower+1
        print("相位中心和同一对角输入的首个本征值",center[0].str(7),scalar.str(10))

        for normalized in (False,True):
            positive = []
            for i,u in enumerate(grams[0],start=-1):
                for j,d in enumerate(grams[1],start=-1):
                    if i==j==0:
                        continue
                    gu,gd = u,d
                    if normalized:
                        gu,gd = u*(trace(up).real/trace(u).real),d*(trace(down).real/trace(d).real)
                    source = sum((c*w for c,w in zip(coeff,words(gu,gd,quadratic),strict=True)),acb_mat(3,3))
                    levels = source.eig(algorithm="rump")
                    assert all(z.imag.contains(0) for z in levels)
                    if root==2 and j==0:
                        assert all(z.real > 3 for z in levels)
                        positive.append((i,j))
                    else:
                        assert any(z.real < arb("-.47") for z in levels)
            print("单位迹延拓、正定的上/下编号偏移",normalized,positive)


if __name__ == "__main__":
    verify()
