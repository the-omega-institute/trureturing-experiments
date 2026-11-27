"""检验精确二次振幅的局部维数、相邻零点读出秩与有限格点正定性。"""
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, arb_mat


def load(name):
    return run_path(str(Path(__file__).with_name(f"ch39-{name}.py")))


def real_columns(vectors):
    return arb_mat([list(row) for row in zip(*vectors, strict=True)])


def row_norm(matrix):
    return max(sum(abs(matrix[i,j]).upper() for j in range(matrix.ncols()))
               for i in range(matrix.nrows()))


def cayley_directions():
    imaginary = acb(0,1)

    def off_diagonal(i,j,z):
        direction = acb_mat(3,3)
        direction[i,j] = z
        direction[j,i] = -acb(z).conjugate()
        return direction

    central = [off_diagonal(2,j,z) for j in (0,1) for z in (1,imaginary)]
    extra = [
        acb_mat([[imaginary,0,0],[0,-imaginary,0],[0,0,0]]),
        acb_mat([[imaginary,0,0],[0,imaginary,0],[0,0,-2*imaginary]]),
        off_diagonal(1,0,1),
        off_diagonal(1,0,imaginary),
    ]
    return central,extra


def amplitude_and_tangents(up,down,k0,identity,amplitude,quadratic):
    trace,bracket = quadratic["trace"],quadratic["bracket"]
    quadratic["normal_frame"](up,down)
    commutator = bracket(up,down)
    normals = [trace(commutator*commutator).real*g - trace(commutator*g).real*commutator
               for g in (bracket(up*up,down),bracket(up,down*down))]

    def constraints(direction):
        values = [trace(normal*k0*direction) for normal in normals]
        return [z.real for z in values] + [z.imag for z in values]

    central,extra = cayley_directions()

    def generator(t):
        return sum((v*e for v,e in zip(t,central,strict=True)), acb_mat(3,3))

    def numerator_jacobian(t):
        h = generator(t)
        return real_columns(constraints(2*v*identity+2*(e+h*e+e*h))
                            for v,e in zip(t,central,strict=True))

    center = list(map(arb, ("1.031651697225703253","-.705825251234304207",
                            "-.408172603011559012",".759489497492028476")))
    original = list(map(arb, ("1.0316517","-.70582525","-.4081726",".7594895")))
    radius = arb("1e-14")
    assert all(abs(a-b)+radius < arb(".001") for a,b in zip(center,original,strict=True))
    box = [x+arb(0,radius.upper()) for x in center]
    jac = numerator_jacobian(center)
    preconditioner = arb_mat([[jac[i,j].mid() for j in range(4)] for i in range(4)]).inv()
    identity4 = arb_mat([[int(i==j) for j in range(4)] for i in range(4)])
    residual = real_columns([constraints(amplitude["cayley_numerator"](center,identity))])
    shift = row_norm(preconditioner*residual)
    slope = row_norm(identity4-preconditioner*numerator_jacobian(box))
    assert slope < 1 and shift+radius*slope < radius

    resolvent = (identity-generator(box)).inv()
    derivatives = [2*resolvent*e*resolvent for e in central+extra]
    jt = real_columns(map(constraints, derivatives[:4]))
    jy = real_columns(map(constraints, derivatives[4:]))
    assert arb("2.38e-28") < jt.det() < arb("2.39e-28")
    correction = -jt.inv()*jy
    changes = [k0*(derivatives[4+j] + sum((derivatives[i]*correction[i,j] for i in range(4)),
                                         acb_mat(3,3))) for j in range(4)]
    print("小盒中心位移、导数误差、四个 Cayley 方向的行列式",
          shift.str(8),slope.str(8),jt.det().str(11))
    return k0*(2*resolvent-identity),changes


def words(up,down,identity):
    return (identity,up,down,up*up,down*down,up*down,down*up)


def infer_coefficients(basis,directions):
    positions = ((0,0),(0,1),(0,2),(1,2),(2,0),(2,1),(2,2))
    minor = acb_mat([[word[i,j] for word in basis] for i,j in positions])
    assert abs(minor.det()) > arb("5.60e-12")
    data = acb_mat([[direction[i,j] for direction in directions] for i,j in positions])
    solution = minor.inv()*data
    return [[solution[i,j] for i in range(7)] for j in range(solution.ncols())]


def verify():
    amplitude = load("single-quadratic-flavor-amplitude")
    quadratic = load("two-zero-quadratic-source")
    gram = load("quadratic-gram-source")
    model,up_box,down_box = load("two-zero-orbit-condition")["calibration"]()
    identity = acb_mat([[int(i==j) for j in range(3)] for i in range(3)])
    amplitude["certify_exact"](model,up_box,down_box,identity)
    adj,trace = model["adj"],quadratic["trace"]

    def raw_gram(charge,index):
        site = gram["at_zero"](model,charge,index)
        return adj(site)*site

    sectors = [[raw_gram(charge,index+i) for i in range(-64,65)]
               for charge,index in ((-4,2345666),(2,964890))]
    # 冻结中央迹。邻号始终使用同一组原始 Gram 多项式系数。
    ups,downs = [[g/trace(sector[64]).real for g in sector] for sector in sectors]
    up,down = ups[64],downs[64]
    x,tangents = amplitude_and_tangents(up,down,model["K"],identity,amplitude,quadratic)
    coefficients,*changes = infer_coefficients(words(up,down,identity),[x,*tangents])

    neighbor = words(ups[65],down,identity)

    def evaluate(c):
        return sum((word*z for word,z in zip(neighbor,c,strict=True)), acb_mat(3,3))

    y = evaluate(coefficients)
    partials = [dy*adj(y) + y*adj(dy) for dy in map(evaluate,changes)]
    minor = real_columns([s[0,0].real,s[2,2].real,s[0,2].imag,s[1,2].real] for s in partials).det()
    assert arb("-3.10e-6") < minor < arb("-2.90e-6"), minor
    c0,cu,cd,cuu,cdd,cud,cdu = coefficients
    up_parts = [c0*identity+cu*u+cuu*(u*u) for u in ups]
    down_parts = [cd*d+cdd*(d*d) for d in downs]
    floor = min(abs((pu+pd+cud*(u*d)+cdu*(d*u)).det()).lower()
                for u,pu in zip(ups,up_parts,strict=True)
                for d,pd in zip(downs,down_parts,strict=True))
    assert floor > arb("7.90e-6")
    print("上型加一号的协方差雅可比子式、16641 个精确振幅行列式的统一下界",
          minor.str(8),floor.str(10))


if __name__ == "__main__":
    verify()
