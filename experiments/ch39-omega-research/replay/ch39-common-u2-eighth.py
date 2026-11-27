"""将第五壳联合根盒与第八壳热矩转换为整数包围，覆盖共用 U(2) 的四个角。"""
from pathlib import Path
from runpy import run_path
from tempfile import TemporaryDirectory
from subprocess import run

from flint import arb, acb

Q = 1 << 24
PERIOD = 1 << 16


def prepare(root, fifth, columns, levels):
    lib = lambda name: run_path(str(root/("ch39-"+name+".py")))
    data = lib("fifth-seventh-shortest-site-gates")["geometry"](
        root, fifth, columns, levels, number=2_306_888_100_093)
    base, _, _, basis, station, weights, heat, field = data
    def floor(v):
        return int((Q*v).lower().floor().unique_fmpz())
    def ceil(v):
        return int((Q*v).upper().ceil().unique_fmpz())
    def mid(v):
        return int((Q*v.mid()).floor().unique_fmpz())
    mass = lib("fifth-seventh-cutoff-separation")["moment_errors"](
        field,heat,weights,[".75553"]*19)
    slopes = (sum(a for _,a,_ in field)/(2*arb.pi()),-sum(g for g,_,_ in field)/(8*arb.pi()))
    assert levels[2] < arb(".75553")
    lines = [[Q,PERIOD], [ceil(levels[2])]]
    for t, slope in zip(mass,slopes,strict=True):
        guard = arb(".01")/abs(slope)
        lines.append([floor(t-guard),ceil(t+guard)])
    for value in (levels[0],levels[1]-levels[0],levels[2]-levels[0]):
        lines.append([floor(value),ceil(value)])
    for weight in weights:
        lines.append([floor(weight),ceil(weight)])
    near,far = (6,7,8),(11,12,18)
    for i in range(19):
        for j in range(2):
            def block(read, write):
                return sum(basis[a,i].conjugate()*station[b,j] for a,b in zip(read,write,strict=True))
            c,d = block(near,near),block(far,far)
            terms = [base[i,j]-c-d,c,d,block(far,near),-block(near,far)]
            values = [(mid(t.real),mid(t.imag)) for t in terms]
            radii = [abs(t-acb(arb(a)/Q,arb(b)/Q)) for t,(a,b) in zip(terms,values,strict=True)]
            assert all(abs(t) < 2 for t in terms)
            size_d = sum(abs(t) for t in terms[1:3])
            size_o = sum(abs(t) for t in terms[3:])
            # 每项中心三角系数的复误差小于 5/Q；最终复投影截断的误差小于 2/Q。
            error = ceil(sum(radii) + (5*(size_d+size_o)+2)/Q)
            lines.append([n for pair in values for n in pair]+[ceil(size_d),ceil(size_o),error])
    # 查表的中心三角值分别向下量化；用对称性扩展至全周期。
    lines.append([0])
    for i in range(1,PERIOD//4):
        value = (2*arb.pi()*i/PERIOD).sin()
        tick = floor(value)
        assert arb(tick)/Q <= value and value < arb(tick+1)/Q
        lines.append([tick])
    lines.append([Q])
    return "\n".join(" ".join(map(str,line)) for line in lines)+"\n"


def experiment(root, fifth, columns, levels):
    text = prepare(root, fifth, columns, levels)
    with TemporaryDirectory() as temp:
        exe = Path(temp)/"cover"
        run(["c++","-std=c++17","-O3",str(root/"ch39-common-u2-eighth.cpp"),"-o",str(exe)],check=True)
        run([str(exe)],input=text,text=True,check=True)


def verify(root):
    joint = run_path(str(root/"ch39-fifth-full-positive-joint.py"))
    experiment(root,*joint["verify"](root,return_operator="spectral"))


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)
