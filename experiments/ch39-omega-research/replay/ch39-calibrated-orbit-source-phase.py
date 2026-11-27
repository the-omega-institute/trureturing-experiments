"""固定两个实际零点的十项见证读数后，检查整个源协方差轨道的闭路相位。"""
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb


def verify():
    model = run_path(str(Path(__file__).with_name("ch39-two-zero-quark-window.py")))
    adj = model["adj"]
    up,down = [model["prepare"](q,arb(x),arb(y)) for q,x,y in (
        (-4,"-0.680061861483031","0.873223965353125"),
        (2,"0.565247837594779","-0.630343671648397"),
    )]
    relative = up*down.inv()
    _,frame = model["eig"](adj(relative)*relative)
    k = model["K"]
    forward, backward = [], []
    for i in (0,1):
        line = acb_mat([[frame[j,i]] for j in range(3)])
        forward.append((down.inv()*line)*(adj(line)*down*k))
        backward.append((k.inv()*down.inv()*line)*(adj(line)*down))
    forward.append(k-forward[0]-forward[1])
    backward.append(k.inv()-backward[0]-backward[1])
    identity = acb_mat([[int(i==j) for j in range(3)] for i in range(3)])
    shell = model["shell"]
    ports = [shell["soft"]["pole"],*shell["boundary_sites"]()]
    gram = model["load"]("native-rh-ray-pvm")["gram_cholesky"](Path(__file__).resolve().parent)
    raw = acb_mat([[gram[j][ports.index(i)] if i in ports else 0 for j in (1,2,3)]
                   for i in range(19)])
    for basis, floor in ((identity,arb("3e6")), (adj(model["sources"])*raw,arb("6e10"))):
        f = [basis.inv()*x for x in forward]
        g = [x*basis for x in backward]
        gram_entry, inverse_entry = {}, {}
        for i in range(3):
            for j in range(3):
                key = (int(i==0)-int(j==0),int(i==1)-int(j==1))
                gram_entry[key] = gram_entry.get(key,0)+(f[i]*adj(f[j]))[2,0]
                inverse_entry[key] = inverse_entry.get(key,0)+(adj(g[i])*g[j])[0,2]
        coeff = {}
        for (i,j),a in gram_entry.items():
            for (u,v),b in inverse_entry.items():
                key = (i+u,j+v)
                coeff[key] = coeff.get(key,0)+a*b
        # Im sum c_k e^{ikθ} 用共轭配对的实三角多项式表示。
        harmonics = {key:(c-coeff[(-key[0],-key[1])].conjugate())/acb(0,2)
                     for key,c in coeff.items()}
        constant = harmonics.pop((0,0)).real
        bound = constant-sum(abs(z) for z in harmonics.values())
        assert bound > floor
        print("整个两相位轨道的有向闭路读数下界",bound.lower().str(14))


if __name__ == "__main__":
    verify()
