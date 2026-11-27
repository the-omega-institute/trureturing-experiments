"""坐标置换等变的六参数门若匹配两壳四式，短边必须混合两个交换模。"""
from pathlib import Path
from runpy import run_path


def verify(root):
    def lib(name):
        return run_path(str(root/("ch39-"+name+".py")))

    operator = lib("fifth-full-positive-joint")["verify"](root,return_operator="spectral")
    torus = lib("two-bond-phase-torus")
    for bound, reverse in ((".035",False),(".05",True)):
        data = torus["prepare"](root,*operator,leakage=bound,antidiagonal=reverse)
        print("长度 12 门的", "对角元" if reverse else "非对角元", "模至多", bound, flush=True)
        torus["cover"](root,data)


if __name__=="__main__":
    verify(Path(__file__).resolve().parent)
