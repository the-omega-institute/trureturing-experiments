"""预登记的四壳两字符共同十二式见证在第十一壳上的三项。"""
from pathlib import Path
from runpy import run_path

from flint import arb


def experiment(root, fifth, columns, levels):
    readout = run_path(str(root/"ch39-cyclic-nine-tenth.py"))["readout"]
    values = readout(root,fifth,columns,levels,"four-cyclic-twelve",2_306_888_100_096)
    for value,(lo,hi) in zip(values,((".88974",".89000"),("-.0020140","-.0020074"),
                                    ("-.07425657725","-.07425657714")),strict=True):
        assert arb(lo)<value<arb(hi)
    return values


if __name__=="__main__":
    root = Path(__file__).resolve().parent
    joint = run_path(str(root/"ch39-fifth-full-positive-joint.py"))
    experiment(root,*joint["verify"](root,return_operator="spectral"))
