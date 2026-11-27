"""第七至十一壳十四个端点的三循环等变共同共轭。"""
from runpy import run_path
from pathlib import Path

PARAMETERS = (
    "0.375", "1.5", "-0.375", "-0.125",
    "0.375", "0.125", "-1.21312523282322", "-1.25",
    "2.75", "0.5", "0.75", "-0.375",
    "-1.375", "-2.125", "0.375", "-1.125",
    "-1", "-0.125", "-1.875", "0.5",
    "0.559972819428683", "-0.125", "-4.5", "-0.625",
    "3.125", "-1.625", "-0.25", "-0",
    "-0.0297188491912269", "-1.5", "6.5", "-1.375",
    "-0.412667890498728", "-2.125", "1.375", "-2",
    "2.5494898741957", "-0.375", "-0.125", "1.5",
    "0.0340681844957389", "0.375", "2.875", "-0.25",
    "2.59934940245681", "-1.125", "-1.14648350303967", "1.125",
    "-2.5", "-0.75", "-0.245308857462684", "-1.25",
)
ACTIVE = (6, 20, 28, 32, 36, 40, 44, 46, 50)
LINES = (0,0,0,15,0)
ENERGIES = ("522.398135900","534.208610459","544.705683199","841.107771361","533.193826337")


def experiment(root,fifth,columns,levels,*,return_chart=False):
    base = run_path(str(root/"ch39-cyclic-two-characters-nine.py"))
    return base["experiment"](root,fifth,columns,levels,return_chart=return_chart,
                              parameters=PARAMETERS,active=ACTIVE,lines=LINES,energies=ENERGIES,
                              uniform_grid=1,ratio_bound=".590",shared_source_error=True,
                              gravity_rows=(0,2,3,4,5,6,7,8,9),align_source=True,
                              gravity_readout=((1,".0480875",".0481259"),))


if __name__=="__main__":
    root=Path(__file__).resolve().parent
    joint=run_path(str(root/"ch39-fifth-full-positive-joint.py"))
    experiment(root,*joint["verify"](root,return_operator="spectral"))
