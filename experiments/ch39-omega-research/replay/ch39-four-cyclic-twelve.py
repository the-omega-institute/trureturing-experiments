"""第七至十壳十二个端点的三循环等变共同共轭。"""
from runpy import run_path
from pathlib import Path

PARAMETERS = (
    "-0.25", "0.875", "-1.78805671755734", "-0",
    "0.625", "-0.25", "-0.625", "1.375",
    "0.75", "-0.375", "-0.75", "1.875",
    "-1.375", "0.375", "-0.125", "-0.375",
    "-0.979785019934898", "-0.875", "-2", "-0.375",
    "-1.125", "0", "0.454987521807512", "-2.125",
    "0.75", "-1.125", "-1", "-0.875",
    "0.60417940291988", "0.375", "1.63273913903233", "-0.125",
    "1.07363550852213", "1", "0.25", "1.125",
    "-1.375", "-1.375", "0.75", "-1.625",
    "-0.0685724458512146", "1.25", "-0.5", "1.25",
    "0.57569470991152", "0.125", "-1.5", "0.625",
    "-0.125", "-1.75", "-1.125", "1.125",
)
ACTIVE = (2, 16, 22, 28, 30, 32, 40, 44)
LINES = (0,0,15,0)
ENERGIES = ("522.398135900","534.208610459","835.615125179","544.755583188")


def experiment(root,fifth,columns,levels,*,return_chart=False):
    base = run_path(str(root/"ch39-cyclic-two-characters-nine.py"))
    return base["experiment"](root,fifth,columns,levels,return_chart=return_chart,
                              parameters=PARAMETERS,active=ACTIVE,lines=LINES,energies=ENERGIES,
                              uniform_grid=1,ratio_bound=".319",shared_source_error=True)


if __name__=="__main__":
    root=Path(__file__).resolve().parent
    joint=run_path(str(root/"ch39-fifth-full-positive-joint.py"))
    experiment(root,*joint["verify"](root,return_operator="spectral"))
