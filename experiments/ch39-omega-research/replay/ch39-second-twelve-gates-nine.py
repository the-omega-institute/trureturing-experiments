"""以第二组独立短边门角重放三壳的九个精确端点。"""
from pathlib import Path
from runpy import run_path


PARAMETERS = (
    "-1.4349903468467", "2.75", "2.9154755206039", "-1.5520978771484",
    "-2.875", "-0.375", "-1.625", "2.3125",
    "0.3125", "1.125", "-2.4375", "0.5625",
    "0.3820798996223", "-2.8015400065673", "-2.3125", "-1.3096027541418",
    "-1.9851707485628", "1", "-0.76155775030621", "1.875",
    "1.375", "0.375", "1.75", "-1.4046933887995",
    "-2.2640348308338", "2.625", "2.6875", "1.4375",
    "0.9375", "0.3125", "-0.875", "-2.0164442068931",
    "-1.7467230661947", "-2.875", "1.7921607656711", "1.4375",
    "2.3984378778098", "-0.67028400535944", "2.1875", "-3.125",
    "0.9375", "-2.9375", "-1.7571023991474", "0.75",
    "2.1788893332621", "1.8929076724496", "0.9375", "-0.3125",
)
ACTIVE = (0, 2, 3, 12, 13, 15, 16, 18, 23, 24, 31, 32, 34, 36, 37, 42, 44, 45)
LINES = (0,0,15)
ENERGIES = ("522.398135900","534.208610459","835.615125179")


def experiment(root, fifth, columns, levels):
    base = run_path(str(root/"ch39-twelve-independent-gates-nine.py"))
    return base["experiment"](root,fifth,columns,levels,parameters=PARAMETERS,active=ACTIVE,
                              lines=LINES,energy_centers=ENERGIES,angle_radius=".0005",
                              ratio_limit=".653",contraction_limit=".333",scaled_widths=True)


if __name__=="__main__":
    root = Path(__file__).resolve().parent
    orig = run_path(str(root/"ch39-fifth-full-positive-joint.py"))
    experiment(root,*orig["verify"](root,return_operator="spectral"))
