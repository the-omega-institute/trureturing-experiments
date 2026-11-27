"""按预登记规则检验第六壳实际荷六枢轴与原能标五式。"""
from pathlib import Path
from runpy import run_path

from flint import arb


def verify(root):
    number=2_306_888_100_091
    verify_exclusions(root,number)
    return verify_third(root,number)


def verify_exclusions(root,number):
    first,second=run_path(str(root/"ch39-fifth-arbitrary-pivot-obstructions.py"))["verify"](root,number=number)
    assert first[0]>arb("429.50")
    assert second[0]>arb("237.67") and second[1]>arb("3699.60")


def verify_third(root,number):
    anchors=(
        "0.00420342851020648430", "1.7276475178979798e-17", "2.7680672715375926e-15",
        "1.3816492240413636e-21", "1.4624879091739624e-31", "0.013824787504628776",
        "3.0002591443844813e-18", "0.1125340074708742", "5.775541216732063e-15",
        "3.91395875246121238e-34", "0.147821318648233560", "0.05957799752803994",
        "1.0426575215595841e-23", "1.005536791229539e-18", "3.670927472477619e-29",
        "6.93801761226990758e-10", "1.417615393727678e-11", "0.00149173655545344035",
        "2.73971088084214272e-20",
    )
    build=run_path(str(root/"ch39-fifth-arbitrary-third-pivot-five.py"))["verify"]
    return build(root,number=number,input_anchors=anchors,active=(18, 15, 17, 10, 0))


if __name__=="__main__":
    verify(Path(__file__).resolve().parent)
