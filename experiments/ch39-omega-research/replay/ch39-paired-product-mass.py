"""Two different closed-fermion parent masses after the joint two-product calibration."""
from fractions import Fraction as Q
from pathlib import Path
from runpy import run_path

root = Path(__file__).resolve().parent
m = run_path(str(root/"ch39-same-offset-square.py"))
base = m["base"]
data = base["data"]


def system(ids):
    pair = data("ch39-paired-products-data.json")
    prefix = data("ch39-tick-separators.json")["prefix"]+data("ch39-adjacent-offset-certificate.json")["tail"]+pair["tail"]
    return m["square_system"](ids,prefix,pair["square_floor"],paired=True)


def main():
    old = set(data("ch39-paired-products-data.json")["case"]["support"])
    for rec,(i,j) in zip(data("ch39-paired-product-mass-data.json"),((30,940),(175,27)),strict=True):
        ids = rec["support"]
        assert set(ids)==(old-{i}|{j})
        w = base["exact_fit"](rec,Q(1,100000),system(ids),50,4691)
        for mass,sign in zip(rec["bracket"],(1,-1),strict=True):
            s,n,_ = base["endpoint"](mass,ids,w)
            assert (s,n)==(sign,144)
        print("Paired observations, 50 positive heights; parent bracket in GeV",rec["bracket"])


if __name__=="__main__":
    main()
