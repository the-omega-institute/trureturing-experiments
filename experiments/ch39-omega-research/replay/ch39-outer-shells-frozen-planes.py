"""在两枚外侧壳上重选两荷谱权，检验训练用过的两个复三平面。"""
from pathlib import Path
from runpy import run_path
from flint import ctx


LEFT = (".0182", ".113", ("-.21668", "-.53516", ".71837", "-.38807"))
CENTRAL = (".00384", ".0384", (".15495", ".07140", ".97665", "-.13053"))


def individual_examples(root, emb):
    examples = (
        (-2, (".367160515", "-.245591947", "-.828036821", ".345300991"),
         ((1, 6, 12, 16), (1, 7, 10, 12))),
        (2, ("-.148712107", "-.078266232", "-.960849725", ".220288247"),
         ((6, 7, 10, 12), (6, 7, 10, 12))),
    )
    for offset, quaternion, chosen in examples:
        node = emb["CENTRAL"] + offset
        emb["verify_nodes"](root, ((node, "0"),), ".1",
                            (".00384", ".0384", quaternion), chosen)


def verify(root):
    emb = run_path(str(root / "ch39-common-embedding-central-right-ckm.py"))
    individual_examples(root, emb)
    cover = run_path(str(root / "ch39-neighbour-all-frame-support.py"))
    for offset in (-2, 2):
        for name, params in (("左侧", LEFT), ("中央与右侧", CENTRAL)):
            ctx.prec = 320
            _, _, _, rotation = emb["common_embedding"](*params)
            print("外侧编号与固定空间", offset, name, flush=True)
            cover["verify"](root, offset, rotation, params[:2])


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)
