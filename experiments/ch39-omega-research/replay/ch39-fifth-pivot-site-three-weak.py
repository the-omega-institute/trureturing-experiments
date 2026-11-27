"""在极点和三处枢轴站之外，任意加一站即可连续三壳消去弱差。"""
from pathlib import Path
from runpy import run_path

from flint import acb_mat, arb, ctx


def experiment(root, fifth, columns, levels):
    def lib(name):
        return run_path(str(root / ("ch39-" + name + ".py")))

    four = lib("four-pivot-sites-weak-horizon")
    graph = lib("three-ray-chiral-boundary")
    native = lib("native-rh-ray-pvm")
    gram = native["gram_cholesky"](root)
    sites = [graph["soft"]["pole"], *graph["boundary_sites"]()]
    assert sites == [14, 8, 7, 6]
    floor, heavy, top = levels
    ctx.prec = 256
    source = fifth * columns * acb_mat([[(heavy-floor).sqrt(), 0], [0, (top-floor).sqrt()]])
    spectra = []
    for j in range(3):
        _, basis = lib("adjacent-pvm-transport")["spectrum"](graph, 6, 2_306_888_100_092+j)
        rows = four["adj"](basis)
        weights = [z.real for z in native["readout"](basis, range(19), sites, gram[3], 3)]
        spectra.append((rows, rows*source, weights))

    planck, mz = arb("1.22089e19"), arb("91.1876")
    centers = [{0: "522.398136", 9: "178.455185", 15: "954.616970"},
               {9: "130.197637", 15: "895.563345"}, {15: "835.615125"}]
    width = arb(".00005")
    groups = [((0, 3, 4), (0, 15, 15)), ((11, 12, 18), (9, 15, 15)),
              ((1, 2, 5, 9, 10, 13, 15, 16, 17), (15, 9, 15))]
    assert sorted(i for group, _ in groups for i in group) == sorted(set(range(19))-set(sites))
    margin = arb(1)
    for group, selected in groups:
        targets = [arb(centers[j][i]) for j, i in enumerate(selected)]
        endpoints = [(center-width, center+width) for center in targets]
        for extra in group:
            keep = sites + [extra]
            local = acb_mat([source.tolist()[i] for i in keep])
            restricted = [(acb_mat([[row[i] for i in keep] for row in rows.tolist()]), old, weights)
                          for rows, old, weights in spectra]
            bound = four["certify_low_lines"](restricted, local, selected, endpoints, floor,
                                               planck, mz, native, lib)
            margin = min(margin, bound)
            print("额外站及六个弱差面的最小余量", extra, bound.str(12), flush=True)
    assert margin > 0


if __name__ == "__main__":
    root = Path(__file__).resolve().parent
    source = run_path(str(root / "ch39-fifth-full-positive-joint.py"))
    experiment(root, *source["verify"](root, return_operator="spectral"))
