"""105 圈联合相位窗内，全部实际零点的正整数质量标尺均唯一。"""
from pathlib import Path
from runpy import run_path
import sys

from flint import arb


def verify(root, library=None):
    library = root if library is None else library
    window = run_path(str(library/'ch39-joint-zero-vacuum-flavor-density.py'))
    integer = run_path(str(library/'ch39-joint-integer-scale-pairs.py'))
    model, initial, joint, bands, density, delta = window['verify'](root)
    kmin, kmax = arb('31.9830'), arb('32.0170')
    assert delta.contains(0)
    integer['verify_scale_layers'](model, initial, joint, bands, kmin, kmax)
    print('整窗共有的全正整数解 (32,724)，零点相对下密度', density.str(18))


if __name__ == '__main__':
    root = Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).resolve().parent
    library = Path(sys.argv[2]) if len(sys.argv)>2 else root
    verify(root, library)
