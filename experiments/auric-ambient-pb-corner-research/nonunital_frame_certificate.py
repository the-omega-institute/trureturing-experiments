"""Adapt the pinned 48-frame source to full H4 input, retaining exact arrays."""
from pathlib import Path
import hashlib
import json
import sys

ROOT = Path(__file__).resolve().parent
REPOSITORY = ROOT.parents[1]
COMMIT = '651b00b5abcfc10de377ae180bcfde38752bdffd'
SOURCE = 'experiments/auric-shared-d-frame-central-pb/frame_kernel.py'
EXPECTED = '8f6224dc5c64959c57d375f82f10f6f6605eec11232b86c0614ae9261b6d0331'
assert sys.version_info >= (3, 12)
raw = (REPOSITORY / SOURCE).read_bytes()
assert hashlib.sha256(raw).hexdigest() == EXPECTED
src = raw.decode()

def replace_once(old, new):
    global src
    assert src.count(old) == 1, old
    src = src.replace(old, new)

# Basis for full input: the original fifteen traceless matrices, followed by I4.
# Output remains the same fifteen-dimensional traceless quotient.
replace_once('wedges=[(i,j) for i in range(15) for j in range(i+1,15)]',
             'wedges=[(i,j) for i in range(16) for j in range(i+1,16)]\n'
             'def input_coord(a):\n'
             '    return np.append(coord(a), int(np.trace(a))*inv4%P)')
replace_once('a=coord(u@c%P@us%P)', 'a=input_coord(u@c%P@us%P)')
replace_once('b=coord(u@t%P@us%P)', 'b=input_coord(u@t%P@us%P)')
replace_once("'nullity':1575-r", "'nullity':1800-r")
replace_once('assert r==1575 and rt==1575 and det!=0 and wrank==105',
             'assert r==1800 and rt==1800 and det!=0 and wrank==120')
saved_argv = sys.argv
sys.argv = [SOURCE, '48']
state = {}
try:
    exec(compile(src, f'{COMMIT}:{SOURCE} [full-input adaptation]', 'exec'), state)
finally:
    sys.argv = saved_argv

np = state['np']
arrays = {
    'nonunital-matrix.u16le': state['mat'],
    'nonunital-pivot-rows.u16le': state['pivots'],
    'nonunital-input-wedges.u16le': state['inputs'],
    'nonunital-denominators.u16le': state['denominators'],
}
records = {}
for name, values in arrays.items():
    array = np.asarray(values, dtype='<u2')
    payload = array.tobytes(order='C')
    (ROOT / name).write_bytes(payload)
    records[name] = {'shape': list(array.shape), 'bytes': len(payload),
                     'sha256': hashlib.sha256(payload).hexdigest()}
result = {
    'status': 'exact modular certificate; research calculation, not independent verification or Lean',
    'source_commit': COMMIT, 'source_path': SOURCE, 'source_sha256': EXPECTED,
    'adapter_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    'input_basis': 'B_0,...,B_14,I4 (original traceless basis unchanged)',
    'input_coord': '(chi(a),Tr(a)/4)',
    'input_wedges': 'lexicographic i<j<16; no 1/2',
    'output_basis': 'original B_0,...,B_14; scalar output removed',
    'row_order': '(frame, elementary leg, offdiagonal output coordinate), lexicographic',
    'column_order': '(input wedge, output basis), output basis fastest',
    'prime': state['P'], 'i_image': state['II'], 'frames': 48,
    'rank_mod_1009': state['r'], 'quotient_nullity': 1800-state['r'],
    'maximal_minor_mod_1009': state['det'],
    'minor_rows': 'leftmost pivot columns of RREF(matrix transpose), increasing',
    'minor_columns': 'all 1800 columns in original order',
    'input_wedge_rank_mod_1009': state['wrank'],
    'full_map_dimension': 1920, 'exact_central_kernel_dimension': 120,
    'denominator_minmax': [min(state['denominators']), max(state['denominators'])],
    'rational_lift': 'Z[1/(2 product_j d_j)]; all d_j between 5 and 26; 1009 prime',
    'array_encoding': 'unsigned 16-bit little-endian, C row-major, zero-based indices',
    'arrays': records,
}
(ROOT / 'nonunital-certificate.json').write_text(json.dumps(result, indent=2, sort_keys=True)+'\n')
print(json.dumps(result, sort_keys=True))
