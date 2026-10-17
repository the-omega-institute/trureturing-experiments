"""Once: exact fixed-basis enclosure and complete signed inertia of original112."""
import datetime
import hashlib
import importlib.util
import json
from pathlib import Path
import resource
import sys
import time
R=Path('/Users/macstudio/omega/research/rh-01a124c9-coordination/signed-source-envelope-20261014')
OUT=R/'signed-eta1-preconditioned112-result-r2512.json'
assert not OUT.exists(), 'Join existing incarnation.'

def load(path,digest,name):
    assert hashlib.sha256(path.read_bytes()).hexdigest()==digest
    spec=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(spec);sys.modules[name]=module
    spec.loader.exec_module(module);return module

core=load(R/'r2380/diagonal_bound_query.py','b278f8a358845ca436166b5595f8043173c7e09678d805d9d2c7e54d9d83d1c7','r2512core')
base=load(R/'r2429/signed_interval_replay.py','8cdfce25b0342e0b0f5b208f186aef1633e69209383a178e42311b76c689256e','r2512base')
small=load(R/'r2490/small_boundary_inertia.py','7f018c43946985a95a15ff5223c49c423ee14a9818c162f9964f693fc72bf7e0','r2512small')
new=load(R/'r2510/fixed_congruence.py','d436cc16ffb63536a04d04ca2e586fdcca8af4b11a93e4f00de91f7c7180e24e','r2512new')
review=json.loads((R/'r2511/root-reviewed-admission.json').read_text())
assert review['admitted'] and review['candidate_sha256']=='d436cc16ffb63536a04d04ca2e586fdcca8af4b11a93e4f00de91f7c7180e24e'
raw=(R/'signed-eta1-border-reduction-result-r2499.json').read_bytes();sha=hashlib.sha256(raw).hexdigest()
assert sha==new.REDUCTION_SHA256
packet=json.loads(raw)
started=time.monotonic_ns();phase='basis proposal'
basis_raw=None;enclosure_raw=None
try:
    basis=new.propose_basis(raw);order=basis['permutation'];blocks=basis['original_label_blocks']
    assert basis['basis_grid_bits']==new.BITS and basis['original_arithmetic_bits']==1024 and basis['dimension']==112
    assert sorted(order)==list(range(112)) and [i for b in blocks for i in b]==order
    small._check_blocks(core,order,blocks,complete=True)
    basis_raw=(json.dumps(basis,indent=2)+'\n').encode()
    with (R/'r2510/actual-fixed-dyadic-basis.json').open('xb') as stream:stream.write(basis_raw)
    phase='one-shot fixed exact congruence enclosure'
    original=small._decode_upper(core,packet['reduced112_upper'],112,1024,'original112')
    permuted=[[original[i][j] for j in order] for i in order]
    matrix=new.enclose_congruence(core,base,permuted,basis['basis_upper_point_integers'],new.BITS)
    enclosure=dict(schema='R2512_FIXED_DYADIC_CONGRUENCE_ENCLOSURE_V1',
        dimension=112,arithmetic_bits=new.BITS,original_arithmetic_bits=1024,
        reduction_result_sha256=sha,basis_sha256=hashlib.sha256(basis_raw).hexdigest(),
        preconditioned_coordinate_labels=order,original_label_blocks_as_coordinate_names=blocks,
        symmetric_encoding=small.ENCODING,matrix_upper=small._export_upper(matrix),
        no_inertia_yet=True,MAIN_or_RH_claim=False)
    enclosure_raw=(json.dumps(enclosure,indent=2)+'\n').encode()
    with (R/'r2510/actual-congruent112-enclosure.json').open('xb') as stream:stream.write(enclosure_raw)
    phase='complete outward signed inertia'
    _,active,inertia,transcript=small._eliminate(core,base,matrix,order,blocks)
    core.require(not active and sum(inertia)==112,'incomplete_coverage','complete transformed112 count')
    core.require(inertia[1]>=56,'inherited_contract_inconsistency','restriction negative index cannot be negative')
    result=dict(schema='R2512_COMPLETE_FIXED_CONGRUENCE_SIGNED_INERTIA_V1',
        status='complete_signed_inertia_conditional_on_inherited_source_contract',
        eta=['1','1'],original_arithmetic_bits=1024,congruence_arithmetic_bits=new.BITS,
        original_reduced_inertia=inertia,border_inertia=[inertia[0],inertia[1]+21,0],
        same_source_restricted_inertia=[1205+inertia[0],inertia[1]-56,0],
        same_source_restriction_positive_definite=inertia==[56,56,0],
        transformed_coordinate_labels=order,pivot_transcript=transcript,
        input_bindings=packet['input_bindings'],padding_account=packet['padding_account'],
        premises=packet['premises']+['The stored exact triangular dyadic basis has positive power-of-two diagonal, so congruence preserves inertia.','All112 preconditioned coordinates were covered by outward signed pivots; midpoint signs are not certificate evidence.'],
        restricted_domination_claim=False,eta_budget_membership_claim=False,
        tau_or_z_upper_claim=False,MAIN_or_RH_claim=False,compiled_or_official_acceptance=False)
    code=0
except Exception as exc:
    known=isinstance(exc,core.BoundaryError);code=2 if known else 1
    result=dict(schema='R2512_FIXED_CONGRUENCE_SIGNED_CERTIFICATE_FAILURE_V1',status='failed',
        failure_phase=phase,failure_code=exc.code if known else type(exc).__name__,
        failure_detail=exc.detail if known else str(exc),original_reduced_inertia=None,
        border_inertia=None,eta_budget_membership_claim=False,MAIN_or_RH_claim=False,
        compiled_or_official_acceptance=False)
result['certificate_bindings']=dict(reduction_result_sha256=sha,
    basis_sha256=hashlib.sha256(basis_raw).hexdigest() if basis_raw else None,
    congruent_enclosure_sha256=hashlib.sha256(enclosure_raw).hexdigest() if enclosure_raw else None,
    congruence_source_sha256=hashlib.sha256(Path(new.__file__).read_bytes()).hexdigest(),
    signed_kernel_source_sha256=hashlib.sha256(Path(small.__file__).read_bytes()).hexdigest())
result['process_receipt']=dict(at=datetime.datetime.now(datetime.timezone.utc).isoformat(),
    intended_process_exit=code,elapsed_ns=time.monotonic_ns()-started,
    maximum_resident_set_size_native=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
    maximum_resident_set_size_unit='bytes' if sys.platform=='darwin' else 'KiB',
    execution_phase=phase)
with OUT.open('x') as stream:json.dump(result,stream,indent=2,sort_keys=True);stream.write('\n')
print(json.dumps(dict(result=str(OUT),status=result['status'],intended_exit=code,process_receipt=result['process_receipt'])))
raise SystemExit(code)
