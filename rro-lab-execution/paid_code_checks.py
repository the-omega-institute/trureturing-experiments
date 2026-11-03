"""Author contracts for the resumable installed-code decoder.

This component does not certify B's surrounding LAB/session controller.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import tempfile
from pathlib import Path

from paid_code import CodeStore, Decoder, decoder_program
from constructor_image import Image, build_image
from constructor_pack import unpack

HERE = Path(__file__).resolve().parent


def check_small(mutation=None):
    import struct
    from nwmg import compile_non_erasing
    from copier import Copier
    source = (((0, 'L', 0), (1, 'R', 1)),)
    raw = b'RROTABLE'+struct.pack('<II',2,1)+b''.join(
        struct.pack('<BbI',w,-1 if d=='L' else 1,n)
        for row in source for w,d,n in row)
    with tempfile.TemporaryDirectory(prefix='paid-code-') as td:
        td = Path(td)
        (td/'binary').write_bytes(raw)
        build_image(td/'binary',td)
        import gzip
        with gzip.open(td/'W_U.block.gz','rb') as f:
            image = f.read()
        store = CodeStore(image)
        ref = Image(_write(td/'raw',image))
        cpu = Decoder(store)
        independently_observed={'read':0,'write':0,'control':0,'refs':0}
        original_rb,original_wb=cpu._rb,cpu._wb
        original_control,original_ref=cpu._control,cpu._ref
        charge_start=dict(cpu.charges)
        def rb(register,index):
            independently_observed['read']+=1
            bit=original_rb(register,index)
            if mutation=='omit-register-reads':cpu.charges['register_bit_reads']-=1
            return bit
        def wb(register,index,value):
            independently_observed['write']+=1
            return original_wb(register,index,value)
        def control(name,value):
            independently_observed['control']+=cpu.CONTROL_WIDTHS[name]
            return original_control(name,value)
        def ref_write(name,value):
            independently_observed['refs']+=cpu.reference_width
            return original_ref(name,value)
        cpu._rb,cpu._wb,cpu._control,cpu._ref=rb,wb,control,ref_write
        expected = compile_non_erasing(Copier(source).rows)
        cases = [0,1,2,3,7,8,12,13,14,len(expected)-2,len(expected)-1,len(expected)]
        outcomes=[]
        for pc in cases:
            before=cpu.census()
            cpu.start(pc)
            saw_busy=False
            while not cpu.ready:
                cpu.tick()
                saw_busy |= cpu.phase != 'idle'
                if cpu.ticks_last % 1024 == 0:cpu.check_census()
            got=cpu.result()
            assert got == ref.fetch(pc) == (expected[pc] if pc<len(expected) else None)
            assert saw_busy
            assert cpu.census()['bits'] == before['bits']
            outcomes.append({'pc':pc,'instruction':got,'ticks':cpu.ticks_last})
        assert cpu.output_acknowledgements == len(cases)
        for quantity,coordinate in [('read','register_bit_reads'),('write','register_bit_writes'),
                                    ('control','control_bit_writes'),('refs','reference_register_write_bits')]:
            assert cpu.charges[coordinate]-charge_start.get(coordinate,0)==independently_observed[quantity], '[FAIL] '+coordinate+'_census'
        assert cpu.charges['directory_reference_read_bits']==cpu.charges['navigation_links']*cpu.reference_width
        assert cpu.charges['page_selection_tests']==15*cpu.charges['code_payload_bit_reads']
        assert cpu.charges['output_bits']==67*len(cases)
        # A successful reply cannot be returned from an unacknowledged bit stream.
        cpu.start(1)
        while cpu.phase != 'ack':
            cpu.tick()
        try:
            cpu.result()
        except RuntimeError:
            pass
        else:
            raise AssertionError('premature decoder return')
        cpu.tick()
        while not cpu.ready:
            cpu.tick()
        assert cpu.result()==expected[1]
        # Audit counters and transcript are never an instruction source.
        first=cpu.result()
        cpu.audit.clear()
        cpu.transcript.clear()
        cpu.start(1)
        cpu.run()
        assert cpu.result()==first
        # Counter corruption is detected by an independent field census.
        cpu.accounted_bits += 1
        try:
            cpu.check_census()
        except AssertionError:
            pass
        else:
            raise AssertionError('census mutation escaped')
        cpu.accounted_bits -= 1
        return {'cases':outcomes,'storage':cpu.census(),
                'store':store.census(),'charges':dict(sorted(cpu.charges.items())),
                'falsifiers':['premature-return','incremental-census-corruption',
                              'audit-erasure-preserves-output','per-primitive-independent-tariffs']}



def _write(path,data):
    path.write_bytes(data)
    return path


def check_full():
    import gzip
    with tempfile.TemporaryDirectory(prefix='paid-W-U-') as td:
        td=Path(td)
        unpack(HERE/'constructor-code','W_U.block.gz',td/'packed')
        with gzip.open(td/'packed','rb') as f:
            data=f.read()
        assert len(data)==590588102
        assert hashlib.sha256(data).hexdigest()=='6eceb93c62b1357ef2faf30debd798e615462dfa93dd90235ed1ed4282003a19'
        store=CodeStore(data)
        ref=Image(_write(td/'raw',data))
        cpu=Decoder(store)
        # Both destination kinds, local template addresses, final M and stopped PC.
        states=[0,1,ref.base(1)-1,ref.base(1),ref.base(ref.n//2),
                ref.base(ref.n-1),ref.nm-1]
        pcs=[13*q+j for q in states for j in (0,1,2,7,8,12)]
        pcs += [13*ref.nm,13*ref.nm+1]
        charges_before=dict(cpu.charges)
        outcomes=[]
        for pc in pcs:
            cpu.start(pc)
            cpu.run()
            got=cpu.result()
            assert got == ref.fetch(pc),(pc,got,ref.fetch(pc))
            cpu.check_census()
            outcomes.append({'pc':pc,'instruction':got,'ticks':cpu.ticks_last})
        return {'raw_bytes':len(data),'raw_sha256':hashlib.sha256(data).hexdigest(),
                'cases':outcomes,'store':store.census(),'storage':cpu.census(),
                'installation_charges':charges_before,
                'cumulative_charges':dict(sorted(cpu.charges.items())),
                'program':decoder_program().manifest(),
                'full_LAB_law_run':False,'B_complete':False}


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--full',action='store_true')
    ap.add_argument('--mutation',choices=['omit-register-reads'])
    ap.add_argument('--output',type=Path)
    ap.add_argument('--check',type=Path)
    a=ap.parse_args()
    result={'small':check_small(a.mutation)}
    if a.full:result['full']=check_full()
    result['source_hashes']={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest()
                             for n in ('paid_code.py','paid_code_checks.py')}
    result['A_interface_sha256']=hashlib.sha256((HERE/'constructor/INTERFACE.json').read_bytes()).hexdigest()
    result['B_complete']=False
    result['C_released']=False
    result['scope']='Code-decoder component and its recorded tariffs; no complete installed ML-I interpreter, LAB/ML-Q session or simultaneous source/observer account.'
    data=json.dumps(result,indent=2,sort_keys=True)+'\n'
    if a.output:a.output.write_text(data)
    if a.check:assert a.check.read_text()==data,'Saved paid-code result differs from exact rerun'
    print(json.dumps(result if not a.full else {
        'checks':'pass','full_cases':len(result['full']['cases']),
        'raw_bytes':result['full']['raw_bytes'],'B_complete':False},sort_keys=True))
    return 0

if __name__=='__main__':
    raise SystemExit(main())
