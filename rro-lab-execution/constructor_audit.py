"""Full closed-image/encoding audit, with explicit reference-machine scope."""
from __future__ import annotations
import argparse,ast,gzip,hashlib,json,struct,tempfile
from pathlib import Path
from collections import Counter
from constructor_pack import digest,unpack
from constructor_checks import graph_bytes,primitive_checks
from evaluator_program import build,incidence
from constructor_compiler import generic_checks,initial_binary
from constructor_image import Image,BinaryTable,image_checks
from copier import initial_word
from nwmg import encode
from framing import parse_arithmetic,parse_arithmetic_output
from pins import suppliers
HERE=Path(__file__).resolve().parent


def audit(output,full=True):
    p=build();raw=graph_bytes(p);assert digest(HERE/'constructor/evaluator.stack.bin')=={'bytes':len(raw),'sha256':hashlib.sha256(raw).hexdigest()}
    assert json.loads((HERE/'constructor/evaluator.stack.json').read_text())==json.loads(json.dumps(p))
    assert set(x[0] for x in p['code'])=={'goto','push','pop','peek','halt'}
    n=len(p['code']);nr=len(p['registers'])
    for x in p['code']:
        if x[0]=='goto':assert 0<=x[1]<n
        if x[0]=='push':assert 0<=x[1]<nr and x[2] in (0,1) and 0<=x[3]<n
        if x[0] in ('pop','peek'):assert 0<=x[1]<nr and all(0<=q<n for q in x[2:])
    for path in HERE.glob('*.py'):ast.parse(path.read_bytes(),filename=str(path))
    checks={'graph_states':n,'graph_registers':nr,'graph_sha256':hashlib.sha256(raw).hexdigest(),'graph_operations':dict(Counter(x[0] for x in p['code'])),'primitive_cases':primitive_checks(),'compiler_cases':generic_checks(),'W_U_image_cases':image_checks()}
    # Source-independent correspondence of the public line/point constants.
    old,_,_=suppliers();literal=old.literal(0);pointsets=[tuple(sorted({old.RIGHT[t][0] for t in row})) for row in literal]
    assert len(set(pointsets))==36 and set(Counter(pointsets).values())=={2}
    from itertools import product
    pts=list(product(range(3),repeat=3));lines=[]
    for d in ((0,0,1),(1,0,0),(1,2,1),(1,1,1)):
        choices={line for line in pointsets if tuple(sorted(tuple((line[0][j]+k*d[j])%3 for j in range(3)) for k in range(3)))==line}
        assert len(choices)==9;lines.extend(sorted(choices,key=lambda line:line[0]))
    expected=[int(point in line) for line in lines for point in pts]
    assert incidence()==expected
    literal_degrees=[]
    for source in (0,1):
        rows=old.literal(source);ld=sorted(set(map(len,rows)));rd=sorted({sum(t in row for row in rows) for t in range(243)})
        assert ld==[27] and rd==[8];literal_degrees.append({'left':ld,'right':rd})
    checks['independent_literal_correspondence']={'entries':972,'ones':sum(expected),'type_capacities':[2,9],'literal_degrees':literal_degrees}
    references=json.loads((HERE/'constructor/constructor-results.json').read_text());runs=[]
    for r in references['reference_runs']:
        name=r['name'];ip=HERE/'constructor'/(name+'.input.bits');op=HERE/'constructor'/(name+'.output.bits');bits=ip.read_text()
        assert digest(ip)=={'bytes':r['input_bits'],'sha256':r['input_sha256']}
        assert digest(op)=={'bytes':r['output_bits'],'sha256':r['output_sha256']}
        assert r['reference']['success'];parse_arithmetic(bits);values,ratios=parse_arithmetic_output(op.read_text())
        assert [str(x) for x in values]==r['values']
        assert (None if ratios is None else [str(x) for x in ratios])==r['ratios']
        binary=initial_binary(bits,nr);copied=initial_word(binary);binput=encode(copied)
        for suffix,word in [('source-binary',binary),('copier-input',copied),('B-input',binput)]:
            path=HERE/'constructor'/(name+'.'+suffix+'.bits');path.write_text(word)
        runs.append({'name':name,'accepted_ML_F_bits':len(bits),'four_input_symbols':len(binary)//2,'binary_input_bits':len(binary),'non_erasing_initial_input_bits':len(copied),'B_input_bits':len(binput),'primitive_reference_steps':r['reference']['steps'],'output_bits':r['output_bits'],'reference_machine':'Literal five-op finite Boolean-stack graph; complete full-law run.','full_literal_binary_run':False})
    checks['input_output_encodings']=runs
    if full:
        with tempfile.TemporaryDirectory(prefix='constructor-audit-',dir=HERE) as td:
            td=Path(td);bp=td/'binary';ip=td/'W.gz';flat=td/'image'
            unpack(HERE/'constructor-tables','evaluator.binary.bin',bp)
            unpack(HERE/'constructor-code','W_U.block.gz',ip)
            with gzip.open(ip,'rb') as f,flat.open('wb') as out:
                import shutil
                shutil.copyfileobj(f,out)
            expected_image=json.loads((HERE/'constructor-code/image-results.json').read_text())
            assert digest(flat)=={'bytes':expected_image['uncompressed_image_bytes'],'sha256':expected_image['uncompressed_image_sha256']}
            im=Image(flat);source=BinaryTable(bp);assert source.n==im.n
            # All bindings and bases, checked against the materialized table.
            for k,(at,size) in enumerate(im.ts):
                for local in range(size):
                    for bit in (0,1):
                        w,d,dest=struct.unpack_from('<Bbq',im.data,at+20*local+10*bit)
                        assert w in (0,1) and w>=bit and d in (-1,1) and (0<=dest<size or dest in (-1,-2))
            for q in range(source.n):
                k,d0,d1=struct.unpack_from('<BQQ',im.data,im.binding_start+17*q)
                x,y=source.row(q);assert source.kind(q)==k
                assert im.base(q+1)-im.base(q)==im.ts[k][1]
                assert d0==im.base(x[2]) and d1==im.base(y[2])
            assert im.base(0)==0 and im.base(im.n)==im.nm
            checks['full_image_closure']={'resolved_bindings':source.n,'template_transition_cases':sum(2*n for _,n in im.ts),'non_erasing_states_covered':im.nm,'non_erasing_rules_covered':2*im.nm,'W_U_instructions_covered':13*im.nm+1,'coverage_method':'Every resolved source binding and every literal template transition, with checked contiguous bases; Cartesian instantiation covers all expanded rows.','flat_W_U_materialized':False}
    checks['scope']='Author finite checks and ordinary all-input proof in theory §§10–13; no Lean, independent review, paid runtime or complete LAB execution.'
    output.write_text(json.dumps(checks,indent=2,sort_keys=True)+'\n');return checks

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,default=HERE/'constructor/audit-results.json');ap.add_argument('--small',action='store_true');a=ap.parse_args()
    print(json.dumps(audit(a.output,not a.small),sort_keys=True),flush=True)
