"""Analytical allocation/material obstruction; never allocates a C reserve."""
from __future__ import annotations
import argparse
import ast
import datetime
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
from packed_constructor import PackedStorage,PackedSource
HERE=Path(__file__).resolve().parent
assert sys.version_info>=(3,12) and __debug__

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def call(args):
    fds=tuple(int(x) for x in os.environ.get('WORKTREE_SCOPE_FDS','').split(',') if x)
    return subprocess.check_output(args,pass_fds=fds,text=True).strip()
def census(k):
    o=8109329+3*k;r=32574125+13*k;b=9546573383+16*k;w=(o+62).bit_length()
    return dict(O=o,R=r,B=b,w=w,S=b+w*r+o)
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);ap.add_argument('--theory-root',type=Path,required=True);args=ap.parse_args()
    s=PackedStorage();block=s.block([None]*4,bytearray(1),8);source=PackedSource(s,0,block)
    slot=s.state.slot(block._handle);ds=s.state.slot(source._handle)
    sizes={name:sys.getsizeof(value) for name,value in {'PackedBlock':block,'PackedSource':source,'PackedSlot':slot,'RefView':block.refs,'PackedRefs':slot.refs,'empty_bytearray':bytearray()}.items()}
    hb=sizes['PackedBlock']+sizes['PackedSlot']+sizes['RefView']+sizes['PackedRefs']+2*sizes['empty_bytearray']
    hd=sizes['PackedSource']+sizes['PackedSlot']+sizes['PackedRefs']+2*sizes['empty_bytearray']
    hc=2*hb+hd;rows=[]
    for name in ('positive-full','positive-partial','null-nine-distinct-left'):
        f=HERE/'constructor'/(name+'.input.bits');raw=f.read_bytes();L=len(raw.strip());m=32*(L+363);b=5*m//8;d=(13*m*m+m)//16+751;W=d+29*m//8+6
        old=HERE/('full-execution-'+name+'.json');old_result=json.loads(old.read_bytes());old_loader=old_result['stages']['converted']['required_complete_loader']
        assert (m,b,d,W)==tuple(old_loader[x] for x in ('B_input_bits','marked_input_bits','D_route','work_objects_at_complete_load'))
        rows.append(dict(case=name,input_bits=L,input_sha256=sha(f),converted_bits=m,marked_bits=b,D=d,W_loader_only=W,
            installer_at_K_equal_W=census(W),host_allocation_lower_bytes=W*hc,
            loader_comparison_source={'artifact':old.name,'sha256':sha(old),'scope':'prior measured loader prefix/count supplier; not rerun or full C'},
            sufficient_paid_reserve=False))
    # An AST audit supplies exact reuse boundaries for the unchanged consumers.
    parent='737f3e0b2ccb84e29f365be6d12ff519c6890582'
    checks={'native_engine.py':{'NativeEngine':['microtick','tick','run','resume','bw','census','check_census'],'Core':None},
       'native_install.py':{'BitKernel':None},'service_native.py':{'ServiceEngine':None,'OriginalConnection':['send','invoke','ack'],'ServiceCore':None},
       'session.py':{'OriginalPort':None},'lab.py':{'Buffer':None,'Lab':None}}
    audited=[]
    for file,classes in checks.items():
        old=ast.parse(call(['git','-C',str(HERE.parent),'show',parent+':rro-lab-execution/'+file]));new=ast.parse((HERE/file).read_bytes())
        for cls,names in classes.items():
            oc=next(x for x in old.body if isinstance(x,ast.ClassDef) and x.name==cls);nc=next(x for x in new.body if isinstance(x,ast.ClassDef) and x.name==cls)
            if names is None:assert ast.dump(oc)==ast.dump(nc)
            else:
                for name in names:
                    om=next(x for x in oc.body if isinstance(x,ast.FunctionDef) and x.name==name);nm=next(x for x in nc.body if isinstance(x,ast.FunctionDef) and x.name==name)
                    assert ast.dump(om)==ast.dump(nm),(file,cls,name)
            audited.append(dict(file=file,cls=cls,methods=names or 'entire class'))
    unchanged=('field_graph.py','primitive_engine.py','paid_code.py','field_decoder.py','field_io.py','field_service.py','framing.py','constructor_compiler.py','constructor_image.py','constructor_ir.py','evaluator_program.py','copier.py','nwmg.py','native_bits.c')
    for name in unchanged:
        old=call(['git','-C',str(HERE.parent),'show',parent+':rro-lab-execution/'+name])
        # Trailing LF is ignored only for this AST/text comparison; raw hashes remain in result.
        assert (HERE/name).read_text().strip()==old.strip(),name
    m=32*363;d=(13*m*m+m)//16+751
    free=shutil.disk_usage(HERE)
    result=dict(schema=1,scope='source-derived lower bounds and concrete Python-backend allocation obstruction; no full C run',
        sampled_at_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),python=sys.version,shallow_bytes=sizes,
        min_host_bytes_per_source_block=hb,min_host_bytes_per_descriptor=hd,min_host_bytes_per_installed_reserve_cell=hc,
        relaxed_L0_W=d+29*m//8+6,original_cases=rows,
        host={'RAM_bytes':int(call(['sysctl','-n','hw.memsize'])),'logical_cpus':int(call(['sysctl','-n','hw.logicalcpu'])),
              'memory_pressure':call(['memory_pressure']),'disk_total_bytes':free.total,'disk_free_bytes':free.free},
        paid_plan={'status':'incomplete','K200000_loads_no_accepted_payload':200000<d+29*m//8+6,
            'missing_sufficient_bound':['finite common-run reserve for acquisition/parser/input/cache/receipts/connection/provider buffers',
                'post-loader complete copier/B/LAB retained growth and pending scratch','returned-head decoder/output/cleanup and old jobs',
                'owner widths and same-cut overlap for all these retained objects','backed host reservation including all metadata and preparation overlap'],
            'precise_obstruction':'K>=W is necessary only. No upper bound on the full-session reserve is established. This backend additionally retains two proxy block records and one descriptor record per reserve cell; their nonshared shallow allocations alone exceed available RAM plus currently free disk for each original loader lower bound. OS compression/swap and larger future backing are not certified by this inequality.',
            'not_inferred':['full-run RSS/time from bounded probes','sufficient K from loader W','feasibility from installed RAM','full C completion from stack reference or component results']},
        unchanged_transition_AST=audited,unchanged_files=list(unchanged),reference_head=parent,
        current_source_hashes={p.name:sha(p) for p in sorted(HERE.glob('*.py'))},native_bits_sha256=sha(HERE/'native_bits.c'),
        theory_sha256=sha(args.theory_root/'docs/develop/theory/RECURSIVE_RELATIONAL_OBSERVATION_MONOTONE_LAB_INDEXED_EXECUTION.md'),
        C_released=False,C_complete=False,D_released=False,owned_child_processes_joined=True)
    for row in rows:
        assert row['host_allocation_lower_bytes']>result['host']['RAM_bytes']+free.free
    tmp=args.output.with_suffix('.json.tmp');tmp.write_text(json.dumps(result,sort_keys=True,indent=2)+'\n');tmp.replace(args.output)
    print(json.dumps({'result':str(args.output),'sha256':sha(args.output),'host_bytes_per_cell_lower':hc,'sufficient_plan':False}))
if __name__=='__main__':main()
