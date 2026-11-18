"""Reproducible post-loader necessary bound and conditional census.

No old or replacement full C run, large reserve, oracle, or law computation.
Expected invariant: actual adapter/copier records and provider census agree
with the symbolic formulas. Any differing count or retained field falsifies
that transfer. Final capacity remains a necessary bound, never sufficient.
"""
from __future__ import annotations
import ast,hashlib,json,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent.parent))
from copier import Copier,initial_word,decode_records
from nwmg import compile_non_erasing,encode
from constructor_compiler import initial_binary,binary_row
from lab import Arena,Obj
from session import prepare_source,OriginalPort
from service_native import OriginalConnection
from native_install import Block
from framing import parse_arithmetic
from packed_constructor import PackedStorage,PackedSource

HERE=Path(__file__).resolve().parent;OLD=HERE.parent
assert sys.version_info >= (3,12) and __debug__

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def census(k):
    o=8109329+3*k;r=32574125+13*k;b=9546573383+16*k;w=(o+62).bit_length()
    return dict(O=o,R=r,B=b,w=w,S=b+w*r+o)
def lower(a,s):
    q=16*(s+1)*(a+s+2)-16
    return dict(a=a,s=s,span=16*(s+1)*(a+s+2),q=q,K_min=1+q*(q+1)//2)
def source_census(k):
    d=(5*k*k+13*k+6)//2;v=8*k+8;h=3*k+3
    o=v+h+d;r=16+4*v+2*h+d;w=o.bit_length()
    return dict(O=o,R=r,B=d,w=w,S=d+w*r+o)
def connection_census(q,a):
    f=int(q>0);o=1+2*a+5*q+f;r=10+2*a+7*q+2*f;b=2*a+3*q
    w=max(1,o.bit_length());return dict(O=o,R=r,B=b,w=w,S=b+w*r+o)
def same(c,expected):
    assert (c['objects'],c['reference_slots'],c['bit_fields'],c['reference_width'],c['bits'])==tuple(expected[x] for x in ('O','R','B','w','S'))

def tiny_copies():
    rows=[]
    copier=Copier([((1,'R',0),(1,'L',0))]);code=compile_non_erasing(copier.rows)
    for word in ('','0','1','01','101'):
        inp=encode(initial_word(word));tape={i:int(b) for i,b in enumerate(inp) if int(b)}
        pc=h=steps=copies=0;high=len(inp)-1
        while copies<3:
            assert steps<150000
            op=code[pc];steps+=1
            if op=='M':tape[h]=1;pc+=1
            elif op=='L':h-=1;pc+=1
            elif op=='R':h+=1;pc+=1
            else:pc=op[1] if tape.get(h,0) else pc+1
            high=max(high,h)
            if pc%13==0 and pc//13 in copier.boundaries:
                copies+=1
                d=decode_records(lambda i:tape.get(2*i+1,0),h//2)
                v=lower(max(1,len(word)),copies)
                assert 2*(d['end']+8)==v['span'] and d['copies']==copies
                assert high+1>=v['q']
                rows.append(dict(word=word,copies=copies,span=v['span'],
                                 conservative_positions=v['q'],instantiated_positions=high+1,B_steps=steps))
    return rows

def adapter_checks():
    # Every nonterminal four-symbol transition takes exactly four binary rules.
    count=0
    for write in range(4):
        for move in (-1,0,1):
            four=[[(write,move,0)]*4]
            for read in range(4):
                tape={0:read//2,1:read%2};h=q=steps=0
                while steps==0 or q!=0:
                    assert steps<5
                    w,d,q=binary_row(four,q//11,q%11)[tape.get(h,0)]
                    tape[h]=w;h+=1 if d=='R' else -1;steps+=1
                assert steps==4 and h==2*move
                assert (tape.get(0,0),tape.get(1,0))==(write//2,write%2)
                count+=1
    return count

def invoke(connection,word):
    # Finite external peer fixture for the ORIGINAL provider, not an evaluator.
    head=tail=None
    for ch in word:
        cell=Block([None]*4,bytearray(1),8);cell.put(0,int(ch=='R'))
        if tail is None:head=cell
        else:tail.refs[1]=cell
        tail=cell
    header=Block([None]*4,bytearray(1),8);header.refs[1]=head
    target=Block([None]*4,bytearray(1),8);target.refs[1]=header
    connection.invoke(target);connection.ack(target)

def provider_checks():
    rows=[]
    for k in (0,1,4,60,135):
        source=prepare_source([0]*k);same(source.a.census(),source_census(k))
        port=OriginalPort(source);c=OriginalConnection(port)
        same(c.a.census(),connection_census(0,0))
        a=0
        for q,word in enumerate(('','LL','LL','R'*k+'LL','L'*130,'R'*(k+10)),1):
            invoke(c,word);a+=len(word)
            same(c.a.census(),connection_census(q,a))
            assert c.a.peak_bits<=connection_census(q,a)['S']+connection_census(q,a)['w']*connection_census(q,a)['R']
        assert source.root is port.root
        rows.append(dict(k=k,source=source_census(k),calls=q,A=a,connection=connection_census(q,a)))
    # Literal semantic code paths never allocate in the installed owner, apart
    # from the provider's explicitly separate connection Arena.
    alloc={}
    for filename,cls in [('native_engine.py','NativeEngine'),('service_native.py','ServiceEngine')]:
        tree=ast.parse((OLD/filename).read_bytes())
        node=next(n for n in tree.body if isinstance(n,ast.ClassDef) and n.name==cls)
        fn=next(n for n in node.body if isinstance(n,ast.FunctionDef) and n.name=='microtick')
        hits=[n.func.attr for n in ast.walk(fn) if isinstance(n,ast.Call) and isinstance(n.func,ast.Attribute)
              and n.func.attr in ('alloc','allocate','node','block')]
        # CPU/control helper names are not reserve allocations; unexpected hits
        # deliberately fail this bounded source-path check.
        assert not hits,(filename,hits);alloc[filename]=hits
    return dict(rows=rows,direct_allocation_calls=alloc,
                indirect_provider='ServiceEngine.sourcecall -> OriginalConnection.invoke -> OriginalPort.call',
                excluded_host_allocations='observer/output Python list appends are physical host allocations, outside this modeled-owner census')

def main():
    p=PackedStorage();b=p.block([None]*4,bytearray(1),8);d=PackedSource(p,0,b)
    slot=p.state.slot(b._handle)
    sizes={k:sys.getsizeof(v) for k,v in dict(PackedBlock=b,PackedSource=d,PackedSlot=slot,
              RefView=b.refs,PackedRefs=slot.refs,empty_bytearray=bytearray()).items()}
    hc=2*(sizes['PackedBlock']+sizes['PackedSlot']+sizes['RefView']+sizes['PackedRefs']+2*sizes['empty_bytearray'])
    hc+=sizes['PackedSource']+sizes['PackedSlot']+sizes['PackedRefs']+2*sizes['empty_bytearray']
    old=json.loads((OLD/'constructor/constructor-results.json').read_bytes());cases=[]
    for name in ('positive-full','positive-partial','null-nine-distinct-left'):
        path=OLD/'constructor'/(name+'.input.bits');raw=path.read_text();n,theta,cache,cyl=parse_arithmetic(raw)
        a=2*(len(raw)+362);assert len(initial_binary(raw,361))==a
        lo=lower(a,60*n);k=lo['K_min'];q=15*n+14 if name!='positive-partial' else 15*n+13
        bits=[dict(row).get(h,0) for row in cache for h in range(15)]
        schedule=['','LL','LRR','LLR','LL','LR','LLL','LRL']
        schedule+=['R'*i+'LL' for i in range(15*n) if (i%15,bits[i]) in cache[i//15]]
        schedule+=['LL','R'*15+'LL','R'*45+'LL','L'*130,'R'*(15*n+10),'LL']
        assert q==len(schedule);A=sum(map(len,schedule))
        sc=source_census(15*n);cc=connection_census(q,A)
        extra=sc['S']+cc['S']+max(sc['w']*sc['R'],cc['w']*cc['R'])
        ref=next(r for r in old['reference_runs'] if r['name']==name)
        assert ref['input_sha256']==sha(path)
        t=len(raw)+ref['reference']['pushes']
        # FS1's explicit persistent layout only, conditional on the old counted
        # StackProgram supplier. No inferred full native execution or timings.
        generation_sum=sum((i.bit_length()-1) for i in range(1,t+1)) if t<10000 else None
        # Closed formula for sum floor(log2 i), avoiding a run/large enumeration.
        ell=t.bit_length()-1
        g=sum(j*(1<<j) for j in range(ell))+ell*(t-(1<<ell)+1)
        if generation_sum is not None:assert generation_sum==g
        Dlo=3*g+2*t # each FIFO directory edge has depth >= generation.
        Olo=4+10*t+Dlo # b=0 relaxation, includes stack frames/queue, excludes register/code etc.
        Rlo=12+30*t+Dlo
        wlo=Olo.bit_length();bitslo=Olo+wlo*Rlo+Dlo
        cases.append(dict(case=name,N=n,L=len(raw),input_sha256=sha(path),lower=lo,
                    native_installer_at_necessary_K=census(k),host_shallow_lower_bytes=hc*k,
                    schedule=dict(q=q,A=A),source=sc,connection=cc,K_conditional_extra_bits=extra,
                    FS1_reference_conditional=dict(host_shallow_object_bytes=sys.getsizeof(Obj('FS1')),
                        explicit_host_object_lower_bytes=Olo*sys.getsizeof(Obj('FS1')),pushes=ref['reference']['pushes'],allocations=t,
                        stack_steps=ref['reference']['steps'],D_lower=Dlo,O_lower=Olo,R_lower=Rlo,
                        width_lower=wlo,retained_bits_lower=bitslo,
                        scope='conditional on prior counted stack supplier; not rerun or native work certificate')))
    files=['copier.py','nwmg.py','constructor_compiler.py','evaluator_program.py','constructor_ir.py',
           'lab.py','session.py','service_native.py','native_engine.py','native_install.py','packed_constructor.py',
           'constructor/constructor-results.json','constructor/evaluator.stack.bin']
    result=dict(schema=1,scope='necessary post-loader bounds and K-conditional modeled census; no full C run',
       old_C_complete=False,sufficient_K=False,host_reservation=False,complete_operation_vector=False,
       python=sys.version,host_shallow_sizes=sizes,host_shallow_bytes_per_reserve_cell=hc,
       adapter_cases=adapter_checks(),copier=tiny_copies(),provider=provider_checks(),cases=cases,
       census_scope='sequential same-cut dominating envelope; not attained peak, extraction bound, time or host bound',
       source_hashes={name:sha(OLD/name) for name in files},program_sha256=sha(Path(__file__)))
    out=HERE/'post-loader-results.json';out.write_text(json.dumps(result,sort_keys=True,indent=2)+'\n')
    print(json.dumps(dict(result=out.name,sha256=sha(out),exit=0,sufficient_K=False)))
if __name__=='__main__':main()
