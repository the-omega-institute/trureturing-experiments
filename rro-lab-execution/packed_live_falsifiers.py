"""Executed packed-field/code mutations; no evaluator or original C run."""
from __future__ import annotations
import argparse,hashlib,json,sys,time,resource
from pathlib import Path
from packed_constructor import PackedStorage,PackedBlock
from field_graph import Graph,Cell,OPS
from native_engine import NativeEngine

def main():
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);args=p.parse_args()
    assert sys.version_info >= (3,12) and __debug__
    start=time.perf_counter()
    g=Graph();g.label('again');g.emit('wait')
    g.get('picked','x',0);g.read('v','picked',0)
    g.emit('emit','v');g.emit('end');g.emit('ack');g.emit('return');g.j('again')
    program=g.seal()
    x=Cell([None]*4);y=Cell([None]*4,1);z=Cell([None]*4,0);x.refs[0]=y
    storage=PackedStorage();e=NativeEngine(program,1,bindings={'x':x},storage=storage,entry='again')
    e.cells.extend((x,y,z))
    try:
        e.seal();e.run()
        px,py,pz=(e.forward_observer[c] for c in (x,y,z))
        assert all(isinstance(c,PackedBlock) for c in e.cells)
        def frame():
            before=len(e.output);e.resume();e.run()
            out=e.output[before:];assert len(out)==1
            return out[0]
        baseline=frame();assert baseline==1
        px.refs[0]=pz
        aliased=frame();assert aliased==0
        pz.put(0,1)
        value=frame();assert value==1
        pz.put(0,0)
        row=e.code[next(i for i,n in enumerate(program.nodes) if n.op=='read')]
        original=row.bits;row.bits=(original&~31)|OPS.index('goto')
        # The changed encoded dispatch skips the actual zero read; the prior
        # retained bit is physically printed. A graph/reference fallback would
        # still read zero, so it fails this mutation witness.
        dispatch=frame();assert dispatch==1
        row.bits=original
        restored=frame();assert restored==0
        n=e.installation.accounted['objects'];e.installation.accounted['objects']=n-1
        detected=False
        try:e.census()
        except AssertionError:detected=True
        finally:e.installation.accounted['objects']=n
        assert detected;e.census()
        px.grant=0;grant_detected=False
        try:e.check_census()
        except AssertionError:grant_detected=True
        finally:px.grant=1
        assert grant_detected;e.check_census()
        e.charges['native_controller_selections']-=1
        assert e.native_steps!=e.charges['native_controller_selections']
        e.charges['native_controller_selections']+=1
        assert e.native_steps==sum(e.native_events.values())==e.charges['native_controller_selections']
        assert e.charges['output_ends']==e.charges['acknowledgements']==e.charges['returns']==5
        result={'schema':1,'checks':'pass','scope':'actual installed NativeEngine; five legal continued frames with intervening test mutations; no evaluator answer',
                'physical_output':e.output,'executed_mutations':['alias-successor','bit-value','encoded-dispatch'],
                'rejected_account_mutations':['census','grant','controller-charge'],'no_reference_fallback_witness':True,
                'native_steps':e.native_steps,'installation_steps':e.installation.steps,'modeled_census':e.census(),
                'storage':storage.metrics(),'seconds':time.perf_counter()-start,
                'host_ru_maxrss':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                'C_released':False,'C_complete':False,'D_released':False,
                'source_hashes':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(Path(__file__).parent.glob('*.py'))}}
        tmp=args.output.with_suffix('.json.tmp');tmp.write_text(json.dumps(result,sort_keys=True,indent=2)+'\n');tmp.replace(args.output)
        print(json.dumps({'checks':'pass','physical_output':e.output,'result':str(args.output),'sha256':hashlib.sha256(args.output.read_bytes()).hexdigest()}),flush=True)
    finally:e.close()
if __name__=='__main__':main()
