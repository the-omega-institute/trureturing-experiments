"""A DP alias falsifier and invocation portability checks for the constructor."""
import json,os,subprocess,sys,tempfile
from pathlib import Path
from constructor_ir import Assembler,primitives
from evaluator_program import define_evaluator,encode_numbers
from constructor_checks import graph_bytes
HERE=Path(__file__).resolve().parent

def main():
    # All four original slots, two proposed endpoint blocks, and two seven-bit
    # labels both forced to zero. Injective label assignment is impossible.
    def program(mutant):
        a=Assembler();primitives(a);define_evaluator(a)
        original=a.call
        def call(name,outs,*args):
            if mutant and a.context=='countlabels' and name=='get' and outs==['countlabels.term']:
                assert args[0]=='countlabels.old';args=('countlabels.new',args[1])
            return original(name,outs,*args)
        a.call=call
        def entry(a):
            data=[1 if h<7 else 0 for j in range(4) for h in range(15)]
            a.raw('p',encode_numbers([0,1,0,1]));a.raw('d',encode_numbers(data));a.op('countlabels','input','p','d',4,7,0,128)
        return a.compile(entry)
    descriptors=tuple(int(x) for x in os.environ.get('WORKTREE_SCOPE_FDS','').split(',') if x.isdigit())
    with tempfile.TemporaryDirectory(prefix='constructor path ',dir=HERE) as td:
        td=Path(td);exe=td/'reference machine'
        subprocess.run(['cc','-O3','-std=c11',str(HERE/'constructor_reference.c'),'-o',str(exe)],check=True,pass_fds=descriptors)
        inp=td/'input';inp.write_text('');runs=[]
        for mutant in (False,True):
            graph=td/'graph';out=td/'output';graph.write_bytes(graph_bytes(program(mutant)))
            result=subprocess.run([str(exe),str(graph),str(inp),str(out),'1000000000'],capture_output=True,text=True,check=True,pass_fds=descriptors)
            value=sum(int(b)<<j for j,b in enumerate(out.read_text()))
            assert value==(2 if mutant else 0),(mutant,value)
            stats=json.loads(result.stdout);stats.pop('cpu_seconds');runs.append({'old_array_replaced_by_new':mutant,'label_count':value,'literal_reference':stats})
        result=subprocess.run([sys.executable,'-B',str(HERE/'constructor_compiler.py'),'--checks-only'],cwd=td,env={'PATH':'/usr/bin:/bin:/usr/sbin:/sbin','WORKTREE_SCOPE_FDS':os.environ.get('WORKTREE_SCOPE_FDS','')},capture_output=True,text=True,pass_fds=descriptors)
        assert result.returncode==0,(result.returncode,result.stderr)
        portability={'empty_environment_except_declared_PATH_and_scope_FDs':True,'cwd_contains_spaces':True,'exit_code':result.returncode,'compiled_cases':json.loads(result.stdout)}
    r={'DP_old_new_alias_falsifier':runs,'portability':portability};(HERE/'constructor/regression-results.json').write_text(json.dumps(r,indent=2,sort_keys=True)+'\n');print(json.dumps(r,sort_keys=True))
if __name__=='__main__':main()
