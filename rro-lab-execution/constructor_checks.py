"""Author checks for the constructor, distinct from original full acceptance."""
from __future__ import annotations
import argparse,hashlib,json,struct,subprocess,tempfile
from pathlib import Path
from constructor_ir import Assembler,primitives,run_graph
from evaluator_program import build,encode_numbers
from framing import u,v,mask_word,cylinder_word,parse_arithmetic_output,arithmetic_input
from fractions import Fraction

HERE=Path(__file__).resolve().parent
OP={'halt':0,'goto':1,'push':2,'pop':3,'peek':4}
def graph_bytes(p):
    out=bytearray(b'RROSTACK'+struct.pack('<III',len(p['registers']),len(p['code']),p['start']))
    for x in p['code']:
        k=x[0]
        row=[OP[k],0,0,0,0,0]
        if k=='halt':row[2]=int(x[1])
        elif k=='goto':row[2]=x[1]
        else:row[1:len(x)]=x[1:]
        out.extend(struct.pack('<6I',*row))
    return bytes(out)

def decode_nat(r):return sum(b<<i for i,b in enumerate(reversed(r)))
def primitive_checks():
    from random import Random
    rng=Random(112)
    checked=0
    for name in ['add','sub','lt','eq','mul','div','gcd','reduce','qadd','encu','encq','readu','readv','get','put','popword','prepend']:
        a=Assembler();primitives(a)
        # One graph per tested routine, argument words actually received in stack0.
        arity=len(a.functions[name][0]);outputs=a.functions[name][1]
        def main(a):
            for i in range(arity):
                a.call('popword',['arg'+str(i),'input'],'input')
            a.call(name,['r.'+x for x in outputs],*['arg'+str(i) for i in range(arity)])
        p=a.compile(main)
        for rep in range(12):
            x,y=rng.randrange(128),rng.randrange(1,128)
            args=[x,y]
            if name=='reduce':args=[x,y]
            elif name=='qadd':args=[x,y,y,x+1]
            elif name in ['encu']:args=[x]
            elif name=='encq':args=[x,y]
            elif name=='readu':args=[list(map(int,u(x)))]
            elif name=='readv':args=[list(map(int,v('01'*rep)))]
            elif name in ['get','put']:args=[encode_numbers([0,7,12,3]),rep%4]+([x] if name=='put' else [])
            elif name=='popword':args=[encode_numbers([x,y])]
            elif name=='prepend':args=[x,encode_numbers([y])]
            def word(n):return ([] if n==0 else list(map(int,bin(n)[:1:-1]))) if isinstance(n,int) else n
            bits=''.join(''.join(map(str,[z for b in word(n) for z in (1,b)]+[0])) for n in args)
            r=run_graph(p,bits,20_000_000);assert r['success'],(name,args,r['state'])
            got={s:list(reversed(r['registers'][p['registers']['r.'+s]])) for s in outputs}
            ns={s:sum(b<<i for i,b in enumerate(w)) for s,w in got.items()}
            if name=='add':assert ns['r']==x+y
            elif name=='sub':assert ns['borrow']==int(x<y) and (x<y or ns['r']==x-y)
            elif name=='lt':assert ns['r']==int(x<y)
            elif name=='eq':assert ns['r']==int(x==y)
            elif name=='mul':assert ns['r']==x*y
            elif name=='div':assert (ns['q'],ns['r'])==divmod(x,y)
            elif name=='gcd':
                from math import gcd
                assert ns['r']==gcd(x,y)
            elif name in ['reduce','qadd']:
                expected=Fraction(x,y) if name=='reduce' else Fraction(x,y)+Fraction(y,x+1)
                assert (ns['n'],ns['d'])==(expected.numerator,expected.denominator)
            elif name=='encu':assert ''.join(map(str,got['r']))==u(x)
            elif name=='encq':
                from framing import rational
                assert ''.join(map(str,got['r']))==rational(Fraction(x,y))
            elif name=='readu':assert ns['n']==x and got['s']==[]
            elif name=='readv':assert got['w']==list(map(int,'01'*rep)) and got['s']==[]
            elif name=='get':assert ns['r']==[0,7,12,3][rep%4]
            elif name=='put':
                arr=[0,7,12,3];arr[rep%4]=x;assert got['seq']==encode_numbers(arr)
            elif name=='popword':assert ns['word']==x and got['seq']==encode_numbers([y])
            elif name=='prepend':assert got['seq']==encode_numbers([x,y])
            checked+=1
    return checked


def cases():
    complete=tuple(tuple((h,0) for h in range(15)) for _ in range(4))
    null=tuple(tuple((h,int(b)) for h,b in enumerate(format(j,'07b')+'00000000')) for j in range(9))
    partial=complete[:-1]+(tuple((h,0) for h in range(14)),)
    yield 'positive-full',arithmetic_input(4,Fraction(1,3),complete,[complete,None]), complete,[complete,None]
    yield 'positive-partial',arithmetic_input(4,Fraction(1,3),partial,[complete,None]),partial,[complete,None]
    yield 'null-nine-distinct-left',arithmetic_input(9,Fraction(1,3),null,[null,None]),null,[null,None]
    # Original valid unreduced prior bits, identical program and accepted domain.
    yield 'positive-unreduced-endpoint',u(4)+'0'+u(0)+u(6)+mask_word(complete)+u(0),complete,[]


def check_references(out,run_partial=False):
    out=Path(out);out.mkdir(exist_ok=True)
    p=build();raw=graph_bytes(p);(out/'evaluator.stack.bin').write_bytes(raw)
    (out/'evaluator.stack.json').write_text(json.dumps(p,separators=(',',':'))+'\n')
    completed=[]
    with tempfile.TemporaryDirectory(prefix='rro-constructor-') as td:
        exe=Path(td)/'reference';subprocess.run(['cc','-O3','-std=c11',str(HERE/'constructor_reference.c'),'-o',str(exe)],check=True)
        for name,bits,cache,cylinders in cases():
            if name=='positive-partial' and not run_partial:continue
            ip=out/(name+'.input.bits');op=out/(name+'.output.bits');ip.write_text(bits)
            r=subprocess.run([str(exe),str(out/'evaluator.stack.bin'),str(ip),str(op),'100000000000'],text=True,capture_output=True)
            if r.returncode:raise AssertionError((name,r.returncode,r.stdout,r.stderr))
            stats=json.loads(r.stdout);stats.pop('cpu_seconds')
            values,ratios=parse_arithmetic_output(op.read_text())
            if name.startswith('positive'):
                # Independent UJB edge argument: a fixed repeated full pair has
                # probability 1/(1944^(N-1)*128*256) in either source.
                z=Fraction(1,1944**3*128*256)
                if name=='positive-partial':
                    # The three full reads force one literal edge. The last
                    # partial right word permits its label or one other label.
                    # UJB source symmetry / independent <=4 supplier checks
                    # below supplies this case, not this repeated-pair formula.
                    from pins import suppliers
                    old,_,_=suppliers();from laws import Evaluator
                    expected=Evaluator(old,None).evaluate(4,Fraction(1,3),cache,cylinders)
                else:expected=([z,z,z]+([z,z,z,0,0,0] if cylinders else []),[Fraction(1),Fraction(0)] if cylinders else [])
            else:expected=([Fraction(0)]*9,None)
            assert (values,ratios)==expected,(name,values,ratios,expected)
            completed.append({'name':name,'input_bits':len(bits),'output_bits':op.stat().st_size,'input_sha256':hashlib.sha256(bits.encode()).hexdigest(),'output_sha256':hashlib.sha256(op.read_bytes()).hexdigest(),'values':[str(x) for x in values],'ratios':None if ratios is None else [str(x) for x in ratios],'reference':stats,'expected_source':'UJB16.2 repeated-edge/full-injection argument; null degree bound; partial <=4-edge inversion supplier'})
    return {'graph_states':len(p['code']),'registers':len(p['registers']),'graph_bytes':len(raw),'graph_sha256':hashlib.sha256(raw).hexdigest(),'reference_runs':completed,'reference_scope_complete':True,'acceptance_complete':False,'scope':'Primitive reference graph; no LAB execution or complete resource account.'}

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,default=HERE/'constructor');ap.add_argument('--partial',action='store_true');ap.add_argument('--primitives-only',action='store_true');args=ap.parse_args()
    n=primitive_checks();print('primitive cases',n,flush=True)
    if not args.primitives_only:
        result=check_references(args.output,args.partial);result['primitive_cases']=n;(args.output/'constructor-results.json').write_text(json.dumps(result,indent=2,sort_keys=True)+'\n');print(json.dumps(result),flush=True)
