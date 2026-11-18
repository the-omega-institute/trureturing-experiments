"""Finite Boolean-stack assembler; numbers/arrays are data, never host opcodes.

The emitted graph has only StackProgram's five constructors. Python executes
assembly and the reference interpreter, not an arithmetic runtime callback.
"""
from __future__ import annotations

class Graph:
    def __init__(self):
        self.code=[]; self.labels={}; self.registers={'input':0}; self.serial=0
    def reg(self,s):
        if s not in self.registers: self.registers[s]=len(self.registers)
        return self.registers[s]
    def new(self,s='L'):
        self.serial+=1; return f'{s}:{self.serial}'
    def mark(self,s):
        assert s not in self.labels,s
        self.labels[s]=len(self.code)
    def emit(self,*x): self.code.append(x)
    def go(self,s): self.emit('goto',s)
    def push(self,r,b): self.emit('push',self.reg(r),int(b))
    def pop(self,r,e,z,o): self.emit('pop',self.reg(r),e,z,o)
    def peek(self,r,e,z,o): self.emit('peek',self.reg(r),e,z,o)
    def clear(self,r):
        l=self.new(); end=self.new(); self.mark(l); self.pop(r,end,l,l); self.mark(end)
    def move(self,s,d):
        assert s!=d
        l,z,o,end=[self.new() for _ in range(4)]
        self.mark(l); self.pop(s,end,z,o)
        self.mark(z); self.push(d,0); self.go(l)
        self.mark(o); self.push(d,1); self.go(l); self.mark(end)
    def copy(self,s,d):
        if s==d:return
        self.clear(d); self.clear('__copy'); self.move(s,'__copy')
        l,z,o,end=[self.new() for _ in range(4)]
        self.mark(l);self.pop('__copy',end,z,o)
        self.mark(z);self.push(s,0);self.push(d,0);self.go(l)
        self.mark(o);self.push(s,1);self.push(d,1);self.go(l);self.mark(end)
    def literal(self,r,bits):
        self.clear(r)
        for b in reversed(list(bits)):self.push(r,int(b))
    def const(self,r,n):
        assert isinstance(n,int) and n>=0
        self.literal(r,[] if n==0 else map(int,bin(n)[:1:-1]))
    def finish(self):
        out=[]
        for i,x in enumerate(self.code):
            k=x[0]
            if k=='halt':out.append(x)
            elif k=='goto':out.append((k,self.labels[x[1]]))
            elif k=='push':out.append((*x,i+1))
            else:out.append((k,x[1],*(self.labels[a] for a in x[2:])))
        assert all((x[0]=='halt' or all(0<=q<len(out) for q in (x[1:] if x[0]=='goto' else x[3:] if x[0]=='push' else x[2:]))) for x in out)
        return {'registers':self.registers,'labels':self.labels,'start':0,'code':out}

class Assembler:
    """Acyclic shared subroutines, with an explicit eight-bit return stack each.

    Calls copy actual words into private argument stacks. No procedure invokes
    itself, directly or indirectly. Return dispatch is a finite binary tree.
    """
    def __init__(self):
        self.g=Graph();self.functions={};self.calls={};self.context='main';self.control=[]
    def var(self,s):return self.context+'.'+s
    def set(self,d,x):
        if isinstance(x,int):self.g.const(d,x)
        else:self.g.copy(x,d)
    def raw(self,d,bits):self.g.literal(d,bits)
    def define(self,name,args,outs,body):self.functions[name]=(args,outs,body)
    def call(self,name,outs,*args):
        ins,res,_=self.functions[name]
        assert len(args)==len(ins) and len(outs)==len(res),(name,args,outs)
        for p,a in zip(ins,args):self.set(name+'.'+p,a)
        target=self.g.new('return.'+name); calls=self.calls.setdefault(name,[])
        number=len(calls);assert number<256,(name,number)
        calls.append(target)
        for bit in reversed([(number>>i)&1 for i in range(8)]):self.g.push(name+'.ret',bit)
        self.g.go('fn.'+name);self.g.mark(target)
        for p,o in zip(res,outs):self.g.copy(name+'.'+p,o)
    def op(self,name,d,*a):self.call(name,[d],*a)
    def ifzero(self,x):
        yes,no,end=[self.g.new() for _ in range(3)]
        self.g.peek(x,yes,no,no);self.g.mark(yes);self.control.append(('if',no,end))
    def otherwise(self):
        _,no,end=self.control[-1];self.g.go(end);self.g.mark(no)
        self.control[-1]=('else',None,end)
    def endif(self):
        k,no,end=self.control.pop()
        if k=='if':self.g.mark(no)
        self.g.mark(end)
    def loop(self):
        start,end=self.g.new(),self.g.new();self.g.mark(start)
        self.control.append(('loop',start,end))
    def breakzero(self,x):
        loop=next(v for v in reversed(self.control) if v[0]=='loop')
        nxt=self.g.new();self.g.peek(x,loop[2],nxt,nxt);self.g.mark(nxt)
    def repeat(self):
        k,start,end=self.control.pop();assert k=='loop';self.g.go(start);self.g.mark(end)
    def ret(self):self.g.go('ret.'+self.context)
    def compile(self,main):
        main(self);assert not self.control;self.g.emit('halt',True)
        # Compilation may discover calls from routine bodies; all are already declared.
        for name,(ins,outs,body) in self.functions.items():
            self.context=name;self.g.mark('fn.'+name);body(self);assert not self.control
            self.g.mark('ret.'+name);self.g.go('dispatch.'+name)
        for name,calls in self.calls.items():
            self.g.mark('dispatch.'+name)
            def tree(ids,depth):
                if depth==8:
                    self.g.go(calls[ids[0]] if ids else 'bad');return
                z,o=self.g.new(),self.g.new();self.g.pop(name+'.ret','bad',z,o)
                self.g.mark(z);tree([i for i in ids if not ((i>>depth)&1)],depth+1) if any(not ((i>>depth)&1) for i in ids) else self.g.go('bad')
                self.g.mark(o);tree([i for i in ids if (i>>depth)&1],depth+1) if any((i>>depth)&1 for i in ids) else self.g.go('bad')
            tree(list(range(len(calls))),0)
        for name in self.functions:
            if name not in self.calls:self.g.mark('dispatch.'+name);self.g.go('bad')
        self.g.mark('bad');self.g.emit('halt',False)
        return self.g.finish()

def norm(g,src,dst,tmp):
    """src is high-bit-first, trim high zeros, then reverse into dst."""
    l,z,o,end=[g.new() for _ in range(4)]
    g.mark(l);g.peek(src,end,z,o);g.mark(z);g.pop(src,'bad',l,l)
    g.mark(o);g.go(end);g.mark(end);g.clear(dst);g.move(src,dst)

def primitives(a):
    def arithmetic(a,subtract=False):
        g=a.g;p=a.context+'.';x,y,r,t=p+'a',p+'b',p+'r',p+'temp'
        g.clear(t)
        entries=[g.new(),g.new()];done=[g.new(),g.new()];g.go(entries[0])
        for carry in (0,1):
            g.mark(entries[carry]); xe,xz,xo=[g.new() for _ in range(3)];g.pop(x,xe,xz,xo)
            for xv,lab in [(None,xe),(0,xz),(1,xo)]:
                g.mark(lab);ye,yz,yo=[g.new() for _ in range(3)];g.pop(y,ye,yz,yo)
                for yv,yl in [(None,ye),(0,yz),(1,yo)]:
                    g.mark(yl)
                    if xv is None and yv is None:
                        if carry and not subtract:g.push(t,1)
                        g.go(done[carry]);continue
                    v=(xv or 0)-(yv or 0)-carry if subtract else (xv or 0)+(yv or 0)+carry
                    g.push(t,v%2);g.go(entries[int(v<0) if subtract else v//2])
        for c in (0,1):
            g.mark(done[c]);g.const(p+'borrow',c);g.go(p+'normalize')
        g.mark(p+'normalize');norm(g,t,r,p+'unused')
    a.define('add',['a','b'],['r'],lambda a:arithmetic(a))
    a.define('sub',['a','b'],['r','borrow'],lambda a:arithmetic(a,True))
    def cmp(a):
        v=a.var;a.call('sub',[v('delta'),v('r')],v('a'),v('b'))
    a.define('lt',['a','b'],['r'],cmp)
    def eq(a):
        v=a.var;a.call('sub',[v('delta'),v('less')],v('a'),v('b'))
        a.ifzero(v('less'));a.ifzero(v('delta'));a.set(v('r'),1);a.otherwise();a.set(v('r'),0);a.endif();a.otherwise();a.set(v('r'),0);a.endif()
    a.define('eq',['a','b'],['r'],eq)
    def mul(a):
        g=a.g;v=a.var;a.set(v('r'),0);a.loop();a.breakzero(v('b'))
        z,o,end=[g.new() for _ in range(3)];g.pop(v('b'),'bad',z,o)
        g.mark(o);a.op('add',v('r'),v('r'),v('a'));g.go(end);g.mark(z);g.mark(end)
        a.ifzero(v('a'));a.otherwise();g.push(v('a'),0);a.endif();a.repeat()
    a.define('mul',['a','b'],['r'],mul)
    def div(a):
        g=a.g;v=a.var;a.set(v('q'),0);a.set(v('r'),0);g.clear(v('rev'));g.move(v('a'),v('rev'))
        a.loop();a.breakzero(v('rev'))
        z,o,end=[g.new() for _ in range(3)];g.pop(v('rev'),'bad',z,o)
        g.mark(z);a.set(v('bit'),0);g.go(end);g.mark(o);a.set(v('bit'),1);g.mark(end)
        a.op('add',v('r'),v('r'),v('r'));a.op('add',v('r'),v('r'),v('bit'))
        a.op('lt',v('small'),v('r'),v('b'));a.ifzero(v('small'))
        a.call('sub',[v('r'),v('junk')],v('r'),v('b'));g.push(v('q'),1)
        a.otherwise();g.push(v('q'),0);a.endif();a.repeat()
        g.clear(v('tmp'));g.move(v('q'),v('tmp'));norm(g,v('tmp'),v('q'),v('unused'))
    a.define('div',['a','b'],['q','r'],div)
    def gcd(a):
        v=a.var;a.loop();a.breakzero(v('b'));a.call('div',[v('q'),v('rem')],v('a'),v('b'))
        a.set(v('a'),v('b'));a.set(v('b'),v('rem'));a.repeat();a.set(v('r'),v('a'))
    a.define('gcd',['a','b'],['r'],gcd)
    def reduce(a):
        v=a.var;a.op('gcd',v('g'),v('n'),v('d'));a.call('div',[v('n'),v('junk')],v('n'),v('g'));a.call('div',[v('d'),v('junk')],v('d'),v('g'))
    a.define('reduce',['n','d'],['n','d'],reduce)
    def qadd(a):
        v=a.var;a.op('mul',v('x'),v('an'),v('bd'));a.op('mul',v('y'),v('bn'),v('ad'))
        a.op('add',v('n'),v('x'),v('y'));a.op('mul',v('d'),v('ad'),v('bd'));a.call('reduce',[v('n'),v('d')],v('n'),v('d'))
    a.define('qadd',['an','ad','bn','bd'],['n','d'],qadd)
    def popword(a):
        g=a.g;v=a.var;g.clear(v('tmp'));g.clear(v('word'))
        l,z,o,end=[g.new() for _ in range(4)];g.mark(l);g.pop(v('seq'),'bad',end,o)
        g.mark(o);z1,o1=[g.new() for _ in range(2)];g.pop(v('seq'),'bad',z1,o1)
        g.mark(z1);g.push(v('tmp'),0);g.go(l);g.mark(o1);g.push(v('tmp'),1);g.go(l)
        g.mark(end);g.move(v('tmp'),v('word'))
    a.define('popword',['seq'],['word','seq'],popword)
    def prepend(a):
        g=a.g;v=a.var;g.clear(v('rev'));g.move(v('word'),v('rev'));g.push(v('seq'),0)
        l,z,o,end=[g.new() for _ in range(4)];g.mark(l);g.pop(v('rev'),end,z,o)
        g.mark(z);g.push(v('seq'),0);g.push(v('seq'),1);g.go(l)
        g.mark(o);g.push(v('seq'),1);g.push(v('seq'),1);g.go(l);g.mark(end)
    a.define('prepend',['word','seq'],['seq'],prepend)
    def get(a):
        v=a.var;a.loop();a.breakzero(v('i'));a.call('popword',[v('r'),v('seq')],v('seq'));a.call('sub',[v('i'),v('junk')],v('i'),1);a.repeat()
        a.call('popword',[v('r'),v('seq')],v('seq'))
    a.define('get',['seq','i'],['r'],get)
    def put(a):
        v=a.var;a.set(v('prefix'),0);a.loop();a.breakzero(v('i'))
        a.call('popword',[v('old'),v('seq')],v('seq'));a.op('prepend',v('prefix'),v('old'),v('prefix'));a.call('sub',[v('i'),v('junk')],v('i'),1);a.repeat()
        a.call('popword',[v('old'),v('seq')],v('seq'));a.op('prepend',v('seq'),v('word'),v('seq'))
        a.loop();a.breakzero(v('prefix'));a.call('popword',[v('old'),v('prefix')],v('prefix'));a.op('prepend',v('seq'),v('old'),v('seq'));a.repeat()
    a.define('put',['seq','i','word'],['seq'],put)
    def zeros(a):
        v=a.var;a.set(v('seq'),0);a.loop();a.breakzero(v('n'));a.g.push(v('seq'),0);a.call('sub',[v('n'),v('junk')],v('n'),1);a.repeat()
    a.define('zeros',['n'],['seq'],zeros)
    def concat(a):
        g=a.g;v=a.var;g.clear(v('tmp'));g.move(v('a'),v('tmp'));g.move(v('tmp'),v('b'))
    a.define('concat',['a','b'],['b'],concat)
    def encv(a):
        g=a.g;v=a.var;g.clear(v('len'));g.clear(v('tmp'))
        l,z,o,end=[g.new() for _ in range(4)];g.mark(l);g.pop(v('w'),end,z,o)
        g.mark(z);g.push(v('tmp'),0);g.push(v('len'),1);g.go(l)
        g.mark(o);g.push(v('tmp'),1);g.push(v('len'),1);g.go(l)
        g.mark(end);g.move(v('tmp'),v('w'));g.push(v('w'),0);g.move(v('len'),v('w'))
    a.define('encv',['w'],['w'],encv)
    def encu(a):
        g=a.g;v=a.var;g.clear(v('rev'));a.ifzero(v('n'));g.push(v('n'),0);a.endif()
        g.move(v('n'),v('rev'));a.op('encv',v('r'),v('rev'))
    a.define('encu',['n'],['r'],encu)
    def encq(a):
        v=a.var;a.call('reduce',[v('n'),v('d')],v('n'),v('d'));a.op('encu',v('x'),v('n'));a.op('encu',v('y'),v('d'));a.op('concat',v('r'),v('x'),v('y'));a.g.push(v('r'),0)
    a.define('encq',['n','d'],['r'],encq)
    def readv(a):
        g=a.g;v=a.var;g.clear(v('length'));g.clear(v('temp'));g.clear(v('w'))
        l,o,end=[g.new() for _ in range(3)];g.mark(l);g.pop(v('s'),'bad',end,o);g.mark(o);g.push(v('length'),1);g.go(l);g.mark(end)
        l,z,o,end=[g.new() for _ in range(4)];g.mark(l);g.pop(v('length'),end,z,z);g.mark(z)
        bz,bo=[g.new(),g.new()];g.pop(v('s'),'bad',bz,bo);g.mark(bz);g.push(v('temp'),0);g.go(l);g.mark(bo);g.push(v('temp'),1);g.go(l);g.mark(end);g.move(v('temp'),v('w'))
    a.define('readv',['s'],['w','s'],readv)
    def readu(a):
        v=a.var;a.call('readv',[v('w'),v('s')],v('s'));a.g.clear(v('n'));a.g.move(v('w'),v('n'))
        # Valid ML-F unsigned words are canonical MSB first. Zero -> empty.
        a.g.clear(v('tmp'));a.g.move(v('n'),v('tmp'));norm(a.g,v('tmp'),v('n'),v('unused'))
    a.define('readu',['s'],['n','s'],readu)


def run_graph(program,bits,fuel=1_000_000_000):
    """Literal counted primitive interpreter. Lists here implement stack storage."""
    regs=[[] for _ in program['registers']];regs[0]=list(map(int,reversed(bits)))
    code=program['code'];q=program['start'];steps=0;peak=len(bits);total=len(bits)
    while steps<fuel:
        x=code[q];k=x[0];steps+=1
        if k=='halt':return {'success':x[1],'steps':steps,'registers':regs,'state':q,'peak_stack_bits':peak}
        if k=='goto':q=x[1];continue
        r=regs[x[1]]
        if k=='push':r.append(x[2]);q=x[3];total+=1
        else:
            b=r[-1] if r else None
            if k=='pop' and r:r.pop();total-=1
            q=x[2 if b is None else 3+b]
        peak=max(peak,total)
    raise RuntimeError(('primitive fuel exhausted',steps,q))
