"""The one fixed indexed-law program. All routines assemble to Boolean-stack ops.

No Fraction, supplier import, input-specific compilation or host runtime law
operation occurs here. The only assembly-time table is public UJB15.2 geometry.
"""
from constructor_ir import Assembler,primitives
from itertools import product


def fori(a,name,bound,body,start=0):
    i=a.var(name); test=a.var(name+'.test');a.set(i,start);a.loop();a.op('lt',test,i,bound);a.breakzero(test);body(i);a.op('add',i,i,1);a.repeat()

def dec(a,x):a.call('sub',[x,a.var('discard')],x,1)

def bit(a,src,dst):
    g=a.g;z,o,end=[g.new() for _ in range(3)];g.pop(src,'bad',z,o)
    g.mark(z);a.set(dst,0);g.go(end);g.mark(o);a.set(dst,1);g.mark(end)

def encode_numbers(items):
    out=[]
    for n in items:
        for b in ([] if not n else map(int,bin(n)[:1:-1])):out.extend((1,b))
        out.append(0)
    return out

def incidence():
    pts=list(product(range(3),repeat=3));lines=[]
    for d in ((0,0,1),(1,0,0),(1,2,1),(1,1,1)):
        choices={tuple(sorted(tuple((p[j]+k*d[j])%3 for j in range(3)) for k in range(3))) for p in pts}
        lines.extend(sorted(choices,key=lambda line:line[0]))
    assert len(lines)==36 and sum(p in line for line in lines for p in pts)==108
    return [int(p in line) for line in lines for p in pts]


def define_evaluator(a):
    def readcache(a):
        v=a.var;a.op('mul',v('size'),v('n'),15);a.op('zeros',v('d'),v('size'))
        def mask(i):
            bit(a,v('s'),v('value'));a.op('put',v('d'),v('d'),i,v('value'))
        fori(a,'i',v('size'),mask)
        def value(i):
            a.op('get',v('value'),v('d'),i);a.ifzero(v('value'));a.otherwise();bit(a,v('s'),v('value'));a.op('add',v('value'),v('value'),1);a.op('put',v('d'),v('d'),i,v('value'));a.endif()
        fori(a,'i',v('size'),value)
    a.define('readcache',['s','n'],['d','s'],readcache)
    def conjoin(a):
        v=a.var;a.set(v('ok'),1);a.op('mul',v('size'),v('n'),15)
        def pos(i):
            a.op('get',v('x'),v('d'),i);a.op('get',v('y'),v('e'),i)
            a.ifzero(v('x'));a.op('put',v('d'),v('d'),i,v('y'))
            a.otherwise();a.ifzero(v('y'));a.otherwise();a.op('eq',v('same'),v('x'),v('y'));a.ifzero(v('same'));a.set(v('ok'),0);a.endif();a.endif();a.endif()
        fori(a,'i',v('size'),pos)
    a.define('conjoin',['d','e','n'],['d','ok'],conjoin)
    def blocks(a):
        v=a.var;a.set(v('r'),0)
        def pos(i):
            a.op('get',v('b'),v('p'),i);a.op('add',v('b'),v('b'),1);a.op('lt',v('less'),v('r'),v('b'));a.ifzero(v('less'));a.otherwise();a.set(v('r'),v('b'));a.endif()
        fori(a,'i',v('n'),pos)
    a.define('blocks',['p','n'],['r'],blocks)
    def nextpart(a):
        v=a.var;a.set(v('more'),0);a.set(v('i'),v('n'));dec(a,v('i'));a.loop();a.breakzero(v('i'));a.set(v('m'),0)
        def prefix(j):
            a.op('get',v('x'),v('p'),j);a.op('lt',v('less'),v('m'),v('x'));a.ifzero(v('less'));a.otherwise();a.set(v('m'),v('x'));a.endif()
        fori(a,'j',v('i'),prefix);a.op('add',v('m'),v('m'),1);a.op('get',v('x'),v('p'),v('i'));a.op('lt',v('less'),v('x'),v('m'))
        a.ifzero(v('less'));a.otherwise();a.op('add',v('x'),v('x'),1);a.op('put',v('p'),v('p'),v('i'),v('x'));a.set(v('more'),1)
        a.op('add',v('start'),v('i'),1)
        fori(a,'j',v('n'),lambda j:a.op('put',v('p'),v('p'),j,0),v('start'));a.ret();a.endif();dec(a,v('i'));a.repeat()
    a.define('nextpart',['p','n'],['p','more'],nextpart)
    def side(a):
        # Fully observed labels force the sole possible equality partition.
        # Otherwise return the initial RGS; countlabels handles the general case.
        v=a.var;a.set(v('full'),1);a.op('zeros',v('p'),v('n'));a.set(v('r'),0)
        def slot(j):
            def coord(h):
                a.op('mul',v('index'),j,15);a.op('add',v('index'),v('index'),v('off'));a.op('add',v('index'),v('index'),h);a.op('get',v('x'),v('d'),v('index'));a.ifzero(v('x'));a.set(v('full'),0);a.endif()
            fori(a,'h',v('width'),coord)
        fori(a,'j',v('n'),slot)
        a.ifzero(v('full'));a.ret();a.endif()
        def assign(j):
            a.set(v('found'),0);a.set(v('chosen'),0)
            def previous(k):
                a.set(v('same'),1)
                def coord(h):
                    a.op('mul',v('index'),j,15);a.op('add',v('index'),v('index'),v('off'));a.op('add',v('index'),v('index'),h);a.op('get',v('x'),v('d'),v('index'))
                    a.op('mul',v('index'),k,15);a.op('add',v('index'),v('index'),v('off'));a.op('add',v('index'),v('index'),h);a.op('get',v('y'),v('d'),v('index'));a.op('eq',v('test'),v('x'),v('y'));a.ifzero(v('test'));a.set(v('same'),0);a.endif()
                fori(a,'h',v('width'),coord)
                a.ifzero(v('same'));a.otherwise();a.op('get',v('chosen'),v('p'),k);a.set(v('found'),1);a.endif()
            fori(a,'k',j,previous)
            a.ifzero(v('found'));a.set(v('chosen'),v('r'));a.op('add',v('r'),v('r'),1);a.endif();a.op('put',v('p'),v('p'),j,v('chosen'))
        fori(a,'j',v('n'),assign)
    a.define('side',['d','n','width','off'],['p','full'],side)
    def countlabels(a):
        v=a.var;a.call('blocks',[v('r')],v('p'),v('n'));a.set(v('count'),0)
        a.op('mul',v('size'),v('r'),v('width'));a.op('zeros',v('requirements'),v('size'));a.set(v('ok'),1)
        def sample(j):
            a.op('get',v('block'),v('p'),j)
            def coord(h):
                a.op('mul',v('di'),j,15);a.op('add',v('di'),v('di'),v('off'));a.op('add',v('di'),v('di'),h);a.op('get',v('value'),v('d'),v('di'))
                a.ifzero(v('value'));a.otherwise();a.op('mul',v('ri'),v('block'),v('width'));a.op('add',v('ri'),v('ri'),h);a.op('get',v('oldvalue'),v('requirements'),v('ri'))
                a.ifzero(v('oldvalue'));a.op('put',v('requirements'),v('requirements'),v('ri'),v('value'))
                a.otherwise();a.op('eq',v('same'),v('oldvalue'),v('value'));a.ifzero(v('same'));a.set(v('ok'),0);a.endif();a.endif();a.endif()
            fori(a,'h',v('width'),coord)
        fori(a,'j',v('n'),sample);a.ifzero(v('ok'));a.ret();a.endif()
        a.set(v('size'),1);fori(a,'b',v('r'),lambda b:a.op('add',v('size'),v('size'),v('size')))
        a.op('zeros',v('old'),v('size'));a.op('put',v('old'),v('old'),0,1)
        def label(l):
            a.op('zeros',v('bits'),v('width'));a.set(v('work'),l)
            # h runs in MSB coordinate order after filling the bits backwards.
            def extract(k):
                a.call('div',[v('work'),v('value')],v('work'),2);a.call('sub',[v('h'),v('discard')],v('width'),k);dec(a,v('h'));a.op('add',v('value'),v('value'),1);a.op('put',v('bits'),v('bits'),v('h'),v('value'))
            fori(a,'k',v('width'),extract)
            a.op('zeros',v('allowed'),v('r'));a.set(v('any'),0)
            def block(b):
                a.set(v('ok'),1)
                def check(h):
                    a.op('mul',v('ri'),b,v('width'));a.op('add',v('ri'),v('ri'),h);a.op('get',v('value'),v('requirements'),v('ri'))
                    a.ifzero(v('value'));a.otherwise();a.op('get',v('lb'),v('bits'),h);a.op('eq',v('same'),v('value'),v('lb'));a.ifzero(v('same'));a.set(v('ok'),0);a.endif();a.endif()
                fori(a,'h',v('width'),check);a.op('put',v('allowed'),v('allowed'),b,v('ok'));a.op('add',v('any'),v('any'),v('ok'))
            fori(a,'b',v('r'),block);a.ifzero(v('any'));a.otherwise()
            a.set(v('new'),v('old'))
            def subset(s):
                a.op('get',v('value'),v('old'),s);a.set(v('power'),1)
                def addblock(b):
                    a.call('div',[v('quotient'),v('discard')],s,v('power'));a.call('div',[v('discard'),v('member')],v('quotient'),2)
                    a.ifzero(v('member'));a.otherwise();a.op('get',v('ok'),v('allowed'),b);a.ifzero(v('ok'));a.otherwise();a.call('sub',[v('other'),v('discard')],s,v('power'));a.op('get',v('term'),v('old'),v('other'));a.op('add',v('value'),v('value'),v('term'));a.endif();a.endif();a.op('add',v('power'),v('power'),v('power'))
                fori(a,'b',v('r'),addblock);a.op('put',v('new'),v('new'),s,v('value'))
            fori(a,'s',v('size'),subset,1);a.set(v('old'),v('new'));a.endif()
        fori(a,'label',v('alphabet'),label);a.set(v('last'),v('size'));dec(a,v('last'));a.op('get',v('count'),v('old'),v('last'))
    a.define('countlabels',['p','d','n','width','off','alphabet'],['count'],countlabels)
    def falling(a):
        v=a.var;a.set(v('r'),1)
        def factor(i):
            a.call('sub',[v('f'),v('negative')],v('c'),i);a.ifzero(v('negative'));a.op('mul',v('r'),v('r'),v('f'));a.otherwise();a.set(v('r'),0);a.endif()
        fori(a,'i',v('n'),factor)
    a.define('falling',['c','n'],['r'],falling)
    def odometer(a):
        v=a.var;a.set(v('more'),0);a.set(v('i'),v('n'));a.loop();a.breakzero(v('i'));dec(a,v('i'));a.op('get',v('x'),v('map'),v('i'));a.op('add',v('x'),v('x'),1)
        a.op('lt',v('left'),v('i'),v('r'));a.ifzero(v('left'));a.set(v('radix'),v('rr'));a.otherwise();a.set(v('radix'),v('lr'));a.endif()
        a.op('lt',v('less'),v('x'),v('radix'));a.ifzero(v('less'));a.op('put',v('map'),v('map'),v('i'),0)
        a.otherwise();a.op('put',v('map'),v('map'),v('i'),v('x'));a.set(v('more'),1);a.ret();a.endif();a.repeat()
    a.define('odometer',['map','n','r','lr','rr'],['map','more'],odometer)
    def incident(a):
        # Public 36 x 27 Boolean relation lowered into ordinary bit branches.
        # Missing high bits of canonical naturals are zero, as in scalar ALU.
        g=a.g;v=a.var;table=incidence()
        zero,one=g.new(),g.new()
        def tree(depth,l,r):
            if depth==11:
                g.go(one if l<36 and r<27 and table[27*l+r] else zero);return
            reg=v('l') if depth<6 else v('r')
            z,o=g.new(),g.new();g.pop(reg,z,z,o)
            g.mark(z);tree(depth+1,l,r)
            g.mark(o);tree(depth+1,l+(1<<depth) if depth<6 else l,r if depth<6 else r+(1<<(depth-6)))
        tree(0,0,0)
        g.mark(zero);a.set(v('ok'),0);a.ret()
        g.mark(one);a.set(v('ok'),1)
    a.define('incident',['l','r'],['ok'],incident)
    def inject(a):
        v=a.var;a.set(v('J'),0)
        # Every literal left has 27 neighbors; every literal right has 8.
        # Necessary degree tests preserve all zero terms before enumerating maps.
        a.op('mul',v('seenlength'),v('r'),v('s'));a.op('zeros',v('seen'),v('seenlength'));a.op('zeros',v('ldeg'),v('r'));a.op('zeros',v('rdeg'),v('s'));a.set(v('ok'),1)
        def occurrence(j):
            a.op('get',v('l'),v('p'),j);a.op('get',v('rright'),v('q'),j);a.op('mul',v('edge'),v('l'),v('s'));a.op('add',v('edge'),v('edge'),v('rright'));a.op('get',v('old'),v('seen'),v('edge'))
            a.ifzero(v('old'));a.op('put',v('seen'),v('seen'),v('edge'),1)
            for side,idx,limit in [('ldeg','l',27),('rdeg','rright',8)]:
                a.op('get',v('deg'),v(side),v(idx));a.op('add',v('deg'),v('deg'),1);a.op('put',v(side),v(side),v(idx),v('deg'));a.op('lt',v('excess'),limit,v('deg'));a.ifzero(v('excess'));a.otherwise();a.set(v('ok'),0);a.endif()
            a.endif()
        fori(a,'j',v('n'),occurrence);a.ifzero(v('ok'));a.ret();a.endif()
        a.ifzero(v('source'));a.set(v('lr'),36);a.set(v('rr'),27);a.set(v('lc'),2);a.set(v('rc'),9)
        a.otherwise();a.set(v('lr'),9);a.set(v('rr'),9);a.set(v('lc'),8);a.set(v('rc'),27);a.endif()
        a.op('add',v('mapsize'),v('r'),v('s'));a.op('zeros',v('map'),v('mapsize'));a.set(v('more'),1)
        a.loop();a.breakzero(v('more'));a.op('zeros',v('lo'),v('lr'));a.op('zeros',v('ro'),v('rr'));a.set(v('weight'),1)
        def vertex(k):
            a.op('get',v('ty'),v('map'),k);a.op('lt',v('left'),k,v('r'))
            for isleft in [False,True]:
                if not isleft:a.ifzero(v('left'))
                else:a.otherwise()
                occ,cap=('lo','lc') if isleft else ('ro','rc')
                a.op('get',v('used'),v(occ),v('ty'));a.call('sub',[v('available'),v('negative')],v(cap),v('used'));a.ifzero(v('negative'));a.op('mul',v('weight'),v('weight'),v('available'));a.otherwise();a.set(v('weight'),0);a.endif();a.op('add',v('used'),v('used'),1);a.op('put',v(occ),v(occ),v('ty'),v('used'))
            a.endif()
        fori(a,'k',v('mapsize'),vertex)
        def edge(j):
            a.op('get',v('l'),v('p'),j);a.op('get',v('rright'),v('q'),j);a.op('add',v('rright'),v('rright'),v('r'));a.op('get',v('lt'),v('map'),v('l'));a.op('get',v('rt'),v('map'),v('rright'))
            a.ifzero(v('source'));a.op('incident',v('incident'),v('lt'),v('rt'))
            a.otherwise();a.op('eq',v('incident'),v('lt'),v('rt'));a.endif();a.ifzero(v('incident'));a.set(v('weight'),0);a.endif()
        fori(a,'j',v('n'),edge);a.op('add',v('J'),v('J'),v('weight'));a.call('odometer',[v('map'),v('more')],v('map'),v('mapsize'),v('r'),v('lr'),v('rr'));a.repeat()
    a.define('inject',['p','q','n','r','s','source'],['J'],inject)
    def normalizer(a):
        v=a.var
        for s in ['z0n','z1n']:a.set(v(s),0)
        for s in ['z0d','z1d']:a.set(v(s),1)
        a.ifzero(v('valid'));a.ret();a.endif()
        a.set(v('power'),1);fori(a,'pow',v('n'),lambda k:a.op('mul',v('power'),v('power'),1944))
        a.call('side',[v('p'),v('fullp')],v('d'),v('n'),7,0);a.set(v('morep'),1)
        a.loop();a.breakzero(v('morep'));a.op('blocks',v('r'),v('p'),v('n'));a.op('lt',v('excess'),72,v('r'));a.ifzero(v('excess'))
        a.ifzero(v('fullp'));a.op('countlabels',v('L'),v('p'),v('d'),v('n'),7,0,128);a.otherwise();a.set(v('L'),1);a.endif()
        a.ifzero(v('L'));a.otherwise();a.call('side',[v('q'),v('fullq')],v('d'),v('n'),8,7);a.set(v('moreq'),1)
        a.loop();a.breakzero(v('moreq'));a.op('blocks',v('s'),v('q'),v('n'));a.op('lt',v('excess'),243,v('s'));a.ifzero(v('excess'))
        a.ifzero(v('fullq'));a.op('countlabels',v('R'),v('q'),v('d'),v('n'),8,7,256);a.otherwise();a.set(v('R'),1);a.endif()
        a.ifzero(v('R'));a.otherwise();a.op('falling',v('lf'),128,v('r'));a.op('falling',v('rf'),256,v('s'));a.op('mul',v('den'),v('power'),v('lf'));a.op('mul',v('den'),v('den'),v('rf'));a.op('mul',v('LR'),v('L'),v('R'))
        for source in (0,1):
            a.op('inject',v('J'),v('p'),v('q'),v('n'),v('r'),v('s'),source);a.op('mul',v('num'),v('J'),v('LR'))
            a.ifzero(v('num'));a.otherwise();a.call('qadd',[v(f'z{source}n'),v(f'z{source}d')],v(f'z{source}n'),v(f'z{source}d'),v('num'),v('den'));a.endif()
        a.endif();a.endif()
        a.ifzero(v('fullq'));a.call('nextpart',[v('q'),v('moreq')],v('q'),v('n'));a.otherwise();a.set(v('moreq'),0);a.endif();a.repeat()
        a.endif();a.endif()
        a.ifzero(v('fullp'));a.call('nextpart',[v('p'),v('morep')],v('p'),v('n'));a.otherwise();a.set(v('morep'),0);a.endif();a.repeat()
    a.define('normalizer',['d','n','valid'],['z0n','z0d','z1n','z1d'],normalizer)
    def mixture(a):
        v=a.var;a.call('sub',[v('complement'),v('discard')],v('b'),v('a'));a.op('mul',v('an'),v('n0'),v('complement'));a.op('mul',v('ad'),v('d0'),v('b'));a.op('mul',v('bn'),v('n1'),v('a'));a.op('mul',v('bd'),v('d1'),v('b'));a.call('qadd',[v('n'),v('d')],v('an'),v('ad'),v('bn'),v('bd'))
    a.define('mixture',['n0','d0','n1','d1','a','b'],['n','d'],mixture)
    def addword(a):
        v=a.var;a.op('encq',v('word'),v('n'),v('d'));a.op('encv',v('word'),v('word'));a.op('concat',v('out'),v('out'),v('word'))
    a.define('addword',['out','n','d'],['out'],addword)


def main(a):
    v=a.var;g=a.g
    a.call('readu',[v('n'),'input'],'input');bit(a,'input',v('sign'));a.call('readu',[v('a'),'input'],'input');a.call('readu',[v('b'),'input'],'input')
    a.call('readcache',[v('cache'),'input'],'input',v('n'));a.call('readu',[v('t'),'input'],'input')
    a.call('normalizer',[v('z0n'),v('z0d'),v('z1n'),v('z1d')],v('cache'),v('n'),1)
    a.call('mixture',[v('zn'),v('zd')],v('z0n'),v('z0d'),v('z1n'),v('z1d'),v('a'),v('b'))
    a.set(v('out'),0);a.set(v('ratios'),0)
    for n,d in [('z0n','z0d'),('z1n','z1d'),('zn','zd')]:a.op('addword',v('out'),v('out'),v(n),v(d))
    def cylinder(k):
        a.call('readv',[v('item'),'input'],'input');bit(a,v('item'),v('conflict'));a.ifzero(v('conflict'))
        a.call('readcache',[v('e'),v('item')],v('item'),v('n'));a.call('conjoin',[v('data'),v('valid')],v('cache'),v('e'),v('n'))
        a.otherwise();a.set(v('valid'),0);a.set(v('data'),v('cache'));a.endif()
        a.call('normalizer',[v('m0n'),v('m0d'),v('m1n'),v('m1d')],v('data'),v('n'),v('valid'))
        a.call('mixture',[v('mn'),v('md')],v('m0n'),v('m0d'),v('m1n'),v('m1d'),v('a'),v('b'))
        for n,d in [('m0n','m0d'),('m1n','m1d'),('mn','md')]:a.op('addword',v('out'),v('out'),v(n),v(d))
        a.ifzero(v('zn'));a.otherwise();a.op('mul',v('rn'),v('mn'),v('zd'));a.op('mul',v('rd'),v('md'),v('zn'));a.op('addword',v('ratios'),v('ratios'),v('rn'),v('rd'));a.endif()
    fori(a,'k',v('t'),cylinder)
    a.op('mul',v('count'),v('t'),3);a.op('add',v('count'),v('count'),3);a.op('encu',v('header'),v('count'));a.op('concat',v('out'),v('header'),v('out'))
    a.ifzero(v('zn'));a.raw(v('tail'),[0]);a.otherwise();a.op('encu',v('header'),v('t'));a.op('concat',v('tail'),v('header'),v('ratios'));g.push(v('tail'),1);a.endif()
    a.op('concat','input',v('out'),v('tail'))


def build():
    a=Assembler();primitives(a);define_evaluator(a);return a.compile(main)

if __name__=='__main__':
    import json,sys
    p=build();json.dump(p,sys.stdout,separators=(',',':'));print()
