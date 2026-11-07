"""Static lowering of the published register decoder to direct field graphs.

Neither Decoder.tick/result nor constructor_image.Image selects execution.
Linking is public preparation. Engine installation remains a separate boundary.
"""
from field_graph import Graph, Engine, loader, driver, growth
from paid_code import decoder_program


def page_selector(g):
    # Physical field selection is an actual complete finite binary decision tree.
    # The leaf's index is a literal in code, never a computed runtime integer.
    def tree(depth,index):
        if depth==15:
            g.ref_names.add('page.node');g.bit_names.add('page.value')
            g.emit('page','page.value','page.node',index)
            g.ret('page.return');return
        one=g.fresh('page.right')
        g.test(f'page.offset.{14-depth}',one)
        tree(depth+1,2*index);g.label(one);tree(depth+1,2*index+1)
    g.label('page.select');tree(0,0)


def load(g,dst,src,width,height):
    g.wcopy('load.address',src)
    for bit in range(width):
        g.cp('page.node','image.root')
        for level in range(height):
            right,done=g.fresh('dir.right'),g.fresh('dir.done')
            g.test(f'load.address.{12+height-1-level}',right)
            g.emit('edge','page.node','page.node',0);g.j(done)
            g.label(right);g.emit('edge','page.node','page.node',1);g.label(done)
        for j in range(12):g.bitcopy(f'page.offset.{j+3}',f'load.address.{j}')
        for j in range(3):g.put(f'page.offset.{j}',((bit%8)>>j)&1)
        resume=g.fresh('page.resume')
        g.bind('page.return',resume);g.j('page.select');g.label(resume)
        g.bitcopy(f'{dst}.{bit}','page.value')
        if bit%8==7 and bit+1<width:
            # Carry addition remains ordinary graph bits and tests.
            g.constant('load.one',1);g.alu('add','load.address','load.address','load.one')
    for bit in range(width,64):g.put(f'{dst}.{bit}',0)


def lower(g,register_graph,height,boot_done,finish):
    g.label('fetch');g.j('decoder.'+str(register_graph.labels['decode']))
    for i,(op,dst,a,b,imm,nxt,alt) in enumerate(register_graph.rows):
        g.label('decoder.'+str(i))
        if op=='const':g.constant(dst,imm)
        elif op in ('copy','shl','shr','add','sub','eq','lt'):g.alu(op,dst,a,b or None,imm)
        elif op.startswith('load'):load(g,dst,a,int(op[4:]),height)
        elif op=='branch':g.test(f'{a}.0','decoder.'+str(alt))
        elif op=='jump':g.j('decoder.'+str(alt))
        elif op=='wait':g.j(boot_done)
        elif op=='finish':g.j(finish)
        elif op=='fault':g.j('fault')
        else:raise AssertionError(op)
        g.j('decoder.'+str(nxt))


def decoder_graph(height,graph=None):
    graph=graph or decoder_program()
    g=Graph();g.label('start');g.j('decoder.'+str(graph.labels['boot']))
    g.label('decoder.boot.done');g.emit('wait')
    g.j('fetch')
    lower(g,graph,height,'decoder.boot.done','decoder.finished')
    g.label('decoder.finished')
    for bit in range(3):g.emit('emit',f'out_op.{bit}')
    for bit in range(64):g.emit('emit',f'out_dest.{bit}')
    g.emit('end');g.emit('ack');g.emit('return');g.emit('wait');g.j('fetch')
    page_selector(g)
    g.label('fault');g.emit('fault','invalid code image/address')
    return g.seal()


def image_b_graph(height):
    graph=decoder_program();g=Graph()
    g.label('start');g.constant('bpc',0);g.j('decoder.'+str(graph.labels['boot']))
    lower(g,graph,height,'load','image.fetched')
    g.label('image.fetched');g.ret('fetch.return')
    loader(g);driver(g,image=True);growth(g);page_selector(g)
    g.label('fault');g.emit('fault','invalid image or LAB phase')
    g.label('exhausted');g.emit('fault','finite supplied pool exhausted')
    return g.seal()


def request(e,pc):
    # Test harness supplies 64 actual incoming bits. This adapter has no paid
    # stream/validation correspondence yet, and is never a LAB instruction source.
    assert e.waiting and 0<=pc<1<<64
    for i in range(64):e.bw(f'input_pc.{i}',(pc>>i)&1)
    e.resume();e.run()
    return e.observe_word('out_op'),e.observe_word('out_dest')


def execution_graph(height,nregisters=361):
    """Composition for conversion, initialized B prefixes and actual output scan.

    Its source/session installer and original-port service are still absent.
    This entry is author preparation, not C release or integrated acceptance.
    """
    from field_io import conversions, output_decoder, print_buffer
    graph=decoder_program();g=Graph()
    g.label('start');g.constant('bpc',0);g.j('decoder.'+str(graph.labels['boot']))
    lower(g,graph,height,'convert','image.fetched')
    g.label('image.fetched');g.ret('fetch.return')
    conversions(g,nregisters)
    g.label('converted.done');g.get('input','converted',1);g.j('load')
    loader(g);driver(g,image=True,halt_continuation='decode.output')
    output_decoder(g)
    g.label('decoded.done');print_buffer(g,'decoded')
    g.label('printed');g.emit('wait');g.j('printed')
    growth(g);page_selector(g)
    g.label('fault');g.emit('fault','invalid image, LAB phase, or returned output')
    g.label('exhausted');g.emit('fault','finite supplied pool exhausted')
    return g.seal()
