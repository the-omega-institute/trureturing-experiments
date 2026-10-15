"""Prepare a bounded actual four-child assembly probe, never a global certificate.

Every context carries its literal depth/code/index. A suffix is decided at
that full actual prefix; every internal node earns its entire split proposal.
Native data only choose the proof shape. Lean remains responsible for every
node, every numerical suffix and the final actual-context assembly.
"""
import argparse,hashlib,json,pathlib,struct
ap=argparse.ArgumentParser();ap.add_argument('--cap',type=int,default=64)
a=ap.parse_args();assert 1<=a.cap<=256
r=pathlib.Path(__file__).resolve().parent.parent
cp=pathlib.Path('/Users/macstudio/omega/.formalization-goal-01a12631/checkpoints/ejm-mdd-global/existing-cover-original-suffix-alignment-20261014/original-source-identity-stream')
raw=(cp/'dag-87458265.bin').read_bytes()
assert hashlib.sha256(raw).hexdigest()=='bf188b5ef8d359fc12872297ba3c3db9b633debdec7db6eddf924265c880bf59'
magic,size=struct.unpack_from('<II',raw);assert (magic,size)==(0x454a4d32,1558777)
def row(j): return struct.unpack_from('<6I',raw,8+24*j)
def count(j):
 kind,_,*children=row(j)
 return 1+sum(count(c) for c in children) if kind==2 else 1
# The independently bound original pilot has a real numeric subtree.
assert count(55)==201
ns='D5.S3.Quantum.Entanglement.EJM'
name='TransientEJMOriginalAdaptivePilot'
lines=[f'import {ns}.EJMOriginalReplayAssembly','set_option autoImplicit false',
 'set_option relaxedAutoImplicit false','set_option maxHeartbeats 100000000',
 'set_option maxRecDepth 65536',f'namespace {ns}.{name}',
 f'open {ns}.EJMOriginalReplayAssembly',
 'open D5.S0.Certificates.FiniteCachedSymmetryDAGReplay (checkSplit checked_split_node)']
node_proofs={};contexts=[]
def fin(n):return f'⟨{n}, by decide⟩'
def emit(j,depth,code):
 assert 0<=depth<=20 and 0<=code<4**depth
 assert count(j)<=201
 label=f'accepted_{depth}_{code}'
 kind,coord,*children=row(j)
 contexts.append(dict(index=j,depth=depth,code=code,contexts=count(j)))
 if count(j)<=a.cap:
  lines.extend([f'private theorem {label} : Accepted {depth} {code} {fin(j)} := by',
   '  decide +kernel',''])
  return label
 assert kind==2 and coord==depth and depth<20 and all(c<j for c in children)
 child_labels=[emit(c,depth+1,4*code+x) for x,c in enumerate(children)]
 key=(depth,j)
 child_vector='!['+','.join(fin(c) for c in children)+']'
 if key not in node_proofs:
  node_label=f'node_{depth}_{j}';node_proofs[key]=node_label
  lines.extend([f'private theorem {node_label} :',
   f'    EJMOriginalGlobalDAGData.nodes {21-depth} {fin(j)} =',
   f'      .split (coordinate {depth}) ({child_vector}) := by',
   '  apply checked_split_node', '  decide +kernel',''])
 lines.extend([f'private theorem {label} : Accepted {depth} {code} {fin(j)} := by',
  f'  apply assemble {depth} {code} {fin(j)} (by decide) ({child_vector}) {node_proofs[key]}',
  '  intro a','  fin_cases a'])
 lines.extend('  · exact '+c for c in child_labels);lines.append('')
 return label
label=emit(55,15,0)
lines.extend([f'#print axioms {label}',f'end {ns}.{name}',''])
p=r/'D5/S3/Quantum/Entanglement/EJM'/f'{name}.lean';p.write_text('\n'.join(lines))
record=dict(status='prepared bounded actual adaptive assembly probe, uncompiled',
 kernel_accepted=False,full_global_accepted=False,root=55,depth=15,code=0,
 suffix_cap=a.cap,full_contexts=201,contexts=contexts,node_equalities=len(node_proofs),
 source=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest())
(r/'.ejm-original-global-transfer1/adaptive-pilot.json').write_text(json.dumps(record,indent=2)+'\n')
print(json.dumps({k:v for k,v in record.items() if k!='contexts'}))
