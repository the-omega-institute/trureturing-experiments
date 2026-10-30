"""Reusable capacity-coloring checker extracted from the completed probe.

Retention modification: wall-time/process wrappers were removed; a node limit
raises NodeLimit (unknown). The branching, cache and isolate filling are retained.
Only the pre-existing small controls are replayed by verify_retained.py.
"""
from packing_core import NodeLimit, independent

def colorpack(a,ks,max_nodes=1000000):
 """A full packing iff every vertex can be allocated to bounded independent color classes plus unrestricted slack color."""
 n=len(a);ks=tuple(ks);slack=n-sum(ks);assert slack>=0
 ell=len(ks);caps=list(ks)+[slack];masks=[0]*(ell+1);nodes=0
 # Isolates can fill every remaining capacity after core assignment, including slack.
 core=sum(1<<v for v in range(n) if a[v]);iso=((1<<n)-1)^core
 failed=set()
 def rec(rem):
  nonlocal nodes
  nodes+=1
  if nodes>max_nodes:raise NodeLimit()
  if not rem:return True
  key=(rem,tuple(sorted((ks[c],masks[c]) for c in range(ell))),masks[ell])
  if key in failed:return False
  # Minimum remaining color-domain, then largest degree in remaining core.
  best=None;chosen=None;opts=None
  rr=rem
  while rr:
   b=rr&-rr;rr-=b;v=b.bit_length()-1
   allowed=[c for c in range(ell) if caps[c] and not a[v]&masks[c]]
   if caps[ell]:allowed.append(ell)
   if not allowed:failed.add(key);return False
   rank=(len(allowed),-(a[v]&rem).bit_count(),-a[v].bit_count(),v)
   if best is None or rank<best:best=rank;chosen=b;opts=allowed
  v=chosen.bit_length()-1
  # Prefer tightly constrained filled classes; leave unrestricted slack last.
  opts.sort(key=lambda c:(c==ell,-(a[v]&rem).bit_count() if c<ell else 0,-caps[c],c))
  empty_seen=set()
  for c in opts:
   if c<ell and not masks[c]:
    if ks[c] in empty_seen:continue
    empty_seen.add(ks[c])
   caps[c]-=1;masks[c]|=chosen
   if rec(rem^chosen):return True
   masks[c]^=chosen;caps[c]+=1
  failed.add(key);return False
 if not rec(core):return None,{'nodes':nodes,'failed_states':len(failed)}
 iso_vs=[v for v in range(n) if (iso>>v)&1];j=0
 for c in range(ell+1):
  for _ in range(caps[c]):masks[c]|=1<<iso_vs[j];j+=1
 assert j==len(iso_vs)
 # Input order retained here; ordinary witness validation independent of search.
 used=0
 for k,m in zip(ks,masks[:ell]):
  assert m.bit_count()==k and not (used&m) and independent(a,m);used|=m
 return tuple(masks[:ell]),{'nodes':nodes,'failed_states':len(failed)}
