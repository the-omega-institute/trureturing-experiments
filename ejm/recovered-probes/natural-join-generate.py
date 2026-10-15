import argparse,itertools,json,pathlib
parser=argparse.ArgumentParser()
parser.add_argument('--input-directory',required=True,type=pathlib.Path)
parser.add_argument('--worktree',required=True,type=pathlib.Path)
args=parser.parse_args()
p=args.input_directory
inp=json.loads((p/'support-orbit-join/input.json').read_text())
perms=json.loads((p/'local-perms.json').read_text())['perms']
support=[set(),set()]
for row in inp['source_rows']:
    a,b=row['assignment'],row['block']
    support[b].update(tuple(v[a[g[j]]] for j in range(len(a)))
        for g in perms[b] for v in itertools.permutations(range(4)))
assert list(map(len,support))==[132,168]
def rows(name,b):
    vals=',\n  '.join('['+', '.join(map(str,a))+']' for a in sorted(support[b]))
    return f'def {name} : List (List (Fin 4)) :=\n  [\n  {vals}\n  ]\n'
prefix='''/- GID: D5/S3/Quantum/Entanglement/EJM/EJMOriginalSupportConflict
   generality: I
   mirror-B: D5/B/S3/Quantum/Entanglement/EJM/EJMOriginalSupportConflict
   mirror-E: none(waiver:original-support-natural-join-refutation)
   anchors: [D5/S0/Certificates/FiniteNaturalJoin,
     D5/S3/Quantum/Entanglement/EJMInflationScore]
   utility: kind=certified-instance; basis=refutes=gid:D5/S3/Quantum/Entanglement/EJM/EJMOriginalSupportConflict.claim; result=D5/S3/Quantum/Entanglement/EJM/EJMOriginalSupportConflict.result; claim=D5/S3/Quantum/Entanglement/EJM/EJMOriginalSupportConflict.claim
   digest: No original twenty-coordinate event simultaneously occupies the full outcome orbits at singleton8 and triangles2 and3. -/

import D5.S0.Certificates.FiniteNaturalJoin
import D5.S3.Quantum.Entanglement.EJMInflationScore

set_option autoImplicit false
set_option relaxedAutoImplicit false
set_option maxRecDepth 100000
set_option linter.style.longLine false

namespace D5.S3.Quantum.Entanglement.EJM.EJMOriginalSupportConflict

open D5.S0.Certificates.FiniteNaturalJoin
open D5.S3.Quantum.Entanglement.ElegantJointMeasurementInflation (Node)
open D5.S3.Quantum.Entanglement.EJMInflationScore (lhsTriangle lhsSingleton)

/- These literal supports are the unions of the twelve local-coordinate and
all twenty-four simultaneous outcome relabelings of original triangle rows
104,259,332,464 and singleton rows220,252. They include all outcomes, impose
no fixed coordinate, and use the actual ordered original marginal maps. -/

'''
suffix='''
def triangleRelation (i : Fin 16) : List (Row Node (Fin 4)) :=
  triangleSupport.map fun row => (List.ofFn (lhsTriangle i)).zip row

def singletonRelation (i : Fin 16) : List (Row Node (Fin 4)) :=
  singletonSupport.map fun row => (List.ofFn (lhsSingleton i)).zip row

/-- A full event in all three literal original support relations. -/
def claim : Prop :=
  ∃ e : Node → Fin 4, Covers (singletonRelation 8) e ∧
    Covers (triangleRelation 2) e ∧ Covers (triangleRelation 3) e

set_option maxHeartbeats 4000000 in
private theorem checked_empty :
    join (join (singletonRelation 8) (triangleRelation 2))
      (triangleRelation 3) = [] := by decide

/-- The actual empty join refutes simultaneous membership, without assuming
an outcome normalization, a global replay receipt or a score inequality. -/
theorem result : ¬ claim := by
  rintro ⟨e, hs, ht2, ht3⟩
  have h := join_covers _ _ e (join_covers _ _ e hs ht2) ht3
  rw [checked_empty] at h
  simpa [Covers] using h

#print axioms result

end D5.S3.Quantum.Entanglement.EJM.EJMOriginalSupportConflict
'''
path=args.worktree/'D5/S3/Quantum/Entanglement/EJM/EJMOriginalSupportConflict.lean'
path.write_text(prefix+rows('triangleSupport',0)+'\n'+rows('singletonSupport',1)+suffix)
print(path.name, len(path.read_bytes()), 'bytes; full support sizes',list(map(len,support)))
