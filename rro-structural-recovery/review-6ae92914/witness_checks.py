"""Actual-source evaluator checks for positive closure and second mixed formula."""
from fractions import Fraction as F
from pathlib import Path
import json,sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from exact_evaluator import normalizer
from exact_checks import complete,falling,M
# The second mixed prediction: compute the two full joint evidence masses using
# the unrestricted evaluator, rather than only dividing literal graph counts.
five=tuple(complete(a,t) for a,t in ((0,0),(0,1),(1,0),(1,1),(2,0)))
six=five+(complete(2,1),)
z5=normalizer(five); z6=normalizer(six)
assert z5==tuple(F(x,M**5*falling(128,3)*falling(256,2)) for x in (863136,2122848))
assert z6==tuple(F(x,M**6*falling(128,3)*falling(256,2)) for x in (653184,2122848))
predictions={}
for t in (F(0),F(1,7),F(1,3),F(1,2),F(1)):
 ratio=((1-t)*z6[0]+t*z6[1])/((1-t)*z5[0]+t*z5[1])
 assert ratio==F(28+63*t,M*(37+54*t))
 predictions[str(t)]=str(ratio)
# The seven-edge positive witness has a new forced eighth incidence under each
# source. Its full joint likelihood ratio must be exactly 1/M, separately.
seven=tuple(complete(a,t) for a,t in ((0,0),(0,1),(1,0),(1,1),(2,0),(2,1),(3,0)))
eight=seven+(complete(3,1),)
z7=normalizer(seven); z8=normalizer(eight)
assert all(x>0 for x in z7)
assert tuple(z8[i]/z7[i] for i in (0,1))==(F(1,M),)*2
result={'status':'passed','second_mixed_prediction':predictions,'Z5':list(map(str,z5)),'Z6':list(map(str,z6)),'Z7':list(map(str,z7)),'Z8':list(map(str,z8)),'forcing_ratio_both_sources':str(F(1,M)),'scope':'Bounded exact witnesses, not an all-cache proof.'}
Path(__file__).with_name('witness-results.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,sort_keys=True))
