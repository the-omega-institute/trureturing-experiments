import json,subprocess
polys={}
for p in ['polynomials3.txt','polynomials4.txt','polynomials5.txt']:
 for l in open(p):
  b,c=l.strip().split(';');b=tuple(map(int,b.rstrip(',').split(',')));c=tuple(map(int,c.rstrip(',').split(',')));polys.setdefault(c,b)
for r in json.load(open('targeted-checked-results.json'))['results']:polys.setdefault(tuple(r['last_coefficients']),tuple(r['bound']))
with open('product-input.txt','w') as f:
 f.write(str(len(polys))+'\n')
 for c,b in polys.items():f.write(str(len(c))+' '+' '.join(map(str,c))+' '+str(len(b))+' '+' '.join(map(str,b))+'\n')
print('unique_polynomials',len(polys),flush=True)
