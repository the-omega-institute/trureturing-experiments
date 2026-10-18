#!/usr/bin/env python3
"""Stream-check the full existing f,g table with standard-library exact arithmetic."""
import csv,gzip,json,math,resource,sys,time
from fractions import Fraction as F
from pathlib import Path
ROOT=Path(__file__).resolve().parent
if not __debug__:
 raise SystemExit('Validation requires assertions: rerun Python without -O or -OO.')
if hasattr(sys,'set_int_max_str_digits'):sys.set_int_max_str_digits(0)
start=time.monotonic();expected_n=1;C=1;exceptions=[];first_squareful=None;maximum=(0,0);first=[]
with gzip.open(ROOT/'bounds_1_10000.csv.gz','rt') as handle:
 reader=csv.DictReader(handle)
 assert reader.fieldnames==['N','f_num','f_den','g_num','g_den','ratio_num','ratio_den','c'],reader.fieldnames
 for m in range(1,5001):
  # Independent central-binomial recurrence, no source generator imported.
  C,remainder=divmod(C*2*(2*m-1),m)
  assert remainder==0,(m,remainder)
  v=(C&-C).bit_length()-1
  b=C>>v;G=(4**m-C)>>v;c=math.gcd(m,G)
  assert b%2==G%2==1 and math.gcd(b,G)==1
  f=F(m*C,2**(2*m-1));g=F(4**m-C,2**(2*m-1));q=f/g
  assert g.numerator==G and q.denominator==G//c and q.numerator==m*b//c
  if m<=100 or m%97==0 or m==5000:assert C==math.comb(2*m,m)
  if c>maximum[1]:maximum=(m,c)
  if c>1:
   exceptions.extend([(2*m-1,m,c),(2*m,m,c)])
   if first_squareful is None and any(c%(p*p)==0 for p in range(3,math.isqrt(c)+1,2)):first_squareful=(m,c)
  for n in (2*m-1,2*m):
   row=next(reader);assert int(row['N'])==n==expected_n
   actual=tuple(int(row[k]) for k in ['f_num','f_den','g_num','g_den','ratio_num','ratio_den','c'])
   expected=(f.numerator,f.denominator,g.numerator,g.denominator,q.numerator,q.denominator,c)
   assert actual==expected,(n,actual,expected)
   if n<=36:first.append(q.numerator)
   expected_n+=1
 assert next(reader,None) is None
assert expected_n==10001 and exceptions[0]==(19,10,5) and len(exceptions)==396
assert first_squareful==(2466,9) and maximum==(871,871)
# Read every displayed term from the fresh public internal page export.
text=(ROOT.parent/'sources/entry-internal.txt').read_text()
terms=[]
for line in text.splitlines():
 if line.startswith(('%S ','%T ','%U ')):terms.extend(int(x) for x in line[3:].split(',') if x.strip())
assert len(terms)==36 and terms==first,(len(terms),terms)
with (ROOT/'cancellation-factors.csv').open('w',newline='') as f:
 w=csv.writer(f);w.writerow(['n','m','c']);w.writerows(exceptions)
rss=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
summary={'checked_n':[1,10000],'table_rows':10000,'entry_terms_checked':len(terms),'exception_count':len(exceptions),'pair_count':len(exceptions)//2,'first_exception':exceptions[0],'first_squareful_m_c':first_squareful,'max_m_c':maximum,'elapsed_seconds':round(time.monotonic()-start,3),'peak_rss_bytes':rss if sys.platform=='darwin' else rss*1024,'single_process':True}
(ROOT/'verification.json').write_text(json.dumps(summary,indent=2))
print(json.dumps(summary,indent=2))
