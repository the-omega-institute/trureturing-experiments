#!/usr/bin/env python3
"""Exact enlarged-patch cover; an external computational proof, not Lean admission.

Reuses historical dyadic core and Fraction kernel. The core differs only in its
uniqueness-box radius: 2^-12 rather than 2^-16. Canonical shift fixation is
transported from independently revalidated centre branches by equivariance,
continuity, local uniqueness, and connectedness, not by an unchecked label.
"""
from fractions import Fraction as F
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import argparse,itertools,json,subprocess,hashlib
import check_real_x_global_cover as original

SCALE=1<<40
DELTA=F(1,1<<24)
RADIUS=F(1,1<<12)
I,C=original.I,original.C
require=original.require

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def make_bounds(out):
 H=original.old.parameter_matrix([F(-1,5),F(0)],DELTA)
 witnesses=original.old.squared_minor_separation(H)
 require(len(witnesses)==30,'missing strict-X witnesses')
 lines=['40']
 for row in H:
  for z in row:
   v=[]
   for x in (z.re,z.im):v.extend([x.lo.numerator*SCALE//x.lo.denominator,-((-x.hi.numerator*SCALE)//x.hi.denominator)])
   lines.append(' '.join(map(str,v)))
 out.write_text('\n'.join(lines)+'\n')
 return witnesses

def graph_audit(data,enclosures):
 lines=enclosures.read_text().splitlines();require(len(lines)==61 and lines[0]=='60','malformed enclosures')
 rays={}
 for line in lines[1:]:
  a=list(map(int,line.split()));require(len(a)==12,'malformed root');idx,mask=a[:2]
  require(0<=idx<60 and idx not in rays and 0<=mask<32,'bad root label')
  X=[I(F(a[2+2*j],SCALE),F(a[3+2*j],SCALE))for j in range(5)]
  rays[idx]=original.base.phase_map(X,[0]+[2*((mask>>j)&1)for j in range(5)])
 edges={tuple(e)for e in data['allowed_edges']};can=set(data['canonical_vertices']);colors=data['noncanonical_colors']
 require(len(edges)==114 and len(can)==6 and len(colors)==60,'bad supergraph')
 require(all(0<=i<j<60 for i,j in edges)and all(0<=i<60 for i in can),'bad endpoints')
 require(all(colors[i]in(0,1)for i in range(60)if i not in can),'bad colors')
 require(all((i in can and j in can)or(i not in can and j not in can and colors[i]!=colors[j])for i,j in edges),'noncanonical triangle may be allowed')
 count=0;minimum=None
 for i,j in itertools.combinations(range(60),2):
  ip=sum((a.conj()*b for a,b in zip(rays[i],rays[j])),C.point())/6;n=ip.normsq()
  require(n.hi<1,f'projective distinctness failed {i},{j}')
  if(i,j)not in edges:
   require(n.lo>F(1,500000),f'nonedge margin failed {i},{j}')
   minimum=n.lo if minimum is None else min(minimum,n.lo);count+=1
 require(count==1656,'missing nonedge')
 return {'distinct_rays':60,'pair_checks':1770,'certified_nonedges':1656,'nonedge_squared_modulus_strict_lower':'1/500000','canonical_vertices':sorted(can),'canonical_shift_fixed_by':'continuous simple-root branch propagation from independently revalidated centre','all_orthogonal_triples_canonical':True}

def run(input_path,output,jobs,cap):
 require(__debug__,'assertions must be enabled');require(jobs>0 and cap>0,'bad resource limit')
 output.mkdir(parents=True,exist_ok=False)
 data=json.loads(input_path.read_text());here=Path(__file__).parent
 old_core=here/'check_real_x_global_cover.cpp';new_core=here/'cover_larger_patch.cpp'
 old_text=old_core.read_text();require(old_text.count('1LL<<24')==2,'historical core changed')
 require(new_core.read_text()==old_text.replace('1LL<<24','1LL<<28'),'unexpected core changes')
 pre=output/'centre_preflight';pre.mkdir()
 centres,old_bounds,_=original.prepare(data,pre)
 old_binary=pre/'old_core';new_binary=output/'cover_core'
 for source,binary in [(old_core,old_binary),(new_core,new_binary)]:
  subprocess.run(['g++','-O3','-std=c++17','-Wall','-Wextra','-Werror',str(source),'-o',str(binary)],check=True)
 # Node cap 1 intentionally leaves phase space uncovered. All 60 root inclusion
 # validations and enclosure emissions run before that traversal begins.
 p=subprocess.run([str(old_binary),str(centres),'0','1',str(pre/'chart0.json'),str(old_bounds),str(pre/'root_enclosures.txt')],capture_output=True,text=True)
 require(p.returncode==2,'centre preflight did not reject incomplete cover')
 q=json.loads((pre/'chart0.json').read_text());require(q['status']=='INCOMPLETE'and q['pending']>0,'bad incomplete-cover control')
 centre_graph=original.graph_audit(data,pre/'root_enclosures.txt',centres)
 require(centre_graph['canonical_shift_fixed'],'centre symmetry not certified')
 require(RADIUS>original.RADIUS,'old uniqueness boxes not contained in new boxes')
 # Exactly the same centre file supplies both cores, so every old uniqueness
 # box is inside its newly validated box. At alpha=-1/5 the unique roots agree.
 bounds=output/'matrix_bounds.txt';witnesses=make_bounds(bounds)
 def chart(mask):
  report=output/f'chart_{mask:02d}.json'
  args=[str(new_binary),str(centres),str(mask),str(cap),str(report),str(bounds)]
  if mask==0:args.append(str(output/'root_enclosures.txt'))
  p=subprocess.run(args,capture_output=True,text=True)
  require(p.returncode==0,f'chart {mask} failed: {p.stderr}{p.stdout}')
  r=json.loads(report.read_text());require(r['chart']==mask and r['status']=='COVERED'and r['pending']==r['unresolved']==0,'incomplete chart')
  require(r['nodes']==2*(r['excluded']+r['known_root_leaves'])+r['contracted']-1,'invalid complete search-tree counts')
  return r
 with ThreadPoolExecutor(max_workers=jobs)as ex:charts=list(ex.map(chart,range(32)))
 require([r['chart']for r in charts]==list(range(32)),'missing chart')
 graph=graph_audit(data,output/'root_enclosures.txt')
 report={'status':'PASS','parameter_center':['-1/5','0'],'parameter_radius':str(DELTA),'uniqueness_radius':str(RADIUS),'all_32_charts_covered':True,'nodes':sum(r['nodes']for r in charts),'pending':0,'unresolved':0,'global_root_count':60,'orthogonal_triple_collision':'3','all_actual_completion_affinities':'2','graph':graph,'centre_identification':{'same_centres_file':True,'old_radius':str(original.RADIUS),'new_radius':str(RADIUS),'old_box_contained_in_new':True,'fresh_centre_shift_fixation':centre_graph},'charts':charts,'fourier_separation_witnesses':witnesses,'lean_kernel_verified':False,'scope':'Only the explicitly constructed closed two-real-parameter rectangle; no full strict-X claim','hashes':{k:sha(p)for k,p in [('driver',Path(__file__)),('historical_core',old_core),('enlarged_core',new_core),('historical_driver',Path(original.__file__)),('parameter_generator',Path(original.old.__file__)),('fraction_kernel',Path(original.base.__file__)),('input',input_path),('centres',centres),('matrix_bounds',bounds),('root_enclosures',output/'root_enclosures.txt')]}}
 (output/'verification.json').write_text(json.dumps(report,indent=2)+'\n')
 return report
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('input',type=Path);p.add_argument('--output',type=Path,required=True);p.add_argument('--jobs',type=int,default=4);p.add_argument('--max-nodes-per-chart',type=int,default=1000000);a=p.parse_args();r=run(a.input,a.output,a.jobs,a.max_nodes_per_chart);print(json.dumps({k:v for k,v in r.items()if k not in('charts','fourier_separation_witnesses')},indent=2))
