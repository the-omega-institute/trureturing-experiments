"""Exploratory floating point diagnostic, not certificate evidence."""
from pathlib import Path
import json,numpy as np
import urllib.request,zipfile,io
url='https://static.renyi.hu/ai-shared/daniel/witsenhausen/double_cap.zip'
with urllib.request.urlopen(url,timeout=90) as response:
 archive=zipfile.ZipFile(io.BytesIO(response.read()))
rows=[]
for N in [14,20,33]:
 d=f'sphere_erdos_best_{N}_data'
 j=json.loads(archive.read(d+'/dual.json'))
 pairs=list(j['distance_dictionary'].values())
 vertices=np.load(io.BytesIO(archive.read(d+'/vertices_p2zr5.npy')))
 v=vertices[:,:,0]+np.sqrt(5)*vertices[:,:,1]
 v=v/np.linalg.norm(v,axis=1)[:,None]
 t=np.array([np.dot(v[a[0][0]],v[a[0][1]]) for a in pairs])
 w=np.array(j['w_d']);assert len(t)==len(w)
 # P_even(-t)=P_even(t), so no arbitrary antipodal signs enter this diagnostic.
 # One recurrence sweep avoids quadratic repeated special-function evaluation.
 grid=np.concatenate(([0.0],t));previous=np.ones(len(grid));current=grid.copy()
 beta=np.zeros(65537);beta[0]=j['w_0']+np.sum(w)
 for k in range(2,65537):
  nxt=((2*k-1)*grid*current-(k-1)*previous)/k
  previous,current=current,nxt
  if k%2==0:beta[k]=j['w_0']*current[0]+current[1:]@w
 blocks=[]
 for start in [32,128,512,2048,8192,32768]:
  defects=beta[start:2*start+1:2]
  blocks.append({'even_degree_block':[start,2*start],'max_abs_beta_minus_w1':float(np.max(np.abs(defects))),'min_beta':float(j['w_1']+np.min(defects))})
 row={'pairs':N,'source_q1':j['q_1'],'diagonal_limit_w1':j['w_1'],'distinct_off_diagonal_classes':len(t),'blocks':blocks};rows.append(row)
 print(json.dumps(row))
print('DIAGNOSTIC_ONLY: floating point; no all-degree bound or Delta optimum certified')

