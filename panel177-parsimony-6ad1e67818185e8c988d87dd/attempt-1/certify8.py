import json,itertools,math,sys,pathlib,io,zipfile
from fractions import Fraction as F
import numpy as np
from full_lp import model
ROOT=pathlib.Path(__file__).parent
B=[F(i,64) for i in [43,34,34,29,27,25,22,20]]
CAP=F(2867,11856)

def load_saved(path):
    if path.exists():return np.load(path)
    with zipfile.ZipFile(ROOT/"numerical-data.zip") as archive:
        return np.load(io.BytesIO(archive.read(path.name)))


def rational_allocation(L,path,b,cap):
    m=model(L,'O');raw=load_saved(path)['x'];xx=[F(float(x)).limit_denominator(10**7) for x in raw[:-1]]
    schedules=[];total=F(0)
    for s,st,nc in m['groups']:
        cells=[xx[st+c*(L+1):st+(c+1)*(L+1)] for c in range(nc)];mass=sum(cells[0]);assert all(x>=0 for cell in cells for x in cell);assert all(sum(cell)==mass for cell in cells)
        total+=mass
        if mass:schedules.append({'schedule':[int(i)+1 for i in s],'mass':str(mass),'cells':[[str(x) for x in cell] for cell in cells]})
    assert total==1
    laws=[];mins={i:dict(A=F(1),C=F(1),P=F(1)) for i in range(1,L+1)};slack=F(1)
    for wi,w in enumerate(m['words']):
        law=[F(0)]*(L+1)
        for row in schedules:
            s=[i-1 for i in row['schedule']];c=0 if not s else int(w[s[0]]) if len(s)==1 else 0 if w[s[0]]==0 else 1+int(w[s[1]])
            law=[x+F(y) for x,y in zip(law,row['cells'][c])]
        assert sum(law)==1;assert all(x>=0 for x in law);n=int(m['n'][wi])
        if n:
            demands=[1-b[n-1]-cap,b[n-1]-cap,1-cap];actual=[law[0],law[n],law[0]+law[n]]
            for typ,a,d in zip(['A','C','P'],actual,demands):assert a>=d;mins[n][typ]=min(mins[n][typ],a);slack=min(slack,a-d)
        else:assert law[0]>=1-cap;positive=law;slack=min(slack,law[0]-(1-cap))
        laws.append([str(x) for x in law])
    return {'L':L,'profile':[str(x) for x in b],'cap':str(cap),'allocation':schedules,'minima':{str(i):{k:str(x) for k,x in d.items()} for i,d in mins.items()},'positive_law':[str(x) for x in positive],'minimum_row_slack':str(slack),'word_laws':laws}


def exact_score():
    m=model(8,'N');raw=-load_saved(ROOT/'varied8-N.npz')['dual'];weights=[F(float(x)).limit_denominator(1000000)*1482 for x in raw];assert all(x.denominator==1 and x>=0 for x in weights);weights=[int(x) for x in weights]
    f=[[0]*9 for _ in m['words']];tokens=[];beta=[0]*8;K=sum(weights);C0=0
    for a,(wi,t) in zip(weights,m['rows']):
        if not a:continue
        n=int(m['n'][wi]);typ='T' if not n else 'ACP'[t];word=''.join(map(str,m['words'][wi]));tokens.append([typ,word,a]);E=[0] if typ in ['A','T'] else [n] if typ=='C' else [0,n]
        for z in E:f[wi][z]+=a
        if typ=='A':beta[n-1]-=a
        if typ=='C':beta[n-1]+=a
        if typ in ['A','P','T']:C0+=a
    supports={}
    for kind in ['N','O']:
        mm=model(8,kind);data=[]
        for s,st,nc in mm['groups']:
            cells=[[0]*9 for _ in range(nc)]
            for w,rewards in zip(mm['words'],f):
                c=0 if not s else int(w[s[0]]) if len(s)==1 else (0 if w[s[0]]==0 else 1+int(w[s[1]])) if kind=='O' else 2*int(w[s[0]])+int(w[s[1]])
                cells[c]=[x+y for x,y in zip(cells[c],rewards)]
            data.append({'schedule':[int(i)+1 for i in s],'cell_maxima':[max(v) for v in cells]})
        supports[kind]={'B':max(sum(d['cell_maxima']) for d in data),'all_schedules':data}
    W=C0+beta[0]+beta[1];assert K==1482 and beta==[0,26,128,168,158,150,158,150] and C0==544 and W==570 and supports['N']['B']==570 and supports['O']['B']==570
    lower=(sum(x*y for x,y in zip(beta,B))+C0-supports['N']['B'])/K;assert lower==CAP
    return {'tokens':tokens,'K':K,'beta':beta,'C0':C0,'W':W,'supports':supports,'lower_at_profile':str(lower),'Delta_N':0,'Delta_O':0}

if __name__=='__main__':
    result={'eight':rational_allocation(8,ROOT/'varied8-O.npz',B,CAP),'score':exact_score(),'deletions':[]}
    for j in range(8):
        bb=B[:j]+B[j+1:];d=rational_allocation(7,ROOT/f'delete-{j+1}-O.npz',bb,F(15,64));d['deleted_original_coordinate']=j+1;result['deletions'].append(d)
    # Same eight-allocation works for the full h interval by endpoint checking.
    low=B.copy();low[0]=F(34,64);rational_allocation(8,ROOT/'varied8-O.npz',low,CAP)
    result['head_interval']=['17/32','43/64'];result['seven_projection_cap']='15/64';result['projection_gap']=str(CAP-F(15,64));result['verification']='All nonnegative cell entries/equal owning masses/global mass/all literal A,C,P,T rows checked with fractions.Fraction; all N/O support maxima over every label and schedule checked as integers.'
    (ROOT/'exact-eight-certificate.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['eight','score','deletions']},indent=2));print('eight minima',json.dumps(result['eight']['minima']));print('score',result['score']['beta'],result['score']['K'],result['score']['supports']['N']['B']);print('positive',result['eight']['positive_law'])
