from joint_search import *
rng=np.random.default_rng(6767002)
data=json.loads((HERE/'parent_curves.json').read_text());cat=json.loads((HERE/'pair_catalogue.json').read_text())
jseeds=[mk for mk in []];mseeds=[]
def mkpair(ww):
    (r,b),(s,c)=ww;return [(r+s)/2,(b+c)/2,(c-b)/(s-r),(s-r)/2]
for name,term,collection in [('J',2,jseeds),('M',3,mseeds)]:
    collection.extend(mkpair([q[:2] for q in d['endpoints']]) for d in data[name][::20])
    for d in cat[name]:
        if max(abs(d['responses'][1][i]-d['responses'][0][i]) for i in [0,1,term])<1e-5 and min(q[1] for q in d['endpoints'])>1e-5:
            collection.append(mkpair(d['endpoints']))
# Draw more initializations close to the positive J-arc seed than the exhausted negative J tube.
positive=[z for z in jseeds if source(*worlds(np.r_[z,mseeds[0]])[0])[0]>0]
print('seeds',len(jseeds),len(mseeds),'positive',len(positive),flush=True)
out={'status':'numerical candidates only','starts':0,'roots':[],'best':[]};best=1e10
for i in range(3000):
    if i<2000:
        j=(positive if i%3 else jseeds)[rng.integers(len(positive if i%3 else jseeds))];m=mseeds[rng.integers(len(mseeds))]
        z=np.r_[j,m];z=z*(1+rng.normal(0,[.07,.3,.3,.15,.07,.3,.3,.15]))
    else:
        def randpair():
            R=rng.uniform(.05,.95);h=rng.uniform(.005,min(R,1-R)*.95);B=10**rng.uniform(-5,2);k=B/h*rng.uniform(-.99,.99);return [R,B,k,h]
        z=np.r_[randpair(),randpair()]
    try:
        sol=root(equations,z,jac=lambda zz:jac(equations,zz),tol=1e-10,options={'maxfev':500});e=float(np.max(np.abs(equations(sol.x))))
        if valid(sol.x) and e<best:
            best=e;row={'start':i,'z':sol.x.tolist(),'error_scaled':e,'worlds':[list(w) for w in worlds(sol.x)],'responses':[source(*w).tolist() for w in worlds(sol.x)]};out['best'].append(row);print('BEST',row,flush=True)
        if valid(sol.x) and e<1e-8 and abs(source(*worlds(sol.x)[0])[0]-24)>1e-4:
            w=np.array(worlds(sol.x))
            if not any(np.linalg.norm(w-np.array(q['worlds']))<1e-3 for q in out['roots']):out['roots'].append(row);print('ROOT',row,flush=True)
    except (ValueError,ZeroDivisionError,OverflowError):pass
    out['starts']=i+1
    if (i+1)%500==0: print('COUNT',i+1,flush=True);(HERE/'joint_catalogue_candidates.json').write_text(json.dumps(out,indent=2)+'\n')
    if len(out['roots'])>=8:break
(HERE/'joint_catalogue_candidates.json').write_text(json.dumps(out,indent=2)+'\n')
