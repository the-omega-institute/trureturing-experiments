#include <algorithm>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>
using namespace std;
using U=uint64_t;
using V=vector<int>;
vector<V> points;
V bound;
void gen(V &v,int i){if(i==(int)bound.size()){points.push_back(v);return;}for(v[i]=0;v[i]<=bound[i];++v[i])gen(v,i+1);}
bool literal(const V&a,const V&b){int n=a.size();for(int i=0;i<n;++i)if(a[i]>b[i])return false;for(int j=1;j<n;++j)for(int i=j-1;i>=0&&a[i]>a[j];--i)if(a[i]-a[j]>b[i]-b[j])return false;return true;}
V par(const V&a){int n=a.size();V p(n,-1);for(int i=0;i<n;++i)for(int j=i+1;j<n;++j)if(a[j]<a[i]){p[i]=j;break;}return p;}
int forest(const V&a,const V&b,const V&p){int n=a.size(),h=0;for(int i=0;i<n;++i){int d=b[i]-a[i];if(d<0)return -1;if(p[i]>=0)d-=b[p[i]]-a[p[i]];if(d<0)return -1;h+=d;}return h;}
vector<int> peaks(const vector<U>&c){vector<int>r;for(int m=0;m<(int)c.size();++m){bool ok=true;for(int j=0;j<m;++j)if(c[j]>c[j+1])ok=false;for(int j=m;j+1<(int)c.size();++j)if(c[j+1]>c[j])ok=false;if(ok)r.push_back(m);}return r;}
void pv(const V&v){cout<<'[';for(int i=0;i<(int)v.size();++i){if(i)cout<<',';cout<<v[i];}cout<<']';}
void pc(const vector<U>&v){cout<<'[';for(int i=0;i<(int)v.size();++i){if(i)cout<<',';cout<<v[i];}cout<<']';}
int main(int argc,char**argv){auto begin=chrono::steady_clock::now();string mode=argv[1];for(int i=2;i<argc;++i)bound.push_back(stoi(argv[i]));int n=bound.size(),D=accumulate(bound.begin(),bound.end(),0);V v(n);gen(v,0);int box=points.size();if((mode=="single"||mode=="singlechecked")){vector<V>s;for(const V&a:points)if(literal(a,bound))s.push_back(a);points.swap(s);}int N=points.size();if(N>35000){cerr<<"IDEAL_LIMIT "<<N<<'\n';return 4;}vector<V>parents;for(const V&a:points)parents.push_back(par(a));vector<vector<int>>down(N);vector<vector<U>>f(N,vector<U>(D+1));U comparable=0,control=0;for(int q=0;q<N;++q)for(int p=0;p<=q;++p){int h=forest(points[p],points[q],parents[p]);bool le=h>=0;if((mode=="control"||mode=="checked"||mode=="singlechecked") && literal(points[p],points[q])!=le){cerr<<"ORDER_MISMATCH";return 2;}if(le){down[q].push_back(p);++f[q][h];++comparable;}}
if(mode=="control"){for(int p=0;p<N;++p){vector<int>dist(N,-1);dist[p]=0;for(int q=p+1;q<N;++q){int mx=-1;for(int r:down[q])if(r<q&&dist[r]>=0)mx=max(mx,dist[r]+1);dist[q]=mx;}for(int q=p;q<N;++q){int h=forest(points[p],points[q],parents[p]);if(dist[q]!=h){cerr<<"HEIGHT_MISMATCH ";pv(points[p]);pv(points[q]);return 3;}if(h>=0)++control;}}}
ofstream dump;if(mode=="dump")dump.open("polynomial-dump.txt");int checked=0,bad=0;vector<U>last;V lastb;U triplecount=0;bool first=true;
cout<<"{\"mode\":\""<<mode<<"\",\"bound\":";pv(bound);cout<<",\"coordinate_box\":"<<box<<",\"vertices\":"<<N<<",\"comparable_pairs\":"<<comparable<<",\"dp_height_controls\":"<<control<<",\"failures\":[";
for(int b=0;b<N;++b){if((mode=="single"||mode=="singlechecked")&&points[b]!=bound)continue;int db=accumulate(points[b].begin(),points[b].end(),0);vector<U>F(db+1);for(int q:down[b])for(int d=0;d<=db;++d)F[d]+=f[q][d];if(F[0]!=down[b].size()||F[db]!=1){cerr<<"ENDPOINT_MISMATCH";return 5;}if(mode=="dump"){for(int z:points[b])dump<<z<<",";dump<<";";for(U z:F)dump<<z<<",";dump<<"\n";}auto pk=peaks(F);++checked;triplecount+=accumulate(F.begin(),F.end(),U(0));last=F;lastb=points[b];if(pk.empty()){if(!first)cout<<',';first=false;cout<<"{\"b\":";pv(points[b]);cout<<",\"coefficients\":";pc(F);cout<<'}';++bad;}}
cout<<"],\"checked_b\":"<<checked<<",\"interval_triples\":"<<triplecount<<",\"last_b\":";pv(lastb);cout<<",\"last_coefficients\":";pc(last);cout<<",\"top_f_coefficients\":";pc(f.back());cout<<",\"elapsed_seconds\":"<<chrono::duration<double>(chrono::steady_clock::now()-begin).count()<<"}\n";return 0;}
