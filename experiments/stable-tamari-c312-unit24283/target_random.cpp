#include <algorithm>
#include <chrono>
#include <cstdint>
#include <functional>
#include <iostream>
#include <random>
#include <vector>
using namespace std; using U=uint64_t;
bool le(const vector<int>& a,const vector<int>& b){int n=a.size();for(int i=0;i<n;i++)if(a[i]>b[i])return false;for(int j=0;j<n;j++)for(int i=0;i<j;i++){bool pos=true;for(int k=i;k<j;k++)if(a[k]<=a[j]){pos=false;break;}if(pos&&a[i]-a[j]>b[i]-b[j])return false;}return true;}
int ht(const vector<int>&a,const vector<int>&b){int s=0;for(int i=0;i<a.size();i++){int p=-1;for(int j=i+1;j<a.size();j++)if(a[j]<a[i]){p=j;break;}int d=b[i]-a[i];if(p>=0)d-=b[p]-a[p];if(d<0)return -1;s+=d;}return s;}
vector<U> poly(const vector<int>&b){int D=0;for(int x:b)D+=x;vector<U> f(D+1),a(b.size());function<void(int)> go=[&](int k){if(k==b.size()){vector<int> ai(a.begin(),a.end());if(le(ai,b))f[ht(ai,b)]++;return;}for(a[k]=0;a[k]<=unsigned(b[k]);a[k]++)go(k+1);};go(0);while(f.size()>1&&f.back()==0)f.pop_back();return f;}
int main(){mt19937 rng(20261030);uniform_int_distribution<int> dn(3,8),dv(0,20);long long accepted=0;for(int z=0;z<20000;z++){int n=dn(rng);vector<int>b(n);long long box=1;for(int&i:b){i=dv(rng);box*=i+1;}if(box>100000){z--;continue;}auto f=poly(b);auto r=b;reverse(r.begin(),r.end());auto g=poly(r);accepted++;if(f!=g){cerr<<"COUNTEREXAMPLE";for(int x:b)cerr<<" "<<x;cerr<<"\n";return 1;}}cout<<"RANDOM_NO_COUNTEREXAMPLE seed=20261030 samples="<<accepted<<" max_coordinate=20 box_limit=100000\n";}
