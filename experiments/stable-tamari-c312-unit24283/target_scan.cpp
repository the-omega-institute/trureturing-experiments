#include <algorithm>
#include <chrono>
#include <cstdint>
#include <functional>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>
using namespace std; using U=uint64_t;
bool le(const vector<int>& a,const vector<int>& b){int n=a.size();for(int i=0;i<n;i++)if(a[i]>b[i])return false;for(int j=0;j<n;j++)for(int i=0;i<j;i++){bool pos=true;for(int k=i;k<j;k++)if(a[k]<=a[j]){pos=false;break;}if(pos && a[i]-a[j]>b[i]-b[j])return false;}return true;}
int height(const vector<int>& a,const vector<int>& b){int n=a.size(),s=0;for(int i=0;i<n;i++){int p=-1;for(int j=i+1;j<n;j++)if(a[j]<a[i]){p=j;break;}int d=b[i]-a[i];if(p>=0)d-=b[p]-a[p]; if(d<0)return -1;s+=d;}return s;}
vector<U> poly(const vector<int>& b){int n=b.size(),D=0;for(int x:b)D+=x;vector<U> f(D+1);vector<int>a(n);function<void(int)> rec=[&](int k){if(k==n){if(le(a,b)){int h=height(a,b);if(h<0){cerr<<"NEGHEIGHT\n";exit(3);}f[h]++;}return;}for(a[k]=0;a[k]<=b[k];a[k]++)rec(k+1);};rec(0);while(f.size()>1&&f.back()==0)f.pop_back();return f;}
string key(const vector<int>&v){string s;for(int x:v){s+=to_string(x);s+=',';}return s;}
void gen(int n,int M,vector<int>&b,int k,vector<vector<int>>&all){if(k==n){all.push_back(b);return;}for(b[k]=0;b[k]<=M;b[k]++)gen(n,M,b,k+1,all);}
void pv(const vector<int>&v){cout<<"[";for(int i=0;i<v.size();i++){if(i)cout<<",";cout<<v[i];}cout<<"]";}
void pc(const vector<U>&v){cout<<"[";for(int i=0;i<v.size();i++){if(i)cout<<",";cout<<v[i];}cout<<"]";}
int main(int ac,char**av){if(ac<3)return 2;int n=stoi(av[1]),M=stoi(av[2]);auto st=chrono::steady_clock::now();vector<int>b(n);vector<vector<int>>all;gen(n,M,b,0,all);unordered_map<string,vector<U>> mp;mp.reserve(all.size()*2);long long total=0;for(size_t z=0;z<all.size();z++){auto f=poly(all[z]);total+=1;mp.emplace(key(all[z]),move(f));}for(auto &bv:all){auto rv=bv;reverse(rv.begin(),rv.end());auto f=mp[key(bv)],g=mp[key(rv)];if(f!=g){cout<<"COUNTEREXAMPLE b=";pv(bv);cout<<" rev=";pv(rv);cout<<" f=";pc(f);cout<<" fr=";pc(g);cout<<"\n";return 1;}}
cout<<"NO_COUNTEREXAMPLE n="<<n<<" max="<<M<<" candidates="<<all.size()<<" elapsed="<<chrono::duration<double>(chrono::steady_clock::now()-st).count()<<"\n";return 0;}
