#include <algorithm>
#include <cstdint>
#include <iostream>
#include <vector>
using namespace std;using U=__uint128_t;
struct P {vector<U>c;vector<int>b;};
void pu(U x){if(x>=10)pu(x/10);cout<<char('0'+x%10);}
int main(){int n;cin>>n;vector<P>a(n);for(auto &p:a){int k;cin>>k;p.c.resize(k);for(auto &x:p.c){uint64_t t;cin>>t;x=t;}cin>>k;p.b.resize(k);for(int&x:p.b)cin>>x;}uint64_t count=0;for(int i=0;i<n;i++)for(int j=i;j<n;j++){vector<U>c(a[i].c.size()+a[j].c.size()-1);for(int x=0;x<a[i].c.size();x++)for(int y=0;y<a[j].c.size();y++)c[x+y]+=a[i].c[x]*a[j].c[y];bool down=false,bad=false;for(int k=0;k+1<c.size();k++){if(c[k]>c[k+1])down=true;if(down&&c[k]<c[k+1])bad=true;}count++;if(bad){cout<<"{\"status\":\"candidate\",\"b\":[";bool first=true;for(int x:a[i].b){if(!first)cout<<',';first=false;cout<<x;}cout<<",0";for(int x:a[j].b)cout<<','<<x;cout<<"],\"coefficients\":[";for(int k=0;k<c.size();k++){if(k)cout<<',';pu(c[k]);}cout<<"],\"products_checked\":"<<count<<"}\n";return 0;}}cout<<"{\"status\":\"no-failure\",\"unique_input_polynomials\":"<<n<<",\"products_checked\":"<<count<<"}\n";}
