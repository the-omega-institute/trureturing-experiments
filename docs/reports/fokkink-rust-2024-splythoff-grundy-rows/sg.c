#include <stdio.h>
#include <stdlib.h>
#include <string.h>
int main(int argc,char**argv){
  int N=atoi(argv[1]); int X=atoi(argv[2]); int GM=atoi(argv[3]);
  int *G=malloc(sizeof(int)*(N+1)*(N+1)); char *seen=malloc(4*N+10);
  #define A(x,y) G[(x)*(N+1)+(y)]
  for(int s=0;s<=2*N;s++) for(int x=(s>N?s-N:0);x<=(s<N?s:N);x++){int y=s-x;
    memset(seen,0,4*N+10);
    for(int k=1;k<=x;k++) seen[A(x-k,y)]=1;
    for(int k=1;k<=y;k++) seen[A(x,y-k)]=1;
    int m=x<y?x:y; for(int k=1;k<=m;k++) seen[A(x-k,y-k)]=1;
    if(m>0){int r=(x>y?x:y)-m; for(int c=0;c<=r;c++) seen[A(c,r-c)]=1;}
    int g=0; while(seen[g]) g++; A(x,y)=g;}
  long bad=0, checked=0; int maxratio_num=0;
  /* rows x<=X, values g<=GM: exactly one y in [0,N] with G=g, and it satisfies y<=5x+2g+3; also symmetry */
  for(int x=0;x<=X;x++){
    for(int g=0;g<=GM;g++){int cnt=0,pos=-1; for(int y=0;y<=N;y++) if(A(x,y)==g){cnt++;pos=y;}
      checked++; int bound=5*x+2*g+3;
      if(bound>N){ printf("bound exceeds N at x=%d g=%d\n",x,g); bad++; continue; }
      if(cnt!=1||pos>bound){ bad++; if(bad<10) printf("FAIL x=%d g=%d cnt=%d pos=%d bound=%d\n",x,g,cnt,pos,bound);} }
  }
  for(int x=0;x<=N;x++) for(int y=0;y<=N;y++) if(A(x,y)!=A(y,x)) {bad++; break;}
  /* injectivity in every row up to N */
  for(int x=0;x<=N;x++){ char *u=calloc(4*N+10,1); for(int y=0;y<=N;y++){ if(u[A(x,y)]){bad++;break;} u[A(x,y)]=1;} free(u);}
  printf("N=%d rows<=%d values<=%d checked=%ld bad=%ld\n",N,X,GM,checked,bad);
  return bad?1:0;
}
