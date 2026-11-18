/* Literal interpreter for the five Boolean-stack instructions.
 * This is a reference machine, not the paid LAB runtime or an evaluator oracle.
 * Build: cc -O3 -std=c11 constructor_reference.c -o <temporary executable>
 */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <string.h>
#include <time.h>
typedef struct {unsigned char *p; size_t n,cap;} Word;
static void fail(const char*s){fprintf(stderr,"%s\n",s);exit(1);}
static void push(Word*w,unsigned char b){if(w->n==w->cap){w->cap=w->cap?2*w->cap:16;w->p=realloc(w->p,w->cap);if(!w->p)fail("memory");}w->p[w->n++]=b;}
static uint32_t read32(FILE*f){unsigned char b[4];if(fread(b,1,4,f)!=4)fail("short code");return (uint32_t)b[0]|(uint32_t)b[1]<<8|(uint32_t)b[2]<<16|(uint32_t)b[3]<<24;}
int main(int argc,char**argv){
 if(argc!=5){fprintf(stderr,"usage: reference graph.bin input.bits output.bits primitive-fuel\n");return 1;}
 FILE*f=fopen(argv[1],"rb");if(!f)fail("graph open");char magic[8];if(fread(magic,1,8,f)!=8||memcmp(magic,"RROSTACK",8))fail("magic");
 uint32_t nr=read32(f),nc=read32(f),q=read32(f);uint32_t(*code)[6]=calloc(nc,sizeof(*code));Word*r=calloc(nr,sizeof(*r));if(!code||!r||!nc||!nr||q>=nc)fail("header");
 for(uint32_t i=0;i<nc;i++)for(int j=0;j<6;j++)code[i][j]=read32(f);if(fgetc(f)!=EOF)fail("extra code");fclose(f);
 for(uint32_t i=0;i<nc;i++){uint32_t*x=code[i];if(x[0]>4)fail("opcode");if(x[0]>=2&&x[1]>=nr)fail("register");if(x[0]==1&&x[2]>=nc)fail("goto");if(x[0]==2&&(x[2]>1||x[3]>=nc))fail("push");if(x[0]>=3)for(int j=2;j<5;j++)if(x[j]>=nc)fail("branch");}
 f=fopen(argv[2],"rb");if(!f)fail("input open");int c;while((c=fgetc(f))!=EOF){if(c!='0'&&c!='1')fail("input alphabet");push(&r[0],c-'0');}fclose(f);
 for(size_t i=0;i<r[0].n/2;i++){unsigned char b=r[0].p[i];r[0].p[i]=r[0].p[r[0].n-1-i];r[0].p[r[0].n-1-i]=b;}
 uint64_t fuel=strtoull(argv[4],0,10),steps=0,total=r[0].n,peak=total,pushes=0,pops=0,peeks=0,jumps=0;int success=-1;clock_t begin=clock();
 while(steps<fuel){uint32_t*x=code[q];steps++;Word*w=&r[x[1]];switch(x[0]){
 case 0:success=x[2];goto done;
 case 1:q=x[2];jumps++;break;
 case 2:push(w,x[2]);total++;if(total>peak)peak=total;pushes++;q=x[3];break;
 case 3:case 4:{uint32_t b=w->n?w->p[w->n-1]+3:2;if(x[0]==3){pops++;if(w->n){w->n--;total--;}}else peeks++;q=x[b];break;}
 default:fail("opcode");}}
 done:
 if(success<0){fprintf(stderr,"fuel exhausted steps=%"PRIu64" state=%u\n",steps,q);return 2;}
 f=fopen(argv[3],"wb");if(!f)fail("output open");for(size_t i=r[0].n;i>0;i--)fputc('0'+r[0].p[i-1],f);fclose(f);
 printf("{\"success\":%s,\"steps\":%"PRIu64",\"pushes\":%"PRIu64",\"pop_instructions\":%"PRIu64",\"peeks\":%"PRIu64",\"gotos\":%"PRIu64",\"peak_stack_bits\":%"PRIu64",\"retained_stack_bits\":%"PRIu64",\"output_bits\":%zu,\"halt_state\":%u,\"cpu_seconds\":%.6f}\n",success?"true":"false",steps,pushes,pops,peeks,jumps,peak,total,r[0].n,q,(double)(clock()-begin)/CLOCKS_PER_SEC);
 return success?0:1;
}
