/// implement testsuit
#include<cstdlib>
#include<cstdio>
#include<cstring>
#include<time.h>
#include<algorithm>
#include<stdint.h>
#include<cmath>
#include<map>
using std::pair;
using std::sort;
using std::min;
using std::max;
typedef uint64_t u64;
typedef uint32_t u32;
typedef int64_t i64;
typedef int32_t i32;
namespace RNG{
	struct xoshiro{
		inline u64 rotl(const u64 x, int k) {return (x << k) | (x >> (64 - k));}
		u64 s[4];
		inline u64 operator()(void) {
			const u64 result_starstar = rotl(s[1] * 5, 7) * 9;
			const u64 t = s[1] << 17;
			s[2] ^= s[0];
			s[3] ^= s[1];
			s[1] ^= s[2];
			s[0] ^= s[3];
			s[2] ^= t;
			s[3] = rotl(s[3], 45);
			return result_starstar;
		}
		void jump(void) {
			static const uint64_t JUMP[] = { 0x180ec6d33cfd0aba, 0xd5a61266f0c9392c, 0xa9582618e03fc9aa, 0x39abdc4529b1661c };
			uint64_t s0 = 0;
			uint64_t s1 = 0;
			uint64_t s2 = 0;
			uint64_t s3 = 0;
			for(int i = 0; i < sizeof JUMP / sizeof *JUMP; i++)
				for(int b = 0; b < 64; b++) {
					if (JUMP[i] & UINT64_C(1) << b) {
						s0 ^= s[0];
						s1 ^= s[1];
						s2 ^= s[2];
						s3 ^= s[3];
					}
					this->operator()();
				}
				
			s[0] = s0;
			s[1] = s1;
			s[2] = s2;
			s[3] = s3;
		}
		inline void seed(void*t){memcpy(s,t,32);}
		void vrseed(){
		}
		xoshiro(){
		}
	};
}
namespace Timer{
	u64 _rdtsc(){
		register u64 rr;
		asm(
			"rdtsc\n"
			"shlq $32,%%rdx\n"
			"addq %%rdx,%%rax\n"
			:"=a"(rr)::"rdx"
		);
		return rr;
	}
}
using namespace std;
//#include "amazn.h"
/// fast board
namespace fb2{
	// typedef u64 clock_t;
	clock_t clockst(){return clock();}
	RNG::xoshiro rnd;
	void srandit(){
		rnd.vrseed();
	}
	const char Initial[8][3]={{0,2,1},{2,0,1},{5,0,1},{7,2,1},{0,5,2},{2,7,2},{5,7,2},{7,5,2}};
	const int d[8][2]={{0,1},{1,1},{1,0},{1,-1},{0,-1},{-1,-1},{-1,0},{-1,1}};
	typedef uint8_t u8;
	const u8 type_move=0x01;
	const u8 type_shoot=0x02;
	struct state{u8 st,ch,fr,to;};
	typedef state pstate;
	struct board{
		u64 state;
		u8 pl[8];
		inline void apply(pstate sta){
			if(sta.st&type_move){
				state|=1ull<<sta.to;
				state^=1ull<<sta.fr;
				pl[sta.ch]=sta.to;
			}else{
				state|=1ull<<sta.to;
			}
		}
		inline void initialize(){
			for(int i=0;i<8;++i){
				pl[i]=Initial[i][0]*8+Initial[i][1];
				state|=1ull<<pl[i];
			}
		}
		inline int request(int t){
			for(int i=0;i<8;++i)if(pl[i]==t)return i;
			return -1;
		}
	};
	inline int popcnt(u64 t){return __builtin_popcountll(t);}
	inline int ctz(u64 t){return __builtin_ctzll(t);}
	struct node{
		u32 cnt;state sta;
		double total;
		node*st,*ed;
		u64 untr[4];// 8+8+8+8+32=64 byte.
		inline int unmoved(){return popcnt(untr[0])+popcnt(untr[1])+popcnt(untr[2])+popcnt(untr[3]);}
		inline pair<int,node*> wmax(){
			pair<int,node*>p(0,NULL);
			for(node*i=st;i;i=i->ed)p=max(p,pair<int,node*>(i->cnt,i));
			return p;
		}
	};
	u64 find(u64 brd,char pos){
		u64 ret=0;
		for(int i=0;i<8;++i){
			int dx=d[i][0],dy=d[i][1];
			int x=(pos>>3)+dx,y=(pos&7)+dy;
			while((x&7)==x&&(y&7)==y){
				int z=x*8+y;
				if((brd>>z)&1)break;
				ret|=1ull<<z;
				x+=dx,y+=dy;
			}
		}
		return ret;
	}
	inline void bfs(char*dist,board brd,int step,int ply,int plz){
		char q[64];int ql=0,qr=0;
		memset(dist,-1,64);
		for(int i=ply;i<plz;++i)dist[brd.pl[i]]=0,q[qr++]=brd.pl[i];
		while(ql!=qr){
			int v=q[ql++];
			for(int i=0;i<8;i+=step){
				int dx=d[i][0],dy=d[i][1];
				int x=(v>>3)+dx,y=(v&7)+dy;
				while((x&7)==x&&(y&7)==y){
					int z=x*8+y;
					if((brd.state>>z)&1)break;
					if(dist[z]==-1)dist[z]=dist[v]+1,q[qr++]=z;
					x+=dx,y+=dy;
				}
			}
		}
	}
	inline void bfs(char*adist,char*bdist,board brd,int step){
		bfs(adist,brd,step,0,4);
		bfs(bdist,brd,step,4,8);
	}
	inline double que(char a,char b){
		if(a==b)return 0.;
		if(a==-1)return -1;else if(b==-1)return 1;
		return a<b?1:-1;
	}
	void init_lib(char*li){
		memset(li,8,64);
		for(int i=0;i<8;++i)li[i*8]-=3,li[i*8+7]-=3;
		for(int i=0;i<8;++i)li[i]-=3,li[56+i]-=3;
		li[0]++,li[8]++,li[56]++,li[63]++;
	}
	inline double sigmoid(double t){
		if(abs(5*t)>40)return t>0.?1.:-1.;
		double u=exp(-5*t);
		return (1-u)/(1+u);
	}
	#define sam cnt
	#define win total
	char vmin(char a,char b){
		return (unsigned char)a<(unsigned char)b?a:b;
	}
	inline double eval(board brd,int g=0){
		char qu[2][64],ki[2][64];
//		char ken[8][64];
		u64 qen[8];
		bfs(qu[0],qu[1],brd,1);bfs(ki[0],ki[1],brd,2);
//		memset(qu,-1,sizeof qu);
//		memset(ki,-1,sizeof ki);
		for(int i=0;i<8;++i){
//			/*bfs(qen[i],brd,1,i,i+1),*/bfs(ken[i],brd,2,i,i+1);
			qen[i]=find(brd.state,brd.pl[i]);
/*			for(int j=0;j<64;++j){
//				qu[i>>2][j]=vmin(qu[i>>2][j],qen[i][j]);
				ki[i>>2][j]=vmin(ki[i>>2][j],ken[i][j]);
			}*/
		}
		double cqu=0.,cki=0.;
//		double mob[8];
		double pw2[64];
		pw2[0]=0.,pw2[1]=1.;
		for(int i=2;i<64;++i)pw2[i]=pw2[i-1]/2.;
//		char lib[64];
//		init_lib(lib);
//		memset(mob,0,sizeof mob);
		double ome=0.;
		for(int i=0;i<64;++i)if(!((brd.state>>i)&1)){
			cqu+=que(qu[0][i],qu[1][i]);
			cki+=que(ki[0][i],ki[1][i]);
/*			for(int j=0;j<8;++j)
				if(((qen[j]>>i)&1)&&~qu[1^(j>>2)][i])mob[j]+=pw2[ken[j][i]+1]*lib[i];*/
			if(~qu[0][i]&&~qu[1][i])ome+=pw2[abs(qu[0][i]-qu[1][i])+1];
		}
		double cmo=0.;
		for(int i=0;i<8;++i)cmo+=(i>3?-1:1)*popcnt(qen[i]);
		double ux=cqu+cki*(ome>=40?1:ome/40)+cmo*(ome>=80?.8:.8*cmo/80);
		return sigmoid(ux);
	}
	template<int size=64>
	struct Rave{
		double sum[size];
		int cnt[size];
		int tot;
		void init(){
			memset(sum,0,sizeof sum);
			memset(cnt,0,sizeof cnt);
			tot=0;
		}
		inline void upd(int t,double val){++cnt[t],++tot,sum[t]+=val;}
		inline double cal(int t,double sig){
			return sig*sum[t]/cnt[t]+sqrt(log(tot)/sum[t])*0.7;
		}
	};
	const int max_uct_size=4*1000000;
	const double W=15.;
	struct UCT{
		clock_t global_timer;
		void reset_timer(){global_timer=clockst();}
		double durat(){clock_t nowti=clockst();return double(nowti-global_timer)/CLOCKS_PER_SEC;}
		static constexpr double TimeLimit=5.00;// 0.95 seconds max.
		bool timeout(double tt=1.){return durat()>=TimeLimit*tt;}// timer
		node p   [max_uct_size];
		node*pst;
		void rest(){
			pst=p;
			reset_timer();
		}
		node*getn(){
			node*ret=pst++;
			ret->st=ret;
			ret->ed=0;
			return ret;
		}
		node*gennode(board bd,state sta,bool isroot=0){
			node*ret=getn();
			ret->sta=sta;
			if(!isroot)bd.apply(sta);
			int pcnt=0;
			int tr=(popcnt(bd.state)&1)*4;
			if((!isroot)&&(sta.st&type_move)){
				ret->untr[0]=find(bd.state,bd.pl[sta.ch]);
				pcnt=popcnt(ret->untr[0]);
				ret->untr[1]=0;
				ret->untr[2]=0;
				ret->untr[3]=0;
			}else{
				ret->untr[0]=find(bd.state,bd.pl[0|tr]);
				ret->untr[1]=find(bd.state,bd.pl[1|tr]);
				ret->untr[2]=find(bd.state,bd.pl[2|tr]);
				ret->untr[3]=find(bd.state,bd.pl[3|tr]);
				pcnt=ret->unmoved();
			}
			if(pcnt==0)ret->total=(sta.st&type_move)?-1.:1.,ret->cnt=1;
			else ret->total=0,ret->cnt=0,ret->st=0;
			return ret;
		}
		node*sele(node*t,board bd){
			node*r=t->st;
			if(r==t)return NULL;
			double ee=1.0*sqrt(log(t->cnt));
			auto ret=pair<double,node*>((r->total)/(r->cnt)+ee/sqrt(r->cnt),t->st);
			for(r=r->ed;r;r=r->ed){
				ret=max(ret,pair<double,node*>(
					(r->total)/(r->cnt)+ee/sqrt(r->cnt),r
				));
			}
			return ret.second;
		}
		inline int randr(u64*s,u64*t){
			int stl=0,sta[256];
			for(int g=0;s!=t;g+=64,++s){
				u64 u=*s;
				while(u){
					int v=ctz(u);
					u^=1ull<<v;
					sta[stl++]=g|v;
				}
			}
			return stl?sta[rnd()%stl]:-1;
		}
		double moneval(board bd,int dep){
			int r=(popcnt(bd.state)&1)*4;
			u64 mv[4];
			for(int i=0;i<4;++i)mv[i]=find(bd.state,bd.pl[i+r]);
			int mov=randr(mv,mv+4);
			if(mov==-1){
				return -1.;
			}
			if(dep==0){
				return r?-eval(bd):eval(bd);
			}
			int pl=(mov>>6)+r;mov&=63;
			int rp=bd.pl[pl];
			bd.apply((state){type_move,pl,rp,mov});
			mv[0]=find(bd.state,mov);
			int mov2=randr(mv,mv+1);
			bd.apply((state){type_shoot,pl,mov,mov2});
			double evalz=-moneval(bd,dep-1);
			return evalz;
		}
		double mcts_done(board bd,node*nod,int dep){
			if(nod->st==nod){return nod->total;}
			state u=(state){0,0,0,0};
			if(nod->sta.st&type_move){
				int mov=randr(nod->untr,nod->untr+1);
				u=(state){type_shoot,nod->sta.ch,bd.pl[nod->sta.ch],mov};
				bd.apply(u);
				++dep;
			}
			dep=dep+1>>1;
			dep=(dep&1)?3:2;
			double qu=-moneval(bd,dep);
			nod->sam++;
			nod->total+=qu;
			return qu;
		}
		double mcts_eval(board bd,node*nod,int dep){
			int rq=(popcnt(bd.state)&1)*4;
			double u=0.;
			int type=nod->sta.st&type_move;
			state r;node*nx;
			if(nod->st!=nod){
				if(nod->unmoved()){
					int mov=randr(nod->untr,nod->untr+(type?1:4));
					int ux=type?nod->sta.ch:((mov>>6)+rq);
					nod->untr[mov>>6]^=1ull<<(mov&63);
					r=(state){type?type_shoot:type_move,ux,bd.pl[ux],mov&63};
					nx=gennode(bd,r);
					nx->ed=nod->st;
					nod->st=nx;
					bd.apply(r);
					u=mcts_done(bd,nx,dep+1);
				}else{
					nx=sele(nod,bd);
					r=nx->sta;
					bd.apply(r);
					u=mcts_eval(bd,nx,dep+1);
				}
				u*=type?1:-1;
				nod->cnt++;
				nod->total+=u;
			}else u=nod->total;
			return u;
		}
		struct mover{char ch,fr,to,ar;};
		mover MCTS(board bd){
			rest();
			int cnt=0;
			node*root=gennode(bd,(state){type_shoot,0,0,0},1);
			if(root->st==root)return (mover){-1,-1,-1,-1};
			double r=popcnt(bd.state)<10?2:1;
			while(!timeout(r)){
				++cnt;if(cnt>max_uct_size-100)break;
				double r=mcts_eval(bd,root,0);
			}
			pair<int,node*>p(0,NULL);
			node*fa=NULL;
			for(node*i=root->st;i;i=i->ed){
				auto q=i->wmax();
				if(q>p)p=q,fa=i;
			}
			if(!fa)return (mover){-1,-1,-1,-1};
//			printf("iter: %d\n",cnt);
			return (mover){fa->sta.ch,fa->sta.fr,fa->sta.to,p.second->sta.to};
		}
	};
}
u32 _call(u64 board,u64 players,char*p){
	fb2::board t;
	fb2::rnd.seed(p);
	t.state=board;
	for(int i=0;i<8;++i){
		t.pl[i]=players&255;
		players>>=8;
	}
	fb2::UCT* gt=new fb2::UCT;
	fb2::UCT::mover g=gt->MCTS(t);
	delete gt;
	return (u32(g.fr)<<16)|(u32(g.to)<<8)|(u32(g.ar));
}
extern "C"{
	u32 mov(u64 board,u64 players,char*p){
		return _call(board,players,p);
	}
}
