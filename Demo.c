// Prisoners Dilemma game on a small-world graph constructed from a square lattice
// Some players are blocked to give their strategy (other players cannot adopt their strategy)
//--------------鍚?  姝?-----------------
// standard include
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>

// define priority classes
#define NORMAL_PRIORITY_CLASS       0x00000020
#define IDLE_PRIORITY_CLASS         0x00000040
#define HIGH_PRIORITY_CLASS         0x00000080
#define REALTIME_PRIORITY_CLASS     0x00000100

// define parameters
#define L           100      /* lattice size                   */
#define SIZE        (L*L)    /* number of sites                */
#define MC_STEPS    50000   /* run-time in MCS     */
#define K           0.1     /* temperature */
#define Q           0      /* Q portion of links are rewired */
#define REPUTATION  100
#define ALPHA       1
#define NAMEOUT     "K4b075r5Q2"
#define RANDOMIZE   3145215
#define IN 4
#define ZH 70

#define derta 0.5
#define Derta 1

int defector,cooperator;
double r;
int steps;

double pc;	// initial proportion of cooperators(C)
double pa;	// initial proportion of A-type players
double w;	// transfer probability: 1 for A and w(<1.0) for B
double fact,cost;


typedef int       tomb1[SIZE];
typedef long int  tomb3[SIZE][IN];
typedef double    tomb6[SIZE];
typedef double    tomb9[MC_STEPS];

tomb1 player_s1;          /* matrix, containing player's strategies: 0 (C) & 1(D) */
tomb1 player_s1m;
tomb1 player_r1;          //reputation
tomb1 player_r1m;
tomb3 player_n1;

tomb6 Si;
tomb6 Sim;
tomb9 each_p;
/* matrix, containing players neighbours */

void prodgraph(void);      	/* creates host graph */
void initial();				/* initial state */
void game(void);
void tongji(void);
void c_initial(void);
double c_stra(int,int);

FILE *outfile1;
FILE *outfile2;

//鈥溾€樑撀犅€脢陋藱聽藵鈮に欌€λ櫰捖Ｃ该埪Ｂㄢ墹陋鈥濃垰蟺鈥姑€赂,梅卤惟鈥濃€濃垰忙脮鈥撯€撀Ｂㄢ€濃垰randf()酶鈥︹€溾€樏仿蔽┾€濃墹藱鈥λ?-1卢藱鈼娾€灻λ樷€樎烩垜梅鈮ぢ郝灯捗€脢陋藱聽藵拢篓randi(x),鈮に欌€λ?---x-1碌茠脌脢陋藱鈥櫵毬犓?/******************************************** RNG procedures ***********************************************************/
#define NN 624
#define MM 397
#define MATRIX_A 0x9908b0df   /* constant vector a */
#define UPPER_MASK 0x80000000 /* most significant w-r bits */
#define LOWER_MASK 0x7fffffff /* least significant r bits */
#define TEMPERING_MASK_B 0x9d2c5680
#define TEMPERING_MASK_C 0xefc60000
#define TEMPERING_SHIFT_U(y)  (y >> 11)
#define TEMPERING_SHIFT_S(y)  (y << 7)
#define TEMPERING_SHIFT_T(y)  (y << 15)
#define TEMPERING_SHIFT_L(y)  (y >> 18)

static unsigned long mt[NN]; /* the array for the state vector  */
static int mti=NN+1; /* mti==NN+1 means mt[NN] is not initialized */
void sgenrand(unsigned long seed)
{
    int i;
    for (i=0;i<NN;i++)
    {
        mt[i] = seed & 0xffff0000;
        seed = 69069 * seed + 1;
        mt[i] |= (seed & 0xffff0000) >> 16;
        seed = 69069 * seed + 1;
    }
    mti = NN;
}

void lsgenrand(unsigned long seed_array[])
{
    int i;
    for (i=0;i<NN;i++)
        mt[i] = seed_array[i];
    mti=NN;
}

double genrand()
{
    unsigned long y;
    static unsigned long mag01[2]={0x0, MATRIX_A};
    if (mti >= NN)
    {
        int kk;
        if (mti == NN+1) sgenrand(4357);
        for (kk=0;kk<NN-MM;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+MM] ^ (y >> 1) ^ mag01[y & 0x1];
        }
        for (;kk<NN-1;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+(MM-NN)] ^ (y >> 1) ^ mag01[y & 0x1];
        }
        y = (mt[NN-1]&UPPER_MASK)|(mt[0]&LOWER_MASK);
        mt[NN-1] = mt[MM-1] ^ (y >> 1) ^ mag01[y & 0x1];
        mti = 0;
    }
    y = mt[mti++]; y ^= TEMPERING_SHIFT_U(y); y ^= TEMPERING_SHIFT_S(y) & TEMPERING_MASK_B;
    y ^= TEMPERING_SHIFT_T(y) & TEMPERING_MASK_C; y ^= TEMPERING_SHIFT_L(y);
    return y;
}

double randf(){ return ( (double)genrand() * 2.3283064370807974e-10 ); }
long randi(unsigned long LIM){ return((unsigned long)genrand() % LIM); }
/******************************************** END of RNG **************************************************/


void initial(void)
{
    int i,j;
    double ran_p;

    for(i=0;i<SIZE;i++)
    {
        Si[i]=0.5;
        player_r1[i] = randi(REPUTATION) + 1;
        ran_p=randf();
        if (ran_p <= Si[i])
        {
            player_s1[i] = 0;
            //cooperator++;
        } else {
            player_s1[i] = 1;
            //defector++;
        }

    }

}

void c_initial(void)
{
    int i;
    double ran_p;

    for(i=0;i<SIZE;i++)
    {
        ran_p=randf();
        if (ran_p <= Si[i])
        {
            player_s1[i] = 0;
        } else {
            player_s1[i] = 1;
        }
    }
}


// creates first a square grid graph and then rewires Q links
void prodgraph(void)
{
    int iu, ju;
    long int player1,player2;
    int i,j;
    int ii, jj, k;
    for(i=0; i<L; i++)           // 4 neighborhood
        for(j=0; j<L; j++)
        {
            // the first player
            player1 = L * j + i;

            k=0;
            for(ii=-1; ii<=1; ii++)
                for(jj=-1; jj<=1; jj++)
                {
                    if((ii!=0||jj!=0)&&(ii==0||jj==0))
                    {
                        iu = (i+ii+L)%L;
                        ju = (j+jj+L)%L;
                        player2 = L*ju + iu;
                        player_n1[player1][k]=player2;
                        k++;
                    }
                }
            //printf("k=%d, produce!\n", k);
        }
}


double calc_payoff(int player1)
{
    int i;
    int player2;
    double payoff;
    payoff=0.0;

    for(i=0; i<IN; i++)
    {
        player2 = player_n1[player1][i];
        if(player_s1[player1]==0 && player_s1[player2]==0)
            payoff += 1;
        else if(player_s1[player1]==0 && player_s1[player2]==1)
            payoff += -r;
        else if(player_s1[player1]==1 && player_s1[player2]==0)
            payoff += 1+r;
        else
            payoff += 0;
    }

    return payoff;

}

double c_stra(int player1,int player2)
{

    //----------Sx------------
    double Sx;
    if(player_r1[player2]>=ZH)
    {
        Sx=derta;
    }else{
        Sx=0.5*derta;
    }
    //------------------------

    if (player_s1[player2] == 0)
    {
        Si[player1] += Sx;
    }
    else
    {
        Si[player1] -= Sx;
    }

    if(Si[player1]>=1){
        Si[player1]=1;
    }else if(Si[player1]<=0){
        Si[player1]=0;
    }else{
        Si[player1]=Si[player1];
    }
    return Si[player1];
}

void game(void)     //脮炉卢脕梅鈥撁€藰鈥濃€撯垙藛脙脗惟炉鈥撯€撯墹铿偮€樷€斆熍撯垶
{
    int i;
    int strat1,strat2;
    double Ux,Uy;
    int reputation1,reputation2;
    int type1,type2;
    int player1,player2;
    int suiji;            //脌脢陋藱聽藵

    double p,dP;
    double qqqq;
    double wy;


    for(i=0;i<SIZE;i++)
    {
        player1 = (int) randi(SIZE);
        strat1 = player_s1[player1];
        reputation1=player_r1[player1];

        Ux=calc_payoff(player1);

        suiji= (int) randi(IN);
        player2 = player_n1[player1][suiji];
        strat2 = player_s1[player2];
        reputation2=player_r1[player2];

        Uy=calc_payoff(player2);

        //if(strat1!=strat2)
        if(Ux<=Uy)
        {
            dP=Ux-Uy;
            wy=pow(((reputation2+0.0)/REPUTATION),ALPHA);

            p=wy/(1+exp(dP/K));
            qqqq=randf();
            if(qqqq<p){
                player_s1[player1]=strat2;
                Si[player1]=c_stra(player1,player2);
            }
        }
        strat1=player_s1[player1];
        if(strat1==0 && reputation1<100){
            player_r1[player1] += Derta;
            if(player_r1[player1]>100){
                player_r1[player1]=100;
            }
        }
        else if(strat1==1 && reputation1>1){
            player_r1[player1] -= Derta;
            if(player_r1[player1]<1){
                player_r1[player1]=1;
            }
        }
    }

}

void tongji(void)
{
    int i;
    cooperator=0;
    defector=0;
    for(i=0;i<SIZE;i++)
    {
        if(player_s1[i]==0){
            cooperator++;
        }
        else{
            defector++;
        }
    }

}

void each(void)
{
    int i;
    for(i=0;i<MC_STEPS;i++)
    {
        each_p[i]=0.0;
    }
}


// the main program
int main()
{

    int i;
    double aa,x,XX,sum,ave_p;

    char na[25];
    char fn[85];

    int run;

    sgenrand(time(NULL));
    prodgraph();
    each();

    printf("=============start===============\n");
    printf("d=%.1f_D=%d_ALPHA=%d\n",derta, Derta, ALPHA);

    strcpy(fn, "RT");
    sprintf(na, "_d=%.1f_D=%d_ALPHA=%d.txt",derta, Derta, ALPHA);
    strcat(fn, na);
    outfile2=fopen(fn,"w+");

    for(run=1; run<=1; run++)        // 10 independent runs
    {
        printf("run=%d\n",run);
        if(outfile2==NULL)
        {
            printf("can not open the file for writing!");
            abort();
        }

        initial();
        r=0.3;


        for (steps=0; steps<MC_STEPS; steps++)
        {
            tongji();
            game();
            c_initial();
            if(steps%1==0)
            {
                x=(double)cooperator/SIZE;
                each_p[steps] += x;
                printf("%d\t %f\n", steps, x);
            }

        }


    }
    printf("-----------Average-------------\n");
    for (i=0;i<MC_STEPS;i++)
    {
        ave_p=each_p[i]/1;
        fprintf(outfile2, "%d\t %f\n", i+1, ave_p);
        printf("%d\t %f\n", i, ave_p);
    }

    fclose(outfile2);
    fn[0]='\0';

    printf("=============end===============\n");

    return 0;
}
