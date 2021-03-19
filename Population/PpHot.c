//
// Created by 11577 on 2021/3/19.
//
//
// standard include
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>

// define parameters
#define L           100      /* lattice size                   */
#define SIZE        (L*L)    /* number of sites                */
#define MC_STEPS    50000   /* run-time in MCS     */
#define K           0.1     /* temperature */
#define RUN    1
#define IN     4
#define alpha1 2.0
#define alpha2 5.0


int steps;
int defector, cooperator;

double b;
double u;

typedef int       tomb1[SIZE];
typedef long int  tomb3[SIZE][IN];

typedef double    tomb6[SIZE];
typedef double    tomb9[MC_STEPS];

tomb1 player_s1;            /* matrix, containing player's strategies: 0 (C) & 1(D) */
tomb3 player_n1;            /* matrix, containing players neighbours */
tomb1 player_tp;


tomb9 each_p;


void prodgraph(void);
void initial(void);
void game(void);
void tongji(void);


FILE *outfile2;

/******************************************** RNG procedures ***********************************************************/
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
    int i;
    double ran_p,pc;

    for(i = 0; i < SIZE; i++)
    {
        pc=0.5;
        ran_p=randf();
        if (ran_p <= pc)
        {
            player_s1[i] = 0;
        } else {
            player_s1[i] = 1;
        }
    }
    for(i = 0; i < SIZE; i++)
    {
        player_tp[i]=(int)randi(2);		// 分成两个种群0、1
        //player_tp[i]=1;   //验证
    }
}



// creates first a square grid graph and then rewires Q links
void prodgraph(void)
{
    int iu, ju;
    long int player1,player2;
    int i,j;

    //int ii, jj, k;

    // set up an initial square lattice, 4 neighborhood
    for(i=0; i<L; i++)
    {
        for(j=0; j<L; j++)
        {
            // the first player
            player1 = L * j + i;

            // and its four nearest neighbors
            iu = i + 1;
            ju = j;
            if (iu==L) iu = 0;
            player2 = L * ju + iu;
            player_n1[player1][0] = player2;

            iu = i;
            ju = j + 1;
            if (ju==L) ju = 0;
            player2 = L * ju + iu;
            player_n1[player1][1] = player2;

            iu = i - 1;
            ju = j;
            if (i==0) iu = L - 1;
            player2 = L * ju + iu;
            player_n1[player1][2] = player2;

            iu = i;
            ju = j - 1;
            if (j==0) ju = L - 1;
            player2 = L * ju + iu;
            player_n1[player1][3] = player2;
        }
    }
}


double calc_payoff(int player1)
{
    int i;
    int player2;
    double payoff,bp;
    payoff=0.0;

    for(i=0; i<IN; i++)
    {
        player2 = player_n1[player1][i];

        if(player_tp[player1] == player_tp[player2]) {
            if(player_s1[player1]==0 && player_s1[player2]==0)
                payoff += 1;
            else if(player_s1[player1]==0 && player_s1[player2]==1)
                payoff += 0;
            else if(player_s1[player1]==1 && player_s1[player2]==0)
                payoff += b;
            else
                payoff += 0;
        }else{
            payoff += 0;
        }
    }

    return payoff;

}



void game(void)
{
    int i, j,k,m1,m2;
    int strat1,strat2;
    double U1,U2;
    int player1,player2,type1,type2;
    int suiji;

    double p,dP,dP_h;
    double ran_p;
    double alpha;

    for(i=0;i<SIZE;i++)
    {
        player1 = (int) randi(SIZE);
        strat1 = player_s1[player1];
        U1=calc_payoff(player1);
        type1 = player_tp[player1];

        suiji= (int) randi(IN);
        player2 = player_n1[player1][suiji];
        strat2 = player_s1[player2];
        U2=calc_payoff(player2);
        type2 = player_tp[player2];

        if(strat1!=strat2)   //考虑理性  非理性
            //if(U1<=U2)
        {
            dP = U1 - U2;

            alpha = type1 == 0 ? alpha1 : alpha2;
            dP_h = U1 - alpha;

            p = (double)((1-u)/(1 + exp(dP/K)) + u/(1 + exp(dP_h/K)));
            //p = 1/(1 + exp(dP/K));
            ran_p=randf();
            if(ran_p<=p)
            {
                player_s1[player1]=strat2;
                player_tp[player1] = type2;
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



int main()
{

    double aa,x,XX;

    char na[25];
    char fn[85];

    int MC;

    MC=MC_STEPS-2001;

    sgenrand(time(NULL));
    prodgraph();

    printf("=============start===============\n");
    printf("Hope=%.1f~%.1f\n",alpha1,alpha2);
    strcpy(fn, "Hot");
    sprintf(na, "_Hope=%.1f~%.1f\n",alpha1,alpha2);
    strcat(fn, na);
    outfile2=fopen(fn,"w+");

    for(b=1.0; b<=2.01; b=b+0.05)
    {


        if(outfile2==NULL)
        {
            printf("can not open the file for writing!");
            abort();
        }

        for(u=0.0; u<=1.01; u=u+0.05)        // 10 independent runs
        {
            aa=0;

            initial();

            for (steps=0; steps<MC_STEPS; steps++)
            {
                game();
                tongji();
                if(steps>MC)
                {
                    x=(double)cooperator/SIZE;
                    aa+=x;
                }
            }
            XX=aa/2000;
            fprintf(outfile2, "%.2f\t %.2f\t %f\n",b, u,XX);
            printf("%.2f\t %.2f\t %f\n",b, u, XX);
        }

    }
    fclose(outfile2);
    fn[0]='\0';

    printf("=============end===============\n");

    return 0;
}





