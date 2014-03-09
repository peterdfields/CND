#include <stdio.h>
#include <string.h>
#include <stddef.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>

#include "chisqr.h"
#include "utils.h"


#define INTERVALS 100
#define REPS 10000
#define DX 0.1

#define NUCLOCI 20
#define CYTLOCI 3
#define INDIVID 300
#define MAXNAME 64
#define MAXALLELES 40
#define SSLB 10
#define Zalpha 1.96
#define Zbeta 1.28
#define MAXEXACTSS 1000

extern int debug;
extern double margthresh;


struct estimators {
  double PM;
  double PA;
  double PAM;
  double PAA;
  double PAAM;
  double PAAm;
  double PAa;
  double PAaM;
  double PAam;
  double Paa;
  double PaaM;
  double Paam;
  double rho1;
  double rho2;
  int eN;
  double epower;
  double esize;
  int nAAM;
  int nAaM;
  int naaM;
  int nAAm;
  int nAam;
  int naam;
  int samplesize;
};

typedef struct estimators anestim;

/*  These are in Dstats.c  */
double *calc_prs_X();
void write_anest();
double *calc_the_vars();
int roundit();
int *calc_ss();
double *calc_diseq();
double *calc_corrs();
double *calc_sdevs();
double *calc_prs();
void calc_diseq_prime();
int samplesize();
void DoTheStats(char *outfile, char *errfile, char *inputf, char *prog,  anestim   *counts, int numsets, char **nnames,char **cnames, int texflag);

/*  These are in Ddatain.c */
int     get_the_data();
anestim  *anestvector();
void    free_anestvector();
int calc_the_est();
void clear_anestim();

/* In chisqr.c*/
double  lngamma();
double  chiprob();
double  chiprobexact();

/*   

void create_sample();
int calc_stats();
double *calc_prs_sim();
double *doxsim();
void    write_sum_file();
void    write_chr_file();

void shuffle_matrix();
double prob_sample();
double *do_monte_carlo();
*/

/*
  Subroutines to do exact sample sizes:  Sanchez et al.
*/

typedef struct enode {
  double nullprob;
  double condprob;
  struct enode *next;    
  struct enode *prev;    
} exactnode;


void  DrivePower(anestim *counts, char *outfile, char *errfile );
void DriveSanchez(anestim *counts , char *outfile, char *errfile);   
void   SanchezBanner(FILE *outf);
int    SanchezSampleSize(anestim *theest, anestim *workest, double alpha, double beta,char *outfile);
exactnode *SanchezPower(char *outfile, anestim *workest, exactnode *first, double size, double rho1, double rho2, double *thepower);
int AsympSampleSize(anestim workest,double alpha, double beta) ;
int AsympSS(double paam, double paa, double pm, double zalpha2, double zbeta);
double AsympPower(  anestim *workest , double zalpha2);
double asympow(double paam, double paa, double pm, int nn, double zalpha2);
double SanchezC( anestim workest, int iNN, int in1d, int ind1, int ind2);
double SanchezNullProb( anestim workest  );
double SanchezMaxMargProb( anestim workest, int iNN, int in1d, int ind1, int ind2);
double SanchezMargProb( anestim workest, int iNN, int in1d, int ind1, int ind2, double rho1, double rho2 );
double SanchezCubic(double dv, double c2, double c3, double c4);
double SanchezSampCondProb( anestim workest, double rho1, double rho2 , double nullprob);
void   SanchezSolveDiseq(double rho1, double rho2, double p1d, double pd1, double pd2, double *d1, double *d2);
double SanchezCubicDerivative(double dv, double c2, double c3, double c4);
void   SanchezCalcModelD( anestim *workest, double p1d, double pd1, double pd2,  double d1, double d2 );
void   SanchezShowModel(FILE *outf, anestim *workest , double d1, double d2);
void   SanchezNewtonRafeson(double mind, double maxd, double c2, double c3, double c4, double *d1, double dthresh);
struct enode *alloc_enode();
exactnode  *clean_enodes(exactnode *anode) ;
void dealloc_enodes(struct enode *anode);
void push_enode(exactnode *first, double nullprob, double condprob);
void pull_enode(exactnode *first );
double qgaus(double lb, double ub);
double standardnormal( double z );
void  SanchezTable5(anestim *counts, char *outfile, char *errfile );
void UnPrepNormFacts(int nfacts);
void PrepNormFacts(int nfacts);
void  CalcSanchezN(anestim *counts, char *outfile, char *errfile );
double Phi(double pr);
double PhiProb(double z) ;
