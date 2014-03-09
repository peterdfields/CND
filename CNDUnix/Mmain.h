/*  
          Mmain.h,  header for Mmain.c
          
Copyright (C) 1996 Christopher J. Basten.

This file is part of CNDm. CNDmiseq is free software; you
can redistribute it and/or modify it under the terms of the GNU  General
Public License as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

CNdiseq is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
more details.

You should have received a copy of the GNU General Public License along with
CNDm; see the file COPYING.  If not, write to the Free Software
Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
*/
#include <stddef.h>
#include "utils.h"
#include "chisqr.h"

#define INTERVALS 100
#define REPS 1000
#define DX 0.1

#define NUCLOCI 20
#define CYTLOCI 10
#define INDIVID 500
#define MAXNAME 64
#define MAXALLELES 40
#define SSLB 100
#define Zalpha 1.96
#define Zalpha2 2.576
#define Zbeta 1.28

struct dataset {
  int kk;              /* Number of nuclear alleles, 2 by default */
  int mm;              /* Number of cytotypes, 2 by default */
  int nn;              /* Sample size */
  int **counts;        /* counts = imatrix(0,kk*(kk+1)/2,0,mm) */
  int **kmcounts;        /* counts = imatrix(0,kk,0,mm) */
  int *acounts;        /* acounts = ivector(1,kk) are the nuclear allelic counts... */
  char cnames[MAXNAME];        /* name of the cytoplasmic system */
  char nnames[MAXNAME];        /* name of the nuclear system */
                       /* The following are all dmatrix(0,kk*(kk+1)/2,0,mm) */
  double **diseq_ests; /* Estimates of the disequilibria */
  double **ndiseq_ests;/* Estimates of the normalized disequilibria */
  double **freq_ests;  /* Estimates of the joint and marginal frequencies */
  double *pstar_ests;  /* Estimates of P*, = 1 - sumj Pij */
  double **var0_ests;  /* Estimates of asymptotic variances under H0 */
  double **var1_ests;  /* Estimates of asymptotic variances under H1 */
  double **stat_ests;  /* Estimates of asymptotic test statistics */
  double **pr_ests;    /* Estimates of probabilities of data */
  double MonteCarlo;   /* MonteCarlo probability of whole data set */
  double Markov;       /* Markov probability of the whole data set */
};


typedef struct dataset adataset;

int colindex();
void  joint_counts();
void allele_counts();
void write_samplesizes();
void write_samplesizes_tex();
void shuffle_matrix();
double prob_sample();
void  calc_Exact_prs();
void calc_stdev();
void write_thematrix();
void write_thematrix_tex();
void write_counts();
void write_counts_tex();
adataset *datavector();
void free_datavector();
adataset *get_data();
int get_numsets();
void calc_afreq();
void calc_adiseq();
void calc_ndiseq();
void calc_avars();
void calc_astat();
void calc_aprs();
void write_asymp_results();
void write_asymp_results_tex();
void calc_MonteCarlo_prs();
void pickNswitch();
void  calc_Markov_prs();
void   pick_ints();
void switch_it();
double themin();
void write_exact_results();

int roundit();
int samplesize();
int newsamplesize();
void clear_dataset();
 
double  lngamma();
double  chiprob();
double  chiprobexact();

