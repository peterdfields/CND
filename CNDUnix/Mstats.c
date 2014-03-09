/*  
          Mstats.c,  subroutines to calculate asymptotic statistics for
          CNDm, a program to calculate
          cytonuclear disequilibria for a multiallelic system.
          
Copyright (C) 1996 Christopher J. Basten.

This file is part of CNDm. CNDm is free software; you
can redistribute it and/or modify it under the terms of the GNU  General
Public License as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

CNDm is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
more details.

You should have received a copy of the GNU General Public License along with
CNDm; see the file COPYING.  If not, write to the Free Software
Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
*/
#include "Mmain.h"
void check_counts();
char Exact(int [][2], float *both_tail, float *closest);

/*
  Compute the allelic counts from the joint cytonuclear genotypic array
*/
void  allele_counts(thedata,numsets,errorfile)
  adataset *thedata;
  int numsets;
  char *errorfile;
{
  int iset,ii,jj,rows,cols,colctr;
  FILE *fileptr;
  for ( iset = 1 ; iset <= numsets ; iset++ ) {
    colctr = 0;
    rows = (thedata+iset)->mm;
    cols = (thedata+iset)->kk * ( (thedata+iset)->kk + 1 ) / 2 ;
    (thedata+iset)->acounts = ivector(1,(thedata+iset)->kk);
    for ( ii = 1 ; ii <= (thedata+iset)->kk ; ii++ )
      *((thedata+iset)->acounts+ii) = 0;
    if ( (thedata+iset)->nn > 0 )
      for ( ii = 1 ; ii <= (thedata+iset)->kk ; ii++ )
        for ( jj = ii ; jj <= (thedata+iset)->kk ; jj++ ) {
          colctr = colctr+1;
          *( (thedata+iset)->acounts + ii )  =  *( (thedata+iset)->acounts+ii)+*(*( (thedata+iset)->counts +0)+colctr);
          *( (thedata+iset)->acounts + jj )  =  *( (thedata+iset)->acounts+jj)+*(*( (thedata+iset)->counts +0)+colctr);
        }
    else {
      fileptr = fileopen(errorfile,"a");
      fprintf(fileptr,"\n\nallele_counts:  Sample size for dataset %3d is <0...\n",iset);
      fileclose(errorfile,fileptr);
    }
  } 
}

/*
  Compute the counts of nik, ni, and 2nk from the joint cytonuclear genotypic array.
*/
void  joint_counts(thedata,numsets,errorfile)
  adataset *thedata;
  int numsets;
  char *errorfile;
{
  int iset,ii,jj,kk,cols,colctr;
  FILE *fileptr;
  for ( iset = 1 ; iset <= numsets ; iset++ ) {
    (thedata+iset)->kmcounts = imatrix(0,(thedata+iset)->mm,0,(thedata+iset)->kk);
    (thedata+iset)->pstar_ests = dvector(1,(thedata+iset)->kk);
    for ( ii = 0 ; ii <= (thedata+iset)->mm ; ii++ )
      for ( jj = 0 ; jj <= (thedata+iset)->kk ; jj++ )
        (thedata+iset)->kmcounts[ii][jj] = 0;
    if ( (thedata+iset)->nn > 0 ) 
      for ( kk = 1 ; kk <= (thedata+iset)->mm ; kk++ )  {
        colctr = 0;
        for ( ii = 1 ; ii <= (thedata+iset)->kk ; ii++ )
          for ( jj = ii ; jj <= (thedata+iset)->kk ; jj++ ) {
            colctr = colctr+1;
            (thedata+iset)->kmcounts[0][0] = (thedata+iset)->kmcounts[0][0] + 2 * (thedata+iset)->counts[kk][colctr];
            (thedata+iset)->kmcounts[kk][ii] = (thedata+iset)->kmcounts[kk][ii] + (thedata+iset)->counts[kk][colctr];
            (thedata+iset)->kmcounts[kk][jj] = (thedata+iset)->kmcounts[kk][jj] + (thedata+iset)->counts[kk][colctr];
            (thedata+iset)->kmcounts[0][ii] = (thedata+iset)->kmcounts[0][ii] + (thedata+iset)->counts[kk][colctr];
            (thedata+iset)->kmcounts[0][jj] = (thedata+iset)->kmcounts[0][jj] + (thedata+iset)->counts[kk][colctr];
            (thedata+iset)->kmcounts[kk][0] = (thedata+iset)->kmcounts[kk][0] + 2 * (thedata+iset)->counts[kk][colctr];
            (thedata+iset)->pstar_ests[ii]  = (thedata+iset)->pstar_ests[ii] + (double) (thedata+iset)->counts[kk][colctr] / (double) (thedata+iset)->nn ;
          }
       }
    else {
      fileptr = fileopen(errorfile,"a");
      fprintf(fileptr,"\n\njoint_counts:  Sample size for dataset %3d is <0...\n",iset);
      fileclose(errorfile,fileptr);
    }   
  } 
}

/*
  Calculate the asymptotic frequencies.  If ag == 0, then do it for the 
  genotypic values (Pijk), otherwise do it for the allelic values (Pik).
*/
void  calc_afreq(thedata,numsets,errorfile,ag)
  adataset *thedata;
  int numsets;
  char *errorfile;
  int ag;
{
  int iset,ii,jj,rows,cols,**mptr;
  double ss;
  FILE *fileptr;
  for ( iset = 1 ; iset <= numsets ; iset++ ) {
    rows = (thedata+iset)->mm;
    if ( ag == 0 ) {
      ss = (thedata+iset)->nn ;
      cols = (thedata+iset)->kk * ( (thedata+iset)->kk + 1 ) / 2 ;
      mptr = (thedata+iset)->counts;
    }
    else {
      ss = 2 * (thedata+iset)->nn ;
      cols = (thedata+iset)->kk ;
      mptr = (thedata+iset)->kmcounts;
    }
    if ( (thedata+iset)->freq_ests == NULL )
      (thedata+iset)->freq_ests = dmatrix(0,rows,0,(thedata+iset)->kk * ( (thedata+iset)->kk + 1 ) / 2 );
    if ( (thedata+iset)->nn > 0 )
      for ( ii = 0 ; ii <= rows ; ii++ )
        for ( jj = 0 ; jj <= cols ; jj++ )
          *(*( (thedata+iset)->freq_ests + ii ) + jj) = (double) *(*(mptr+ii)+jj) / ss ;
    else {
      fileptr = fileopen(errorfile,"a");
      fprintf(fileptr,"\n\ncalc_afreq:  Sample size for dataset %3d is<0...\n",iset);
      fileclose(errorfile,fileptr);
    }
  } 
}


/*
  Calculate the asymptotic disequilibria.  If ag == 0, then do it for the 
  genotypic values (Pijk), otherwise do it for the allelic values (Pik).
*/

void  calc_adiseq(thedata,numsets,errorfile,ag)
  adataset *thedata;
  int numsets;
  char *errorfile;
  int ag;
{
  int iset,ii,jj,rows,cols,colctr;
  double ss,twonn;
  FILE *fileptr;
  for ( iset = 1 ; iset <= numsets ; iset++ ) {
    colctr = 0;
    twonn = 2.0 *  (double) *(*((thedata+iset)->counts+0)+0);
    twonn = 1.0/twonn;
    rows = (thedata+iset)->mm;
    if ( ag == 0 )
      cols = (thedata+iset)->kk * ( (thedata+iset)->kk + 1 ) / 2 ;
    else
      cols = (thedata+iset)->kk;
    if ( (thedata+iset)->diseq_ests == NULL )
      (thedata+iset)->diseq_ests = dmatrix(0,rows,0, (thedata+iset)->kk * ( (thedata+iset)->kk + 1 ) / 2);
    else
      for ( ii = 0; ii <= rows ; ii++ )
        for ( jj = 0 ; jj <= (thedata+iset)->kk * ( (thedata+iset)->kk + 1 ) / 2 ; jj++ )
           (thedata+iset)->diseq_ests[ii][jj] = 0.0;
    if ( (thedata+iset)->nn > 0 ) {
      for ( ii = 1 ; ii <= rows ; ii++ )
        for ( jj = 1 ; jj <= cols ; jj++ )
           (thedata+iset)->diseq_ests[ii][jj] =  (thedata+iset)->freq_ests[ii][jj]  -   (thedata+iset)->freq_ests[0][jj] * (thedata+iset)->freq_ests[ii][0]  ;
      if ( ag == 0 ) 
        for ( ii = 1 ; ii <=  (thedata+iset)->kk ; ii++ )
          for ( jj = ii ; jj <= (thedata+iset)->kk ; jj++ ) {
            colctr = colctr+1;
            if ( ii != jj )
              *(*( (thedata+iset)->diseq_ests + 0)+colctr) =(double)*((thedata+iset)->acounts+ii) * twonn*(double)*((thedata+iset)->acounts+jj)*twonn  - *(*((thedata+iset)->freq_ests+0)+colctr)/ 2.0 ;
            else
              *(*( (thedata+iset)->diseq_ests + 0)+colctr) = -(double)*((thedata+iset)->acounts+ii) * twonn *(double)*((thedata+iset)->acounts+jj) *twonn  + *(*((thedata+iset)->freq_ests+0)+colctr);
          }
    }
    else {
      fileptr = fileopen(errorfile,"a");
      fprintf(fileptr,"\n\ncalc_adiseq: Sample size for dataset %3d is < 0...\n",iset);
      fileclose(errorfile,fileptr);
    }
  } 
}

/*
  Get the column index from the alleles.  Assumes that
  11 12 13 .... 1r
     22 23 .... 2r
        33 .... 3r
           .....rr
  is mapped onto a vector in the order
  
    11 12 13 ... 1r 22 23 .... 2r 33 .... 3r ... rr

  col goes from 1 to r(r+1)/2
  
  ij is the genotype.
  r is the number of alleles.

*/
int colindex(i,j,r)
  int i,j,r;
{
  int col;
  col = j-i + 1 + (i-1)*(r+1) - i*(i-1)/2;
  return(col);
}

/*
  Calculate the normalized asymptotic disequilibria.  If ag == 0, then do it for the 
  genotypic values (Dijk), otherwise do it for the allelic values (Dik).
*/

void  calc_ndiseq(thedata,numsets,errorfile,ag)
  adataset *thedata;
  int numsets;
  char *errorfile;
  int ag;
{
  int iset,ii,jj,rows,cols,colctr;
  double twonn,ub1,ub2,lb1,lb2,ub3,lb3,pij,pk,pijk,dijk,pii,pstar;
  FILE *fileptr;
  for ( iset = 1 ; iset <= numsets ; iset++ ) {
    colctr = 0;
    twonn = 2.0 *  (double) *(*((thedata+iset)->counts+0)+0);
    twonn = 1.0/twonn;
    rows = (thedata+iset)->mm;
    if ( ag == 0 )
      cols = (thedata+iset)->kk * ( (thedata+iset)->kk + 1 ) / 2 ;
    else
      cols = (thedata+iset)->kk;
    if ( (thedata+iset)->ndiseq_ests == NULL )
      (thedata+iset)->ndiseq_ests = dmatrix(0,rows,0,(thedata+iset)->kk * ( (thedata+iset)->kk + 1 ) / 2);
    if ( (thedata+iset)->nn > 0 ) {
      for ( ii = 1 ; ii <= rows ; ii++ )
        for ( jj = 1 ; jj <= cols ; jj++ ) {
           pij = *(*( (thedata+iset)->freq_ests+0)+jj); 
           pk = *(*( (thedata+iset)->freq_ests+ii)+0);
           lb1 = pij*pk;
           lb2 = (1-pij)*(1-pk);
           ub1 = pij*(1-pk);
           ub2 = (1-pij)*pk;
           if ( lb2 < lb1 )
             lb1 = lb2;
           if ( ub2 < ub1 )
             ub1 = ub2;
           if ( ag == 1 ) {
           /*we will need pii and pstar = 1 - sum pij */
             colctr = colindex(jj,jj,(thedata+iset)->kk);
             pii = (double) (thedata+iset)->counts[0][colctr] /  (double)  (thedata+iset)->nn ;
             pstar = 1 - (thedata+iset)->pstar_ests[jj];
             lb3 = 0.5 * ( pii*pk + pstar * (1.0-pk) );
             ub3 = 0.5 * ( pii*(1.0-pk) + pstar * pk );
             if ( lb3 < lb1 )
               lb1 = lb3;
             if ( ub3 < ub1 )
               ub1 = ub3;
           }
           lb1 = -lb1;
           pijk = (thedata+iset)->freq_ests[ii][jj];
           dijk = (thedata+iset)->diseq_ests[ii][jj];
           if ( dijk < 0.0  && lb1 != 0.0 )
             dijk = -dijk / lb1;
           if ( dijk > 0.0 && ub1 != 0.0 )
             dijk = dijk / ub1;
           (thedata+iset)->ndiseq_ests[ii][jj] = dijk ;
            
      }
      colctr = 0;
      if ( ag == 0 )
        for ( ii = 1 ; ii <=  (thedata+iset)->kk ; ii++ )
          for ( jj = ii ; jj <= (thedata+iset)->kk ; jj++ ) {
            colctr = colctr+1;
            if ( ii != jj )
              *(*( (thedata+iset)->diseq_ests + 0)+colctr) = (double) *((thedata+iset)->acounts+ii) * twonn * (double)*((thedata+iset)->acounts+jj)*twonn  - *(*( (thedata+iset)->freq_ests+0)+colctr)/ 2.0 ;
            else
              *(*( (thedata+iset)->diseq_ests + 0)+colctr) = - (double)*((thedata+iset)->acounts+ii) * twonn * (double)*((thedata+iset)->acounts+jj)*twonn  + *(*( (thedata+iset)->freq_ests+0)+colctr);
          }
    }
    else {
      fileptr = fileopen(errorfile,"a");
      fprintf(fileptr,"\n\ncalc_ndiseq: Sample size for dataset %3d is < 0...\n",iset);
      fileclose(errorfile,fileptr);
    }
  } 
}


/*
  Calculate the normalized asymptotic disequilibria.  If ag == 0, then do it for the 
  genotypic values (Dijk), otherwise do it for the allelic values (Dik).
  
  This requires that the genotypic work was done first.
*/
void  calc_avars(thedata,numsets,errorfile,ag)
  adataset *thedata;
  int numsets;
  char *errorfile;
  int ag;
{
  int iset,ii,jj,rows,cols,colctr;
  double m1,m2,d1,nn,piik,pii,resid,residt;
  FILE *fileptr;
  for ( iset = 1 ; iset <= numsets ; iset++ ) {
    rows = (thedata+iset)->mm;
    if ( ag == 0 )
      cols = (thedata+iset)->kk * ( (thedata+iset)->kk + 1 ) / 2 ;
    else 
      cols = (thedata+iset)->kk;
    if ( (thedata+iset)->var0_ests == NULL ) {
      (thedata+iset)->var0_ests = dmatrix(0,rows,0,(thedata+iset)->kk * ( (thedata+iset)->kk + 1 ) / 2);
      (thedata+iset)->var1_ests = dmatrix(0,rows,0,(thedata+iset)->kk * ( (thedata+iset)->kk + 1 ) / 2);
    }
    for ( ii = 0 ; ii <= rows ; ii++ ) {
      *(*( (thedata+iset)->var0_ests + ii)+0) = 0.0;
      *(*( (thedata+iset)->var1_ests + ii)+0) = 0.0;
    }
    for ( ii = 0 ; ii <= cols ; ii++ ) {
      *(*( (thedata+iset)->var0_ests + 0)+ii) = 0.0;
      *(*( (thedata+iset)->var1_ests + 0)+ii) = 0.0;
    }
    if ( (thedata+iset)->nn > 0 ) {
      if ( ag == 0 )
          nn = (double) (thedata+iset)->nn;
      else 
          nn = 2.0 * (double) (thedata+iset)->nn;      
      for ( ii = 1 ; ii <= rows ; ii++ )
        for ( jj = 1 ; jj <= cols ; jj++ ) {
          m1 =  *(*( (thedata+iset)->freq_ests+0)+jj);
          m2 =  *(*( (thedata+iset)->freq_ests+ii)+0);
          d1 =  *(*( (thedata+iset)->diseq_ests+ii)+jj);
          if ( ag == 0 ) {
             resid  =  0;   
             residt =  d1*(1.0-2.0*m1)*(1.0-2.0*m2) - d1*d1;
          }
          else {
             colctr = colindex(jj,jj,(thedata+iset)->kk);
             pii = (double) (thedata+iset)->counts[0][colctr] /  (double)  (thedata+iset)->nn ;
             piik = (double) (thedata+iset)->counts[ii][colctr] /  (double)  (thedata+iset)->nn ;
             if ( (thedata+iset)->pr_ests[ii][colctr] < 0.05 )   
               resid  =  m2*(1-m2)*(pii-m1*m1) ; 
             else       
               resid  =  m2*(1-m2)*(pii-m1*m1) + (1-2.0*m2)*(piik-pii*m2); 
             residt =  d1*(1.0-2.0*m1)*(1.0-2.0*m2) - d1*d1;
          }
          *(*( (thedata+iset)->var0_ests + ii)+jj) =  (m1 * (1.0-m1) * m2*(1.0-m2) + resid)/ nn ;
          *(*( (thedata+iset)->var1_ests + ii)+jj) =  (m1 * (1.0-m1) * m2*(1.0-m2) + resid + residt) /  nn ; 
        }
    }
    else {
      fileptr = fileopen(errorfile,"a");
      fprintf(fileptr,"\n\ncalc_avars: Sample size for dataset %3d is < 0...\n",iset);
      fileclose(errorfile,fileptr);
    }
  } 
}

/*
  Calculate the  asymptotic test statistics for the disequilibria.  If ag == 0, then do it for the 
  genotypic values (Dijk), otherwise do it for the allelic values (Dik).
*/
void  calc_astat(thedata,numsets,errorfile,ag)
  adataset *thedata;
  int numsets;
  char *errorfile;
  int ag;
{
  int iset,ii,jj,rows,cols;
  double m1,m2,d1,nn;
  FILE *fileptr;
  for ( iset = 1 ; iset <= numsets ; iset++ ) {
    rows = (thedata+iset)->mm;
    if ( ag == 0 )
      cols = (thedata+iset)->kk * ( (thedata+iset)->kk + 1 ) / 2 ;
    else 
      cols = (thedata+iset)->kk;
    if ( (thedata+iset)->stat_ests == NULL ) 
      (thedata+iset)->stat_ests = dmatrix(0,rows,0,(thedata+iset)->kk * ( (thedata+iset)->kk + 1 ) / 2);
    for ( ii = 0 ; ii <= rows ; ii++ ) 
      *(*( (thedata+iset)->stat_ests + ii)+0) = 0.0;
    for ( ii = 0 ; ii <= cols ; ii++ ) 
      *(*( (thedata+iset)->stat_ests + 0)+ii) = 0.0;
    if ( (thedata+iset)->nn > 0 )
      for ( ii = 1 ; ii <= rows ; ii++ )
        for ( jj = 1 ; jj <= cols ; jj++ ) {
          d1 =  *(*( (thedata+iset)->diseq_ests+ii)+jj);
          *(*( (thedata+iset)->stat_ests + ii)+jj) =  d1*d1 / *(*((thedata+iset)->var0_ests + ii)+jj); 
        }
    else {
      fileptr = fileopen(errorfile,"a");
      fprintf(fileptr,"\n\ncalc_astat: Sample size for dataset %3d is < 0...\n",iset);
      fileclose(errorfile,fileptr);
    }
  } 
}


/*
  Calculate the  asymptotic test statistics p values for the disequilibria. 
  If ag == 0, then do it for the 
  genotypic values (Dijk), otherwise do it for the allelic values (Dik).
*/

void  calc_aprs(thedata,numsets,errorfile,ag)
  adataset *thedata;
  int numsets;
  char *errorfile;
{
  int iset,ii,jj,rows,cols;
  double m1,m2,d1,nn;
  FILE *fileptr;
  for ( iset = 1 ; iset <= numsets ; iset++ ) {
    rows = (thedata+iset)->mm;
    if ( ag == 0 )
      cols = (thedata+iset)->kk * ( (thedata+iset)->kk + 1 ) / 2 ;
    else 
      cols = (thedata+iset)->kk;
    if ( (thedata+iset)->pr_ests == NULL ) 
      (thedata+iset)->pr_ests = dmatrix(0,rows,0,(thedata+iset)->kk * ( (thedata+iset)->kk + 1 ) / 2);
    for ( ii = 0 ; ii <= rows ; ii++ ) 
      *(*( (thedata+iset)->pr_ests + ii)+0) = 0.0;
    for ( ii = 0 ; ii <= cols ; ii++ ) 
      *(*( (thedata+iset)->pr_ests + 0)+ii) = 0.0;
    if ( (thedata+iset)->nn > 0 )
      for ( ii = 1 ; ii <= rows ; ii++ )
        for ( jj = 1 ; jj <= cols ; jj++ ) {
          *(*( (thedata+iset)->pr_ests + ii)+jj) = chiprob(1,*(*((thedata+iset)->stat_ests + ii)+jj) ); 
        }
    else {
      fileptr = fileopen(errorfile,"a");
      fprintf(fileptr,"\n\ncalc_astat: Sample size for dataset %3d is < 0...\n",iset);
      fileclose(errorfile,fileptr);
    }
  } 
}

/*
  Write  asymptotic results. If ag == 0, then do it for the 
  genotypic values (Dijk), otherwise do it for the allelic values (Dik).
*/

void  write_asymp_results(thedata,numsets,outfile,errorfile,ag)
  adataset *thedata;
  int numsets;
  char *outfile,*errorfile;
  int ag;
{
  int iset,ii,jj,rows,cols;
  FILE *fileptr;
  fileptr = fileopen(outfile,"a");
  for ( ii = 1 ; ii <= numsets ; ii++ )  {
    for ( jj = 1 ; jj <= 80 ; jj++ ) 
      fprintf(fileptr,"=");
    fprintf(fileptr,"\n");
    fprintf(fileptr,"\n\n\tAnalysis of nuclear system %s by cytotype %s\n\n",(thedata+ii)->nnames,(thedata+ii)->cnames);


    write_counts((thedata+ii),fileptr,ag);
    fprintf(fileptr,"\t\tFrequencies...\n");   
   
    write_thematrix(fileptr,(thedata+ii)->freq_ests,(thedata+ii)->mm,(thedata+ii)->kk,ag);
    fprintf(fileptr,"\t\tDisequilibria...\n");
   
    write_thematrix(fileptr,(thedata+ii)->diseq_ests,(thedata+ii)->mm,(thedata+ii)->kk,ag);
    fprintf(fileptr,"\t\tNormalized Disequilibria for Dijk (not Dij)...\n");
   
    write_thematrix(fileptr,(thedata+ii)->ndiseq_ests,(thedata+ii)->mm,(thedata+ii)->kk,ag);
    fprintf(fileptr,"\t\tVariances under H0...\n");
   
    write_thematrix(fileptr,(thedata+ii)->var0_ests,(thedata+ii)->mm,(thedata+ii)->kk,ag);
    fprintf(fileptr,"\t\tVariances under H1...\n");
   
    write_thematrix(fileptr,(thedata+ii)->var1_ests,(thedata+ii)->mm,(thedata+ii)->kk,ag);


    calc_stdev((thedata+ii)->var0_ests,(thedata+ii)->mm,(thedata+ii)->kk,ag);
    calc_stdev((thedata+ii)->var1_ests,(thedata+ii)->mm,(thedata+ii)->kk,ag);
    fprintf(fileptr,"\t\tSample standard deviations under H0...\n");
   
    write_thematrix(fileptr,(thedata+ii)->var0_ests,(thedata+ii)->mm,(thedata+ii)->kk,ag);
    fprintf(fileptr,"\t\tSample standard deviations under H1...\n");
   
    write_thematrix(fileptr,(thedata+ii)->var1_ests,(thedata+ii)->mm,(thedata+ii)->kk,ag);

    fprintf(fileptr,"\tFor alpha = 0.05\n");

    write_samplesizes(fileptr,(thedata+ii),ag, (double) Zalpha);
    fprintf(fileptr,"\tFor alpha = 0.01\n");

    write_samplesizes(fileptr,(thedata+ii),ag, (double) Zalpha2);

    fprintf(fileptr,"\t\tTest statistics...\n");
   
    write_thematrix(fileptr,(thedata+ii)->stat_ests,(thedata+ii)->mm,(thedata+ii)->kk,ag);
    fprintf(fileptr,"\t\tProbabilities of test statistics under X(1)...\n");
   
    write_thematrix(fileptr,(thedata+ii)->pr_ests,(thedata+ii)->mm,(thedata+ii)->kk,ag);
    for ( jj = 1 ; jj <= 80 ; jj++ ) 
      fprintf(fileptr,"=");
    fprintf(fileptr,"\n");
  }
  fileclose(outfile,fileptr);


}

void  write_asymp_results_tex(thedata,numsets,outfile,errorfile,ag)
  adataset *thedata;
  int numsets;
  char *outfile,*errorfile;
  int ag;
{
  int iset,ii,jj,rows,cols;
  FILE *fileptr;
  fileptr = fileopen(outfile,"a");
  for ( ii = 1 ; ii <= numsets ; ii++ ) {
    fprintf(fileptr,"\n\\section{Analysis of nuclear system %s by cytotype %s}\n\n",(thedata+ii)->nnames,(thedata+ii)->cnames);
    write_counts_tex((thedata+ii),fileptr,ag);
    fprintf(fileptr,"\\subsection{Frequencies}\n");   
    write_thematrix_tex(fileptr,(thedata+ii)->freq_ests,(thedata+ii)->mm,(thedata+ii)->kk,ag);
    fprintf(fileptr,"\\subsection{Disequilibria}\n");
    write_thematrix_tex(fileptr,(thedata+ii)->diseq_ests,(thedata+ii)->mm,(thedata+ii)->kk,ag);
    fprintf(fileptr,"\\subsection{Normalized Disequilibria}\n");
    write_thematrix_tex(fileptr,(thedata+ii)->ndiseq_ests,(thedata+ii)->mm,(thedata+ii)->kk,ag);
    fprintf(fileptr,"\\subsection{Variances under $H_0$}\n");
    write_thematrix_tex(fileptr,(thedata+ii)->var0_ests,(thedata+ii)->mm,(thedata+ii)->kk,ag);
    fprintf(fileptr,"\\subsection{Variances under $H_1$}\n");
    write_thematrix_tex(fileptr,(thedata+ii)->var1_ests,(thedata+ii)->mm,(thedata+ii)->kk,ag);
    calc_stdev((thedata+ii)->var0_ests,(thedata+ii)->mm,(thedata+ii)->kk,ag);
    calc_stdev((thedata+ii)->var1_ests,(thedata+ii)->mm,(thedata+ii)->kk,ag);
    fprintf(fileptr,"\\subsection{Sample standard deviations under $H_0$}\n");
    write_thematrix_tex(fileptr,(thedata+ii)->var0_ests,(thedata+ii)->mm,(thedata+ii)->kk,ag);
    fprintf(fileptr,"\\subsection{Sample standard deviations under $H_1$}\n");
    write_thematrix_tex(fileptr,(thedata+ii)->var1_ests,(thedata+ii)->mm,(thedata+ii)->kk,ag);
    write_samplesizes_tex(fileptr,(thedata+ii),ag);
    fprintf(fileptr,"\\subsection{Asymptotic Test Statistics}\n");
    write_thematrix_tex(fileptr,(thedata+ii)->stat_ests,(thedata+ii)->mm,(thedata+ii)->kk,ag);
    fprintf(fileptr,"\\subsection{Probabilities of test statistics under $\\chi^2_1$}\n");
    write_thematrix_tex(fileptr,(thedata+ii)->pr_ests,(thedata+ii)->mm,(thedata+ii)->kk,ag);
    fprintf(fileptr,"\n");
  }
  fileclose(outfile,fileptr);
}

void write_samplesizes(fileptr,thedata,ag,zalph)
FILE *fileptr;
adataset *thedata;
int ag;
float zalph;
{
  int rows,cols,iii,ii,jj,jjj,nn;
  double d0,d1,dalt,ss;
  if ( *(*(thedata->counts+0)+0) > 0 )
    ss = sqrt( (double) *(*(thedata->counts+0)+0) );
  else
    return;
  fprintf(fileptr,"\n");
  for ( ii = 1 ; ii <= 80 ; ii++ ) 
    fprintf(fileptr,"*");
  fprintf(fileptr,"\n");
  cols = thedata->mm;
  if ( ag == 0 )
    rows = thedata->kk*(thedata->kk+1)/2;
  else
    rows = thedata->kk ;
  fprintf(fileptr,"\n\n\tSample sizes to detect specified disequilibria...Power = 90%%");
  if ( ag == 0 )
    fprintf(fileptr,"\n\t\tNuclear Genotypes (rows) by Cytotypes (columns)\n   |");
  else
    fprintf(fileptr,"\n\t\tNuclear Alleles (rows) by Cytotypes (columns)\n   |");
  for ( ii = 1 ; ii <= cols ; ii++ )
    fprintf(fileptr," %4d",ii);
  fprintf(fileptr," | Margin \n");
  for ( ii = 1 ; ii <= 5*(cols+1)+10 ; ii++ )
    fprintf(fileptr,"-");
  fprintf(fileptr,"\n");
  if ( ag == 0 ) {
	  iii = 0;
	  for ( ii = 1 ; ii <= thedata->kk ; ii++ )
	    for ( jj = ii ; jj <= thedata->kk ; jj++ ) {
	      iii = iii+1;
	      fprintf(fileptr," %2d,%-2d |",ii,jj);
	      for ( jjj = 1 ; jjj <= cols ; jjj++ ) {
	        d1 = ss * *(*(thedata->var1_ests+jjj)+iii);
	        d0 = ss * *(*(thedata->var0_ests+jjj)+iii);
	        dalt =  *(*(thedata->diseq_ests+jjj)+iii);
	        nn = newsamplesize(d0,d1,dalt,zalph);
	        if ( nn != 9999 )
	          fprintf(fileptr," %4d",nn );
	        else
	          fprintf(fileptr,"     ");
	      }
	      nn = 0;
	      fprintf(fileptr," | %4d\n", nn);
	    }
  }
  else {
    for ( jj = 1 ; jj <= thedata->kk ; jj++ ) {
      fprintf(fileptr," %2d    |", jj);
      for ( jjj = 1 ; jjj <= cols ; jjj++ ) {
        d1 = ss * *(*(thedata->var1_ests+jjj)+jj);
        d0 = ss * *(*(thedata->var0_ests+jjj)+jj);
        dalt =  *(*(thedata->diseq_ests+jjj)+jj);
        nn = newsamplesize(d0,d1,dalt,zalph);
        if ( nn != 9999 )
          fprintf(fileptr," %4d",nn );
        else
          fprintf(fileptr,"     ");
      }
      nn = 0;
      fprintf(fileptr," | %4d\n", nn);
    }
  }
  for ( ii = 1 ; ii <= 5*(cols+1)+10 ; ii++ )
    fprintf(fileptr,"-");
  fprintf(fileptr,"\n");
  for ( ii = 1 ; ii <= 80 ; ii++ ) 
    fprintf(fileptr,"*");
  fprintf(fileptr,"\n");
  d1 = 0.0;

  fprintf(fileptr,"\n\n\tSample sizes to detect specified disequilibria...Power = 50%%");
  if ( ag == 0 )
    fprintf(fileptr,"\n\t\tNuclear Genotypes (rows) by Cytotypes (columns)\n   |");
  else
    fprintf(fileptr,"\n\t\tNuclear Alleles (rows) by Cytotypes (columns)\n   |");
  for ( ii = 1 ; ii <= cols ; ii++ )
    fprintf(fileptr," %4d",ii);
  fprintf(fileptr," | Margin \n");
  for ( ii = 1 ; ii <= 5*(cols+1)+10 ; ii++ )
    fprintf(fileptr,"-");
  fprintf(fileptr,"\n");
  if ( ag == 0 ) {
	  iii = 0;
	  for ( ii = 1 ; ii <= thedata->kk ; ii++ )
	    for ( jj = ii ; jj <= thedata->kk ; jj++ ) {
	      iii = iii+1;
	      fprintf(fileptr," %2d,%-2d |",ii,jj);
	      for ( jjj = 1 ; jjj <= cols ; jjj++ ) {
	        d0 = ss * *(*(thedata->var0_ests+jjj)+iii);
	        dalt =  *(*(thedata->diseq_ests+jjj)+iii);
	        nn = newsamplesize(d0,d1,dalt,zalph);
	        if ( nn != 9999 )
	          fprintf(fileptr," %4d",nn );
	        else 
	          fprintf(fileptr,"     ");
	      }
	      nn = 0;
	      fprintf(fileptr," | %4d\n", nn);
	    }
  }
  else {
    for ( jj = 1 ; jj <= thedata->kk ; jj++ ) {
      fprintf(fileptr," %2d    |", jj);
      for ( jjj = 1 ; jjj <= cols ; jjj++ ) {
        d0 = ss * *(*(thedata->var0_ests+jjj)+jj);
        dalt =  *(*(thedata->diseq_ests+jjj)+jj);
        nn = newsamplesize(d0,d1,dalt,zalph);
        if ( nn != 9999 )
          fprintf(fileptr," %4d",nn );
        else 
          fprintf(fileptr,"     ");
      }
      nn = 0;
      fprintf(fileptr," | %4d\n", nn);
    }
  }
  for ( ii = 1 ; ii <= 5*(cols+1)+10 ; ii++ )
    fprintf(fileptr,"-");
  fprintf(fileptr,"\n");
  for ( ii = 1 ; ii <= 80 ; ii++ ) 
    fprintf(fileptr,"*");
  fprintf(fileptr,"\n");
}

void write_samplesizes_tex(fileptr,thedata,ag)
FILE *fileptr;
adataset *thedata;
{
  int rows,cols,iii,ii,jj,jjj,nn;
  double d0,d1,dalt,ss;
  if ( *(*(thedata->counts+0)+0) > 0 )
    ss = sqrt( (double) *(*(thedata->counts+0)+0) );
  else
    return;
  fprintf(fileptr,"\n");
  cols = thedata->mm;
  if ( ag == 0 )
    rows = thedata->kk*(thedata->kk+1)/2;
  else
    rows = thedata->kk ;
  fprintf(fileptr,"\\subsection{Minimum sample sizes}\n");
  fprintf(fileptr,"\nThese are the sample sizes to detect the estimated disequilibria in this data set with power  90\\%%.");
  if ( ag == 0 )
    fprintf(fileptr,"\nRows are nuclear genotypes  and columns are cytotypes.\n");
  else
    fprintf(fileptr,"\nRows are nuclear alleles  and columns are cytotypes.\n");
  fprintf(fileptr,"\n\\begin{center}\n\\begin{tabular}{");
  for ( ii = 0 ; ii <= cols+1 ; ii++ )
    fprintf(fileptr,"c");
  fprintf(fileptr,"}\\hline\\hline\n");
  for ( ii = 1 ; ii <= cols ; ii++ )
    fprintf(fileptr,"      & %4d",ii);
  fprintf(fileptr," & Margin \\\\ \\hline\n");
  if ( ag == 0 ) {
  iii = 0;
  for ( ii = 1 ; ii <= thedata->kk ; ii++ )
    for ( jj = ii ; jj <= thedata->kk ; jj++ ) {
      iii = iii+1;
      fprintf(fileptr," %2d,%-2d  ",ii,jj);
      for ( jjj = 1 ; jjj <= cols ; jjj++ ) {
        d1 = ss * *(*(thedata->var1_ests+jjj)+iii);
        d0 = ss * *(*(thedata->var0_ests+jjj)+iii);
        dalt =  *(*(thedata->diseq_ests+jjj)+iii);
        nn = samplesize(d0,d1,dalt);
        if ( nn != 9999 )
          fprintf(fileptr," & %4d",nn );
        else
          fprintf(fileptr," &    ");
      }
      nn = 0;
      fprintf(fileptr," & %4d \\\\ \n", nn);
    }
  }
  else {
    for ( jj = 1 ; jj <= thedata->kk ; jj++ ) {
      fprintf(fileptr," %2d   "  ,jj);
      for ( jjj = 1 ; jjj <= cols ; jjj++ ) {
        d1 = ss * *(*(thedata->var1_ests+jjj)+jj);
        d0 = ss * *(*(thedata->var0_ests+jjj)+jj);
        dalt =  *(*(thedata->diseq_ests+jjj)+jj);
        nn = samplesize(d0,d1,dalt);
        if ( nn != 9999 )
          fprintf(fileptr," & %4d",nn );
        else
          fprintf(fileptr," &    ");
      }
      nn = 0;
      fprintf(fileptr," & %4d \\\\ \n", nn);
    }
  }
  fprintf(fileptr,"\\hline\\hline\n\\end{tabular}\n\\end{center}\n\n");

  d1 = 0.0;
  fprintf(fileptr,"\nThese are the sample sizes to detect the estimated disequilibria in this data set with power  50\\%%.");
  if ( ag == 0 )
    fprintf(fileptr,"\nRows are nuclear genotypes  and columns are cytotypes.\n");
  else
    fprintf(fileptr,"\nRows are nuclear alleles  and columns are cytotypes.\n");
  fprintf(fileptr,"\n\\begin{center}\n\\begin{tabular}{");
  for ( ii = 0 ; ii <= cols+1 ; ii++ )
    fprintf(fileptr,"c");
  fprintf(fileptr,"}\\hline\\hline\n");
  for ( ii = 1 ; ii <= cols ; ii++ )
    fprintf(fileptr," & %4d",ii);
  fprintf(fileptr," & Margin \\\\ \\hline\n");

  if ( ag == 0 ) {
  iii = 0;
  for ( ii = 1 ; ii <= thedata->kk ; ii++ )
    for ( jj = ii ; jj <= thedata->kk ; jj++ ) {
      iii = iii+1;
      fprintf(fileptr," %2d,%-2d  ",ii,jj);
      for ( jjj = 1 ; jjj <= cols ; jjj++ ) {
        d0 = ss * *(*(thedata->var0_ests+jjj)+iii);
        dalt =  *(*(thedata->diseq_ests+jjj)+iii);
        nn = samplesize(d0,d1,dalt);
        if ( nn != 9999 )
          fprintf(fileptr," & %4d",nn );
        else
          fprintf(fileptr," &    ");
      }
      nn = 0;
      fprintf(fileptr," & %4d \\\\ \n", nn);
    }
  }
  else {
    for ( jj = 1 ; jj <= thedata->kk ; jj++ ) {
      fprintf(fileptr," %2d   "  ,jj);
      for ( jjj = 1 ; jjj <= cols ; jjj++ ) {
        d0 = ss * *(*(thedata->var0_ests+jjj)+jj);
        dalt =  *(*(thedata->diseq_ests+jjj)+jj);
        nn = samplesize(d0,d1,dalt);
        if ( nn != 9999 )
          fprintf(fileptr," & %4d",nn );
        else
          fprintf(fileptr," &    ");
      }
      nn = 0;
      fprintf(fileptr," & %4d \\\\ \n", nn);
    }
  }
  fprintf(fileptr,"\\hline\\hline\n\\end{tabular}\n\\end{center}\n\n");


}


/* Use Dimitry's subroutines for the exact test for a 2 x 2 table */
void  calc_Exact_prs(thedata,numsets,outfile,errorfile,tex_flag,ag)  
  adataset *thedata;
  int numsets;
  char *outfile,*errorfile; 
  int tex_flag,ag;
{
  int table[2][2],ii,jj,iset,rows,cols,**mptr;
  float bothtails, closesttail;
  FILE *fileptr;
  double mi,mj,nn;
  
  for ( iset = 1 ; iset <= numsets ; iset++ ) {
    rows = (thedata+iset)->mm;
    if ( ag == 0 ) {
      mptr = (thedata+iset)->counts ;
      cols = (thedata+iset)->kk * ( (thedata+iset)->kk + 1 ) / 2;
    }
    else {
      mptr = (thedata+iset)->kmcounts ;
      cols = (thedata+iset)->kk ;
    }
    /* fill in the table */
    for ( ii = 1 ; ii <= rows ; ii++ ) 
      for ( jj = 1 ; jj <= cols ; jj++ ) {
        table[0][0] = *(*( mptr+ii)+jj);
        table[0][1] = *(*( mptr+ii)+0) - *(*(mptr+ii)+jj);
        table[1][0] =  *(*( mptr+0)+jj) - *(*(mptr+ii)+jj);  
        table[1][1] = *(*( mptr+0)+0) - table[0][0]-table[0][1]- table[1][0];
        mi = *(*( (thedata+iset)->freq_ests+ii)+0);
        mj = *(*( (thedata+iset)->freq_ests+0)+jj);
/*        if ( !( mi > 0.2 && mj > 0.2 && mi < 0.8 && mj < 0.8 && *(*((thedata+iset)->counts+0)+0) > 100 ) ) {*/          
          if (Exact(table, &bothtails, &closesttail) != -1) {
            *(*( (thedata+iset)->pr_ests+ii)+jj) = bothtails;
          }
          else {
            *(*( (thedata+iset)->pr_ests+ii)+jj) = -1.0;
            fileptr = fileopen(errorfile,"a");
            fprintf(fileptr,"\n\ncalc_Exact_prs: Problem with dataset %3d...\n",iset);
            fileclose(errorfile,fileptr);
          }
      }

    fileptr = fileopen(outfile,"a");
    fprintf(fileptr,"\n");
    if ( tex_flag == 0 ) {
      fprintf(fileptr,"\n\n\tAnalysis of nuclear system %s by cytotype%s\n\n",(thedata+iset)->nnames,(thedata+iset)->cnames);
      write_counts((thedata+iset),fileptr,ag);
      fprintf(fileptr,"\t\tProbabilities from Fisher's exact test...\n");   
      write_thematrix(fileptr,(thedata+iset)->pr_ests,(thedata+iset)->mm,(thedata+iset)->kk,ag);
      for ( jj = 1 ; jj <= 80 ; jj++ ) 
        fprintf(fileptr,"=");
      fprintf(fileptr,"\n");
    }
    else {
      fprintf(fileptr,"\n\\subsection{Fisher's exact test for nuclear system %s by cytotype %s}\n",(thedata+iset)->nnames,(thedata+iset)->cnames);
      write_counts_tex((thedata+iset),fileptr,ag);
      write_thematrix_tex(fileptr,(thedata+iset)->pr_ests,(thedata+iset)->mm,(thedata+iset)->kk,ag);
      fprintf(fileptr,"\n");
    }
    fileclose(outfile,fileptr);
  }  
}


void  calc_Markov_prs(thedata,numsets,outfile,errorfile,reps,bb,cc,verbosity,tex_flag,ag) 
  adataset *thedata;
  int numsets;
  char *outfile,*errorfile;
  long reps,bb,cc;
  int verbosity,tex_flag,ag;
{
  int  ii,jj,iset,rows,cols, **counts,**mptr;
  long cntr,prelim,lii,lxcntr;
  FILE *fileptr;
  double rho,tot,var,tot2,hh;
  ldiv_t xx; 
  prelim = 1000L;
  var = 0.0;
  for ( iset = 1 ; iset <= numsets ; iset++ ) {
    if ( verbosity == 1 )
          printf("\n Now doing Markov Chain test for Data Set %d\n\twith %ld reps in %ld batches of %ld observations\n",iset,reps,bb,cc);
    rows = (thedata+iset)->mm;
    if ( ag == 0 ) {
      mptr = (thedata+iset)->counts;
      cols = (thedata+iset)->kk * ((thedata+iset)->kk+1) / 2;
    }
    else {
      mptr = (thedata+iset)->kmcounts;
      cols = (thedata+iset)->kk  ;
    }
    counts = imatrix(1,rows,1,cols);
    for ( ii = 1 ; ii <= rows ; ii++ )
      for ( jj = 1 ; jj <= cols ; jj++ )
        *(*(counts+ii)+jj) = *(*( mptr+ii)+jj);
    rho = 0.0;
    lxcntr = cntr = 0L;
    tot = tot2 = 0.0;
    for ( lii = 1 ; lii <= prelim ; lii++ ) 
      pickNswitch(counts,rows,cols,&rho);
    check_counts(counts,(thedata+iset)->counts,rows,cols,errorfile);
    for ( lii = 1 ; lii <= reps ; lii++ ) {
      pickNswitch(counts,rows,cols,&rho);
      if ( rho <= 0.0 ) {
        cntr += 1L;
        lxcntr += 1L;
      }
      xx = ldiv(lii,cc);
      if ( xx.rem == 0 ) {
        hh = (double) lxcntr/ (double) cc ;
        tot += hh;
        tot2 += hh*hh;
        lxcntr = 0L;
        if ( verbosity == 1 ) {
          printf(".");
          fflush(stdout);
        }
      }      
    }
    if ( verbosity == 1 )
          printf("\n" );
    var = (tot2 -  (tot*tot / ((double) bb)) )/( (double) (bb*(bb-1L)) )  ;
    if ( var > 0.0 )
      var = sqrt(var);
    if ( verbosity == 1 )
          printf("\n " );
    check_counts(counts,mptr,rows,cols,errorfile);
    (thedata+iset)->Markov = (double) cntr / (double) reps ;
    fileptr = fileopen(outfile,"a");
    if ( tex_flag == 0 ) {
	    fprintf(fileptr,"\n\n");
	    for ( jj = 1 ; jj <= 80 ; jj++ ) 
	      fprintf(fileptr,"=");
	    fprintf(fileptr,"\n");
	    if ( ag == 0 )
	      fprintf(fileptr,"\n\n\tMarkov Chain Monte Carlo analysis of nuclear system %s by cytotype %s : Genotypic Disequilibria\n\n",(thedata+iset)->nnames,(thedata+iset)->cnames);
	    else
	      fprintf(fileptr,"\n\n\tMarkov Chain Monte Carlo analysis of nuclear system %s by cytotype %s : Allelic Disequilibria\n\n",(thedata+iset)->nnames,(thedata+iset)->cnames);
	    write_counts((thedata+iset),fileptr,ag);
	    if ( var > 0.0 )
	      fprintf(fileptr,"\n\nMarkov test pr =  %12.8f with sample standard deviation %12.8f",(thedata+iset)->Markov,var);
	    else
	      fprintf(fileptr,"\n\nMarkov test pr =  %12.8f with sample variance %12.8f",(thedata+iset)->Markov,var);
	    fprintf(fileptr,"\nFor %ld / %ld\n\n",cntr,reps);
	    for ( jj = 1 ; jj <= 80 ; jj++ ) 
	      fprintf(fileptr,"=");
	    fprintf(fileptr,"\n");
	}
	else {
	  if ( ag == 0 )
	    fprintf(fileptr,"\n\n\\subsection{Markov Chain Monte Carlo Genotypic Results}  ");
	  else
	    fprintf(fileptr,"\n\n\\subsection{Markov Chain Monte Carlo Allelic Results}  ");
	    fprintf(fileptr," This is the analysis of nuclear system %s by cytotype %s.",(thedata+iset)->nnames,(thedata+iset)->cnames);
	    write_counts_tex((thedata+iset),fileptr,ag);
	    if ( var > 0.0 )
	      fprintf(fileptr," The Markov test probability was  %f\nand had a sample standard deviation of %f.",(thedata+iset)->Markov,var);
	    else
	      fprintf(fileptr,"\n\nThe Markov test probability was  %f\nand had a sample variance of %f.",(thedata+iset)->Markov,var);
	    fprintf(fileptr,"\nThere were %ld observations as rare as the data out of %ld repetitions.  ",cntr,reps);
	    fprintf(fileptr,"\nThe standard error (sampling variance) was based on %ld batches of %ld observations.",bb,cc);

	}
    fileclose(outfile,fileptr);
    free_imatrix(counts,1,rows,1,cols);
  }
}

void check_counts(counts1,counts,rows,cols,errorfile)
  int **counts1,**counts,rows,cols;
  char *errorfile;
{
  int ii,jj,rt,ct;
  FILE *fileptr;
    
  fileptr = fileopen(errorfile,"a");
  fprintf(fileptr,"\n\n\tChecking counts after Markov...\n");

  for ( ii = 1 ; ii <= rows ; ii++ )
    for ( jj = 1 ; jj <= cols ; jj++ )
      if ( *(*(counts1+ii)+jj) < 0 )
        fprintf(fileptr,"\ncount %3d %3d <0",ii,jj);
  for ( ii = 1 ; ii <= rows ; ii++ ) {
    rt = 0;
    for ( jj = 1 ; jj <= cols ; jj++ )
      rt = rt + *(*(counts1+ii)+jj);
      if ( rt != *(*(counts+ii)+0) )
        fprintf(fileptr,"\nrow margin %3d changed",ii);
  }
  for ( jj = 1 ; jj <= cols ; jj++ ) {
    ct = 0;
    for ( ii = 1 ; ii <= rows ; ii++ )
      ct = ct + *(*(counts1+ii)+jj);
      if ( ct != *(*(counts+0)+jj) )
        fprintf(fileptr,"\ncol margin %3d changed",jj);
  }
  fprintf(fileptr,"\n");
  fileclose(errorfile,fileptr);
}

void pickNswitch(counts,rows,cols,rho)
  int **counts,rows,cols;
  double *rho;
{
  int   i1,i2,j1,j2,R,D;
  double pR,pD,qR,qD,prb, aa,bb,cc,dd;
    pick_ints(&i1,&i2,&j1,&j2,rows,cols);
    R = D = 0;
    pR = pD = qR = qD = 0.0;
    aa = (double)  *(*(counts+i1)+j1);
    bb = (double)  *(*(counts+i1)+j2);
    cc = (double)  *(*(counts+i2)+j1);
    dd = (double)  *(*(counts+i2)+j2);
    if ( *(*(counts+i1)+j2) > 0 &&  *(*(counts+i2)+j1) > 0 ) {
      R = 1;
      pR =   bb*cc/( (aa+1.0) * (dd+1.0) );
    }
    if  ( *(*(counts+i1)+j1) > 0 &&  *(*(counts+i2)+j2) > 0 ) {
      D = 1;
      pD =  aa*dd/( (cc+1.0) * (bb+1.0) );
    }

    if ( R==1 & D==1 ) {
      qR = 0.5*themin(1.0,pR);
      qD = 0.5*themin(1.0,pD);
      prb=ranf(1);
      if ( prb  <= qR ) {
        *rho = *rho + log(pR) ;
        switch_it(counts,i1,i2,j1,j2);
      }
      else if ( prb <= qR+qD ) {
        *rho = *rho + log(pD) ;
        switch_it(counts,i1,i2,j2,j1);
      }
    }
    else if ( R==1  ) {
      qR = 0.5*themin(1.0,pR);
      if ( (prb=ranf(1)) <= qR ) {
        *rho = *rho + log(pR) ;
        switch_it(counts,i1,i2,j1,j2);

      }
    }
    else if ( D==1 ) {
       qD =  0.5*themin(1.0,pD);
       if ( (prb=ranf(1)) <= qD ) {
        *rho = *rho + log(pD) ;
        switch_it(counts,i1,i2,j2,j1);
      }
    }


}



void   pick_ints(i1,i2,j1,j2,rows,cols)
  int *i1,*i2,*j1,*j2,rows,cols;
{
  int tmpi; 
  *i2 = *i1 = 0;
  while ( *i1 > rows || *i1 < 1 )
    *i1 = (int) ceil( (double) rows * ranf(1) ) ;    
  while ( *i2 == *i1 || *i2 < 1 || *i2 > rows ) 
    *i2 = (int) ceil((double) rows * ranf(1)) ;   
 /*    
  if ( *i2 < *i1 ) { / *switch them if i2 is smaller than i1* /
    tmpi = *i1;
    *i1 = *i2;
    *i2 = tmpi;
  }   
*/
  *j2 = *j1 = 0;
  while ( *j1 > cols || *j1 < 1 )
    *j1 = (int) ceil((double) cols * ranf(1))  ;    
  while ( *j2 == *j1 || *j2 < 1 || *j2 > cols ) 
    *j2 = (int) ceil((double) cols * ranf(1)) ;
/*   
  if ( *j2 < *j1 ) { / *switch them if j2 is smaller than j1* /
      tmpi = *j1;
      *j1 = *j2;
      *j2 = tmpi;
  } 
*/
}

void switch_it(counts,i1,i2,j1,j2)
  int **counts,i1,i2,j1,j2;
{
   *(*(counts+i1)+j1) += 1;
   *(*(counts+i2)+j2) += 1;
   *(*(counts+i1)+j2) -= 1;
   *(*(counts+i2)+j1) -= 1;
}

double themin(x,y)
  double x,y;
{
 if ( x < y ) return(x);
 else return(y);
}



void calc_stdev(dblptr,mm,kk,ag)
  double **dblptr;
  int mm,kk,ag;
{
  int ii,jj,rows,cols;
  rows = mm;
  if ( ag == 0 )
    cols = kk*(kk+1)/2;
  else
    cols = kk ;
  for ( ii = 0 ; ii <= rows ; ii++ )
    for ( jj = 0 ; jj <= cols ; jj++ )
      if ( *(*(dblptr+ii)+jj) > 0.0 ) 
        *(*(dblptr+ii)+jj) = sqrt(*(*(dblptr+ii)+jj));
}

void write_thematrix(fileptr,dblptr,mm,kk,ag)
  FILE *fileptr;
  double **dblptr;
  int mm,kk,ag;
{
  int rows,cols,iii,ii,jj,jjj;

  cols = mm;
  if ( ag == 0 )
    rows = kk*(kk+1)/2;
  else
    rows = kk ;

  if ( ag == 0 )
    fprintf(fileptr,"\nNuclear Genotypes (rows) by Cytotypes (columns)\n");
  else
    fprintf(fileptr,"\nNuclear Alleles (rows) by Cytotypes (columns)\n");
  for ( ii = 1 ; ii <= 10*(cols+1)+10 ; ii++ )
    fprintf(fileptr,"=");
  fprintf(fileptr,"\n       |");
  for ( ii = 1 ; ii <= cols ; ii++ )
    fprintf(fileptr,"    %2d    ",ii);
  fprintf(fileptr," | Margin \n");
  for ( ii = 1 ; ii <= 10*(cols+1)+10 ; ii++ )
    fprintf(fileptr,"-");
  fprintf(fileptr,"\n");
  if ( ag == 0 ) {
  iii = 0;
  for ( ii = 1 ; ii <= kk ; ii++ )
    for ( jj = ii ; jj <= kk ; jj++ ) {
      iii = iii+1;
      fprintf(fileptr," %2d,%-2d |",ii,jj);
      for ( jjj = 1 ; jjj <= cols ; jjj++ ) 
        fprintf(fileptr," %9.4f",*(*(dblptr+jjj)+iii) );
      fprintf(fileptr," | %9.4f\n",*(*(dblptr+0)+iii) );
    }
  }
  else { 
    for ( jj = 1 ; jj <= kk ; jj++ ) {
      fprintf(fileptr," %2d    |", jj);
      for ( jjj = 1 ; jjj <= cols ; jjj++ ) 
        fprintf(fileptr," %9.4f",*(*(dblptr+jjj)+jj) );
      fprintf(fileptr," | %9.4f\n",*(*(dblptr+0)+jj) );
    }
  }
  for ( ii = 1 ; ii <= 10*(cols+1)+10 ; ii++ )
    fprintf(fileptr,"-");
  fprintf(fileptr,"\n");
  fprintf(fileptr,"Margins|");
  for ( ii = 1 ; ii <= cols ; ii++ )
    fprintf(fileptr," %9.4f",*(*(dblptr+ii)+0));
  fprintf(fileptr," | %9.4f\n\n",*(*(dblptr+0)+0));

  for ( ii = 1 ; ii <= 10*(cols+1)+10 ; ii++ )
    fprintf(fileptr,"=");
  fprintf(fileptr,"\n\n");


}

void write_thematrix_tex(fileptr,dblptr,mm,kk,ag)
  FILE *fileptr;
  double **dblptr;
  int mm,kk;
{
  int rows,cols,iii,ii,jj,jjj;

  cols = mm;
  if ( ag == 0 )
    rows = kk*(kk+1)/2;
  else
    rows = kk ;

  if ( ag == 0 )
    fprintf(fileptr,"\nThe rows are nuclear genotypes and the columns are cytotypes.\n");
  else
    fprintf(fileptr,"\nThe rows are nuclear alleles and the columns are cytotypes.\n");


  fprintf(fileptr,"\n\\begin{center}\n\\begin{tabular}{");
  for ( ii = 0 ; ii <= cols+1 ; ii++ )
    fprintf(fileptr,"c");
  fprintf(fileptr,"}\\hline\\hline\n");
  for ( ii = 1 ; ii <= cols ; ii++ )
    fprintf(fileptr,"  &  %2d    ",ii);
  fprintf(fileptr," & Margins \\\\ \\hline \n");
  if ( ag == 0 ) {
  iii = 0;
  for ( ii = 1 ; ii <= kk ; ii++ )
    for ( jj = ii ; jj <= kk ; jj++ ) {
      iii = iii+1;
      fprintf(fileptr," %2d,%-2d ",ii,jj);
      for ( jjj = 1 ; jjj <= cols ; jjj++ ) 
        fprintf(fileptr," & %9.4f",*(*(dblptr+jjj)+iii) );
      fprintf(fileptr," & %9.4f \\\\ \n",*(*(dblptr+0)+iii) );
    }
  }
  else {
    for ( jj = 1 ; jj <= kk ; jj++ ) {
      fprintf(fileptr," %2d    ", jj);
      for ( jjj = 1 ; jjj <= cols ; jjj++ ) 
        fprintf(fileptr," & %9.4f",*(*(dblptr+jjj)+jj) );
      fprintf(fileptr," & %9.4f \\\\ \n",*(*(dblptr+0)+jj) );
    }
  }
  fprintf(fileptr,"\\hline\n");
  fprintf(fileptr,"Margins ");
  for ( ii = 1 ; ii <= cols ; ii++ )
    fprintf(fileptr," & %9.4f",*(*(dblptr+ii)+0));
  fprintf(fileptr," & %9.4f \\\\ \n",*(*(dblptr+0)+0));
  fprintf(fileptr,"\\hline\\hline\n\\end{tabular}\n\\end{center}\n\n");


}

void  calc_MonteCarlo_prs(thedata,numsets,outfile,errorfile,reps,bb,cc,verbosity,tex_flag,ag)  
  adataset *thedata;
  int numsets;
  char *outfile,*errorfile;
  long reps,bb,cc;
  int verbosity,tex_flag,ag; 
{
  int **table,ii,jj,iset,rows,cols,*cmargins,*rmargins,**mptr;
  long extremectr,lii,lxcntr;
  FILE *fileptr;
  double shuffle_pr,data_pr,hh,tot,tot2,var;
  ldiv_t xx;
  var = 0.0;
  for ( iset = 1 ; iset <= numsets ; iset++ ) {
    tot = tot2 = 0.0;
    if ( verbosity == 1 )
          printf("\n Now doing Monte Carlo Shuffle test for Data Set %d\n\twith %ld reps in %ld batches of %ld observations\n",iset,reps,bb,cc);
    rows = (thedata+iset)->mm;
    if ( ag == 0 ) {
      mptr = (thedata+iset)->counts;
      cols = (thedata+iset)->kk * ( (thedata+iset)->kk + 1 ) / 2;
    }
    else{
      mptr = (thedata+iset)->kmcounts;
      cols = (thedata+iset)->kk  ;
    }
    cmargins = ivector(1,cols);
    rmargins = ivector(1,rows);
    table = imatrix(0,rows,0,cols);
    /* fill in the table */
    for ( ii = 0 ; ii <= rows ; ii++ ) 
      for ( jj = 0 ; jj <= cols ; jj++ ) 
        *(*(table+ii)+jj) = *(*( mptr+ii)+jj);
    data_pr = prob_sample(table,rows,cols);
    extremectr = 0L; lxcntr = 0L;
    for ( lii = 1 ; lii <= reps ; lii++ ) {
      shuffle_matrix(table,rows,cols,cmargins,rmargins);
      shuffle_pr =  prob_sample(table,rows,cols);
      if ( shuffle_pr <= data_pr ) {
        extremectr += 1L;
        lxcntr += 1L;
      }
      xx = ldiv(lii,cc);
      if ( xx.rem == 0 ) {
        hh = (double) lxcntr/ (double) cc ;
        tot = tot + hh;
        tot2 = tot2 + hh*hh;
        lxcntr = 0L;
        if ( verbosity == 1 ) {
          printf(".");
          fflush(stdout);
        }
      }      
    }
    if ( verbosity == 1 )
          printf("\n" );
    var = (tot2 -  (tot*tot / ((double) bb)) )/( (double) (bb*(bb-1L)) )  ;
    if ( var > 0.0 )
      var = sqrt(var);
    free_ivector(cmargins,1,cols);
    free_ivector(rmargins,1,rows);
    free_imatrix(table,0,rows,0,cols);
    (thedata+iset)->MonteCarlo = (double) extremectr / (double) reps ;
    
    fileptr = fileopen(outfile,"a");
    if ( tex_flag == 0 ) {
	    fprintf(fileptr,"\n\n");
	    for ( jj = 1 ; jj <= 80 ; jj++ ) 
	      fprintf(fileptr,"=");
	    fprintf(fileptr,"\n");
	    if ( ag == 0 )
	      fprintf(fileptr,"\n\n\tAnalysis of nuclear system %s by cytotype %s : Genotypic Disequilibria\n\n",(thedata+iset)->nnames,(thedata+iset)->cnames);
	    else
	      fprintf(fileptr,"\n\n\tAnalysis of nuclear system %s by cytotype %s : Allelic Disequilibria\n\n",(thedata+iset)->nnames,(thedata+iset)->cnames);
	    write_counts((thedata+iset),fileptr,ag);
	    if ( var > 0.0 )
	      fprintf(fileptr,"\n\nMonte Carlo Shuffle test pr = %12.8f with sample standard deviation %12.8f",(thedata+iset)->MonteCarlo,var);
	    else
	      fprintf(fileptr,"\n\nMonte Carlo Shuffle test pr = %12.8f with sample variance %12.8f",(thedata+iset)->MonteCarlo,var);
	    fprintf(fileptr,"\nFor %ld / %ld\n\n",extremectr,reps);
	    for ( jj = 1 ; jj <= 80 ; jj++ ) 
	      fprintf(fileptr,"=");
	    fprintf(fileptr,"\n");
	}
	else {
	  if ( ag == 0 )
	    fprintf(fileptr,"\n\n\\subsection{Monte Carlo Simulation for Genotypic Disequilibria}  ");
	  else
	    fprintf(fileptr,"\n\n\\subsection{Monte Carlo Simulation for Allelic Disequilibria}  ");
	    fprintf(fileptr," This is the analysis of nuclear system %s by cytotype %s.",(thedata+iset)->nnames,(thedata+iset)->cnames);
	    write_counts_tex((thedata+iset),fileptr,ag);
	    if ( var > 0.0 )
	      fprintf(fileptr," The Monte Carlo Shuffle test probability was  %f\nwith sample standard deviation of %f.",(thedata+iset)->MonteCarlo,var);
	    else
	      fprintf(fileptr,"\n\nThe Monte Carlo Shuffle test probability was  %f\nwith sample variance of %f.",(thedata+iset)->MonteCarlo,var);
	    fprintf(fileptr,"\nThere were %ld observations as rare as the data out of %ld repetitions.",extremectr,reps);
	    fprintf(fileptr,"\nThe standard error (sampling variance) was based on %ld batches of %ld observations.",bb,cc);

	}
    fileclose(outfile,fileptr);

  }
}



/*
  This  returns the  log of the  exact 
  probability of the sample in the array mm, constrained by it's 
  margins.  

  mm = imatrix(0,rows,0,cols)

  The 0th row from col 1 to cols contains the marginal counts of
  the columns.


  The 0th column from row 1 to rows contains the marginal counts of
  the rows.
  
  *(*(mm+0)+0) contains the sample size, n

  The log of the exact probability is num-den, where 

  num = Sum ln(marginals!)

  and 

  den = ln( n ! ) + Sum ln( counts! )

*/

double prob_sample(mm,rows,cols)
  int **mm,rows,cols;
{
  int ii,jj;
  double den,num;
  num  = 0.0;
  for ( ii = 1 ; ii <= rows ; ii++ ) 
    num = num + lnfactrl(*(*(mm+ii)+0) );
  for ( ii = 1 ; ii <= cols ; ii++ ) 
    num = num + lnfactrl(*(*(mm+0)+ii) ) ;
  den = lnfactrl( *(*(mm+0)+0)) ;
  for ( ii = 1 ; ii <= rows ; ii++ ) 
    for ( jj = 1 ; jj <= cols ; jj++ ) 
      den = den + lnfactrl( *(*(mm+ii)+jj)  );
  return(  num-den  );
}




/*
  This shuffles the matrix mm constrained by the margins.
  mm = imatrix(0,rows,0,cols)
  the 0th column from 1 to rows has the marginal counts for the rows.
  the 0th row from 1 to cols has the marginal counts for the columns.
  element *(*(mm+0)+0) has the sample size.
*/
void shuffle_matrix(mm,rows,cols,cmargins,rmargins)
  int **mm,rows,cols,*cmargins,*rmargins;
{
  int  nn,ii,jj;
  int arow,acol,total;
  double tmp;
  nn = *(*(mm+0)+0);
  for ( ii = 1 ; ii <= rows ; ii++ ) 
    *(rmargins+ii) = *(*(mm+ii)+0);
  for ( ii = 1 ; ii <= cols ; ii++ ) 
    *(cmargins+ii) = *(*(mm+0)+ii);
  for ( ii = 1 ; ii <= rows ; ii++ )
    for ( jj = 1 ; jj <= cols ; jj++ )
      *(*(mm+ii)+jj) = 0;

  for ( ii = nn ; ii > 0 ; ii-- ) {
    tmp =  ranf(ii);
    arow = (int) ceil( (double) ii * tmp ) ;
    tmp =  ranf(ii);
    acol = (int) ceil( (double) ii * tmp) ;
    total = 0;
    for ( jj = 1 ; jj <= rows && arow > total ; jj++ )
      total = total + *(rmargins+jj);
    arow = jj - 1;
    *(rmargins+arow) -= 1;
    total = 0;
    for ( jj = 1 ; jj <= cols && acol > total ; jj++ )
      total = total + *(cmargins+jj);
    acol = jj - 1;
    *(cmargins+acol) -= 1;
    *(*(mm+arow)+acol) += 1;
  }
}

int samplesize(d0,d1,dalt)
  double d0,d1,dalt;
{
  double temp;
  int back;
  back = 1;
  if (dalt == 0 ) 
    back =  -7;
  else 
    temp =   ( (double) Zalpha * d0 + (double) Zbeta * d1 )  / dalt;
  if ( temp > 99.0 || temp < -99.0 ) /* 99^2 < 10,000 = 100^2 */
    back =  9999;
  else if ( back > 0 ) {
    temp = temp*temp;
    back = roundit(temp);
  }
  return back;
}


int newsamplesize(d0,d1,dalt,zalph)
  double d0,d1,dalt,zalph;
{
  double temp;
  int back;
  back = 1;
  if (dalt == 0 ) 
    back =  -7;
  else 
    temp =   ( zalph * d0 + (double) Zbeta * d1 )  / dalt;
  if ( temp > 99.0 || temp < -99.0 ) /* 99^2 < 10,000 = 100^2 */
    back =  9999;
  else if ( back > 0 ) {
    temp = temp*temp;
    back = roundit(temp);
  }
  return back;
}



int roundit(value)
  double value;
{
  int retval;
  double t;
  retval = (int) value;
  t = (double) retval;
  if ( value-t > 0.5 ) 
    retval = retval+1;
  return retval;
}

