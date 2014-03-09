/*  
          Mdatain.c,  subroutines to read datafiles for CNDm.
          
Copyright (C) 1996 Christopher J. Basten.

This file is part of CNDm. CNDm is free software; you
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
#include "Mmain.h"


int     get_numsets(inputf,errfile)
  char   *inputf,*errfile;
{
  FILE   *errorf, *infile;
  int     ii,jj, xnumsets;
  char buffer[MAXLINE],ch,xtemp[MAXNAME];
  infile = fileopen(inputf, "r");
  if (infile == NULL) return( -1);
  do {
    for ( ii = 0 ; ii < MAXLINE ; ii++ ) *(buffer+ii) = '\0';
    for ( ii = 0 ;  ((ch=fgetc(infile)) != EOF) && (ch != '\n') ; ii++ )  *(buffer+ii) = ch;
    if ( *(buffer+0) == '-' && *(buffer+1) == 'l' ) {
      fileclose(inputf,infile);
      get_field(2,xtemp,buffer);
      xnumsets = atoi(xtemp);
          errorf = fileopen(errfile, "a");
          fprintf(errorf,"\n%d set(s) of joint cytonuclear counts\n",xnumsets);
          fileclose(errfile,errorf);
      return(xnumsets);
    }
  } while ( ch != EOF );
  fileclose(inputf,infile); 
  errorf = fileopen(errfile,"a");
  fprintf(errorf,"\nCouldn't get the number of datasets...\n");
  fileclose(errfile,errorf); 
  return(-1);
}


adataset  *datavector(nl, nh)
  int     nl, nh;
{
  adataset  *v;
  int ii,jj;
  v = (adataset *) malloc((unsigned) (nh - nl + 1) * sizeof(adataset));
  if (!v)
    nrerror("allocation failure in datavector()");
  for ( ii = 0 ; ii <= nh-nl ; ii ++ ) {
    (v+ii)->kk = (v+ii)->mm = 2;
    (v+ii)->nn = 0;
    (v+ii)->MonteCarlo = (v+ii)->Markov = -1.0;
    (v+ii)->acounts = NULL;
    (v+ii)->kmcounts = (v+ii)->counts = NULL;
    (v+ii)->diseq_ests = (v+ii)->ndiseq_ests = (v+ii)->freq_ests = (v+ii)->var0_ests = NULL;  
    (v+ii)->var1_ests = (v+ii)->stat_ests = (v+ii)->pr_ests = NULL; 
    (v+ii)->pstar_ests = NULL;       
    for ( jj = 0 ; jj < MAXNAME ; jj++ ) 
      *((v+ii)->cnames+jj) = *((v+ii)->nnames+jj) = '\0';
  }

  return v - nl;
}

void    free_datavector(v, nl, nh)
  adataset  *v;
  int     nl, nh;
{
  int ii,jj,kk,mm,rows,cols;
  for ( ii = nl ; ii <= nh ; ii ++ ) {
    kk = (v+ii)->kk;
    rows = mm = (v+ii)->mm;
    cols = kk*(kk+1)/2;
    if ( (v+ii)->freq_ests != NULL )
      free_dmatrix((v+ii)->freq_ests,0,rows,0,cols);
    if ( (v+ii)->var0_ests != NULL )
      free_dmatrix((v+ii)->var0_ests,0,rows,0,cols);
    if ( (v+ii)->var1_ests != NULL )
      free_dmatrix((v+ii)->var1_ests,0,rows,0,cols);
    if ( (v+ii)->stat_ests != NULL )
      free_dmatrix((v+ii)->stat_ests,0,rows,0,cols);
    if ( (v+ii)->pr_ests != NULL )
      free_dmatrix((v+ii)->pr_ests,0,rows,0,cols);
    if ( (v+ii)->diseq_ests != NULL )
      free_dmatrix((v+ii)->diseq_ests,0,rows,0,cols);
    if ( (v+ii)->ndiseq_ests != NULL )
      free_dmatrix((v+ii)->ndiseq_ests,0,rows,0,cols);
    if ( (v+ii)->counts != NULL )
      free_imatrix((v+ii)->counts,0,rows,0,cols);   
    if ( (v+ii)->kmcounts != NULL )
      free_imatrix((v+ii)->kmcounts,0,mm,0,kk);   
    if ( (v+ii)->acounts != NULL )
      free_ivector((v+ii)->acounts,1,kk);
    if ( (v+ii)->pstar_ests != NULL )
      free_dvector((v+ii)->pstar_ests,1,kk);
  }
  free((char *) (v + nl));
}



adataset *get_data(inputf,errfile,numsets)
  char   *inputf,*errfile;
  int   numsets;
{
  adataset *thedata;
  FILE   *errorf, *infile;
  int     ii,jj, xnumsets,total,row,col,cols,field;
  char buffer[MAXLINE],ch,xtemp[MAXNAME];
  long pos;

  thedata = datavector(1,numsets);
  infile = fileopen(inputf, "r");
  if (infile == NULL) return(NULL);
  jj = 0;
  do {
    for ( ii = 0 ; ii < MAXLINE ; ii++ ) *(buffer+ii) = '\0';
    for ( ii = 0 ;  ((ch=fgetc(infile)) != EOF) && (ch != '\n') ; ii++ )  *(buffer+ii) = ch;
    if ( *(buffer+0) == '-' ) 
      switch ( *(buffer+1) ) {
        case 'C' : 
          get_field(2,xtemp,buffer);
          strcpy((thedata+jj)->cnames, xtemp);  
          break;
        case 'N' : 
          get_field(2,xtemp,buffer);
          strcpy((thedata+jj)->nnames, xtemp); 
          break;
        case 'H' :
        case 'h' :
/*
          pos = ftell(infile);
          fileclose(inputf,infile);
*/
          cols = (thedata+jj)->kk * ((thedata+jj)->kk + 1) / 2 ;
          (thedata+jj)->counts = imatrix(0,(thedata+jj)->mm,0,cols);
/*
          infile = fileopen(inputf,"r");
          fseek(infile,pos,0);
*/
          break;
        case 'M' :
        case 'm' :
        case '0' : 
        case '1' : 
        case '2' : 
        case '3' : 
        case '4' : 
        case '5' : 
        case '6' : 
        case '7' : 
        case '8' : 
        case '9' :        
          get_field(1,xtemp,buffer);
          *(xtemp+0) = ' ';
          if ( *(xtemp+1) == 'M' )
            row = 1;
          else if ( *(xtemp+1) == 'm' )
            row = 2;
          else
            row = atoi(xtemp);
          for ( field = 1 ; field <= cols ; field++ ) {
            get_field(field+1,xtemp,buffer);
            *(*( (thedata+jj)->counts+row)+field) = atoi(xtemp);
          }
          break;
        case 'e' :
          for ( row = 1 ; row <= (thedata+jj)->mm ; row++) {
            total = 0;
            for ( col = 1 ; col <= cols ; col++ )
              total = total + *(*( (thedata+jj)->counts+row)+col) ;
            *(*( (thedata+jj)->counts+row)+0) = total;
          }
          for ( col = 0 ; col <= cols ; col++) {
            total = 0;
            for ( row = 1 ; row <= (thedata+jj)->mm ; row++)
              total = total + *(*( (thedata+jj)->counts+row)+col) ;
            *(*( (thedata+jj)->counts+0)+col) = total;
          }
          (thedata+jj)->nn = *(*( (thedata+jj)->counts+0)+0) ;
          break;
        case 'p' :  
          write_counts((thedata+jj),stdout);
          break;
        case 'a' :
          get_field(2,xtemp,buffer);
          (thedata+jj)->kk = atoi(xtemp);
          break;
        case 'q' :
          fileclose(inputf,infile); 
          errorf = fileopen(errfile, "a");
          fprintf(errorf,"\nRead in the data AOK...\n");
          fileclose(errfile,errorf);
          if ( jj == numsets ) return(thedata);
          else return(NULL);
          break;
        case 'r' :
          get_field(2,xtemp,buffer);
          (thedata+jj)->mm = atoi(xtemp);
          break;
        case 'b' : 
          jj = jj + 1;
          break;
        default :
          break;
       }
  } while ( ch != EOF );
  fileclose(inputf,infile); 
  if ( ch == EOF ) {
          errorf = fileopen(errfile, "a");
          fprintf(errorf,"\nTime to close the input file after an EOF...\n");
          fileclose(errfile,errorf);
  }
  return(NULL);
}


void write_counts(thedata,fileptr,ag)
adataset *thedata;
FILE *fileptr;
int ag;
{
  int rows,cols,iii,ii,jj,jjj;
  int **mptr;
  fprintf(fileptr,"\n");
  for ( ii = 1 ; ii <= 80 ; ii++ ) 
    fprintf(fileptr,"*");
  fprintf(fileptr,"\n");
  cols = thedata->mm;
  if ( ag == 0 ) {
    mptr = thedata->counts;
    rows = thedata->kk*(thedata->kk+1)/2;
    fprintf(fileptr,"\n\n\tJoint cytonuclear counts...");
    fprintf(fileptr,"\n\t\tNuclear Genotypes (rows) by Cytotypes (columns)\n       |");
  }
  else {
    mptr = thedata->kmcounts;
    rows = thedata->kk ;
    fprintf(fileptr,"\n\n\tJoint cytonuclear counts...");
    fprintf(fileptr,"\n\t\tNuclear Alleles (rows) by Cytotypes (columns)\n       |");
  }
  for ( ii = 1 ; ii <= cols ; ii++ )
    fprintf(fileptr," %4d",ii);
  fprintf(fileptr," | Margin \n");
  for ( ii = 1 ; ii <= 5*(cols+1)+10 ; ii++ )
    fprintf(fileptr,"-");
  fprintf(fileptr,"\n");
  iii = 0;
  if ( ag == 0 ) {
    for ( ii = 1 ; ii <= thedata->kk ; ii++ )
      for ( jj = ii ; jj <= thedata->kk ; jj++ ) {
        iii = iii+1;
        fprintf(fileptr," %2d,%-2d |",ii,jj);
        for ( jjj = 1 ; jjj <= cols ; jjj++ ) 
          fprintf(fileptr," %4d",*(*(mptr+jjj)+iii) );
        fprintf(fileptr," | %4d\n",*(*(mptr+0)+iii) );
      }
  }
  else {
      for ( jj = 1 ; jj <= thedata->kk ; jj++ ) {
        fprintf(fileptr," %2d    |",jj );
        for ( jjj = 1 ; jjj <= cols ; jjj++ ) 
          fprintf(fileptr," %4d",*(*(mptr+jjj)+jj) );
        fprintf(fileptr," | %4d\n",*(*(mptr+0)+jj) );
      }
  }
  for ( ii = 1 ; ii <= 5*(cols+1)+10 ; ii++ )
    fprintf(fileptr,"-");
  fprintf(fileptr,"\n");
  fprintf(fileptr,"Margins|");
  for ( ii = 1 ; ii <= cols ; ii++ )
    fprintf(fileptr," %4d",*(*(mptr+ii)+0));
  fprintf(fileptr," | %4d\n\n",*(*(mptr+0)+0));
  for ( ii = 1 ; ii <= 80 ; ii++ ) 
    fprintf(fileptr,"*");
  fprintf(fileptr,"\n");
}

void write_counts_tex(thedata,fileptr,ag)
adataset *thedata;
FILE *fileptr;
int ag;
{
  int rows,cols,iii,ii,jj,jjj,**mptr;
  fprintf(fileptr,"\n");

  cols = thedata->mm;
  if ( ag == 0 ) {
    mptr = thedata->counts;
    rows = thedata->kk*(thedata->kk+1)/2;
    fprintf(fileptr,"\n Joint cytonuclear counts where");
    fprintf(fileptr,"\n the nuclear genotypes are in rows and cytotypes are in columns. \n ");
    fprintf(fileptr,"\n\\begin{center}\n\\begin{tabular}{");
  }
  else {
    mptr = thedata->kmcounts;
    rows = thedata->kk ;
    fprintf(fileptr,"\n Joint cytonuclear counts where");
    fprintf(fileptr,"\n the nuclear alleles are in rows and cytotypes are in columns. \n ");
    fprintf(fileptr,"\n\\begin{center}\n\\begin{tabular}{");
  }
  for ( ii = 0 ; ii <= cols+1 ; ii++ )
    fprintf(fileptr,"c");
  fprintf(fileptr,"}\\hline\\hline\n");
  for ( ii = 1 ; ii <= cols ; ii++ )
    fprintf(fileptr," & %4d",ii);
  fprintf(fileptr," & Margin \\\\ \\hline \n");
  iii = 0;
  if ( ag == 0 ) {
    for ( ii = 1 ; ii <= thedata->kk ; ii++ )
      for ( jj = ii ; jj <= thedata->kk ; jj++ ) {
        iii = iii+1;
        fprintf(fileptr," %2d,%-2d &",ii,jj);
        for ( jjj = 1 ; jjj <= cols ; jjj++ ) 
          fprintf(fileptr," %4d &",*(*(mptr+jjj)+iii) );
        fprintf(fileptr,"   %4d \\\\\n",*(*(mptr+0)+iii) );
      }
  }
  else {
    for ( jj = 1 ; jj <= thedata->kk ; jj++ ) {
      fprintf(fileptr," %2d    &",  jj);
      for ( jjj = 1 ; jjj <= cols ; jjj++ ) 
        fprintf(fileptr," %4d &",*(*(mptr+jjj)+jj) );
      fprintf(fileptr,"   %4d \\\\\n",*(*(mptr+0)+jj) );
    }
  }
  fprintf(fileptr,"\\hline\n");
  fprintf(fileptr,"Margins");
  for ( ii = 1 ; ii <= cols ; ii++ )
    fprintf(fileptr,"& %4d",*(*(mptr+ii)+0));
  fprintf(fileptr," & %4d \\\\ ",*(*(mptr+0)+0));
  fprintf(fileptr,"\\hline\\hline\n\\end{tabular}\n\\end{center}\n\n");
}

