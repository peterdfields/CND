#include "Dmain.h"

anestim  *anestvector(nl, nh)
  int     nl, nh;
{
  anestim  *v;
  v = (anestim *) malloc((unsigned) (nh - nl + 1) * sizeof(anestim));
  if (!v)
    nrerror("allocation failure in vector()");
  return v - nl;
}

void    free_anestvector(v, nl, nh)
  anestim  *v;
  int     nl, nh;
{
  free((char *) (v + nl));
}



int     get_the_data(counts,inputf,errfile,cnames,nnames,numsets)
  anestim *counts;
  char   *inputf,*errfile,**cnames,**nnames;
  int   numsets;
{
  FILE   *errorf, *infile;
  int     ii,jj, xnumsets;
  char buffer[MAXLINE],ch,xtemp[MAXNAME];
  infile = fileopen(inputf, "r");
  if (infile == NULL) return -1;
  jj = 0;
  do {
    for ( ii = 0 ; ii < MAXLINE ; ii++ ) *(buffer+ii) = '\0';
    for ( ii = 0 ;  ((ch=fgetc(infile)) != EOF) && (ch != '\n') ; ii++ )  *(buffer+ii) = ch;
/*    *(buffer+ii) = ' ';*/
    if ( *(buffer+0) == '-' ) 
      switch ( *(buffer+1) ) {
        case 'l' :
          if ( numsets == 0 ) {
            fileclose(inputf,infile);
            get_field(2,xtemp,buffer);
            xnumsets = atoi(xtemp);
            printf("\n%d set(s) of joint cytonuclear counts\n",xnumsets);
            return xnumsets;
          }
          break;
        case 'C' : 
          get_field(2,xtemp,buffer);
          strcpy(*(cnames+jj), xtemp);  
          break;
        case 'N' : 
          get_field(2,xtemp,buffer);
          strcpy(*(nnames+jj), xtemp);  
          break;
        case 'M' : 
          get_field(2,xtemp,buffer);
          (counts+jj)->nAAM = atoi(xtemp);
          get_field(3,xtemp,buffer);
          (counts+jj)->nAaM = atoi(xtemp);
          get_field(4,xtemp,buffer);
          (counts+jj)->naaM = atoi(xtemp);
          break;
        case 'm' :  
          get_field(2,xtemp,buffer);
          (counts+jj)->nAAm = atoi(xtemp);
          get_field(3,xtemp,buffer);
          (counts+jj)->nAam = atoi(xtemp);
          get_field(4,xtemp,buffer);
          (counts+jj)->naam = atoi(xtemp);
          break;
        case 'q' :
          fileclose(inputf,infile); 
          printf("\nRead in the data...\n");
          if ( jj == numsets ) return 0;
          else return -1;
          break;
        case 'p' : 
          fprintf(stderr,"\nData set number: %d\nMitochondrial locus: %s\nNuclear locus: %s",jj,*(cnames+jj),*(nnames+jj));
          fprintf(stderr,"\n\n    AA   Aa   aa");
          fprintf(stderr,"\nM %4d %4d %4d",(counts+jj)->nAAM,(counts+jj)->nAaM,(counts+jj)->naaM);
          fprintf(stderr,"\nm %4d %4d %4d\n",(counts+jj)->nAAm,(counts+jj)->nAam,(counts+jj)->naam);
          break;
        case 'b' : 
          jj = jj + 1;
          clear_anestim((counts+jj));
          break;
        default :
          break;
       }
  } while ( ch != EOF );
  if ( ch == EOF ) printf("\nTime to close the input file after an EOF...\n");
  fileclose(inputf,infile); 
  return -1;
}

int calc_the_est(counts,numsets)
  anestim *counts;
  int numsets;
{
  int ii;
  for ( ii = 1 ; ii <= numsets ; ii++ ) {
    counts[ii].samplesize = counts[ii].nAAM + counts[ii].nAaM + counts[ii].naaM +
    counts[ii].nAAm +counts[ii].nAam + counts[ii].naam;
    counts[ii].PAAM = (double) counts[ii].nAAM / (double) counts[ii].samplesize;
    counts[ii].PAaM = (double) counts[ii].nAaM / (double) counts[ii].samplesize;
    counts[ii].PaaM = (double) counts[ii].naaM / (double) counts[ii].samplesize;
    counts[ii].PAAm = (double) counts[ii].nAAm / (double) counts[ii].samplesize;
    counts[ii].PAam = (double) counts[ii].nAam / (double) counts[ii].samplesize;
    counts[ii].Paam = (double) counts[ii].naam / (double) counts[ii].samplesize;
    counts[ii].PAA = counts[ii].PAAM + counts[ii].PAAm;
    counts[ii].PAa = counts[ii].PAaM + counts[ii].PAam;
    counts[ii].Paa = counts[ii].PaaM + counts[ii].Paam;
    counts[ii].PA = counts[ii].PAA + counts[ii].PAa/2;
    counts[ii].PM = counts[ii].PAAM + counts[ii].PAaM + counts[ii].PaaM;
    counts[ii].PAM = counts[ii].PAAM + counts[ii].PAaM/2;
    counts[ii].rho1 = counts[ii].PAAM * counts[ii].PAam / (counts[ii].PAaM * counts[ii].PAAm ) ;
    counts[ii].rho2 = counts[ii].PAaM * counts[ii].Paam / (counts[ii].PaaM * counts[ii].PAam ) ; 


  }
  return 0;
}
  
void clear_anestim(theest)
  anestim *theest;
{
    theest->samplesize = 0;
    theest->PAAM = 0;
    theest->PAaM = 0;
    theest->PaaM = 0;
    theest->PAAm = 0;
    theest->PAam = 0;
    theest->Paam = 0;
    theest->PAA = 0;
    theest->PAa = 0;
    theest->Paa = 0;
    theest->PA = 0;
    theest->PM = 0;
    theest->PAM = 0;
}
