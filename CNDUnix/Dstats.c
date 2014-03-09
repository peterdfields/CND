#include "Dmain.h"

#ifndef M_PI
#define  M_PI    3.14159265358979323846
#endif

double *factorials, **normals; 

char Exact(int [][2], float *both_tail, float *closest);

double lappr_fct (double i);
double lfct( double x);

/*

   Drive the data analysis for real data.  
*/
void DoTheStats(char *outfile, char *errfile, char *inputf, char *prog,  anestim   *counts, int numsets, char **nnames,char **cnames, int tex_flag) {
  FILE *outf;
  char *tptr;    
  anestim theest ;
  int     ii,l1,indicator;
  int    *samples,*sample;
  double *vars,*diseq,*corrs,*prcorrs,*sdevs,*sdevs1,*vars1,*prexact;
    
  outf = fileopen(outfile, "a");
  if ( outf == NULL ) outf = fileopen(outfile, "w");
  if (outf != NULL) {
    tptr = asctime2();
    if ( tex_flag == 0 )
      fprintf(outf,"\n\n\tIt is %s\n\n",tptr);
    else {
      fprintf(outf,"\\documentclass[12pt]{article}");
	  fprintf(outf,"\n\\setlength{\\textwidth}{6.5in}");
	  fprintf(outf,"\n\\setlength{\\textheight}{8.0in}");
	  fprintf(outf,"\n\\setlength{\\topmargin}{0.0in}");
	  fprintf(outf,"\n\\setlength{\\oddsidemargin}{0in}"); 
	  fprintf(outf,"\n\\begin{document}");
	  fprintf(outf,"\n\\centerline{{\\large {\\bf CNDd Analysis}}}\n\\vspace{3ex}\n");	  
	  fprintf(outf,"\n The analysis of program %s",prog);
      fprintf(outf, " was run at %s\n",tptr);
      fprintf(outf,"\nThe output was written to %s,",outfile);
      fprintf(outf,"errors were logged to %s, and ",errfile);
      fprintf(outf,"the data came from %s.\n",inputf);
    }
    
    for ( l1 = 1 ; l1 <= numsets ; l1++ ) {
	      theest = *(counts+l1);
          diseq = calc_diseq(theest);
          indicator = 5;
          vars = calc_the_vars(theest,indicator);
          for ( ii = 1 ; ii <= 5 ; ii++ )
            if ( *(vars+ii) < 0.0 )
              *(vars+ii) = 0.0;
          corrs = calc_corrs(diseq,vars);
          prcorrs = calc_prs(theest,corrs);

          if ( *(prcorrs+3) >= 0.05  ) {
            indicator = 4;
            free_dvector(vars,1,5);
            vars = calc_the_vars(theest,indicator);
            for ( ii = 1 ; ii <= 5 ; ii++ )
              if ( *(vars+ii) < 0.0 )
                *(vars+ii) = 0.0;
            free_dvector(corrs,1,5);
            corrs = calc_corrs(diseq,vars);
            free_dvector(prcorrs,1,5);
            prcorrs = calc_prs(theest,corrs);
          }
          prexact = calc_prs_X(theest);

          sdevs = calc_sdevs(vars);
          vars1 = calc_the_vars(theest,indicator+10);
          for ( ii = 1 ; ii <= 5 ; ii++ )
            if ( *(vars1+ii) < 0.0 )
              *(vars1+ii) = 0.0;
          sdevs1 = calc_sdevs(vars1);
          samples = calc_ss(theest,diseq,sdevs,sdevs1,1.0);
          sample = calc_ss(theest,diseq,sdevs,sdevs1,0.0);
          for ( ii = 1 ; ii < 6 ; ii++ ) if ( *(diseq+ii) == 0 ) {
            *(sdevs+ii) = *(vars1+ii) = *(vars+ii) = *(sdevs1+ii) = corrs[ii] = *(prcorrs+ii) = 0.0;
            *(sample+ii) = *(samples+ii) = 0;
          }
        if ( tex_flag == 0 ) {
	      fprintf(outf,"\n\t %s by %s\n",*(nnames+l1),*(cnames+l1));
          fprintf(outf,"\tWhere pA = %4.2f, pM = %4.2f and the sample size is %3d\n",theest.PA,theest.PM,theest.samplesize);
          for (ii=1;ii<76;ii++) fprintf(outf,"-");
          fprintf(outf,"+\n               |     DA    |     DAM   |    DAAM   |    DAaM   |    DaaM   |\n");
          for (ii=1;ii<76;ii++) fprintf(outf,"-");
          fprintf(outf,"+\nEstimators     | %9.4g | %9.4g | %9.4g | %9.4g | %9.4g |",*(diseq+1),*(diseq+2),*(diseq+3),*(diseq+4),*(diseq+5));

          calc_diseq_prime(theest,diseq);
          fprintf(outf,"\nNorm. Est.     | %9.4g | %9.4g | %9.4g | %9.4g | %9.4g |",*(diseq+1),*(diseq+2),*(diseq+3),*(diseq+4),*(diseq+5));
          fprintf(outf,"\nVariances      | %9.4g | %9.4g | %9.4g | %9.4g | %9.4g |",*(vars+1),*(vars+2),*(vars+3),*(vars+4),*(vars+5));
          fprintf(outf,"\nVariances (T)  | %9.4g | %9.4g | %9.4g | %9.4g | %9.4g |",*(vars1+1),*(vars1+2),*(vars1+3),*(vars1+4),*(vars1+5));
          fprintf(outf,"\nSt. Errors     | %9.4g | %9.4g | %9.4g | %9.4g | %9.4g |",*(sdevs+1),*(sdevs+2),*(sdevs+3),*(sdevs+4),*(sdevs+5));
          fprintf(outf,"\nSt. Errors (T) | %9.4g | %9.4g | %9.4g | %9.4g | %9.4g |",*(sdevs1+1),*(sdevs1+2),*(sdevs1+3),*(sdevs1+4),*(sdevs1+5));
          fprintf(outf,"\nTest Statistic | %9.4g | %9.4g | %9.4g | %9.4g | %9.4g |",*(corrs+1),*(corrs+2),*(corrs+3),*(corrs+4),*(corrs+5));
          fprintf(outf,"\nPr(Sample)     | %9.4g | %9.4g | %9.4g | %9.4g | %9.4g |",*(prcorrs+1),*(prcorrs+2),*(prcorrs+3),*(prcorrs+4),*(prcorrs+5));
          fprintf(outf,"\nPr(Exact )     | %9.4g | %9.4g | %9.4g | %9.4g | %9.4g |",*(prexact+1),*(prexact+2),*(prexact+3),*(prexact+4),*(prexact+5));
          fprintf(outf,"\nSample Sizes 90| %9d | %9d | %9d | %9d | %9d |",*(samples+1),*(samples+2),*(samples+3),*(samples+4),*(samples+5));
          fprintf(outf,"\nSample Sizes 50| %9d | %9d | %9d | %9d | %9d |",*(sample+1),*(sample+2),*(sample+3),*(sample+4),*(sample+5));
          fprintf(outf,"\n");
          for (ii=1;ii<76;ii++) fprintf(outf,"-");
          fprintf(outf,"+\nThe (counts) vector:   (nAAM, nAaM, naaM, nAAm, nAam, naam)");
          fprintf(outf,"\nThe (counts) values:   (%4d, %4d, %4d, %4d, %4d, %4d)\n",theest.nAAM,theest.nAaM,theest.naaM,theest.nAAm,theest.nAam,theest.naam);
          fprintf(outf,"\n\n");
        }
        else {
	      fprintf(outf,"\n\\section{%s by %s}\n\n",*(nnames+l1),*(cnames+l1));
          fprintf(outf,"Here $p^A = %4.2f$, $p_M = %4.2f$ and the sample size is %3d.\n",theest.PA,theest.PM,theest.samplesize);
          fprintf(outf,"\n\\begin{center}\n\\begin{tabular}{lccccc}\\hline\\hline");
          fprintf(outf,"\n               & $D^A$ & $D^A_M$ & $D^{AA}_M$ & $D^{Aa}_M$   & $D^{aa}_M$ \\\\ \\hline");
          fprintf(outf,"\nEstimators $\\tilde{D}$ & %9.4g & %9.4g & %9.4g & %9.4g & %9.4g \\\\",*(diseq+1),*(diseq+2),*(diseq+3),*(diseq+4),*(diseq+5));
          calc_diseq_prime(theest,diseq);
          fprintf(outf,"\nNorm. Est. $\\tilde{D^\\prime}$ & %9.4g & %9.4g & %9.4g & %9.4g & %9.4g \\\\",*(diseq+1),*(diseq+2),*(diseq+3),*(diseq+4),*(diseq+5));
          fprintf(outf,"\nVariances $\\sigma^2_0$     & %9.4g & %9.4g & %9.4g & %9.4g & %9.4g \\\\",*(vars+1),*(vars+2),*(vars+3),*(vars+4),*(vars+5));
          fprintf(outf,"\nVariances (T)  $\\sigma^2_1$ & %9.4g & %9.4g & %9.4g & %9.4g & %9.4g \\\\",*(vars1+1),*(vars1+2),*(vars1+3),*(vars1+4),*(vars1+5));
          fprintf(outf,"\nSt. Errors  $\\sigma_0$   & %9.4g & %9.4g & %9.4g & %9.4g & %9.4g \\\\",*(sdevs+1),*(sdevs+2),*(sdevs+3),*(sdevs+4),*(sdevs+5));
          fprintf(outf,"\nSt. Errors (T) $\\sigma_1$ & %9.4g & %9.4g & %9.4g & %9.4g & %9.4g \\\\",*(sdevs1+1),*(sdevs1+2),*(sdevs1+3),*(sdevs1+4),*(sdevs1+5));
          fprintf(outf,"\nTest Statistic $n \\tilde{r}^2$ & %9.4g & %9.4g & %9.4g & %9.4g & %9.4g \\\\",*(corrs+1),*(corrs+2),*(corrs+3),*(corrs+4),*(corrs+5));
          fprintf(outf,"\nPr(Sample)  $\\chi^2_1$   & %9.4g & %9.4g & %9.4g & %9.4g & %9.4g \\\\",*(prcorrs+1),*(prcorrs+2),*(prcorrs+3),*(prcorrs+4),*(prcorrs+5));
          fprintf(outf,"\nPr(Exact)      & %9.4g & %9.4g & %9.4g & %9.4g & %9.4g \\\\",*(prexact+1),*(prexact+2),*(prexact+3),*(prexact+4),*(prexact+5));
          fprintf(outf,"\nSample Sizes 90\\%% & %9d & %9d & %9d & %9d & %9d \\\\",*(samples+1),*(samples+2),*(samples+3),*(samples+4),*(samples+5));
          fprintf(outf,"\nSample Sizes 50\\%% & %9d & %9d & %9d & %9d & %9d \\\\",*(sample+1),*(sample+2),*(sample+3),*(sample+4),*(sample+5));
          fprintf(outf,"\\hline\\hline\n\\end{tabular}\n\\end{center}\n\n");
          fprintf(outf,"\nThe counts vector is   $(n^{AA}_M,\\,\\, n^{Aa}_M,\\,\\, n^{aa}_M,\\,\\, n^{AA}_m,\\,\\, n^{Aa}_m,\\,\\, n^{aa}_m) =$");
          fprintf(outf," $(%4d, %4d, %4d, %4d, %4d, %4d)$\n",theest.nAAM,theest.nAaM,theest.naaM,theest.nAAm,theest.nAam,theest.naam);
          fprintf(outf,"\n\n");
        }
          free_dvector(sdevs,1,5);
          free_dvector(sdevs1,1,5);
          free_dvector(prcorrs,1,5);
          free_dvector(prexact,1,5);
          free_dvector(diseq,1,5);
          free_dvector(vars,1,5);
          free_dvector(vars1,1,5);
          free_dvector(corrs,1,5);
          free_ivector(samples,1,5);
          free_ivector(sample,1,5);
 	}
      }
      tptr = asctime2();
      fprintf(outf, "\n\nThis analysis finished at %s\n",tptr);

       if ( tex_flag == 0 )
         fprintf(outf,"\n\n");
       else    
         fprintf(outf,"\n\\end{document}\n");  

         
  fileclose(outfile, outf);
}



void write_anest(theest)
  anestim theest;
{
  fprintf(stdout,"\n\n");
  fprintf(stdout,"\npA   = %10.5f         pM   = %10.5f         pAM  = %10.5f",theest.PA,theest.PM,theest.PAM);
  fprintf(stdout,"\npAA  = %10.5f       pAa  = %10.5f       paa  = %10.5f",theest.PAA,theest.PAa,theest.Paa);
  fprintf(stdout,"\npAAM = %10.5f       pAaM = %10.5f       paaM = %10.5f",theest.PAAM,theest.PAaM,theest.PaaM);
  fprintf(stdout,"\nnAAM = %10d       nAaM = %10d       naaM = %10d",theest.nAAM,theest.nAaM,theest.naaM);
  fprintf(stdout,"\nnAAm = %10d       nAam = %10d       naam = %10d",theest.nAAm,theest.nAam,theest.naam);
  fprintf(stdout,"\n%5d\n",theest.samplesize);
}

void calc_diseq_prime(theest,values)
  anestim theest;
  double *values;
{
  double   const1, const2, PAA, PAa, Paa;
  PAA = theest.PAAM + theest.PAAm;
  PAa = theest.PAaM + theest.PAam;
  Paa = theest.PaaM + theest.Paam;
/*  DA prime */
  const1 = 1.0;
  if ( *(values+1) < 0.0 ) {
    const1 = theest.PA * theest.PA ;
    const2 = (1.0 - theest.PA) * (1.0 - theest.PA) ;
    if ( const2 < const1 )
      const1 = const2;
  }
  else if ( *(values+1) > 0.0 ) 
    const1 = theest.PA * (1.0 - theest.PA);
  if ( const1 != 0.0 ) 
    *(values+1) = *(values+1) / const1 ;
/*  DAM prime */
  const1 = 1.0;
  if ( *(values+2) < 0.0 ) {
    const1 = theest.PA * theest.PM;
    const2 = (1.0 - theest.PA) * (1.0 - theest.PM);
    if ( const2 < const1 ) 
      const1 = const2;
    const2 = 0.5 * ( PAA * theest.PM + Paa * (1.0 - theest.PM) );
    if ( const2 < const1 ) 
      const1 = const2;
  }
  else if ( *(values+2) > 0.0 ) {
    const1 = theest.PA *(1.0 - theest.PM) ;
    const2 = (1.0 - theest.PA) * theest.PM ;
    if ( const2 < const1 ) 
      const1 = const2;
    const2 = 0.5 * ( PAA * (1.0 - theest.PM) + Paa * theest.PM );
    if ( const2 < const1 ) 
      const1 = const2;
  }
  if ( const1 != 0.0 ) 
    *(values+2) = *(values+2) / const1 ;
/*  DAAM prime */
  const1 = 1.0;
  if ( *(values+3) < 0.0 ) {
    const1 = PAA * theest.PM;
    const2 = (1.0 - PAA) * (1.0 - theest.PM);
    if ( const2 < const1 ) 
      const1 = const2;
  }
  else if ( *(values+3) > 0.0 ) {
    const1 = PAA *(1.0 - theest.PM) ;
    const2 = (1.0 - PAA) * theest.PM ;
    if ( const2 < const1 ) 
      const1 = const2;
  }
  if ( const1 != 0.0 ) 
    *(values+3) = *(values+3) / const1 ;
/*  DAaM prime */
  const1 = 1.0;
  if ( *(values+4) < 0.0 ) {
    const1 = PAa * theest.PM;
    const2 = (1.0 - PAa) * (1.0 - theest.PM);
    if ( const2 < const1 ) 
      const1 = const2;
  }
  else if ( *(values+4) > 0.0 ) {
    const1 = PAa *(1.0 - theest.PM) ;
    const2 = (1.0 - PAa) * theest.PM ;
    if ( const2 < const1 ) 
      const1 = const2;
  }
  if ( const1 != 0.0 ) 
    *(values+4) = *(values+4) / const1 ;
/*  DaaM prime */
  const1 = 1.0;
  if ( *(values+5) < 0.0 ) {
    const1 = Paa * theest.PM;
    const2 = (1.0 - Paa) * (1.0 - theest.PM);
    if ( const2 < const1 ) 
      const1 = const2;
  }
  else if ( *(values+5) > 0.0 ) {
    const1 = Paa *(1.0 - theest.PM) ;
    const2 = (1.0 - Paa) * theest.PM ;
    if ( const2 < const1 ) 
      const1 = const2;
  }
  if ( const1 != 0.0 ) 
    *(values+5) = *(values+5) / const1 ;
}


double *calc_diseq(theest)
  anestim theest;
{
  double *values;
  values = dvector(1,5);
  *(values+1) = theest.PAA - theest.PA*theest.PA;
  *(values+2) = theest.PAM - theest.PA*theest.PM;
  *(values+3) = theest.PAAM - theest.PAA*theest.PM;
  *(values+4) = theest.PAaM - theest.PAa*theest.PM;
  *(values+5) = theest.PaaM - theest.Paa*theest.PM;
  return values;
}


double *calc_the_vars(theest,ind)
  anestim theest;
  int ind;
{
  double *values,rem,dA,dAM,dAAM,dAaM,daaM;
  double pA,pM,pAA,pAa,paa;
  double ipAA,ipAa,ipaa,idA,idAAM,idAaM,idaaM;
  int nn;
  values = dvector(1,5);
  nn = theest.samplesize;
  pA = theest.PA;
  pM = theest.PM;
  pAA = theest.PAA;
  pAa = theest.PAa;
  paa = theest.Paa;
  dA = pAA-pA*pA;
  dAAM =  theest.PAAM - pAA*pM;
  dAaM =  theest.PAaM - pAa*pM;
  daaM =  theest.PaaM - paa*pM;
  dAM = theest.PAM - pA*pM;
  if ( ind == 1 || ind > 10 ) idA = dA;
  else idA = 0.0;
  *(values+1) = (pA * ( 1.0 - pA ) * pA * ( 1.0 - pA ) + idA*(1-2.0*pA)*(1-2.0*pA)- idA*idA)/(double) nn;

  if ( ind == 0 || ind == 5 || ind == 15 ) 
    rem = dA*pM*(1-pM) + dAAM*(1-2.0*pM);  
  else if ( ind == 1 || ind == 11 ) 
    rem = dA*pM*(1-pM) + dAAM*(1-2.0*pM) + dAM*( 1.0 - 4.0*pA )*(1-2.0*pM) - 2.0*dAM*dAM;  
  else if ( ind == 2 || ind == 12 ) 
    rem = 0.0;  
  else if ( ind == 3 || ind == 13 ) 
    rem = dAAM*(1-2.0*pM);  
  else if ( ind == 4 || ind == 14 ) 
    rem = dA*pM*(1-pM); 
  if ( ind > 11 ) rem = rem + dAM*( 1.0 - 4.0*pA )*(1-2.0*pM) - 2*dAM*dAM; 

  *(values+2) = 0.5 * (pA * ( 1.0 - pA ) * pM * ( 1.0 - pM ) + rem)/(double) nn;
 
  if ( ind == 1 || ind > 10 ) {
    idAAM = dAAM;
    idAaM = dAaM;
    idaaM = daaM;
  }
  else {
    idAAM = 0.0;
    idAaM = 0.0;
    idaaM = 0.0;
  }
  ipAA = pAA;
  ipAa = pAa;
  ipaa = paa;

  *(values+3) = (ipAA * ( 1.0 - ipAA ) * pM * ( 1.0 - pM ) + idAAM*( 1.0 - 2.0* ipAA )*(1-2.0*pM) - idAAM*idAAM)/(double) nn;
  *(values+4) = (ipAa * ( 1.0 - ipAa ) * pM * ( 1.0 - pM ) + idAaM*( 1.0 - 2.0* ipAa )*(1-2.0*pM) - idAaM*idAaM)/(double) nn;
  *(values+5) = (ipaa * ( 1.0 - ipaa ) * pM * ( 1.0 - pM ) + idaaM*( 1.0 - 2.0* ipaa )*(1-2.0*pM) - idaaM*idaaM)/(double) nn;
  return values;
}

int *calc_ss(theest,diseq,sdevs0,sdevs1,fact)
  anestim theest;
  double *diseq,*sdevs0,*sdevs1,fact;
{
/* fact determines whether beta is 0.9 (fact = 1.0) or 
                                   0.5 (fact = 0.0)
*/
   double del1,del0;
   int *values;
   values = ivector(1,5);

   del0 = *(sdevs0+1) * sqrt(theest.samplesize);
   del1 = *(sdevs1+1) * sqrt(theest.samplesize) * fact;
   *(values+1) = samplesize(del0,del1,*(diseq+1));

   del0 = *(sdevs0+2) * sqrt(theest.samplesize);
   del1 = *(sdevs1+2) * sqrt(theest.samplesize) * fact;
   *(values+2) = samplesize(del0,del1, *(diseq+2));

   del0 = *(sdevs0+3) * sqrt(theest.samplesize);
   del1 = *(sdevs1+3) * sqrt(theest.samplesize) * fact;
   *(values+3) = samplesize(del0,del1, *(diseq+3));

   del0 = *(sdevs0+4) * sqrt(theest.samplesize);
   del1 = *(sdevs1+4) * sqrt(theest.samplesize) * fact;
   *(values+4) = samplesize(del0,del1, *(diseq+4));

   del0 = *(sdevs0+5) * sqrt(theest.samplesize);
   del1 = *(sdevs1+5) * sqrt(theest.samplesize) * fact;
   *(values+5) = samplesize(del0,del1, *(diseq+5));

  return values;
}

int samplesize(d0,d1,dalt)
  double d0,d1,dalt;
{
  double temp;
  int back;
  back = 1;
  if (dalt == 0 ) back =  -7;
  else temp =   ( (double) Zalpha * d0 + (double) Zbeta * d1 )  / dalt;
  if ( temp > 1000 ) back =  -2;
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
  if ( value-t > 0.5 ) retval = retval+1;
  return retval;
}


double *calc_corrs(diseq,vars)
  double *diseq,*vars;
{
  int ii;
  double *values;
  values = dvector(1,5);
  for ( ii = 1 ; ii <= 5 ; ii++)
    if ( *(vars+ii) > 0.0 ) values[ii] = *(diseq+ii) * *(diseq+ii) / *(vars+ii);
    else values[ii] = 99.99;
  return values;
}

double *calc_sdevs(vars)
  double *vars;
{
  int ii;
  double *values;
  values = dvector(1,5);
  for ( ii = 1 ; ii <= 5 ; ii++ )
    if ( *(vars+ii) > 0 )  values[ii] = sqrt(*(vars+ii));
    else if ( *(vars+ii) == 0.0 )
      values[ii] = 0.0;
    else
      values[ii] = -.01;

  return values;
}

double *calc_prs(theest,corrs)
  anestim theest;
  double *corrs;
{
  int ii;
  double *values,pa,pm;
  pa = theest.PA;
  pm = theest.PM; /*
  if ( 0.2 < pa && pa < 0.8 && 0.15 < pm && pm < 0.85 ) {

*/

    values = dvector(1,5);
    for ( ii = 1 ; ii <= 5 ; ii++ ) 
      values[ii] = chiprob(1,corrs[ii]);
/*  }
  else

    values = calc_prs_X(theest);*/
  return values;
}

double *calc_prs_X(theest)
  anestim theest;
{
  int table[2][2] ;
  float bothtails, closesttail;
  double *values;

  values = dvector(1,5);

  *(values+1) = -0.01;

  table[0][0] = 2.0*theest.nAAM +     theest.nAaM;
  table[0][1] =     theest.nAaM + 2.0*theest.naaM;
  table[1][0] = 2.0*theest.nAAm +     theest.nAam;
  table[1][1] =     theest.nAam + 2.0*theest.naam;
  if (Exact(table, &bothtails, &closesttail) != -1)
    *(values+2) = bothtails;
  else
    *(values+2) = -0.01;

  table[0][0] = theest.nAAM;
  table[0][1] = theest.nAaM + theest.naaM;
  table[1][0] = theest.nAAm;
  table[1][1] = theest.nAam + theest.naam;
  if (Exact(table, &bothtails, &closesttail) != -1)
    *(values+3) = bothtails;
  else
    *(values+3) = -0.01;

  table[0][0] = theest.nAaM;
  table[0][1] = theest.nAAM + theest.naaM;
  table[1][0] = theest.nAam;
  table[1][1] = theest.nAAm + theest.naam;
  if (Exact(table, &bothtails, &closesttail) != -1)
    *(values+4) = bothtails;
  else
    *(values+4) = -0.01;

  table[0][0] = theest.naaM;
  table[0][1] = theest.nAAM + theest.nAaM;
  table[1][0] = theest.naam;
  table[1][1] = theest.nAAm + theest.nAam;
  if (Exact(table, &bothtails, &closesttail) != -1)
    *(values+5) = bothtails;
  else
    *(values+5) = -0.01;

  return values;
}


/*
  Minimum sample size required to detect any disequilibrium based on
  the exact test for 2 x 3 tables.  From
    
     M. S. Sanchez, C. J. Basten, A. M. Ferrenberg, M. A. Asmussen and
     J. Arnold, 2004.  Exact sample sizes needed to detect dependence 
     in 2 x 3 Tables.  TPB, submitted.

The raw data are sample counts for each of the six cells:  (nAAM, nAaM, naaM, nAAm, nAam, naam)
theest has estimates from a data set of the probability model and the margins.
We calculate rho1, rho2 and rho3 as well as the disequilibria.  
We want to know the minimum sample size to detect the observed disequilibria with
power 1-beta and size alpha.     


Assume that *factorials and *normals have been created.  
*/
int SanchezSampleSize(anestim *theest, anestim *workest, double alpha, double beta, char *outfile) {
  int minss, lbn, ubn, midn , asympn;
  double  DAAM, DAaM, DaaM,midpower,lpower,upower,lsize,usize,msize,adjust;
  exactnode *first;
  FILE *outf;
  char *theoutfile;

  first = NULL;
  if ( debug == 0 ) 
    theoutfile = NULL;
  else 
    theoutfile = outfile;
  adjust = 0.2;
  minss = 0;
  workest->PAAM = theest->PAAM;
  workest->PAaM = theest->PAaM;
  workest->PaaM = theest->PaaM;
  workest->PAAm = theest->PAAm;
  workest->PAam = theest->PAam;
  workest->Paam = theest->Paam;
  workest->PAA = theest->PAA;
  workest->PA = theest->PA;
  workest->PM = theest->PM;
  workest->PAM = theest->PAM;
  workest->PAa = theest->PAa;
  workest->Paa = theest->Paa;
  workest->rho1 = theest->PAAM * theest->PAam / (theest->PAaM * theest->PAAm ) ; 
  workest->rho2 = theest->PAaM * theest->Paam / (theest->PaaM * theest->PAam ) ; 

  DAAM = theest->PAAM - theest->PAA*theest->PM;
  DAaM = theest->PAaM - theest->PAa*theest->PM;
  DaaM = theest->PaaM - theest->Paa*theest->PM;

  asympn = workest->samplesize = midn =lbn = ubn =   AsympSampleSize(*workest,alpha, beta);
  
  if ( midn > (int) MAXEXACTSS )
    minss = -midn;
  else {        
	  first = alloc_enode();
	  first = SanchezPower( theoutfile, workest, first, alpha, workest->rho1, workest->rho2, &midpower);
	  msize = workest->esize;
	  
	  if ( midpower > 1.0-beta && midn < (int) MAXEXACTSS) { /* try a smaller sample size */
	    do {
	      workest->samplesize = lbn = (int) ((1.0 - adjust) * (double) lbn) ;
	      first = SanchezPower( theoutfile, workest, first, alpha, workest->rho1, workest->rho2, &lpower);
	      lsize = workest->esize;
	    } while ( lpower > 1.0-beta );
	  }
	  else if ( midpower < 1.0-beta && midn < (int) MAXEXACTSS) { /* try a larger sample size */
	    do {
	      workest->samplesize = ubn =  (int) ((1.0 + adjust) * (double) ubn) ;
	      if (ubn <= (int) MAXEXACTSS) 
	        first = SanchezPower( theoutfile, workest, first, alpha, workest->rho1, workest->rho2, &upower);
	      else if ( midn <= (int) MAXEXACTSS ) {
	        workest->samplesize = ubn = (int) MAXEXACTSS;
	        first = SanchezPower( theoutfile, workest, first, alpha, workest->rho1, workest->rho2, &upower); 
	        if ( upower < 1.0-beta) {  /* we only need do this one time.  If not big enough, return the asymp. ss. */
	          minss = -midn;     
	          upower = 2.0-beta;
	        }
	      }
	      else { /*  we get here if midn is larger than MAXEXACTSS and still not large enough.  */
	          minss = -midn;     
	          upower = 2.0-beta;
	      }
	      usize = workest->esize;
	    } while ( upower < 1.0-beta );
	  }
	  else /*  we get here if midn is larger than MAXEXACTSS and still not large enough.  */
	    minss = -midn;
  }
/*  if minss == 0 then we have succeeded in bounding the sample size and can narrow in on it. 
    if not, then we skip this part.    */
  if ( minss == 0 ) {
    do {
	    if ( lbn == midn ) {
	      lpower = midpower;
	      lsize = msize;
	      midn = lbn + (int) ( (double)(ubn-lbn) * (1.0-beta-lpower)/(upower-lpower) );
	      if ( midn == lbn )
	        midn = midn+1;
	      else if (midn == ubn)
	        midn = ubn - 1; 
	      workest->samplesize = midn ;
	      first = SanchezPower( theoutfile, workest, first, alpha, workest->rho1, workest->rho2, &midpower);	 
	      msize = workest->esize;     
	    }
	    else {
	      upower = midpower;
	      usize = msize;
	      midn = lbn + (int) ( (double)(ubn-lbn) * (1.0-beta-lpower)/(upower-lpower) );
	      if ( midn == lbn )
	        midn = midn+1;
	      else if (midn == ubn)
	        midn = ubn - 1; 
	      workest->samplesize = midn ;
	      first = SanchezPower( theoutfile, workest, first, alpha, workest->rho1, workest->rho2, &midpower);     
	      msize = workest->esize; 
	    }
	    if ( midpower < 1.0-beta) {
	        lpower = midpower;
	        lsize = msize;
	        lbn=midn;
	    }
	    else {
	        upower = midpower;
	        usize = msize;
	        ubn=midn;	
	    }        
    } while ( (ubn-lbn) > 1 );
    if ( fabs(upower-1.0+beta) < fabs(lpower-1.0+beta) ) {
      minss = ubn;
      workest->eN = minss;
      workest->epower = upower;
      workest->esize = usize;
    }
    else {
      minss = lbn;
      workest->eN = minss;
      workest->epower = lpower;
      workest->esize = lsize;
    }
  }
  if ( first != NULL )
      dealloc_enodes(first);
  
  
  if ( theoutfile != NULL ) {  
        outf = fileopen(theoutfile, "a");
        SanchezShowModel(outf, workest, DAAM, DAaM );
        fprintf(outf,"\n   Exact sample size:  %d\n#\n#", minss );
        fileclose(theoutfile, outf);  
  }
  
  workest->samplesize = asympn;
  if ( minss < 0 ) {
    workest->epower = 1.0-beta;
    workest->esize = alpha;
  }
  return(minss);
}

/*
  Asymptotic Sample size calcuated from a model in a workest.  
  
*/
int AsympSampleSize(anestim workest,double alpha, double beta) {
   int i,sAAM,sAaM,saaM;
  double  thez, zalpha2,zbeta;
  
  zalpha2 = Phi( 1.0 - alpha/2.0 ); /* two sided test of size alpha*/
 
  zbeta = Phi( 1.0 - beta );   /* power is beta. */
 
  sAAM = AsympSS( workest.PAAM, workest.PAA, workest.PM, zalpha2, zbeta);
  sAaM = AsympSS( workest.PAaM, workest.PAa, workest.PM, zalpha2, zbeta);
  saaM = AsympSS( workest.PaaM, workest.Paa, workest.PM, zalpha2, zbeta);
  if ( sAaM < sAAM )
    sAAM = sAaM;
  if ( saaM < sAAM )
    sAAM = saaM;
  
  return(sAAM);
}

/*
   Calculate the asymptotic sample size needed for size alpha, power 1-beta
   2-sided test.   The form of the variance for each of DAAM, DAaM and DaaM are
   the same.   
*/
int AsympSS(double paam, double paa, double pm, double zalpha2, double zbeta) {
  double D, delta0, delta1,ss;
  
  D = paam - paa * pm ;
  delta0 = sqrt( paa * (1.0 -paa) * pm * (1.0 - pm) );
  delta1 = sqrt( delta0*delta0 + D*(1.0-2.0*paa)*(1.0 - 2.0*pm) - D*D );  
  if ( D != 0.0 )
    ss = (zbeta*delta1 + zalpha2*delta0)/ D;
  else
    ss = pow(2.0, 15.0);
  ss = ss*ss;
  return( (int) ss ); 
}

/*

  Power for a given sample size required to detect any disequilibrium based on
  the exact test for 2 x 3 tables.  From
    
     M. S. Sanchez, A. M. Ferrenberg, M. A. Asmussen, C. J. Basten and
     J. Arnold, 2004.  Exact sample sizes needed to detect dependence 
     in 2 x 3 Tables.  TPB, submitted.



Given a sample size, a model (contained in workest along with rho1 and rho2), calculate the power.    


*/
exactnode *SanchezPower(char *outfile, anestim *workest, exactnode *first, double size, double rho1, double rho2, double *thepower) {
  int i,j,k,lbnAAM,ubnAAM,lbnAaM, ubnAaM;
  double thesize, maxp,thisp,sampprob, margprob, apower,zalpha2,thez,  truemargsum ;
  FILE *outf;
  long start,finish;

  zalpha2 = Phi(1.0-size/2.0);
  start = marktime();
  thesize = *thepower = 0.0;
  for ( i=0; i<=workest->samplesize; i++ )  /*  Loop over all margins.  */
    for ( j=0; j<=workest->samplesize; j++ )
      for ( k=0; k<=workest->samplesize-j; k++ ) {
        maxp = SanchezMaxMargProb(*workest,workest->samplesize,i,j,k);
        if ( maxp > margthresh ) {  /* This margin has sufficient probability to consider. */
          first = clean_enodes(first);
          margprob=SanchezMargProb(*workest, workest->samplesize,i,j,k,rho1,rho2);
          lbnAAM = i+j-workest->samplesize;          
          if ( lbnAAM < 0 )
            lbnAAM = 0;
          ubnAAM = i;
          if (ubnAAM > j)
            ubnAAM = j;
          truemargsum = 0.0;
          for ( workest->nAAM=lbnAAM; workest->nAAM<= ubnAAM ; workest->nAAM++ ) {
            lbnAaM = i+j+k-workest->samplesize-workest->nAAM;
            if ( lbnAaM<0)
              lbnAaM = 0;
            ubnAaM=i-workest->nAAM;
            if ( ubnAaM>k)
              ubnAaM=k;
            for ( workest->nAaM=lbnAaM; workest->nAaM<=ubnAaM ; workest->nAaM++ ) {
              workest->naaM = i - workest->nAAM - workest->nAaM;
              workest->nAAm = j - workest->nAAM;
              workest->nAam = k - workest->nAaM;
              workest->naam = workest->samplesize - j - k - workest->naaM;           
              thisp = SanchezNullProb(*workest);
              sampprob = SanchezSampCondProb(*workest, rho1, rho2, thisp);
              truemargsum += sampprob;
              if (thisp <= size ) {		                
                push_enode(first, thisp, sampprob);
                if ( first->nullprob > size )
                  pull_enode(first); /* get rid of the biggest one. */
              }
              else if (  thisp > 1.0000000001 && debug != 0) 
                printf("\n %2d %2d %2d %2d %2d %2d %2d %2d %2d %2d %f",workest->nAAM, workest->nAaM, workest->naaM, workest->nAAm, workest->nAam, workest->naam, i, j, k, workest->samplesize, thisp);
            }
                    
          }
          *thepower =  *thepower + first->condprob * margprob;
          thesize = thesize + truemargsum * margprob * first->nullprob ;  
         }   /* End of the if maxp > margthresh */
      }
  apower = AsympPower( workest, zalpha2);
  workest->epower = *thepower;
  workest->esize = thesize;
  workest->eN = workest->samplesize;
  if ( outfile != NULL ) {
    finish = marktime();
    outf = fileopen(outfile, "a");
    fprintf(outf, "\n-n  %4d -epower  %6.4f  -esize   %6.4f -apower   %6.4f  -etime %ld", workest->eN ,  *thepower,  thesize, apower, finish-start);
    fileclose(outfile, outf);
  }
  first = clean_enodes(first);
  return(first);
}

/*
   calculate the power based on asymptotic results.   Do for each of DAAM, DAaM and DaaM and return the
   largest.
   See Asmussen and Basten (1994) Table 5 and equation 26.  
   
   delta0 and delta1 are defined in Table 5.
   the power defined by 
   
   
   phi( (sqrt(n)*D - z_a/2 * delta0)  / delta1 ) = 1 - beta
   
   
   and we return the biggest of the three.  
   
   
   
*/
double AsympPower(  anestim *workest , double zalpha2) {
  double  powerAAM, powerAaM, poweraaM;

  powerAAM = asympow(workest->PAAM,workest->PAA, workest->PM, workest->samplesize, zalpha2);
  powerAaM = asympow(workest->PAaM,workest->PAa, workest->PM, workest->samplesize, zalpha2);
  poweraaM = asympow(workest->PaaM,workest->Paa, workest->PM, workest->samplesize, zalpha2);
  

  if ( powerAaM > powerAAM )
    powerAAM = powerAaM;
  if ( poweraaM > powerAAM )
    powerAAM = poweraaM;  

  return(powerAAM);
  
  
}

/*

Let z be z of alpha/2 (right tail of standard normal with alpha/2 area)
d0 and d1 are delta0 and delta1

we fail to reject in [ -z d0/sqrt(n) , z d0/sqrt(n) ]

D is the estimate of D

we reject if D  < -z d0/sqrt(n) or  D > z d0/sqrt(n)

D1 is the hypothosized alternate value of D, and has s.e. d1/sqrt(n)

1-beta = P( D  < -z d0/sqrt(n) | D1) + P( D > z d0/sqrt(n) | D1) 
      
       = P( Z < (-z d0 - D1 sqrt(n)) / d1 )  + P(Z > (z d0 - D1 sqrt(n))/d1 )
       = P( Z < (-z d0 - D1 sqrt(n)) / d1 )  + 1.0 - P(Z < (z d0 - D1 sqrt(n))/d1 )
       = Phi((-z d0 - D1 sqrt(n)) / d1 ) +1.0 -  Phi((z d0 - D1 sqrt(n))/d1)
       
The trick is to work with the positive version of D, since they sum to zero 
in columns.   If D<0, then use paa-paam, paa, 1-pm and -D.  
       
       
*/
double asympow(double paam, double paa, double pm, int nn, double zalpha2) {
  int i;
  double D, delta0, delta1, thezl,thezu, power,pl,pu;
  double pAAM, pAA, pM;
  pAAM = paam;
  pAA = paa;
  pM = pm;
    
  D = pAAM - pAA * pM ;
  if ( D < 0 )  {
    pAAM = paa - paam;
    pM = 1.0 - pm;
    D =  pAAM - pAA * pM ;
  }
    
  
  delta0 = sqrt( pAA * (1.0 -pAA) * pM * (1.0 - pM) );
  delta1 = sqrt( delta0*delta0 + D*(1.0-2.0*pAA)*(1.0 - 2.0*pM) - D*D );  
  
  thezl = ( -zalpha2 * delta0 - D * sqrt( (double) nn ) )/delta1 ;
  thezu = (  zalpha2 * delta0 - D * sqrt( (double) nn ) )/delta1 ;
  pl = PhiProb(thezl);
  pu = PhiProb(thezu);
  power = 1.0 - pu + pl  ;   
  return(power);
}

/*
   Calculate the constant C from the paper.   return log(C)  
   
   can calculate from all counts or just marginal counts.   Use margins if
   iNN != 0 
*/
double SanchezC( anestim workest, int iNN, int in1d, int ind1, int ind2 ) {
  double cc;
  double n1d, nd1, nd2, NN;

	if ( iNN == 0 ) {  
	  n1d = (double) (workest.nAAM + workest.nAaM + workest.naaM) ;
	  nd1 = (double) (workest.nAAM + workest.nAAm ) ;  
	  nd2 = (double) (workest.nAaM + workest.nAam ) ;  
	  NN =(double) (workest.samplesize); 
	  
	} else {
	  n1d = (double) in1d ;
	  nd1 = (double) ind1;  
	  nd2 = (double) ind2 ;  
	  NN =(double) iNN; 

	}  
	cc =      (2.0*n1d+3.0*nd1-NN)*log(workest.PAAM) + (2.0*n1d+3.0*nd2-NN)*log(workest.PAaM) ;
	cc = cc + (2.0*n1d-3.0*nd1-3.0*nd2+2.0*NN)*log(workest.PaaM) + (-2.0*n1d+3.0*nd1+NN)*log(workest.PAAm) ;
	cc = cc + (-2.0*n1d+3.0*nd2+NN)*log(workest.PAam) + (-2.0*n1d-3.0*nd1-3.0*nd2+4.0*NN)*log(workest.Paam) ;
	cc = cc/6.0;
    return( cc );
   
}



/*

   Probability of a sample conditional on the margins.  In Eqn (10), this is the term in
   the double sum at the right.  
*/
double SanchezSampCondProb( anestim workest, double rho1, double rho2, double nullprob ) {

  double rp,   nd1,nd2,nd3, n11, n12, n13  ;

  nd1 = factorials[ (workest.nAAM + workest.nAAm ) ] ;
  nd2 = factorials[ (workest.nAaM + workest.nAam ) ];
  nd3 = factorials[ (workest.naaM + workest.naam ) ];
  n11 = factorials[ workest.nAAM ] ;
  n12 = factorials[ workest.nAaM ];
  n13 = factorials[ workest.naaM ];
    rp =  nd1 + nd2 + nd3 - n11 - n12 - n13  ;

  n11 = factorials[ workest.nAAm ] ;
  n12 = factorials[ workest.nAam ];
  n13 = factorials[ workest.naam ];
    rp = rp - n11 - n12 - n13;
  
  n11 = (double) workest.nAAM;
  n12 = (double) workest.nAaM;
    rp = rp + n11*log(rho1) + (n11+n12)*log(rho2);

  return( exp( rp ) );

}




/*
   Calculate maximum contribution of the margin to the power.  
   If iNN = 0, then assume that workest has data in n..
   otherwise, assume that NN and margins in1d, ind1, ind2 are given.
   
   workest must have a model.  
   
   This is equation 7 in the paper.
   
*/
double SanchezMaxMargProb( anestim workest, int iNN, int in1d, int ind1, int ind2 ) {

  double maxp,cc, n1d, nd1, nd2,nd3, n2d, NN;

	if ( iNN == 0 )  {
	  n1d = factorials[(workest.nAAM + workest.nAaM + workest.naaM)];
	  n2d = factorials[(workest.nAAm + workest.nAam + workest.naam)];
	  nd1 = factorials[(workest.nAAM + workest.nAAm)] ;
	  nd2 = factorials[(workest.nAaM + workest.nAam)];  
	  nd3 = factorials[(workest.naaM + workest.naam)];  
	  NN  = factorials[(workest.samplesize)]; 
	}
	else {  
	  n1d = factorials[in1d];
	  n2d = factorials[(iNN - in1d)];
	  nd1 = factorials[ind1] ;
	  nd2 = factorials[ind2];  
	  nd3 = factorials[(iNN - ind1 - ind2)];  
	  NN  = factorials[iNN]; 
	}
    cc = SanchezC( workest, iNN, in1d, ind1, ind2); 
    maxp =  cc + NN - n1d - n2d + NN - nd1 - nd2 - nd3  ;  
  
  
    return( exp( maxp ) );
  
}


/*
  Calculate the probability of a Margin.  Actually, this is the part of equation 10 to the left of
  the double sum.
*/
double SanchezMargProb( anestim workest, int iNN, int in1d, int ind1, int ind2, double rho1, double rho2 ) {

  double rp,cc, n1d, nd1, nd2,nd3, n2d, NN;


	  n1d = factorials[in1d];
	  n2d = factorials[(iNN - in1d)];
	  nd1 = factorials[ind1];
	  nd2 = factorials[ind2];  
	  nd3 = factorials[(iNN - ind1 - ind2)];  
	  NN  = factorials[iNN]; 

  cc = SanchezC( workest, iNN, in1d, ind1, ind2); 
  rp = cc + NN  - nd1 - nd2 - nd3 + ( (double)iNN / 6.0 - (double)in1d / 3.0 - (double)ind1/2.0) * log( rho1 ) ; 
  rp = rp + ( (double)iNN / 3.0 - 2.0 * (double)in1d / 3.0 - (double)(ind1+ind2) /2.0) * log( rho2 )  ;  
  
  
  return( exp( rp ) );
  
}

/*
   Calculate the probability of a sample { nij } given the margins { n.1, n.2, n1., N }
   under the null hypothesis.  

   workest must counts  { nij }   Equation 9
*/
double SanchezNullProb( anestim workest  ) {

  double rp, n1d, nd1, nd2,nd3,   n2d, NN;
  n1d = factorials[(workest.nAAM + workest.nAaM + workest.naaM)];
  n2d = factorials[(workest.nAAm + workest.nAam + workest.naam)];
  nd1 = factorials[(workest.nAAM + workest.nAAm )];
  nd2 = factorials[(workest.nAaM + workest.nAam )];  
  nd3 = factorials[(workest.naaM + workest.naam )];  
  NN  = factorials[(workest.samplesize) ]; 
  rp = nd1 + nd2 + nd3 + n1d + n2d - NN;
  rp = rp - factorials[workest.nAAM] - factorials[workest.nAaM] - factorials[workest.naaM];
  rp = rp - factorials[workest.nAAm] - factorials[workest.nAam] - factorials[workest.naam];  
  return( exp( rp ) );
  
}
/*
   Given the margins and the rho values, solve for D1 and D2 (DAAM and DAaM).
*/
void SanchezSolveDiseq(double rho1, double rho2, double p1d, double pd1, double pd2, double *d1, double *d2) {
  double mtemp, mind, maxd,dthresh ;
  double c1, c2, c3,c4,rho3,pd3,p2d ;   
  rho3 = rho1 * rho2;
  pd3 = 1.0-(pd1+pd2);
  p2d = 1.0-p1d;
  dthresh = 10E-15;
/*First, iterate on DAAm.   It is bounded by -min[p2. p.1, (1-p2.)(1-p.1)] and min[(1-p2.)p.1, p2. (1-p.1)]*/  
  mind =  p2d * pd1;
  mtemp = (1-p2d) * (1-pd1);
  if (mtemp < mind )
    mind = mtemp;
  mind = -mind; 
  maxd =  (1-p2d) * pd1;
  mtemp = p2d * (1-pd1);
  if (mtemp < maxd )
    maxd = mtemp;
/* original: c1 D^3 + c2 D^2 + c3 D + c4 */
  c1= (1.0-rho1)*(1.0-rho3);
  c2= -( (1.0-rho1)*( (p1d+rho3*p2d)*pd1+(p2d+rho3*p1d)*pd3 )  + (1.0-rho3)*( (p1d+rho1*p2d)*pd1+(p2d+rho1*p1d)*pd2 )  );
  c3= ((p1d+rho1*p2d)*(p1d+rho3*p2d)*pd1+(p1d+rho3*p2d)*(p2d+rho1*p1d)*pd2+(p1d+rho1*p2d)*(p2d+rho3*p1d)*pd3+(1.0-rho1)*(1.0-rho3)*p1d*p2d*(1.0-pd1))*pd1;
  c4= -((1.0-rho1)*(p1d+rho3*p2d)*pd2+(1.0-rho3)*(p1d+rho1*p2d)*pd3)*p1d*p2d*pd1*pd1;
/*  Divide all by c1.  */  
  c2=c2/c1;
  c3=c3/c1;
  c4=c4/c1;    /*  now, d^3 + c2 d^2 + c3 d + c4 = 0 */
  SanchezNewtonRafeson(mind, maxd, c2, c3, c4, d1, dthresh);
/* Next, iterate on DAam.   It is bounded by -min[p2. p.2, (1-p2.)(1-p.2)] and min[(1-p2.)p.2, p2. (1-p.2)] */  
  mind =  p2d * pd2;
  mtemp = (1-p2d) * (1-pd2);
  if (mtemp < mind )
    mind = mtemp;
  mind = -mind; 
  maxd =  (1-p2d) * pd2;
  mtemp = p2d * (1-pd2);
  if (mtemp < maxd )
    maxd = mtemp;
  c1= (1.0-rho1)*(1.0-rho2);
  c2= -(  (1.0-rho1)*( (p1d+rho2*p2d)*pd2+(p2d+rho2*p1d)*pd3 )  - (1.0-rho2)*( (p1d+rho1*p2d)*pd1+(p2d+rho1*p1d)*pd2) )    ;
  c3= -((p1d+rho1*p2d)*(p1d+rho2*p2d)*pd1+(p1d+rho2*p2d)*(p2d+rho1*p1d)*pd2+(p2d+rho1*p1d)*(p2d+rho2*p1d)*pd3-(1.0-rho1)*(1.0-rho2)*p1d*p2d*(1.0-pd2))*pd2;
  c4= -((1.0-rho1)*(p1d+rho2*p2d)*pd1-(1.0-rho2)*(p2d+rho1*p1d)*pd3)*p1d*p2d*pd2*pd2;
/*  Divide all by c1.  */  
  c2=c2/c1;
  c3=c3/c1;
  c4=c4/c1;    /*  now, d^3 + c2 d^2 + c3 d + c4 = 0 */
  SanchezNewtonRafeson(mind, maxd, c2, c3, c4, d2, dthresh);
}

/* 
Safe Newton-Raphson root finder:  p 274 of Num. Rec. C. Yes, it's misspelled. 
*/    
void SanchezNewtonRafeson(double mind, double maxd, double c2, double c3, double c4, double *d1, double dthresh) {
  double mtemp, lv, ld, uv, ud, xl, xh, dxold, mv, mf, dx;
  int i,maxit;
  maxit = 100;
  lv=SanchezCubic(mind,c2,c3,c4);
  ld=SanchezCubicDerivative(mind,c2,c3,c4);
  uv=SanchezCubic(maxd,c2,c3,c4);
  ud=SanchezCubicDerivative(maxd,c2,c3,c4);
  if ( lv*uv < 0.0 )  {
    if ( lv < 0.0 ) {
      xl = mind;
      xh = maxd; 
    } else {
      xh = mind;
      xl = maxd;
    }
    *d1 = 0.5*(maxd+mind);
    dxold=fabs(maxd-mind);
    dx = dxold;
    mv = SanchezCubic(*d1,c2,c3,c4);
    mf = SanchezCubicDerivative(*d1,c2,c3,c4);
    for ( i=1; i<=maxit; i++ ) {
      if ( ( ((*d1-xh)*mf-mv) * ((*d1-xl)*mf-mv) >= 0.0 ) || (fabs(2.0*mv) > fabs(dxold*mf)))  {
        dxold = dx;
        dx=0.5*(xh-xl);
        *d1=xl+dx;
        if (xl == *d1) 
          i = maxit + 1;
      } else {
        dxold=dx;
        dx=mv/mf;
        mtemp=*d1;
        *d1 -= dx;
        if ( mtemp == *d1 )
          i= maxit + 1;
      }
      if ( fabs(dx) < dthresh ) 
        i=maxit + 1;
      else {
        mv = SanchezCubic(*d1,c2,c3,c4);
        mf = SanchezCubicDerivative(*d1,c2,c3,c4);
        if ( mv <0.0)
          xl = *d1;
        else
          xh = *d1;
      }
       
    }
  }
  else
    *d1 = -10.0;
}

/*
Calculate the model form the margins and disequilibria.  Table 2 in the paper. 
d1 = DAAM
d2 = DAaM

*/
void SanchezCalcModelD( anestim *workest, double p1d, double pd1, double pd2,  double d1, double d2 ) {

    workest->PAA = pd1;
    workest->PAa = pd2;
    workest->Paa = 1.0 - pd1 - pd2;
    workest->PM = p1d;
    workest->PAAM = d1 + pd1 * p1d ;
    workest->PAaM = d2 + pd2 * p1d ;
    workest->PaaM = -(d1+d2) + workest->Paa * workest->PM ;
    workest->PAAm = workest->PAA - workest->PAAM;
    workest->PAam = workest->PAa - workest->PAaM;
    workest->Paam = workest->Paa - workest->PaaM;
    workest->PAM = workest->PAAM + 0.5 * workest->PAaM;
    workest->PA = pd1 + 0.5*pd2; 
    workest->rho1 = workest->PAAM * workest->PAam / (workest->PAaM * workest->PAAm ) ;
    workest->rho2 = workest->PAaM * workest->Paam / (workest->PaaM * workest->PAam ) ; 

}
/*
  Header for the output file when doing exact sample sizes.  
*/
void SanchezBanner(FILE *outf) {
  char *tptr;
        tptr = asctime2();

      fprintf(outf,"#\n#\n#\tExact sample sizes needed to detect dependence in 2 x 3 tables.");
      fprintf(outf,"\n#\n#\n#");
      fprintf(outf,"\n#\tAuthors:  M. S. Sanchez, C. J. Basten, A. M. Ferrenberg, M. A. Asmussen and J. Arnold");
      fprintf(outf,"\n#\n#\tSubmitted to Theoretical Population Biology, 2004");
      fprintf(outf,"\n#\n#\tThis particular implementation by Chris Basten (basten@statgen.ncsu.edu).\n#");
      fprintf(outf,"\n#  It is currently %s#\n",tptr);
}
/*

    Show the model.  This will show the cell and margin probabilities as well as the
    input and calculated rho1, rho2 values and the disequilibria.  
*/ 
void SanchezShowModel(FILE *outf, anestim *workest, double d1, double d2 ) {
  double erho1, erho2;
  erho1 = workest->PAAM *workest->PAam /(workest->PAaM *  workest->PAAm);
  erho2 = workest->PAaM *workest->Paam /(workest->PaaM *  workest->PAam);
    fprintf(outf, "\n#\n#  Here is the model.  \n#");
    fprintf(outf,"\n------------------------------------------------");
    fprintf(outf,"\n------------------------------------------------");
    fprintf(outf,"\n                    Factor B ");
    fprintf(outf,"\n          --------------------------");
    fprintf(outf,"\nFactor A      B1       B2       B3   +    Total");
    fprintf(outf,"\n-------------------------------------+----------");
    fprintf(outf,"\n  A1 (M)  %7.4f  %7.4f  %7.4f  |  %7.4f",workest->PAAM,workest->PAaM,workest->PaaM,workest->PM);
    fprintf(outf,"\n  A2 (m)  %7.4f  %7.4f  %7.4f  |  %7.4f",workest->PAAm,workest->PAam,workest->Paam,1.0-workest->PM);
    fprintf(outf,"\n-------------------------------------+----------");
    fprintf(outf,"\nTotal     %7.4f  %7.4f  %7.4f  |   1.0",workest->PAA,workest->PAa,workest->Paa );
    fprintf(outf,"\n------------------------------------------------");
    fprintf(outf,"\n------------------------------------------------");
    fprintf(outf,"\n  rho1 = %7.4f,  rho2 = %7.4f",workest->rho1, workest->rho2);
    fprintf(outf,"\n erho1 = %7.4f, erho2 = %7.4f",erho1, erho2);
    if ( d1 != 0.0 && d2 != 0.0 )
      fprintf(outf,"\n    d1 = %7.4f,    d2 = %7.4f",d1,d2);
    fprintf(outf,"\n------------------------------------------------");
    fprintf(outf,"\n");
    fprintf(outf,"\n  Script readable parameters");
    fprintf(outf,"\n-p1.  %f",workest->PM);
    fprintf(outf,"\n-p.1  %f",workest->PAA);
    fprintf(outf,"\n-p.2  %f",workest->PAa);
    fprintf(outf,"\n-rho1 %f",workest->rho1);
    fprintf(outf,"\n-rho2 %f",workest->rho2);
    fprintf(outf,"\n-d1   %f",d1);
    fprintf(outf,"\n-d2   %f",d2);
    fprintf(outf,"\n#\n#");
 
}
 
double SanchezCubicDerivative(double dv, double c2, double c3, double c4) {
  return( 3.0*dv*dv + 2.0*c2*dv + c3  );
}
double SanchezCubic(double dv, double c2, double c3, double c4) {

  return( dv*dv*dv + c2*dv*dv + c3*dv + c4 );

}
/*
  return log(x!)
*/
double lfct( double x) {

  int t,tm;
  double rt;
  if (x <=  1.0 ) return(0.0); 
  
  tm = (int)(x + 0.5); /* round, do not truncate */
  rt=0.0;
  if (x > 20.0 ) return (lappr_fct(x)) ;
    for(t=2; t<=tm; t++) {
       rt += log((double)t);
    }
  return rt;
}

/*
  Return the log(i!)  approx.
*/
double lappr_fct (double i)
{
  return(  i*log(i) + 0.5*log(2.0*M_PI*i) - i ); 
/*
  return ( ((double)pow(i,i)*(double)exp(-i)*(double)sqrt(2.*M_PI*i)));
  */
}

/*----------------------------------------------------------------------------
  This section allocates memory to hold the information for a single enode
----------------------------------------------------------------------------*/
struct enode *alloc_enode() {
 struct enode *anode;
 anode = (struct enode *) malloc( sizeof(struct enode) ); 
 anode->next = anode->prev = NULL;
 anode->nullprob = anode->condprob = 0.0; 
 return(anode);
}


void dealloc_enodes(struct enode *anode) {
  struct enode *lptr,*tptr;
  int cntr;
  lptr = anode;
  while ( lptr->next != NULL )
    lptr = lptr->next;
  cntr=0;
  while ( lptr != NULL ) {
      cntr++;
	  tptr = lptr->prev;
      free((char *) lptr);
	  lptr = tptr;
  }
  tptr = NULL;
}

/*
  Set all nodes to have nullprob and condprob 0.0.
  return the rightmost node. 
*/
exactnode  *clean_enodes(exactnode *anode) {
  exactnode *current,*rnode;
  
  for ( rnode=anode; rnode->next != NULL; rnode=rnode->next )
    rnode->nullprob=rnode->condprob=0.0;
  
  rnode->nullprob=rnode->condprob=0.0;
  
  for ( current=anode; current != NULL; current=current->prev ) 
    current->nullprob=current->condprob=0.0;
  
  return(rnode);
}

/*
   Assume that first->nullprob <= size.   This new node may push it over the top.
*/
void push_enode(exactnode *first, double nullprob, double condprob) {

  exactnode *cnode;
  if ( first->prev == NULL ) /* if none available, create one*/
    cnode = alloc_enode();
  else {  /* pull the nearest one to the left. */
    cnode = first->prev;
    first->prev = cnode->prev;
    if ( first->prev != NULL )
      first->prev->next = first;
  }
  /* insert right after first */
    cnode->nullprob=nullprob;
    cnode->condprob=condprob;
    cnode->prev = first;
    cnode->next = first->next;
    if ( cnode->next != NULL )
      cnode->next->prev = cnode;
    first->next = cnode;
    first->nullprob += nullprob;
    first->condprob += condprob;
}

/*
   Assume that first->nullprob > size, and that deleting the largest will 
   make it less than or equal to size.
*/
void pull_enode(exactnode *first ) {

  exactnode *cnode,*maxnode;
  double maxvalue;
  maxvalue = 0.0;
  maxnode=NULL;

  for ( cnode=first->next; cnode != NULL; cnode=cnode->next) 
    if ( maxvalue < cnode->nullprob ) {
      maxvalue = cnode->nullprob;
      maxnode = cnode;
    }
    
  if ( maxnode != NULL ) {/*  we found a node */
    first->nullprob = first->nullprob - maxnode->nullprob;
    first->condprob = first->condprob - maxnode->condprob;
    maxnode->prev->next = maxnode->next;
    if (maxnode->next != NULL )
      maxnode->next->prev = maxnode->prev;
    maxnode->next = first;
    maxnode->prev = first->prev;
    if (maxnode->prev != NULL )
      maxnode->prev->next = maxnode;
    first->prev = maxnode;
    maxnode->condprob = maxnode->nullprob = 0;
  
  }
}

/*
  This function drives the Sanchez analysis.   This mainly calculates powers using the
  exact test.   It will calculate the power for N = Nl, Nu incremented by Ni.  
     

*/
void  DriveSanchez(anestim *counts, char *outfile, char *errfile ) {
  FILE *outf;
  anestim workest;
  double d1, d2, thepower,pd1,p1d,pd2,rho1,rho2,size;
  exactnode *first; 
  char *tptr;
  PrepNormFacts(counts[1].nAaM); 


  outf = fileopen(outfile, "a");
  if ( outf == NULL ) 
    outf = fileopen(outfile, "w");
  if (outf != NULL) {
      SanchezBanner(outf);
      first = alloc_enode();
	  pd1 = counts[1].PAA;
	  pd2 = counts[1].PAa;
	  p1d = counts[1].PM;
	  rho1 = counts[1].rho1;
	  rho2 = counts[1].rho2;
	  size = counts[1].PA;
	  workest = *(counts+0);
	  if ( counts[1].nAam == 0 )
	    SanchezSolveDiseq( rho1, rho2, p1d, pd1, pd2, &d1, &d2 );
	  else {
	    d1 = counts[1].PAAM;
	    d2 = counts[1].PAaM;
	  }
	  SanchezCalcModelD( &workest, p1d, pd1, pd2, -d1, -d2 );
      SanchezShowModel(outf,&workest,-d1,-d2);	  

      /*        For each sample size, estimate the power.  */
      fprintf(outf,"\n-start powercurve");
      fileclose(outfile, outf);
      
      thepower = 0.0;

      for ( workest.samplesize = counts[1].nAAM; workest.samplesize <= counts[1].nAaM  && thepower < counts[1].Paa ; workest.samplesize  += counts[1].naaM ) 
        first = SanchezPower(outfile, &workest, first, size, rho1, rho2, &thepower);

      

      outf = fileopen(outfile, "a");
      fprintf(outf,"\n-end powercurve");

      tptr=asctime2();
      fprintf(outf,"\n#\n#  It is currently %s#\n",tptr);
      fileclose(outfile, outf);
  }
  UnPrepNormFacts(counts[1].nAaM);	

}   


/*
     This mainly calculates powers using asymptotic results.   
  It was mainly used for debugging and to create power curves.
 There will be an indicator
  in counts[1].nAam that says whether to calculate D1 and D2 from the margins and
  rho1, rho2 or to use the values of D1 and D2 directly.   */
void  DrivePower(anestim *counts, char *outfile, char *errfile ) {
  FILE *outf;
  anestim workest;
  double d1, d2, thepower,pd1,p1d,pd2,rho1,rho2,size;
  exactnode *first; 
  char *tptr;
  PrepNormFacts(counts[1].nAaM); 


  outf = fileopen(outfile, "a");
  if ( outf == NULL ) 
    outf = fileopen(outfile, "w");
  if (outf != NULL) {
      SanchezBanner(outf);
      first = alloc_enode();
	  pd1 = counts[1].PAA;
	  pd2 = counts[1].PAa;
	  p1d = counts[1].PM;
	  rho1 = counts[1].rho1;
	  rho2 = counts[1].rho2;
	  size = counts[1].PA;
	  workest = *(counts+0);
	  if ( counts[1].nAam == 0 )
	    SanchezSolveDiseq( rho1, rho2, p1d, pd1, pd2, &d1, &d2 );
	  else {
	    d1 = counts[1].PAAM;
	    d2 = counts[1].PAaM;
	  }
	  SanchezCalcModelD( &workest, p1d, pd1, pd2, -d1, -d2 );
      SanchezShowModel(outf,&workest,-d1,-d2);	  
/*  This was in for a test of AsympPower.         
	  for ( workest.samplesize = 50; workest.samplesize <= 1000; workest.samplesize  += 10 ) {
        apower = AsympPower( &workest, 1.96);   
        fprintf(outf,"\n %4d  %6.4f ",workest.samplesize,  apower);   
      } 
*/
      /*        For each sample size, estimate the power.  */
      fprintf(outf,"\n-start powercurve");
      
      for ( workest.samplesize = counts[1].nAAM; workest.samplesize <= counts[1].nAaM; workest.samplesize  += counts[1].naaM ) {
        thepower = AsympPower(&workest,1.96);
        fprintf(outf, "\n  -n  %d   -p  %f ", workest.samplesize, thepower);
      }
      fileclose(outfile, outf);


 
      outf = fileopen(outfile, "a");
      fprintf(outf,"\n-end powercurve");

      tptr=asctime2();
      fprintf(outf,"\n#\n#  It is currently %s#\n",tptr);
      fileclose(outfile, outf);
  }
  UnPrepNormFacts(counts[1].nAaM);	

}   
/*
  This function drives the Sanchez analysis.   This calculate the sample size for a single
  set of parameters contained within the counts[1] structure.   There will be an indicator
  in counts[1].nAam that says whether to calculate D1 and D2 from the margins and
  rho1, rho2 or to use the values of D1 and D2 directly.   

*/
void  CalcSanchezN(anestim *counts, char *outfile, char *errfile ) {
  FILE *outf;
  int k,asympn;
  double  d1, d2,pd1,p1d,pd2,rho1,rho2,size,beta;
  char *tptr;
  long start,finish;
  

  PrepNormFacts((int) MAXEXACTSS); 


  outf = fileopen(outfile, "a");
  if ( outf == NULL ) 
    outf = fileopen(outfile, "w");
  if (outf != NULL) {
      SanchezBanner(outf);
	  pd1 = counts[1].PAA;
	  pd2 = counts[1].PAa;
	  p1d = counts[1].PM;
	  rho1 = counts[1].rho1;
	  rho2 = counts[1].rho2;
	  size = counts[1].PA;
	  beta = counts[1].Paa;

	  if ( counts[1].nAam == 0 )
	    SanchezSolveDiseq( rho1, rho2, p1d, pd1, pd2, &d1, &d2 );
	  else {
	    d1 = counts[1].PAAM;
	    d2 = counts[1].PAaM;
	  }

	  SanchezCalcModelD( &counts[0], p1d, pd1, pd2, -d1, -d2 );
      SanchezShowModel(outf,&counts[0],-d1,-d2);	 
      fileclose(outfile, outf);
      asympn = AsympSampleSize( counts[0], size, beta );
      start = marktime();
      k = SanchezSampleSize( &counts[0], &counts[2], size,beta,outfile);
      finish = marktime();

      outf = fileopen(outfile, "a");

      fprintf(outf, "\n-Nexact %d\n-Nasymp %d\n-alpha %f\n-beta %f\n-power %f\n# ",counts[2].eN,counts[2].samplesize, counts[2].esize, 1.0-counts[2].epower,counts[2].epower);
      fprintf(outf, "\n  Minimum sample size from exact test : %d and asymptotic test: %d",k,asympn);
      fprintf(outf, "\n  Minimum sample size from exact test : %d and asymptotic test: %d   (check)",counts[2].eN,counts[2].samplesize);
	  fprintf(outf, "\n  P(type I) = %6.4f,  P(type II) = %6.4f, Power = %6.4f    as given originally",size,beta, 1.0-beta);	
	  fprintf(outf, "\n  P(type I) = %6.4f,  P(type II) = %6.4f, Power = %6.4f    as achieved\n#", counts[2].esize,1.0- counts[2].epower, counts[2].epower);	
      fprintf(outf, "\n-time %d", finish-start);
      tptr=asctime2();
      fprintf(outf,"\n#\n#  It is currently %s#\n",tptr);
      fileclose(outfile, outf);
  }
  UnPrepNormFacts((int) MAXEXACTSS);	

}   

/*
  This function is a check on table 5.     It checks whether the D11, D12 from the original table
  yield rho1=rho2=0.5 (no).   Also, recalculated D11, D12 with rho1=rho2=0.5 and checks if they
  yield a model with rho1=rho2=0.5 (yes).   
  
  
  table[i] is the i th row of Table 5.  Columns are
  
  1.  p1.   PM
  2.  p.1   PAA
  3.  p.2   PAa                           The first six are given
  4.  D11   DAAM
  5.  D12   DAaM
  6.  N     exact test sample size
  
  7.  Na    asymptotic sample size         Calc. based on the given disequlibria
  8.  exact power                          
  9.  rho1  from the given D11, D12
  10. rho2  from the given D11, D12
  
  11. PAAM
  12. PAaM
  13. PaaM
  14. PAAm
  15. PAam
  16. Paam
  17. DaaM
  18. p.3    Paa
  19. p2.    Pm  
  20. asymptotic power

  21. est D11 from rho's estimated
  22. est D12 from rho's estimated
  
  23. est D11 from rho=0.5                 Set rho1=rho2=0.5 and recalc D11, D12, powers for given SS, 
  24. est D12 from rho=0.5                 Asymp. sample size.   
  25. estimate rho1 from 23, 24 and margins model :   check that 23, 24 really give the rho's input
  26. estimate rho2 from 23, 24 and margins model
  27. asymp. sample size for rho1=rho2=0.5, beta = 0.5
  28. Exact sample size using Basten params and estimators (for D1, D2)
  29. Exact power for 28
  30. Exact size achieved for 28
  
  31. PAAM  when rho1=rho2=0.5 and we estimate D11, D12 here
  32. PAaM
  33. PaaM
  34. PAAm
  35. PAam
  36. Paam
  37. PM
  38. PAA
  39. PAa
  40. Time to calc N for exact test with D1, D2 of 44, 45 below.
  
  41.  p1.   PM
  42.  p.1   PAA
  43.  p.2   PAa                          
  44.  D1  with rho1 = rho2 = 0.5 solved here
  45.  D2   "   "
  46.  Exact sample size using 41-45
  47.  Exact power for 46
  48.  Exact size achieved for 46
  49.  Asymptotic sample size for disequilibria in 44, 45
  50.  Asymptotic power for 49.
  
*/
void  SanchezTable5(anestim *counts, char *outfile, char *errfile ) {
  FILE *outf;
  anestim workest;
  int i,j;
  double  d1, d2, thepower,pd1,p1d,pd2,rho1,rho2,size,beta;
  exactnode *first; 
  double  **table, zalpha2;
  char *tptr;
  long start,finish;
  PrepNormFacts((int) MAXEXACTSS ); 
  

  table=dmatrix(1,20,1,50);
  i=1;
  table[i][1] = 0.1; table[i][2]=0.2; table[i][3]=0.2; table[i][4]=-0.0146; table[i][5]=-0.0025; table[i][6]=282; i++;
  table[i][1] = 0.1; table[i][2]=0.2; table[i][3]=0.3; table[i][4]=-0.0133; table[i][5]=-0.0024; table[i][6]=320; i++;
  table[i][1] = 0.1; table[i][2]=0.3; table[i][3]=0.4; table[i][4]=-0.0142; table[i][5]= 0.0021; table[i][6]=390; i++;
  
  table[i][1] = 0.3; table[i][2]=0.1; table[i][3]=0.1; table[i][4]=-0.0196; table[i][5]=-0.0074; table[i][6]=196; i++;
  table[i][1] = 0.3; table[i][2]=0.1; table[i][3]=0.3; table[i][4]=-0.0171; table[i][5]=-0.0151; table[i][6]=192; i++;
  table[i][1] = 0.3; table[i][2]=0.1; table[i][3]=0.5; table[i][4]=-0.0146; table[i][5]=-0.0140; table[i][6]=225; i++;
  table[i][1] = 0.3; table[i][2]=0.1; table[i][3]=0.7; table[i][4]=-0.0123; table[i][5]=-0.0053; table[i][6]=336; i++;
  table[i][1] = 0.3; table[i][2]=0.2; table[i][3]=0.2; table[i][4]=-0.0312; table[i][5]=-0.0074; table[i][6]=138; i++;
  table[i][1] = 0.3; table[i][2]=0.2; table[i][3]=0.3; table[i][4]=-0.0287; table[i][5]=-0.0078; table[i][6]=149; i++;
  table[i][1] = 0.3; table[i][2]=0.2; table[i][3]=0.4; table[i][4]=-0.0264; table[i][5]=-0.0062; table[i][6]=168; i++;
  table[i][1] = 0.3; table[i][2]=0.3; table[i][3]=0.4; table[i][4]=-0.0320; table[i][5]= 0.0025; table[i][6]=166; i++;
  
  table[i][1] = 0.5; table[i][2]=0.1; table[i][3]=0.1; table[i][4]=-0.0202; table[i][5]=-0.0088; table[i][6]=211; i++;
  table[i][1] = 0.5; table[i][2]=0.1; table[i][3]=0.3; table[i][4]=-0.0180; table[i][5]=-0.0189; table[i][6]=187; i++;
  table[i][1] = 0.5; table[i][2]=0.1; table[i][3]=0.5; table[i][4]=-0.0157; table[i][5]=-0.0189; table[i][6]=202; i++;
  table[i][1] = 0.5; table[i][2]=0.1; table[i][3]=0.7; table[i][4]=-0.0134; table[i][5]=-0.0087; table[i][6]=286; i++;
  table[i][1] = 0.5; table[i][2]=0.2; table[i][3]=0.2; table[i][4]=-0.0338; table[i][5]=-0.0102; table[i][6]=132; i++;
  table[i][1] = 0.5; table[i][2]=0.2; table[i][3]=0.3; table[i][4]=-0.0315; table[i][5]=-0.0114; table[i][6]=137; i++;
  table[i][1] = 0.5; table[i][2]=0.2; table[i][3]=0.4; table[i][4]=-0.0292; table[i][5]=-0.0101; table[i][6]=149; i++;
  table[i][1] = 0.5; table[i][2]=0.3; table[i][3]=0.4; table[i][4]=-0.0367; table[i][5]= 0.0000; table[i][6]=139; i++;
  table[i][1] = 0.1; table[i][2]=0.1; table[i][3]=0.1; table[i][4]=-0.0099; table[i][5]=-0.0031; table[i][6]=340;

  
  workest = *(counts+0);
  size = counts[1].PA ;
  beta = counts[1].Paa;
  zalpha2 = Phi(1.0-size/2.0);

  for ( i=1; i<=20; i++ ) {
    table[i][17] = -( table[i][4] + table[i][5] );
    table[i][18] = 1.0 - table[i][2] - table[i][3];
    table[i][19] = 1.0 - table[i][1];
    table[i][11] = table[i][4] + table[i][1] * table[i][2] ;
    table[i][12] = table[i][5] + table[i][1] * table[i][3] ;
    table[i][13] = table[i][17] + table[i][1] * table[i][18] ;
    table[i][14] = table[i][2] - table[i][11];
    table[i][15] = table[i][3] - table[i][12];
    table[i][16] = table[i][18] - table[i][13];
    table[i][9]  = table[i][11] * table[i][15] / ( table[i][12] * table[i][14] ) ;
    table[i][10] = table[i][12] * table[i][16] / ( table[i][13] * table[i][15] ) ;
    
	table[i][42] = table[i][2];
	table[i][43] = table[i][3];
	table[i][41] = table[i][1];
    
    workest.rho1=workest.rho2=0.5;
    SanchezSolveDiseq( workest.rho1, workest.rho2, table[i][41],table[i][42],table[i][43], &d1, &d2 );
    table[i][44] = -d1;
    table[i][45] = -d2;
	SanchezCalcModelD( &workest, table[i][41],table[i][42],table[i][43], table[i][44], table[i][45] );
	table[i][49] =  AsympSampleSize( workest,size, 1.0-beta ) ;
	workest.samplesize = table[i][49];
    table[i][50] =  AsympPower( &workest, zalpha2);
    if ( debug == 2 )
      printf("\n%3.2f %3.2f %3.2f %7.4f %7.4f %5.0f %5.4f",table[i][41],table[i][42],table[i][43], table[i][44], table[i][45], table[i][49], table[i][50]);
  }
  if ( debug < 2 ) {
	  outf = fileopen(outfile, "a");
	  if ( outf == NULL ) 
	    outf = fileopen(outfile, "w");
	  if (outf != NULL) {
	      SanchezBanner(outf);
	      fprintf(outf,"\n\n   Size:  %f,  beta:  %f" , size, beta);
	      fileclose(outfile, outf);
	      first = alloc_enode();
	      for ( i=1; i<=20; i++ ) {
		    pd1 = table[i][2];
		    pd2 = table[i][3];
		    p1d = table[i][1];
		    rho1 = table[i][9];
		    rho2 = table[i][10];
		    workest.samplesize = (int)  table[i][6];
		    SanchezCalcModelD( &workest, p1d, pd1, pd2, table[i][4], table[i][5] );
		    table[i][7] =  AsympSampleSize( workest,size, 1.0-beta) ;
		    if ( debug == 0 ) {
	          first = SanchezPower(outfile, &workest, first, size, rho1, rho2, &thepower);
	          table[i][8] = thepower; 
	        }
	        table[i][20] =   AsympPower( &workest, zalpha2);

		    SanchezSolveDiseq( rho1, rho2, p1d, pd1, pd2, &d1, &d2 );
		    table[i][21] = -d1;
		    table[i][22] = -d2;
		    SanchezSolveDiseq( 0.5, 0.5, p1d, pd1, pd2, &d1, &d2 );
		    table[i][23] = -d1;
		    table[i][24] = -d2;        
		    SanchezCalcModelD( &workest, p1d, pd1, pd2, table[i][23], table[i][24] );
	        table[i][25] = workest.PAAM * workest.PAam / ( workest.PAAm * workest.PAaM);
	        table[i][26] = workest.PAaM * workest.Paam / ( workest.PAam * workest.PaaM);


		    workest.samplesize =   AsympSampleSize( workest,size, 1.0-beta) ;
	        table[i][27] = (double) workest.samplesize;
	        if ( debug == 0 ) {
	          start = marktime();
	          table[i][46] = table[i][28] = (double) SanchezSampleSize( &workest, &counts[1], size,beta,outfile);
	          table[i][47] = table[i][29] = counts[1].epower;
	          table[i][48] = table[i][30] = counts[1].esize;
              finish = marktime();
              table[i][40] = 0.0001 * (double) (finish-start) ;
	        }

	        
	        
	        table[i][31] = workest.PAAM; 
	        table[i][32] = workest.PAaM;
	        table[i][33] = workest.PaaM;
	        table[i][34] = workest.PAAm; 
	        table[i][35] = workest.PAam;
	        table[i][36] = workest.Paam;
	        table[i][37] = workest.PAAM + workest.PAaM + workest.PaaM;
	        table[i][38] = workest.PAAM + workest.PAAm;
	        table[i][39] = workest.PAaM + workest.PAam;
	      }
	      outf = fileopen(outfile, "a");
	 
	 
	/* ************ *************** */ 
	      fprintf(outf,"\n   The first six columns are from your Table 5.  Na is the sample size corresponding to");
	      fprintf(outf,"\n   N using the asymptotic results.  EP is the power achieved by N given the first 5 ");
	      fprintf(outf,"\n   columns.   rho1 and rho2 are estimated from the first five columns.  Note that they");
	      fprintf(outf,"\n   are not 0.5.");
	 
	      fprintf(outf,"\n-------------------------------------------------------------------------------------------");
	      fprintf(outf,"\n-------------------------------------------------------------------------------------------");
	      fprintf(outf,"\n      p1.      p.1      p.2      D11      D12        N       Na       EP     rho1     rho2");
	      fprintf(outf,"\n-------------------------------------------------------------------------------------------");



	      for ( i=1; i<=20; i++ ) {
	      
	        fprintf(outf,"\n");
	        for (j=1; j<=10; j++ )
	          fprintf(outf," %8.4f",table[i][j]);
	      }
	      fprintf(outf,"\n-------------------------------------------------------------------------------------------");
	      fprintf(outf,"\n-------------------------------------------------------------------------------------------");

	/* ************ *************** */ 
	fprintf(outf,"\nThe first five columns above are used to calculate the first nine columns here.  AP is the ");
	fprintf(outf,"\npower of N above from the asymptotic equations.  ");

	      fprintf(outf,"\n-------------------------------------------------------------------------------------------");
	      fprintf(outf,"\n-------------------------------------------------------------------------------------------");
	      fprintf(outf,"\n     PAAM     PAaM     PaaM     PAAm     PAam     Paam     DaaM      p.3      p2.       AP");
	      fprintf(outf,"\n-------------------------------------------------------------------------------------------");



	      for ( i=1; i<=20; i++ ) {
	      
	        fprintf(outf,"\n");
	        for (j=11; j<=20; j++ )
	          fprintf(outf," %8.4f",table[i][j]);
	      }
	      fprintf(outf,"\n-------------------------------------------------------------------------------------------");
	      fprintf(outf,"\n-------------------------------------------------------------------------------------------");
	 
	 
	/* ************ *************** */ 
	      fprintf(outf,"\nThe sD11 and sD12 are calculated from the first three columns in the first table");
	      fprintf(outf,"\nalong with the last two columns (rho1, rho2) and they recover the 6th and 7th columns");
	      fprintf(outf,"\nin the first table.   The bD11 and bD12 are the disequilibria calculated from rho1 = rho2 =");
	      fprintf(outf,"\n0.5.   bNa is the asymptotic sample size and  using bD11 and bD12  and has asymptotic power given");
	      fprintf(outf,"\nin as(1-b) exact power ex(1-b).  bN is the exact sample size for power %4.2f. The rho's are ", 1.0-beta );
	      fprintf(outf,"\nrecalculated from themargins and the bD11, bD12 values.   ");

	      fprintf(outf,"\n-------------------------------------------------------------------------------------------");
	      fprintf(outf,"\n-------------------------------------------------------------------------------------------");
	      fprintf(outf,"\n     sD11     sD12     bD11     bD12      rho1     rho2     bNa       bN    epower    esize");
	      fprintf(outf,"\n-------------------------------------------------------------------------------------------");



	      for ( i=1; i<=20; i++ ) {
	      
	        fprintf(outf,"\n");
	        for (j=21; j<=30; j++ )
	          fprintf(outf," %8.4f",table[i][j]);
	      }
	      fprintf(outf,"\n-------------------------------------------------------------------------------------------");
	      fprintf(outf,"\n-------------------------------------------------------------------------------------------");
	 
	 /* ************ *************** */ 
	fprintf(outf,"\n These are the six classes estimated with rho1=rho2=0.5 and the \nfirst three columns of the first table.");
	      fprintf(outf,"\n-------------------------------------------------------------------------------------------");
	      fprintf(outf,"\n-------------------------------------------------------------------------------------------");
	      fprintf(outf,"\n     PAAM     PAaM     PaaM     PAAm    PAam     Paam       PM      PAA      PAa      t E-4");
	      fprintf(outf,"\n-------------------------------------------------------------------------------------------");



	      for ( i=1; i<=20; i++ ) {
	      
	        fprintf(outf,"\n");
	        for (j=31; j<=40; j++ )
	          fprintf(outf," %8.4f",table[i][j]);
	      }
	      fprintf(outf,"\n-------------------------------------------------------------------------------------------");
	      fprintf(outf,"\n-------------------------------------------------------------------------------------------");
	 
	 /* ************ *************** */ 
	fprintf(outf,"\n Final table, with D11, D12 estimated assuming rho1=rho2=0.5.  ");
	      fprintf(outf,"\n-------------------------------------------------------------------------------------------");
	      fprintf(outf,"\n-------------------------------------------------------------------------------------------");
	      fprintf(outf,"\n      p1.      p.1      p.2      D11      D12        N     epower    esize     Na   apower ");
	      fprintf(outf,"\n-------------------------------------------------------------------------------------------");



	      for ( i=1; i<=20; i++ ) {
	      
	        fprintf(outf,"\n");
	        for (j=41; j<=50; j++ )
	          fprintf(outf," %8.4f",table[i][j]);
	      }
	      fprintf(outf,"\n-------------------------------------------------------------------------------------------");
	      fprintf(outf,"\n-------------------------------------------------------------------------------------------");
	 
	      tptr=asctime2();
	      fprintf(outf,"\n#\n#  It is currently %s#\n",tptr);


	  }	
  }
  free_dmatrix( table,1,20,1,50);
  UnPrepNormFacts((int) MAXEXACTSS); 

}   

double qgaus(double lb, double ub) {


  int j;
  double xr, xm, dx, s;
  static double x[] = {0.0, 0.1488743389, 0.4333953941, 0.6794095682, 0.8650633666, 0.9739065285};
  static double w[] = {0.0, 0.2955242247, 0.2692667193, 0.2190863625, 0.1494513491, 0.0666713443};

  xm = 0.5*(ub+lb);
  xr = 0.5*(ub-lb);
  s=0.0;
  for ( j=1; j<=5; j++ ){
    dx = xr*x[j];
    s += w[j]*( standardnormal(xm+dx)  + standardnormal(xm-dx) );
  
  }
  return  s *= xr  ;

}


double standardnormal( double z ) {
 double v;
 return  v = exp( -0.5*z*z)/ sqrt( 2.0*M_PI ) ;

}


#define NNORMAL 4000

/*
Two functions to create the table of normals deviates and
log of factorials.  
*/
void UnPrepNormFacts(int nfacts) {
  free_dvector(factorials, 0, nfacts);
  free_dmatrix(normals,0,1,0,(int) NNORMAL);
}

void PrepNormFacts(int nfacts) {
int i;
double lb,ub,thisvalue,lbp,ubp;
/* Create a table to hold normal deviates and probabilities.  */  
    factorials = dvector(0,nfacts);
    normals = dmatrix( 0, 1, 0, (int) NNORMAL );

  
  normals[0][0] = 0.0;
  normals[1][0] = 0.5;
  
  for ( i=1; i<=(int) NNORMAL; i++ ) {
    normals[0][i] = normals[0][i-1] + 0.001;
    thisvalue = qgaus(normals[0][i-1],normals[0][i]);
    normals[1][i] = normals[1][i-1] + thisvalue;
  }
  factorials[0]=factorials[1]=0.0;
  for ( i=2; i<=nfacts; i++ )
    factorials[i]=factorials[i-1]+log((double)i);
    
    
    
/*  Debugging code.
  for ( i=1; i<=100; i++ )  {
    thisvalue = (double) i * 0.03; 
    
    lbp = PhiProb( -thisvalue );
    ubp = PhiProb( thisvalue );   
    lb = Phi( lbp);
    ub = Phi( ubp);
    printf( "\n  %7.3f %7.3f %7.3f %7.3f %7.3f", thisvalue, lbp, lb, ubp, ub);
  }  
*/
}

/*
        Phi and PhiProb yield the standard normal values and probabilities.
        if Z ~ N(0,1), then   
        
          p = P( Z <= z ) = Phi(z)
          
          Phi takes a probability p and gives back z
          PhiProb takes a z and gives back a probability
*/
double Phi(double pr) {
  int i;
  double zp;
  
  if ( pr == 0.5 )
    return(0.0);

  if ( pr < 0.5 ) 
    zp = pr + 0.5; 
  else 
    zp = pr;
    
  for ( i=1; i<=(int) NNORMAL; i++ ) 
    if ( zp > normals[1][i-1] && zp <= normals[1][i]) {
      if ( pr < 0.5 )
        return( -normals[0][i]);
      else
        return( normals[0][i]) ;
     }  
}

double PhiProb(double z) {
  int i;
  double zp;
  
  if ( z == 0.0 )
    return(0.5);
  else if ( z > normals[0][(int) NNORMAL] )
    return( 1.0 );
  else  if ( -z > normals[0][(int) NNORMAL] )
    return( 0.0 );

  if ( z < 0.0 ) 
    zp = -z; 
  else
    zp = z; 
  

    
  for ( i=1; i<=(int) NNORMAL; i++ ) 
    if ( zp > normals[0][i-1] && zp <= normals[0][i] ) {
      if ( z < 0.0 )
        return( 1.0 - normals[1][i] );
      else
        return( normals[1][i]) ;
     }  
}
#undef NNORMAL
