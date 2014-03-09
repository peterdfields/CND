#include "Dmain.h"

/* 
    USAGE:
  
  cnib -o outputfile -e errorfile -i inputfile 

If you use an option, be sure to give a value (eg cnisim -e -o will give nonsense)

*/

#ifdef __MACOS__
#include <console.h>
#endif

int debug;
double margthresh;

/*
  Notes on the program.
*/

int     main(argc, argv)
  int     argc;
  char   *argv[];
{
  char    *inputf,*ipfptr;
  char    *outfile, *errfile,**nnames,**cnames;
  int     errflag,ii,jj,tex_flag,Nl,Nu,Ni;
  int    numsets,Sanchez,dswitch;
  double tmp,pd1,pd2,p1d,r1,r2,size,beta,pu,d1,d2,dthresh;
  long seed;
  anestim  *counts,*work;
  FILE *outf;
  
/*
  Here is some code needed to get a console window on the Macintosh. 
  I used to use Think C but have switched to Code Warrior.   
*/
#ifdef __MACOS__
 /* make environ safe, and give the console window a new
    title */
  unsigned char abuffer[25] = "\p CNDd \0";
/*  console_options.title = abuffer;  Not used with CodeWarrior*/
/* Simulate the UNIX shell with ccommand, this should also
    be the place to specify redirection of input and
    output. */
  argc = ccommand(&argv);
#endif
  debug = 0;
  margthresh = 10E-20; 
  dswitch = Sanchez = 0;
  d1=d2=0.0;
  Nl=50;
  Nu=100;
  Ni=25;
  beta = 0.5;
  pd1 = pd2 = p1d = 0.1;
  r1 = r2 = 0.5;
  pu = 0.99;
  size = 0.05;
  counts = NULL;
  cnames = nnames = NULL;
  seed = -get_a_seed();
  tmp = ranf(seed);
  tex_flag = 0;
  inputf = cvector(0,MAXNAME);
  outfile = cvector(0,MAXNAME);
  errfile = cvector(0,MAXNAME);
  for ( ii = 0 ; ii <= MAXNAME ; ii++ ) *(outfile+ii) = *(errfile+ii) = *(inputf+ii) = '\0';
  for ( ii = 1 ; ii < argc ; ii++ ) 
    if ( *(*(argv+ii)+0) == '-' )
      switch ( *(*(argv+ii)+1) ) {
        case 'h' :
          printf("\nUSAGE:\n");
          printf("\n   CNDd [-o outfile] [-e errorfile] [-i inputfile] -t\n");
          printf("\nDEFAULTS:\n");
          printf("\n\t [-o CNDd.out]");
          printf("\n\t [-e CNDd.err]");
          printf("\n\t [-i CNDd.dat]");
          printf("\n\t [-t ] Means output data analysis in LaTeX2e format ");
          printf("\n\t [-a]  Set alpha = p(Type I error) for exact sample size (%4.2f)",size);
          printf("\n\t [-b]  Set beta = p(Type II error) for same (%3.1f)",beta);
          printf("\n\t [-S ] Set the Sanchez flag for exact sample sizes [0,1,2,3,4,5]");
          
          printf("\n-S with greater than 0 may require a long time. ");
          printf("\nThe following options are for theoretical calculations.");
          printf("\n\t [-p1] Set p.1 (%3.1f)",pd1);
          printf("\n\t [-p2] Set p.2 (%3.1f)",pd2);
          printf("\n\t [-pM] Set p1. (%3.1f)",p1d);
          printf("\n\t [-r1] Set rho1 (%3.1f)",r1);
          printf("\n\t [-r2] Set rho2 (%3.1f)",r2);
          printf("\n\t [-d1] Set D1 (%3.1f) and ignore rho1",d1);
          printf("\n\t [-d2] Set D2 (%3.1f) and ignore rho2",d2);
          printf("\n\t [-pu] Set upper bound for power in power curve",pu);
          printf("\n\t [-Nl] Set lower bound for sample size (%2d)",Nl);
          printf("\n\t [-Nu] Set upper bound for sample size (%3d)",Nu);
          printf("\n\t [-Ni] Set increment on sample size (%2d)",Ni);
          printf("\n\n\tNote that -p1, -p2 and -pM require real valued arguments in (0,1) ");
          printf("\n\twith 0 < p1+p2 < 1.  -Nl, -Nu and -Ni set the lower bound, upper bound ");
          printf("\n\tand increment on the sample size.  -S with value > 1 mean ignore data files. ");
          printf("\n\tThe given values are the defaults. ");
          printf("\n\n");
          exit(0);
        case 'e' :
          strcpy(errfile, argv[ii+1]);  break;
        case 'o' : 
          strcpy(outfile, argv[ii+1]);  break;
        case 't' : 
          tex_flag = 1;
          break;
        case 'i' : 
          strcpy(inputf, argv[ii+1]); 
          break;
        case 'd' : 
          dswitch = 1;
          if ( argv[ii][2] == '1' )
            d1 = atof( argv[ii+1] );
          else if ( argv[ii][2] == '2' )
            d2 = atof( argv[ii+1] );
          else
            strcpy(inputf, argv[ii+1]); 
          break;
        case 'S' : 
          Sanchez = atoi( argv[ii+1] ) ;
          break;
        case 'D' : 
          debug = atoi( argv[ii+1] ) ;
          break;
        case 'p' : 
          if ( argv[ii][2] == '1' )
            pd1 = atof( argv[ii+1] );
          else if ( argv[ii][2] == '2' )
            pd2 = atof( argv[ii+1] );
          else if ( argv[ii][2] == 'M' )
            p1d = atof( argv[ii+1 ] );
          else
            pu = atof( argv[ii+1 ] );
          break;
        case 'N' : 
          if ( argv[ii][2] == 'l' )
            Nl = atoi( argv[ii+1] );
          else if ( argv[ii][2] == 'u' )
            Nu = atoi( argv[ii+1] );
          else
            Ni = atoi( argv[ii+1 ] );
          break;
        case 'r' : 
          if ( argv[ii][2] == '1' )
            r1 = atof( argv[ii+1] );
          else if ( argv[ii][2] == '2' )
            r2 = atof( argv[ii+1] );
          break;
        case 'b' :            
          beta = atof( argv[ii+1] );
          break;
        case 'a' :            
          size = atof( argv[ii+1] );
          break;
        case 'T' :            
          dthresh = (double) atoi( argv[ii+1] );
          if ( dthresh > 0.0 )
            margthresh = pow( 10.0, -dthresh );
          else
            margthresh = 0.0; 
          break;
       }
  if (*(outfile+0) == '\0' )
    strcpy(outfile, "CNDd.out");
  if (*(errfile+0) == '\0' )
    strcpy(errfile, "CNDd.err");
  if (*(inputf+0) == '\0' )
    strcpy(inputf,"CNDd.dat");
  ipfptr = inputf;
  if ( Sanchez == 0 || Sanchez == 1 || Sanchez == 6) /*  We are doing data */
    numsets =  get_the_data(counts,ipfptr,errfile,cnames,nnames,0);
  else
    numsets = 2;
  counts = anestvector(0,numsets);
  for ( ii = 1 ; ii <= numsets ; ii++ ) 
    clear_anestim((counts+ii));

/*  Check parameters and reset if necessary. */    
  if ( p1d < 0.0 || p1d > 1.0 )
    p1d = 0.1;  
  if ( pd1 < 0.0 || pd1 > 1.0 )
    pd1 = 0.1;  
  if ( pd2 < 0.0 || pd2 > 1.0 )
    pd2 = 0.1;  
  if ( pu < 0.0 || pu > 1.0 )
    pu = 0.99;  
  if ( beta < 0.0 || beta > 1.0 )
    beta = 0.5;  
  if ( size < 0.0 || size > 1.0 )
    size = 0.05;  
  if ( r1 <= 0.0 )
    r1 = 0.5;
  if ( r2 <= 0.0 )
    r2 = 0.5;
  if ( Nl < 1 )
    Nl = 50;
  if (Ni < 1 )
    Ni = 25;
  if ( Nu < Nl ) 
    Nu = Nl + Ni;
/*   Decide what to do and then do it. */
  if ( Sanchez == 0 || Sanchez == 1 || Sanchez == 6) {/* Analyze a data set. */
    cnames = cmatrix(0,numsets,0,MAXNAME);
    nnames = cmatrix(0,numsets,0,MAXNAME);
    for ( ii = 0 ; ii <= numsets ; ii ++ )
      for ( jj = 0 ; jj <= MAXNAME ; jj++ ) *(*(cnames+ii)+jj) = *(*(nnames+ii)+jj) = '\0';
    errflag = get_the_data(counts,inputf,errfile,cnames,nnames,numsets);
    if ( errflag < 0 ) { 
      printf("\nHad trouble reading in the data from %s...\n",inputf);
      exit(-1);
    }
    else 
      errflag = calc_the_est(counts,numsets);
    DoTheStats( outfile,  errfile, inputf, argv[0],  counts,   numsets, nnames, cnames,tex_flag) ;
    if (Sanchez == 1 ) {
      outf = fileopen(outfile, "a");
      SanchezBanner(outf);
      fileclose(outfile, outf);
      
      PrepNormFacts((int) MAXEXACTSS ); 
      for ( ii = 1 ; ii <= numsets ; ii ++ ) {
        jj=SanchezSampleSize( &counts[ii], &counts[0], size, beta ,NULL);
        outf = fileopen(outfile, "a");
        d1 = counts[0].PAAM - counts[0].PAA * counts[0].PM;
        d2 = counts[0].PAaM - counts[0].PAa * counts[0].PM;
	    fprintf(outf,"\t %s by %s\n",nnames[ii],cnames[ii]);        
        SanchezShowModel(outf, &counts[0], d1, d2 );
        fprintf(outf, "\n-exactN  %4d    The exact sample size\n-epower  %6.4f  The achieved power\n-esize   %6.4f  The achieved size\n-asympN   %d", counts[0].eN ,  counts[0].epower,  counts[0].esize,counts[0].samplesize );
        fileclose(outfile, outf);          
      }
      UnPrepNormFacts((int) MAXEXACTSS ); 
    }
    else if ( Sanchez == 6 ) {
      PrepNormFacts((int) MAXEXACTSS ); 
      work = anestvector(0,2);
      for ( ii=1; ii<=numsets; ii++ ) {
	    work[1].PAA = counts[ii].PAA;
	    work[1].PAa = counts[ii].PAa;
	    work[1].Paa = pu;
 	    work[1].PM = counts[ii].PM;
	    work[1].PA =  size;
	    work[1].rho1 = counts[ii].rho1;
	    work[1].rho2 = counts[ii].rho2;
	    work[1].nAAM = Nl;
	    work[1].nAaM = Nu;
	    work[1].naaM = Ni;
	    work[1].nAAm = tex_flag;
	    work[1].nAam = 0;
	    DriveSanchez( work,outfile,  errfile );   
      }
      free_anestvector(work,0,2);
      UnPrepNormFacts((int) MAXEXACTSS );     
    }
    printf("\nOutput to %s",outfile);
    printf("\n\n\t....finished analysis.\n");
    free_cmatrix(cnames,0,numsets,0,MAXNAME);
    free_cmatrix(nnames,0,numsets,0,MAXNAME);
  }
  else if ( Sanchez == 2 )  { /* Do a power curve*/
    counts[1].PAA = pd1;
    counts[1].PAa = pd2;
    counts[1].Paa = pu;
    counts[1].PAAM = d1;
    counts[1].PAaM = d2;
    counts[1].PM = p1d;
    counts[1].PA = size;
    counts[1].rho1 = r1;
    counts[1].rho2 = r2;
    counts[1].nAAM = Nl;
    counts[1].nAaM = Nu;
    counts[1].naaM = Ni;
    counts[1].nAAm = tex_flag;
    counts[1].nAam = dswitch;
    DriveSanchez( counts,outfile,  errfile );   
  }
  else if ( Sanchez == 3 ) { /* Do Table 5*/
    counts[1].Paa = beta;
    counts[1].PA = size;
    SanchezTable5( counts,  outfile,  errfile );
  }
  else if ( Sanchez == 4 ) { /* Calculate the exact sample size for the margins */
 
    counts[1].PAA = pd1;
    counts[1].PAa = pd2;
    counts[1].Paa = beta;
    counts[1].PM = p1d;
    counts[1].PA = size;
    counts[1].rho1 = r1;
    counts[1].rho2 = r2;
    counts[1].PAAM = d1;
    counts[1].PAaM = d2;
    counts[1].nAam = dswitch;
    counts[1].nAAm = tex_flag;
    CalcSanchezN( counts,  outfile,  errfile );
  }
  else if ( Sanchez == 5 )  { /* Calculate asymptotic power for the margins  */
 
    counts[1].PAA = pd1;
    counts[1].PAa = pd2;
    counts[1].PM = p1d;
    counts[1].PA = size;
    counts[1].rho1 = r1;
    counts[1].rho2 = r2;
    counts[1].PAAM = d1;
    counts[1].PAaM = d2;
    counts[1].nAAM = Nl;
    counts[1].nAaM = Nu;
    counts[1].naaM = Ni;
    counts[1].nAam = dswitch;
    counts[1].nAAm = tex_flag;
    DrivePower( counts,outfile,  errfile );   
  }

  free_anestvector(counts,0,numsets);
  free_cvector(inputf,0,MAXNAME);
  free_cvector(outfile,0,MAXNAME);
  free_cvector(errfile,0,MAXNAME);
  return(0);
}







 
  /*
  for ( i=1; i<=numsets; i++ ) {
    theest = *(counts+i);
    minss = SanchezSampleSize(theest, workest, 0.05, 0.5);

    maxp = SanchezMaxMargProb(theest,0,0,0,0);
    printf("\n\n   Max prob of the sample with specified margins data set %d:  %f\n\n",i, maxp);
    rho1 = theest.PAAM *  theest.PAam / (theest.PAAm *  theest.PAaM );
    rho2 = theest.PAaM *  theest.Paam / (theest.PaaM *  theest.PAam );
    p1d = theest.PM;
    pd1 =  theest.PAA;
    pd2 =  theest.PAa;
    SanchezSolveDiseq(rho1, rho2, p1d, pd1, pd2, &d1, &d2);
  }
  
  rho1 = 0.5;  rho2 = 0.5; 
  printf("\n\n PM  PAA PAa     DAAM      DAaM   |  check of rho1, rho2 should be 0.5");
  for ( i=1; i<=9; i++ )
    for ( j=1; j<=7; j++ )
      for ( k=1; k<9-j; k++ ) {
        p1d = (double) i * 0.1;
        pd1 = (double) k * 0.1 ;
        pd2 = (double) j * 0.1;

        SanchezSolveDiseq( rho1, rho2, p1d, pd1, pd2, &d1, &d2 );
        SanchezCalcModelD( &workest, p1d, pd1, pd2, -d1, -d2 );
        printf("\n %3.1f %3.1f %3.1f %9.5f %9.5f  |   %5.2f  %5.2f ", p1d,pd1,pd2, -d1, -d2, workest.rho1, workest.rho2 );
      }
  */
