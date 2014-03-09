/*  
          Mmain.c,  main driver for CNDm, a program to calculate
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



/* 
    USAGE:  
  CNdiseq -o outputfile -e errorfile -i inputfile 

If you use an option, be sure to give a value (eg cnisim -e -o will
give nonsense)

  This program will read in the joint cytonuclear counts for
  multiallelic systems.  It will then calculate the frequencies,
  disequilibria and asymptotic variances, as well as the chi sqared
  values for the test statistics.  If wanted, it will also calculate
  approximations to Fisher's exact test using both Monte Carlo and
  Markov chain methods.

*/

#ifdef __MACOS__
#include <console.h>
#endif


int     main(argc, argv)
  int     argc;
  char   *argv[];
{
  FILE   *errorf, *outf;
  char    *inputf,*nowtime;
  char    *outfile, *errfile,**nnames,**cnames;
  int    do_exact,ii,numsets,errflag,verbosity,tex_flag;
  long    seed,reps,MarkovReps,MonteReps,MB,MC,mB,mC;
  double  *tmp,ss;
  adataset *thedata;
/*
  Here is some code needed to get a console window on the Macintosh. 
  I used to use Think C but have switched to Code Warrior.   
*/
#ifdef __MACOS__
 /* make environ safe, and give the console window a new
    title */
  unsigned char abuffer[25] = "\pCNDd\0";
/*  console_options.title = abuffer;  Not used with CodeWarrior*/
/* Simulate the UNIX shell with ccommand, this should also
    be the place to specify redirection of input and
    output. */
  argc = ccommand(&argv);
#endif
/*
  Set some default values for parameters and set the seed for the 
  random number generator.  
*/
  verbosity = 1;
  tex_flag = 0;
  MB = MC = mB = mC = MonteReps = MarkovReps = reps = 0L;
  seed = - get_a_seed();
  ss = ranf(seed);
  do_exact = 1;
/*
  Allocate space for filenames.  Give them default names.  
*/
  inputf = cvector(0,MAXNAME);
  outfile = cvector(0,MAXNAME);
  errfile = cvector(0,MAXNAME);
  strcpy(outfile, "CNDm.out");
  strcpy(errfile, "CNDm.err");
  strcpy(inputf,"CNDm.dat");
/*
 Process the command line parameters.  
*/
  for ( ii = 1 ; ii < argc ; ii++ ) 
    if ( *(*(argv+ii)+0) == '-' )
      switch ( *(*(argv+ii)+1) ) {
        case 'h' :
          printf("\nUSAGE:\n");
          printf("\n   CNDm [-i input] [-o output] [-e logfile] [-r Reps] [-M MonteReps] [-m MarkovReps] [-X ] [-V]\n");
          printf("\nDEFAULTS:\n");
          printf("\n\t [-o %s]",outfile);
          printf("\n\t [-e %s]",errfile);
          printf("\n\t [-i %s]",inputf);
          printf("\n\t [-r %ld ] [-rB %ld ] [-rC %ld ]",reps,MB,MC);
          printf("\n\t [-M %ld ] [-MB %ld ] [-MC %ld ]", MonteReps,MB,MC);            
          printf("\n\t [-m %ld ] [-mB %ld ] [-mC %ld ]", MarkovReps,mB,mC);
          printf("\n\t [-X ] Turns off exact test calculations.");
          printf("\n\t [-V ] Turns off verbosity mode.");
          printf("\n\t [-t ] Output in LaTeX2e format.");
          printf("\n\nNotes: ");
          printf("\n\t1. Setting Reps sets both MonteReps and MarkovReps to Reps.");
          printf("\n\t2. rB is the number of batches, and sets MB and mB.");
          printf("\n\t3. rC is the number of observations per batch, and sets MC and mC.");
          printf("\n\t4. Reps, if positive, will be reset to rB * rC. These default to 10 and 100.");
          printf("\n\t5. Markov and Monte B and C values are set similar to Reps.");
          printf("\n\t6. If you are on a Mac or PC, you need to quit now.\n");
          exit(0);
        case 'e' :
          strcpy(errfile, argv[ii+1]);  break;
        case 'o' : 
          strcpy(outfile, argv[ii+1]);  break;
        case 'X' : 
          do_exact = 0; break;
        case 'r' :
          if (  *(*(argv+ii)+2) == 'B' )
            mB = MB = atol(*(argv+ii+1));
          else if (  *(*(argv+ii)+2) == 'C' )
            mC = MC = atol(*(argv+ii+1));
          else
            MonteReps = MarkovReps = reps = atol(*(argv+ii+1)); 
          break;
        case 'M' :
          if (  *(*(argv+ii)+2) == 'B' )
            MB = atol(*(argv+ii+1));
          else if (  *(*(argv+ii)+2) == 'C' )
            MC = atol(*(argv+ii+1));
          else
            MonteReps = atol(*(argv+ii+1)); 
          break;
        case 'm' :
          if (  *(*(argv+ii)+2) == 'B' )
            mB = atol(*(argv+ii+1));
          else if (  *(*(argv+ii)+2) == 'C' )
            mC = atol(*(argv+ii+1));
          else
            MarkovReps = atol(*(argv+ii+1)); 
          break;
        case 'V' :
          verbosity = 0; break;
        case 't' :
          tex_flag = 1; break;
        case 'i' : 
        case 'd' : 
          strcpy(inputf, argv[ii+1]);  break;
       }
  nowtime = asctime2();

/*
  Check on the values for the Monte Carlo and Markov Chain runs.
  Need to set MB, MC, mB, mC if they aren't set already.
  
  Defaults:  0 for all if nothing set.
  If TotalReps are set, then 10 batches and TotalReps/10 reps in each batch.
  If TotalReps were set and batches were set, then 100 reps in each batch.
  At end, reset TotalReps = batches * reps
*/
  if ( MonteReps > 0L || MB > 0L || MC > 0L) {
    if ( MB == 0L ) {
      MB =  10L ;
      MC = MonteReps / MB;
    }  
    if ( MC == 0L )
      MC = 100L;
    MonteReps = MB * MC;
  }
  if ( MarkovReps > 0L  || mB > 0L || mC > 0L) {
    if ( mB == 0L ) {
      mB =  10L ;
      mC = MarkovReps / mB;
    }  
    if ( mC == 0L )
      mC = 100L;
    MarkovReps = mB * mC;
  }


/*
  Open up the output file and print out all the parameters used, as well
  as the time run.
*/
  outf = fileopen(outfile, "a");
  if ( outf == NULL ) outf = fileopen(outfile, "w");
  if (outf != NULL) {
    if ( tex_flag == 0 ) {
	  fprintf(outf,"\n Analysis of program %s\n",argv[0]);
          fprintf(outf, " This program was run at %s\n",nowtime);
          fprintf(outf,"\n\t Output file: [-o %s]",outfile);
          fprintf(outf,"\n\t Error file:  [-e %s]",errfile);
          fprintf(outf,"\n\t Input file:  [-i %s]",inputf);
          if ( do_exact == 1 )
            fprintf(outf,"\n\t Exact tests for 2 by 2 tables will be calculated: no [-X ]");
          if ( MonteReps > 0L ) 
            fprintf(outf,"\n\t Reps, B and C for Monte Carlo test: [-M %ld ] [-MB %ld ] [-MC %ld ]\n", MonteReps,MB,MC);            
          if ( MarkovReps > 0L ) 
            fprintf(outf,"\n\t Reps, B and C for Markov Chain test: [-m %ld ] [-mB %ld ] [-mC %ld ]\n", MarkovReps,mB,mC);
          if ( MonteReps == 0L && MarkovReps == 0L )
            fprintf(outf,"\n\t No Monte Carlo or Markov approximation for the exact test...\n\n");
      }
   else {
	  fprintf(outf,"\\documentclass[12pt]{article}");
	  fprintf(outf,"\n\\setlength{\\textwidth}{6.5in}");
	  fprintf(outf,"\n\\setlength{\\textheight}{8.0in}");
	  fprintf(outf,"\n\\setlength{\\topmargin}{0.0in}");
	  fprintf(outf,"\n\\setlength{\\oddsidemargin}{0in}"); 
	  fprintf(outf,"\n\\begin{document}");
	  fprintf(outf,"\n\\centerline{{\\large {\\bf CNDm Analysis}}}\n\\vspace{3ex}\n");	  
	  fprintf(outf,"\n The analysis of program %s",argv[0]);
      fprintf(outf, " was run at %s\n",nowtime);
      fprintf(outf,"\nThe output was written to %s,",outfile);
      fprintf(outf,"errors were logged to %s, and ",errfile);
      fprintf(outf,"the data came from %s.\n",inputf);
      if ( do_exact == 1 )
        fprintf(outf,"  Exact tests for 2 by 2 tables were calculated.");
      if ( MonteReps > 0L ) 
        fprintf(outf,"  A Monte Carlo approximation with  %ld reps in  %ld batches of %ld observations per batch was preformed.", MonteReps,MB,MC);            
      if ( MarkovReps > 0L ) 
        fprintf(outf,"  A Markov Chain test with %ld reps in %ld batches of %ld observations per batch was performed.", MarkovReps,mB,mC);
      if ( MonteReps == 0L && MarkovReps == 0L )
        fprintf(outf," Neither a Monte Carlo nor a Markov approximation was performed.\n");    
    }
  }
  fileclose(outfile, outf);
/*  If verbosity = 1, then print appropriate info to the screen... */
  if (verbosity == 1 ) {
	      printf("\n Analysis of program %s\n",argv[0]);
          printf(" It is presently %s\n",nowtime);
          printf("\n\t Output file: [-o %s]",outfile);
          printf("\n\t Error file:  [-e %s]",errfile);
          printf("\n\t Input file:  [-i %s]",inputf);
          if ( do_exact == 1 )
            printf( "\n\t Exact tests for 2 by 2 tables will be calculated: no [-X ]");
          if ( MonteReps > 0L ) 
            printf("\n\t Reps, B and C for Monte Carlo test: [-M %ld ] [-MB %ld ] [-MC %ld ]\n", MonteReps,MB,MC);            
          if ( MarkovReps > 0L ) 
            printf("\n\t Reps, B and C for Markov Chain test: [-m %ld ] [-mB %ld ] [-mC %ld ]\n", MarkovReps,mB,mC);
          if ( MonteReps == 0L && MarkovReps == 0L )
            printf( "\n\t No Monte Carlo or Markov approximation for the exact test...\n\n");
      }
/*
  Open up the log file and print out all the parameters used, as well as the time run.
*/
  outf = fileopen(errfile, "a");
  if ( outf == NULL ) outf = fileopen(errfile, "w");
  if (outf != NULL) {
	  fprintf(outf,"\n  Error messages from program %s\n",argv[0]);
          fprintf(outf, " This program was run at %s\n",nowtime);
          fprintf(outf,"\n\t Output file: [-o %s]",outfile);
          fprintf(outf,"\n\t Error file:  [-e %s]",errfile);
          fprintf(outf,"\n\t Input file:  [-i %s]",inputf);
          if ( do_exact == 1 )
            fprintf(outf,"\n\t Exact tests for 2 by 2 tables will be calculated: no [-X ]");
          if ( MonteReps > 0L ) 
            fprintf(outf,"\n\t Reps, B and C for Monte Carlo test: [-M %ld ] [-MB %ld ] [-MC %ld ]\n", MonteReps,MB,MC);            
          if ( MarkovReps > 0L ) 
            fprintf(outf,"\n\t Reps, B and C for Markov Chain test: [-m %ld ] [-mB %ld ] [-mC %ld ]\n", MarkovReps,mB,mC);
          if ( MonteReps == 0L && MarkovReps == 0L )
            fprintf(outf,"\n\t No Monte Carlo or Markov approximation for the exact test...\n\n");
      }
  fileclose(errfile, outf);
/*
  Get the data from the specified file.  If it is not read correctly,
  then print a message to the log file and bail out.
*/
  numsets =  get_numsets(inputf,errfile);
  thedata = get_data(inputf,errfile,numsets);
  if ( thedata == NULL ) { 
    outf = fileopen(errfile, "a");
    fprintf(outf,"\nHad trouble reading in the data from %s...\n",inputf);
    fileclose(errfile, outf);
    exit(-1);
  }  
  if (verbosity == 1) 
	  printf("\t About to do asymptotic analysis\n" );
  allele_counts(thedata,numsets,errfile);
  joint_counts(thedata,numsets,errfile);
/*
  Calculate the asymptotic statistics for the genotypic disequilibria
  and print them to the output file.
*/  
  calc_afreq(thedata,numsets,errfile,0);
  calc_adiseq(thedata,numsets,errfile,0);
  calc_ndiseq(thedata,numsets,errfile,0);
  calc_avars(thedata,numsets,errfile,0);
  calc_astat(thedata,numsets,errfile,0);
  calc_aprs(thedata,numsets,errfile,0);
  if ( tex_flag == 0 )
    write_asymp_results(thedata,numsets,outfile,errfile,0 );
  else
    write_asymp_results_tex(thedata,numsets,outfile,errfile,0 );
  if (verbosity == 1) 
	  printf("\t Finished asymptotic analysis\n" );
/*
  Calculate the asymptotic statistics for the allelic disequilibria 
  and print them to the output file.
*/
  calc_afreq(thedata,numsets,errfile,1);
  calc_adiseq(thedata,numsets,errfile,1);
  calc_ndiseq(thedata,numsets,errfile,1);
  calc_avars(thedata,numsets,errfile,1);
  calc_astat(thedata,numsets,errfile,1);
  calc_aprs(thedata,numsets,errfile,1);
  if ( tex_flag == 0 )
    write_asymp_results(thedata,numsets,outfile,errfile,1 );
  else
    write_asymp_results_tex(thedata,numsets,outfile,errfile,1 );  
/*
  If do_exact is equal to 1, then calculate exact probabilities
  using Dmitri Zaykin's subroutines.  If MonteReps or MarkovReps are
  greater than 0, then do a Monte Carlo or 
  Markov chain approximation to the disequilibrium for the entire
  data set.
*/
  if (do_exact == 1 && verbosity == 1) 
	  printf("\t About to do exact tests... \n" );
  if ( do_exact == 1 )
    calc_Exact_prs(thedata,numsets,outfile,errfile,tex_flag,0);
  if ( do_exact == 1 )
    calc_Exact_prs(thedata,numsets,outfile,errfile,tex_flag,1);
  if (do_exact == 1 && verbosity == 1) 
	  printf("\t Finished exact tests... \n" );

  if (MonteReps > 0L && verbosity == 1) 
	  printf("\t About to do Monte Carlo Shuffle Test... \n" );
  if ( MonteReps > 0L ) {
    calc_MonteCarlo_prs(thedata,numsets,outfile,errfile,MonteReps,MB,MC,verbosity,tex_flag,0);
    calc_MonteCarlo_prs(thedata,numsets,outfile,errfile,MonteReps,MB,MC,verbosity,tex_flag,1);
  }
  if (MonteReps > 0L && verbosity == 1) 
	  printf("\t Finished Monte Carlo Shuffle Test... \n" );
  if (MarkovReps > 0L && verbosity == 1) 
	  printf("\t About to do Markov Chain Monte Carlo Test... \n" );
  if ( MarkovReps > 0L) {
    calc_Markov_prs(thedata,numsets,outfile,errfile,MarkovReps,mB,mC,verbosity,tex_flag,0);
    calc_Markov_prs(thedata,numsets,outfile,errfile,MarkovReps,mB,mC,verbosity,tex_flag,1);
  }
  if (MarkovReps > 0L && verbosity == 1) 
	  printf("\t Finished Markov Chain Monte Carlo Test... \n" );

/*
  Write the time and the name of the output file to the
  log file.
*/
  nowtime = asctime2();
  outf = fileopen(errfile, "a");
  fprintf(outf, "\n\n\t It is %s\n",nowtime);
  fprintf(outf,"\nOutput to %s",outfile);
  fprintf(outf,"\n\n\t....finished analysis.\n");
  fileclose(errfile, outf);
  outf = fileopen(outfile, "a");
  if ( tex_flag == 0 ) {
    fprintf(outf, "\n\n\t It is %s\n",nowtime);
    fprintf(outf,"\nOutput to %s",outfile);
    fprintf(outf,"\n\n\t....finished analysis.\n");
  }
  else {
    fprintf(outf, "\n\n{\\bf Note:} Analysis finished at %s\n",nowtime);
    fprintf(outf,"\n\\end{document}\n");  
  }
  fileclose(outfile, outf);

  if (  verbosity == 1) {
    printf("\n\n\t It is %s\n",nowtime);
    printf("\n\t Output to %s",outfile);
    printf("\n\t ....finished analysis, deallocating memory.\n");
  }
/*
  Time to deallocate the memory used.
*/
  free_datavector(thedata,1,numsets); 
  free_cvector(inputf,0,MAXNAME);
  free_cvector(outfile,0,MAXNAME);
  free_cvector(errfile,0,MAXNAME);
  if (  verbosity == 1) 
    printf("\t Exiting... Quit window if Macintosh\n\n");
  exit(0);
}
