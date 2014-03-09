/*
 * These are some functions that I wrote that are useful
 * from time to time. Some are available in some
 * implementations.
 *
 *
 */


#include "utils.h"

#include <time.h>
long ix;


char get_next_token(xtemp,nn,fileptr)
  char *xtemp;
  int nn;
  FILE *fileptr;
{
  int ii;
  char ch;
  for ( ii = 0 ; ii < nn ; ii++ ) 
    *(xtemp+ii) = '\0';
  ch = fgetc(fileptr);
  while( isspace(ch) )
    ch = fgetc(fileptr);
  for ( ii = 0 ; ii < nn && !isspace(ch) && ch != EOF ; ii++ ) {
    *(xtemp+ii) = ch;
    ch = fgetc(fileptr);
  }
  return(ch);
}



long get_a_seed()
{
  time_t tptr;
  time(&tptr);
#if defined(__MACOS__)
  tptr = tptr / 2L;
#endif
  return ((long) tptr);
}

long marktime() {
  time_t  mark;
  time(&mark);
  return( (long) mark );
}

char *asctime2()
{
  static char time_buffer[MAXNAME];
  time_t tptr;
  struct tm *tms;
  size_t len;
  time(&tptr);
  tms = localtime(&tptr);
  len = strftime( time_buffer, MAXNAME, "%H:%M:%S on %A, %d %B %Y\n", tms);
  if ( len == 0 ) 
    return NULL;
  else
    return time_buffer;
}


/*
 * dtranspose does an arbitrary transpose of one matrix
 * onto another. You must have allocated the memory that
 * mm1 and mm2 have pointed to.  lr, and lc are the initial
 * row and column while ur and uc are the final row and
 * column for the transpositon.  Note that mm1 =
 * dmatrix(lr,ur,lc,uc); mm2 = dmatrix(lc,uc,lr,ur); at
 * least.  These can be bigger, so that one can transpose
 * an arbitrary part of one matrix onto another.
 *
 * By Chris Basten, January 1994
 */
int
dtranspose(mm1, mm2, lr, lc, ur, uc)
  double        **mm1, **mm2;
  int             lr, lc, ur, uc;
{
  int             ii, jj;
  if (lr >= ur || lc >= uc)
    return (1);
  for (ii = lr; ii <= ur; ii++)
    for (jj = lc; jj <= uc; jj++)
      *(*(mm2 + jj) + ii) = *(*(mm1 + ii) + jj);
  return (0);
}







#ifdef DIVT

div_t
div(numer, denom)
  int             numer;
  int             denom;
{
  div_t           x;
  x.quot = 0;
  x.rem = 0;
  if (denom != 0) {
    x.quot = numer / denom;
    x.rem = numer - x.quot * denom;
  }
  else
    fprintf(stderr, "\n\nSorry, divide by zero in _div...\n");
  return x;
}

ldiv_t
ldiv(numer, denom)
  long            numer;
  long            denom;
{
  ldiv_t          x;
  x.quot = 0;
  x.rem = 0;
  if (denom != 0) {
    x.quot = numer / denom;
    x.rem = numer - x.quot * denom;
  }
  else
    fprintf(stderr, "\n\nSorry, divide by zero in _ldiv...\n");
  return x;
}
#endif

FILE           *
fileopen(name, mode)
  char           *name;
  char           *mode;
{
  /* Open files, writing an error message if it fails. */
  FILE           *fp = fopen(name, mode);
  if (fp == NULL) {
    printf("Can't open %s for ", name);
    switch (mode[0]) {
    case 'r':
      printf("reading\n");
      break;
    case 'w':
      printf("writing\n");
      break;
    case 'a':
      printf("appending\n");
      break;
    default:
      printf("some strange mode\n");
      break;
    }
  }
  return fp;
}

void
fileclose(name, fp)
  char           *name;
  FILE           *fp;
{
  /* Close files, writing an error message if it fails. */
  if (fp != NULL && fclose(fp) == EOF)
    printf("Error closing %s\n", name);
}

void
pause()
{
  char            answer;
  printf("\nType a <CR> to go on...\n");
  while ((answer = getchar()) != '\n');
}


void
get_field(xfield, xtemp, xbuffer)
  int             xfield;
  char           *xtemp, *xbuffer;
{
  int             i, j, field, look;
  for (i = 0; i < MAXNAME; i++)
    *(xtemp + i) = '\0';
  for (i = 0; *(xbuffer + i) != '\0'; i++);
  *(xbuffer + i) = ' ';
  field = 0;
  look = 1;
  for (i = 0; field < xfield; i++) {
    if (*(xbuffer + i) != ' ' && look == 1) {
      field = field + 1;
      look = 0;
    }
    if (*(xbuffer + i) == ' ' && look == 0)
      look = 1;
  }
  for (j = i - 1; !isspace(*(xbuffer + j)); j++)
    *(xtemp + j - i + 1) = *(xbuffer + j);
}


#define  a 16807
#define  b15 32768
#define  b16 65536
#define  p 2147483647
#define  xnorm 4.656612875E-10
/*  Generates a random number in the interval [0,1] */
double ranf(inix)
long inix;
{
  long xhi, xalo, leftlo, fhi, k;
  double xx;
  if (inix < 0)
    ix = -1 * inix;
  xhi = ix / b16;
  xalo = (ix - xhi * b16) * a;
  leftlo = xalo / b16;
  fhi = xhi * a + leftlo;
  k = fhi / b15;
  ix = (((xalo - leftlo * b16) - p) + (fhi - k * b15) * b16) + k;
  if (ix < 0)
    ix = ix + p;
  xx = (double) ix *xnorm;
  return xx;
}
#undef  a  
#undef  b15  
#undef  b16  
#undef  p  
#undef  xnorm  




/*
 * The following come from Numerical Recipes in C, and
 * handle the allocation of memory.  The character
 * allocators were modified by me.
 */


void
nrerror(error_text)
  char            error_text[];
{
  fprintf(stderr, "Numerical Recipes run-time error...\n");
  fprintf(stderr, "%s\n", error_text);
  fprintf(stderr, "...now exiting to system...\n");
  exit(1);
}

float          *
vector(nl, nh)
  int             nl, nh;
{
  float          *v;
  v = (float *) malloc((unsigned) (nh - nl + 1) * sizeof(float));
  if (!v)
    nrerror("allocation failure in vector()");
  return v - nl;
}

double         *
dvector(nl, nh)
  int             nl, nh;
{
  double         *v;
  v = (double *) malloc((unsigned) (nh - nl + 1) * sizeof(double));
  if (!v)
    nrerror("allocation failure in dvector()");
  return v - nl;
}

int            *
ivector(nl, nh)
  int             nl, nh;
{
  int            *v;
  v = (int *) malloc((unsigned) (nh - nl + 1) * sizeof(int));
  if (!v)
    nrerror("allocation failure in ivector()");
  return v - nl;
}
char            *
cvector(nl, nh)
  int             nl, nh;
{
  char           *v;
  v = (char *) malloc((unsigned) (nh - nl + 1) * sizeof(char));
  if (!v)
    nrerror("allocation failure in cvector()");
  return v - nl;
}

long           *
lvector(nl, nh)
  int             nl, nh;
{
  long           *v;
  v = (long *) malloc((unsigned) (nh - nl + 1) * sizeof(long));
  if (!v)
    nrerror("allocation failure in lvector()");
  return v - nl;
}

char          **
cmatrix(nrl, nrh, ncl, nch)
  int             nrl, nrh, ncl, nch;
{
  int             i;
  char          **m;
  m = (char **) malloc((unsigned) (nrh - nrl + 1) * sizeof(char *));
  if (!m)
    nrerror("allocation failure 1 in matrix()");
  m -= nrl;
  for (i = nrl; i <= nrh; i++) {
    m[i] = (char *) malloc((unsigned) (nch - ncl + 1) * sizeof(char));
    if (!m[i])
      nrerror("allocation failure 2 in cmatrix()");
    m[i] -= ncl;
  }
  return m;
}


int           **
imatrix(nrl, nrh, ncl, nch)
  int             nrl, nrh, ncl, nch;
{
  int             i, **m;
  m = (int **) malloc((unsigned) (nrh - nrl + 1) * sizeof(int *));
  if (!m)
    nrerror("allocation failure 1 in imatrix()");
  m -= nrl;

  for (i = nrl; i <= nrh; i++) {
    m[i] = (int *) malloc((unsigned) (nch - ncl + 1) * sizeof(int));
    if (!m[i])
      nrerror("allocation failure 2 in imatrix()");
    m[i] -= ncl;
  }
  return m;
}

float         **
matrix(nrl, nrh, ncl, nch)
  int             nrl, nrh, ncl, nch;
{
  int             i;
  float         **m;
  m = (float **) malloc((unsigned) (nrh - nrl + 1) * sizeof(float *));
  if (!m)
    nrerror("allocation failure 1 in matrix()");
  m -= nrl;

  for (i = nrl; i <= nrh; i++) {
    m[i] = (float *) malloc((unsigned) (nch - ncl + 1) * sizeof(float));
    if (!m[i])
      nrerror("allocation failure 2 in matrix()");
    m[i] -= ncl;
  }
  return m;
}



short         **
smatrix(nrl, nrh, ncl, nch)
  int             nrl, nrh, ncl, nch;
{
  int             i;
  short         **m;
  m = (short **) malloc((unsigned) (nrh - nrl + 1) * sizeof(short *));
  if (!m)
    nrerror("allocation failure 1 in imatrix()");
  m -= nrl;

  for (i = nrl; i <= nrh; i++) {
    m[i] = (short *) malloc((unsigned) (nch - ncl + 1) * sizeof(short));
    if (!m[i])
      nrerror("allocation failure 2 in smatrix()");
    m[i] -= ncl;
  }
  return m;
}

void
free_vector(v, nl, nh)
  float          *v;
  int             nl, nh;
{
  free((char *) (v + nl));
}

void
free_dvector(v, nl, nh)
  double         *v;
  int             nl, nh;
{
  free((char *) (v + nl));
}

void
free_ivector(v, nl, nh)
  int            *v, nl, nh;
{
  free((char *) (v + nl));
}
void
free_cvector(v, nl, nh)
  char            *v;
  int  nl, nh;
{
  free((char *) (v + nl));
}

void
free_lvector(v, nl, nh)
  long           *v;
  int             nl, nh;
{
  free((char *) (v + nl));
}


void
free_cmatrix(m, nrl, nrh, ncl, nch)
  char          **m;
  int             nrl, nrh, ncl, nch;
{
  int             i;
  for (i = nrh; i >= nrl; i--)
    free((char *) (m[i] + ncl));
  free((char *) (m + nrl));
}

void
free_imatrix(m, nrl, nrh, ncl, nch)
  int           **m;
  int             nrl, nrh, ncl, nch;
{
  int             i;
  for (i = nrh; i >= nrl; i--)
    free((char *) (m[i] + ncl));
  free((char *) (m + nrl));
}
void 
free_matrix(m, nrl, nrh, ncl, nch)
  float         **m;
  int             nrl, nrh, ncl, nch;
{
  int             i;
  for (i = nrh; i >= nrl; i--)
    free((char *) (m[i] + ncl));
  free((char *) (m + nrl));
}


void 
free_smatrix(m, nrl, nrh, ncl, nch)
  short         **m;
  int             nrl, nrh, ncl, nch;
{
  int             i;
  for (i = nrh; i >= nrl; i--)
    free((char *) (m[i] + ncl));
  free((char *) (m + nrl));
}


double        **
dmatrix(nrl, nrh, ncl, nch)
  int             nrl, nrh, ncl, nch;
{
  int             i,j;
  double        **m;
  m = (double **) malloc((unsigned) (nrh - nrl + 1) * sizeof(double *));
  if (!m)
    nrerror("allocation failure 1 in dmatrix()");
  m -= nrl;

  for (i = nrl; i <= nrh; i++) {
    m[i] = (double *) malloc((unsigned) (nch - ncl + 1) * sizeof(double));
    if (!m[i])
      nrerror("allocation failure 2 in dmatrix()");
    m[i] -= ncl;
  }
  for ( i = nrl ; i <= nrh ; i++ )
    for ( j = ncl ; j <= nch ; j++ )
      *(*(m+i)+j) = 0.0;
  return m;
}


float         **
submatrix(aa, oldrl, oldrh, oldcl, oldch, newrl, newcl)
  float         **aa;
  int             oldrl, oldrh, oldcl, oldch, newrl, newcl;
{
  int             i, j;
  float         **m;
  m = (float **) malloc((unsigned) (oldrh - oldrl + 1) * sizeof(float *));
  if (!m)
    nrerror("allocation failure in submatrix()");
  m -= newrl;

  for (i = oldrl, j = newrl; i <= oldrh; i++, j++)
    m[j] = aa[i] + oldcl - newcl;

  return m;
}

void
free_dmatrix(m, nrl, nrh, ncl, nch)
  double        **m;
  int             nrl, nrh, ncl, nch;
{
  int             i;
  for (i = nrh; i >= nrl; i--)
    free((char *) (m[i] + ncl));
  free((char *) (m + nrl));
}


void
free_submatrix(bb, nrl, nrh, ncl, nch)
  float         **bb;
  int             nrl, nrh, ncl, nch;
{
  free((char *) (bb + nrl));
}



float         **
convert_matrix(aa, nrl, nrh, ncl, nch)
  float          *aa;
  int             nrl, nrh, ncl, nch;
{
  int             i, j, nrow, ncol;
  float         **m;
  nrow = nrh - nrl + 1;
  ncol = nch - ncl + 1;
  m = (float **) malloc((unsigned) (nrow) * sizeof(float *));
  if (!m)
    nrerror("allocation failure in convert_matrix()");
  m -= nrl;
  for (i = 0, j = nrl; i <= nrow - 1; i++, j++)
    m[j] = aa + ncol * i - ncl;
  return m;
}



void
free_convert_matrix(bb, nrl, nrh, ncl, nch)
  float         **bb;
  int             nrl, nrh, ncl, nch;
{
  free((char *) (bb + nrl));
}

/*****************************************************************************
 *        This is a general Gamma variate generator with the shape
 *     parameter beta > 0.25. Coded from the Algorithm GBH of Cheng and Feast
 *     (1980 Communications of the ACM 23:389-394) by Zhao-Bang Zeng on Jan. 30,
 *     1992.
 *     This was translated into c by Chris Basten, 19 Jan. 1994.
 *****************************************************************************/
double
gamgbh(beta, xix)
  double          beta;
  long            xix;
{
  double          aa, bb, c, d, t, h1, h2, u, u1, u2, w, tmp;

  if (beta <= 0.25)
    nrerror("Error in gamgbh:  beta <= 0.25...");
  aa = beta - 0.25;
  bb = beta / aa;
  c = 2.0 / aa;
  d = c + 2.0;
  t = 1.0 / sqrt(beta);
  h1 = (0.4417 + 0.0245 * t / beta) * t;
  h2 = (0.222 - 0.043 * t) * t;
  do {
    u = ranf(xix);
    u1 = ranf(xix);
    u2 = u1 + h1 * u - h2;
    if (u2 <= 0.0 || u2 >= 1.0)
      tmp = 1.0;
    else {
      w = bb * pow((u1 / u2), 4.0);
      if (c * u2 - d + w + 1.0 / w <= 0.0)
	return (aa * w);
      else
	tmp = c * log(u2) - log(w) + w - 1.0;
    }
  } while (tmp >= 0.0);
  return (aa * w);
}

/****************************************************************************
      This gamma generator generates gamma variates by composition-
   rejection from the Weibull density for the shape parameter beta
   smaller than 1 (0 < beta < 1). Coded from I. Vaduva (1977 Math.
   Operationsforsch. Statist., Ser. Statistics, Vol. 8:545-576) by
   Zhao-Bang Zeng on Feb. 5 1992.   (Best one of this kind)

   Translated into c by Chris Basten on 19 January 1994.
 *****************************************************************************/
double
gamnl1(beta, aa, bb, pp, xix)
  double          beta, aa, bb, pp;
  long            xix;
{
  double          gamnl1r, u1, u2, u0;

  if (beta >= 1.0 || beta <= 0.0)
    nrerror("beta in gamnl1 is not in (0,1)");


  if (ranf(xix) > pp) {
    u0 = ranf(xix);
    gamnl1r = u0 = pow(u0, aa);
    u1 = ranf(xix);
    while (u0 >= u1) {
      u2 = ranf(xix);
      if (u1 < u2) {
	u0 = ranf(xix);
	gamnl1r = u0 = pow(u0, aa);
      }
      else
	u0 = u2;
      u1 = ranf(xix);
    }
  }
  else
    do {
      u0 = ranf(xix);
      u0 = pow(u0, bb);
      gamnl1r = 1.0 - log(ranf(xix));
    } while (gamnl1r >= u0);
  return (gamnl1r);
}


/*
  Return a random integer from 1 to in, using ranf(ix) as the
  source of uniform deviates.
     Fortran source from John Monahan,
     and translated into c by Chris Basten 19 January 1994.
*/
long
iran(xix, in)
  long            xix, in;
{
  double          rv;
  rv = ranf(xix);
  return ((long) (rv * (double) in) + 1);
}

