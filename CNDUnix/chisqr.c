#include "chisqr.h"
#include "utils.h"

#define chinn 500

/***********************************************************************/
/*           Function total = lngamma(num,den)                         */
/*                                                                     */
/*           Calculates the value of ln{(num/den)!}                    */
/*        num is a positive integer and den is 1 or 2.                 */
/***********************************************************************/

double  lngamma(num, den)
  int     num, den;
{
  double xx;
  xx = ( (double) num / (double) den ) + 1.0 ;
  xx = gammln(xx);
  return( xx );
 
}


/* 
This is a modified Num. Rec. in C. function to calculate the log of the factorial.
*/
double lnfactrl(n)
int n;
{
        double a ;
        int j;
        double gammln();
        void nrerror();

        if (n < 0) nrerror("Negative factorial in routine FACTRL");
        if (n > 32) return( gammln(n+1.0) );
        a = 1.0;
        for ( j = 2 ; j<= n ; j++ )
          a = a * (double) j;
        return(log(a));
}



/*
This is a Numerical Recipes in C function  
*/
double gammq(a,x)
double a,x;
{
	double gamser,gammcf,gln;
	void gcf(),gser(),nrerror();

	if (x < 0.0 || a <= 0.0) nrerror("Invalid arguments in routine GAMMQ");
	if (x < (a+1.0)) {
		gser(&gamser,a,x,&gln);
		return 1.0-gamser;
	} else {
		gcf(&gammcf,a,x,&gln);
		return gammcf;
	}
}




#define ITMAX 1000
#define EPS 3.0e-7
/*
This is a Numerical Recipes in C function  
*/
void gcf(gammcf,a,x,gln)
double a,x,*gammcf,*gln;
{
	int n;
	double gold=0.0,g,fac=1.0,b1=1.0;
	double b0=0.0,anf,ana,an,a1,a0=1.0;
	double gammln();
	void nrerror();

	*gln=gammln(a);
	a1=x;
	for (n=1;n<=ITMAX;n++) {
		an=(double) n;
		ana=an-a;
		a0=(a1+a0*ana)*fac;
		b0=(b1+b0*ana)*fac;
		anf=an*fac;
		a1=x*a0+anf*a1;
		b1=x*b0+anf*b1;
		if (a1) {
			fac=1.0/a1;
			g=b1*fac;
			if (fabs((g-gold)/g) < EPS) {
				*gammcf=exp(-x+a*log(x)-(*gln))*g;
				return;
			}
			gold=g;
		}
	}
	nrerror("a too large, ITMAX too small in routine GCF");
}


/*
This is a Numerical Recipes in C function  
*/
void gser(gamser,a,x,gln)
double a,x,*gamser,*gln;
{
	int n;
	double sum,del,ap;
	double gammln();
	void nrerror();

	*gln=gammln(a);
	if (x <= 0.0) {
		if (x < 0.0) nrerror("x less than 0 in routine GSER");
		*gamser=0.0;
		return;
	} else {
		ap=a;
		del=sum=1.0/a;
		for (n=1;n<=ITMAX;n++) {
			ap += 1.0;
			del *= x/ap;
			sum += del;
			if (fabs(del) < fabs(sum)*EPS) {
				*gamser=sum*exp(-x+a*log(x)-(*gln));
				return;
			}
		}
		nrerror("a too large, ITMAX too small in routine GSER");
		return;
	}
}

#undef ITMAX
#undef EPS


/*
This is a Numerical Recipes in C function to calculate ln(Gamma(xx)) for xx > 0.
*/
double gammln(xx)
double xx;
{
	double x,tmp,ser;
	static double cof[6]={76.18009173,-86.50532033,24.01409822,
		-1.231739516,0.120858003e-2,-0.536382e-5};
	int j;

	x=xx-1.0;
	tmp=x+5.5;
	tmp -= (x+0.5)*log(tmp);
	ser=1.0;
	for (j=0;j<=5;j++) {
		x += 1.0;
		ser += cof[j]/x;
	}
	return -tmp+log(2.50662827465*ser);
}



/**************************************************************************/
/*       Function  prob = chiprob(dof,value)                              */
/*                                                                        */
/*  Calculates the chi-squared probability for 'dof' degrees of           */
/* freedom for the 'value' of the X-square.              */
/**************************************************************************/

double  chiprob(dof, value)
  int     dof;
  double  value;
{
  double ddof,prob;
  ddof = (double) dof;
  if ( value > 0.0 )
    prob =  gammq(ddof,value/2.0);
  else
    prob = 1.0;

  return prob;
}
/**************************************************************************/
/*       Function  prob = chiprobexact(locus)                             */
/*                                                                        */
/*  Calculates the chi-squared probability for 'dof' degrees of           */
/* freedom for the 'value' of the X-square.  Uses 'lngamma'.              */
/* n1 is the number of A alleles, n2 the number of a alleles and en the   */
/* sample size. Obs is the observed number of Aa individuals.             */
/**************************************************************************/

double  chiprobexact(counts)
  int    *counts;
{
  int     eins, zwei, drei, vier, funf, sechs, sieben, lindex, ex, na, odd, top;
  int     n1, n2, obs, en;
  double  temp, cumprob, perm, perm2;
  double *px;
  div_t   x;
  px = dvector(0, chinn);
  n1 = *(counts + 0);
  n2 = *(counts + 1);
  obs = *(counts + 3);
  en = (*(counts + 0) + *(counts + 1)) / 2;
  if (n1 > n2)
    na = n2;
  else
    na = n1;
  x = div(en, 2);
  odd = x.rem;
  eins = en + 1;
  zwei = na + 1;
  drei = 2 * en - na + 1;
  sieben = 2 * en + 1;
  perm2 = lngamma(sieben, 1);
  top = en / 2;
  temp = lngamma(eins, 1);
  temp = temp + lngamma(zwei, 1);
  perm = temp + lngamma(drei, 1);
  for (lindex = 0; lindex <= top; lindex++) {
    temp = perm;
    ex = 2 * lindex + odd;
    vier = ((na - ex) / 2) + 1;
    temp = temp + ((double) ex) * log(2.0) - lngamma(vier, 1);
    funf = ex + 1;
    temp = temp - lngamma(funf, 1);
    sechs = en + 1 - ((na + ex) / 2);
    temp = temp - lngamma(sechs, 1);
    temp = temp - perm2;
    px[ex] = exp(temp);
  }
  cumprob = 0;
  for (lindex = 0; lindex <= top; lindex++) {
    ex = 2 * lindex + odd;
    if (px[ex] <= px[obs])
      cumprob = cumprob + px[ex];
  }
  free_dvector(px, 0, chinn);
  return cumprob;
}
