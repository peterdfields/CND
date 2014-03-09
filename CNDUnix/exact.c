/* 
  This was programmed by Dimitry Zaykin of North Carolina State University.
  He sent it to Chris Basten on 10 May 1994.
*/
/*  Fisher's exact test */
#include <stdio.h>
#include <math.h>
#ifndef M_PI
#define  M_PI    3.14159265358979323846
#endif
#define DMY     (-999.999)
#define MIN(a,b)    (((a) < (b)) ? (a) : (b))
char exactfit;

  struct tw {
    int   x1,x2,x3,x4,     /* entries                                */
	  m1,m2,m3,m4,     /* marginals: m1=1*, m2=2*, m3=*1, m4=*2. */
	  n;               /* sum                                    */
    float e1,e2,e3,e4,     /* expected                               */
	  pf,              /* closest tail probability               */
	  pb;              /* both tails probability of              */
  };

  struct entry {
	int   x1,x2,x3,x4, /* entries                                */
	m1,m2,m3,m4,       /* marginals: m1=1*, m2=2*, m3=*1, m4=*2. */
	n;                 /* sum                                    */
  };


char Exact(int [][2], float *both_tail, float *closest);
static  char    minx(struct entry *),
		twotail_exact(struct tw *, double *, double *);
static          struct entry copy_entry(struct tw *);
static  double  hypergeom(struct entry *);
static  void    GetExpec(struct tw *),
		getparam(int [][2], struct tw *);

static  long
double  fct(long double x), appr_fct (long double) ;
static  int     ABS( int ), search(int *, int );

int search( int *x, int num)
{  /*  returns the index of the maximal element from an array of size num  */
   int  i,
	t=0, mx=x[0];

      for(i=1; i<num; i++) {
	if( x[i] > mx ) {
	    t=i;
	    mx=x[i];
	}
      }
   return t;
}

char minx (struct entry *a)
{
 return( ( a->x1 * a->x4 - a->x2 * a->x3  ) > 0 ? 2:1) ;
}


 char Exact(int cells[][2], float *both_tail, float *closest)
 {
  struct tw it;
  double s_tai, prob;
    getparam(cells, &it);
    GetExpec( &it);
    if( twotail_exact(&it, &prob, &s_tai) == 1) {
     *both_tail = (float)prob;
     *closest =   (float)s_tai;
     return 1;
    }
    return -1;
 }


 void getparam(int cells[][2], struct tw *i)
 {
  i->x1 = cells[0][0]; i->x2 = cells[0][1];
  i->x3 = cells[1][0]; i->x4 = cells[1][1];
  return;
 }

    void GetExpec(struct tw *mm)
    {
      mm->m1 = mm->x1 + mm->x2;
      mm->m2 = mm->x3 + mm->x4;
      mm->m3 = mm->x1 + mm->x3;
      mm->m4 = mm->x2 + mm->x4;
      mm->n  = mm->x1 + mm->x2 + mm->x3 + mm->x4;
      mm->e1 = (float)((double)mm->m1 * (double)mm->m3 / (double)mm->n);
      mm->e2 = (float)((double)mm->m1 * (double)mm->m4 / (double)mm->n);
      mm->e3 = (float)((double)mm->m2 * (double)mm->m3 / (double)mm->n);
      mm->e4 = (float)((double)mm->m2 * (double)mm->m4 / (double)mm->n);
      return;
    }


char twotail_exact(struct tw *m, double *cp, double *sp)
#define PR (*cp)
#define FT (*sp)
  {
  char k;
  struct entry em;
  double sec_prob, check, orig_p;
  int orig_cross;
        em=copy_entry(m);
        if( !(em.x1 * em.x4 - em.x2 * em.x3) ) {
	  printf("exact fit -> det = 0\n");
          *cp=1.;
          *sp=0.5;
          exactfit=1;
          return 1;
        }

                            /* probability of the original table */
        if((PR=hypergeom(&em)) < 0.0) {
	   PR=FT=DMY;
	   return -1;
        }

        k=minx(&em);
    orig_p=PR;
    switch (k) {   /* first tail */

      case 1:
        while( em.x1 && em.x4 ) {
           --em.x1;       /* decrement 11-22 diagonal */
           --em.x4;
           ++em.x2;
           ++em.x3;
           if((check=hypergeom(&em)) < 0.0) {
                        PR=FT=DMY;
			return -1;
           }
           PR+=check;
        }
      break;

      case 2:
        while( em.x2 && em.x3 ) {
           --em.x2;        /* decrement 12-21 diagonal */
           --em.x3;
           ++em.x1;
           ++em.x4;
           if((check=hypergeom(&em)) < 0.0) {
			PR=FT=DMY;
                        return -1;
           }
           PR+=check;
        }
      break;
    }
#define X1 em.x1
#define X2 em.x2
#define X3 em.x3
#define X4 em.x4


    /* second tail */

    em=copy_entry(m);
    orig_cross = ABS( X1*X4 - X2*X3 );
    sec_prob=0.0;
    switch (k) {

        case 1:
        while( em.x2 && em.x3 ) {
           --em.x2;        /* decrement 12-21 diagonal */
           --em.x3;
           ++em.x1;
           ++em.x4;
	   if ( ABS( X1*X4 - X2*X3 ) < orig_cross ) continue;
           if((check=hypergeom(&em)) < 0.0) {
                        PR=FT=DMY;
                        return -1;
           }
           if(check > orig_p) continue;
           sec_prob+=check;
        }
        break;

        case 2:
        while( em.x1 && em.x4 ) {
           --em.x1;        /* decrement  11-22 diagonal */
	   --em.x4;
           ++em.x2;
           ++em.x3;
           if ( ABS( X1*X4 - X2*X3 ) < orig_cross  ) continue;
           if((check=hypergeom(&em)) < 0.0) {
                        PR=FT=DMY;
                        return -1;
           }
           if(check > orig_p) continue;
           sec_prob+=check;
        }
        break;
    }

    /* probabilities  */
    FT=PR;         /* first tail (closest) */
    PR+=sec_prob;  /* both tails */
    return 1;

    #undef PR
    #undef FT
    #undef X1
    #undef X2
    #undef X3
    #undef X4
  }


  int ABS ( int x) { return (x < 0 ?  -x:x); }

  struct entry copy_entry(struct tw *a)
  {
     struct entry e;
     e.x1=a->x1; e.x2=a->x2; e.x3=a->x3;  e.x4=a->x4;
     e.m1=a->m1; e.m2=a->m2; e.m3=a->m3;  e.m4=a->m4;
     e.n=a->n;
     return e;
  }


  double hypergeom(struct entry *a)
  {
  long double PROB;
  int           x[4], y[5],i,j,v,
                kx,ky,max_xy,min_xy;
  long double   ux,uy;

  x[0]=a->m1,   x[1]=a->m2,   /* x[] - nominator terms */
  x[2]=a->m3,   x[3]=a->m4;

  y[0]=a->x1,   y[1]=a->x2,   /* y[] - denominator terms */
  y[2]=a->x3,   y[3]=a->x4,
  y[4]=a->n;

         PROB=ux=uy=1.;
         /*_control87(EM_INEXACT|EM_OVERFLOW|EM_UNDERFLOW,MASKALL);*/
         for(;;) {
           i=search(x,4);
	   j=search(y,5);
                  if(!x[i]) break;
           kx=x[i];
           ky=y[j];
           x[i]=y[j]=0;
                if(kx == ky) continue; /* this pair is canceled out */
           min_xy=MIN(kx,ky);
           max_xy = (min_xy == kx) ? ky : kx;

             if (max_xy == ky)
              for(v = min_xy+1; v <= max_xy; ++v) { uy*=v; }
             else
             if (max_xy == kx)
              for(v = min_xy+1; v <= max_xy; ++v) { ux*=v; }
	   PROB*= ux/uy;
           ux=uy=1.;
         }

         j=search(y,5);
         if( y[j] )  PROB/= fct( (long double)y[j]);
         return (double)PROB;
}

long double fct( long double x) {

  int t,
  tm=(int)(x + 0.5); /* round, do not truncate */
  long double rt=1.;
  if (x > 13.0 ) return (appr_fct(x)) ;
    for(t=2; t<=tm; t++) {
       rt*=(long double)t;
    }
  return rt;
}

long double appr_fct (long double i)
{
  return ( ((long double)pow(i,i)*(long double)exp(-i)*(long double)sqrt(2.*M_PI*i)));
}




