/*  
          utils.h, header for utils.c.
          
Copyright (C) 1996 Christopher J. Basten.

This file is part of CND. CND is free software; you
can redistribute it and/or modify it under the terms of the GNU  General
Public License as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

CND is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
more details.

You should have received a copy of the GNU General Public License along with
CND; see the file COPYING.  If not, write to the Free Software
Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
*/
#include <string.h>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>


#if defined(BORLAND)
  #include <windows.h>
  void WYield();
#endif



extern long     ix;

#define MAXNAME 64
#define MAXLINE 1024
#define SEED 11221960
#define BIG 2147483647.0


int is_pinteger();
int myfgets();

void shuffle_ivector();

int isfile();
int get_int(); 
char get_next_line();
char get_next_token();
int dtranspose();      /*  Transpose a matrix of doubles 
                           **mm1, **mm2, lr,lc,ur,uc     */
long get_a_seed();
char *asctime2(); 
long marktime();

#ifdef DIVT
typedef struct {
  int     quot;
  int     rem;
}       div_t;

typedef struct {
  long     quot;
  long     rem;
}       ldiv_t;
div_t   div();
ldiv_t  ldiv();
#endif


extern FILE *fileopen();  /* open a file with error messages */
extern void fileclose();  /* close a file with error messages */

void pause();             /* read characters until a <CR> */
void get_field();         /* gets the specified field in a string */


double ranf();
 
int      *ivector();
long     *lvector();
float    *vector();
double   *dvector();
char     *cvector();
int     **imatrix();
short   **smatrix();
char    **cmatrix();
float   **convert_matrix();
float   **matrix();
float   **submatrix();
double  **dmatrix();
void    free_ivector();
void    free_lvector();
void    free_vector();
void    free_dvector();
void    free_cvector();
void    free_imatrix();
void    free_smatrix();
void    free_cmatrix();
void    free_dmatrix();
void    free_matrix();
void    free_submatrix();
void    free_convert_matrix();
void    nrerror();

long iran();
double gamnl1(); 
double gamgbh();

