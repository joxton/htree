#include <R.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>



void revsort(double *a, int *ib, int n)
{
/* Sort a[] into descending order by "heapsort";
 * sort ib[] alongside;
 * if initially, ib[] = 1...n, it will contain the permutation finally
 */

    int l, j, ir, i;
    double ra;
    int ii;

    if (n <= 1) return;

    a--; ib--;

    l = (n >> 1) + 1;
    ir = n;

    for (;;) {
	if (l > 1) {
	    l = l - 1;
	    ra = a[l];
	    ii = ib[l];
	}
	else {
	    ra = a[ir];
	    ii = ib[ir];
	    a[ir] = a[1];
	    ib[ir] = ib[1];
	    if (--ir == 1) {
		a[1] = ra;
		ib[1] = ii;
		return;
	    }
	}
	i = l;
	j = l << 1;
	while (j <= ir) {
	    if (j < ir && a[j] > a[j + 1]) ++j;
	    if (ra > a[j]) {
		a[i] = a[j];
		ib[i] = ib[j];
		j += (i = j);
	    }
	    else
		j = ir + 1;
	}
	a[i] = ra;
	ib[i] = ii;
    }
}




void sort_test(int *niter,double *a,int *n,double *res,int *b)
{

  // void revsort(double *a, int *ib, int n)

  double *ax;
  int *ib;
  int k,i;

  ib=(int*) malloc(sizeof(int)*n[0]);
  ax=(double*) malloc(sizeof(double)*n[0]);

  for(k=0;k<niter[0];k++)
    {

      for(i=0;i<n[0];i++)
	{
	  ax[i]=a[i];
	  ib[i]=i;
	}
      revsort(ax,ib,n[0]);

    }

  for(i=0;i<n[0];i++)
    {
      res[i]=ax[i];
      b[i]=ib[i];
    }

  
}







