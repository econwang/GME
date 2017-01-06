#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "stplugin.h"
extern void estimate_(double [],double (*)[],double *, double *,double *,double (*)[],double [],int *,int *,int *,int *);

STDLL stata_call(int argc, char *argv[])
{
	ST_int i,j,m,n,k,dz1,dz2,dv;
	ST_retcode rc;
	ST_double z;
	double h;
	int xcol;
	char buf[80];
	double **data;
	double *y;
	double *x;
	char entropy[] = "me";
	
	
/* import STATA DATA as C Matrix */
	n = 0;
	k = SF_nvars();
	for (j=SF_in1();j<=SF_in2(); j++) {
		if(SF_ifobs(j)) {
			n++;
		}
	}                                      /* no. of observations in the data */
	/*
	sprintf(buf,"\nThe num. of observations is: %d\n",n);
	SF_display(buf);
	*/
	data = (double **) calloc(k,sizeof(double*));     /* data matrix, with k variables and n observations */
	for (i=0;i<k;i++) {
		data[i] = (double *) calloc(n,sizeof(double));
	}
                   
	for(i=0; i<=k-1;i++) {
		m = 0;
		for(j=SF_in1(); j<=SF_in2(); j++) {
			if (SF_ifobs(j)) {
				if (rc = SF_vdata(i+1,j,&z)) return(rc);
				data[i][m] = z;
				m++;
			}
		}
	}
	
/* Import Z Matrix */
	dz1 = SF_row(argv[0]);
	dz2 = SF_col(argv[0]);
	double Zentropy[dz2][dz1];
	for (i=0;i<=dz1-1;i++) {
		for (j=0;j<=dz2-1;j++) {
		SF_mat_el(argv[0],i+1,j+1,&z);
		Zentropy[j][i] = z;
		}
	}
   
/* Import V Matrix */
	dv = SF_row(argv[1]);
	double Ventropy[dv];
	for (i=0;i<=dv-1;i++) {
		SF_mat_el(argv[1],i+1,1,&z);
		Ventropy[i] = z;
	}

/* Y data */
   y = (double *) calloc(n,sizeof(double));
   for (i=0;i<n;i++) {
		y[i] = data[0][i];
   }
   
/* X data */
   xcol = k-1;
   x = (double *) calloc(xcol*n,sizeof(double));
	for (i=0;i<xcol;i++) {
			for (j=0;j<n;j++) {
				x[i*n+j] = data[i+1][j];
			}
	}

   free(data);
/* call the fortran subroutine */
	double beta[xcol];
	double vcov[xcol][xcol];
	estimate_(beta,vcov,&h,y,x,Zentropy,Ventropy,&n,&xcol,&dz2,&dv);
	free(y);
	free(x);
/* Export Results to STATA Matrix */
	for (i=0;i<=xcol-1;i++) {
		SF_mat_store(argv[2],i+1,1,beta[i]);
		}
    for (i=0;i<xcol;i++) {
		for (j=0;j<xcol;j++) {
			SF_mat_store(argv[3],j+1,i+1,vcov[i][j]);
			sprintf(buf,"\nThe num. of observations is: %1f\n",vcov[i][j]);
			SF_display(buf);
		}
	}
	SF_scal_save(entropy,h);
	return(0);

}

