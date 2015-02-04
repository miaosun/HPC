/*
 * Jacobi.c
 *
 *  Created on: Feb 1, 2015
 *      Author: miaosun
 */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#define PI 3.14159265

#define max 20                         /* number of grid points */
#define TOLERANCE 0.00000001
#define ITERATION 100000000

double average(double p[max][max])
{
	double res = 0;
	for(int i=0; i<max; i++)
		for(int j=0; j<max; j++)
			res += p[i][j];

	return res/(max*max);
}

double rms_norm (double a[max][max] )
{
	int i;
	int j;
	double v;

	v = 0.0;

	for ( j = 0; j < max; j++ )
	{
		for ( i = 0; i < max; i++ )
		{
			v = v + a[i][j] * a[i][j];
		}
	}
	v = sqrt ( v / ( double ) ( max * max )  );

	return v;
}

double analytical_norm(int iter)
{
	double a[max][max];

	double h = 1.0 / (max-1);
	double res = 0.0;
	long double b = 0.0;
	long double c = 0.0;
	for(int i=0; i<max; i++)
	{
		for(int j=0; j<max; j++)
		{
			double x = i*h;
			double y = j*h;
			for(int k=1; k<=iter; k++)
			{
				if(k != 2)
				{
					b = 4*(cos(k*PI)-1)*(1/tanh(k*PI))*sin(k*PI*y)*(tanh(k*PI)*cosh(k*PI*x)-sinh(k*PI*x));
					c = (PI*(pow(k,3)-4*k));

					res += b/c;
				}
			}

			a[i][j] = res;
		}
	}

	return rms_norm(a);
}

int main()
{
	double p[max][max];
	int i, j, iter;
	time_t time1, time2;
	double seconds;

	FILE *output;			/* save data in laplace.dat */
	output = fopen("laplace.dat","w");

	for(i=0; i<max; i++)                 /* clear the array  */
	{
		for (j=0; j<max; j++)
			p[i][j] = 0;
	}

	for(i=1; i<max-1; i++)
	{
		p[i][0] = sin(PI*i/max) * sin(PI*i/max);
	}

	double old_norm = 0.0;
	double norm = 0.0;
	double analytical = 0.0;

	time(&time1);

	for(iter=1; iter<ITERATION; iter++)               /* iterations */
	{
		old_norm = norm;
		for(i=1; i<(max-1); i++)                  /* x-direction */
		{
			for(j=1; j<(max-1); j++)               /* y-direction */
			{
				p[i][j] = 0.25*(p[i+1][j]+p[i-1][j]+p[i][j+1]+p[i][j-1]);
			}

		}
		for(int i=0; i<max; i++)
		{
			for(int j=0; j<max; j++)
				printf("%.2f\t", p[i][j]);
			printf("\n");
		}
		norm = rms_norm(p);

		if(fabs(norm-old_norm) < TOLERANCE)
		{
			printf("\nConverged after %d iterations\n", iter);

			time(&time2);
			seconds = difftime(time2, time1);
			printf("Time spent: %.2f seconds\n", seconds);

			analytical = analytical_norm(iter);
			printf("\nError: %f\n", fabs(analytical-norm));
			break;
		}
		printf("\n");
	}



	for (i=0; i<max ; i++)         /* write data gnuplot 3D format */
	{
		for (j=0; j<max; j++)
		{
			fprintf(output, "%.2f\t",p[i][j]);
		}
		fprintf(output, "\n");	  /* empty line for gnuplot */
	}
	printf("data stored in laplace.dat\n");
	fclose(output);
}
