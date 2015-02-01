/*
 * Jacobi.c
 *
 *  Created on: Feb 1, 2015
 *      Author: miaosun
 */


#include <stdio.h>
#define PI 3.14159265

#define max 20                         /* number of grid points */

double average(double p[max][max])
{
	double res = 0;
	for(int i=0; i<max; i++)
		for(int j=0; j<max; j++)
			res += p[i][j];

	return res/(max*max);
}

main()
{
	double x, p[max][max];
	int i, j, iter, y;

	FILE *output;			/* save data in laplace.dat */
	output = fopen("laplace.dat","w");

	for(i=0; i<max; i++)                 /* clear the array  */
	{
		for (j=0; j<max; j++)
			p[i][j] = 0;
	}

	for(i=1; i<max-1; i++)
		p[i][0] = sin(PI*1/i) * sin(PI*1/i);

	double old_avg = 0.0;
	double avg = 0.0;
	for(iter=0; iter<10000; iter++)               /* iterations */
	{
		old_avg = avg;
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
		avg = average(p);

		if(avg-old_avg < 0.00000001)
		{
			printf("\nConverged after %d iterations\n", iter);
			exit(0);
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
