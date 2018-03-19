#include "param.h"
#include "fileio.h"
#include <stdio.h>
#include <stdlib.h>

void savefile(cellData* A, int id, int* dims, int val, int iter) //sauvegarde de data au format data_mpi[n].txt
{	
	FILE* f;
	char name[20];
	
	if(val == 1)
		sprintf(name, "matrice%d.out", id+1);
	if(val == 2)
		sprintf(name, "hole%d.txt", id+1);
	if(val == 3)
		sprintf(name, "hole_etendu%d.txt", id+1);
	if(val == 4)
		sprintf(name, "matrice_wall%d.txt", id+1);
	if(val == 5)
		sprintf(name, "matrice%d_iter%d.data", id+1, iter);
	
	f=fopen(name, "w");
	
	int x, y;
	x = NX/dims[0]+2;
	y = NY/dims[1]+2;
	
	int i, j;
	for(i=1; i<x-1; i++)
	{
		for(j=1; j<y-1; j++)
		{
			if(val == 1 || val == 5)
				fprintf(f, "%d ", A[i*y+j].matrice);
			if(val == 2)
				fprintf(f, "%d ", A[i*y+j].hole);
			if(val == 3)
				fprintf(f, "%d ", A[i*y+j].hole_etendu);
			if(val == 4)
				fprintf(f, "%d ", A[i*y+j].matrice_wall);
		}
		fprintf(f, "\n");
	}
	
	fclose(f);
}

void readfile(cellData* A, int id, int* dims) //sauvegarde de data au format data_mpi[n].txt
{	
	FILE* f;
	char name[20];
	
	sprintf(name, "matrice%d.out", id+1);
	f=fopen(name, "r");
	
	int x, y;
	x = NX/dims[0]+2;
	y = NY/dims[1]+2;
	
	int i, j;
	for(i=0; i<x; i++)
	{
		for(j=0; j<y; j++)
		{
			A[i*y+j].matrice = 1;
			
			if(i>0 && i<x-1 && j>0 && j<y-1)
				fscanf(f, "%d", &A[i*y+j].matrice);	
		}
	}
	
	fclose(f);
}
