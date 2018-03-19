#include "param.h"
#include "mpi.h"
#include "hex.h"
#include "fmatlab.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void initializationHex(cellData* A, cellData* B, int id, int* dims, int init) //creation des pointeurs sur voisins NW NE E SE SW S
{	
	int x, y;
	x = NX/dims[0]+2;
	y = NY/dims[1]+2;
	
	int i, j;
	for(i=1; i<x-1; i++)
	{
		for(j=1; j<y-1; j++)
		{		
			if((j+1)%2 != 0) //y odd
			{
				C[i*y+j].neighbor[0] = &A[i*y+j+1];
				C[i*y+j].neighbor[1] = &A[(i+1)*y+j+1];
				C[i*y+j].neighbor[2] = &A[(i+1)*y+j];
				C[i*y+j].neighbor[3] = &A[(i+1)*y+j-1];
				C[i*y+j].neighbor[4] = &A[i*y+j-1];
				C[i*y+j].neighbor[5] = &A[(i-1)*y+j];
			}
			else //y even
			{
				C[i*y+j].neighbor[0] = &A[(i-1)*y+j+1];
				C[i*y+j].neighbor[1] = &A[i*y+j+1];
				C[i*y+j].neighbor[2] = &A[(i+1)*y+j];
				C[i*y+j].neighbor[3] = &A[i*y+j-1];
				C[i*y+j].neighbor[4] = &A[(i-1)*y+j-1];
				C[i*y+j].neighbor[5] = &A[(i-1)*y+j];
			}
		}
	}
}

void topologyHex(int* direction, int id, int* dims) //direction = [west east north south northwest northeast southwest southeast]
{
		int nbSlaves;
		nbSlaves = dims[0]*dims[1]-1;
		
		direction[0] = id-1;
		direction[1] = id+1;
		direction[2] = id-dims[1];
		direction[3] = id+dims[1];
		direction[4] = direction[2]-1;
		direction[5] = direction[2]+1;
		direction[6] = direction[3]-1;
		direction[7] = direction[3]+1;
		
		if(direction[0] < 0 || direction[0] > nbSlaves || (id+1)%dims[1] == 1 || dims[1] == 1)
			direction[0] = MPI_PROC_NULL;
		if(direction[1] < 0 || direction[1] > nbSlaves || (id+1)%dims[1] == 0)
			direction[1] = MPI_PROC_NULL;
		if(direction[2] < 0 || direction[2] > nbSlaves)
			direction[2] = MPI_PROC_NULL;
		if(direction[3] < 0 || direction[3] > nbSlaves)
			direction[3] = MPI_PROC_NULL;
		if(direction[4] < 0 || direction[4] > nbSlaves || (id+1)%dims[1] == 1)
			direction[4] = MPI_PROC_NULL;
		if(direction[5] < 0 || direction[5] > nbSlaves || (id+1)%dims[1] == 0)
			direction[5] = MPI_PROC_NULL;
		if(direction[6] < 0 || direction[6] > nbSlaves || (id+1)%dims[1] == 1 || dims[1] == 1)
			direction[6] = MPI_PROC_NULL;
		if(direction[7] < 0 || direction[7] > nbSlaves || (id+1)%dims[1] == 0)
			direction[7] = MPI_PROC_NULL;
}

void updateHex(int i, int j, cellData* A, cellData* B, int id, int* dims, int op) //update de B a partir des valeurs de A
{
	int x, y;
	x = NX/dims[0]+2;
	y = NY/dims[1]+2;
	
	if(op == 0)
	{
		B[(i-1)*(y-2)+(j-1)].matrice = A[i*y+j].matrice;
	}
	
	if(op == 1)
	{
		double delta;
		delta = 0;
		
		int n;
		for(n=0; n<6; n++)
			delta = delta+C[i*y+j].neighbor[n]->wa;
		
		B[(i-1)*(y-2)+(j-1)].wa = (1-A[i*y+j].matrice)*(0.5*delta/6.0+0.5*A[i*y+j].wa);
		
		if(B[(i-1)*(y-2)+(j-1)].wa < 1)
			B[(i-1)*(y-2)+(j-1)].wa = 0;
	}
	
	if(op == 2)
	{
		double delta;
		delta = 0;
		
		int n;
		for(n=0; n<6; n++)
			delta = delta+C[i*y+j].neighbor[n]->wa;
		
		B[(i-1)*(y-2)+(j-1)].wa = (1-A[i*y+j].matrice)*(0.5*delta/6.0+0.5*A[i*y+j].wa);
		
		if(B[(i-1)*(y-2)+(j-1)].wa > 0)
			B[(i-1)*(y-2)+(j-1)].wa = 1;
	}
	
	if(op == 3)
	{
		double delta;
		delta = 0;
		
		int n;
		for(n=0; n<6; n++)
			delta = delta+C[i*y+j].neighbor[n]->wa;
		
		B[(i-1)*(y-2)+(j-1)].wa = 0.5*delta/6.0+0.5*A[i*y+j].wa;
		
		if(B[(i-1)*(y-2)+(j-1)].wa > 0)
			B[(i-1)*(y-2)+(j-1)].wa = 1;
	}
	
	if(op == 4)
	{	
		if(A[i*y+j].hole == 0)
		{
			double delta;
			delta = 0;
			
			int n;
			for(n=0; n<6; n++)
				delta = delta+C[i*y+j].neighbor[n]->gf;
			
			B[(i-1)*(y-2)+(j-1)].gf = 0.5*delta/6.0+0.5*A[i*y+j].gf;
		}
		else
		{
			B[(i-1)*(y-2)+(j-1)].gf = A[i*y+j].gf;
		}
	}
	
	if(op == 5)
	{
		B[(i-1)*(y-2)+(j-1)].matrice_wall = 0;
		
		int n;
		for(n=0; n<6; n++)
		{
			if(C[i*y+j].neighbor[n]->matrice == 0 && A[i*y+j].matrice == 1)
				B[(i-1)*(y-2)+(j-1)].matrice_wall = 1;
		}
	}
	
	if(op == 6)
	{
		double test;
		test = (double)rand()/(double)RAND_MAX;
		
		if(A[i*y+j].matrice_wall == 1 && test < P_MITOSIS)
		{
			double liste_dr[6];
			int liste_direction[6];
			
			int already_moved;
			
			int n;
			for(n=0; n<6; n++)
				liste_dr[n] = (double)rand()/(double)RAND_MAX;
			
			sort(liste_dr, liste_direction, 6);
			
			already_moved = 0;
			
			for(n=0; n<6; n++)
			{
				if(C[i*y+j].neighbor[liste_direction[n]]->matrice == 0 && already_moved == 0)
				{
					C[i*y+j].neighbor[liste_direction[n]]->matrice = 1;
					already_moved = 1;
				}
			}
		}
	}
	
	if(op == 7)
	{
		double delta;
		delta = 0;
		
		int n;
		for(n=0; n<6; n++)
			delta = delta+C[i*y+j].neighbor[n]->mf;
		
		B[(i-1)*(y-2)+(j-1)].mf = 0.5*delta/6.0+0.5*A[i*y+j].mf;
	}
}

void upgradeHex(cellData* A, cellData* B, int id, int* dims, int op) //recopie de B dans A
{
	int x, y;
	x = NX/dims[0]+2;
	y = NY/dims[1]+2;
	
	int i, j;
	for(i=1; i<x-1; i++) //parcours de B pour remplir A
	{
		for(j=1; j<y-1; j++) //largeur B = largeur A -2
		{
			if(op == 0)
				A[i*y+j].matrice = B[(i-1)*(y-2)+(j-1)].matrice;
			if(op == 1 || op == 2 || op == 3)
				A[i*y+j].wa = B[(i-1)*(y-2)+(j-1)].wa;
			if(op == 4)
				A[i*y+j].gf = B[(i-1)*(y-2)+(j-1)].gf;
			if(op == 5)
				A[i*y+j].matrice_wall = B[(i-1)*(y-2)+(j-1)].matrice_wall;
			if(op == 7)
				A[i*y+j].mf = B[(i-1)*(y-2)+(j-1)].mf;
			if(op == 8)
				A[i*y+j].matrice = A[i*y+j].matrice+B[(i-1)*(y-2)+(j-1)].matrice;
		}
	}
}

void updateHexBords(cellData* A, cellData* B, int id, int* dims, int op)
{
	int x, y;
	x = NX/dims[0]+2;
	y = NY/dims[1]+2;
	
	int i, j;
	for(i=1; i<x-1; i++) //calcul bords de B
	{
		j = 1;
		updateHex(i, j, A, B, id, dims, op);
		j = y-2;
		updateHex(i, j, A, B, id, dims, op);
	}
	
	for(j=1+1; j<y-1-1; j++) //calcul bords de B
	{
		i = 1;
		updateHex(i, j, A, B, id, dims, op);
		i = x-2;
		updateHex(i, j, A, B, id, dims, op);
	}
}

void updateHexInt(cellData* A, cellData* B, int id, int* dims, int op)
{
	int x, y;
	x = NX/dims[0]+2;
	y = NY/dims[1]+2;
	
	int i, j;
	for(i=1+1; i<x-1-1; i++) //calcul interieur de B
	{
		for(j=1+1; j<y-1-1; j++)
		{
			updateHex(i, j, A, B, id, dims, op);
		}
	}
}

void globalindHex(int i, int j, int id, int* dims, int* globalind)
{
	int x, y;
	x = NX/dims[0];
	y = NY/dims[1];
	
	globalind[0] = i+floor(id/dims[1])*x;
	globalind[1] = j+(id%dims[1])*y;
}
