#include "fmatlab.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int contour(int* tab, int nx, int ny, int* indices)
{
	int i, j;
	
	int nbindices;
	int* tmp;
	
	tmp = malloc(nx*ny*sizeof(int)); //indices non tries
	
	int tmpval;
	
	nbindices = 0;
	for(j=1; j<ny-1; j++)
	{
		for(i=nx-2; i>0; i--)
		{
			if(tab[i*ny+j] == 1)
			{
				tmpval = tab[(i-1)*ny+(j-1)]+tab[(i-1)*ny+j]+tab[(i-1)*ny+(j+1)]+tab[i*ny+(j-1)]+tab[i*ny+(j+1)]+tab[(i+1)*ny+(j-1)]+tab[(i+1)*ny+j]+tab[(i+1)*ny+(j+1)];
				
				if(tmpval < 8)
				{
					tmp[nbindices] = sub2ind(nx, ny, i, j);
					nbindices = nbindices+1;
				}
			}
		}
	}
	
	//variables tri+origine
	int i_r, nbtmp;
	double dist;
	int val_r[2];
	int val_t[2];
	
	//origine = A CHANGER
	int tmp_dist, orientation; //orientation = 0 si origine a gauche, > 0 sinon
	orientation = 0;
	ind2sub(nx, ny, tmp[0], val_t);
	tmp_dist = val_t[1];
	for(i=0; i<nbindices; i++) //cherche si origine a gauche ou a droite
	{
		ind2sub(nx, ny, tmp[i], val_t);
		if(val_t[1] < tmp_dist) //bord gauche
		{
			tmp_dist = val_t[1];
			orientation = i;
		}
		if((ny-1-val_t[1]) < tmp_dist) //bord droit
		{
			tmp_dist = ny-1-val_t[1];
			orientation = i;
		}
	}
	
	if(orientation != 0) //si origine a droite, reparcourt tab
	{
		nbindices = 0;
		for(j=ny-2; j>0; j--)
		{
			for(i=nx-2; i>0; i--)
			{	
				if(tab[i*ny+j] == 1)
				{
					tmpval = tab[(i-1)*ny+(j-1)]+tab[(i-1)*ny+j]+tab[(i-1)*ny+(j+1)]+tab[i*ny+(j-1)]+tab[i*ny+(j+1)]+tab[(i+1)*ny+(j-1)]+tab[(i+1)*ny+j]+tab[(i+1)*ny+(j+1)];
					
					if(tmpval < 8)
					{
						tmp[nbindices] = sub2ind(nx, ny, i, j);
						nbindices = nbindices+1;
					}
				}
			}
		}
	}
	
	//tri
	indices[0] = tmp[0];
	remElt(tmp, nbindices, 0);
	nbtmp = nbindices-1;
	
	while(nbtmp > 0)
	{
		dist = nx*ny;
		ind2sub(nx, ny, indices[nbindices-nbtmp-1], val_r);
		
		for(i=0; i<=nbtmp; i++)
		{
			ind2sub(nx, ny, tmp[i], val_t);
			
			if((pow(val_t[0]-val_r[0], 2)+pow(val_t[1]-val_r[1], 2)) < dist){
				indices[nbindices-nbtmp] = tmp[i];
				i_r = i;
				dist = pow(val_t[0]-val_r[0], 2)+pow(val_t[1]-val_r[1], 2);
			}
		}
		
		remElt(tmp, nbindices, i_r);
		nbtmp = nbtmp-1;
	}
	
	//index to indices i,j
	for(j=0; j<nbindices; j++)
	{
		ind2sub(nx, ny, indices[j], val_t);
		i = 0;
		indices[i*nx*ny+j] = val_t[0]+1;
		i = 1;
		indices[i*nx*ny+j] = val_t[1]+1;
	}
	
	//enleve les outliers
	tmpval = 1;
	for(i=0; i<nbindices-1; i++)
	{
		if(abs(indices[i]-indices[i+1]) > 2 || abs(indices[nx*ny+i]-indices[nx*ny+i+1]) > 2)
			break;
		
		tmpval = tmpval+1;
	}
	nbindices = tmpval;
	
	//orientation
	if(orientation == 0) //cas_gauche
	{
		if(atan2(indices[1]-indices[nbindices-1], indices[nx*ny+1]-indices[nx*ny+nbindices-1]) < 0)
		{
			nbtmp = floor((nbindices-1-1)/2);
			for(i=1; i<nbtmp; i++)
			{
				tmpval = indices[i];
				indices[i] = indices[nbindices-1-i+1];
				indices[nbindices-1-i+1] = tmpval;
				
				tmpval = indices[nx*ny+i];
				indices[nx*ny+i] = indices[nx*ny+nbindices-1-i+1];
				indices[nx*ny+nbindices-1-i+1] = tmpval;
			}
		}
	}
	else //cas_droite
	{
		if(atan2(indices[1]-indices[nbindices-1], indices[nx*ny+1]-indices[nx*ny+nbindices-1]) > 0)
		{
			nbtmp = floor((nbindices-1-1)/2);
			for(i=1; i<nbtmp; i++)
			{
				tmpval = indices[i];
				indices[i] = indices[nbindices-1-i+1];
				indices[nbindices-1-i+1] = tmpval;
				
				tmpval = indices[nx*ny+i];
				indices[nx*ny+i] = indices[nx*ny+nbindices-1-i+1];
				indices[nx*ny+nbindices-1-i+1] = tmpval;
			}
		}
	}
	
	//ajout du premier element a la fin
	i = 0;
	indices[i*nx*ny+nbindices] = indices[0];
	i = 1;
	indices[i*nx*ny+nbindices] = indices[i*nx*ny];
	
	nbindices = nbindices+1;
	
	free(tmp);
	
	return(nbindices);
}

int sub2ind(int nx, int ny, int i, int j)
{
	return(i*ny+j);
}

void ind2sub(int nx, int ny, int ind, int* val)
{
	val[0] = floor(ind/ny);
	val[1] = ind-val[0]*ny;
}

void remElt(int* tab, int size, int i2r)
{
	int i;
	for(i=i2r; i<size-1; i++)
		tab[i] = tab[i+1];
}



void interp1(double* x, double* y, double* x1, double* y1, int nx, int nx1)
{
	double alpha;
	
	int i, j;
	for(i=0; i<nx1; i++)
	{
		j = 0;
		while(1)
		{
			if(x1[i] < x[0])
			{
				y1[i] = y[0];
				break;
			}
			
			if(x1[i] <= x[j])
			{
				alpha = (double)(x1[i]-x[j-1])/(double)(x[j]-x[j-1]);
				y1[i] = (1.0-alpha)*y[j-1]+alpha*y[j];
				break;
			}
			else
				j = j+1;
			
			if(j >= nx)
			{
				y1[i] = y[nx-1];
				break;
			}
		}
	}
}



void sort(double* list, int* ordre, int size)
{
	double tmp;
	int tmpp;
	
	int i, j;
	
	for(i=0; i<size; i++)
		ordre[i] = i;
	
	for (i=size-1; i>0; i--)
	{		
		for (j=1; j<=i; j++)
		{
			if (list[j-1] > list[j])
			{
				tmp = list[j-1];
				list[j-1] = list[j];
				list[j] = tmp;
				
				tmpp = ordre[j-1];
				ordre[j-1] = ordre[j];
				ordre[j] = tmpp;
			}
		}
	}
}
