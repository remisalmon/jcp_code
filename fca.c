#include "param.h"
#include "fileio.h"
#include "fca.h"
#include "fmatlab.h"
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>

void CAinitialisation(cellData* A, cellData* B, int id, int* dims, int init)
{
	int x, y;
	x = NX/dims[0]+2;
	y = NY/dims[1]+2;
	
	int i, j;
	if(init == 1)
	{
		for(i=0; i<x; i++)
		{
			for(j=0; j<y; j++)
			{
				A[i*y+j].wa = 0;
				A[i*y+j].oldwa = 0;
				A[i*y+j].gf = 0;
				A[i*y+j].matrice_mitosis = 0;
				A[i*y+j].mf = 0;
				A[i*y+j].matrice = 1;
				A[i*y+j].matrice_wall = 0;
				A[i*y+j].hole = 0;
				A[i*y+j].hole_etendu = 0;
				A[i*y+j].mi = 0;
				
				if(i>0 && i<x-1 && j>0 && j<y-1)
				{
					double dx, dy;
					dx = (double)LX/(NX-1);
					dy = (double)LY/(NY-1);
					
					int globalind[2];
					globalindHex(i, j, id, dims, globalind);
					
					double rotation[4];
					rotation[0] = cos(PHI_ELL*PI/180);
					rotation[1] = -sin(PHI_ELL*PI/180);
					rotation[2] = sin(PHI_ELL*PI/180);
					rotation[3] = cos(PHI_ELL*PI/180);
					
					double xg, yg;
					xg = globalind[0]*dx-LX/2.0;
					yg = globalind[1]*dy-LY/2.0;
					
					double theta;
					theta = atan2(yg, xg);
					
					double distance;
					distance = sqrt(pow(xg, 2)+pow(yg, 2));
					
					double xt, yt;
					xt = cos(theta);
					yt = sin(theta);
					
					double cr[2];
					cr[0] = rotation[0]*xt+rotation[1]*yt;
					cr[1] = rotation[2]*xt+rotation[3]*yt;
					
					double cost, sint;
					cost = cos(theta+PHI_ELL*PI/180);
					sint = sin(theta+PHI_ELL*PI/180);
					
					double r;
					if(WOUND_OPTION_ELL == 1)
						r = (A_ELL*B_ELL)/sqrt(pow(B_ELL*cost, 2)+pow(A_ELL*sint, 2));
					else
						r = A_ELL*fabs(cr[0])+B_ELL*fabs(cr[1]);
					
					if(distance > r)
						A[i*y+j].matrice = 1;
					else
						A[i*y+j].matrice = 0;
				}
			}
		}
	}
	else
	{
		readfile(A, id, dims);
		
		MPI_Barrier(MPI_COMM_WORLD);
		
		updateHexBords(A, B, id, dims, 0);
		
		sendMPI(B, direction, request, dims);
		
		updateHexInt(A, B, id, dims, 0);
		
		recvMPI(A, direction, request, id, dims);
		
		upgradeHex(A, B, id, dims, 0);
	}
	
	initializationHex(A, B, id, dims, init);
	
	CAfind_wounded_area(A, B, id, dims);
}

void CAextraire_contour(cellData* A, cellData* B, int id, int* dims)
{
	CAfind_wounded_area(A, B, id, dims);
	
	savefile(A, id, dims, 2, -1); //ecrit hole*.txt
	
	int x, y;
	x = NX/dims[0]+2;
	y = NY/dims[1]+2;
	
	int i, j;
	for(i=0; i<x; i++)
	{
		for(j=0; j<y; j++)
		{
			A[i*y+j].wa = A[i*y+j].hole;
		}
	}
	
	int kt;
	
	for(kt=1; kt<=round(MAX_DIFFUSION_GROWTH_FACTOR/2); kt++)
	{
		updateHexBords(A, B, id, dims, 3);
		
		sendMPI(B, direction, request, dims);
		
		updateHexInt(A, B, id, dims, 3);
		
		recvMPI(A, direction, request, id, dims);
		
		upgradeHex(A, B, id, dims, 3);
	}
	
	for(i=0; i<x; i++)
	{
		for(j=0; j<y; j++)
		{
			A[i*y+j].hole_etendu = (int)A[i*y+j].wa;
		}
	}
	
	savefile(A, id, dims, 3, -1); //ecrit hole_etendu*.txt
	
	MPI_Barrier(MPI_COMM_WORLD); //debut code sequentiel
	
	double tstart_seq, tend_seq;
	if(id == PRINTID && DEBUG == 1)
		tstart_seq = MPI_Wtime();
	
	if(id == 0)
	{
		if(SICORTEX == 0)
		{
			system("cd /home/remi/Desktop/CA/code/code_seq && ./code_seq"); //execute code sequentiel
		}
		
		if(SICORTEX == 1)
		{
			system("tar -cf data_hole.tar hole*.txt && scp -q data_hole.tar mpi_param.data remi:/home/remi/Desktop/CA/code/"); //envoie donnees
			system("ssh remi 'cd /home/remi/Desktop/CA/code/ && tar -xf data_hole.tar && rm data_hole.tar && cd code_seq/ && ./code_seq'"); //execute code sequentiel
			system("rm hole*.txt data_hole.tar 2>/dev/null"); //donnees inutiles
			system("scp -q remi:/home/remi/Desktop/CA/code/*.txt ."); //recupere coordonnees ANSYS
			system("scp -q remi:/home/remi/Desktop/CA/code/trace_E.out ."); //recupere trace_E
		}
	}
	
	MPI_Barrier(MPI_COMM_WORLD); //fin code sequentiel
	
	if(id == PRINTID && DEBUG == 1)
	{
		tend_seq = MPI_Wtime();
		printf("temps seq = %.02f sec\n", tend_seq-tstart_seq);
	}
}

void CAfind_wounded_area(cellData* A, cellData* B, int id, int* dims)
{	
	int x, y;
	x = NX/dims[0]+2;
	y = NY/dims[1]+2;
	
	int i, j;
	for(i=0; i<x; i++)
	{
		for(j=0; j<y; j++)
		{
			A[i*y+j].wa = 1-A[i*y+j].matrice;
		}
	}
	
	int kt;
	for(kt=1; kt<=1; kt++)
	{
		updateHexBords(A, B, id, dims, 1);
		
		sendMPI(B, direction, request, dims);
		
		updateHexInt(A, B, id, dims, 1);
		
		recvMPI(A, direction, request, id, dims);
		
		upgradeHex(A, B, id, dims, 1);
	}
	
	for(i=0; i<x; i++)
	{
		for(j=0; j<y; j++)
		{
			A[i*y+j].oldwa = A[i*y+j].wa;	
		}
	}
	
	int test, test_local, test_global, test_unif;
	double max_local;
	
	test = 1;
	
	while(test > 0)
	{
		updateHexBords(A, B, id, dims, 2);
		
		sendMPI(B, direction, request, dims);
		
		updateHexInt(A, B, id, dims, 2);
		
		recvMPI(A, direction, request, id, dims);
		
		upgradeHex(A, B, id, dims, 2);
		
		test_unif = CAtest_uniformite(A, dims);
		
		if(test_unif == 0)
		{
			max_local = 0;
			for(i=1; i<x-1; i++)
			{
				for(j=1; j<y-1; j++)
				{
					if(fabs(A[i*y+j].wa-A[i*y+j].oldwa) > max_local)
						max_local = fabs(A[i*y+j].wa-A[i*y+j].oldwa);
				}
			}
			if(max_local == 0)
				test_local = 0;
			else
				test_local = 1;
		}
		else
		{
			test_local = 0;
		}
		
		for(i=0; i<x; i++)
		{
			for(j=0; j<y; j++)
			{
				A[i*y+j].oldwa = A[i*y+j].wa;
			}
		}
		
		test_global = 0;
		
		MPI_Allreduce(&test_local, &test_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
		
		if(test_global == 0)
			test = 0;
	}
	
	for(i=0; i<x; i++)
	{
		for(j=0; j<y; j++)
		{
			A[i*y+j].hole = (int)A[i*y+j].wa;
		}
	}
}

void CAcompute_mitosis_probability(cellData* A, cellData* B, int id, int* dims)
{
	int x, y;
	x = NX/dims[0]+2;
	y = NY/dims[1]+2;

	int i, j;	
	for(i=0; i<x; i++)
	{
		for(j=0; j<y; j++)
		{
			A[i*y+j].gf = A[i*y+j].hole;
		}
	}
	
	int kt_diffusion;
	for(kt_diffusion=0; kt_diffusion<MAX_DIFFUSION_GROWTH_FACTOR; kt_diffusion++)
	{
		updateHexBords(A, B, id, dims, 4);
		
		sendMPI(B, direction, request, dims);
		
		updateHexInt(A, B, id, dims, 4);
		
		recvMPI(A, direction, request, id, dims);
		
		upgradeHex(A, B, id, dims, 4);
	}
	
	if(RATE_OF_DECAY_OF_GF == 0)
	{
		for(i=0; i<x; i++)
		{
			for(j=0; j<y; j++)
			{
				if(A[i*y+j].gf > 0)
					A[i*y+j].gf = 1;
			}
		}
	}
	
	for(i=0; i<x; i++)
	{
		for(j=0; j<y; j++)
		{
			A[i*y+j].matrice_wall = 0;
		}
	}
	
	updateHexBords(A, B, id, dims, 5);
	
	sendMPI(B, direction, request, dims);
	
	updateHexInt(A, B, id, dims, 5);
	
	recvMPI(A, direction, request, id, dims);
	
	upgradeHex(A, B, id, dims, 5);
	
	double dx, dy;
	dx = (double)LX/(NX-1);
	dy = (double)LY/(NY-1);
	
	int globalind[2];
	double xg, yg;
	
	int jE, JE;
	double distance_E, E_theta;
	
	for(i=0; i<x; i++)
	{
		for(j=0; j<y; j++)
		{
			A[i*y+j].matrice_mitosis = 0;
			
			if((A[i*y+j].hole+A[i*y+j].matrice_wall) == 0)
			{
				if(A[i*y+j].gf > 0)
				{
					globalindHex(i, j, id, dims, globalind);
					xg = globalind[0]*dx-LX/2.0;
					yg = globalind[1]*dy-LY/2.0;
					
					xg = ANSYS_SCALING*xg;
					yg = ANSYS_SCALING*yg;
					
					distance_E = pow(xg-X_ansys_E[0], 2)+pow(yg-Y_ansys_E[0], 2);
					JE = 0;
					
					for(jE=0; jE<N_WOUND-2; jE++)
					{
						if((pow(xg-X_ansys_E[jE], 2)+pow(yg-Y_ansys_E[jE], 2)) < distance_E)
						{
							distance_E = pow(xg-X_ansys_E[jE], 2)+pow(yg-Y_ansys_E[jE], 2);
							JE = jE;
						}
					}
					
					E_theta = Energy[JE];
					
					A[i*y+j].matrice_mitosis = P_MITOSIS_INSIDE*(COEF_E0+COEF_E1*E_theta/Energy_0)*A[i*y+j].gf;
				}
			}
		}
	}
}

void CAinvade_wound(cellData* A, cellData* B, int id, int* dims)
{
	int x, y;
	x = NX/dims[0]+2;
	y = NY/dims[1]+2;
	
	double liste_ri[x-2];
	double liste_rj[y-2];
	int liste_i[x-2];
	int liste_j[y-2];
	
	int i, j, i2, j2;
	for(i=0; i<x; i++)
	{
		for(j=0; j<y; j++)
		{
			A[i*y+j].matrice_wall = 0;
		}
	}
	
	updateHexBords(A, B, id, dims, 5);
	
	sendMPI(B, direction, request, dims);
	
	updateHexInt(A, B, id, dims, 5);
	
	recvMPI(A, direction, request, id, dims);
	
	upgradeHex(A, B, id, dims, 5);
	
	/*for(j=0; j<y-2; j++)
		liste_rj[j] = (double)rand()/(double)RAND_MAX;
	
	sort(liste_rj, liste_j, y-2);
	
	for(j2=1; j2<y-1; j2++) //calcul bords de B
	{
		j = liste_j[j2-1]+1;
		
		i = 1;
		updateHex(i, j, A, B, id, dims, 6);
		i = x-2;
		updateHex(i, j, A, B, id, dims, 6);
	}
	
	for(i=0; i<x-4; i++)
		liste_ri[i] = (double)rand()/(double)RAND_MAX;
	
	sort(liste_ri, liste_i, x-4);
	
	for(i2=1+1; i2<x-1-1; i2++) //calcul bords de B
	{
		i = liste_i[i2-2]+2;
		
		j = 1;
		updateHex(i, j, A, B, id, dims, 6);
		j = y-2;
		updateHex(i, j, A, B, id, dims, 6);
	}
	
	updateHexBords(A, B, id, dims, 0);
	
	sendMPI(B, direction, request, dims);
	
	for(i=0; i<x-4; i++)
		liste_ri[i] = (double)rand()/(double)RAND_MAX;
	
	sort(liste_ri, liste_i, x-4);
	
	for(i2=1+1; i2<x-1-1; i2++) //calcul interieur de B
	{
		i = liste_i[i2-2]+2;
		
		for(j=0; j<y-4; j++)
			liste_rj[j] = (double)rand()/(double)RAND_MAX;
		
		sort(liste_rj, liste_j, y-4);
		
		for(j2=1+1; j2<y-1-1; j2++)
		{
			j = liste_j[j2-2]+2;
			
			updateHex(i, j, A, B, id, dims, 6);
		}
	}
	
	recvMPI(A, direction, request, id, dims);
	
	updateHexInt(A, B, id, dims, 0);
	
	upgradeHex(A, B, id, dims, 0);*/
	
	for(i=1; i<x-1; i++)
		liste_ri[i] = (double)rand()/(double)RAND_MAX;
	
	sort(liste_ri, liste_i, x-2);
	
	for(i2=1; i2<x-1; i2++)
	{
		i = liste_i[i2-1]+1;
		
		for(j=1; j<y-1; j++)
			liste_rj[j] = (double)rand()/(double)RAND_MAX;
		
		sort(liste_rj, liste_j, y-2);
		
		for(j2=1; j2<y-1; j2++)
		{
			j = liste_j[j2-1]+1;
			
			updateHex(i, j, A, B, id, dims, 6);
		}
	}
}

void CAcell_mobility(cellData* A, cellData* B, int id, int* dims)
{	
	int x, y;
	x = NX/dims[0]+2;
	y = NY/dims[1]+2;
	
	int i, j;
	for(i=0; i<x; i++)
	{
		for(j=0; j<y; j++)
		{
			A[i*y+j].mf = A[i*y+j].matrice;
		}
	}
	
	int kt_diffusion;
	for(kt_diffusion=0; kt_diffusion<MAX_DIFFUSION_MOBILITY; kt_diffusion++)
	{
		updateHexBords(A, B, id, dims, 7);
		
		sendMPI(B, direction, request, dims);
		
		updateHexInt(A, B, id, dims, 7);
		
		recvMPI(A, direction, request, id, dims);
		
		upgradeHex(A, B, id, dims, 7);
	}
	
	int masse0, masse0_local, masse, masse_local, test_dichotomie, kd;
	double tolA, tolB, tol;
	
	masse0 = 0;
	masse0_local = 0;
	
	for(i=1; i<x-1; i++)
	{
		for(j=1; j<y-1; j++)
		{
			masse0_local = masse0_local+A[i*y+j].matrice;
		}
	}
	
	MPI_Allreduce(&masse0_local, &masse0, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	
	tolA = 0;
	tolB = 1;
	
	test_dichotomie = NX*NY;
	kd = 0;
	
	while(test_dichotomie != 0 && kd < 100)
	{
		kd = kd+1;
		
		tol = 0.5*(tolA+tolB);
		
		for(i=1; i<x-1; i++)
		{
			for(j=1; j<y-1; j++)
			{				
				if(A[i*y+j].mf >= tol)
					B[(i-1)*(y-2)+(j-1)].matrice = 1;
				else
					B[(i-1)*(y-2)+(j-1)].matrice = 0;
			}
		}
		
		masse = 0;
		masse_local = 0;
		for(i=1; i<x-1; i++)
		{
			for(j=1; j<y-1; j++)
			{
				masse_local = masse_local+B[(i-1)*(y-2)+(j-1)].matrice;
			}
		}
		
		MPI_Allreduce(&masse_local, &masse, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
		
		if(masse > masse0)
			tolA = tol;
		else
			tolB = tol;
		
		test_dichotomie = masse-masse0;
	}
	
	upgradeHex(A, B, id, dims, 0);
}

void CAcell_division_inside_tissue(cellData* A, cellData* B, int id, int* dims)
{
	int x, y;
	x = NX/dims[0]+2;
	y = NY/dims[1]+2;
	
	int i, j;
	for(i=0; i<x; i++)
	{
		for(j=0; j<y; j++)
		{
			A[i*y+j].mi = 0;
			A[i*y+j].mf = 0;
		}
	}
	
	double p_bump;
	
	for(i=1; i<x-1; i++)
	{
		for(j=1; j<y-1; j++)
		{
			if((A[i*y+j].hole+A[i*y+j].matrice_wall) == 0)
			{
				p_bump = (double)rand()/(double)RAND_MAX;
				
				if(p_bump < A[i*y+j].matrice_mitosis)
				{
					A[i*y+j].mi = 1;
					A[i*y+j].mf = 1;
				}
			}
		}
	}
	
	int kt_diffusion;
	for(kt_diffusion=0; kt_diffusion<MIN(NX, NY)/2; kt_diffusion++)
	{
		updateHexBords(A, B, id, dims, 7);
		
		sendMPI(B, direction, request, dims);
		
		updateHexInt(A, B, id, dims, 7);
		
		recvMPI(A, direction, request, id, dims);
		
		upgradeHex(A, B, id, dims, 7);
	}
	
	for(i=0; i<x; i++)
	{
		for(j=0; j<y; j++)
		{
			A[i*y+j].mf = A[i*y+j].mf*A[i*y+j].hole;
		}
	}
	
	double maxMF_local, maxMF;
	
	maxMF_local = 0;
	
	for(i=1; i<x-1; i++)
	{
		for(j=1; j<y-1; j++)
		{
			if(A[i*y+j].mf > maxMF_local)
				maxMF_local = A[i*y+j].mf;
		}
	}
	
	MPI_Allreduce(&maxMF_local, &maxMF, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
	
	for(i=0; i<x; i++)
	{
		for(j=0; j<y; j++)
		{
			A[i*y+j].mf = A[i*y+j].mf/maxMF+0.001*((double)rand()/(double)RAND_MAX)*A[i*y+j].hole;
		}
	}
	
	maxMF_local = 0;
	
	for(i=1; i<x-1; i++)
	{
		for(j=1; j<y-1; j++)
		{
			if(A[i*y+j].mf > maxMF_local)
				maxMF_local = A[i*y+j].mf;
		}
	}
	
	MPI_Allreduce(&maxMF_local, &maxMF, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
	
	int masse0, masse0_local, masse, masse_local, test_dichotomie, kd;
	double tolA, tolB, tol;
	
	masse0 = 0;
	masse0_local = 0;
	
	for(i=1; i<x-1; i++)
	{
		for(j=1; j<y-1; j++)
		{
			masse0_local = masse0_local+A[i*y+j].mi;
		}
	}
	
	MPI_Allreduce(&masse0_local, &masse0, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	
	tolA = 0;
	tolB = 1.1*maxMF;
	
	test_dichotomie = NX*NY;
	kd = 0;
	
	while(test_dichotomie != 0 && kd < 100)
	{
		kd = kd+1;
		
		tol = 0.5*(tolA+tolB);
		
		for(i=1; i<x-1; i++)
		{
			for(j=1; j<y-1; j++)
			{		
				if(A[i*y+j].mf >= tol)
					B[(i-1)*(y-2)+(j-1)].matrice = 1;
				else
					B[(i-1)*(y-2)+(j-1)].matrice = 0;
			}
		}
		
		masse = 0;
		masse_local = 0;
		for(i=1; i<x-1; i++)
		{
			for(j=1; j<y-1; j++)
			{
				masse_local = masse_local+B[(i-1)*(y-2)+(j-1)].matrice;
			}
		}
		
		MPI_Allreduce(&masse_local, &masse, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
		
		if(masse > masse0)
			tolA = tol;
		else
			tolB = tol;
		
		test_dichotomie = masse-masse0;
	}
	
	upgradeHex(A, B, id, dims, 8);
}

int CAtest_uniformite(cellData* A, int* dims)
{	
	int x, y;
	x = NX/dims[0]+2;
	y = NY/dims[1]+2;
	
	int test;
	test = 0;
	
	int i, j;
	for(i=1; i<x-1; i++)
	{
		for(j=1; j<y-1; j++)
		{
			test = test+A[i*y+j].matrice;
		}
	}
	
	if(test == (x-2)*(y-2) || test == 0)
		test = 1;
	else
		test = 0;
	
	return(test);
}
