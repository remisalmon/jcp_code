#include "param.h"
#include "hex.h"
#include "fca.h"
#include "fmpi.h"
#include "fileio.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

int loadParam(int id);
void saveParam(int id, int* dims);
void loadTrace_E(int id);
//double rand_debug(void);

int main(int argc, char *argv[])
{	
	int id, nproc, nx, ny, source, dest, rc, t, op;
	double tstart, tend;
	
	srand((unsigned int)time(NULL)); //initialisation de rand()
	/*index_rand = -2;
	index_rand_tot = 0;*/
	
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);
	MPI_Comm_rank(MPI_COMM_WORLD, &id);
	
	int dims[2];
	dims[0] = (int)sqrt(nproc); //topologie carree
	dims[1] = dims[0];
	
	if(nproc < 1 || NX%dims[0] != 0 || NY%dims[1] != 0) //test topologie
	{
		MPI_Abort(MPI_COMM_WORLD, rc);
		exit(1);
	}
	
	if(id == PRINTID)
		printf("MPI %d.%d initialise avec %d proc, topologie %dx%d, tableaux %dx%d\n", MPI_VERSION, MPI_SUBVERSION, nproc, dims[0], dims[1], NX/dims[0], NY/dims[1]);
	
	if(id == PRINTID && DEBUG == 1)
		tstart = MPI_Wtime();
	
	nx = NX/dims[0]; //nbe de lignes par processeurs
	ny = NY/dims[1]; //nbe de colones par processeurs
	
	int blockcounts[2];
	blockcounts[0] = 5;
	blockcounts[1] = 5;
	MPI_Aint offsets[2], extent;
	MPI_Type_extent(MPI_DOUBLE, &extent);
	offsets[0] = 0;
	offsets[1] = 5*extent;
	MPI_Datatype oldtypes[2];
	oldtypes[0] = MPI_DOUBLE;
	oldtypes[1] = MPI_INT;
	MPI_Type_struct(2, blockcounts, offsets, oldtypes, &cellDataType);
	MPI_Type_commit(&cellDataType); //MPI data type pour structure cellData{}, definit dans hex.h
	
	MPI_Type_vector(ny, 1, 1, cellDataType, &rowB);
	MPI_Type_commit(&rowB);
	MPI_Type_vector(nx, 1, ny, cellDataType, &colB);
	MPI_Type_commit(&colB);
	MPI_Type_vector(ny, 1, 1, cellDataType, &rowA);
	MPI_Type_commit(&rowA);
	MPI_Type_vector(nx, 1, ny+2, cellDataType, &colA);
	MPI_Type_commit(&colA);
	
	cellData* A;
	cellData* B;
	A = (cellData*) malloc((nx+2)*(ny+2)*sizeof(cellData)); //tableau de cells + ghost cells pour recevoir les donnees
	B = (cellData*) malloc(nx*ny*sizeof(cellData)); //tableau pour stocker le calcul
	C = (cell*) malloc((nx+2)*(ny+2)*sizeof(cell)); //tableau de pointeurs sur voisins (dans A), definit dans hex.h
	
	topologyHex(direction, id, dims); // initialisation topologie des processus
	
	saveParam(id, dims); //ecrit mpi_param.data pour code sequentiel
	
	int init;
	init = loadParam(id); //lit CAinput.data et Energy.data
	
	// --- DEBUT CODE HEALING ---
	
	CAinitialisation(A, B, id, dims, init);
	
	if(init == 0)
	{
		loadTrace_E(id);
		CAcompute_mitosis_probability(A, B, id, dims);
	}
	
	for(t=1; t<=(1-init)*MAXT; t++) //init == 1 -> pas d'iterations
	{
		if(id == PRINTID && SICORTEX == 0)
			printf("iteration %d\n", t);
		
		CAinvade_wound(A, B, id, dims);
		
		CAcell_mobility(A, B, id, dims);
		
		if(P_MITOSIS_INSIDE > 0)
		{
			CAfind_wounded_area(A, B, id, dims);
			
			CAcompute_mitosis_probability(A, B, id, dims);
			
			CAcell_division_inside_tissue(A, B, id, dims);
		}
	}
	
	savefile(A, id, dims, 1, -1); //ecrit matrice*.out
	
	CAextraire_contour(A, B, id, dims);
	
	// --- FIN CODE HEALING ---
	
	free(C);
	free(B);
	free(A);
	
	if(id == PRINTID && DEBUG == 1)
	{
		tend = MPI_Wtime();
		printf("temps tot = %.02f sec\n", tend-tstart);
	}
	
	MPI_Finalize();
	
	return(0);
}

int loadParam(int id)
{
	int init;
	
	if(id == 0)
	{
		FILE* f;
		f=fopen("CAinput.data", "r");
		
		fscanf(f, "%d %lf", &init, &Energy_0);
		
		fclose(f);
		
		f=fopen("Energy.data", "r");
		
		int i;
		for(i=0; i<N_WOUND-2; i++)
			fscanf(f, "%lf\n", &Energy[i]);
		
		fclose(f);
	}
	
	MPI_Bcast(&init, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&Energy_0, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(Energy, N_WOUND-2, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	
	return(init);
}

void saveParam(int id, int* dims)
{
	if(id == 0)
	{
		FILE* f;
		f=fopen("mpi_param.data", "w");
		
		fprintf(f, "%d %d %d %d %d", NX, NY, dims[0], dims[1], DSTEP);
		
		fclose(f);
	}
}

void loadTrace_E(int id)
{	
	if(id == 0)
	{
		FILE* f;
		f=fopen("trace_E.out", "r");
		
		int i;
		for(i=0; i<N_WOUND-2; i++)
		{
			fscanf(f, "%lf %lf\n", &X_ansys_E[i], &Y_ansys_E[i]);
		}
		
		fclose(f);
	}
	
	MPI_Bcast(X_ansys_E, N_WOUND-2, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(Y_ansys_E, N_WOUND-2, MPI_DOUBLE, 0, MPI_COMM_WORLD);
}

/*double rand_debug(void)
{
	if(index_rand == -2)
	{
		FILE* f;
		f=fopen("random.data", "r");
		
		int i;
		for(i=0; i<10000; i++)
			fscanf(f, "%lf\n", &tab_rand[i]);
		
		fclose(f);
		
		index_rand = -1;
	}
	
	index_rand = index_rand+1;
	
	index_rand_tot = index_rand_tot+1;
	
	if(index_rand == 10000)
		index_rand = 0;
	
	return(tab_rand[index_rand]);
}*/
