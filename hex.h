#include "mpi.h"

#ifndef HEX_H
#define HEX_H

typedef struct cellData //structure pour MPI (A, B)
{
	double wa;
	double oldwa;
	double gf;
	double matrice_mitosis;
	double mf;
	int matrice;
	int matrice_wall;
	int hole;
	int hole_etendu;
	int mi;
} cellData;

typedef struct cell //structure cellule de la grille (C)
{
	struct cellData* neighbor[6]; //pointeur sur voisins NW NE E SE SW W
} cell;

cell* C;

MPI_Datatype cellDataType;
MPI_Datatype rowB;
MPI_Datatype colB;
MPI_Datatype rowA;
MPI_Datatype colA;

MPI_Request request[8];
int direction[8]; //direction = [west east north south northwest northeast southwest southeast]

#endif

void initializationHex(cellData* A, cellData* B, int id, int* dims, int init); //initialisation de la grille
void topologyHex(int* direction, int id, int* dims); //direction = [west east north south northwest northeast southwest southeast]
void updateHex(int i, int j, cellData* A, cellData* B, int id, int* dims, int op); //update de B a partir des valeurs de A
void upgradeHex(cellData* A, cellData* B, int id, int* dims, int op); //recopie de B dans A
void updateHexBords(cellData* A, cellData* B, int id, int* dims, int op);
void updateHexInt(cellData* A, cellData* B, int id, int* dims, int op);
void globalindHex(int i, int j, int id, int* dims, int* globalind);
