#include "param.h"
#include "fmpi.h"

void sendMPI(cellData* B, int* direction, MPI_Request* request, int* dims)
{
	int x, y;
	x = NX/dims[0];
	y = NY/dims[1];
	
	MPI_Isend(&B[0], 1, colB, direction[0], 1, MPI_COMM_WORLD, &request[0]); //envoie W
	
	MPI_Isend(&B[y-1], 1, colB, direction[1], 1, MPI_COMM_WORLD, &request[1]); //envoie E
	
	MPI_Isend(&B[0], 1, rowB, direction[2], 1, MPI_COMM_WORLD, &request[2]); //envoie N
	
	MPI_Isend(&B[(x-1)*y], 1, rowB, direction[3], 1, MPI_COMM_WORLD, &request[3]); //envoie S
	
	MPI_Isend(&B[0], 1, cellDataType, direction[4], 1, MPI_COMM_WORLD, &request[4]); //envoie NW
	
	MPI_Isend(&B[y-1], 1, cellDataType, direction[5], 1, MPI_COMM_WORLD, &request[5]); //envoie NE
	
	MPI_Isend(&B[(x-1)*y], 1, cellDataType, direction[6], 1, MPI_COMM_WORLD, &request[6]); //envoie SW
	
	MPI_Isend(&B[x*y-1], 1, cellDataType, direction[7], 1, MPI_COMM_WORLD, &request[7]); //envoie SE
}

void recvMPI(cellData* A, int* direction, MPI_Request* request, int id, int* dims)
{
	int x, y;
	x = NX/dims[0];
	y = NY/dims[1];
	
	MPI_Recv(&A[1*(y+2)], 1, colA, direction[0], 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE); //reception W
	
	MPI_Recv(&A[1*(y+2)+y+2-1], 1, colA, direction[1], 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE); //reception E
	
	MPI_Recv(&A[1], 1, rowA, direction[2], 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE); //reception N
	
	MPI_Recv(&A[(x+2-1)*(y+2)+1], 1, rowA, direction[3], 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE); //reception S
	
	MPI_Recv(&A[0], 1, cellDataType, direction[4], 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE); //reception NW
	
	MPI_Recv(&A[y+2-1], 1, cellDataType, direction[5], 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE); //reception NE
	
	MPI_Recv(&A[(x+2-1)*(y+2)], 1, cellDataType, direction[6], 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE); //reception SW
	
	MPI_Recv(&A[(y+2)*(x+2)-1], 1, cellDataType, direction[7], 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE); //reception SE
	
	MPI_Waitall(8, request, MPI_STATUSES_IGNORE);
}
