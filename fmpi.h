#include "mpi.h"
#include "hex.h"

void sendMPI(cellData* B, int* direction, MPI_Request* request, int* dims);
void recvMPI(cellData* A, int* direction, MPI_Request* request, int id, int* dims);
