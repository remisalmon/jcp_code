#include "hex.h"
#include "fmpi.h"

void CAinitialisation(cellData* A, cellData* B, int id, int* dims, int init);
void CAextraire_contour(cellData* A, cellData* B, int id, int* dims);
void CAfind_wounded_area(cellData* A, cellData* B, int id, int* dims);
void CAcompute_mitosis_probability(cellData* A, cellData* B, int id, int* dims);
void CAinvade_wound(cellData* A, cellData* B, int id, int* dims);
void CAcell_mobility(cellData* A, cellData* B, int id, int* dims);
void CAcell_division_inside_tissue(cellData* A, cellData* B, int id, int* dims);
int CAtest_uniformite(cellData* A, int* dims);
