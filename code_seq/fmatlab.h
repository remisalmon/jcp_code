int contour(int* tab, int nx, int ny, int* indices);
int sub2ind(int nx, int ny, int i, int j);
void ind2sub(int nx, int ny, int ind, int* val);
void remElt(int* tab, int size, int i2r);

void interp1(double* x, double* y, double* x1, double* y1, int nx, int nx1);

void sort(double* list, int* ordre, int size);
