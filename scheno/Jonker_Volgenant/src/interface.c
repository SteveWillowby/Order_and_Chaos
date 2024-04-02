#include "assignAlgs2D.h"
#include <time.h>
#include <stdlib.h>
#include <stddef.h>
#include <stdio.h>

bool bipartite_assignment_with_duals(const bool maximize, double *  C, double *gain, int * colForRow, int *  rowForCol, double *  u, double *  v, const size_t numRow, const size_t numCol) {
    ptrdiff_t * col4row;
    ptrdiff_t * row4col;
    void *tempBuffer;

    tempBuffer=malloc(assign2DSimpBufferSize(numRow,numCol)+(numRow+numCol)*sizeof(ptrdiff_t));
    col4row=(ptrdiff_t*)tempBuffer;
    row4col=col4row+numRow;
    void* bufferStart=(void*)(row4col+numCol);

    for(size_t i=0; i<numCol; ++i){
        u[i] = 0.0;
    }
    for(size_t i=0; i<numRow; ++i){
        v[i] = 0.0;
    }

    bool feasible = assign2DC(maximize, C, gain, col4row, row4col, bufferStart, u, v, numRow, numCol);

    for (size_t i=0; i<numRow; ++i) {
        colForRow[i] = col4row[i];
    }
    for (size_t i=0; i<numCol; ++i) {
        rowForCol[i] = row4col[i];
    }
    free(tempBuffer);
    return(feasible);
}

bool bipartite_assignment(const bool maximize, double *  C, double *gain, int * colForRow, int *  rowForCol, const size_t numRow, const size_t numCol) {
    ptrdiff_t * col4row;
    ptrdiff_t * row4col;
    double *  u;
    double *  v;
    void *tempBuffer;

    tempBuffer=malloc(assign2DSimpBufferSize(numRow,numCol)+(numRow+numCol)*sizeof(ptrdiff_t));
    col4row=(ptrdiff_t*)tempBuffer;
    row4col=col4row+numRow;
    void* bufferStart=(void*)(row4col+numCol);
    u = (double *) malloc(numCol * sizeof(double));
    v = (double *) malloc(numRow * sizeof(double));
    for(size_t i=0; i<numCol; ++i){
        u[i] = 0.0;
    }
    for(size_t i=0; i<numRow; ++i){
        v[i] = 0.0;
    }

    bool feasible = assign2DC(maximize, C, gain, col4row, row4col, bufferStart, u, v, numRow, numCol);

    for (size_t i=0; i<numRow; ++i) {
        colForRow[i] = col4row[i];
    }
    for (size_t i=0; i<numCol; ++i) {
        rowForCol[i] = row4col[i];
    }
    free(tempBuffer);
    free(u);
    free(v);
    return(feasible);
}

// int main() {
//     srand(time(NULL));   // Initialization, should only be called once.
//     double *costMatrix;
    
//     int * col4row;
//     int * row4col;
//     // double *  u;
//     // double *  v;
//     // void *tempBuffer;
    
//     double gain;

//     size_t numRow = 3;
//     size_t numCol = 3;
//     costMatrix = (double *) malloc(numRow * numCol * sizeof(double));
//     col4row = (int *) malloc(numCol * sizeof(int));
//     row4col = (int *) malloc(numRow * sizeof(int));
//     // Assign costs to the costMatrix
//     for(int i=0; i<numRow; ++i){
//         for(int j=0; j<numCol; ++j){
//             costMatrix[i+numRow*j]  =  ((double)rand() / RAND_MAX);      // Returns a pseudo-random integer between 0 and RAND_MAX.
//         }
//     }
//     //  double assign2DCBasic(const double *C, ptrdiff_t * col4row, ptrdiff_t *  row4col, void *tempBuffer, double *  u, double *  v, const size_t numRow, const size_t numCol)       
//     // double totalCost = assign2DCBasic(costMatrix, col4row, row4col, bufferStart, u, v, numRow, numCol);

//     // bool assign2DC(const bool maximize, double *  C, double *  gain, ptrdiff_t *  col4row, ptrdiff_t *  row4col, void *tempBuffer, double *  u, double *  v, const size_t numRow, const size_t numCol)

//     bool feasible = bipartite_assignment(false, costMatrix, &gain, col4row, row4col, numRow, numCol);
    
//     printf("Total cost: %f.\n", gain);

//     for(int i=0; i<numRow; ++i){
//         for(int j=0; j<numCol; ++j){
//             printf("%f2.3 ", costMatrix[i+numRow*j]);      // Returns a pseudo-random integer between 0 and RAND_MAX.
//         }
//         printf("\n");
//     }
//     printf("\n");

//     for (int i=0; i<numRow; ++i) {
//         printf("%d -> %d (%f2.3)\n", i, col4row[i], costMatrix[i+numRow*col4row[i]]);
//     }
//     free(costMatrix);
//     free(col4row);
//     free(row4col);
//     // free(tempBuffer);
//     // free(u);
//     // free(v);
//     return 0;
// }
