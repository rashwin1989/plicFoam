#ifndef  _MY_UMFPACK_H_
#define  _MY_UMFPACK_H_

/*   How to use UMFPACK

1) to customize two following functions to setup matrix shape
     void memoryMatVec(LPT_UMFPACK umf, int n);
     void compRow_setRowCol(LPT_UMFPACK umf, int n) 

2) to initialize data by calling
     initialUmfpack(&umf);
     memoryMatVec(&umf, n);
     compRow_setRowCol(&umf, n);

3) to edit matrix and RHS vector using
     iz = searchIndex(&umf, row_idx, col_idx);
     umf.AE[iz] = val;
     umf.rhs[row_idx] = f;

4) to solve linear equations
     solveUmfpack(&umf);
   and solutions are avaible in
     umf.x[row_idx]

5) finally, release the memory if UMFPACK is not needed
     free_memory(&umf);
 */


#include "umfpack.h"

typedef struct {
  int rows, cols;   //number of rows and columns 
  int nz;           //number of nonzero entries

  double *rhs, *x;    //righr hand side and solution vectors

  //---umfpack data structure
  int  *Arow, *Acol, *Acol_s;
  double *AE;

  //---compressed-row matrix
  int  *nRowPtr, *nColIdx;
  double *rVal;

  // data of information and control
  double Info[ UMFPACK_INFO ], 
      Control[ UMFPACK_CONTROL ];
  int status;
  void *Symbolic, *Numeric ;

}T_UMFPACK, *LPT_UMFPACK;

#define _ALLOCATE_MEM_  \
  umf->rhs=(double *)calloc(rows,sizeof(double)); \
  umf->x=(double *)calloc(rows,sizeof(double));   \
  umf->AE=(double *)calloc(nz,sizeof(double));    \
  umf->Acol=(int *)calloc(cols+1,sizeof(int));    \
  umf->Acol_s= (int *)calloc(cols+1,sizeof(int)); \
  umf->Arow=(int *)calloc(nz,sizeof(int));        \
  umf->rVal=(double *)calloc(nz,sizeof(double));  \
  umf->nRowPtr=(int *)calloc(rows+1,sizeof(int)); \
  umf->nColIdx=(int *)calloc(nz,sizeof(int));     

// Initialize
// ---------------------------------------------
void initialUmfpack(LPT_UMFPACK umf);
// ---------------------------------------------
// include:
void constructUmfpack(LPT_UMFPACK umf);
// set all value be 0 in the matrix
void setValue0(LPT_UMFPACK umf);
void setControl(LPT_UMFPACK umf);
// ---------------------------------------------
void memoryMatVec(LPT_UMFPACK umf, int n);
//
// Related to transform from 
// CompRow Matrix to CompCol Matrix
// ---------------------------------------------
// ---------------------------------------------
void compRow_setRowCol(LPT_UMFPACK umf, int n);
void compCol_setRowCol(LPT_UMFPACK umf);
void copyValRow_to_Col(LPT_UMFPACK umf);
int  searchIndex(LPT_UMFPACK umf, int nr, int nc);
// ---------------------------------------------
// ---------------------------------------------

// decompose
void decomposeUmfpack(LPT_UMFPACK umf);
// ---------------------------------------------
// include
void symbolicDecomposeUmfpack(LPT_UMFPACK umf);
void numericDecomposeUmfpack(LPT_UMFPACK umf);
// ---------------------------------------------

// Solve
// ---------------------------------------------
void solveUmfpack(LPT_UMFPACK umf);
// ---------------------------------------------

// Free
// ---------------------------------------------
void free_memory(LPT_UMFPACK umf);
// ---------------------------------------------

// A Maxwell-Stefan solver's Matrix structure
void flux_memoryMatVec(LPT_UMFPACK umf, int n);
void flux_compRow_setRowCol(LPT_UMFPACK umf, int n);

void Maxwell_Stefan_mass_flux(
       double Z, // compressibility factor
       int n, double *MW, double *x, double *Dij, double *rhs, 
       double *flux, LPT_UMFPACK umf);

#endif
