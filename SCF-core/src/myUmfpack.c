#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "MACROS.H"
#include "myUmfpack.h"


void	initialUmfpack(LPT_UMFPACK umf)
{
  constructUmfpack(umf);
  setControl(umf);
}

void	constructUmfpack(LPT_UMFPACK umf){
  umf->Arow=0;
  umf->Acol=0; 
  umf->AE=0;
  umf->rhs=0;
  umf->x=0;
  umf->Symbolic=0; 
  umf->Numeric=0;

  umf->rVal=0;
  umf->nColIdx=0; 
  umf->Acol_s=0; 
  umf->nRowPtr=0;
}

void	memoryMatVec(LPT_UMFPACK umf, int n)
{
  int	i;
  int	rows, cols, nz;

  //rows and cols
  rows = n*2 + 1;
  cols = rows;
  umf->rows = rows;  
  umf->cols = cols;

  nz = 0;
  for (i=0; i<n; i++)  // eq. marked as i=0,n-1
    nz += n*2;
  for (i=1; i<n; i++)  // eq. marked as i=1,n-1
    nz += 3;
  nz += n*2;           // two eq.: sum_xA_i=0, sum_xB_i=0
  umf->nz=nz;

  _ALLOCATE_MEM_;
  setValue0(umf);
}

void  setValue0(LPT_UMFPACK umf){
  int i;
  int nz=umf->nz;

  for ( i=0; i<nz; i++ ){
    umf->AE[i] = 1.0;
    umf->rVal[i] = 1.0;
  }
}

void	setControl(LPT_UMFPACK umf){
  umfpack_di_defaults (umf->Control) ;
  //umf->Control [UMFPACK_DENSE_ROW] = 0.4 ;
  //umf->Control [UMFPACK_DENSE_COL] = 0.4 ;
}

int	searchIndex(LPT_UMFPACK umf, int nr, int nc){
  int i, idx[2];
  for (i=0; i<2; i++)
    idx[i] = umf->Acol[nc+i];
  for (i=idx[0]; i<idx[1]; i++)
    if ( umf->Arow[i] == nr ){
      return i;
    }
  return -1;
}

void    compRow_setRowCol(LPT_UMFPACK umf, int n) 
{
  int i, j, iz, row;

  umf->nRowPtr[0] = 0;
  iz = 0;

  for (i=0; i<n; i++){ /* rows: 0, n-1 */
    umf->nRowPtr[i+1] = umf->nRowPtr[i] + n*2;
    for (j=0; j<n*2; j++){ /* cols: 0, 2*n-1 */
      umf->rVal[ iz ] = 1;
      umf->nColIdx[ iz++ ] = j;
    }                    
  }                     

  row = n;
  for (i=0; i<n-1; i++){ /* rows: row+0, row+n-1 */
    umf->nRowPtr[row+i+1] = umf->nRowPtr[row+i] + 3;
    umf->rVal[ iz ] = 1; umf->nColIdx[ iz++ ] =  i+1;         //   i+1;
    umf->rVal[ iz ] = 1; umf->nColIdx[ iz++ ] = (i+1)+n;      //  (i+1)+n;
    umf->rVal[ iz ] = 1; umf->nColIdx[ iz++ ] =  2*n;         //   2*n;
  }                     

  row += (n-1);
  umf->nRowPtr[row+1] = umf->nRowPtr[row] + n;
  for (j=0; j<n; j++){ /* cols: 0, n-1 */
    umf->rVal[ iz ] = 1;
    umf->nColIdx[ iz++ ] = j;
  }

  row++;
  umf->nRowPtr[row+1] = umf->nRowPtr[row] + n;
  for (j=0; j<n; j++){ /* cols: 0, n-1 */
    umf->rVal[ iz ] = 1;
    umf->nColIdx[ iz++ ] = j+n;
  }

  printf("\numfpack non-zero number is %d, & accumulative counting number is %d\n\n", umf->nz, iz);
  //umf->nz = iz;
  compCol_setRowCol(umf);
}

void    compCol_setRowCol(LPT_UMFPACK umf)
{
  int i, j, k, count;
#ifdef		_DEBUG_
  FILE *f;
#endif

  for (i=0; i<=umf->cols; i++)
    umf->Acol_s[i]=0;

  for (i=0; i<umf->nz; i++)
    umf->Acol_s[ umf->nColIdx[i] ] ++ ;

  umf->Acol[0] = 0 ;
  for (i=0; i<umf->cols; i++)
    umf->Acol[i+1] = umf->Acol[i] + umf->Acol_s[i];

  for (i=0; i<=umf->cols; i++)
    umf->Acol_s[i] = umf->Acol[i];

  count = 0;

  for (i=1; i<=umf->rows; i++)
    for (j=count; j<umf->nRowPtr[i]; j++){
      k=umf->Acol_s[umf->nColIdx[j]];

      umf->AE[k] = umf->rVal[j];
      umf->Arow[k] = i-1;

      umf->Acol_s[ umf->nColIdx[j] ]++;
      count++;
    }

#ifdef		_DEBUG_
  f = fopen("compCol.dat", "w");
  fprintf(f, "%d %d %d\n", 
      umf->rows,
      umf->cols,
      umf->nz);
  for (i=0; i<=umf->cols; i++){
    fprintf(f, "%d ", umf->Acol[i]);
  }
  fprintf(f, "\n");
  for (i=0; i<umf->nz; i++){
    fprintf(f, "%d ", umf->Arow[i]);
  }
  fprintf(f, "\n");
  for (i=0; i<umf->nz; i++){
    fprintf(f, "%le ", umf->AE[i]);
  }
  fprintf(f, "\n");
#endif
}

void    copyValRow_to_Col(LPT_UMFPACK umf)
{
  int i, j, count;

  for (i=0; i<=umf->cols; i++)
    umf->Acol_s[i] = umf->Acol[i];

  count = 0;

  for (i=1; i<=umf->rows; i++)
    for (j=count; j<umf->nRowPtr[i]; j++){
      umf->AE[ umf->Acol_s[
        umf->nColIdx[j] ] ] = umf->rVal[j];

      umf->Acol_s[ umf->nColIdx[j] ]++;
      count++;
    }
}

void	free_memory(LPT_UMFPACK umf){
  //---clear memory
  if(umf->x!=0){
    free(umf->x); 
    free(umf->rhs); 
    free(umf->AE); 
    free(umf->Acol); 
    free(umf->Arow);

    umf->Arow=0;
    umf->Acol=0; 
    umf->AE=0;
    umf->rhs=0;
    umf->x=0;
  }

  if (!umf->nColIdx) free(umf->nColIdx);
  if (!umf->Acol_s) free(umf->Acol_s);
  if (!umf->nRowPtr) free(umf->nRowPtr);
  if (!umf->rVal) free(umf->rVal);

  if(umf->Symbolic!=0) {
    umfpack_di_free_symbolic (&umf->Symbolic) ; 
    umf->Symbolic=0; 
  }
  if(umf->Numeric!=0) {
    umfpack_di_free_numeric (&umf->Numeric) ; 
    umf->Numeric=0; 
  }
}

//---decompose
void	decomposeUmfpack(LPT_UMFPACK umf){
  symbolicDecomposeUmfpack(umf);
  numericDecomposeUmfpack(umf);
}
void	symbolicDecomposeUmfpack(LPT_UMFPACK umf){
  int ierr;
  if(umf->Symbolic!=0) {
    umfpack_di_free_symbolic (&umf->Symbolic) ; 
  }
  ierr = umfpack_di_symbolic (umf->rows, umf->cols, 
      umf->Acol, umf->Arow, umf->AE, 
      &umf->Symbolic, umf->Control, umf->Info) ;
  if (ierr>0)
    printf("symbolic decompose return %d\n", ierr);
}
void	numericDecomposeUmfpack(LPT_UMFPACK umf){
  int ierr;
  if(umf->Numeric!=0) {
    umfpack_di_free_numeric (&umf->Numeric) ; 
  }
  ierr = umfpack_di_numeric (umf->Acol, umf->Arow, 
      umf->AE, umf->Symbolic, 
      &umf->Numeric, umf->Control, umf->Info) ;
  if (ierr>0)
    printf("numeric decompose return %d\n", ierr);
}

//---solve
void	solveUmfpack(LPT_UMFPACK umf){
  int ierr;
  /*
     double det, ex, mx;
     ierr = umfpack_di_get_determinant(&mx, &ex, 
     umf->Numeric, umf->Info);
     det = mx * pow (10.0, ex) ;
     printf("return value is %d det is %le ", ierr, det);
   */

  ierr = umfpack_di_solve (UMFPACK_A, umf->Acol, 
      umf->Arow, umf->AE, 
      umf->x, umf->rhs, 
      umf->Numeric, umf->Control, umf->Info);
  if (ierr>0)
  {
    printf("solver return %d\n", ierr);
    printf("condition number is %le \n", umf->Info[UMFPACK_RCOND]);
  }
}

//----------------------------------------------------------------------------
//---MAXWELL STEFAN SOLVER FOR MULTICOMPONENT MASS TRANSFER-------------------
//----------------------------------------------------------------------------
void	flux_memoryMatVec(LPT_UMFPACK umf, int n)
{
  int	i;
  int	rows, cols, nz;

  //rows and cols
  rows = n;
  cols = rows;
  umf->rows = rows;  
  umf->cols = cols;

  nz = 0;
  for (i=0; i<n; i++)  // eq. marked as i=0,n-1
    nz += n;
  umf->nz=nz;

  _ALLOCATE_MEM_;
  setValue0(umf);
}

void    flux_compRow_setRowCol(LPT_UMFPACK umf, int n) 
{
  int i, j, iz;

  umf->nRowPtr[0] = 0;
  iz = 0;

  for (i=0; i<n; i++){ /* rows: 0, n-1 */
    umf->nRowPtr[i+1] = umf->nRowPtr[i] + n;
    for (j=0; j<n; j++){ /* cols: 0, n-1 */
      umf->rVal[ iz ] = 1;
      umf->nColIdx[ iz++ ] = j;
    }                    
  }                     

  //printf("\nMaxwell-Stefan: umfpack non-zero number is %d, & accumulative counting number is %d\n\n", umf->nz, iz);
  //umf->nz = iz;
  compCol_setRowCol(umf);
}

void Maxwell_Stefan_mass_flux
(
    double Z, // compressibility factor
    int n, double *MW, double *x, double *Dij, double *rhs, 
    double *flux, LPT_UMFPACK umf
)
{
    int i, j, idx, iz;
    double small = 1e-8;
    //double sum, rhs_sum;
    double *y = (double*)malloc(sizeof(double)*n);
    int i0 = 0;

    //sum = 0;
    //for (i=0; i<n; i++) sum += rhs[i];
    //rhs_sum = sum;

    if (Z<1e-6) Z=1e-6;

    for (i=0; i<n; i++) 
    {
        //printf("i = %d\n", i);
        //printf("x[i] = %9.21e  small = %9.21e\n", x[i], small);
        if (i==i0) continue;

        double sum = 0.0;

        for (j=0; j<n; j++) 
        {
            //printf("j = %d\n", j);
            if (j==i) continue;

            idx = i + j*n;
            iz = searchIndex(umf, i, j);
            //printf("idx = %d  iz = %d\n", idx, iz);

            if (x[i]>small) 
            {
                //umf->AE[iz] = - x[i]/Dij[idx];
                //sum += x[j]/Dij[idx];
                umf->AE[iz] = - x[i]/(Dij[idx]*MW[j])*1e3;
                sum += x[j]/(Dij[idx]*MW[i])*1e3;
            }
            else
            {
                umf->AE[iz] = 0.0;
            }
            //printf("Dij[idx] = %9.21e  MW[j] = %9.21e\n", Dij[idx], MW[j]);
            //printf("AE[iz] = %9.21e  sum[i] = %9.21e\n", umf->AE[iz], sum);
        }

        iz = searchIndex(umf, i, i);
        umf->AE[iz] = sum;
        umf->rhs[i] = rhs[i];
        //umf->rhs[i] = rhs[i] - rhs_sum*y[i];
        //umf->rhs[i] = rhs[i] - rhs_sum*x[i];
    
        if (x[i]<=small) 
        {
            umf->AE[iz] = 1.0;
            umf->rhs[i] = 0.0;
        }
    }

    i=i0;
    for (j=0; j<n; j++) 
    {
        iz = searchIndex(umf, i, j);
        umf->AE[iz] = 1.0;
    }
    umf->rhs[i] = 0.0;
  
    //printf("My %9.2le\n", sqrt(MW[i]*y[i]));
    /*
    printf("\tMaxwell-Stefan matrix:\n");
    for (i=0; i<n; i++) 
    {
        printf("\t\t");
        for (j=0; j<n; j++) 
        {
            iz = searchIndex(umf, i, j);
            printf("%9.2le ", umf->AE[iz]);
        }
        printf("rhs %9.2le ", umf->rhs[i]);
        printf("x %9.2le\n", umf->x[i]);
    }
        */
    decomposeUmfpack(umf);
    solveUmfpack(umf);
    /*
    printf("\tMaxwell-Stefan matrix:\n");
    for (i=0; i<n; i++) 
    {
        printf("\t\t");
        for (j=0; j<n; j++) 
        {
            iz = searchIndex(umf, i, j);
            printf("%9.2le ", umf->AE[iz]);
        }
        printf("rhs %9.2le ", umf->rhs[i]);
        printf("x %9.2le\n", umf->x[i]);
    }
        */
    //printf("g_dx=%lf\n", g_T_s);
    /*
    printf("\tMaxwell-Stefan matrix:\n");
    for (i=0; i<n; i++) {
        printf("\t\t");
        for (j=0; j<n; j++) {
            iz = searchIndex(umf, i, j);
            printf("%9.2le ", umf->AE[iz]);
        }
        printf("rhs %9.2le ", umf->rhs[i]);
        printf("x %9.2le\n", umf->x[i]);
        //printf("My %9.2le\n", sqrt(MW[i]*y[i]));
    }
        */

    // convert molar flux into mass flux
    //sum = 0;
    //for (i=0; i<n; i++) sum += umf->x[i]*MW[i]*1e-3;

    // copy solutions out
    //for (i=0; i<n; i++) flux[i] = umf->x[i]*MW[i]*1e-3 - sum*y[i];
    for (i=0; i<n; i++) flux[i] = umf->x[i];

    free(y);

    /*
    printf("\tDij matrix:\n");
    for (i=0; i<n; i++) 
    {
        printf("\t\t");
        for (j=0; j<n; j++) 
        {
            int idx = i + j*n;
            printf("%9.2le ", Dij[idx]);
        }
        printf("\n");
    }
        */
}

