#include "transport_LLE2.h"
#include "gsl_optimizations2.h"
#include "MACROS.H"

/* code adopted from examples of GSL
 * http://www.gnu.org/software/gsl/manual/html_node/Multimin-Examples.html#Multimin-Examples
 */

#define     _LLE_INPUT_ \
    optim->n,  \
    optim->Pc, \
    optim->Tc, \
    optim->Vc, \
    optim->w,  \
    optim->MW, \    
    optim->Tb, \
    optim->SG, \
    optim->H8, \
    optim->kij, \
    \
    optim->n_flux_type,   \
    optim->dL, optim->dR, \
    optim->P,   \
    optim->T_s, \
    \
    optim->x_L,    \
    optim->x_R,    \
    optim->x_1,    \
    optim->x_2,    \
    optim->lnphi_1, optim->lnphi_2,\
    optim->Dij_1,   optim->Dij_2,  \
    optim->H_1,     optim->H_2,    \
    optim->flux_m_1,\
    optim->flux_m_2,\
    optim->flux_umf

#define     _ASSIGN_PARAMETERS_ \
    optim.n  = n;  \
    optim.Pc = Pc; \
    optim.Tc = Tc; \
    optim.Vc = Vc; \
    optim.w  = w;  \
    optim.MW = MW; \
    optim.Tb = Tb; \
    optim.SG = SG; \
    optim.H8 = H8; \
    optim.kij = kij; \
    \
    optim.n_flux_type = n_flux_type;   \
    optim.dL = dL; optim.dR = dR; \
    optim.P  = P;   \
    optim.T_s= T_s; \
    \
    optim.x_1 = x_1;     optim.x_L = x_L;    \
    optim.x_2 = x_2;     optim.x_R = x_R;    \
    optim.lnphi_1 = lnphi_1; optim.lnphi_2 = lnphi_2;\
    optim.Dij_1   = Dij_1;   optim.Dij_2   = Dij_2;  \
    optim.H_1     = H_1;     optim.H_2     = H_2;    \
    optim.flux_m_1 = flux_m_1;\
    optim.flux_m_2 = flux_m_2;\
    optim.flux_umf = flux_umf;

double transport_LLE2_f_(const gsl_vector *v, void *params)
{
  int j, n;

  LPT_GSL_OPTIM optim;
  double *x_1, *x_2;

  if (v==NULL || params==NULL) return GSL_NAN;

  optim = (LPT_GSL_OPTIM) params;
  x_1 = optim->x_1;
  x_2 = optim->x_2;
  n = optim->n;

  for (j=2; j<n; j++) {
    double _val = gsl_vector_get(v, j-2);
    if (j!=n-1) {
      x_1[j] = _val;
    }else{
      x_2[j] = _val;
    }
  }
  
  return transport_LLE2_eval_func(_LLE_INPUT_);
}

int my2_gsl_find_transport_LLE2(_PARAMETERS_)
{
  int j;
  double  tol;

  const gsl_multimin_fminimizer_type * gsl_T = 
        gsl_multimin_fminimizer_nmsimplex2;
  gsl_multimin_fminimizer *gsl_s = NULL;
  gsl_vector *gsl_ss, *gsl_x;
  gsl_multimin_function minex_func;

  T_GSL_OPTIM optim;
  void * gsl_par = &optim;

  size_t gsl_iter = 0;
  int gsl_status;
  double gsl_size;

  {_ASSIGN_PARAMETERS_}

  /* Starting point */
  gsl_x = gsl_vector_alloc (n-2);
  gsl_ss = gsl_vector_alloc (n-2);
  gsl_vector_set_all (gsl_ss, 1e-8);
  for (j=2;j<n;j++) {
    if (j!=n-1) {
      gsl_vector_set (gsl_x, j-2, x_1[j]);

    }else{
      gsl_vector_set (gsl_x, j-2, x_2[j]);

    }
  }

  /* Initialize method and iterate */
  minex_func.n = n-2;
  minex_func.f = transport_LLE2_f_;
  minex_func.params = gsl_par;

  gsl_s = gsl_multimin_fminimizer_alloc (gsl_T, n-2);
  gsl_multimin_fminimizer_set (gsl_s, &minex_func, gsl_x, gsl_ss);

  gsl_set_error_handler_off ();

  do
  {
    gsl_iter++;
    gsl_status = gsl_multimin_fminimizer_iterate(gsl_s);

    if (gsl_status) break;

    gsl_size = gsl_multimin_fminimizer_size (gsl_s);
    gsl_status = gsl_multimin_test_size (gsl_size, 1e-12);

    /*
    if (gsl_status == GSL_SUCCESS)
    {
      printf ("converged to minimum at\n");
    }

    printf ("%5d ", (int)gsl_iter);
    for (j=0;j<n-2;j++) printf("%10.6le ", gsl_vector_get (gsl_s->x, j)); 
    printf ("%le %le\n", gsl_s->fval, gsl_size);
    */
  }
  while (gsl_status == GSL_CONTINUE && gsl_iter < 1000);

  //transport_LLE2_eval_func(_LLE_INPUT_);

  gsl_status = gsl_iter;

  tol = 1e-3;
  if (gsl_s->fval>tol) gsl_status = -1;
  else if (fabs(x_1[n-1]-x_2[n-1])<2e-3) {
    if (fabs(x_L[n-1]-x_R[n-1])<1e-2) {
      gsl_status = -2;
      printf("x-: ");for (j=0;j<n;j++) printf("%9.10le ", x_L[j]);printf("\n");
      printf("x1: ");for (j=0;j<n;j++) printf("%9.10le ", x_1[j]);printf("\n");
      printf("x2: ");for (j=0;j<n;j++) printf("%9.10le ", x_2[j]);printf("\n");
      printf("x+: ");for (j=0;j<n;j++) printf("%9.10le ", x_R[j]);printf("\n");
    }else{
      gsl_status = -1;
    }
  }

  for (j=0;j<n;j++) {
      if (x_1[j]<0 || x_2[j]<0 
       || x_1[j]>1 || x_2[j]>1) 
      {
          gsl_status = -1;
          printf("X1: ");for (j=0;j<n;j++) printf("%9.10le ", x_1[j]);printf("\n");
          printf("X2: ");for (j=0;j<n;j++) printf("%9.10le ", x_2[j]);printf("\n");
      }
  }

  //MY_PRINTF("\tGSL2 returns %d, f0=%le\n", gsl_status, gsl_s->fval);

  if (gsl_status==-1) {
    _COPY_VECTOR_(x_L, x_1, n);
    _COPY_VECTOR_(x_R, x_2, n);
  }

  gsl_vector_free(gsl_x);
  gsl_vector_free(gsl_ss);
  gsl_multimin_fminimizer_free (gsl_s);

  return gsl_status;
}

double transport_LLE_f_(const gsl_vector *v, void *params)
{
  int j, k, id, n;

  LPT_GSL_OPTIM optim;
  double *x_1;

  if (v==NULL || params==NULL) return GSL_NAN;

  optim = (LPT_GSL_OPTIM) params;
  x_1 = optim->x_1;
  n = optim->n;

  // set k = 0 by default
  // or k can be calculated as the component with the maximum mole fraction in x_L[0<=i<=n-2]
  k = 0;

  id = 0;
  for (j=0; j<n-1; j++) {
    if (j!=k) {
      x_1[j] = gsl_vector_get(v, id);
      id++;
    }
  }
  
  return transport_LLE_eval_func(k, _LLE_INPUT_);
}

int my2_gsl_find_transport_LLE(_PARAMETERS_)
{
  int j, k, id;
  double  tol;

  const gsl_multimin_fminimizer_type * gsl_T = 
        gsl_multimin_fminimizer_nmsimplex2;
  gsl_multimin_fminimizer *gsl_s = NULL;
  gsl_vector *gsl_ss, *gsl_x;
  gsl_multimin_function minex_func;

  size_t gsl_iter = 0;
  int gsl_status;
  double gsl_size;

  T_GSL_OPTIM optim;
  void * gsl_par = &optim;  

  {_ASSIGN_PARAMETERS_}

  // set k = 0 by default
  // or k can be calculated as the component with the maximum mole fraction in x_1[0<=i<=n-2]
  k = 0;

  /* Starting point */
  gsl_x = gsl_vector_alloc (n-2);
  gsl_ss = gsl_vector_alloc (n-2);
  gsl_vector_set_all (gsl_ss, 1e-8);

  id = 0;
  for (j=0;j<n-1;j++) {
    if (j!=k) {
      gsl_vector_set (gsl_x, id, x_1[j]);

      id++;
    }
  }

  /* Initialize method and iterate */
  minex_func.n = n-2;
  minex_func.f = transport_LLE_f_;
  minex_func.params = gsl_par;

  gsl_s = gsl_multimin_fminimizer_alloc (gsl_T, n-2);
  gsl_multimin_fminimizer_set (gsl_s, &minex_func, gsl_x, gsl_ss);

  gsl_set_error_handler_off ();

  do
  {
    gsl_iter++;
    gsl_status = gsl_multimin_fminimizer_iterate(gsl_s);

    if (gsl_status) break;

    gsl_size = gsl_multimin_fminimizer_size (gsl_s);
    gsl_status = gsl_multimin_test_size (gsl_size, 1e-9);

    /*
    if (gsl_status == GSL_SUCCESS)
    {
      printf ("converged to minimum at\n");
    }

    printf ("%5d ", (int)gsl_iter);
    for (j=0;j<n-2;j++) printf("%10.6le ", gsl_vector_get (gsl_s->x, j)); 
    printf ("%le %le\n", gsl_s->fval, gsl_size);
    */
  }
  while (gsl_status == GSL_CONTINUE && gsl_iter < 1000);

  gsl_status = gsl_iter;

  tol = 1e-3;
  if (gsl_s->fval>tol) gsl_status = -1;
  else if (fabs(x_1[n-1]-x_2[n-1])<2e-3) {
    if (fabs(x_L[n-1]-x_R[n-1])<1e-2) {
      gsl_status = -2;
      printf("x1: ");for (j=0;j<n;j++) printf("%9.10le ", x_1[j]);printf("\n");
      printf("x2: ");for (j=0;j<n;j++) printf("%9.10le ", x_2[j]);printf("\n");
    }else{
      gsl_status = -1;
    }
  }

  for (j=0;j<n;j++) {
    if (x_1[j]<0 || x_2[j]<0) {
      gsl_status = -1;
    }
  }

  if (gsl_status<0) {
    printf("\tGSL steps %d\n", gsl_iter);
    printf("x1: ");for (j=0;j<n;j++) printf("%9.10le ", x_1[j]);printf("\n");
    printf("x2: ");for (j=0;j<n;j++) printf("%9.10le ", x_2[j]);printf("\n");
  }
  //MY_PRINTF("\tGSL returns %d, f0=%le\n", gsl_status, gsl_s->fval);

  if (gsl_status==-1) {
    _COPY_VECTOR_(x_L, x_1, n);
    _COPY_VECTOR_(x_R, x_2, n);
  }

  gsl_vector_free(gsl_x);
  gsl_vector_free(gsl_ss);
  gsl_multimin_fminimizer_free (gsl_s);

  return gsl_status;
}

