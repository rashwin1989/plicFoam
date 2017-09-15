#ifndef     _MACROS_H_
#define     _MACROS_H_

#include <stdlib.h>


// I. new and delete
// 1. one dimensional vextor var[size]
#define     _NEW_(var,type,size)        var = (type*) malloc(sizeof(type)*(size));
#define     _DELETE_(var)               free(var);

// 2. two dimensional matrix var[size1][size2]
#define     _NEW2_(var,type,size1, size2)  {\
         int _i; \
         _NEW_(var,type*,size1);\
         for (_i=0; _i<size1; _i++) _NEW_(var[_i],type,size2);\
}
#define     _DELETE2_(var,size1)  {\
         int _i; \
         for (_i=0; _i<size1; _i++) _DELETE_(var[_i]);\
         _DELETE_(var); \
}

// II. copy grid data
// 1. one dimensional vector var[size] to var_n[size]
#define     _COPY_VECTOR_(var,var_n,size) {\
         int _i; \
         for (_i=0; _i<size; _i++) var_n[_i]=var[_i];\
}
// 2. two dimensional matrix var[size1][size2] to var_n[size1][size2]
#define     _COPY_VECTOR2_(var,var_n,size1,size2) {\
         int _i, _j; \
         for (_i=0; _i<size1; _i++) \
           for (_j=0; _j<size2; _j++) \
             var_n[_i][_j]=var[_i][_j];\
}

#endif
