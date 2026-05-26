/*! \file
Copyright (c) 2003, The Regents of the University of California, through
Lawrence Berkeley National Laboratory (subject to receipt of any required
approvals from U.S. Dept. of Energy)

All rights reserved.

The source code is distributed under BSD license, see the file License.txt
at the top-level directory.
*/

#ifndef __SLU_MT_DDDEFS
#define __SLU_MT_DDDEFS

#include "slu_mt_ddefs.h"

#ifdef __cplusplus
#include <cstddef>
#include <type_traits>
#include "doubledouble.h"
using superlu_dd_real = doubledouble::DoubleDouble;
static_assert(std::is_standard_layout<superlu_dd_real>::value,
              "superlu_dd_real must have a stable C-compatible field layout");
static_assert(sizeof(superlu_dd_real) == 2 * sizeof(double),
              "superlu_dd_real must contain exactly two doubles");
static_assert(offsetof(superlu_dd_real, upper) == 0,
              "superlu_dd_real upper field must be first");
static_assert(offsetof(superlu_dd_real, lower) == sizeof(double),
              "superlu_dd_real lower field must follow upper");
#else
typedef struct {
    double upper;
    double lower;
} superlu_dd_real;
#endif

#ifdef __cplusplus
extern "C" {
#endif

extern void
pddgssv(int_t, SuperMatrix *, int_t *, int_t *, SuperMatrix *, SuperMatrix *,
        SuperMatrix *, int_t *);

extern void
ddCreate_CompCol_Matrix(SuperMatrix *, int_t, int_t, int_t, superlu_dd_real *,
                        int_t *, int_t *, Stype_t, Dtype_t, Mtype_t);

extern void
ddCreate_CompCol_Permuted(SuperMatrix *, int_t, int_t, int_t, superlu_dd_real *,
                          int_t *, int_t *, int_t *, Stype_t, Dtype_t, Mtype_t);

extern void
ddCreate_Dense_Matrix(SuperMatrix *, int_t, int_t, superlu_dd_real *, int_t,
                      Stype_t, Dtype_t, Mtype_t);

#ifdef __cplusplus
}
#endif

#endif
