#include "downstream_consumer.h"

#include <stdexcept>

#include "slu_mt_ddefs.h"
#include "slu_mt_dddefs.h"

namespace downstreamConsumerProg {

namespace {

double solve_with_double()
{
    auto *colptr = static_cast<int_t *>(SUPERLU_MALLOC(3 * sizeof(int_t)));
    auto *rowind = static_cast<int_t *>(SUPERLU_MALLOC(2 * sizeof(int_t)));
    auto *avals = static_cast<double *>(SUPERLU_MALLOC(2 * sizeof(double)));
    auto *rhsvals = static_cast<double *>(SUPERLU_MALLOC(2 * sizeof(double)));
    if (!colptr || !rowind || !avals || !rhsvals) {
        throw std::runtime_error("allocation failed");
    }
    colptr[0] = 0;
    colptr[1] = 1;
    colptr[2] = 2;
    rowind[0] = 0;
    rowind[1] = 1;
    avals[0] = 2.0;
    avals[1] = 4.0;
    rhsvals[0] = 6.0;
    rhsvals[1] = 20.0;
    int_t perm_c[] = {0, 1};
    int_t perm_r[] = {0, 0};
    int_t info = 0;

    SuperMatrix A, B, L, U;
    dCreate_CompCol_Matrix(&A, 2, 2, 2, avals, rowind, colptr, SLU_NC, SLU_D, SLU_GE);
    dCreate_Dense_Matrix(&B, 2, 1, rhsvals, 2, SLU_DN, SLU_D, SLU_GE);

    pdgssv(1, &A, perm_c, perm_r, &L, &U, &B, &info);
    if (info != 0) {
        throw std::runtime_error("pdgssv failed");
    }

    auto *x = static_cast<double *>(((DNformat *) B.Store)->nzval);
    const double result = x[0] + x[1];

    Destroy_CompCol_Matrix(&A);
    Destroy_SuperMatrix_Store(&B);
    SUPERLU_FREE(rhsvals);
    Destroy_SuperNode_SCP(&L);
    Destroy_CompCol_NCP(&U);
    return result;
}

double solve_with_doubledouble()
{
    auto *colptr = static_cast<int_t *>(SUPERLU_MALLOC(3 * sizeof(int_t)));
    auto *rowind = static_cast<int_t *>(SUPERLU_MALLOC(2 * sizeof(int_t)));
    auto *avals = static_cast<superlu_dd_real *>(SUPERLU_MALLOC(2 * sizeof(superlu_dd_real)));
    auto *rhsvals = static_cast<superlu_dd_real *>(SUPERLU_MALLOC(2 * sizeof(superlu_dd_real)));
    if (!colptr || !rowind || !avals || !rhsvals) {
        throw std::runtime_error("allocation failed");
    }
    colptr[0] = 0;
    colptr[1] = 1;
    colptr[2] = 2;
    rowind[0] = 0;
    rowind[1] = 1;
    avals[0] = superlu_dd_real(2.0);
    avals[1] = superlu_dd_real(4.0);
    rhsvals[0] = superlu_dd_real(6.0);
    rhsvals[1] = superlu_dd_real(20.0);
    int_t perm_c[] = {0, 1};
    int_t perm_r[] = {0, 0};
    int_t info = 0;

    SuperMatrix A, B, L, U;
    ddCreate_CompCol_Matrix(&A, 2, 2, 2, avals, rowind, colptr, SLU_NC, SLU_D, SLU_GE);
    ddCreate_Dense_Matrix(&B, 2, 1, rhsvals, 2, SLU_DN, SLU_D, SLU_GE);

    pddgssv(1, &A, perm_c, perm_r, &L, &U, &B, &info);
    if (info != 0) {
        throw std::runtime_error("pddgssv failed");
    }

    auto *x = static_cast<superlu_dd_real *>(((DNformat *) B.Store)->nzval);
    const double result = x[0].upper + x[1].upper;

    Destroy_CompCol_Matrix(&A);
    Destroy_SuperMatrix_Store(&B);
    SUPERLU_FREE(rhsvals);
    Destroy_SuperNode_SCP(&L);
    Destroy_CompCol_NCP(&U);
    return result;
}

}  // namespace

double solve_with_runtime_dispatch(Precision precision)
{
    return precision == Precision::Double ? solve_with_double() : solve_with_doubledouble();
}

}  // namespace downstreamConsumerProg
