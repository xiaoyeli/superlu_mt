#include "slu_mt_dddefs.h"

extern "C" {

void ddCreate_CompCol_Matrix(SuperMatrix *A, int_t m, int_t n, int_t nnz,
                             superlu_dd_real *nzval, int_t *rowind,
                             int_t *colptr, Stype_t stype, Dtype_t dtype,
                             Mtype_t mtype)
{
    NCformat *Astore;

    A->Stype = stype;
    A->Dtype = dtype;
    A->Mtype = mtype;
    A->nrow = m;
    A->ncol = n;
    A->Store = (void *) SUPERLU_MALLOC(sizeof(NCformat));
    Astore = (NCformat *) A->Store;
    Astore->nnz = nnz;
    Astore->nzval = nzval;
    Astore->rowind = rowind;
    Astore->colptr = colptr;
}

void ddCreate_CompCol_Permuted(SuperMatrix *A, int_t m, int_t n, int_t nnz,
                               superlu_dd_real *nzval, int_t *rowind,
                               int_t *colbeg, int_t *colend, Stype_t stype,
                               Dtype_t dtype, Mtype_t mtype)
{
    NCPformat *Astore;

    A->Stype = stype;
    A->Dtype = dtype;
    A->Mtype = mtype;
    A->nrow = m;
    A->ncol = n;
    A->Store = (void *) SUPERLU_MALLOC(sizeof(NCPformat));
    Astore = (NCPformat *) A->Store;
    Astore->nnz = nnz;
    Astore->nzval = nzval;
    Astore->rowind = rowind;
    Astore->colbeg = colbeg;
    Astore->colend = colend;
}

void ddCreate_Dense_Matrix(SuperMatrix *X, int_t m, int_t n,
                           superlu_dd_real *x, int_t ldx, Stype_t stype,
                           Dtype_t dtype, Mtype_t mtype)
{
    DNformat *Xstore;

    X->Stype = stype;
    X->Dtype = dtype;
    X->Mtype = mtype;
    X->nrow = m;
    X->ncol = n;
    X->Store = (void *) SUPERLU_MALLOC(sizeof(DNformat));
    Xstore = (DNformat *) X->Store;
    Xstore->lda = ldx;
    Xstore->nzval = x;
}

void ddCreate_SuperNode_Permuted(SuperMatrix *L, int_t m, int_t n, int_t nnz,
                                 superlu_dd_real *nzval, int_t *nzval_colbeg,
                                 int_t *nzval_colend, int_t *rowind,
                                 int_t *rowind_colbeg, int_t *rowind_colend,
                                 int_t *col_to_sup, int_t *sup_to_colbeg,
                                 int_t *sup_to_colend, Stype_t stype,
                                 Dtype_t dtype, Mtype_t mtype)
{
    SCPformat *Lstore;

    L->Stype = stype;
    L->Dtype = dtype;
    L->Mtype = mtype;
    L->nrow = m;
    L->ncol = n;
    L->Store = (void *) SUPERLU_MALLOC(sizeof(SCPformat));
    Lstore = (SCPformat *) L->Store;
    Lstore->nnz = nnz;
    Lstore->nsuper = col_to_sup[n];
    Lstore->nzval = nzval;
    Lstore->nzval_colbeg = nzval_colbeg;
    Lstore->nzval_colend = nzval_colend;
    Lstore->rowind = rowind;
    Lstore->rowind_colbeg = rowind_colbeg;
    Lstore->rowind_colend = rowind_colend;
    Lstore->col_to_sup = col_to_sup;
    Lstore->sup_to_colbeg = sup_to_colbeg;
    Lstore->sup_to_colend = sup_to_colend;
}

void ddfill(superlu_dd_real *a, int_t alen, superlu_dd_real dval)
{
    int_t i;
    for (i = 0; i < alen; ++i) a[i] = dval;
}

}
