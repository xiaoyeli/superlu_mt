/*
 * -- SuperLU MT routine (version 1.0) --
 * Univ. of California Berkeley, Xerox Palo Alto Research Center,
 * and Lawrence Berkeley National Lab.
 * August 15, 1997
 *
 * Modified by Daniel C. Simkins, Jr.
 * University of South Florida
 * July, 2010
 */
#include <stdlib.h>
#include <stdio.h>
#include "pdsp_defs.h"
#include "slu_mt_util.h"
#include "slu_mt_Cnames.h"

typedef long long fptr;  /* 64-bit by default */

typedef struct {
    SuperMatrix *L;
    SuperMatrix *U;
    int *perm_c;
    int *perm_r;
} factors_t;

void
c_bridge_pdgssv_(int *nprocs, int *iopt, int *n, int *nnz, int *nrhs, double *values,
		 int *rowind, int *colptr, double *b, int *ldb, 
		 fptr *f_factors, /* a handle containing the address
				     pointing to the factored matrices */
		 int *info)
{
    SuperMatrix A, B, *L, *U;
    SCPformat *Lstore;
    NCPformat *Ustore;
    int      *perm_r; /* row permutations from partial pivoting */
    int      *perm_c; /* column permutation vector */
    int      panel_size, permc_spec, i, relax;
    superlu_memusage_t superlu_memusage;
    Gstat_t  Gstat;
    factors_t *LUfactors;
    trans_t  trans;

    trans = NOTRANS;

    if(*iopt == 1) /* factor and solve */
      {
	/* Adjust to 0-based indexing */
	for (i = 0; i < *nnz; ++i) --rowind[i];
	for (i = 0; i <= *n; ++i) --colptr[i];

	dCreate_CompCol_Matrix(&A, *n, *n, *nnz, values, rowind, colptr, 
			       SLU_NC, SLU_D, SLU_GE);
	dCreate_Dense_Matrix(&B, *n, *nrhs, b, *ldb, SLU_DN, SLU_D, SLU_GE);

	L = (SuperMatrix *) SUPERLU_MALLOC( sizeof(SuperMatrix) );
	U = (SuperMatrix *) SUPERLU_MALLOC( sizeof(SuperMatrix) );
	if ( !(perm_r = intMalloc(*n)) ) USER_ABORT("Malloc fails for perm_r[].");
	if ( !(perm_c = intMalloc(*n)) ) USER_ABORT("Malloc fails for perm_c[].");

	/*
	 * Get column permutation vector perm_c[], according to permc_spec:
	 *   permc_spec = 0: natural ordering 
	 *   permc_spec = 1: minimum degree ordering on structure of A'*A
	 *   permc_spec = 2: minimum degree ordering on structure of A'+A
	 *   permc_spec = 3: approximate minimum degree for unsymmetric matrices
	 */    	
	permc_spec = 1;
	get_perm_c(permc_spec, &A, perm_c);

	panel_size = sp_ienv(1);
    
	pdgssv(*nprocs, &A, perm_c, perm_r, L, U, &B, info);

	if ( *info == 0 ) {

	  Lstore =  L->Store;
	  Ustore =  U->Store;
	  printf("#NZ in factor L = %d\n", Lstore->nnz);
	  printf("#NZ in factor U = %d\n", Ustore->nnz);
	  printf("#NZ in L+U = %d\n", Lstore->nnz + Ustore->nnz - L->ncol);
	
	  superlu_dQuerySpace(*nprocs, L, U, panel_size, &superlu_memusage);
	  printf("L\\U MB %.3f\ttotal MB needed %.3f\texpansions %d\n",
		 superlu_memusage.for_lu/1e6, superlu_memusage.total_needed/1e6,
		 superlu_memusage.expansions);
	
	} else {
	  printf("dgssv() error returns INFO= %d\n", *info);
	  if ( info <= n ) { /* factorization completes */
	    superlu_dQuerySpace(*nprocs, L, U, panel_size, &superlu_memusage);
	    printf("L\\U MB %.3f\ttotal MB needed %.3f\texpansions %d\n",
		   superlu_memusage.for_lu/1e6, 
		   superlu_memusage.total_needed/1e6,
		   superlu_memusage.expansions);
	  }
	}
	/* Restore to 1-based indexing */
	for (i = 0; i < *nnz; ++i) ++rowind[i];
	for (i = 0; i <= *n; ++i) ++colptr[i];

	/* Save the LU factors in the factors handle */
	LUfactors = (factors_t*) SUPERLU_MALLOC(sizeof(factors_t));
	LUfactors->L = L;
	LUfactors->U = U;
	LUfactors->perm_c = perm_c;
	LUfactors->perm_r = perm_r;
	*f_factors = (fptr) LUfactors;
	
	Destroy_SuperMatrix_Store(&A);
	Destroy_SuperMatrix_Store(&B);
      }
    else if( *iopt == 2) /* do another solve */
      {
	/* Initialize the statistics variables. */
	/* Extract the LU factors in the factors handle */
	LUfactors = (factors_t*) *f_factors;
	L = (LUfactors->L);
	U = (LUfactors->U);
	perm_c = LUfactors->perm_c;
	perm_r = LUfactors->perm_r;

	panel_size         = sp_ienv(1);
	relax              = sp_ienv(2);
 
	/* 	n = A.ncol; */
	StatAlloc(*n, *nprocs, panel_size, relax, &Gstat);
	StatInit(*n, *nprocs, &Gstat);

	dCreate_Dense_Matrix(&B, *n, *nrhs, b, *ldb, SLU_DN, SLU_D, SLU_GE);

        /* Solve the system A*X=B, overwriting B with X. */

/*       NOTE: The order of the arguments perm_r and perm_c are
	 opposite for this call than for the call in the single
	 threaded version of SuperLU. */

        dgstrs (trans, L, U, perm_r, perm_c, &B, &Gstat, info);

	Destroy_SuperMatrix_Store(&B);

	PrintStat(&Gstat);
	StatFree(&Gstat);

      }
    else if( *iopt == 3)
      {
	/*     SUPERLU_FREE (perm_r); */
	/*     SUPERLU_FREE (perm_c); */
	/*     Destroy_SuperNode_Matrix(&L); */
	/*     Destroy_CompCol_Matrix(&U); */
	/* Free the LU factors in the factors handle */
	LUfactors = (factors_t*) *f_factors;

	/* 	printf("f_factors=%ld\n",*f_factors) ; */

	if(LUfactors)
	  {
	    if(LUfactors->perm_r)SUPERLU_FREE (LUfactors->perm_r);
	    SUPERLU_FREE (LUfactors->perm_c);
	    Destroy_SuperNode_Matrix(LUfactors->L);
	    Destroy_CompCol_Matrix(LUfactors->U);
	    SUPERLU_FREE (LUfactors->L);
	    SUPERLU_FREE (LUfactors->U);
	    SUPERLU_FREE (LUfactors);
	  }
      }
    else
      {
	fprintf(stderr,"Invalid iopt=%d passed to c_bridge_pdgssv()\n",*iopt);
	exit(-1);
      }
}

