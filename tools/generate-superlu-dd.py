#!/usr/bin/env python3

from __future__ import annotations

import argparse
import re
from pathlib import Path


SOURCES = [
    "dgstrs.c",
    "dmyblas2.c",
    "dsp_blas2.c",
    "pdmemory.c",
    "pdgssv.c",
    "pdgstrf.c",
    "pdgstrf_bmod1D.c",
    "pdgstrf_bmod1D_mv2.c",
    "pdgstrf_bmod2D.c",
    "pdgstrf_bmod2D_mv2.c",
    "pdgstrf_column_bmod.c",
    "pdgstrf_column_dfs.c",
    "pdgstrf_copy_to_ucol.c",
    "pdgstrf_factor_snode.c",
    "pdgstrf_init.c",
    "pdgstrf_panel_bmod.c",
    "pdgstrf_panel_dfs.c",
    "pdgstrf_pivotL.c",
    "pdgstrf_snode_bmod.c",
    "pdgstrf_snode_dfs.c",
    "pdgstrf_thread.c",
    "pdgstrf_thread_finalize.c",
    "pdgstrf_thread_init.c",
]

WHOLE_DOUBLE_REWRITE = {
    "dgstrs.c",
    "dmyblas2.c",
    "dsp_blas2.c",
    "pdmemory.c",
    "pdgstrf_bmod1D.c",
    "pdgstrf_bmod1D_mv2.c",
    "pdgstrf_bmod2D.c",
    "pdgstrf_bmod2D_mv2.c",
    "pdgstrf_column_bmod.c",
    "pdgstrf_column_dfs.c",
    "pdgstrf_copy_to_ucol.c",
    "pdgstrf_factor_snode.c",
    "pdgstrf_panel_bmod.c",
    "pdgstrf_panel_dfs.c",
    "pdgstrf_pivotL.c",
    "pdgstrf_snode_bmod.c",
    "pdgstrf_snode_dfs.c",
    "pdgstrf_thread.c",
}

REPLACEMENTS = [
    (r'\bpdgstrf_', 'pddgstrf_'),
    (r'\bpdgstrf\b', 'pddgstrf'),
    (r'\bpdgssv\b', 'pddgssv'),
    (r'\bdgstrs\b', 'ddgstrs'),
    (r'\bdprint_soln\b', 'ddprint_soln'),
    (r'\bsp_dtrsv\b', 'sp_ddtrsv'),
    (r'\bdlsolve\b', 'ddlsolve'),
    (r'\bdusolve\b', 'ddusolve'),
    (r'\bdmatvec2\b', 'ddmatvec2'),
    (r'\bdmatvec\b', 'ddmatvec'),
    (r'\bdCreate_', 'ddCreate_'),
    (r'\bdoubleMalloc\b', 'ddRealMalloc'),
    (r'\bdoubleCalloc\b', 'ddRealCalloc'),
    (r'\bdallocateA\b', 'ddallocateA'),
    (r'\bdPresetMap\b', 'ddPresetMap'),
    (r'\bsuperlu_dQuerySpace\b', 'superlu_ddQuerySpace'),
    (r'\bcopy_mem_double\b', 'copy_mem_dd'),
    (r'\bdexpanders\b', 'ddexpanders'),
    (r'\bdfill\b', 'ddfill'),
]


def generate_internal_header(src_dir: Path, out_dir: Path) -> None:
    text = (src_dir / "slu_mt_ddefs.h").read_text()
    text = text.replace("__SLU_MT_DDEFS", "__SLU_MT_DD_INTERNAL")
    text = text.replace("slu_mt_ddefs.h", "slu_mt_dd_internal.h")
    insert_after = '#include "pxgstrf_synch.h"\n'
    snippet = """
#ifdef __cplusplus
#include <cmath>
#include "doubledouble.h"
using superlu_dd_real = doubledouble::DoubleDouble;
inline double superlu_dd_abs(const superlu_dd_real &x) { return x.abs().upper; }
inline double superlu_dd_abs(double x) { return std::fabs(x); }
#else
#error "Double-double support requires C++ compilation"
#endif

#ifdef USE_VENDOR_BLAS
#undef USE_VENDOR_BLAS
#endif

#define fabs superlu_dd_abs
"""
    text = text.replace(insert_after, insert_after + snippet, 1)
    for pattern, replacement in REPLACEMENTS:
        text = re.sub(pattern, replacement, text)
    text = re.sub(r'\bdouble\s*\*', 'superlu_dd_real *', text)
    text = text.replace("    double  *lusup;   /* L supernodes */",
                        "    superlu_dd_real *lusup;   /* L supernodes */")
    text = text.replace("    double  *ucol;    /* U columns */",
                        "    superlu_dd_real *ucol;    /* U columns */")
    text = text.replace("extern void    ddfill(superlu_dd_real *, int_t, double);",
                        "extern void    ddfill(superlu_dd_real *, int_t, superlu_dd_real);")
    (out_dir / "slu_mt_dd_internal.h").write_text(text)


def patch_common(text: str) -> str:
    text = text.replace('#include "slu_mt_ddefs.h"', '#include "slu_mt_dd_internal.h"')
    for pattern, replacement in REPLACEMENTS:
        text = re.sub(pattern, replacement, text)
    return text


def patch_double_rewrite(filename: str, text: str) -> str:
    if filename in WHOLE_DOUBLE_REWRITE:
        text = re.sub(r'\bdouble\b', 'superlu_dd_real', text)

    if filename == "pdgstrf_bmod1D.c":
        text = text.replace("superlu_dd_real *utime = Gstat->utime;", "double *utime = Gstat->utime;")
        text = text.replace("superlu_dd_real f_time;", "double f_time;")
    elif filename == "pdgstrf_bmod1D_mv2.c":
        text = text.replace("superlu_dd_real *utime = Gstat->utime;", "double *utime = Gstat->utime;")
        text = text.replace("superlu_dd_real f_time;", "double f_time;")
    elif filename == "pdgstrf_bmod2D.c":
        text = text.replace("superlu_dd_real *utime = Gstat->utime;", "double *utime = Gstat->utime;")
        text = text.replace("superlu_dd_real f_time;", "double f_time;")
    elif filename == "pdgstrf_bmod2D_mv2.c":
        text = text.replace("superlu_dd_real *utime = Gstat->utime;", "double *utime = Gstat->utime;")
        text = text.replace("superlu_dd_real f_time;", "double f_time;")
    elif filename == "pdgstrf_panel_bmod.c":
        text = text.replace("superlu_dd_real   t1, t2; /* temporary time */",
                            "double   t1, t2; /* temporary time */")
    elif filename == "pdgstrf_thread.c":
        text = text.replace("superlu_dd_real *utime = Gstat->utime;", "double *utime = Gstat->utime;")
        text = text.replace("superlu_dd_real t1, t2, t, stime;", "double t1, t2, t, stime;")
        text = text.replace("superlu_dd_real     diag_pivot_thresh = superlumt_options->diag_pivot_thresh;",
                            "double     diag_pivot_thresh = superlumt_options->diag_pivot_thresh;")
    elif filename == "pdgstrf.c":
        text = text.replace('#include "slu_mt_dd_internal.h"\n',
                            '#include "slu_mt_dd_internal.h"\nextern "C" double usertimer_();\n',
                            1)
        text = text.replace("    double    usertimer_();\n", "")
        text = text.replace("    superlu_dd_real    usertimer_();\n", "")
    elif filename == "pdgstrf_pivotL.c":
        text = text.replace("const superlu_dd_real u,", "const double u,")
        text = text.replace("register superlu_dd_real pivmax, rtemp, thresh;",
                            "register double pivmax, rtemp, thresh;")
    elif filename == "pdgstrf_factor_snode.c":
        text = text.replace("const superlu_dd_real diag_pivot_thresh,",
                            "const double diag_pivot_thresh,")
    elif filename == "dgstrs.c":
        text = text.replace('printf("\\t" IFMT ": %.10f\\n", i, soln[i]);',
                            'printf("\\t" IFMT ": %.10f\\n", i, static_cast<double>(soln[i].upper));')
    return text


def patch_cpp_casts(filename: str, text: str) -> str:
    file_specific = {
        "dgstrs.c": [
            ("Bstore = B->Store;", "Bstore = (DNformat *) B->Store;"),
            ("Bmat = Bstore->nzval;", "Bmat = (superlu_dd_real *) Bstore->nzval;"),
            ("Lstore = L->Store;", "Lstore = (SCPformat *) L->Store;"),
            ("Lval = Lstore->nzval;", "Lval = (superlu_dd_real *) Lstore->nzval;"),
            ("Ustore = U->Store;", "Ustore = (NCPformat *) U->Store;"),
            ("Uval = Ustore->nzval;", "Uval = (superlu_dd_real *) Ustore->nzval;"),
        ],
        "dsp_blas2.c": [
            ("Lstore = (SCPformat*) L->Store;", "Lstore = (SCPformat*) L->Store;"),
            ("Ustore = (NCPformat*) U->Store;", "Ustore = (NCPformat*) U->Store;"),
            ("Lval = (superlu_dd_real*) Lstore->nzval;", "Lval = (superlu_dd_real*) Lstore->nzval;"),
            ("Uval = (superlu_dd_real*) Ustore->nzval;", "Uval = (superlu_dd_real*) Ustore->nzval;"),
            ("Astore = A->Store;", "Astore = (NCformat *) A->Store;"),
            ("Aval = Astore->nzval;", "Aval = (superlu_dd_real *) Astore->nzval;"),
        ],
        "pdgssv.c": [
            ("Astore = A->Store;", "Astore = (NCformat *) A->Store;"),
            ("Bstore = B->Store;", "Bstore = (DNformat *) B->Store;"),
            ("NRformat *Astore = (NCformat *) A->Store;", "NRformat *Astore = (NRformat *) A->Store;"),
        ],
        "pdgstrf_thread_init.c": [
            ("Astore     = A->Store;", "Astore     = (NCPformat *) A->Store;"),
        ],
        "pdgstrf_panel_dfs.c": [
            ("Astore     = A->Store;", "Astore     = (NCPformat *) A->Store;"),
            ("a          = Astore->nzval;", "a          = (superlu_dd_real *) Astore->nzval;"),
        ],
        "pdgstrf_factor_snode.c": [
            ("Astore   = A->Store;", "Astore   = (NCPformat *) A->Store;"),
            ("a        = Astore->nzval;", "a        = (superlu_dd_real *) Astore->nzval;"),
        ],
        "pdmemory.c": [
            ("extern void    copy_mem_int    (int_t, void *, void *);",
             'extern "C" void copy_mem_int(int_t, void *, void *);'),
            ("extern void    user_bcopy      (char *, char *, int_t);",
             'extern "C" void user_bcopy(char *, char *, int_t);'),
            ("Lstore = L->Store;", "Lstore = (SCPformat *) L->Store;"),
            ("Ustore = U->Store;", "Ustore = (NCPformat *) U->Store;"),
            ("Lstore   = L->Store;", "Lstore   = (SCPformat *) L->Store;"),
            ("Ustore   = U->Store;", "Ustore   = (NCPformat *) U->Store;"),
            ("lusup = ddexpanders[LUSUP].mem = Lstore->nzval;",
             "lusup = ddexpanders[LUSUP].mem = (superlu_dd_real *) Lstore->nzval;"),
            ("ucol  = ddexpanders[UCOL].mem  = Ustore->nzval;;",
             "ucol  = ddexpanders[UCOL].mem  = (superlu_dd_real *) Ustore->nzval;;"),
            ("Astore   = A->Store;", "Astore   = (NCPformat *) A->Store;"),
            ("void *new)", "void *new_mem)"),
            ("superlu_dd_real *dnew = new;", "superlu_dd_real *dnew = (superlu_dd_real *) new_mem;"),
            ("superlu_dd_real   alpha = EXPAND;", "double   alpha = EXPAND;"),
        ],
    }
    for old, new in file_specific.get(filename, []):
        text = text.replace(old, new)
    if filename == "pdmemory.c":
        text = re.sub(r'copy_mem_dd\(int_t howmany, void \*old, void \*new\)',
                      'copy_mem_dd(int_t howmany, void *old, void *new_mem)', text)
        text = re.sub(r'superlu_dd_real\s+\*dnew = new;',
                      'superlu_dd_real *dnew = (superlu_dd_real *) new_mem;', text)
        text = re.sub(r'superlu_dd_real\s+alpha = EXPAND;',
                      'double alpha = EXPAND;', text)
    return text


def generate_source(src_dir: Path, out_dir: Path, filename: str) -> None:
    text = (src_dir / filename).read_text()
    text = patch_common(text)
    text = patch_double_rewrite(filename, text)
    text = patch_cpp_casts(filename, text)
    out_name = filename[:-2] + "_dd.cpp"
    (out_dir / out_name).write_text(text)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-dir", required=True)
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()

    src_dir = Path(args.source_dir)
    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    generate_internal_header(src_dir, out_dir)
    for filename in SOURCES:
        generate_source(src_dir, out_dir, filename)


if __name__ == "__main__":
    main()
