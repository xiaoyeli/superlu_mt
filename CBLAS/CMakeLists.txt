add_library(blas STATIC
    # single precision real
    isamax.c sasum.c saxpy.c scopy.c sdot.c snrm2.c srot.c sscal.c
    sgemv.c ssymv.c strsv.c sger.c ssyr2.c

    # double precision real
    idamax.c dasum.c daxpy.c dcopy.c ddot.c dnrm2.c drot.c dscal.c
    dgemv.c dsymv.c dtrsv.c dger.c dsyr2.c

    # single precision complex
    icamax.c scasum.c caxpy.c ccopy.c scnrm2.c cscal.c
    cgemv.c chemv.c ctrsv.c cgerc.c cher2.c

    # double precision complex
    izamax.c dzasum.c zaxpy.c zcopy.c dznrm2.c zscal.c dcabs1.c
    zgemv.c zhemv.c ztrsv.c zgerc.c zher2.c
)

set_target_properties(
    blas
    PROPERTIES
        POSITION_INDEPENDENT_CODE ON
        C_STANDARD 99
)

target_include_directories(
    blas
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}
        ${PROJECT_SOURCE_DIR}/SRC
)

install(
    TARGETS blas
    DESTINATION ${CMAKE_INSTALL_LIBDIR}
)
