/** @file pmetis_fwrap.h
 * @brief Fortran77 interface of parmetis library
 *
 * @author Adam Peplinski
 * @date Mar 4, 2016
 *
 * @ingroup nekp4est
 */

#ifndef NEKP4EST_PMETIS_PMETIS_FWRAP_H_
#define NEKP4EST_PMETIS_PMETIS_FWRAP_H_

#if !defined(__parmetis_h__)
#error "pmetis_fwrap.h requires parmetis.h"
#endif

#if !defined(NAME_H)
#error "pmetis_fwrap.h requires name.h"
#endif

/* Function definitions */
#define MAX( a, b ) ( ( a > b) ? a : b )
#define MIN( a, b ) ( ( a > b) ? b : a )

/* FORTRAN interface
 * #define fpmetis_ FORTRAN_NAME(fpmetis_,FPMETIS_)
 */
#define fpmetis_part FORTRAN_NAME(fpmetis_part,FPMETIS_PART)
#define fpmetis_rpart FORTRAN_NAME(fpmetis_rpart,FPMETIS_RPART)
#define fpmetis_refine FORTRAN_NAME(fpmetis_refine,FPMETIS_REFINE)

/**
 *
 * @param fmpicomm
 * @param vtxdist
 * @param xadj
 * @param adjncy
 * @param nparts
 * @param part
 */
void fpmetis_part(MPI_Fint * fmpicomm, int * vtxdist,int * xadj, int * adjncy,
		int * nparts, int * part);

/**
 *
 * @param fmpicomm
 * @param vtxdist
 * @param xadj
 * @param adjncy
 * @param nparts
 * @param part
 * @param itrf
 */
void fpmetis_rpart(MPI_Fint * fmpicomm, int * vtxdist,int * xadj, int * adjncy,
		int * nparts, int * part, double * itrf);

/**
 *
 * @param fmpicomm
 * @param vtxdist
 * @param xadj
 * @param adjncy
 * @param nparts
 * @param part
 */
void fpmetis_refine(MPI_Fint * fmpicomm, int * vtxdist,int * xadj, int * adjncy,
		int * nparts, int * part);

#endif /* NEKP4EST_PMETIS_PMETIS_FWRAP_H_ */
