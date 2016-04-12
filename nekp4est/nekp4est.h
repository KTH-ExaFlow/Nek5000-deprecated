/** @file nekp4est.h
 *  @brief Definitions for C-Fortran77 interface in nekp4est
 *
 * @author Adam Peplinski
 * @date Feb 20, 2016
 *
 * @ingroup nekp4est
 */
/** @defgroup nekp4est Nonconforming Nek5000 with p4est
 *
 * AMR version of nek5000
 */

#ifndef NEKP4EST_H_
#define NEKP4EST_H_

/* check preprocessing flags */
#if !defined(N_DIM)
#error "Undefined N_DIM check makenek options"
#endif

#if N_DIM < 2 || N_DIM > 3
#error "Wrong N_DIM value"
#endif

#if !defined(N_PSCL)
#error "Undefined N_PSCL check makenek options"
#endif

#if N_PSCL < 0
#error "Wrong N_PSCL value"
#endif

/*
 * To define the way string is transferred from f77 to C
 * Important for number of strings transferred to C routine
 */
#undef LEN_FOLLOW_STRING

/** @defgroup refmarkers Refinement markers
 *
 * Numbers designating the type and action for refinement
 *
 * @ingroup nekp4est
 */
/* @{ @ingroup refmarkers */
#define NP4_RM_NONE       0     /**< No refinement action */
#define NP4_RM_H_REF      1     /**< Refine by splitting element */
#define NP4_RM_H_CRS     -1     /**< Coarsen by merging element */
#define NP4_RM_P_REF      2     /**< Refine by rising polynomial order in element */
#define NP4_RM_P_CRS     -2     /**< Coarsen by lowering polynomial order in element */
/* @} */

#endif /* NEKP4EST_H_ */
