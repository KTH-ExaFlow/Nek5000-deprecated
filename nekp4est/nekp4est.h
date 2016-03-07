/** @file nekp4est.h
 *  @brief Definitions for C-Fortran77 interface in nekp4est
 *
 * @author Adam Peplinski
 * @date Feb 20, 2016
 *
 * @ingroup nekp4est
 */
/** @defgroup nekp4est nekp4est
 *
 * AMR version of nek5000
 */

#ifndef NEKP4EST_H_
#define NEKP4EST_H_
/*
 * To define the way string is transferred from f77 to C
 * Important for number of strings transferred to C routine
 */
#undef LEN_FOLLOW_STRING

/* check preprocessing flags */
#if !defined(N_DIM)
#error "Undefined N_DIM check makenek options"
#endif

#if !defined(N_PSCL)
#error "Undefined N_PSCL check makenek options"
#endif

#endif /* NEKP4EST_H_ */
