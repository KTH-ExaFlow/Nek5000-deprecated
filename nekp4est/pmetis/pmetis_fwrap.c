/*
 * pmetis_fwrap.c
 * Fortran interface for parmetis library.
 *
 *  Created on: Mar 4, 2016
 *      Author: Adam Peplinski
 */
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

#include <parmetis.h>

#include "name.h"
#include "pmetis_fwrap.h"

/* Graph partitioning */
void fpmetis_part(MPI_Fint * fmpicomm, int * vtxdist,int * xadj, int * adjncy,
		int * nparts, int * part)
{
  int err, i;
  idx_t wgtflag=0;  // no weights
  idx_t numflag=0;  // C numbering
  idx_t ncon=1;     // number of vertex weights
  idx_t options[3]; // options
  idx_t edgcut;     // number of crossed edges
  real_t tpwgts[1][*nparts];
  real_t ubvec[1];
  real_t rtmp; 
  MPI_Comm mpicomm;

  // fraction of vertex weight
  rtmp = 1.0/(real_t)(*nparts);
  for(i=0;i<*nparts;++i){
    tpwgts[0][i] = rtmp;
  }

  // imbalance tolerance
  ubvec[0] = 1.05;

  // options; default
  options[0] = 0;

  // MPI communicator
  mpicomm = MPI_Comm_f2c(*fmpicomm);

  err = ParMETIS_V3_PartKway(vtxdist,xadj,adjncy,NULL,NULL,&wgtflag,&numflag,
		  &ncon,nparts,tpwgts,ubvec,options,&edgcut,part,&mpicomm);
}
/* Graph repartitioning */
void fpmetis_rpart(MPI_Fint * fmpicomm, int * vtxdist,int * xadj, int * adjncy,
		int * nparts, int * part, double * itrf)
{
	int err, i;
	idx_t wgtflag=0;   // no weights
	idx_t numflag=0;   // C numbering
	idx_t ncon=1;      // number of vertex weights
	idx_t options[4];  // options
	idx_t edgcut;      // number of crossed edges
	real_t tpwgts[1][*nparts];
	real_t ubvec[1];
	real_t itr;
	real_t rtmp;
	MPI_Comm mpicomm;

	// fraction of vertex weight
	rtmp = 1.0/(real_t)(*nparts);
	for(i=0;i<*nparts;++i){
		tpwgts[0][i] = rtmp;
	}

	// imbalance tolerance
	ubvec[0] = 1.05;

	// itr bounds
	itr = (real_t) *itrf;
	itr = MAX(itr, 1.0E-6);
	itr = MIN(itr, 1.0E6);

	// options; default
	//options[0] = 0;
	options [ 0 ] = 1; // set to zero for default
	options [ 1 ] = 0; // get timings
	options [ 2 ] = 15; // random seed
	options [ 3 ] = PARMETIS_PSR_UNCOUPLED; // sub-domains and processors are un-coupled

	// MPI communicator
	mpicomm = MPI_Comm_f2c(*fmpicomm);

	err = ParMETIS_V3_AdaptiveRepart(vtxdist,xadj,adjncy,NULL,NULL,NULL,&wgtflag,
			&numflag,&ncon,nparts,tpwgts,ubvec,&itr,options,&edgcut,part,&mpicomm);
}

void fpmetis_refine(MPI_Fint * fmpicomm, int * vtxdist,int * xadj, int * adjncy,
		int * nparts, int * part)
{
  int err, i;
  idx_t wgtflag=0;    // no weights
  idx_t numflag=0;    // C numbering
  idx_t ncon=1;       // number of vertex weights
  idx_t options[4];   // options
  idx_t edgcut;       // number of crossed edges
  real_t tpwgts[1][*nparts];
  real_t ubvec[1];
  real_t rtmp;
  MPI_Comm mpicomm;

  // fraction of vertex weight
  rtmp = 1.0/(real_t)(*nparts);
  for(i=0;i<*nparts;++i){
    tpwgts[0][i] = rtmp;
  }

  // imbalance tolerance
  ubvec[0] = 1.05;

  // options; default
  //options[0] = 0;
  options[0] = 1; // set to zero for default
  options[1] = 0; // get timings
  options[2] = 15; // random seed
  options[3] = PARMETIS_PSR_UNCOUPLED; // sub-domains and processors are un-coupled

  // MPI communicator
  mpicomm = MPI_Comm_f2c(*fmpicomm);

  err = ParMETIS_V3_RefineKway(vtxdist,xadj,adjncy,NULL,NULL,&wgtflag,&numflag,
		  &ncon,nparts,tpwgts,ubvec,options,&edgcut,part,&mpicomm);
}
