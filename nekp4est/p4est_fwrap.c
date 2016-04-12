/*
 * p4est_fwrap.c
 * Fortran interface for p4est library.
 *
 *  Created on: Feb 21, 2016
 *      Author: Adam Peplinski
 */
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

#include "nekp4est.h"

#if N_DIM == 2
#undef P4_TO_P8
#else
#include <p4est_to_p8est.h>
#endif

#ifndef P4_TO_P8
#include <p4est_algorithms.h>
#include <p4est_bits.h>
#include <p4est_extended.h>
#include <p4est_ghost.h>
#include <p4est_nodes.h>
#include <p4est_vtk.h>
#include <p4est_lnodes.h>
#include <p4est_mesh.h>
#include <p4est_iterate.h>
#else
#include <p8est_algorithms.h>
#include <p8est_bits.h>
#include <p8est_extended.h>
#include <p8est_ghost.h>
#include <p8est_nodes.h>
#include <p8est_vtk.h>
#include <p8est_lnodes.h>
#include <p8est_mesh.h>
#include <p8est_iterate.h>
#endif

#include "name.h"

#include "p4est_fwrap.h"

/** Global variables */
static p4est_connectivity_t *connect_nek=NULL; /**< Nek5000 connectivity structure */
static p4est_t              *tree_nek=NULL;    /**< Nek5000 tree structure */
static p4est_mesh_t         *mesh_nek=NULL;    /**< Nek5000 mesh structure */
static p4est_ghost_t        *ghost_nek=NULL;   /**< Nek5000 ghost zone structure */
static p4est_nodes_t        *nodes_nek=NULL;   /**< Nek5000 vertex numbering structure */
static p4est_lnodes_t	    *lnodes_nek=NULL;  /**< Nek5000 GLL node numbering structure */
static int                   nellv;            /**< number of local V-mesh elements */
static int                   max_level;        /**< local max quadrant level */

/*------------------------------------------------------------------------------
 * internal subroutines for interaction with Fortran
------------------------------------------------------------------------------*/
/** @brief Initialize element data for 0 level tree
 *
 * @details This is dummy routine required by fp4est_tree_new.
 * I do not initialize tree data here, as there are special routines dedicated
 * to transfer data between Nek5000 and p4est, however it could be used if
 * one would like to generate connectivity and tree directly in Nek5000
 *
 * @param p4est
 * @param which_tree
 * @param quadrant
 */
void
init_msh_dat (p4est_t * p4est, p4est_topidx_t which_tree,
		p4est_quadrant_t * quadrant)
{
}

/* Wrappers */
/* initialize */
void fp4est_init(int * log_threshold)
{
	p4est_init (NULL, *log_threshold);
}

/* Connectivity */
#ifdef P4_TO_P8
void fp4est_cnn_new(int * num_vertices, int * num_trees,
		int * num_edges,
		int * num_corners,
		double *vertices,
		int * tree_to_vertex, int * tree_to_tree,
		int * tree_to_face,
		int * tree_to_edge, int * ett_offset,
		int * edge_to_tree, int * edge_to_edge,
		int * tree_to_corner, int * ctt_offset,
		int * corner_to_tree, int * corner_to_corner)
#else
void fp4est_cnn_new(int * num_vertices, int * num_trees,
		int * num_corners,
		double *vertices,
		int * tree_to_vertex, int * tree_to_tree,
		int * tree_to_face,
		int * tree_to_corner, int * ctt_offset,
		int * corner_to_tree, int * corner_to_corner)
#endif
{
	int il;

	p4est_topidx_t num_verticesl, num_treesl;
#ifdef P4_TO_P8
	p4est_topidx_t num_edgesl;
#endif
	p4est_topidx_t num_cornersl;
	double verticesl[*num_vertices*3];
	p4est_topidx_t tree_to_vertexl[*num_trees*P4EST_CHILDREN];
	p4est_topidx_t tree_to_treel[*num_trees*P4EST_FACES];
	int8_t tree_to_facel[*num_trees*P4EST_FACES];
#ifdef P4_TO_P8
	p4est_topidx_t tree_to_edgel[*num_trees*P8EST_EDGES];
	p4est_topidx_t ett_offsetl[*num_edges+1];
	p4est_topidx_t edge_to_treel[ett_offset[*num_edges]];
	int8_t edge_to_edgel[ett_offset[*num_edges]];
#endif
	p4est_topidx_t tree_to_cornerl[*num_trees*P4EST_CHILDREN];
	p4est_topidx_t ctt_offsetl[*num_corners+1];
	p4est_topidx_t corner_to_treel[ctt_offset[*num_corners]];
	int8_t corner_to_cornerl[ctt_offset[*num_corners]];

	num_verticesl = (p4est_topidx_t) *num_vertices;
	num_treesl    = (p4est_topidx_t) *num_trees;
#ifdef P4_TO_P8
	num_edgesl    = (p4est_topidx_t) *num_edges;
#endif
	num_cornersl  = (p4est_topidx_t) *num_corners;

	for(il=0;il<*num_vertices*3;++il){
		verticesl[il] = (double) vertices[il];
	}

	for(il=0;il<*num_trees*P4EST_CHILDREN;++il){
		tree_to_vertexl[il] = (p4est_topidx_t) tree_to_vertex[il];
	}
	for(il=0;il<*num_trees*P4EST_FACES;++il){
		tree_to_treel[il] = (p4est_topidx_t) tree_to_tree[il];
		tree_to_facel[il] = (int8_t) tree_to_face[il];
	}

#ifdef P4_TO_P8
	for(il=0;il<*num_trees*P8EST_EDGES;++il){
		tree_to_edgel[il] = (p4est_topidx_t) tree_to_edge[il];
	}
	for(il=0;il<= *num_edges;++il){
		ett_offsetl[il] = (p4est_topidx_t) ett_offset[il];
	}
	for(il=0;il<ett_offset[*num_edges];++il){
		edge_to_treel[il] = (p4est_topidx_t) edge_to_tree[il];
		edge_to_edgel[il] = (int8_t) edge_to_edge[il];
	}
#endif

	for(il=0;il<*num_trees*P4EST_CHILDREN;++il){
		tree_to_cornerl[il] = (p4est_topidx_t) tree_to_corner[il];
	}
	for(il=0;il<= *num_corners;++il){
		ctt_offsetl[il] = (p4est_topidx_t) ctt_offset[il];
	}
	for(il=0;il<ctt_offset[*num_corners];++il){
		corner_to_treel[il] = (p4est_topidx_t) corner_to_tree[il];
		corner_to_cornerl[il] = (int8_t) corner_to_corner[il];
	}

#ifdef P4_TO_P8
	connect_nek = p4est_connectivity_new_copy(num_verticesl, num_treesl,
			num_edgesl,
			num_cornersl,
			verticesl, tree_to_vertexl,
			tree_to_treel, tree_to_facel,
			tree_to_edgel, ett_offsetl,
			edge_to_treel, edge_to_edgel,
			tree_to_cornerl, ctt_offsetl,
			corner_to_treel, corner_to_cornerl);
#else
	connect_nek = p4est_connectivity_new_copy(num_verticesl, num_treesl,
			num_cornersl,
			verticesl, tree_to_vertexl,
			tree_to_treel, tree_to_facel,
			tree_to_cornerl, ctt_offsetl,
			corner_to_treel, corner_to_cornerl);
#endif
}

void fp4est_cnn_del()
{
	if (connect_nek) {
		p4est_connectivity_destroy(connect_nek);
	}
}

void fp4est_cnn_attr(int * enable_tree_attr)
{
	p4est_connectivity_set_attr (connect_nek, *enable_tree_attr);
}

void fp4est_cnn_valid(int * is_valid)
{
	*is_valid = p4est_connectivity_is_valid(connect_nek);
}

void fp4est_cnn_complete()
{
	p4est_connectivity_complete(connect_nek);
}

void fp4est_cnn_save(char *filename, int len_f)
{
	p4est_connectivity_save (filename,connect_nek);
}

void fp4est_cnn_load(char * filename, int len_f)
{
	connect_nek = p4est_connectivity_load (filename, NULL);
}

/* tree_ management */
void fp4est_tree_new(MPI_Fint * fmpicomm, int * min_level)
{
	MPI_Comm mpicomm;
	mpicomm = MPI_Comm_f2c(*fmpicomm);
	tree_nek = p4est_new_ext (mpicomm,connect_nek, 0, *min_level, 1,
			sizeof(user_data_t), init_msh_dat, NULL);
}

void fp4est_tree_del()
{
	if (tree_nek) {
		p4est_destroy (tree_nek);
	}
}

void fp4est_tree_valid(int * is_valid)
{
	*is_valid = p4est_is_valid(tree_nek);
}

void fp4est_tree_save(int *save_data, char * filename, int len_f)
{
	p4est_save (filename, tree_nek, *save_data);
}

void fp4est_tree_load(MPI_Fint * fmpicomm, int *load_data, char * filename, int len_f)
{
	MPI_Comm mpicomm;
	mpicomm = MPI_Comm_f2c(*fmpicomm);
	tree_nek = p4est_load (filename, mpicomm, sizeof(user_data_t), *load_data,
			NULL, &connect_nek);
}

/* tree and grid info */
void fp4est_ghost_new()
{
	ghost_nek = p4est_ghost_new(tree_nek,P4EST_CONNECT_FULL);
}

void fp4est_ghost_del()
{
	if (ghost_nek) {
		p4est_ghost_destroy(ghost_nek);
	}
}

void fp4est_mesh_new()
{
	mesh_nek = p4est_mesh_new(tree_nek,ghost_nek,P4EST_CONNECT_FULL);
}

void fp4est_mesh_del()
{
	if (mesh_nek) {
		p4est_mesh_destroy(mesh_nek);
	}
}

void fp4est_nodes_new()
{
	nodes_nek = p4est_nodes_new(tree_nek,ghost_nek);
}

void fp4est_nodes_del()
{
	if (nodes_nek) {
		p4est_nodes_destroy(nodes_nek);
	}
}

void fp4est_lnodes_new(int ldgr)
{
	lnodes_nek = p4est_lnodes_new(tree_nek,ghost_nek,ldgr);
}

void fp4est_lnodes_del()
{
	if (lnodes_nek) {
		p4est_lnodes_destroy(lnodes_nek);
	}
}

/* nekp4est internal load balance */
void fp4est_part(int * partforcoarsen)
{
        p4est_partition (tree_nek,*partforcoarsen,NULL);
}

/* refinement, coarsening, balance */

/* I/O (VTK) */
void fp4est_vtk_write(char * filename, int len_f)
{
	p4est_vtk_write_file (tree_nek, NULL, filename);
}

void fp4est_vtk_iscalar(int * iscalar,int *num,char * filename, int len_f)
{
	int il,jl,err;
	double rscalar[*num*P4EST_CHILDREN];

	for(il=0;il<*num;++il){
		for(jl=0;jl<P4EST_CHILDREN;jl++){
			rscalar[il*P4EST_CHILDREN+jl] = (double) iscalar[il];
		}
	}

	err = p4est_vtk_write_header (tree_nek, NULL, 1.0, 0, 0, 0, 0, "scalar", NULL, filename);
	err = p4est_vtk_write_point_scalar (tree_nek, NULL, filename, "scalar", rscalar);
	err = p4est_vtk_write_footer (tree_nek, filename);
}

/* routines required by Nek5000 for data exchange */

/** @brief Iterate over elements to count V-mesh elements
 *
 * @details Required by fp4est_msh_get_size
 *
 * @param info
 * @param user_data
 */
void count_mshv(p4est_iter_volume_info_t * info, void *user_data)
{
	int loc_level;
	user_data_t  *data = (user_data_t *) info->quad->p.user_data;

	if (data->imsh == 0){
	  nellv = nellv+1;
	}
    // find max local level
    loc_level = (int) info->quad->level;
	max_level = (loc_level > max_level ? loc_level : max_level);
}

/* get mesh size */
void fp4est_msh_get_size(int * nelgt,int * nelgit, int * nelt, int * nelv, int * maxl)
{
  // get global number of quadrants
  *nelgt = (int) tree_nek->global_num_quadrants;
  // zero based global position of local quadrants
  *nelgit = (int)
		tree_nek->global_first_quadrant[tree_nek->mpirank];
  // number of local quadrants
  *nelt = (int)tree_nek->local_num_quadrants;

  // count number of V-mesh elements and find current max level
  nellv = 0;
  max_level = 0;

#ifdef P4_TO_P8
    p4est_iterate(tree_nek,ghost_nek,NULL,count_mshv,NULL,NULL,NULL);
#else
    p4est_iterate(tree_nek,ghost_nek,NULL,count_mshv,NULL,NULL);
#endif

  *nelv = nellv;
  *maxl = max_level;
}

/** @brief Iterate over element volumes to transfer element data
 *
 * @details Required by fp4est_msh_get
 *
 * @param info
 * @param user_data
 */
void iter_msh_dat(p4est_iter_volume_info_t * info, void *user_data)
{
	user_data_t        *data = (user_data_t *) info->quad->p.user_data;

	// which quad (local and global element number)
	p4est_tree_t       *tree;
	p4est_locidx_t      iwl;
	int                 iwlt, iwg, level;

	// to get element vertices
	double              xyz[3*P4EST_CHILDREN];
	int                 ifc, ifl, ic, ref, set;
	int8_t              ifm, ifmt;
	p4est_locidx_t      icl;

	// to correct element connectivity
	// I check only E (element) and P (periodic) b.c.
	p4est_locidx_t     *quad_to_quad = (p4est_locidx_t *) mesh_nek->quad_to_quad;
	int8_t             *quad_to_face = (int8_t *) mesh_nek->quad_to_face;
	// quads on different processor
	p4est_ghost_t      *ghost_layer = (p4est_ghost_t *) info->ghost_layer;
	p4est_quadrant_t   *quad;
	int                *ghost_to_proc = (int *) mesh_nek->ghost_to_proc;

	// external function to collect element data
	extern void nek_get_msh_dat(int * iwt, int * iwq, int * imsh, int * igrp,
					int * level, char (*)[N_PSCL+2][6][3], double (*)[N_PSCL+2][6][5],
				    int (*)[6]);
	// external function to collect mesh refinement history
	extern void nek_get_msh_hst(int * iwq_o, int * iwq_p, int (*) [P4EST_CHILDREN]);

	// common blocks to exchange data with Nek5000
	extern struct
	{
	  int ibc, nfldt;
	}
	nekp4est_cbci;

	// get quad number
	tree = p4est_tree_array_index (info->p4est->trees, info->treeid);
	// local quad number
	iwl = info->quadid + tree->quadrants_offset;
	iwlt = (int) iwl;
	// global quad number
	iwg = (int) tree_nek->global_first_quadrant[tree_nek->mpirank] + (int) iwl;

	// quadrant level
	level = (int) info->quad->level;

	// check connectivity
	for(ifc=0;ifc<P4EST_FACES;++ifc){
	  ic = (int) iwl*P4EST_FACES+ifc;
	  // is it inner or periodic face
	  if (iwl == quad_to_quad[ic] && ifc == (int) quad_to_face[ic]){
		  // not inner or periodic face
		  for (ifl=nekp4est_cbci.ibc;ifl<=nekp4est_cbci.nfldt;++ifl){
			  //check connection type (should not be E or P)
			  if(data->cbc[ifl][ifc][2]==' '&&data->cbc[ifl][ifc][1]==' '&&
					  (data->cbc[ifl][ifc][0]=='E'||data->cbc[ifl][ifc][0]=='P')){
				  printf("Connectivity error: %i %i %i %i %c%c%c\n",
						  tree_nek->mpirank,iwl,ifc,ifl,data->cbc[ifl][ifc][0],
						  data->cbc[ifl][ifc][1],data->cbc[ifl][ifc][2]);
				  printf("internal/external face mismatch.\n");
				  SC_ABORT("Aborting: iter_msh_dat\n");
		      }
		  }
	  }
	  else{
		  // inner or periodic face
		  for (ifl=nekp4est_cbci.ibc;ifl<=nekp4est_cbci.nfldt;++ifl){
			  // check connection type (should be E or P)
			  // this check should eliminate external faces here
			  if(data->cbc[ifl][ifc][2]==' '&&data->cbc[ifl][ifc][1]==' '&&
					  (data->cbc[ifl][ifc][0]=='E'||data->cbc[ifl][ifc][0]=='P')){
				  //set neighbour global number
				  icl = quad_to_quad[ic];
				  if(icl < mesh_nek->local_num_quadrants){
					  // local quad
					  data->bc[ifl][ifc][0] = (double) icl +
							  tree_nek->global_first_quadrant[tree_nek->mpirank]+ 1;
				  }
				  else{
					  // remote quad
					  icl = icl - mesh_nek->local_num_quadrants;
					  quad = p4est_quadrant_array_index (&ghost_layer->ghosts, (size_t) icl);
					  data->bc[ifl][ifc][0] = (double) quad->p.piggy3.local_num +
							  tree_nek->global_first_quadrant[ghost_to_proc[icl]]+ 1;
				  }
				  //set neighbour face
				  ifm = quad_to_face[ic];
#ifdef P4_TO_P8
				  if (ifm>=0&&ifm<=23){
					  // equal face sizes
					  data->bc[ifl][ifc][4] = 0.0;
					  // position
					  data->bc[ifl][ifc][3] = 0.0;
					  // orientation
					  ifmt = ifm / P4EST_FACES; // orientation r
					  ifm = ifm % P4EST_FACES;  // neighbour face
					  // find permutation
					  ref = p8est_face_permutation_refs[ifc][ifm];
					  set = p8est_face_permutation_sets[ref][ifmt];
					  data->bc[ifl][ifc][2] = (double) set;
					  // face number
					  data->bc[ifl][ifc][1] = (double) 1+ ifm;
				  }
				  else if(ifm>=24&&ifm<=119){
					  // Double-sized neighbour; marked 'J  ' in original nek5000
					  // I mark them 'E  ' or 'P  ' setting last parameter in bc array to 1
					  // important for subroutines:
					  //		get_fast_bc (no modification necessary)
					  //		and in dssum.f
					  // double-sized neighbour
					  data->bc[ifl][ifc][4] = 1.0;
					  // position
					  ifmt = ifm / 24;
					  data->bc[ifl][ifc][3] = (double) ifmt - 1;
					  // orientation
					  ifm = ifm - 24*ifmt;
					  ifmt = ifm / P4EST_FACES; // orientation r
					  ifm = ifm % P4EST_FACES;  // neighbour face
					  // find permutation
					  ref = p8est_face_permutation_refs[ifc][ifm];
					  set = p8est_face_permutation_sets[ref][ifmt];
					  data->bc[ifl][ifc][2] = (double) set;
					  // face number
					  data->bc[ifl][ifc][1] = (double) 1+ ifm;
				  }
				  else if(ifm>=-24&&ifm<=-1){
					  // Half-sized neighbour; marked 'SP ' in original nek5000
					  // I mark them 'E  ' or 'P  ' setting last parameter in bc array to 1
					  // important for subroutines:
					  //		get_fast_bc (no modification necessary)
					  //		and in dssum.f
					  // half-sized neighbour
					  data->bc[ifl][ifc][4] = 2.0;
					  // position
					  data->bc[ifl][ifc][3] = 0.0;
					  // orientation
					  ifm = ifm + 24;
					  ifmt = ifm / P4EST_FACES; // orientation r
					  ifm = ifm % P4EST_FACES;  // neighbour face
					  // find permutation
					  ref = p8est_face_permutation_refs[ifc][ifm];
					  set = p8est_face_permutation_sets[ref][ifmt];
					  data->bc[ifl][ifc][2] = (double) set;
					  // face number
					  data->bc[ifl][ifc][1] = (double) 1+ ifm;
					  // element number
					  data->bc[ifl][ifc][0] = 0.0;
				  }
				  else{
					  SC_ABORT("Wrong face number, aborting: iter_mshv\n");
				  }
#else
				  if (ifm>=0&&ifm<=7){
					  // equal face sizes
					  data->bc[ifl][ifc][4] = 0.0;
					  // position
					  data->bc[ifl][ifc][3] = 0.0;
					  // orientation
					  ifmt = ifm / P4EST_FACES;
					  data->bc[ifl][ifc][2] = (double) ifmt;
					  // face number
					  data->bc[ifl][ifc][1] = (double) 1+ ifm % P4EST_FACES;
				  }
				  else if(ifm>=8&&ifm<=23){
					  // Double-sized neighbour; marked 'J  ' in original nek5000
					  // I mark them 'E  ' or 'P  ' setting last parameter in bc array to 1
					  // important for subroutines:
					  //		get_fast_bc (no modification necessary)
					  //		and in dssum.f
					  // double-sized neighbour
					  data->bc[ifl][ifc][4] = 1.0;
					  // position
					  ifmt = ifm / 8;
					  data->bc[ifl][ifc][3] = (double) ifmt -1;
					  // orientation
					  ifm = ifm - 8*ifmt;
					  ifmt = ifm / P4EST_FACES;
					  data->bc[ifl][ifc][2] = (double) ifmt;
					  // face number
					  data->bc[ifl][ifc][1] = (double) 1+ ifm % P4EST_FACES;
				  }
				  else if(ifm>=-8&&ifm<=-1){
					  // Half-sized neighbour; marked 'SP ' in original nek5000
					  // I mark them 'E  ' or 'P  ' setting last parameter in bc array to 1
					  // important for subroutines:
					  //		get_fast_bc (no modification necessary)
					  //		and in dssum.f
					  // half-sized neighbour
					  data->bc[ifl][ifc][4] = 2.0;
					  // position
					  data->bc[ifl][ifc][3] = 0.0;
					  // orientation
					  ifm = ifm + 8;
					  ifmt = ifm / P4EST_FACES;
					  data->bc[ifl][ifc][2] = (double) ifmt;
					  // face number
					  data->bc[ifl][ifc][1] = (double) 1+ ifm % P4EST_FACES;
					  // element number
					  data->bc[ifl][ifc][0] = 0.0;
				  }
				  else{
					  SC_ABORT("Wrong face number, aborting: iter_mshv\n");
				  }
#endif
			  }
			  else{
				  // not corretly marked inner face
				  printf("Connectivity error: %i %i %i %i %c%c%c\n",
						  tree_nek->mpirank,iwl,ifc,ifl,data->cbc[ifl][ifc][0],
						  data->cbc[ifl][ifc][1],data->cbc[ifl][ifc][2]);
				  SC_ABORT("Aborting: iter_msh_dat\n");
			  }
		  } // loop over fluids
	  }
	} // loop over faces

	// transfer mesh data to Nek5000
	nek_get_msh_dat(&iwlt,&iwg,&data->imsh,&data->igrp,&level,
			&data->cbc,&data->bc,&data->crv);

	// transfer refinement history to Nek5000
	nek_get_msh_hst(&data->gln_el,&data->gln_parent,&data->gln_children);

	//	fill current global numbering of elements in nek5000 and reset refinement history for 0 level blocs
	data->gln_el = iwg;
	data->gln_parent = -1;
	for(ifc=0;ifc<P4EST_CHILDREN;++ifc){
		data->gln_children[ifc] = -1;
	}
}

// get mesh data to Nek5000
void fp4est_msh_get_dat()
{
#ifdef P4_TO_P8
  p4est_iterate(tree_nek,ghost_nek,NULL,iter_msh_dat,NULL,NULL,NULL);
#else
  p4est_iterate(tree_nek,ghost_nek,NULL,iter_msh_dat,NULL,NULL);
#endif
}


// get global vertex numbering
void fp4est_msh_get_node(int * lnelt, int * node)
{
  // number of processors and processor id.
  const int num_procs = tree_nek->mpisize;
  const int rank = tree_nek->mpirank;
  // types inside nodes structure
  // independent nodes
  p4est_indep_t *indep;
  //hanging nodes
#ifdef P4_TO_P8
  p8est_hang2_t *hang2;
  p8est_hang4_t *hang4;
#else
  p4est_hang2_t *hang2;
#endif
  // numbers of independent and hanging nodes
  const int ni = (int) nodes_nek->indep_nodes.elem_count;
  const int fi = (int) nodes_nek->face_hangings.elem_count;
#ifdef P4_TO_P8
  const int ei = (int) nodes_nek->edge_hangings.elem_count;
#endif
  // number of local elements
  const int vi = (int) nodes_nek->num_local_quadrants;
  // local node offset
  const int offi = (int) nodes_nek->offset_owned_indeps;
  int il, jl, kl, it, nl, nr;
  p4est_locidx_t  offset;

  // loca number of quadrants
  *lnelt = (int) nodes_nek->num_local_quadrants;

  // quad to vertex global map
  for (il=0;il<vi;++il){
    for(jl=0;jl<P4EST_CHILDREN;++jl){

      nl = (int) nodes_nek->local_nodes[il*P4EST_CHILDREN+jl];
      if (nl < ni){
          // independent node
    	  // do nothing
      }
      else if(nl < ni+fi){
          // face hanging
    	  // find local vertex number on the face
    	  SC_ABORT("Hanging node; not supported yet\n");
      }
#ifdef P4_TO_P8
      else if(nl < ni+fi+ei){
          // edge hanging
    	  // find local vertex number on the edge
    	  SC_ABORT("Hangign node; not supported yet\n");
      }
#endif
      else{
    	  SC_ABORT("Wrong node number, aborting: fp4est_msh_get_node\n");
      }
      // get global number of this node
	  indep = (p4est_indep_t *) sc_array_index (&nodes_nek->indep_nodes,nl);
	  if (nl < offi) {
		  nr = nodes_nek->nonlocal_ranks[nl];
	  }
	  else{
		  nr = rank;
	  }
	  offset = 0;
	  for(kl=0;kl<nr;++kl){
		  offset = offset + nodes_nek->global_owned_indeps[kl];
	  }
	  node[il*(P4EST_CHILDREN+2)+2+jl] = (int) offset + indep->p.piggy3.local_num +1;

    }
  }
}
