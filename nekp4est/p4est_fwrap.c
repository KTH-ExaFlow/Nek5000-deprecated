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
static p4est_connectivity_t *connect_nek = NULL; /**< Nek5000 connectivity structure */
static p4est_t *tree_nek = NULL; /**< Nek5000 tree structure */
static p4est_mesh_t *mesh_nek = NULL; /**< Nek5000 mesh structure */
static p4est_ghost_t *ghost_nek = NULL; /**< Nek5000 ghost zone structure */
static p4est_nodes_t *nodes_nek = NULL; /**< Nek5000 vertex numbering structure */
static p4est_lnodes_t *lnodes_nek = NULL; /**< Nek5000 GLL node numbering structure */
static user_data_t *ghost_data = NULL; /**< user data for ghost cells */

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
void init_msh_dat(p4est_t * p4est, p4est_topidx_t which_tree,
p4est_quadrant_t * quadrant) {
}

/* Wrappers */
/* initialize */
void fp4est_init(int * log_threshold) {
	p4est_init(NULL, *log_threshold);
}

/* Connectivity */
#ifdef P4_TO_P8
void fp4est_cnn_new(int * num_vertices, int * num_trees, int * num_edges,
		int * num_corners, double *vertices, int * tree_to_vertex,
		int * tree_to_tree, int * tree_to_face, int * tree_to_edge,
		int * ett_offset, int * edge_to_tree, int * edge_to_edge,
		int * tree_to_corner, int * ctt_offset, int * corner_to_tree,
		int * corner_to_corner)
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
	double verticesl[*num_vertices * 3];
	p4est_topidx_t tree_to_vertexl[*num_trees * P4EST_CHILDREN];
	p4est_topidx_t tree_to_treel[*num_trees * P4EST_FACES];
	int8_t tree_to_facel[*num_trees * P4EST_FACES];
#ifdef P4_TO_P8
	p4est_topidx_t tree_to_edgel[*num_trees * P8EST_EDGES];
	p4est_topidx_t ett_offsetl[*num_edges + 1];
	p4est_topidx_t edge_to_treel[ett_offset[*num_edges]];
	int8_t edge_to_edgel[ett_offset[*num_edges]];
#endif
	p4est_topidx_t tree_to_cornerl[*num_trees * P4EST_CHILDREN];
	p4est_topidx_t ctt_offsetl[*num_corners + 1];
	p4est_topidx_t corner_to_treel[ctt_offset[*num_corners]];
	int8_t corner_to_cornerl[ctt_offset[*num_corners]];

	num_verticesl = (p4est_topidx_t) *num_vertices;
	num_treesl = (p4est_topidx_t) *num_trees;
#ifdef P4_TO_P8
	num_edgesl = (p4est_topidx_t) *num_edges;
#endif
	num_cornersl = (p4est_topidx_t) *num_corners;

	for (il = 0; il < *num_vertices * 3; ++il) {
		verticesl[il] = (double) vertices[il];
	}

	for (il = 0; il < *num_trees * P4EST_CHILDREN; ++il) {
		tree_to_vertexl[il] = (p4est_topidx_t) tree_to_vertex[il];
	}
	for (il = 0; il < *num_trees * P4EST_FACES; ++il) {
		tree_to_treel[il] = (p4est_topidx_t) tree_to_tree[il];
		tree_to_facel[il] = (int8_t) tree_to_face[il];
	}

#ifdef P4_TO_P8
	for (il = 0; il < *num_trees * P8EST_EDGES; ++il) {
		tree_to_edgel[il] = (p4est_topidx_t) tree_to_edge[il];
	}
	for (il = 0; il <= *num_edges; ++il) {
		ett_offsetl[il] = (p4est_topidx_t) ett_offset[il];
	}
	for (il = 0; il < ett_offset[*num_edges]; ++il) {
		edge_to_treel[il] = (p4est_topidx_t) edge_to_tree[il];
		edge_to_edgel[il] = (int8_t) edge_to_edge[il];
	}
#endif

	for (il = 0; il < *num_trees * P4EST_CHILDREN; ++il) {
		tree_to_cornerl[il] = (p4est_topidx_t) tree_to_corner[il];
	}
	for (il = 0; il <= *num_corners; ++il) {
		ctt_offsetl[il] = (p4est_topidx_t) ctt_offset[il];
	}
	for (il = 0; il < ctt_offset[*num_corners]; ++il) {
		corner_to_treel[il] = (p4est_topidx_t) corner_to_tree[il];
		corner_to_cornerl[il] = (int8_t) corner_to_corner[il];
	}

#ifdef P4_TO_P8
	connect_nek = p4est_connectivity_new_copy(num_verticesl, num_treesl,
			num_edgesl, num_cornersl, verticesl, tree_to_vertexl, tree_to_treel,
			tree_to_facel, tree_to_edgel, ett_offsetl, edge_to_treel,
			edge_to_edgel, tree_to_cornerl, ctt_offsetl, corner_to_treel,
			corner_to_cornerl);
#else
	connect_nek = p4est_connectivity_new_copy(num_verticesl, num_treesl,
			num_cornersl,
			verticesl, tree_to_vertexl,
			tree_to_treel, tree_to_facel,
			tree_to_cornerl, ctt_offsetl,
			corner_to_treel, corner_to_cornerl);
#endif
}

void fp4est_cnn_del() {
	if (connect_nek) p4est_connectivity_destroy(connect_nek);
}

void fp4est_cnn_attr(int * enable_tree_attr) {
	p4est_connectivity_set_attr(connect_nek, *enable_tree_attr);
}

void fp4est_cnn_valid(int * is_valid) {
	*is_valid = p4est_connectivity_is_valid(connect_nek);
}

void fp4est_cnn_complete() {
	p4est_connectivity_complete(connect_nek);
}

void fp4est_cnn_save(char *filename, int len_f) {
	p4est_connectivity_save(filename, connect_nek);
}

void fp4est_cnn_load(char * filename, int len_f) {
	connect_nek = p4est_connectivity_load(filename, NULL);
}

/* tree_ management */
void fp4est_tree_new(MPI_Fint * fmpicomm, int * min_level) {
	MPI_Comm mpicomm;
	mpicomm = MPI_Comm_f2c(*fmpicomm);
	tree_nek = p4est_new_ext(mpicomm, connect_nek, 0, *min_level, 1,
			sizeof(user_data_t), init_msh_dat, NULL);
}

void fp4est_tree_del() {
	if (tree_nek) p4est_destroy(tree_nek);
}

void fp4est_tree_valid(int * is_valid) {
	*is_valid = p4est_is_valid(tree_nek);
}

void fp4est_tree_save(int *save_data, char * filename, int len_f) {
	p4est_save_ext(filename, tree_nek, *save_data,0);
}

void fp4est_tree_load(MPI_Fint * fmpicomm, int *load_data, char * filename,
		int len_f) {
	MPI_Comm mpicomm;
	mpicomm = MPI_Comm_f2c(*fmpicomm);
	tree_nek = p4est_load_ext(filename, mpicomm, sizeof(user_data_t), *load_data,1,0,
	NULL, &connect_nek);
}

/* tree and grid info */
void fp4est_ghost_new() {
	if (ghost_nek) p4est_ghost_destroy(ghost_nek);
	ghost_nek = p4est_ghost_new(tree_nek, P4EST_CONNECT_FULL);
	/* ghost data */
	if (ghost_data) P4EST_FREE(ghost_data);
	ghost_data = P4EST_ALLOC (user_data_t, ghost_nek->ghosts.elem_count);
	p4est_ghost_exchange_data (tree_nek, ghost_nek, ghost_data);
}

void fp4est_ghost_del() {
	if (ghost_nek) p4est_ghost_destroy(ghost_nek);
	if (ghost_data) P4EST_FREE(ghost_data);
}

void fp4est_mesh_new() {
	if (mesh_nek) p4est_mesh_destroy(mesh_nek);
	mesh_nek = p4est_mesh_new(tree_nek, ghost_nek, P4EST_CONNECT_FULL);
}

void fp4est_mesh_del() {
	if (mesh_nek) p4est_mesh_destroy(mesh_nek);
}

void fp4est_nodes_new() {
	if (nodes_nek) p4est_nodes_destroy(nodes_nek);
	nodes_nek = p4est_nodes_new(tree_nek, ghost_nek);
}

void fp4est_nodes_del() {
	if (nodes_nek) p4est_nodes_destroy(nodes_nek);
}

void fp4est_lnodes_new() {
	int ldgr;
	/* I set degree to -N_DIM to be able to distinguish between vertices, edges and faces.
	 * Element interior is discarded. */
	ldgr = -N_DIM;
	if (lnodes_nek) p4est_lnodes_destroy(lnodes_nek);
	lnodes_nek = p4est_lnodes_new(tree_nek, ghost_nek, ldgr);
}

void fp4est_lnodes_del() {
	if (lnodes_nek) p4est_lnodes_destroy(lnodes_nek);
}

/* nekp4est internal load balance */
void fp4est_part(int * partforcoarsen) {
	p4est_partition(tree_nek, *partforcoarsen, NULL);
}

/* refinement, coarsening, balance */

/* I/O (VTK) */
void fp4est_vtk_write(char * filename, int len_f) {
	p4est_vtk_write_file(tree_nek, NULL, filename);
}

void fp4est_vtk_iscalar(int * iscalar, int *num, char * filename, int len_f) {
	int il, jl, err;
	double rscalar[*num * P4EST_CHILDREN];

	for (il = 0; il < *num; ++il) {
		for (jl = 0; jl < P4EST_CHILDREN; jl++) {
			rscalar[il * P4EST_CHILDREN + jl] = (double) iscalar[il];
		}
	}

	err = p4est_vtk_write_header(tree_nek, NULL, 1.0, 0, 0, 0, 0, "scalar",
	NULL, filename);
	err = p4est_vtk_write_point_scalar(tree_nek, NULL, filename, "scalar",
			rscalar);
	err = p4est_vtk_write_footer(tree_nek, filename);
}

/* routines required by Nek5000 for data exchange */

/** @brief Iterate over elements to count V-mesh elements
 *
 * @details Required by fp4est_msh_get_size
 *
 * @param info
 * @param user_data
 */
void count_mshv(p4est_iter_volume_info_t * info, void *user_data) {
	int loc_level;
	user_data_t *data = (user_data_t *) info->quad->p.user_data;
	int *lmax = (int *) user_data;

    // coult V-type elements
	if (data->imsh == 0) {
		lmax[0] = lmax[0] + 1;
	}
	// find max local level
	loc_level = (int) info->quad->level;
	lmax[1] = (loc_level > lmax[1] ? loc_level : lmax[1]);
}

/* get mesh size */
void fp4est_msh_get_size(int * nelgt, int * nelgit, int * nelt,
		int * nelv, int * maxl) {
	int lmax[2];
	// get global number of quadrants
	*nelgt = (int) tree_nek->global_num_quadrants;
	// zero based global position of local quadrants
	*nelgit = (int) tree_nek->global_first_quadrant[tree_nek->mpirank];
	// number of local quadrants
	*nelt = (int) tree_nek->local_num_quadrants;

	// count number of V-mesh elements and find current max level
	lmax[0] = 0;
	lmax[1] = 0;

#ifdef P4_TO_P8
	p4est_iterate(tree_nek, ghost_nek,(void *) &lmax, count_mshv,
			NULL, NULL, NULL);
#else
	p4est_iterate(tree_nek, ghost_nek,(void *) &lmax, count_mshv,
			NULL, NULL);
#endif

	*nelv = lmax[0];
	*maxl = lmax[1];
}

/* data type for mesh data transfer between nek5000 and p4est*/
typedef struct transfer_data_s {
	int  ibc, ebc; /**< start and end for */
	int  nelv; /**< global number of V-type elements */
	int  lelt; /**< array dimension for bc arrays */
    int *level; /**< pointer to element level array */
	int *igrp; /**< pointer to element group array */
	int *crv; /**< curvature data */
	double *bc; /**< boundary condition parameters */
	char *cbc; /**< boundary condition type */

} transfer_data_t;

/** @brief Iterate over element volumes to transfer element mesh data
 *
 * @details Required by fp4est_msh_get_dat
 *
 * @param info
 * @param user_data
 */
void iter_msh_dat(p4est_iter_volume_info_t * info, void *user_data) {
	user_data_t *data = (user_data_t *) info->quad->p.user_data;
	transfer_data_t *trans_data = (transfer_data_t *) user_data;

	// which quad (local and global element number)
	p4est_tree_t *tree;
	p4est_locidx_t iwl;
	int iwlt, iwg;
	int ifc, ifl, ib, ic, il;// loop index

	// get quad number
	tree = p4est_tree_array_index(info->p4est->trees, info->treeid);
	// local quad number
	iwl = info->quadid + tree->quadrants_offset;
	iwlt = (int) iwl;
	// global quad number
	iwg = (int) tree_nek->global_first_quadrant[tree_nek->mpirank] + iwlt;

	// test V versus T meshes; no V mesh beyond nelv
	if (data->imsh && iwg < trans_data->nelv) {
		printf("V/T-type element mismatch: %i, %i\n",iwg,trans_data->nelv);
		SC_ABORT("Aborting: iter_msh_dat\n");
	}

	// quadrant level
	trans_data->level[iwlt] = (int) info->quad->level;

	// element group mark
	trans_data->igrp[iwlt] = data->igrp;

	// curvature data
	for (ic = 0; ic < 6; ic++) {
		trans_data->crv[iwlt*6+ic] = data->crv[ic];
	}

	// boundary condition
	// loop over fluids
	for (ifl = trans_data->ibc; ifl <= trans_data->ebc; ++ifl) {
		// skip velocity bc for T-type mesh
		if (data->imsh && ifl == NP4_VFLD) continue;
		// loop over faces
		for (ifc = 0; ifc < P4EST_FACES; ifc++) {
			ib = ((ifl*trans_data->lelt + iwlt)*6+ifc);
			ic = ib*5;
			for (il = 0; il < 5; il++) {
				trans_data->bc[ic+il] = data->bc[ifl][ifc][il];
			}
			ic = ib*3;
			for (il = 0; il < 3; il++) {
				trans_data->cbc[ic+il] = data->cbc[ifl][ifc][il];
			}
		}
	}
}

// get mesh data to Nek5000
void fp4est_msh_get_dat(int * ibc, int * ebc, int * nelv, int *lelt,
		int * igrp, int * level, int * crv, double * bc, char * cbc) {
	transfer_data_t transfer_data;
	transfer_data.ibc = *ibc;
	transfer_data.ebc = *ebc;
	transfer_data.nelv = *nelv;
	transfer_data.lelt = *lelt;
	transfer_data.igrp = igrp;
	transfer_data.level = level;
	transfer_data.crv = crv;
	transfer_data.bc = bc;
	transfer_data.cbc = cbc;
#ifdef P4_TO_P8
	p4est_iterate(tree_nek, ghost_nek,(void *) &transfer_data,
			iter_msh_dat, NULL, NULL, NULL);
#else
	p4est_iterate(tree_nek, ghost_nek,(void *) &transfer_data,
			iter_msh_dat, NULL, NULL);
#endif
}

/* data type for refinement history data transfer between nek5000 and p4est*/
typedef struct transfer_hst_s {
	int map_nr; /**< local number of unchanged elements */
	int rfn_nr; /**< local number of refined elements */
	int crs_nr; /**< local number of coarsened elements */
    int *glgl_map; /**< element number mapping for unchanged elements */
	int *glgl_rfn; /**< element number mapping for refined elements */
	int *glgl_crs; /**< element number mapping for coarsened elements */
} transfer_hst_t;

#define MIN( a, b ) ( ( a > b) ? b : a )

/** @brief Iterate over element volumes to transfer refinement history data
 *
 * @details Required by fp4est_msh_get_hst
 *
 * @param info
 * @param user_data
 */
void iter_msh_hst(p4est_iter_volume_info_t * info, void *user_data) {
	user_data_t *data = (user_data_t *) info->quad->p.user_data;
	transfer_hst_t *trans_data = (transfer_hst_t *) user_data;

	// which quad (local and global element number)
	p4est_tree_t *tree;
	p4est_locidx_t iwl;
	int iwlt, iwg;
	int ifc, ifl, ic, il;// loop index

	// get quad number
	tree = p4est_tree_array_index(info->p4est->trees, info->treeid);
	// local quad number
	iwl = info->quadid + tree->quadrants_offset;
	iwlt = (int) iwl;
	// global quad number
	iwg = (int) tree_nek->global_first_quadrant[tree_nek->mpirank] + iwlt;

	// check refinement status
    ic = data->gln_children[0];
    for (il=1; il < P4EST_CHILDREN; il++) ic = MIN(ic,data->gln_children[il]);
    if (data->gln_parent == -1 && ic == -1) {
    	// no refinement
    	// count elements
    	trans_data->map_nr = trans_data->map_nr + 1;
    	// set element map
    	trans_data->glgl_map[iwlt] = data->gln_el + 1;
    } else if (data->gln_parent /= -1) {
    	// refinement
    	// count elements
    	trans_data->rfn_nr = trans_data->rfn_nr + 1;
    	// set dummy element map
    	trans_data->glgl_map[iwlt] = 0;
    	ic = (trans_data->rfn_nr-1)*3;
    	// current global element number
    	trans_data->glgl_rfn[ic] = iwg + 1;
    	// old parent element number
    	trans_data->glgl_rfn[ic +1] = data->gln_parent + 1;
    	// child position; numbered 0,..,P4EST_CHILDREN-1
    	trans_data->glgl_rfn[ic +2] = data->gln_el;
    } else {
    	// coarsening
    	// count elements
    	trans_data->crs_nr = trans_data->crs_nr + 1;
    	// set dummy element map
    	trans_data->glgl_map[iwlt] = 0;
    	// new global position
    	ic =(trans_data->crs_nr - 1)*2*P4EST_CHILDREN;
    	trans_data->glgl_crs[ic] = iwg + 1;
    	// old global position
    	trans_data->glgl_crs[ic+1] = data->gln_children[0] + 1;
    	for (il = 1; il < P4EST_CHILDREN; il++) {
    		// new dummy global position
    		trans_data->glgl_crs[ic+il*2] = 0;
    		// old global position
    		trans_data->glgl_crs[ic+il*2+1] = data->gln_children[il] + 1;
    	}
    }

	//	fill current global numbering of elements in nek5000 and reset refinement history
	data->gln_el = iwg;
	data->gln_parent = -1;
	for (ifc = 0; ifc < P4EST_CHILDREN; ++ifc) data->gln_children[ifc] = -1;
}

// get refinement history data to Nek5000
void fp4est_msh_get_hst(int * map_nr, int * rfn_nr, int * crs_nr, int *glgl_map,
		int * glgl_rfn, int * glgl_crs) {
	transfer_hst_t transfer_data;
	transfer_data.map_nr = 0;
	transfer_data.rfn_nr = 0;
	transfer_data.crs_nr = 0;
	transfer_data.glgl_map = glgl_map;
	transfer_data.glgl_rfn = glgl_rfn;
	transfer_data.glgl_crs = glgl_crs;
#ifdef P4_TO_P8
	p4est_iterate(tree_nek, ghost_nek,(void *) &transfer_data, iter_msh_hst,
			NULL, NULL, NULL);
#else
	p4est_iterate(tree_nek, ghost_nek,(void *) &transfer_data, iter_msh_hst,
			NULL, NULL);
#endif
	*map_nr = transfer_data.map_nr;
	*rfn_nr = transfer_data.rfn_nr;
	*crs_nr = transfer_data.crs_nr;
}

// get global vertex numbering
void fp4est_msh_get_node(int * lnelt, int * node) {
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
	if (nodes_nek) {
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
		p4est_locidx_t offset;

		// local number of quadrants
		*lnelt = (int) nodes_nek->num_local_quadrants;

		// quad to vertex global map
		for (il = 0; il < vi; ++il) {
			for (jl = 0; jl < P4EST_CHILDREN; ++jl) {

				nl = (int) nodes_nek->local_nodes[il * P4EST_CHILDREN + jl];
				if (nl < ni) {
					// independent node
					// do nothing
				} else if (nl < ni + fi) {
					// face hanging
					// find local vertex number on the face
					SC_ABORT("Hanging node; not supported yet\n");
				}
#ifdef P4_TO_P8
				else if (nl < ni + fi + ei) {
					// edge hanging
					// find local vertex number on the edge
					SC_ABORT("Hangign node; not supported yet\n");
				}
#endif
				else {
					SC_ABORT(
							"Wrong node number, aborting: fp4est_msh_get_node\n");
				}
				// get global number of this node
				indep = (p4est_indep_t *) sc_array_index(
						&nodes_nek->indep_nodes, nl);
				if (nl < offi) {
					nr = nodes_nek->nonlocal_ranks[nl];
				} else {
					nr = rank;
				}
				offset = 0;
				for (kl = 0; kl < nr; ++kl) {
					offset = offset + nodes_nek->global_owned_indeps[kl];
				}
				node[il * (P4EST_CHILDREN + 2) + 2 + jl] = (int) offset
						+ indep->p.piggy3.local_num + 1;

			}
		}
	} else {
		SC_ABORT("nodes_nek not allocated; aborting: fp4est_msh_get_node\n");
	}
}

/* get GLL node numbering */
void fp4est_msh_get_lnode(int * lnelt, int * lnoden, int * gnoden,
		p4est_gloidx_t * lnodes, int * hang_elm, int * hang_fsc, int * hang_edg) {
	int il, jl, kl;
	int hanging_face[P4EST_FACES];
#ifdef P4_TO_P8
	int hanging_edge[P8EST_EDGES];
#endif
	if (lnodes_nek) {
		const int vnd = (int) lnodes_nek->vnodes;
		const int owned = (int) lnodes_nek->owned_count;
		const int local = (int) lnodes_nek->num_local_nodes;
		const p4est_gloidx_t offset = lnodes_nek->global_offset;
		*lnelt = (int) lnodes_nek->num_local_elements;
		*lnoden = local;
		*gnoden = owned;
		for (il = 0; il < lnodes_nek->num_local_elements; ++il) {
			for (jl = 0; jl < P4EST_FACES; ++jl) {
				hanging_face[jl] = -1;
			}
#ifdef P4_TO_P8
			for (jl = 0; jl < P8EST_EDGES; ++jl) {
				hanging_edge[jl] = -1;
			}
			hang_elm[il] = p4est_lnodes_decode(lnodes_nek->face_code[il],
					hanging_face, hanging_edge);
			for (jl = 0; jl < P4EST_FACES; ++jl) {
				hang_fsc[il * P4EST_FACES + jl] = hanging_face[jl];
			}
			for (jl = 0; jl < P8EST_EDGES; ++jl) {
				hang_edg[il * P8EST_EDGES + jl] = hanging_edge[jl];
			}
#else
			hang_elm[il] = p4est_lnodes_decode(lnodes_nek->face_code[il],hanging_face);
			for (jl=0;jl<P4EST_FACES;++jl) {
				hang_fsc[il*P4EST_FACES +jl] = hanging_face[jl];
			}
#endif
			for (jl = 0; jl < vnd; ++jl) {
				kl = lnodes_nek->element_nodes[il * vnd + jl];
				if (kl < owned) {
					lnodes[il * vnd + jl] = (p4est_gloidx_t) 1 + offset + kl;
				} else if (kl < local) {
					kl =kl - owned;
					lnodes[il * vnd + jl] = (p4est_gloidx_t) 1 +
							lnodes_nek->nonlocal_nodes[kl];
				} else {
					SC_ABORT("Wrong node number; aborting: fp4est_msh_get_lnode\n");
				}
			}
		}
	} else {
		SC_ABORT("lnodes_nek not allocated; aborting: fp4est_msh_get_lnode\n");
	}

}

/** @brief Iterate over faces to get alignment
 *
 * @details Required by fp4est_msh_get_algn
 *
 * @param info
 * @param user_data
 */
void algn_fcs_get(p4est_iter_face_info_t * info, void *user_data) {
	int mpirank = info->p4est->mpirank;
	p4est_gloidx_t gfirst_quad =  info->p4est->global_first_quadrant[mpirank];
	int orient = (int) info->orientation;
	sc_array_t *sides = &(info->sides);
	p4est_iter_face_side_t *side;
	int nside = (int) sides->elem_count;
	int *fcs_arr = (int *) user_data;
	int il, jl;
	int iref, pref, pset, ipos;
	int8_t face[2], ftmp;
	if (info->tree_boundary&&orient){
		/* face is on the outside of the tree; compare orientation of different trees*/
		/* find the reference side; lowest face number */
		P4EST_ASSERT (nside <= 2);
		for (il = 0; il < nside; ++il) {
			side = p4est_iter_fside_array_index_int (sides, il);
			face[il] = side->face;
		}
		if (nside == 1) {
			/* single side no alignment */
			iref = nside +1;
		} else {
			/* 2 sides; find permutation set */
			if (face[0]<face[1]) {
				iref = 0;
#ifdef P4_TO_P8
				pref = p8est_face_permutation_refs[face[0]][face[1]];
				pset = p8est_face_permutation_sets[pref][orient];
#else
				pset = orient;
#endif
			} else {
				iref = 1;
#ifdef P4_TO_P8
				pref = p8est_face_permutation_refs[face[1]][face[0]];
				pset = p8est_face_permutation_sets[pref][orient];
#else
				pset = orient;
#endif
			}
		}

	} else {
		/* face is on the interior of the tree or orientation == 0;
		 * all quads aligned */
		iref = nside +1;
	}
	for (il = 0; il < nside; ++il) {
		side = p4est_iter_fside_array_index_int (sides, il);
		if (il == iref) {
			orient = pset;
		} else {
			orient = 0;
		}
		if (side->is_hanging) {
			/* hanging face */
			for (jl = 0; jl < P4EST_HALF; jl++) {
				if (!side->is.hanging.is_ghost[jl]){
					// local node
					ipos = (int) (side->treeid +
							side->is.hanging.quadid[jl] - gfirst_quad);
					ipos = ipos*P4EST_FACES + (int) side->face;
					fcs_arr[ipos] = orient;
				}
			}
		} else {
			if (!side->is.full.is_ghost) {
				// local node
				ipos = (int) (side->treeid + side->is.full.quadid -
						gfirst_quad);
				ipos = ipos*P4EST_FACES + (int) side->face;
				fcs_arr[ipos] = orient;
			}
		}
	}
}

#ifdef P4_TO_P8
/** @brief Iterate over edges to get alignment
 *
 * @details Required by fp4est_msh_get_algn
 *
 * @param info
 * @param user_data
 */
void algn_edg_get(p8est_iter_edge_info_t * info, void *user_data) {
	int mpirank = info->p4est->mpirank;
	p4est_gloidx_t gfirst_quad =  info->p4est->global_first_quadrant[mpirank];
	sc_array_t *sides = &(info->sides);
	p8est_iter_edge_side_t *side;
	int nside = (int) sides->elem_count;
	int *edg_arr = (int *) user_data;
	int il, jl;
	int ipos, tmin;
	if (info->tree_boundary){
		/* face is on the outside of the tree; compare orientation of
		 * different trees*/
		for (il = 0; il < nside; ++il) {
			side = p8est_iter_eside_array_index_int (sides, il);
			if (side->is_hanging) {
				/* hanging face */
				for (jl = 0; jl < 2; jl++) {
					if (!side->is.hanging.is_ghost[jl]){
						// local node
						ipos = (int) (side->treeid + side->is.hanging.quadid[jl] -
								gfirst_quad);
						ipos = ipos*P8EST_EDGES + (int) side->edge;
						edg_arr[ipos] = (int) side->orientation;
					}
				}
			} else {
				if (!side->is.full.is_ghost) {
					// local node
					ipos = (int) (side->treeid + side->is.full.quadid -
							gfirst_quad);
					ipos = ipos*P8EST_EDGES + (int) side->edge;
					edg_arr[ipos] = (int) side->orientation;
				}
			}
		}
	} else {
		/* edge is on the interior of the tree; all quads are aligned */
		for (il = 0; il < nside; ++il) {
			side = p8est_iter_eside_array_index_int (sides, il);
			if (side->is_hanging) {
				/* hanging face */
				for (jl = 0; jl < 2; jl++) {
					if (!side->is.hanging.is_ghost[jl]){
						// local node
						ipos = (int) (side->treeid +
								side->is.hanging.quadid[jl] - gfirst_quad);
						ipos = ipos*P8EST_EDGES + (int) side->edge;
						edg_arr[ipos] = 0;
					}
				}
			} else {
				if (!side->is.full.is_ghost) {
					// local node
					ipos = (int) (side->treeid + side->is.full.quadid -
							gfirst_quad);
					ipos = ipos*P8EST_EDGES + (int) side->edge;
					edg_arr[ipos] = 0;
				}
			}
		}
	}
}
#endif

/* get face and edge alignment */
void fp4est_msh_get_algn(int * fcs_algn, int *const edg_algn,
		int * lnelt) {
	*lnelt = (int) tree_nek->local_num_quadrants;
#ifdef P4_TO_P8
	p4est_iterate(tree_nek, ghost_nek,(void *) fcs_algn, NULL,
			algn_fcs_get, NULL, NULL);
	p4est_iterate(tree_nek, ghost_nek,(void *) edg_algn, NULL,
			NULL, algn_edg_get, NULL);
#else
	p4est_iterate(tree_nek, ghost_nek,(void *) fcs_algn, NULL,
			algn_fcs_get, NULL);
#endif
}

/* get graph for partitioning */
void fp4est_msh_get_graph(int * node_num, int * graph, int * graph_offset)
{
	// loop indexes
	int il,jl,kl,hl;
	//
	int offset,element,el_count,node, itmp;
	// number of processors, proc. id., local number of elements
	const int num_procs = tree_nek->mpisize;
	const int rank = tree_nek->mpirank;
	const int el_local = (int) mesh_nek->local_num_quadrants;
	// variables to deal with nonuniform connections
	sc_array_t  *halfs;
	p4est_locidx_t *halfentries;
	sc_array_t  *ghosts = &(ghost_nek->ghosts);
	p4est_quadrant_t *lquad;
	// halfs position
	halfs = mesh_nek->quad_to_half;
	// create the graph
	// initial offset end element count
	offset=0;
	el_count = 0;
	graph_offset[0] = offset;
	// in case no elements is counted
	graph_offset[1] = offset;
	// loop over elements
	for(il=0;il<el_local;++il){
		// to count elements
		itmp = offset;
		// loop over faces
		for(jl=0;jl<P4EST_FACES;++jl){
			kl=il*P4EST_FACES +jl;
			// half-size neighbour; multiple entries per face
			if (mesh_nek->quad_to_face[kl]<0){
				halfentries = (p4est_locidx_t *)
						sc_array_index (halfs,mesh_nek->quad_to_quad[kl]);
				for(hl=0;hl<P4EST_HALF;++hl){
					element = (int) halfentries[hl];
					// is the element non local
					if (element>=el_local){
						element = element - el_local;
						lquad = (p4est_quadrant_t *)
								sc_array_index (ghosts,element);
						node = (int) mesh_nek->ghost_to_proc[element];
						element = (int) lquad->p.piggy3.local_num;
						element =  element +
								(int) tree_nek->global_first_quadrant[node];
					}
					else{
						element =  tree_nek->global_first_quadrant[rank] +
								element;
					}
					graph[offset] = element;
					offset = offset+1;
				}
			}
			// single entrance per face
			else{
				element = (int) mesh_nek->quad_to_quad[kl];
				// do I point myself
				if (element!=il){
					// is the element non local
					if (element>=el_local){
						element = element - el_local;
						lquad = (p4est_quadrant_t *) sc_array_index (ghosts,element);
						node = (int) mesh_nek->ghost_to_proc[element];
						element = (int) lquad->p.piggy3.local_num;
						element =  element + (int) tree_nek->global_first_quadrant[node];
					}
					else{
						element =  tree_nek->global_first_quadrant[rank] +  element;
					}
					graph[offset] = element;
					offset = offset+1;
				}
			}
		}
		// set graph offset
		if (itmp<offset){
			el_count = el_count +1;
			graph_offset[el_count] = offset;
		}
	}
	*node_num = el_count;
}

/** @brief Iterate over element volumes to transfer element refinement mark
 *
 * @details Required by fp4est_refm_put
 *
 * @param info
 * @param user_data
 */
void iter_refm(p4est_iter_volume_info_t * info, void *user_data) {
	user_data_t *data = (user_data_t *) info->quad->p.user_data;
	int *ref_mark = (int *) user_data;

	// which quad (local and global element number)
	p4est_tree_t *tree;
	p4est_locidx_t iwl;
	int iwlt, iwg;
	int il;// loop index

	// get quad number
	tree = p4est_tree_array_index(info->p4est->trees, info->treeid);
	// local quad number
	iwl = info->quadid + tree->quadrants_offset;
	iwlt = (int) iwl;
	// global quad number
	iwg = (int) tree_nek->global_first_quadrant[tree_nek->mpirank] + iwlt;

    // correct global element number on nek5000 side if necessary; it can happen
    // if .rea file was used for mesh generation
    if(data->gln_el == -1){
            data->gln_el = iwg;
            data->gln_parent = -1;
            for(il=0;il<P4EST_CHILDREN;il++){
                    data->gln_children[il] = -1;
            }
    }

    // refinement mark for given quad
    data->ref_mark = ref_mark[iwlt];
}

/* fill ref_mark in p4est block */
void fp4est_refm_put(int * ref_mark) {
#ifdef P4_TO_P8
	p4est_iterate(tree_nek, ghost_nek,(void *) ref_mark, iter_refm,
			NULL, NULL, NULL);
#else
	p4est_iterate(tree_nek, ghost_nek,(void *) ref_mark, iter_refm,
			NULL, NULL);
#endif
}

/** @brief Iterate over faces to correct Nek5000 boundary conditions
 *
 * @details Required by fp4est_bc_check
 *
 * @param info
 * @param user_data
 */
void iter_bc_chk(p4est_iter_face_info_t * info, void *user_data) {
	int *fldn = (int *) user_data;
	user_data_t *data;

	int mpirank = info->p4est->mpirank;
	p4est_gloidx_t gfirst_quad =  info->p4est->global_first_quadrant[mpirank];

	sc_array_t *sides = &(info->sides);
	p4est_iter_face_side_t *side;
	int nside = (int) sides->elem_count;


	int il, jl, kl, ifl; // loop index
	int iwl; // global quad number
	int face; // quad face
	int imsh[2]; // mesh type mark

	if (info->tree_boundary) {
		/* Face is on the outside of the tree; it can be any type of boundary
		 * condition. V-type type mesh can have external bc inside the mesh if
		 * a neighbor is T-type quad.
		 */
		P4EST_ASSERT (nside <= 2);
		if (nside == 1) {
			/* external face; no E, P bc */
			side = p4est_iter_fside_array_index_int (sides, 0);
			face = (int) side->face;
			if (side->is_hanging) {
				/* hanging face */
				for (jl = 0; jl < P4EST_HALF; jl++) {
					if (!side->is.hanging.is_ghost[jl]) {
						/* local node */
						data = (user_data_t *) side->is.hanging.quad[jl]->p.user_data;
						/* loop over fluids */
						for (ifl = fldn[0]; ifl <= fldn[1]; ifl++) {
							/* skip velocity bc for T-type mesh */
							if (data->imsh && ifl == NP4_VFLD) continue;
							/* check connection type (should not be E or P) */
							if (data->cbc[ifl][face][2] == ' '
									&& data->cbc[ifl][face][1] == ' '
											&& (data->cbc[ifl][face][0] == 'E'
													|| data->cbc[ifl][face][0] == 'P')) {
								iwl = (int) (gfirst_quad + side->treeid + side->is.hanging.quadid[jl]);
								printf("Connectivity error: %i %i %i %i %c%c%c\n",
										tree_nek->mpirank, iwl, face, ifl,
										data->cbc[ifl][face][0], data->cbc[ifl][face][1],
										data->cbc[ifl][face][2]);
								printf("external face marked as internal.\n");
								SC_ABORT("Aborting: iter_bc_chk\n");
							}
						} // fld loop
					} // is ghost
				} // children loop

			} else {
				if (!side->is.full.is_ghost) {
					/* local node */
					data = (user_data_t *) side->is.full.quad->p.user_data;
					/* loop over fluids */
					for (ifl = fldn[0]; ifl <= fldn[1]; ifl++) {
						/* skip velocity bc for T-type mesh */
						if (data->imsh && ifl == NP4_VFLD) continue;
						/* check connection type (should not be E or P) */
						if (data->cbc[ifl][face][2] == ' '
								&& data->cbc[ifl][face][1] == ' '
										&& (data->cbc[ifl][face][0] == 'E'
												|| data->cbc[ifl][face][0] == 'P')) {
							iwl = (int) (gfirst_quad + side->treeid + side->is.full.quadid);
							printf("Connectivity error: %i %i %i %i %c%c%c\n",
									tree_nek->mpirank, iwl, face, ifl,
									data->cbc[ifl][face][0], data->cbc[ifl][face][1],
									data->cbc[ifl][face][2]);
							printf("external face marked as internal.\n");
							SC_ABORT("Aborting: iter_bc_chk\n");
						}
					} // fld loop
				} // is ghost
			} // hanging

		} else {
			/* internal face; any face type possible */
			/* collect quad type; for different values of imsh V-type
			 * elements should point to external bc
			 */
			imsh[0] = 0;
			imsh[1] = 0;
			for (il = 0; il < nside; ++il) {
				side = p4est_iter_fside_array_index_int (sides, il);
				if (side->is_hanging) {
					/* imsh for all children is the same */
					if (side->is.hanging.is_ghost[0]) {
						data = (user_data_t *) &ghost_data[side->is.hanging.quadid[0]];
					} else {
						data = (user_data_t *) side->is.hanging.quad[0]->p.user_data;
					}
				} else {
					if (side->is.full.is_ghost) {
						data = (user_data_t *) &ghost_data[side->is.full.quadid];
					} else {
						data = (user_data_t *) side->is.full.quad->p.user_data;
					}
				}
				imsh[il] = data->imsh;
			}
			/* test bc */
			for (il = 0; il < nside; ++il) {
				side = p4est_iter_fside_array_index_int (sides, il);
				face = (int) side->face;
				if (side->is_hanging) {
					/* hanging face */
					for (jl = 0; jl < P4EST_HALF; jl++) {
						if (!side->is.hanging.is_ghost[jl]) {
							/* local node */
							data = (user_data_t *) side->is.hanging.quad[jl]->p.user_data;
							/* loop over fluids */
							for (ifl = fldn[0]; ifl <= fldn[1]; ifl++) {
								/* skip velocity bc for T-type mesh */
								if (data->imsh && ifl == NP4_VFLD) continue;
								/* velocity bc for V-type element neighbor to T-type element
								 * should not be E or P
								 */
								if (imsh[0] != imsh[1] && ifl == NP4_VFLD) {
									if (data->cbc[ifl][face][2] == ' '
											&& data->cbc[ifl][face][1] == ' '
													&& (data->cbc[ifl][face][0] == 'E'
															|| data->cbc[ifl][face][0] == 'P')) {
										iwl = (int) (gfirst_quad + side->treeid + side->is.hanging.quadid[jl]);
										printf("Connectivity error: %i %i %i %i %c%c%c\n",
												tree_nek->mpirank, iwl, face, ifl,
												data->cbc[ifl][face][0], data->cbc[ifl][face][1],
												data->cbc[ifl][face][2]);
										printf("velocity external face marked as internal.\n");
										SC_ABORT("Aborting: iter_bc_chk\n");
									}
								} else {
									/* not velocity or not V-T meshes boundary - all
									 * internal elements
									 */
									if (data->cbc[ifl][face][2] == ' '
											&& data->cbc[ifl][face][1] == ' '
													&& (data->cbc[ifl][face][0] == 'E'
															|| data->cbc[ifl][face][0] == 'P')) {
										/* For internal faces I do not store neither global quad number
										 * nor face number of the neighbor, as Nek5000 data structure
										 * does not allow to keep all the information for hanging nodes
										 */
										for (kl=0; kl < 5; kl++) data->bc[ifl][face][kl] = 0.0;
									} else {
										iwl = (int) (gfirst_quad + side->treeid + side->is.hanging.quadid[jl]);
										printf("Connectivity error: %i %i %i %i %c%c%c\n",
												tree_nek->mpirank, iwl, face, ifl,
												data->cbc[ifl][face][0], data->cbc[ifl][face][1],
												data->cbc[ifl][face][2]);
										printf("internal face marked as external.\n");
										SC_ABORT("Aborting: iter_bc_chk\n");
									}
								}
							} // fld loop
						} // is ghost
					} // children loop

				} else {
					if (!side->is.full.is_ghost) {
						/* local node */
						data = (user_data_t *) side->is.full.quad->p.user_data;
						/* loop over fluids */
						for (ifl = fldn[0]; ifl <= fldn[1]; ifl++) {
							/* skip velocity bc for T-type mesh */
							if (data->imsh && ifl == NP4_VFLD) continue;
							/* velocity bc for V-type element neighbor to T-type element
							 * should not be E or P
							 */
							if (imsh[0] != imsh[1] && ifl == NP4_VFLD) {
								if (data->cbc[ifl][face][2] == ' '
										&& data->cbc[ifl][face][1] == ' '
												&& (data->cbc[ifl][face][0] == 'E'
														|| data->cbc[ifl][face][0] == 'P')) {
									iwl = (int) (gfirst_quad + side->treeid + side->is.full.quadid);
									printf("Connectivity error: %i %i %i %i %c%c%c\n",
											tree_nek->mpirank, iwl, face, ifl,
											data->cbc[ifl][face][0], data->cbc[ifl][face][1],
											data->cbc[ifl][face][2]);
									printf("velocity external face marked as internal.\n");
									SC_ABORT("Aborting: iter_bc_chk\n");
								}
							} else {
								/* not velocity or not V-T meshes boundary - all
								 * internal elements
								 */
								if (data->cbc[ifl][face][2] == ' '
										&& data->cbc[ifl][face][1] == ' '
												&& (data->cbc[ifl][face][0] == 'E'
														|| data->cbc[ifl][face][0] == 'P')) {
									/* For internal faces I do not store neither global quad number
									 * nor face number of the neighbor, as Nek5000 data structure
									 * does not allow to keep all the information for hanging nodes
									 */
									for (kl=0; kl < 5; kl++) data->bc[ifl][face][kl] = 0.0;
								} else {
									iwl = (int) (gfirst_quad + side->treeid + side->is.full.quadid);
									printf("Connectivity error: %i %i %i %i %c%c%c\n",
											tree_nek->mpirank, iwl, face, ifl,
											data->cbc[ifl][face][0], data->cbc[ifl][face][1],
											data->cbc[ifl][face][2]);
									printf("internal face marked as external.\n");
									SC_ABORT("Aborting: iter_bc_chk\n");
								}
							}
						} // fld loop
					} // is ghost
				} // hanging
			} // side loop
		} // forest internal/external face
	} else {
		/* face is on the interior of the tree; all faces are E (or J in Nek5000
		 * notation; I do not use J as this could cause problems with cyclic
		 * bc; see Nek500 routines check_cyclic, rotate_cyc)
		 *  as even
		 * for V-type mesh external and periodic boundary can be defined on
		 * tree faces only
		 */
		P4EST_ASSERT (nside == 2);
		for (il = 0; il < nside; ++il) {
			side = p4est_iter_fside_array_index_int (sides, il);
			face = (int) side->face;
			if (side->is_hanging) {
				/* hanging face */
				for (jl = 0; jl < P4EST_HALF; jl++) {
					if (!side->is.hanging.is_ghost[jl]) {
						// local node
						data = (user_data_t *) side->is.hanging.quad[jl]->p.user_data;
						// loop over fluids
						for (ifl = fldn[0]; ifl <= fldn[1]; ifl++) {
							/* skip velocity bc for T-type mesh */
							if (data->imsh && ifl == NP4_VFLD) continue;
							data->cbc[ifl][face][0] = 'E';
							data->cbc[ifl][face][1] = ' ';
							data->cbc[ifl][face][2] = ' ';
							/* For internal faces I do not store neither global quad number
							 * nor face number of the neighbor, as Nek5000 data structure
							 * does not allow to keep all the information for hanging nodes
							 */
							for (kl=0; kl < 5; kl++) data->bc[ifl][face][kl] = 0.0;
						} // fld loop
					} // is ghost
				} // children loop

			} else {
				if (!side->is.full.is_ghost) {
					/* local node */
					data = (user_data_t *) side->is.full.quad->p.user_data;
					/* loop over fluids */
					for (ifl = fldn[0]; ifl <= fldn[1]; ifl++) {
						/* skip velocity bc for T-type mesh */
						if (data->imsh && ifl == NP4_VFLD) continue;
						data->cbc[ifl][face][0] = 'E';
						data->cbc[ifl][face][1] = ' ';
						data->cbc[ifl][face][2] = ' ';
						/* For internal faces I do not store neither global quad number
						 * nor face number of the neighbor, as Nek5000 data structure
						 * does not allow to keep all the information for hanging nodes
						 */
						for (kl=0; kl < 5; kl++) data->bc[ifl][face][kl] = 0.0;
					} // fld loop
				} // is ghost
			} // hanging
		} // side loop
	} // tree internal/external

}

/* Check boundary conditions for V- and T-type mesh */
void fp4est_bc_check(int * ibc, int * ebc) {
	int fldn[2];
	fldn[0] = *ibc;
	fldn[1] = *ebc;
#ifdef P4_TO_P8
	p4est_iterate(tree_nek, ghost_nek, (void *) &fldn, NULL, iter_bc_chk, NULL, NULL);
#else
	p4est_iterate(tree_nek, ghost_nek, (void *) &fldn, NULL, iter_bc_chk, NULL);
#endif
}
