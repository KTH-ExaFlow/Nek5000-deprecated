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

/*------------------------------------------------------------------------------
 * Global varaibles
------------------------------------------------------------------------------*/
static p4est_t            *p4est_nekton;
static p4est_connectivity_t *connectivity_nekton;
static p4est_ghost_t *ghost_nekton;
static p4est_nodes_t *nodes_nekton;
static p4est_lnodes_t	*lnodes_nekton;
static p4est_mesh_t *mesh_nekton;

/*------------------------------------------------------------------------------
 * internal subroutines for interaction with fortran
------------------------------------------------------------------------------*/
// initialise element data reading 0 level tree
// required by fp4est_tree_new but not really used
void
init_mshv (p4est_t * p4est, p4est_topidx_t which_tree,
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
	int i;

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

	for(i=0;i<*num_vertices*3;++i){
		verticesl[i] = (double) vertices[i];
	}

	for(i=0;i<*num_trees*P4EST_CHILDREN;++i){
		tree_to_vertexl[i] = (p4est_topidx_t) tree_to_vertex[i];
	}
	for(i=0;i<*num_trees*P4EST_FACES;++i){
		tree_to_treel[i] = (p4est_topidx_t) tree_to_tree[i];
		tree_to_facel[i] = (int8_t) tree_to_face[i];
	}

#ifdef P4_TO_P8
	for(i=0;i<*num_trees*P8EST_EDGES;++i){
		tree_to_edgel[i] = (p4est_topidx_t) tree_to_edge[i];
	}
	for(i=0;i<= *num_edges;++i){
		ett_offsetl[i] = (p4est_topidx_t) ett_offset[i];
	}
	for(i=0;i<ett_offset[*num_edges];++i){
		edge_to_treel[i] = (p4est_topidx_t) edge_to_tree[i];
		edge_to_edgel[i] = (int8_t) edge_to_edge[i];
	}
#endif

	for(i=0;i<*num_trees*P4EST_CHILDREN;++i){
		tree_to_cornerl[i] = (p4est_topidx_t) tree_to_corner[i];
	}
	for(i=0;i<= *num_corners;++i){
		ctt_offsetl[i] = (p4est_topidx_t) ctt_offset[i];
	}
	for(i=0;i<ctt_offset[*num_corners];++i){
		corner_to_treel[i] = (p4est_topidx_t) corner_to_tree[i];
		corner_to_cornerl[i] = (int8_t) corner_to_corner[i];
	}

#ifdef P4_TO_P8
	connectivity_nekton = p4est_connectivity_new_copy(num_verticesl, num_treesl,
			num_edgesl,
			num_cornersl,
			verticesl, tree_to_vertexl,
			tree_to_treel, tree_to_facel,
			tree_to_edgel, ett_offsetl,
			edge_to_treel, edge_to_edgel,
			tree_to_cornerl, ctt_offsetl,
			corner_to_treel, corner_to_cornerl);
#else
	connectivity_nekton = p4est_connectivity_new_copy(num_verticesl, num_treesl,
			num_cornersl,
			verticesl, tree_to_vertexl,
			tree_to_treel, tree_to_facel,
			tree_to_cornerl, ctt_offsetl,
			corner_to_treel, corner_to_cornerl);
#endif
}

void fp4est_cnn_del()
{
	p4est_connectivity_destroy(connectivity_nekton);
}

void fp4est_cnn_attr(int * enable_tree_attr)
{
	p4est_connectivity_set_attr (connectivity_nekton, *enable_tree_attr);
}

void fp4est_cnn_valid(int * is_valid)
{
	*is_valid = p4est_connectivity_is_valid(connectivity_nekton);
}

void fp4est_cnn_complete()
{
	p4est_connectivity_complete(connectivity_nekton);
}

void fp4est_cnn_save(char *filename, int len_f)
{
	p4est_connectivity_save (filename,connectivity_nekton);
}

void fp4est_cnn_load(char * filename, int len_f)
{
	connectivity_nekton = p4est_connectivity_load (filename, NULL);
}

/* tree_ management */
void fp4est_tree_new(MPI_Fint * fmpicomm, int * min_level)
{
	MPI_Comm mpicomm;
	mpicomm = MPI_Comm_f2c(*fmpicomm);
	p4est_nekton = p4est_new_ext (mpicomm,connectivity_nekton, 0, *min_level, 1,
			sizeof(user_data_t), init_mshv, NULL);
}

void fp4est_tree_del()
{
	p4est_destroy (p4est_nekton);
}

void fp4est_tree_valid(int * is_valid)
{
	*is_valid = p4est_is_valid(p4est_nekton);
}

void fp4est_tree_save(int *save_data, char * filename, int len_f)
{
	p4est_save (filename, p4est_nekton, *save_data);
}

void fp4est_tree_load(MPI_Fint * fmpicomm, int *load_data, char * filename, int len_f)
{
	MPI_Comm mpicomm;
	mpicomm = MPI_Comm_f2c(*fmpicomm);
	p4est_nekton = p4est_load (filename, mpicomm, sizeof(user_data_t), *load_data,
			NULL, &connectivity_nekton);
}

/* tree and grid info */
void fp4est_ghost_new()
{
	ghost_nekton = p4est_ghost_new(p4est_nekton,P4EST_CONNECT_FULL);
}

void fp4est_ghost_del()
{
	p4est_ghost_destroy(ghost_nekton);
}

void fp4est_mesh_new()
{
	mesh_nekton = p4est_mesh_new(p4est_nekton,ghost_nekton,P4EST_CONNECT_FULL);
}

void fp4est_mesh_del()
{
	p4est_mesh_destroy(mesh_nekton);
}

void fp4est_nodes_new()
{
	nodes_nekton = p4est_nodes_new(p4est_nekton,ghost_nekton);
}

void fp4est_nodes_del()
{
	p4est_nodes_destroy(nodes_nekton);
}

void fp4est_lnodes_new(int ldgr)
{
	lnodes_nekton = p4est_lnodes_new(p4est_nekton,ghost_nekton,ldgr);
}

void fp4est_lnodes_del()
{
	p4est_lnodes_destroy(lnodes_nekton);
}

/* nekp4est internal load balance */
void fp4est_part()
{
        /* Partition: The quadrants are redistributed for equal element count.  The
         * partition can optionally be modified such that a family of octants, which
         * are possibly ready for coarsening, are never split between processors. */
        int partforcoarsen= 0;
        p4est_partition (p4est_nekton,partforcoarsen,NULL);
}

/* place for refinement, coarsening, balance */

/* I/O (VTK) */
void fp4est_vtk_write(char * filename, int len_f)
{
	p4est_vtk_write_file (p4est_nekton, NULL, filename);
}

void fp4est_vtk_iscalar(int * iscalar,int *num,char * filename, int len_f)
{
	int i,j,err;
	double rscalar[*num*P4EST_CHILDREN];

	for(i=0;i<*num;++i){
		for(j=0;j<P4EST_CHILDREN;j++){
			rscalar[i*P4EST_CHILDREN+j] = (double) iscalar[i];
		}
	}

	err = p4est_vtk_write_header (p4est_nekton, NULL, 1.0, 0, 0, 0, 0, "scalar", NULL, filename);
	err = p4est_vtk_write_point_scalar (p4est_nekton, NULL, filename, "scalar", rscalar);
	err = p4est_vtk_write_footer (p4est_nekton, filename);
}
