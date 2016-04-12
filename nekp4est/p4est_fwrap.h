/** @file p4est_fwrap.h
 *  @brief Fortran77 interface of p4est library
 *
 * @author Adam Peplinski
 * @date Feb 20, 2016
 *
 * @ingroup nekp4est
 */
#ifndef NEKP4EST_P4EST_FWRAP_H_
#define NEKP4EST_P4EST_FWRAP_H_

#if !defined(NEKP4EST_H_)
#error "p4est_fwrap.h requires nekp4est.h"
#endif

#if !defined(P4EST_H) && !defined(P8EST_H)
#error "p4est_fwrap.h requires p4est.h or p8est.h"
#endif

#if !defined(NAME_H)
#error "p4est_fwrap.h requires name.h"
#endif

/* FORTRAN interface
 * #define fp4est_ FORTRAN_NAME(fp4est_,FP4EST_)
 */
/* wrappers */
/* initialize */
#define fp4est_init         FORTRAN_NAME(fp4est_init,FP4EST_INIT)
/* Connectivity */
#define fp4est_cnn_new      FORTRAN_NAME(fp4est_cnn_new,FP4EST_CNN_NEW)
#define fp4est_cnn_del      FORTRAN_NAME(fp4est_cnn_del,FP4EST_CNN_DEL)
#define fp4est_cnn_attr     FORTRAN_NAME(fp4est_cnn_attr,FP4EST_CNN_ATTR)
#define fp4est_cnn_valid    FORTRAN_NAME(fp4est_cnn_valid,FP4EST_CNN_VALID)
#define fp4est_cnn_compl    FORTRAN_NAME(fp4est_cnn_compl,FP4EST_CNN_COMPL)
#define fp4est_cnn_save     FORTRAN_NAME(fp4est_cnn_save,FP4EST_CNN_SAVE)
#define fp4est_cnn_load     FORTRAN_NAME(fp4est_cnn_load,FP4EST_CNN_LOAD)
/* tree_ management */
#define fp4est_tree_new     FORTRAN_NAME(fp4est_tree_new,FP4EST_TREE_NEW)
#define fp4est_tree_del     FORTRAN_NAME(fp4est_tree_del,FP4EST_TREE_DEL)
#define fp4est_tree_valid   FORTRAN_NAME(fp4est_tree_valid,FP4EST_TREE_VALID)
#define fp4est_tree_save    FORTRAN_NAME(fp4est_tree_save,FP4EST_TREE_SAVE)
#define fp4est_tree_load    FORTRAN_NAME(fp4est_tree_load,FP4EST_TREE_LOAD)
/* tree and grid info */
#define fp4est_ghost_new    FORTRAN_NAME(fp4est_ghost_new,FP4EST_GHOST_NEW)
#define fp4est_ghost_del    FORTRAN_NAME(fp4est_ghost_del,FP4EST_GHOST_DEL)
#define fp4est_mesh_new     FORTRAN_NAME(fp4est_mesh_new,FP4EST_MESH_NEW)
#define fp4est_mesh_del     FORTRAN_NAME(fp4est_mesh_del,FP4EST_MESH_DEL)
#define fp4est_nodes_new    FORTRAN_NAME(fp4est_nodes_new,FP4EST_NODES_NEW)
#define fp4est_nodes_del    FORTRAN_NAME(fp4est_nodes_del,FP4EST_NODES_DEL)
#define fp4est_lnodes_new   FORTRAN_NAME(fp4est_lnodes_new,FP4EST_LNODES_NEW)
#define fp4est_lnodes_del   FORTRAN_NAME(fp4est_lnodes_del,FP4EST_LNODES_DEL)
/* nekp4est internal load balance */
#define fp4est_part         FORTRAN_NAME(fp4est_part,FP4EST_PART)
/* refinement, coarsening, balance */
/* I/O (VTK) */
#define fp4est_vtk_write    FORTRAN_NAME(fp4est_vtk_write,FP4EST_VTK_WRITE)
#define fp4est_vtk_iscalar  FORTRAN_NAME(fp4est_vtk_iscalar,FP4EST_VTK_ISCALAR)
/* routines required by Nek5000 for data exchange */
#define fp4est_msh_get_size  FORTRAN_NAME(fp4est_msh_get_size,FP4EST_MSH_GET_SIZE)
#define fp4est_msh_get_dat   FORTRAN_NAME(fp4est_msh_get_dat,FP4EST_MSH_GET_DAT)
#define fp4est_msh_get_node   FORTRAN_NAME(fp4est_msh_get_node,FP4EST_MSH_GET_NODE)
/* callback routines */
#define nek_init_msh_dat FORTRAN_NAME(nek_init_msh_dat,NEK_INIT_MSH_DAT)
#define nek_get_msh_dat FORTRAN_NAME(nek_get_msh_dat,NEK_GET_MSH_DAT)
#define nek_get_msh_hst FORTRAN_NAME(nek_get_msh_hst,NEK_GET_MSH_HST)
#define nek_refine_mark FORTRAN_NAME(nek_refine_mark,NEK_REFINE_MARK)
/* Fortran common blocks required by p4est */
#define nekp4est_cbci FORTRAN_NAME(nekp4est_cbci,NEKP4EST_CBCI)

/** Data type for user variables; required by p4est */
typedef struct
{
	int imsh; /**< velocity (0) and temperature (1) mesh indicator */
	int igrp; /**< element group */
	int crv[6]; /**< curvature data set by face curvature group the face belongs to:
	0 - no curvature;
	<0 - curved internal face; should be described by GLL points; nothing to do
	>0 - curved external face; keeps the curvature group; correct on refinement */
	char cbc[N_PSCL+2][6][3]; /**< boundary condition data */
	double bc[N_PSCL+2][6][5];
	/* @{ */
	double x[P4EST_CHILDREN],y[P4EST_CHILDREN],z[P4EST_CHILDREN]; /**< vertices of the element */
	/* @} */
	int ref_mark; /**< integer to store refinement mark; 0 - nothing, 1 - refine, -1 coarsen */
	// to keep track of changes of nek5000 global element numbering
	int gln_el; /**< old element global numbering; nek5000 side */
	int gln_parent; /**< old parent global numbering; nek5000 side */
	int gln_children[P4EST_CHILDREN]; /**< old children global numbering; nek5000 side */
}
user_data_t;


/* Wrappers */
/** Initialize p4est package setting log verbosity
 *
 * @param [in] log_threshold  Log level
 */
void fp4est_init(int * log_threshold);

/** Initialize elements connectivity
 *
 * @param num_vertices
 * @param num_trees
 * @param num_edges
 * @param num_corners
 * @param vertices
 * @param tree_to_vertex
 * @param tree_to_tree
 * @param tree_to_face
 * @param tree_to_edge
 * @param ett_offset
 * @param edge_to_tree
 * @param edge_to_edge
 * @param tree_to_corner
 * @param ctt_offset
 * @param corner_to_tree
 * @param corner_to_corner
 */
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
		int * corner_to_tree, int * corner_to_corner);
#else
void fp4est_cnn_new(int * num_vertices, int * num_trees,
		int * num_corners,
		double *vertices,
		int * tree_to_vertex, int * tree_to_tree,
		int * tree_to_face,
		int * tree_to_corner, int * ctt_offset,
		int * corner_to_tree, int * corner_to_corner);
#endif

/** Destroy mesh connectivity */
void fp4est_cnn_del();

/** Allocate or free the attribute fields in a connectivity
 *
 * @param enable_tree_attr
 */
void fp4est_cnn_attr(int * enable_tree_attr);

/** Check connectivity consistency
 *
 * @param [out] is_valid   non zero for correct connectivity
 */
void fp4est_cnn_valid(int * is_valid);

/** Internally connect a connectivity based on tree_to_vertex information.
 * Only internal connectivity; nor mesh periodicity
 */
void fp4est_cnn_complete();

/** Save connectivity in a file
 *
 * @param filename
 * @param len_f
 */
void fp4est_cnn_save(char *filename, int len_f);

/** Load a connectivity structure from a file
 *
 * @param filename
 * @param len_f
 */
void fp4est_cnn_load(char * filename, int len_f);

/** Generate forest
 *
 * @param fmpicomm
 * @param min_level
 */
void fp4est_tree_new(MPI_Fint * fmpicomm, int * min_level);

/** Destroy tree */
void fp4est_tree_del();

/** Check tree consistency
 *
 * @param [out] is_valid     non zero for correct tree
 */
void fp4est_tree_valid(int * is_valid);

/** Save tree to the file
 *
 * @param [in] save_data      if non zero save user data
 * @param filename
 * @param len_f
 */
void fp4est_tree_save(int *save_data, char * filename, int len_f);

/** Load tree from a file
 *
 * @param fmpicomm
 * @param [in] load_data           if non zero read user data
 * @param filename
 * @param len_f
 */
void fp4est_tree_load(MPI_Fint * fmpicomm, int *load_data, char * filename, int len_f);

/** Build ghost layer */
void fp4est_ghost_new();

/** Destroy ghosr layer */
void fp4est_ghost_del();

/** Generate mesh information */
void fp4est_mesh_new();

/** Destroy mesh information */
void fp4est_mesh_del();

/** Generate new node information */
void fp4est_nodes_new();

/** Destroy node information */
void fp4est_nodes_del();

/** Generate new global nodes (GLL points) numbering
 *
 * @param [in] ldgr      polynomial degree +1
 */
void fp4est_lnodes_new(int ldgr);

/** Destroy global node numberring */
void fp4est_lnodes_del();

/** Forest partitioning for p4est
 *
 * @param partforcoarsen   partitioning strategy:
 * 0 - equal element count
 * 1 - octants families prepared for coarsening
 */
void fp4est_part(int * partforcoarsen);


/* place for refinement, coarsening and balancing */

/** Write tree structure to VTK file
 *
 * @param filename
 * @param len_f
 */
void fp4est_vtk_write(char * filename, int len_f);

/** Write integer scalar field to VTK file
 *
 * @param iscalar
 * @param num
 * @param filename
 * @param len_f
 */
void fp4est_vtk_iscalar(int * iscalar,int *num,char * filename, int len_f);

/** Get mesh size information to Nek5000
 *
 * @param [out] nelgt   global element number
 * @param [out] nelgit  element offset (number of elements on lower nid's)
 * @param [out] nelt    number ov T-type elements
 * @param [out] nelv    number ov V-type elements
 * @param [out] maxl    current max level
 */
void fp4est_msh_get_size(int * nelgt,int * nelgit, int * nelt, int * nelv, int * maxl);

/** Get mesh data to Nek5000
 *
 */
void fp4est_msh_get_dat();

/** Get global vertex numbering to Nek5000
 *
 * @param [out] lnelt   number of local elements
 * @param [out] node    array with global vertex numbering
 */
void fp4est_msh_get_node(int * lnelt, int * node);

#endif /* NEKP4EST_P4EST_FWRAP_H_ */
