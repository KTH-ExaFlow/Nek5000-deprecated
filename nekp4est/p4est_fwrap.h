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
#define fp4est_init          FORTRAN_NAME(fp4est_init,FP4EST_INIT)
/* Connectivity */
#define fp4est_cnn_new       FORTRAN_NAME(fp4est_cnn_new,FP4EST_CNN_NEW)
#define fp4est_cnn_del       FORTRAN_NAME(fp4est_cnn_del,FP4EST_CNN_DEL)
#define fp4est_cnn_attr      FORTRAN_NAME(fp4est_cnn_attr,FP4EST_CNN_ATTR)
#define fp4est_cnn_valid     FORTRAN_NAME(fp4est_cnn_valid,FP4EST_CNN_VALID)
#define fp4est_cnn_compl     FORTRAN_NAME(fp4est_cnn_compl,FP4EST_CNN_COMPL)
#define fp4est_cnn_save      FORTRAN_NAME(fp4est_cnn_save,FP4EST_CNN_SAVE)
#define fp4est_cnn_load      FORTRAN_NAME(fp4est_cnn_load,FP4EST_CNN_LOAD)
/* tree_ management */
#define fp4est_tree_new      FORTRAN_NAME(fp4est_tree_new,FP4EST_TREE_NEW)
#define fp4est_tree_del      FORTRAN_NAME(fp4est_tree_del,FP4EST_TREE_DEL)
#define fp4est_tree_valid    FORTRAN_NAME(fp4est_tree_valid,FP4EST_TREE_VALID)
#define fp4est_tree_save     FORTRAN_NAME(fp4est_tree_save,FP4EST_TREE_SAVE)
#define fp4est_tree_load     FORTRAN_NAME(fp4est_tree_load,FP4EST_TREE_LOAD)
/* tree and grid info */
#define fp4est_ghost_new     FORTRAN_NAME(fp4est_ghost_new,FP4EST_GHOST_NEW)
#define fp4est_ghost_del     FORTRAN_NAME(fp4est_ghost_del,FP4EST_GHOST_DEL)
#define fp4est_mesh_new      FORTRAN_NAME(fp4est_mesh_new,FP4EST_MESH_NEW)
#define fp4est_mesh_del      FORTRAN_NAME(fp4est_mesh_del,FP4EST_MESH_DEL)
#define fp4est_nodes_new     FORTRAN_NAME(fp4est_nodes_new,FP4EST_NODES_NEW)
#define fp4est_nodes_del     FORTRAN_NAME(fp4est_nodes_del,FP4EST_NODES_DEL)
#define fp4est_lnodes_new    FORTRAN_NAME(fp4est_lnodes_new,FP4EST_LNODES_NEW)
#define fp4est_lnodes_del    FORTRAN_NAME(fp4est_lnodes_del,FP4EST_LNODES_DEL)
/* nekp4est internal load balance */
#define fp4est_part          FORTRAN_NAME(fp4est_part,FP4EST_PART)
/* refinement, coarsening, balance */
#define fp4est_refine        FORTRAN_NAME(fp4est_refine,FP4EST_REFINE)
#define fp4est_coarsen       FORTRAN_NAME(fp4est_coarsen,FP4EST_COARSEN)
#define fp4est_balance       FORTRAN_NAME(fp4est_balance,FP4EST_BALANCE)
/* I/O (VTK) */
#define fp4est_vtk_write     FORTRAN_NAME(fp4est_vtk_write,FP4EST_VTK_WRITE)
#define fp4est_vtk_iscalar   FORTRAN_NAME(fp4est_vtk_iscalar,FP4EST_VTK_ISCALAR)
/* routines required by Nek5000 for data exchange */
#define fp4est_msh_get_size  FORTRAN_NAME(fp4est_msh_get_size,FP4EST_MSH_GET_SIZE)
#define fp4est_msh_get_dat   FORTRAN_NAME(fp4est_msh_get_dat,FP4EST_MSH_GET_DAT)
#define fp4est_msh_get_hst   FORTRAN_NAME(fp4est_msh_get_hst,FP4EST_MSH_GET_HST)
#define fp4est_msh_get_node  FORTRAN_NAME(fp4est_msh_get_node,FP4EST_MSH_GET_NODE)
#define fp4est_msh_get_lnode FORTRAN_NAME(fp4est_msh_get_lnode,FP4EST_MSH_GET_LNODE)
#define fp4est_msh_get_algn  FORTRAN_NAME(fp4est_msh_get_algn,FP4EST_MSH_GET_ALGN)
#define fp4est_msh_get_graph FORTRAN_NAME(fp4est_msh_get_graph,FP4EST_MSH_GET_GRAPH)
#define fp4est_refm_put      FORTRAN_NAME(fp4est_refm_put,FP4EST_REFM_PUT)
#define fp4est_bc_check      FORTRAN_NAME(fp4est_bc_check,FP4EST_BC_CHECK)

/** Data type for user variables; required by p4est */
typedef struct user_data_s {
	int imsh; /**< velocity (0) and temperature (1) mesh indicator */
	int igrp; /**< element group */
	int crv[6]; /**< curvature data set by face curvature group the face belongs to:
	 0 - no curvature;
	 <0 - curved internal face; should be described by GLL points; nothing to do
	 >0 - curved external face; keeps the curvature group; correct on refinement */
	char cbc[N_PSCL + 2][6][3]; /**< boundary condition data */
	double bc[N_PSCL + 2][6][5];
	/* @{ */
	double x[P4EST_CHILDREN], y[P4EST_CHILDREN],
	       z[P4EST_CHILDREN]; /**< vertices of the element */
	/* @} */
	int ref_mark; /**< integer to store refinement mark; definition in nekp4est.h */
	// to keep track of changes of nek5000 global element numbering
	int gln_el; /**< old element global numbering; nek5000 side */
	int gln_parent; /**< old parent global numbering; nek5000 side */
	int gln_children[P4EST_CHILDREN]; /**< old children global numbering; nek5000 side */
} user_data_t;

/* Wrappers */
/** Initialize p4est package setting log verbosity
 *
 * @param [in] log_threshold  Log level
 */
void fp4est_init(int * log_threshold)
;

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
void fp4est_cnn_del()
;

/** Allocate or free the attribute fields in a connectivity
 *
 * @param enable_tree_attr
 */
void fp4est_cnn_attr(int * enable_tree_attr)
;

/** Check connectivity consistency
 *
 * @param [out] is_valid   non zero for correct connectivity
 */
void fp4est_cnn_valid(int * is_valid)
;

/** Internally connect a connectivity based on tree_to_vertex information.
 * Only internal connectivity; nor mesh periodicity
 */
void fp4est_cnn_complete();

/** Save connectivity in a file
 *
 * @param filename
 * @param len_f
 */
void fp4est_cnn_save(char *filename, int len_f)
;

/** Load a connectivity structure from a file
 *
 * @param filename
 * @param len_f
 */
void fp4est_cnn_load(char * filename, int len_f)
;

/** Generate forest
 *
 * @param fmpicomm
 * @param min_level
 */
void fp4est_tree_new(MPI_Fint * fmpicomm, int * min_level)
;

/** Destroy tree */
void fp4est_tree_del()
;

/** Check tree consistency
 *
 * @param [out] is_valid     non zero for correct tree
 */
void fp4est_tree_valid(int * is_valid)
;

/** Save tree to the file
 *
 * @param [in] save_data      if non zero save user data
 * @param filename
 * @param len_f
 */
void fp4est_tree_save(int *save_data, char * filename, int len_f)
;

/** Load tree from a file
 *
 * @param fmpicomm
 * @param [in] load_data           if non zero read user data
 * @param filename
 * @param len_f
 */
void fp4est_tree_load(
		MPI_Fint * fmpicomm, int *load_data, char * filename, int len_f)
;

/** Build ghost layer */
void fp4est_ghost_new()
;

/** Destroy ghost layer */
void fp4est_ghost_del()
;

/** Generate mesh information */
void fp4est_mesh_new()
;

/** Destroy mesh information */
void fp4est_mesh_del()
;

/** Generate new node information */
void fp4est_nodes_new()
;

/** Destroy node information */
void fp4est_nodes_del()
;

/** Generate new global nodes (GLL points) numbering */
void fp4est_lnodes_new()
;

/** Destroy global node numberring */
void fp4est_lnodes_del()
;

/** Forest partitioning for p4est
 *
 * @param partforcoarsen   partitioning strategy:
 * 0 - equal element count
 * 1 - octants families prepared for coarsening
 */
void fp4est_part(int * partforcoarsen)
;

/** Perform tree refinement
 *
 * @param max_level    max refinement level
 */
void fp4est_refine(int *max_level)
;

/** Perform tree coarsening
 */
void fp4est_coarsen()
;

/** Perform 2:1 tree balancing
 */
void fp4est_balance()
;

/** Write tree structure to VTK file
 *
 * @param filename
 * @param len_f
 */
void fp4est_vtk_write(char * filename, int len_f)
;

/** Write integer scalar field to VTK file
 *
 * @param iscalar
 * @param num
 * @param filename
 * @param len_f
 */
void fp4est_vtk_iscalar(int * iscalar,int *num,char * filename, int len_f)
;

/** Get mesh size information to Nek5000
 *
 * @param [out] nelgt   global element number
 * @param [out] nelgit  element offset (number of elements on lower nid's)
 * @param [out] nelt    number of T-type elements
 * @param [out] nelv    number of V-type elements
 * @param [out] maxl    current max level
 */
void fp4est_msh_get_size(
		int * nelgt,int * nelgit, int * nelt, int * nelv, int * maxl)
;

/** Get mesh data to Nek5000
 *
 * @param[in]  ibc   starting position in cb and cbc arrays
 * @param[in]  ebc   ending position in cb and cbc arrays
 * @param[in]  nelv  global number of V-type elements
 * @param[in]  lelt  array dimension for bc, cbc
 * @param[out] igrp  element group
 * @param[out] level element level
 * @param[out] crv   element curvature
 * @param[out] bc    element boundary condition parameters
 * @param[out] cbc   element boundary condition type
 */
void fp4est_msh_get_dat(int * ibc, int * nfld, int * nelv, int *lelt, int * igrp,
		int * level, int * crv, double * bc, char * cbc)
;

/** Get refinement history data to Nek5000
 *
 * @param map_nr     local number of unchanged elements
 * @param rfn_nr     local number of refined elements
 * @param crs_nr     local number of coarsened elements
 * @param glgl_map   element number mapping for unchanged elements
 * @param glgl_rfn   element number mapping for refined elements
 * @param glgl_crs   element number mapping for coarsened elements
 */
void fp4est_msh_get_hst(int * map_nr, int * rfn_nr, int * crs_nr, int *glgl_map,
		int * glgl_rfn, int * glgl_crs)
;

/** Get global vertex numbering to Nek5000
 *
 * @param [out] lnelt   number of local elements
 * @param [out] node    array with global vertex numbering
 */
void fp4est_msh_get_node(int * lnelt, int * node)
;

/** Get global vertex, face and edge numbering to Nek5000
 *
 * @param lnelt      local number of elements
 * @param lnoden     local number of nodes (vert., fac., edg.)
 * @param gnoden     local number of owned nodes
 * @param lnodes     global node numbering list
 * @param hang_elm   is any face/edge hanging list
 * @param hang_fsc   hanging face list
 * @param hang_edg   hanging edge list
 */
void fp4est_msh_get_lnode(
		int * lnelt, int * lnoden, int * gnoden,p4est_gloidx_t * lnodes,
		int * hang_elm, int * hang_fsc, int * hang_edg)
;

/** Get face and edge orientation
 *
 * @param fcs_algn   face alignment
 * @param edg_algn   edge alignment
 * @param lnelt      local number of elements
 */
void fp4est_msh_get_algn(int * fcs_algn, int *const edg_algn, int * lnelt)
;

/** Get adjacency graph for partitioning
 *
 * @param node_num      node number in the graph (to create vtxdist in ParMetis notation)
 * @param graph         graph - (adjncy in ParMetis notation)
 * @param graph_offset  graph_offset - (xadj in ParMetis notation)
 */
void fp4est_msh_get_graph(int * node_num, int * graph, int * graph_offset)
;

/** Fill ref_mark in p4est block
 *
 * @param ref_mark   refinement mark array
 */
void fp4est_refm_put(int * ref_mark)
;

/** Check boundary conditions for V- and T-type mesh
 *
 * @param ibc    starting position in cb and cbc arrays
 * @param ebc    ending position in cb and cbc arrays
 */
void fp4est_bc_check(int * ibc, int * ebc)
;

#endif /* NEKP4EST_P4EST_FWRAP_H_ */
