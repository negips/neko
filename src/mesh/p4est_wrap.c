/*
 * p4est_wrap.c
 * Wrapper for p4est and sc libraries.
 *
 *  Created on: Aug 19, 2022
 *      Author: Adam Peplinski
 */

// libsc speciffic definitions
#define SC_ENABLE_MPI
#define SC_ENABLE_MPIIO

// p4est speciffic definitions
// This is dimension rleated; for now 3D only
#define N_DIM 3

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

#include <sc.h>

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

#include "./p4est_wrap.h"

/* Global variables
 *  notice I use tree_nek->user_pointer to store ghost quadrants data
 */
static p4est_connectivity_t *connect_neko = NULL; /**< connectivity structure */
static p4est_t *tree_neko = NULL; /**< tree structure */
static p4est_mesh_t *mesh_neko = NULL; /**< mesh structure */
static p4est_ghost_t *ghost_neko = NULL; /**< ghost zone structure */
static p4est_nodes_t *nodes_neko = NULL; /**< vertex numbering structure */
static p4est_lnodes_t *lnodes_neko = NULL; /**< GLL node numbering structure */
static p4est_t *tree_neko_compare = NULL; /**< tree structure to perform comparison */

/* Wrappers */
/* initialise/finalise */
void wp4est_init(MPI_Fint fmpicomm,
		 int catch_signals, int print_backtrace,
		 int log_threshold)
{
  MPI_Comm mpicomm;
  mpicomm = MPI_Comm_f2c(fmpicomm);
  sc_init (mpicomm, catch_signals, print_backtrace,
  	   NULL, log_threshold);
  p4est_init(NULL, log_threshold);
}

void wp4est_finalize(int log_priority)
{
  sc_package_print_summary (log_priority);
  sc_finalize ();
}

/* Connectivity; only required routines */
void wp4est_cnn_del() {
  if (connect_neko) p4est_connectivity_destroy(connect_neko);
  connect_neko =  NULL;
}

void wp4est_cnn_valid(int * is_valid) {
  *is_valid = p4est_connectivity_is_valid(connect_neko);
}

/* tree_ management */
void wp4est_tree_del() {
  if (tree_neko) p4est_destroy(tree_neko);
  tree_neko = NULL;
}

void wp4est_tree_valid(int * is_valid) {
  *is_valid = p4est_is_valid(tree_neko);
}

void wp4est_tree_save(int save_data, char filename[]) {
  p4est_save_ext(filename, tree_neko, save_data,0);
}

void wp4est_tree_load(MPI_Fint fmpicomm, int load_data, char filename[]) {
  MPI_Comm mpicomm;
  mpicomm = MPI_Comm_f2c(fmpicomm);
  tree_neko = p4est_load_ext(filename, mpicomm, sizeof(user_data_t), load_data,1,0,
			     NULL, &connect_neko);
  tree_neko->user_pointer = NULL;
}

/* tree and grid info */
void wp4est_ghost_new() {
  if (ghost_neko) p4est_ghost_destroy(ghost_neko);
  ghost_neko = p4est_ghost_new(tree_neko, P4EST_CONNECT_FULL);
  /* ghost data */
  if (tree_neko->user_pointer) P4EST_FREE(tree_neko->user_pointer);
  tree_neko->user_pointer = (void *) P4EST_ALLOC (user_data_t, ghost_neko->ghosts.elem_count);
  p4est_ghost_exchange_data (tree_neko, ghost_neko, (user_data_t *) tree_neko->user_pointer);
}

void wp4est_ghost_del() {
  if (ghost_neko) p4est_ghost_destroy(ghost_neko);
  ghost_neko = NULL;
  /* ghost data */
  if (tree_neko->user_pointer) P4EST_FREE(tree_neko->user_pointer);
  tree_neko->user_pointer = NULL;
}

void wp4est_mesh_new() {
  if (mesh_neko) p4est_mesh_destroy(mesh_neko);
  mesh_neko = p4est_mesh_new(tree_neko, ghost_neko, P4EST_CONNECT_FULL);
}

void wp4est_mesh_del() {
  if (mesh_neko) p4est_mesh_destroy(mesh_neko);
  mesh_neko = NULL;
}

void wp4est_nodes_new() {
  if (nodes_neko) p4est_nodes_destroy(nodes_neko);
  nodes_neko = p4est_nodes_new(tree_neko, ghost_neko);
}

void wp4est_nodes_del() {
  if (nodes_neko) p4est_nodes_destroy(nodes_neko);
  nodes_neko = NULL;
}

void wp4est_lnodes_new(int degree) {
  if (lnodes_neko) p4est_lnodes_destroy(lnodes_neko);
  lnodes_neko = p4est_lnodes_new(tree_neko, ghost_neko, degree);
}

void wp4est_lnodes_del() {
  if (lnodes_neko) p4est_lnodes_destroy(lnodes_neko);
  lnodes_neko = NULL;
}

/* p4est internal load balance */
void wp4est_part(int partforcoarsen) {
  p4est_partition(tree_neko, partforcoarsen, NULL);
}

/** @brief Iterate over faces to correct Neko boundary conditions
 *
 * @details Required by wp4est_bc_check
 *
 * @param info
 * @param user_data
 */
void iter_bc_chk(p4est_iter_face_info_t * info, void *user_data) {
  user_data_t *data;

  int mpirank = info->p4est->mpirank;
  p4est_gloidx_t gfirst_quad =  info->p4est->global_first_quadrant[mpirank];
  /* ghost data */
  user_data_t *ghost_data = (user_data_t *) info->p4est->user_pointer;

  sc_array_t *sides = &(info->sides);
  p4est_iter_face_side_t *side;
  int nside = (int) sides->elem_count;


  int il, jl, kl; // loop index
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
      /* external face; no E (0), P (-1) bc */
      side = p4est_iter_fside_array_index_int (sides, 0);
      face = (int) side->face;
      if (side->is_hanging) {
	/* hanging face */
	for (jl = 0; jl < P4EST_HALF; jl++) {
	  if (!side->is.hanging.is_ghost[jl]) {
	    /* local node */
	    data = (user_data_t *) side->is.hanging.quad[jl]->p.user_data;
	    /* check connection type (should not be E or P) */
	    if (data->bc[face] == 0 || data->bc[face] == -1) {
	      iwl = (int) (gfirst_quad + side->treeid + side->is.hanging.quadid[jl]);
	      printf("Connectivity error: %i %i %i %i\n",
		     tree_neko->mpirank, iwl, face, data->bc[face]);
	      printf("external face marked as internal.\n");
	      SC_ABORT("Aborting: iter_bc_chk\n");
	    }
	  } // is ghost
	} // children loop

      } else {
	if (!side->is.full.is_ghost) {
	  /* local node */
	  data = (user_data_t *) side->is.full.quad->p.user_data;
	  /* check connection type (should not be E or P) */
	  if (data->bc[face] == 0 || data->bc[face] == -1) {
	    iwl = (int) (gfirst_quad + side->treeid + side->is.full.quadid);
	    printf("Connectivity error: %i %i %i %i\n",
		   tree_neko->mpirank, iwl, face, data->bc[face]);
	    printf("external face marked as internal.\n");
	    SC_ABORT("Aborting: iter_bc_chk\n");
	  }
	} // is ghost
      } // hanging

    } else {
      /* internal face; any face type possible */
      /* collect quad type;
       * for different values of imsh V-type elements should point to external bc
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
	      /* velocity bc for V-type element neighbor to T-type element
	       * should not be E or P
	       */
	      if (imsh[0] != imsh[1]) {
		if (data->bc[face] == 0 || data->bc[face] == -1) {
		  iwl = (int) (gfirst_quad + side->treeid + side->is.hanging.quadid[jl]);
		  printf("Connectivity error: %i %i %i %i\n",
			 tree_neko->mpirank, iwl, face, data->bc[face]);
		  printf("velocity external face marked as internal.\n");
		  SC_ABORT("Aborting: iter_bc_chk\n");
		}
	      } else {
		/* not velocity or not V-T meshes boundary - all
		 * internal elements
		 */
		if (data->bc[face] != 0 && data->bc[face] != -1) {
		  iwl = (int) (gfirst_quad + side->treeid + side->is.hanging.quadid[jl]);
		  printf("Connectivity error: %i %i %i %i\n",
			 tree_neko->mpirank, iwl, face,data->bc[face]);
		  printf("internal face marked as external.\n");
		  SC_ABORT("Aborting: iter_bc_chk\n");
		}
	      }
	    } // is ghost
	  } // children loop

	} else {
	  if (!side->is.full.is_ghost) {
	    /* local node */
	    data = (user_data_t *) side->is.full.quad->p.user_data;
	    /* velocity bc for V-type element neighbor to T-type element
	     * should not be E or P
	     */
	    if (imsh[0] != imsh[1]) {
	      if (data->bc[face] == 0|| data->bc[face] == -1) {
		iwl = (int) (gfirst_quad + side->treeid + side->is.full.quadid);
		printf("Connectivity error: %i %i %i %i\n",
		       tree_neko->mpirank, iwl, face, data->bc[face]);
		printf("velocity external face marked as internal.\n");
		SC_ABORT("Aborting: iter_bc_chk\n");
	      }
	    } else {
	      /* not velocity or not V-T meshes boundary - all
	       * internal elements
	       */
	      if (data->bc[face] != 0 && data->bc[face] != -1) {
		iwl = (int) (gfirst_quad + side->treeid + side->is.full.quadid);
		printf("Connectivity error: %i %i %i %i\n",
		       tree_neko->mpirank, iwl, face, data->bc[face]);
		printf("internal face marked as external.\n");
		SC_ABORT("Aborting: iter_bc_chk\n");
	      }
	    }
	  } // is ghost
	} // hanging
      } // side loop
    } // forest internal/external face
  } else {
    /* face is on the interior of the tree; all faces are E (0)
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
	    if (data->bc[face] != 0 ) {
	      iwl = (int) (gfirst_quad + side->treeid + side->is.hanging.quadid[jl]);
	      printf("Connectivity error: %i %i %i %i\n",
		     tree_neko->mpirank, iwl, face, data->bc[face]);
	      printf("tree internal face marked as external.\n");
	      SC_ABORT("Aborting: iter_bc_chk\n");
	    }
	  } // is ghost
	} // children loop

      } else {
	if (!side->is.full.is_ghost) {
	  /* local node */
	  data = (user_data_t *) side->is.full.quad->p.user_data;
	  if (data->bc[face] != 0) {
	    iwl = (int) (gfirst_quad + side->treeid + side->is.full.quadid);
	    printf("Connectivity error: %i %i %i %i\n",
		   tree_neko->mpirank, iwl, face, data->bc[face]);
	    printf("tree internal face marked as external.\n");
	    SC_ABORT("Aborting: iter_bc_chk\n");
	  }
	} // is ghost
      } // hanging
    } // side loop
  } // tree internal/external

}

/* Check boundary conditions */
void wp4est_bc_check() {
#ifdef P4_TO_P8
  p4est_iterate(tree_neko, ghost_neko, NULL, NULL, iter_bc_chk, NULL, NULL);
#else
  p4est_iterate(tree_neko, ghost_neko, NULL, NULL, iter_bc_chk, NULL);
#endif
}

/* routines for data exchange between Neko and p4est */

/** @brief Iterate over elements to count V-mesh elements
 *
 * @details Required by wp4est_msh_get_size
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
void wp4est_msh_get_size(int * mdim, int64_t * nelgt, int64_t * nelgto,
			 int32_t * nelt, int * nelv, int * maxl) {
  int lmax[2];
  // mesh dimension
  *mdim = (int) P4EST_DIM;
  // get global number of quadrants
  *nelgt = tree_neko->global_num_quadrants;
  // zero based global position of local quadrants
  *nelgto = tree_neko->global_first_quadrant[tree_neko->mpirank];
  // number of local quadrants
  *nelt = tree_neko->local_num_quadrants;

  // count number of V-mesh elements and find current max level
  lmax[0] = 0;
  lmax[1] = 0;

#ifdef P4_TO_P8
  p4est_iterate(tree_neko, ghost_neko,(void *) &lmax, count_mshv,
		NULL, NULL, NULL);
#else
  p4est_iterate(tree_neko, ghost_neko,(void *) &lmax, count_mshv,
		NULL, NULL);
#endif

  *nelv = lmax[0];
  *maxl = lmax[1];
}

/* get node list size */
void wp4est_nds_get_size(int * nowin, int * nowsh, int * oowin,
			 int * nin, int * nhf, int * nhe) {
  if (nodes_neko) {
    // number of owned independent nodes
    *nowin = (int) nodes_neko->num_owned_indeps;
    // number of owned shared
    *nowsh = (int) nodes_neko->num_owned_shared;
    // position of the first independent owned node
    *oowin = (int) nodes_neko->offset_owned_indeps;
    // numbers of local independent and hanging nodes (face and edge)
    *nin = (int) nodes_neko->indep_nodes.elem_count;
    *nhf = (int) nodes_neko->face_hangings.elem_count;
#ifdef P4_TO_P8
    *nhe = (int) nodes_neko->edge_hangings.elem_count;
#endif
  } else {
    SC_ABORT("nodes_neko not allocated; aborting: wp4est_nds_get_size\n");
  }
}

/* independent node list */
void wp4est_nds_get_ind(int64_t * nglid, int * nown, double * ncoord) {
  int il, jl;
  int64_t id;
  p4est_indep_t *node;
  double vxyz[3];
  if (nodes_neko) {
    sc_array_t  *indep_nodes = &(nodes_neko->indep_nodes);
    // numbers of local independent nodes
    const int ni = (int) nodes_neko->indep_nodes.elem_count;
    // independent node offset
    const int oi = (int) nodes_neko->offset_owned_indeps;
    // loop over nodes owned by other mpi rank
    for(il=0;il<oi;++il){
      // extract node
      node = (p4est_indep_t *) sc_array_index (indep_nodes,il);
      // get global id
      id = (int64_t) node->p.piggy3.local_num + 1; // conversion to fortran numberring
      for(jl=0;jl<nodes_neko->nonlocal_ranks[il];++jl) {
	id += (int64_t) nodes_neko->global_owned_indeps[jl];
      }
      nglid[il] = id;
      // node owner
      nown[il] = (int) nodes_neko->nonlocal_ranks[il];
      // get physical coordinates
      p4est_qcoord_to_vertex (connect_neko, node->p.piggy3.which_tree,
			      node->x, node->y,
#ifdef P4_TO_P8
			      node->z,
#endif
			      vxyz);
      // copy coordinates
      for(jl=0;jl<N_DIM;++jl){
	ncoord[il*N_DIM+jl] = vxyz[jl];
      }
    }
    // loop over nodes owned by this mpi rank
    // local id offset
    id = (int64_t) 1; // conversion to fortran numberring
    for(jl=0;jl<tree_neko->mpirank;++jl) {
      id += (int64_t) nodes_neko->global_owned_indeps[jl];
    }
    for(il=oi;il<ni;++il){
      // extract node
      node = (p4est_indep_t *) sc_array_index (indep_nodes,il);
      // get global id
      nglid[il] = id + (int64_t) node->p.piggy3.local_num;
      // node owner
      nown[il] = (int) tree_neko->mpirank;
      // get physical coordinates
      p4est_qcoord_to_vertex (connect_neko, node->p.piggy3.which_tree,
			      node->x, node->y,
#ifdef P4_TO_P8
			      node->z,
#endif
			      vxyz);
      // copy coordinates
      for(jl=0;jl<N_DIM;++jl){
	ncoord[il*N_DIM+jl] = vxyz[jl];
      }
    }
  } else {
    SC_ABORT("nodes_neko not allocated; aborting: wp4est_nds_get_ind\n");
  }
}

/* face hanging node list */
void wp4est_nds_get_hfc(int * depend, double * ncoord) {
  int il, jl;
#ifdef P4_TO_P8
  const int ndep = 4;
  p8est_hang4_t *node;
#else
  const int ndep = 2;
  p4est_hang2_t *node;
#endif
  double vxyz[3];

  if (nodes_neko) {
    sc_array_t  *nodes = &(nodes_neko->face_hangings);
    // numbers of local face hanging nodes
    const int ni = (int) nodes_neko->face_hangings.elem_count;
    // loop over nodes
    for(il=0;il<ni;++il){
      // extract node
#ifdef P4_TO_P8
      node = (p8est_hang4_t *) sc_array_index (nodes,il);
#else
      node = (p4est_hang2_t *) sc_array_index (nodes,il);
#endif
      // get independent nodes mapping
      for(jl=0;jl<ndep;++jl){
        depend[il*ndep+jl] = (int) node->p.piggy.depends[jl] +1; // conversion to fortran numberring
      }
      // get physical coordinates
      p4est_qcoord_to_vertex (connect_neko, node->p.piggy.which_tree,
                              node->x, node->y,
#ifdef P4_TO_P8
                              node->z,
#endif
                              vxyz);
      // copy coordinates
      for(jl=0;jl<N_DIM;++jl){
        ncoord[il*N_DIM+jl] = vxyz[jl];
      }
    }
  }else {
    SC_ABORT("nodes_neko not allocated; aborting: wp4est_nds_get_hfc\n");
  }
}


/* edge hanging node list */
void wp4est_nds_get_hed(int * depend, double * ncoord) {
  int il, jl;
#ifdef P4_TO_P8
  const int ndep = 2;
  p8est_hang2_t *node;
  double vxyz[3];

  if (nodes_neko) {
    sc_array_t  *nodes = &(nodes_neko->edge_hangings);
    // numbers of local face hanging nodes
    const int ni = (int) nodes_neko->edge_hangings.elem_count;
    // loop over nodes
    for(il=0;il<ni;++il){
      // extract node
      node = (p8est_hang2_t *) sc_array_index (nodes,il);
      // get independent nodes mapping
      for(jl=0;jl<ndep;++jl){
        depend[il*ndep+jl] = (int) node->p.piggy.depends[jl] + 1; // conversion to fortran numberring
      }
      // get physical coordinates
      p4est_qcoord_to_vertex (connect_neko, node->p.piggy.which_tree,
                              node->x, node->y,
                              node->z,
                              vxyz);
      // copy coordinates
      for(jl=0;jl<N_DIM;++jl){
        ncoord[il*N_DIM+jl] = vxyz[jl];
      }
    }
  }else {
    SC_ABORT("nodes_neko not allocated; aborting: wp4est_nds_get_hed\n");
  }
#endif
}

/* get vertex to node mapping */
void wp4est_nds_get_vmap(int * vmap) {
  int il, jl;

  if (nodes_neko) {
    // number of local elements
    const int vi = (int) nodes_neko->num_local_quadrants;
    // quad to vertex local map
    for (il = 0; il < vi; ++il) {
      for (jl = 0; jl < P4EST_CHILDREN; ++jl) {
	vmap[il * P4EST_CHILDREN + jl] = (int) nodes_neko->local_nodes[il * P4EST_CHILDREN + jl] + 1; // conversion to fortran numberring
      }
    }
  }else {
    SC_ABORT("nodes_neko not allocated; aborting: wp4est_nds_get_vnmap\n");
  }
}


/* data type for mesh data transfer between nek5000 and p4est*/
typedef struct transfer_data_s {
  int64_t *gidx; /**< pointer to global element index array */
  int *level; /**< pointer to element level array */
  int *igrp; /**< pointer to element group array */
  int *crv; /**< pointer to face projection array */
  int *bc; /**< pointer to boundary condition array */
  double *coord; /**< pointer to approximate vertex coordinates array */
  int *falg; /**< face alignment */
} transfer_data_t;

/** @brief Iterate over element volumes to transfer approximate coordinates of the element vertices
 *
 * @details Required by wp4est_elm_get_dat
 *
 * @param info
 * @param user_data
 */
void iter_datav(p4est_iter_volume_info_t * info, void *user_data) {
  user_data_t *data = (user_data_t *) info->quad->p.user_data;
  transfer_data_t *trans_data = (transfer_data_t *) user_data;

  // which quad (local and global element number)
  p4est_tree_t *tree;
  p4est_locidx_t iwl;
  p4est_gloidx_t iwlt, iwg;
  int iwli, il, jl;
  p4est_quadrant_t node;
  double vxyz[3];

  // get quad number
  tree = p4est_tree_array_index(info->p4est->trees, info->treeid);
  // local quad number
  iwl = info->quadid + tree->quadrants_offset;
  iwli = (int) iwl;
  // global quad number
  iwg = (p4est_gloidx_t) info->p4est->global_first_quadrant[info->p4est->mpirank] + iwl;

  // global quadrant index
  trans_data->gidx[iwli] = iwg + 1; // conversion to fortran numberring

  // quadrant level
  trans_data->level[iwli] = (int) info->quad->level;

  // element group mark
  trans_data->igrp[iwli] = data->igrp;

  // curvature and boundary data
  for (il = 0; il < P4EST_FACES; il++) {
    trans_data->crv[iwli*P4EST_FACES+il] = data->crv[il];
    trans_data->bc[iwli*P4EST_FACES+il] = data->bc[il];
  }

  // get corner coordinates
  for (il=0;il < P4EST_CHILDREN; ++il){
    p4est_quadrant_corner_node (info->quad, il, &node);
    p4est_qcoord_to_vertex (info->p4est->connectivity, info->treeid,
			    node.x, node.y,
#ifdef P4_TO_P8
			    node.z,
#endif
			    vxyz);
    // copy coordinates
    for(jl=0;jl<N_DIM;++jl){
      trans_data->coord[(iwli*P4EST_CHILDREN+il)*N_DIM+jl] = vxyz[jl];
    }
  }
}

/** @brief Iterate over faces to get alignment
 *
 * @details Required by wp4est_elm_get_dat
 *
 * @param info
 * @param user_data
 */
void iter_algf(p4est_iter_face_info_t * info, void *user_data) {
  // which quad (local element number)
  p4est_tree_t *tree;
  p4est_locidx_t iwl;
  int iwlt;
  // orientation and side info
  int orient = (int) info->orientation;
  sc_array_t *sides = &(info->sides);
  p4est_iter_face_side_t *side;
  int nside = (int) sides->elem_count;
  transfer_data_t *trans_data = (transfer_data_t *) user_data;
  //int *fcs_arr = (int *) user_data;
  int il, jl;
  int iref, pref, pset;
  int8_t face[2], ftmp;

  /* compare orientation of different trees*/
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
	  tree = p4est_tree_array_index(info->p4est->trees, side->treeid);
	  // local quad number
	  iwl =  side->is.hanging.quadid[jl] + tree->quadrants_offset;
	  iwlt = (int) iwl;
	  iwlt = iwlt*P4EST_FACES + (int) side->face;
	  trans_data->falg[iwlt] = orient;
	}
      }
    } else {
      if (!side->is.full.is_ghost) {
	// local node
	tree = p4est_tree_array_index(info->p4est->trees, side->treeid);
	// local quad number
	iwl =  side->is.full.quadid + tree->quadrants_offset;
	iwlt = (int) iwl;
	iwlt = iwlt*P4EST_FACES + (int) side->face;
	trans_data->falg[iwlt] = orient;
      }
    }
  }
}

/* get element info */
void wp4est_elm_get_dat(int64_t * gidx, int * level, int * igrp, int * crv,
			int * bc,double * coord, int * falg) {
  transfer_data_t transfer_data;

  transfer_data.gidx = gidx;
  transfer_data.igrp = igrp;
  transfer_data.level = level;
  transfer_data.crv = crv;
  transfer_data.bc = bc;
  transfer_data.coord = coord;
  transfer_data.falg = falg;
#ifdef P4_TO_P8
  p4est_iterate(tree_neko, ghost_neko,(void *) &transfer_data, iter_datav,
		iter_algf, NULL, NULL);
#else
  p4est_iterate(tree_neko, ghost_neko,(void *) &transfer_data, iter_datav,
		iter_algf, NULL);
#endif
}

/* degrees of freedom numbering */
void wp4est_elm_get_lnode(int * lnnum, int * lnown, int64_t * lnoff, int * lnodes) {

  int il, jl, kl;
  if (lnodes_neko) {
    const int vnd = (int) lnodes_neko->vnodes;
    const int owned = (int) lnodes_neko->owned_count;
    const int local = (int) lnodes_neko->num_local_nodes;
    const p4est_gloidx_t offset = lnodes_neko->global_offset;
    *lnnum = (int) local;
    *lnown = (int) owned;
    *lnoff = (int64_t) offset;
    for (il = 0; il < lnodes_neko->num_local_elements*vnd; ++il) {
      lnodes[il] = (int) lnodes_neko->element_nodes[il] + 1; // conversion to fortran numberring
    }
  } else {
    SC_ABORT("lnodes_nek not allocated; aborting: wp4est_elem_get_lnode\n");
  }
}

/* get shares list size */
void wp4est_sharers_get_size(int * nrank, int * nshare) {
  int il, jl;
  p4est_lnodes_rank_t *lnode;
  if (lnodes_neko) {
    // number of ranks in a sharers array
    *nrank = (int) lnodes_neko->sharers->elem_count;
    sc_array_t  *sharers = lnodes_neko->sharers;
    jl = 0;
    // loop over mpi rank sharing the node
    for(il=0;il < lnodes_neko->sharers->elem_count;++il){
      lnode = (p4est_lnodes_rank_t *) sc_array_index (sharers,il);
      // number of nodes shared with given rank
      jl += (int) lnode->shared_nodes.elem_count;
      /*      printf("Sharers: %i %i %i %i %i %i %i %i\n",
	     tree_neko->mpirank, il, lnode->rank, (int) lnode->shared_nodes.elem_count,
	     lnode->shared_mine_offset, lnode->shared_mine_count,
	     lnode->owned_offset, lnode->owned_count);*/
    }
    *nshare = jl;
  } else {
    SC_ABORT("lnodes_neko not allocated; aborting: wp4est_sharers_get_size\n");
  }
}

/* independent lnode list */
void wp4est_sharers_get_ind(int64_t * nglid, int * lrank, int * loff, int * lshare) {
  int il, jl;
  p4est_lnodes_rank_t *lnode;
  p4est_locidx_t *snode;
  if (lnodes_neko) {
    const int owned = (int) lnodes_neko->owned_count;
    const int local = (int) lnodes_neko->num_local_nodes;
    const p4est_gloidx_t offset = lnodes_neko->global_offset;
    // get global index list of local nodes
    // local owned
    for (il = 0; il < owned ; ++il) {
      nglid[il] = (p4est_gloidx_t) 1 + offset + il;
    }
    // local not owned
    for (il = owned; il < local ; ++il) {
      nglid[il] = (p4est_gloidx_t) 1 + lnodes_neko->nonlocal_nodes[il-owned]; // conversion to fortran numberring
    }
    // get node sharing info
    sc_array_t  *sharers = lnodes_neko->sharers;
    // loop over mpi rankssharing the node
    loff[0] = 1;
    for(il=0;il < lnodes_neko->sharers->elem_count;++il){
      lnode = (p4est_lnodes_rank_t *) sc_array_index (sharers,il);
      lrank[il] = (int) lnode->rank;
      loff[il+1] = loff[il] + (int) lnode->shared_nodes.elem_count;
      sc_array_t  *node_list = &(lnode->shared_nodes);
      for(jl=0;jl < lnode->shared_nodes.elem_count;++jl){
	snode = (p4est_locidx_t *) sc_array_index(node_list,jl);
	lshare[loff[il]-1+jl] = 1 + *snode; // conversion to fortran numberring
      } 
    }
  } else {
    SC_ABORT("lnodes_neko not allocated; aborting: wp4est_sharers_get_ind\n");
  }
}


/* get hanging element/face/edge information */
void wp4est_hang_get_info(int * hang_elm, int * hang_fsc, int * hang_edg) {

  int il, jl;
  int hanging_face[P4EST_FACES];
#ifdef P4_TO_P8
  int hanging_edge[P8EST_EDGES];
#endif
  if (lnodes_neko) {
    for (il = 0; il < lnodes_neko->num_local_elements; ++il) {
      for (jl = 0; jl < P4EST_FACES; ++jl) {
	hanging_face[jl] = -1;
      }
#ifdef P4_TO_P8
      for (jl = 0; jl < P8EST_EDGES; ++jl) {
	hanging_edge[jl] = -1;
      }
      hang_elm[il] = p4est_lnodes_decode(lnodes_neko->face_code[il],
					 hanging_face, hanging_edge);
      for (jl = 0; jl < P4EST_FACES; ++jl) {
	hang_fsc[il * P4EST_FACES + jl] = hanging_face[jl];
      }
      for (jl = 0; jl < P8EST_EDGES; ++jl) {
	hang_edg[il * P8EST_EDGES + jl] = hanging_edge[jl];
      }
#else
      hang_elm[il] = p4est_lnodes_decode(lnodes_neko->face_code[il],hanging_face);
      for (jl=0;jl<P4EST_FACES;++jl) {
	hang_fsc[il*P4EST_FACES +jl] = hanging_face[jl];
      }
#endif
    }
  } else {
    SC_ABORT("lnodes_neko not allocated; aborting: wp4est_msh_get_tplg\n");
  }
}

/* extract family information */
void wp4est_fml_get_info(int64_t * family, int * nelf) {
  p4est_topidx_t jt;
  p4est_tree_t *tree;
  sc_array_t *tquadrants;
  int isfamily;
  size_t zz, incount, window;
  p4est_quadrant_t *ci[P4EST_CHILDREN];
  int64_t igs, its, iqs;
  int iqf;

  // initialise number of family quads
  iqf = 0;

  // global quadrants shift
  igs = (int64_t) tree_neko->global_first_quadrant[tree_neko->mpirank] + 1; // conversion to fortran numberring

  /* loop over all local trees */
  for (jt = tree_neko->first_local_tree; jt <= tree_neko->last_local_tree; ++jt) {
    tree = p4est_tree_array_index (tree_neko->trees, jt);
    tquadrants = &tree->quadrants;

    window = 0;  // start position of sliding window in array
    incount = tquadrants->elem_count;  // number of quadrants
    its = (int64_t) tree->quadrants_offset; // local quadrant shift
    while (window + P4EST_CHILDREN <= incount) {
      isfamily = 1;
      for (zz = 0; zz < P4EST_CHILDREN; ++zz) {
	ci[zz] = p4est_quadrant_array_index (tquadrants, window + zz);
	if (zz != (size_t) p4est_quadrant_child_id (ci[zz])) {
	  isfamily = 0;
	  break;
	}
      }
      P4EST_ASSERT (!isfamily || p4est_quadrant_is_familypv (ci));
      if (isfamily) {
	// get global element number in the family
	for (zz = 0; zz < P4EST_CHILDREN; ++zz) {
	  iqs = (int64_t) window + zz;
	  family[iqf] = igs + its + iqs;
	  ++iqf;
	}
	window += P4EST_CHILDREN;
      }
      else {
	++window;
      }
    }
  }

  // copy number of family quads
  *nelf = iqf;
}
