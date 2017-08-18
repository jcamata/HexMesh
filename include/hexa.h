/* 
 * File:   hexa.h
 * Author: camata
 *
 * Created on March 20, 2015, 12:00 PM
 */

#ifndef HEXA_H
#define	HEXA_H

#include <stdio.h>
#include <stdlib.h>
#include <sc.h>
#include <sc_containers.h>
#include <gts.h>
#include <vector>

#include "hilbert.h"

#define HEXA_DEBUG_

typedef struct {
    int32_t id;
    int32_t x,y,z;
    int color;
} octant_node_t;

typedef struct {
    uint64_t id;
    bool ref;
    bitmask_t coord[2];
} octant_edge_t;

typedef struct {
    int32_t    x,y,z;
    int8_t     level;
    int     pad;
    int     tem;
    int8_t     pml_id;
    int        n_mat;
    unsigned int edge_id[12];
    bool edge_ref[12];
    octant_node_t nodes[8];
    octant_edge_t edge[12];
    int32_t    id;
} octant_t;

typedef struct shared_node
{
  int32_t      id;
  int32_t       x, y, z;
  int32_t       listSz;
  int32_t      rankList[8];
} shared_node_t;

typedef struct shared_edge{
  uint64_t     id;
  int32_t      listSz;
  int32_t      rankList[4];
} shared_edge_t;

typedef struct 
{
    int rank;
    sc_array_t idxs;
    //sc_array_t message;
} message_t;

typedef struct
{
    sc_array_t SendTo;
    sc_array_t RecvFrom;
    uint32_t max_recvbuf_size;
    uint32_t max_sendbuf_size;
    int     nrequests;
} comm_map_t;


typedef struct 
{
    GtsSurface *s;      // Surface
    GNode      *bbt;    // Boundary box tree used in 
                        // interceptions operations
    GtsBBox    *bbox;   // Surface boundry box
} GeometryData;


typedef struct {
    
    int64_t total_n_elements;
    int64_t total_n_nodes;
    int64_t total_n_edges;
    int64_t total_n_mat;
    
    int32_t local_n_elements;
    int32_t local_n_nodes;
    int32_t local_n_edges;
    int32_t max_step;
    
    sc_array_t      elements;
    sc_array_t      edges;
    sc_array_t      nodes;
    sc_array_t      shared_nodes; 
    sc_array_t      shared_edges;
    sc_array_t      edges_ref;
    int64_t *global_id;
    //uint64_t *global_edge_id;
    int32_t *part_nodes;
    int32_t ncellx;
    int32_t ncelly;
    int32_t ncellz;
    int32_t max_z;
    int     mpi_rank;
    int     mpi_size;
    int     max_levels;
    int32_t x_start;
    int32_t x_end;
    int32_t y_start;
    int32_t y_end;
    double x_startc;
    double x_endc;
    double y_startc;
    double y_endc;
    int neighbors[9];
    comm_map_t comm_map;
    comm_map_t  comm_map_edge;
#ifdef HEXA_DEBUG_
    FILE* fdbg;
#endif
    GeometryData gdata;

} hexa_tree_t;



void copy_octant(octant_t *orig, octant_t* dest);

void hexa_init(int argc, char* argv[], hexa_tree_t* mesh);

void hexa_finalize(hexa_tree_t* mesh);

void hexa_tree_init(hexa_tree_t* mesh, int max_levels);

void hexa_tree_destroy(hexa_tree_t* mesh);

void hexa_tree_cube(hexa_tree_t* mesh);

int  hexa_tree_write_vtk(hexa_tree_t* mesh,  const char *filename);

void hexa_transition_element(hexa_tree_t* mesh, int i, int j, int k, int step, int level);

void hexa_processors_interval(hexa_tree_t* mesh);

void hexa_mesh(hexa_tree_t* tree);

int  hexa_mesh_write_vtk(hexa_tree_t* mesh,  const char *filename, std::vector<double> *coords);

void hexa_mesh_write_unv(hexa_tree_t* mesh, const char* root_name, std::vector<double> *coords);

void hexa_mesh_write_msh(hexa_tree_t* mesh, const char* root_name, std::vector<double> *coords);

void hexa_mesh_write_h5(hexa_tree_t* mesh, const char* root_name, std::vector<double> coords);

void hexa_debug_face_hanging(hexa_tree_t* mesh);

unsigned processors_hash_fn (const void *v, const void *u);

int processors_equal_fn (const void *v1, const void *v2, const void *u);

void hexa_insert_shared_node(sc_hash_array_t    *shared_nodes, octant_node_t* node, int processor);

unsigned node_hash_fn (const void *v, const void *u);

int node_equal_fn (const void *v1, const void *v2, const void *u);

#endif	/* HEXA_H */

