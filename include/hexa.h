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

typedef struct
{
	int32_t id;
	int32_t x,y,z;
	int color;
	int fixed;
} octant_node_t;

typedef struct
{
	uint64_t id;
	bool ref;
	bitmask_t coord[2];
} octant_edge_t;

typedef struct
{
	uint64_t id;
	int nodes[4];
	bool ref;
} octant_face_t;

typedef struct {
	double coord[3];
	int node_id;
} node_t;

typedef struct
{
	uint64_t id;
	int32_t  list_elem;
	//int32_t  list_edges;
	int32_t  elem[16]; // we must check this values
	//int32_t  edges[3*16]; // we must check this values
} octant_vertex_t;

typedef struct
{
	int32_t    x,y,z;
	int8_t     level;
	int     pad;
	int     cut;
	int     tem;
	int     inipad;
	int     initem;
	int8_t     pml_id;
	int        n_mat;
	unsigned int edge_id[12];
	unsigned int father;
	bool edge_ref[12];
	bool ghost;
	octant_node_t nodes[8];
	octant_edge_t edge[12];
	octant_face_t face[6];
	int64_t    id;
} octant_t;

typedef struct
{
	int64_t id[8];
	bool    cut;
	bool    face[6];
	bool    edge[12];
	int     mat[8];
} octree_t;

typedef struct shared_node
{
	int32_t      id;
	int32_t       x, y, z;
	int32_t       listSz;
	int32_t      rankList[8];
} shared_node_t;

typedef struct shared_edge
{
	uint64_t     id;
	int32_t      listSz;
	int32_t      rankList[4];
} shared_edge_t;

typedef struct
{
	uint64_t     id;
	int32_t      nodes[4];
	int32_t      listSz;
	int32_t      rankList;
} shared_face_t;

typedef struct shared_element
{
	uint64_t     id;
	int32_t       x, y, z;
	int32_t      listSz;
	int32_t      rankList[4];
} shared_element_t;

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
	sc_array_t      faces;
	sc_array_t      vertex;
	sc_array_t		oct;

	sc_array_t      shared_nodes;
	sc_array_t      shared_edges;
	sc_array_t      shared_elements;
	sc_array_t      shared_faces;
	sc_array_t      shared_vertex;

	sc_array_t      edges_ref;
	int64_t *global_id;

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

int const EdgeVerticesMap[12][2] = {
		{0, 1}, // Edge 0
		{1, 2}, // Edge 1
		{2, 3}, // Edge 2
		{3, 0}, // Edge 3
		{0, 4}, // Edge 4
		{1, 5}, // Edge 5
		{2, 6}, // Edge 6
		{3, 7}, // Edge 7
		{4, 5}, // Edge 8
		{5, 6}, // Edge 9
		{6, 7}, // Edge 10
		{7, 4}  // Edge 11
};

int const EdgeVerticesMap_surf_diagonal[12][2] = {
		{0, 5},
		{1, 4},
		{3, 6},
		{2, 7},
		{0, 7},
		{3, 4},
		{1, 6},
		{2, 5},
		{4, 6},
		{5, 7},
		{0, 2},
		{1, 3}
};

int const EdgeVerticesMap_vol_diagonal[4][2] = {
		{0, 6},
		{1, 7},
		{2, 4},
		{3, 5}
};

int const FaceEdgesMap[6][4] = {
		{4, 11, 7, 3},
		{5, 1, 6, 9},
		{0, 5, 8, 4},
		{2, 6, 10, 7},
		{8, 9, 10, 11},
		{0, 1, 2, 3}
};

int const FaceNodesMap[6][4] = {
		{0,4,7,3},
		{1,5,6,2},
		{0,4,5,1},
		{2,6,7,3},
		{4,5,6,7},
		{0,1,2,3}
};

int const VertexEdgeMap[8][3] = {
		{0,3,4}, //Vertex 0
		{0,1,5}, //Vertex 1
		{1,2,6}, //Vertex 2
		{2,3,7}, //Vertex 3
		{4,8,11}, //Vertex 4
		{5,8,9}, //Vertex 5
		{6,9,10}, //Vertex 6
		{7,10,11}  //Vertex 7
};

int const EdgeEdgeMap[12][4] = {
		{3, 4, 1, 5}, // Edge 0
		{0, 5, 2, 6}, // Edge 1
		{1, 6, 3, 7}, //      2
		{2, 7, 0, 4,}, //3
		{0, 3, 8, 11}, //4
		{1, 0, 9, 8}, //5
		{2, 1, 10, 9}, //6
		{3, 2, 11, 10}, //7
		{11, 4, 9, 5}, //8
		{8, 5, 10, 6}, //9
		{9, 6, 11, 7}, //10
		{10, 7, 8, 9}//11
};
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
void GetMeshFromSurface(hexa_tree_t* tree, const char* surface_topo, std::vector<double>& coords);
void GetInterceptedElements(hexa_tree_t* mesh, std::vector<double>& coords, std::vector<int>& elements_ids, const char* surface_bathy);
void Apply_material(hexa_tree_t* mesh, std::vector<double>& coords,std::vector<int>& element_ids, const char* surface_bathy);
void ExtrudePMLElements(hexa_tree_t* mesh, std::vector<double>& coords);
void MovingNodes(hexa_tree_t* mesh, std::vector<double>& coords, std::vector<int>& nodes_b_mat, const char* surface);
void MeshOpt(hexa_tree_t* mesh, std::vector<double>& coords, std::vector<int> material_fixed_nodes);
void Adjust_material(hexa_tree_t *mesh);
void UntagleMesh(hexa_tree_t* mesh, std::vector<double>& coords, std::vector<int> material_fixed_nodes);
void hexa_mesh_write_msh(hexa_tree_t* mesh, const char* root_name, std::vector<double> *coords);
void hexa_mesh_write_h5(hexa_tree_t* mesh, const char* root_name, std::vector<double> coords);
int  hexa_mesh_write_vtk(hexa_tree_t* mesh,  const char *filename, std::vector<double> *coords);
void hexa_debug_face_hanging(hexa_tree_t* mesh);
unsigned processors_hash_fn (const void *v, const void *u);
int processors_equal_fn (const void *v1, const void *v2, const void *u);
unsigned edge_hash_fn(const void *v, const void *u);
int edge_equal_fn(const void *v, const void *u, const void *w);
unsigned vertex_hash_id(const void *v, const void *u);
int vertex_equal_id(const void *v, const void *u, const void *w);
void hexa_insert_shared_node(sc_hash_array_t    *shared_nodes, octant_node_t* node, int processor);
unsigned node_hash_fn (const void *v, const void *u);
int node_equal_fn (const void *v1, const void *v2, const void *u);
void hexa_mesh_destroy(hexa_tree_t* mesh);
int node_comp (const void *v, const void *u);
GtsPoint* SegmentTriangleIntersection(GtsSegment * s, GtsTriangle * t);
GtsPoint* LinearMapHex(const double* cord_in_ref, const double* cord_in_x, const double* cord_in_y, const double* cord_in_z);
GtsSurface* SurfaceRead(const char* fname);
gdouble distance(GtsPoint *p, gpointer bounded);
//void communicate_global_ids(hexa_tree_t* mesh);



#endif	/* HEXA_H */

