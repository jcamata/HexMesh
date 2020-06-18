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

typedef struct pillow
{
	int32_t id;
	int32_t a;
	int32_t b;
	int8_t mata;
	int8_t matb;
	bool pa;
	bool pb;
	//	uint64_t elem_a;
	//	uint64_t elem_b;
	int32_t x,y,z;
	int8_t  list_elem;
	int8_t  list_face[8];
	int32_t  elem[8]; // we must check this values
	//int32_t  mat[8]; // we must check this values
	//uint64_t noderef[8];
	//int32_t  octpos[8]; // we must check this values
	int8_t   face[8][3]; // we must check this values
	//int8_t nfaces;
} pillow_t;

typedef struct octant_node
{
	int32_t id;
	int32_t x,y,z;
	int fixed;
	int color;
} octant_node_t;

typedef struct
{
	uint64_t id;
	bool ref;
	uint64_t coord[2];
} octant_edge_t;

typedef struct octant_surf
{
	bool    ext;
} octant_surf_t;

typedef struct {
	double coord[3];
	int node_id;
	bool flag; //maybe can be removed...
} node_t;

typedef struct octant_vertex
{
	uint64_t id;
	int32_t  list_elem;
	int32_t  elem[16]; // we must check this values
} octant_vertex_t;

typedef struct edge
{
	uint64_t id;
	bool ref;
	int32_t  list_elem;
	int32_t  elem[8]; // we must check this values
} edge_t;

typedef struct octant
{
	int32_t    x,y,z;
	int8_t     level;
	int     pad;
	int8_t     pml_id;
	int        n_mat;
	int tem;
	int father;
	octant_node_t nodes[8];
	octant_edge_t edge[12];
	octant_surf_t surf[6];
	int64_t    id;
	bool       boundary;
} octant_t;

typedef struct octree
{
	int64_t id[8];
	bool    cut;
	bool    face[6];
	bool    edge[12];
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

typedef struct shared_octant
{
	int32_t id;
	octant_node_t nodes[8];
	double coord[8][3];
	int listSz;
	int rankList[4];
	int rank;
} shared_octant_t;

typedef struct 
{
	int rank;
	sc_array_t idxs;
} message_t;

typedef struct message_el
{
	int rank;
	sc_array_t idxs;
	sc_array_t nodes;
	sc_array_t coord;
} message_el_t;

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
	sc_array_t      nodes;
	sc_array_t      edges;
	sc_array_t      vertex;
	sc_array_t		oct;
	sc_array_t      outsurf;

	sc_array_t      shared_nodes;
	sc_array_t      shared_edges;

	int64_t *global_id;
	int64_t *global_edge_id;

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

	int neighbors[9];
	comm_map_t comm_map;
	comm_map_t  comm_map_edge;
#ifdef HEXA_DEBUG_
	FILE* fdbg;
#endif
	FILE* profile;
	GeometryData gdata;
	GeometryData tdata;

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

int const EdgeElemOctMap[12][2] = {
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

int const OctNeighbourMap[8][3] = {
		{1, 3, 4}, //El 0
		{2, 0, 5}, //El 1
		{3, 1, 6}, //El 2
		{0, 2, 7}, //El 3
		{5, 7, 0}, //El 4
		{6, 4, 1}, //El 5
		{7, 5, 2}, //El 6
		{4, 6, 3}  //El 7
};
/*
int const EdgeVerticesMap_surf_diagonal[12][2] = {
		{0, 5}, //y-
		{1, 4}, //y-
		{3, 6}, //y+
		{2, 7}, //y+
		{0, 7}, //x-
		{3, 4}, //x-
		{1, 6}, //x+
		{2, 5}, //x+
		{4, 6}, //z-
		{5, 7}, //z-
		{0, 2}, //z+
		{1, 3}  //z+
};

int const EdgeVerticesMap_vol_diagonal[4][2] = {
		{0, 6},
		{1, 7},
		{2, 4},
		{3, 5}
};
 */

int const FaceEdgesMap[6][4] = {
		{4, 11, 7, 3},
		{5, 1, 6, 9},
		{0, 5, 8, 4},
		{2, 6, 10, 7},
		{8, 9, 10, 11},
		{0, 1, 2, 3}
};

int const FaceNormalEdgesMap[6][4] = {
		{0, 8,10, 2},
		{0, 8,10, 2},
		{1, 9,11, 3},
		{1, 9,11, 3},
		{5, 4, 7, 6},
		{5, 4, 7, 6}
};

int const FaceNodesMap[6][4] = {
		{0,4,7,3},
		{1,5,6,2},
		{1,5,4,0},
		{2,6,7,3},
		{1,0,3,2},
		{5,4,7,6}
};

int const FaceNodesMap_inv[6][4] = {
		{1,5,6,2},
		{0,4,7,3},
		{2,6,7,3},
		{1,5,4,0},
		{5,4,7,6},
		{1,0,3,2}
};

int const FaceNodesMapRef[6][4] = {
		{3,0,4,7},
		{2,1,5,6},
		{0,1,5,4},
		{3,2,6,7},
		{7,4,5,6},
		{3,0,1,2}
};

int const FaceNodesMapRef_inv[6][4] = {
		{2,1,5,6},
		{3,0,4,7},
		{3,2,6,7},
		{0,1,5,4},
		{3,0,1,2},
		{7,4,5,6}
};

int const vert_ord[4][3]= {
		{0, 1, 3},
		{1, 2, 0},
		{2, 3, 1},
		{3, 0, 2}
};

int const VertexEdgeMap[8][3] = {
		{0,3,4}, //Vertex 0
		{1,0,5}, //Vertex 1
		{2,1,6}, //Vertex 2
		{3,2,7}, //Vertex 3
		{8,11,4}, //Vertex 4
		{9,8,5}, //Vertex 5
		{10,9,6}, //Vertex 6
		{11,10,7}  //Vertex 7
};

int const VertexSurfMap[8][3] = {
		{0,2,4}, //Vertex 0
		{1,2,4}, //Vertex 1
		{1,3,4}, //Vertex 2
		{0,3,4}, //Vertex 3
		{0,2,5}, //Vertex 4
		{1,2,5}, //Vertex 5
		{1,3,5}, //Vertex 6
		{0,3,5}  //Vertex 7
};

int const VertexVertexMap[8][3] = {
		{1,3,4},
		{0,2,5},
		{1,3,6},
		{0,2,7},
		{5,7,0},
		{4,6,1},
		{5,7,2},
		{4,6,3}
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

int const coord_ref[8][3] = {
		{-1,-1,-1}, // node 0
		{ 1,-1,-1}, // node 1
		{ 1, 1,-1}, // node 2
		{-1, 1,-1}, //node 1
		{-1,-1, 1}, //node 4
		{ 1,-1, 1}, //node 5
		{ 1, 1, 1}, //node 6
		{-1, 1, 1} //node 7
};

void copy_octant(octant_t *orig, octant_t* dest);
void hexa_init(int argc, char* argv[], hexa_tree_t* mesh);
void hexa_finalize(hexa_tree_t* mesh);
void hexa_tree_init(hexa_tree_t* mesh, int max_levels);
void hexa_tree_destroy(hexa_tree_t* mesh);
void hexa_tree_cube(hexa_tree_t* mesh);
int  hexa_tree_write_vtk(hexa_tree_t* mesh,  const char *filename);
void hexa_transition_element(hexa_tree_t* mesh, int i, int j, int k, int step, int level, int ext);
void hexa_processors_interval(hexa_tree_t* mesh);
void hexa_mesh(hexa_tree_t* tree);
void GetMeshFromSurface(hexa_tree_t* tree, const char* surface_topo, std::vector<double>& coords);
void GetInterceptedElements(hexa_tree_t* mesh, std::vector<double>& coords, std::vector<int>& elements_ids, const char* surface_bathy);
void Apply_material(hexa_tree_t* mesh, std::vector<double>& coords,std::vector<int>& element_ids, const char* surface_bathy);
void ExtrudePMLElements(hexa_tree_t* mesh, std::vector<double>& coords);
void MovingNodes(hexa_tree_t* mesh, std::vector<double>& coords, std::vector<int>& nodes_b_mat, const char* surface);
void MeshOpt(hexa_tree_t* mesh, std::vector<double>& coords, std::vector<int> material_fixed_nodes);
void PillowingInterface(hexa_tree_t *mesh, std::vector<double>& coords, std::vector<int>& nodes_b_mat);
void Adjust_material(hexa_tree_t *mesh);
void ExtrudeToOctree(hexa_tree_t* mesh,std::vector<double>& coords);
void hexa_element_init(octant_t *elem);
void hexa_element_copy(octant_t * elem, octant_t * elemCopy);
void MeshOptimization(hexa_tree_t* mesh, std::vector<double>& coords, std::vector<int> material_fixed_nodes);
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
unsigned el_hash_id(const void *v, const void *u);
int el_equal_id(const void *v, const void *u, const void *w);
void hexa_mesh_destroy(hexa_tree_t* mesh);
int node_comp (const void *v, const void *u);
GtsPoint* SegmentTriangleIntersection(GtsSegment * s, GtsTriangle * t);
GtsPoint* LinearMapHex(const double* cord_in_ref, const double* cord_in_x, const double* cord_in_y, const double* cord_in_z);
GtsSurface* SurfaceRead(const char* fname);
gdouble distance(GtsPoint *p, gpointer bounded);
int AddPoint(hexa_tree_t* mesh, sc_hash_array_t* hash_nodes, GtsPoint *p, std::vector<double> &coords, int x, int y, int z);
void CopyPropEl(hexa_tree_t* mesh, int id, octant_t *elem1);
void SetElemPML(hexa_tree_t* tree, octant_t *elem);
void ApplyElement(hexa_tree_t* mesh, std::vector<double>& coords, int id, int iel, int* id_node,
		double* local_ref_x, double* local_ref_y, double* local_ref_z, std::vector<int>& ord, sc_hash_array_t* hash_nodes);

//void communicate_global_ids(hexa_tree_t* mesh);

#endif	/* HEXA_H */

