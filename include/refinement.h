#include <stdio.h>
#include <stdlib.h>
#include <sc.h>
#include <sc_containers.h>
#include <gts.h>
#include <vector>

void IdentifyTemplate(hexa_tree_t* mesh, const std::vector<double>& coords, std::vector<int>& elements_ids);

GtsPoint* SegmentTriangleIntersection(GtsSegment * s, GtsTriangle * t);

void ApplyOctreeTemplate(hexa_tree_t* mesh, std::vector<double>& coords, std::vector<int>& elements_ids);

GtsPoint* LinearMapHex(const double* cord_in_ref, const double* cord_in_x, const double* cord_in_y, const double* cord_in_z);

int AddPoint(hexa_tree_t* mesh, sc_hash_array_t* hash, GtsPoint *p, std::vector<double> &coords);

void IdentifyTemplate(hexa_tree_t* mesh, const std::vector<double>& coords, std::vector<int>& elements_ids);

GtsSurface* SurfaceRead(const char* fname);

gdouble distance(GtsPoint *p, gpointer bounded);

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

typedef struct {
	double coord[3];
	int node_id;
} node_t;
