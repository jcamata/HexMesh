#include <stdio.h>
#include <stdlib.h>
#include <sc.h>
#include <sc_containers.h>
#include <gts.h>
#include <vector>

void IdentifyTemplate(hexa_tree_t* mesh, const std::vector<double>& coords, std::vector<int>& elements_ids);

GtsPoint* SegmentTriangleIntersection(GtsSegment * s, GtsTriangle * t);

void ApplyOctreeTemplate(hexa_tree_t* mesh, std::vector<double>& coords, std::vector<int>& elements_ids);

//unsigned edge_hash_fn(const void *v, const void *u);

//int edge_equal_fn(const void *v, const void *u, const void *w);

void IdentifyTemplate(hexa_tree_t* mesh, const std::vector<double>& coords, std::vector<int>& elements_ids);

GtsSurface* SurfaceRead(const char* fname);

 int const EdgeVerticesMap[12][2] = {
		{0, 1}, // Edge 0
		{1, 2}, // Edge 1
		{2, 3}, //      2
		{3, 0},
		{0, 4},
		{1, 5},
		{2, 6},
		{3, 7},
		{4, 5},
		{5, 6},
		{6, 7},
		{7, 4}
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
