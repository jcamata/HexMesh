#include <gts.h>
#include <glib.h>
#include <vector>
#include <iostream>
using namespace std;
#include <sc.h>
#include <sc_io.h>
#include <sc_containers.h>

#include "hexa.h"
#include "hilbert.h"

typedef struct {
	bitmask_t coord[2];
	int node_id;
} node_in_edge_t;

unsigned edge_hash_fn1(const void *v, const void *u) {
	const node_in_edge_t *q = (const node_in_edge_t*) v;
	uint32_t a, b, c;

	a = (uint32_t) q->coord[0];
	b = (uint32_t) q->coord[1];
	c = (uint32_t) 0;
	sc_hash_mix(a, b, c);
	sc_hash_final(a, b, c);
	return (unsigned) c;
}

int edge_equal_fn1(const void *v, const void *u, const void *w) {
	const node_in_edge_t *e1 = (const node_in_edge_t*) v;
	const node_in_edge_t *e2 = (const node_in_edge_t*) u;

	return (unsigned) ((e1->coord[0] == e2->coord[0]) && (e2->coord[1] == e2->coord[1]));

}

int AddPointOnEdge(int* nodes, sc_hash_array_t* hash, int &npoints, GtsPoint *p, std::vector<double> &coords) {
	size_t position;
	node_in_edge_t *r;
	node_in_edge_t key;
	key.coord[0] = nodes[0];
	key.coord[1] = nodes[1];

	r = (node_in_edge_t*) sc_hash_array_insert_unique(hash, &key, &position);
	if (r != NULL) {
		r->coord[0] = key.coord[0];
		r->coord[1] = key.coord[1];
		r->node_id = npoints;
		npoints++;
		coords.push_back(p->x);
		coords.push_back(p->y);
		coords.push_back(p->z);
		return r->node_id;
	} else {
		r = (node_in_edge_t*) sc_array_index(&hash->a, position);
		return r->node_id;
	}

}

int EdgeVerticesMap[12][2] = {
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

int EdgeVerticesMap_surf_diagonal[12][2] = {
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

int EdgeVerticesMap_vol_diagonal[4][2] = {
		{0, 6},
		{1, 7},
		{2, 4},
		{3, 5}
};

int FaceEdgesMap[6][4] = {
		{4, 11, 7, 3},
		{5, 1, 6, 9},
		{0, 5, 8, 4},
		{2, 6, 10, 7},
		{8, 9, 10, 11},
		{0, 1, 2, 3}
};

int FaceNodesMap[6][4] = {
		{0,4,7,3},
		{1,5,6,2},
		{0,4,5,1},
		{2,6,7,3},
		{4,5,6,7},
		{0,1,2,3}
};

int EdgeEdgeMap[12][4] = {
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



int CheckTemplate(hexa_tree_t* mesh, const std::vector<double>& coords, std::vector<int>& elements_ids, bool flag) {
	bool face_intecepted[6];
	int el_0 = 0;
	int el_1 = 0;
	int el_2 = 0;
	int el_3 = 0;
	int el_not_handle;
	FILE * fdbg;

	fdbg = fopen("intercepted_faces.dbg", "w");

	for (int iel = 0; iel < elements_ids.size(); ++iel) {

		octant_t *elem = (octant_t*) sc_array_index(&mesh->elements, elements_ids[iel]);

		int Edge2GNode[12][2]={0};
		int Edge2GNode_s[12][2]={0};
		int Edge2GNode_v[4][2]={0};

		GtsSegment * segments[12]={0};
		GtsSegment * segments_s[12]={0};
		GtsSegment * segments_v[4]={0};

		GtsPoint * point[12]={NULL};
		GtsPoint * point_s[12]={NULL};
		GtsPoint * point_v[4]={NULL};
		bool edge_list[12]={false};
		bool edge_list_s[12]={false};
		bool edge_list_v[4]={false};
		int ed_cont = 0;

		for (int edge = 0; edge < 12; ++edge) {
			point[edge] = NULL;
			int node1 = elem->nodes[EdgeVerticesMap[edge][0]].id;
			int node2 = elem->nodes[EdgeVerticesMap[edge][1]].id;

			Edge2GNode[edge][0] = node1 <= node2 ? node1 : node2;
			Edge2GNode[edge][1] = node1 >= node2 ? node1 : node2;

			GtsVertex *v1 = gts_vertex_new(gts_vertex_class(), coords[node1 * 3], coords[node1 * 3 + 1], coords[node1 * 3 + 2]);
			GtsVertex *v2 = gts_vertex_new(gts_vertex_class(), coords[node2 * 3], coords[node2 * 3 + 1], coords[node2 * 3 + 2]);

			segments[edge] = gts_segment_new(gts_segment_class(), v1, v2);
			GtsBBox *sb = gts_bbox_segment(gts_bbox_class(), segments[edge]);
			GSList* list = gts_bb_tree_overlap(mesh->gdata.bbt, sb);
			edge_list[edge] = false;
			if (list == NULL) continue;
			while (list) {
				GtsBBox *b = GTS_BBOX(list->data);
				point[edge] = SegmentTriangleIntersection(segments[edge], GTS_TRIANGLE(b->bounded));
				if (point[edge]) {
					edge_list[edge] = true;
					ed_cont++;
					break;
				}
				list = list->next;
			}
		}

		for (int edge = 0; edge < 12; ++edge) {
			point_s[edge] = NULL;
			int node1 = elem->nodes[EdgeVerticesMap_surf_diagonal[edge][0]].id;
			int node2 = elem->nodes[EdgeVerticesMap_surf_diagonal[edge][1]].id;

			Edge2GNode_s[edge][0] = node1 <= node2 ? node1 : node2;
			Edge2GNode_s[edge][1] = node1 >= node2 ? node1 : node2;

			GtsVertex *v1 = gts_vertex_new(gts_vertex_class(), coords[node1 * 3], coords[node1 * 3 + 1], coords[node1 * 3 + 2]);
			GtsVertex *v2 = gts_vertex_new(gts_vertex_class(), coords[node2 * 3], coords[node2 * 3 + 1], coords[node2 * 3 + 2]);

			segments_s[edge] = gts_segment_new(gts_segment_class(), v1, v2);
			GtsBBox *sb = gts_bbox_segment(gts_bbox_class(), segments_s[edge]);
			GSList* list = gts_bb_tree_overlap(mesh->gdata.bbt, sb);
			edge_list_s[edge]=false;
			if (list == NULL) continue;
			while (list) {
				GtsBBox *b = GTS_BBOX(list->data);
				point_s[edge] = SegmentTriangleIntersection(segments_s[edge], GTS_TRIANGLE(b->bounded));
				if (point_s[edge]) {
					edge_list_s[edge]=true;
					break;
				}
				list = list->next;
			}
		}

		for (int edge = 0; edge < 4; ++edge) {
			point_v[edge] = NULL;
			int node1 = elem->nodes[EdgeVerticesMap_vol_diagonal[edge][0]].id;
			int node2 = elem->nodes[EdgeVerticesMap_vol_diagonal[edge][1]].id;

			Edge2GNode_v[edge][0] = node1 <= node2 ? node1 : node2;
			Edge2GNode_v[edge][1] = node1 >= node2 ? node1 : node2;

			GtsVertex *v1 = gts_vertex_new(gts_vertex_class(), coords[node1 * 3], coords[node1 * 3 + 1], coords[node1 * 3 + 2]);
			GtsVertex *v2 = gts_vertex_new(gts_vertex_class(), coords[node2 * 3], coords[node2 * 3 + 1], coords[node2 * 3 + 2]);

			segments_v[edge] = gts_segment_new(gts_segment_class(), v1, v2);
			GtsBBox *sb = gts_bbox_segment(gts_bbox_class(), segments_v[edge]);
			GSList* list = gts_bb_tree_overlap(mesh->gdata.bbt, sb);
			edge_list_v[edge]=false;
			if (list == NULL) continue;
			while (list) {
				GtsBBox *b = GTS_BBOX(list->data);
				point_v[edge] = SegmentTriangleIntersection(segments_v[edge], GTS_TRIANGLE(b->bounded));
				if (point_v[edge]) {
					edge_list_v[edge]=true;
					break;
				}
				list = list->next;
			}
		}

		//check parallel faces
		for (int face = 0; face < 6; ++face) {
			face_intecepted[face] = false;
			for (int fe = 0; fe < 4; ++fe) {
				int edge = FaceEdgesMap[face][fe];
				if (point[edge] != NULL) {
					face_intecepted[face] = true;
					break;
				}
			}
			if (face_intecepted[face]) continue;
		}

		fprintf(fdbg, "Faces:  Elem. %d: %d %d %d %d %d %d\n", elements_ids[iel], face_intecepted[0],
				face_intecepted[1], face_intecepted[2],
				face_intecepted[3], face_intecepted[4], face_intecepted[5]);
		fprintf(fdbg, "Edges:  Elem. %d: %d %d %d %d %d %d %d %d %d %d %d %d %d\n", elements_ids[iel], edge_list[0],
				edge_list[1], edge_list[2], edge_list[3],
				edge_list[4], edge_list[5], edge_list[6],
				edge_list[7], edge_list[8], edge_list[9],
				edge_list[10], edge_list[11], ed_cont);

		int n_parallel_faces = (face_intecepted[0] && face_intecepted[1]) +
				(face_intecepted[2] && face_intecepted[3]) +
				(face_intecepted[4] && face_intecepted[5]);

		// element verification
		if (n_parallel_faces == 2 && ed_cont == 4) {
			// Check template 1.
			if (    (!edge_list[0]) && (edge_list[1]) && (!edge_list[2]) && (edge_list[3]) &&
					(!edge_list[4]) && (!edge_list[5]) && (!edge_list[6]) && (!edge_list[7]) &&
					(!edge_list[8]) && (edge_list[9]) && (!edge_list[10]) && (edge_list[11]) ) {

				elem->cut = 1;

			} else if (     (edge_list[0]) && (!edge_list[1]) && (edge_list[2]) && (!edge_list[3]) &&
					(!edge_list[4]) && (!edge_list[5]) && (!edge_list[6]) && (!edge_list[7]) &&
					(edge_list[8]) && (!edge_list[9]) && (edge_list[10]) && (!edge_list[11]) ) {

				elem->cut = 1;

			} else if (     (!edge_list[0]) && (!edge_list[1]) && (!edge_list[2]) && (!edge_list[3]) &&
					(edge_list[4]) && (edge_list[5]) && (edge_list[6]) && (edge_list[7]) &&
					(!edge_list[8]) && (!edge_list[9]) && (!edge_list[10]) && (!edge_list[11]) ) {

				elem->cut = 1;

			} else {
				elem->cut = -10;
				fprintf(fdbg, "Warning: Bad Element: %d, case 1\n",elements_ids[iel]);
				fprintf(fdbg,"Edges:  Elem. %d: %d %d %d %d %d %d %d %d %d %d %d %d %d\n", elements_ids[iel], edge_list[0],
						edge_list[1], edge_list[2], edge_list[3],
						edge_list[4], edge_list[5], edge_list[6],
						edge_list[7], edge_list[8], edge_list[9],
						edge_list[10], edge_list[11], ed_cont);
				el_not_handle++;
				continue;
			}
		}

		if (n_parallel_faces == 1 && ed_cont==4) {
			// Check template 2.
			//Edge 0
			if ( (!edge_list[0]) && (edge_list[1]) && (!edge_list[2]) && (edge_list[3]) &&
					(edge_list[4]) && (edge_list[5]) && (!edge_list[6]) && (!edge_list[7]) &&
					(!edge_list[8]) && (!edge_list[9]) && (!edge_list[10]) && (!edge_list[11]) &&
					(edge_list_s[4]) && (edge_list_s[6]) ) {

				elem->cut = 2;

			}//Edge 1
			else if ((edge_list[0]) && (!edge_list[1]) && (edge_list[2]) && (!edge_list[3]) &&
					(!edge_list[4]) && (edge_list[5]) && (edge_list[6]) && (!edge_list[7]) &&
					(!edge_list[8]) && (!edge_list[9]) && (!edge_list[10]) && (!edge_list[11])&&
					(edge_list_s[1]) && (edge_list_s[3]) ) {

				elem->cut = 2;

			}//Edge 2
			else if ((!edge_list[0]) && (edge_list[1]) && (!edge_list[2]) && (edge_list[3]) &&
					(!edge_list[4]) && (!edge_list[5]) && (edge_list[6]) && (edge_list[7]) &&
					(!edge_list[8]) && (!edge_list[9]) && (!edge_list[10]) && (!edge_list[11])&&
					(edge_list_s[5]) && (edge_list_s[7]) ) {

				elem->cut = 2;

			}//Edge 3
			else if ((edge_list[0]) && (!edge_list[1]) && (edge_list[2]) && (!edge_list[3]) &&
					(edge_list[4]) && (!edge_list[5]) && (!edge_list[6]) && (edge_list[7]) &&
					(!edge_list[8]) && (!edge_list[9]) && (!edge_list[10]) && (!edge_list[11])&&
					(edge_list_s[0]) && (edge_list_s[2])) {

				elem->cut = 2;

			}//Edge 4
			else if ((edge_list[0]) && (!edge_list[1]) && (!edge_list[2]) && (edge_list[3]) &&
					(!edge_list[4]) && (!edge_list[5]) && (!edge_list[6]) && (!edge_list[7]) &&
					(edge_list[8]) && (!edge_list[9]) && (!edge_list[10]) && (edge_list[11])&&
					(edge_list_s[8]) && (edge_list_s[10])) {

				elem->cut = 2;

			}//Edge 5
			else if ((edge_list[0]) && (edge_list[1]) && (!edge_list[2]) && (!edge_list[3]) &&
					(!edge_list[4]) && (!edge_list[5]) && (!edge_list[6]) && (!edge_list[7]) &&
					(edge_list[8]) && (edge_list[9]) && (!edge_list[10]) && (!edge_list[11])&&
					(edge_list_s[9]) && (edge_list_s[11])) {

				elem->cut = 2;

			}//Edge 6
			else if ((!edge_list[0]) && (edge_list[1]) && (edge_list[2]) && (!edge_list[3]) &&
					(!edge_list[4]) && (!edge_list[5]) && (!edge_list[6]) && (!edge_list[7]) &&
					(!edge_list[8]) && (edge_list[9]) && (edge_list[10]) && (!edge_list[11])&&
					(edge_list_s[8]) && (edge_list_s[10])) {

				elem->cut = 2;

			}//Edge 7
			else if ((!edge_list[0]) && (!edge_list[1]) && (edge_list[2]) && (edge_list[3]) &&
					(!edge_list[4]) && (!edge_list[5]) && (!edge_list[6]) && (!edge_list[7]) &&
					(!edge_list[8]) && (!edge_list[9]) && (edge_list[10]) && (edge_list[11])&&
					(edge_list_s[9]) && (edge_list_s[11])) {

				elem->cut = 2;

			}//Edge 8
			else if ((!edge_list[0]) && (!edge_list[1]) && (!edge_list[2]) && (!edge_list[3]) &&
					(edge_list[4]) && (edge_list[5]) && (!edge_list[6]) && (!edge_list[7]) &&
					(!edge_list[8]) && (edge_list[9]) && (!edge_list[10]) && (edge_list[11])&&
					(edge_list_s[5]) && (edge_list_s[7])) {

				elem->cut = 2;

			}//Edge 9
			else if ((!edge_list[0]) && (!edge_list[1]) && (!edge_list[2]) && (!edge_list[3]) &&
					(!edge_list[4]) && (edge_list[5]) && (edge_list[6]) && (!edge_list[7]) &&
					(edge_list[8]) && (!edge_list[9]) && (edge_list[10]) && (!edge_list[11])&&
					(edge_list_s[0]) && (edge_list_s[2])) {

				elem->cut = 2;

			}//Edge 10
			else if ((!edge_list[0]) && (!edge_list[1]) && (!edge_list[2]) && (!edge_list[3]) &&
					(!edge_list[4]) && (!edge_list[5]) && (edge_list[6]) && (edge_list[7]) &&
					(!edge_list[8]) && (edge_list[9]) && (!edge_list[10]) && (edge_list[11])&&
					(edge_list_s[4]) && (edge_list_s[6])) {

				elem->cut = 2;

			}//Edge 11
			else if ((!edge_list[0]) && (!edge_list[1]) && (!edge_list[2]) && (!edge_list[3]) &&
					(edge_list[4]) && (!edge_list[5]) && (!edge_list[6]) && (edge_list[7]) &&
					(edge_list[8]) && (!edge_list[9]) && (edge_list[10]) && (!edge_list[11])&&
					(edge_list_s[1]) && (edge_list_s[3])) {

				elem->cut = 2;

			} else {
				elem->cut = -10;
				fprintf(fdbg, "Warning: Bad Element: %d, case 2\n",elements_ids[iel]);
				fprintf(fdbg, "Edges:  Elem. %d: %d %d %d %d %d %d %d %d %d %d %d %d %d\n", elements_ids[iel], edge_list[0],
						edge_list[1], edge_list[2], edge_list[3],
						edge_list[4], edge_list[5], edge_list[6],
						edge_list[7], edge_list[8], edge_list[9],
						edge_list[10], edge_list[11], ed_cont);
				el_not_handle++;
				continue;
			}
		}

		if (n_parallel_faces == 0 && ed_cont == 3) {
			// Check template 3.
			//Corner 0
			if ((edge_list[0]) && (!edge_list[1]) && (!edge_list[2]) && (edge_list[3]) &&
					(edge_list[4]) && (!edge_list[5]) && (!edge_list[6]) && (!edge_list[7]) &&
					(!edge_list[8]) && (!edge_list[9]) && (!edge_list[10]) && (!edge_list[11]) &&
					(edge_list_s[4]) && (edge_list_s[0]) && (edge_list_s[10]) &&
					(edge_list_v[0]) ) {

				elem->cut = 3;

			}//Corner 1
			else if ((edge_list[0]) && (edge_list[1]) && (!edge_list[2]) && (!edge_list[3]) &&
					(!edge_list[4]) && (edge_list[5]) && (!edge_list[6]) && (!edge_list[7]) &&
					(!edge_list[8]) && (!edge_list[9]) && (!edge_list[10]) && (!edge_list[11])&&
					(edge_list_s[1]) && (edge_list_s[6]) && (edge_list_s[11])  &&
					(edge_list_v[1]) ) {

				elem->cut = 3;

			} //Corner 2
			else if ((!edge_list[0]) && (edge_list[1]) && (edge_list[2]) && (!edge_list[3]) &&
					(!edge_list[4]) && (!edge_list[5]) && (edge_list[6]) && (!edge_list[7]) &&
					(!edge_list[8]) && (!edge_list[9]) && (!edge_list[10]) && (!edge_list[11]) &&
					(edge_list_s[3]) && (edge_list_s[7]) && (edge_list_s[10])  &&
					(edge_list_v[2]) ) {

				elem->cut = 3;

			} //Corner 3
			else if ((!edge_list[0]) && (!edge_list[1]) && (edge_list[2]) && (edge_list[3]) &&
					(!edge_list[4]) && (!edge_list[5]) && (!edge_list[6]) && (edge_list[7]) &&
					(!edge_list[8]) && (!edge_list[9]) && (!edge_list[10]) && (!edge_list[11]) &&
					(edge_list_s[5]) && (edge_list_s[2]) && (edge_list_s[11])  &&
					(edge_list_v[3]) ) {

				elem->cut = 3;

			} //Corner 4
			else if ((!edge_list[0]) && (!edge_list[1]) && (!edge_list[2]) && (!edge_list[3]) &&
					(edge_list[4]) && (!edge_list[5]) && (!edge_list[6]) && (!edge_list[7]) &&
					(edge_list[8]) && (!edge_list[9]) && (!edge_list[10]) && (edge_list[11]) &&
					(edge_list_s[1]) && (edge_list_s[5]) && (edge_list_s[8])  &&
					(edge_list_v[2]) ) {

				elem->cut = 3;

			} //Corner 5
			else if ((!edge_list[0]) && (!edge_list[1]) && (!edge_list[2]) && (!edge_list[3]) &&
					(!edge_list[4]) && (edge_list[5]) && (!edge_list[6]) && (!edge_list[7]) &&
					(edge_list[8]) && (edge_list[9]) && (!edge_list[10]) && (!edge_list[11]) &&
					(edge_list_s[7]) && (edge_list_s[0]) && (edge_list_s[9])  &&
					(edge_list_v[3]) ) {

				elem->cut = 3;

			} //Corner 6
			else if ((!edge_list[0]) && (!edge_list[1]) && (!edge_list[2]) && (!edge_list[3]) &&
					(!edge_list[4]) && (!edge_list[5]) && (edge_list[6]) && (!edge_list[7]) &&
					(!edge_list[8]) && (edge_list[9]) && (edge_list[10]) && (!edge_list[11]) &&
					(edge_list_s[2]) && (edge_list_s[6]) && (edge_list_s[8])  &&
					(edge_list_v[0]) ) {

				elem->cut = 3;

			} //Corner 7
			else if ((!edge_list[0]) && (!edge_list[1]) && (!edge_list[2]) && (!edge_list[3]) &&
					(!edge_list[4]) && (!edge_list[5]) && (!edge_list[6]) && (edge_list[7]) &&
					(!edge_list[8]) && (!edge_list[9]) && (edge_list[10]) && (edge_list[11]) &&
					(edge_list_s[4]) && (edge_list_s[3]) && (edge_list_s[9])  &&
					(edge_list_v[1]) ) {
				elem->cut = 3;
			}else {
				elem->cut = -10;
				fprintf(fdbg, "Warning: Bad Element: %d, case 3\n",elements_ids[iel]);
				fprintf(fdbg, "Edges:  Elem. %d: %d %d %d %d %d %d %d %d %d %d %d %d %d\n", elements_ids[iel], edge_list[0],
						edge_list[1], edge_list[2], edge_list[3],
						edge_list[4], edge_list[5], edge_list[6],
						edge_list[7], edge_list[8], edge_list[9],
						edge_list[10], edge_list[11], ed_cont);
				el_not_handle++;
				continue;
			}
		}



		//clean points
		for (int edge = 0; edge < 4; edge++) {
			if (point_v[edge]) gts_object_destroy(GTS_OBJECT(point_v[edge]));
			point_v[edge] = NULL;
		}

		for (int edge = 0; edge < 12; edge++) {
			if (point_s[edge]) gts_object_destroy(GTS_OBJECT(point_s[edge]));
			point_s[edge] = NULL;
		}

		for (int edge = 0; edge < 12; edge++) {
			if (point[edge]) gts_object_destroy(GTS_OBJECT(point[edge]));
			point[edge] = NULL;
		}

		if(elem->cut==0){el_0++;}
		if(elem->cut==1){el_1++;}
		if(elem->cut==2){el_2++;}
		if(elem->cut==3){el_3++;}
		//element not handle
		if(elem->cut==-1){
			elem->cut = -10;
			fprintf(fdbg, "Warning: Bad Element: %d, case -1\n",elements_ids[iel]);
			fprintf(fdbg, "Edges:  Elem. %d: %d %d %d %d %d %d %d %d %d %d %d %d %d\n", elements_ids[iel], edge_list[0],
					edge_list[1], edge_list[2], edge_list[3],
					edge_list[4], edge_list[5], edge_list[6],
					edge_list[7], edge_list[8], edge_list[9],
					edge_list[10], edge_list[11], ed_cont);
			el_not_handle++;
		}
	}

	if(flag){
		printf("el_not_handle: %d\n",el_not_handle);
		printf("case_0: %d\n",el_0);
		printf("case_1: %d\n",el_1);
		printf("case_2: %d\n",el_2);
		printf("case_3: %d\n",el_3);
	}
	fclose(fdbg);
	return el_not_handle;
}

void CutTemplate(hexa_tree_t* mesh, std::vector<double>& coords, std::vector<int>& elements_ids) {
	bool clamped = true;
	bool face_intecepted[6];
	int conn_p[4];
	int original_conn[8];
	FILE * fdbg;

	int Edge2GNode[12][2]={0};
	int Edge2GNode_s[12][2]={0};
	int Edge2GNode_v[4][2]={0};

	GtsSegment * segments[12]={0};
	GtsSegment * segments_s[12]={0};
	GtsSegment * segments_v[4]={0};

	GtsPoint * point[12]={NULL};
	GtsPoint * point_s[12]={NULL};
	GtsPoint * point_v[4]={NULL};
	bool edge_list[12]={false};
	bool edge_list_s[12]={false};
	bool edge_list_v[4]={false};
	int ed_cont = 0;

	fdbg = fopen("intercepted_faces.dbg", "w");

	sc_hash_array_t* hash_nodes = sc_hash_array_new(sizeof (node_in_edge_t), edge_hash_fn1, edge_equal_fn1, &clamped);

	for (int iel = 0; iel < elements_ids.size(); ++iel) {

		octant_t *elem = (octant_t*) sc_array_index(&mesh->elements, elements_ids[iel]);
		for (int i = 0; i < 8; i++) original_conn[i] = elem->nodes[i].id;

		for (int edge = 0; edge < 12; ++edge) {
			point[edge] = NULL;
			int node1 = elem->nodes[EdgeVerticesMap[edge][0]].id;
			int node2 = elem->nodes[EdgeVerticesMap[edge][1]].id;

			Edge2GNode[edge][0] = node1 <= node2 ? node1 : node2;
			Edge2GNode[edge][1] = node1 >= node2 ? node1 : node2;

			GtsVertex *v1 = gts_vertex_new(gts_vertex_class(), coords[node1 * 3], coords[node1 * 3 + 1], coords[node1 * 3 + 2]);
			GtsVertex *v2 = gts_vertex_new(gts_vertex_class(), coords[node2 * 3], coords[node2 * 3 + 1], coords[node2 * 3 + 2]);

			segments[edge] = gts_segment_new(gts_segment_class(), v1, v2);
			GtsBBox *sb = gts_bbox_segment(gts_bbox_class(), segments[edge]);
			GSList* list = gts_bb_tree_overlap(mesh->gdata.bbt, sb);
			edge_list[edge] = false;
			if (list == NULL) continue;
			while (list) {
				GtsBBox *b = GTS_BBOX(list->data);
				point[edge] = SegmentTriangleIntersection(segments[edge], GTS_TRIANGLE(b->bounded));
				if (point[edge]) {
					edge_list[edge] = true;
					ed_cont++;
					break;
				}
				list = list->next;
			}
		}

		//check the diagonals in the surface and find the intersections
		for (int edge = 0; edge < 12; ++edge) {
			point_s[edge] = NULL;
			int node1 = elem->nodes[EdgeVerticesMap_surf_diagonal[edge][0]].id;
			int node2 = elem->nodes[EdgeVerticesMap_surf_diagonal[edge][1]].id;

			Edge2GNode_s[edge][0] = node1 <= node2 ? node1 : node2;
			Edge2GNode_s[edge][1] = node1 >= node2 ? node1 : node2;

			GtsVertex *v1 = gts_vertex_new(gts_vertex_class(), coords[node1 * 3], coords[node1 * 3 + 1], coords[node1 * 3 + 2]);
			GtsVertex *v2 = gts_vertex_new(gts_vertex_class(), coords[node2 * 3], coords[node2 * 3 + 1], coords[node2 * 3 + 2]);

			segments_s[edge] = gts_segment_new(gts_segment_class(), v1, v2);
			GtsBBox *sb = gts_bbox_segment(gts_bbox_class(), segments_s[edge]);
			GSList* list = gts_bb_tree_overlap(mesh->gdata.bbt, sb);
			edge_list_s[edge]=false;
			if (list == NULL) continue;
			while (list) {
				GtsBBox *b = GTS_BBOX(list->data);
				point_s[edge] = SegmentTriangleIntersection(segments_s[edge], GTS_TRIANGLE(b->bounded));
				if (point_s[edge]) {
					edge_list_s[edge]=true;
					break;
				}
				list = list->next;
			}
		}

		//check the diagonals in the volume and find the intersections
		for (int edge = 0; edge < 4; ++edge) {
			point_v[edge] = NULL;
			int node1 = elem->nodes[EdgeVerticesMap_vol_diagonal[edge][0]].id;
			int node2 = elem->nodes[EdgeVerticesMap_vol_diagonal[edge][1]].id;

			Edge2GNode_v[edge][0] = node1 <= node2 ? node1 : node2;
			Edge2GNode_v[edge][1] = node1 >= node2 ? node1 : node2;

			GtsVertex *v1 = gts_vertex_new(gts_vertex_class(), coords[node1 * 3], coords[node1 * 3 + 1], coords[node1 * 3 + 2]);
			GtsVertex *v2 = gts_vertex_new(gts_vertex_class(), coords[node2 * 3], coords[node2 * 3 + 1], coords[node2 * 3 + 2]);

			segments_v[edge] = gts_segment_new(gts_segment_class(), v1, v2);
			GtsBBox *sb = gts_bbox_segment(gts_bbox_class(), segments_v[edge]);
			GSList* list = gts_bb_tree_overlap(mesh->gdata.bbt, sb);
			edge_list_v[edge]=false;
			if (list == NULL) continue;
			while (list) {
				GtsBBox *b = GTS_BBOX(list->data);
				point_v[edge] = SegmentTriangleIntersection(segments_v[edge], GTS_TRIANGLE(b->bounded));
				if (point_v[edge]) {
					edge_list_v[edge]=true;
					break;
				}
				list = list->next;
			}
		}

		//check parallel faces
		for (int face = 0; face < 6; ++face) {
			face_intecepted[face] = false;
			for (int fe = 0; fe < 4; ++fe) {
				int edge = FaceEdgesMap[face][fe];
				if (point[edge] != NULL) {
					face_intecepted[face] = true;
					break;
				}
			}
			if (face_intecepted[face]) continue;
		}

		fprintf(fdbg, "Faces:  Elem. %d: %d %d %d %d %d %d\n", elements_ids[iel], face_intecepted[0],
				face_intecepted[1], face_intecepted[2],
				face_intecepted[3], face_intecepted[4], face_intecepted[5]);
		fprintf(fdbg, "Edges:  Elem. %d: %d %d %d %d %d %d %d %d %d %d %d %d %d\n", elements_ids[iel], edge_list[0],
				edge_list[1], edge_list[2], edge_list[3],
				edge_list[4], edge_list[5], edge_list[6],
				edge_list[7], edge_list[8], edge_list[9],
				edge_list[10], edge_list[11], ed_cont);

		int n_parallel_faces = (face_intecepted[0] && face_intecepted[1]) +
				(face_intecepted[2] && face_intecepted[3]) +
				(face_intecepted[4] && face_intecepted[5]);

		// verificação dos elementos
		if(elem->cut==1){

			// Check template 1.
			int cut_edge[4];
			int connec_order_in[8];
			int connec_order_out1[8];
			int connec_order_out2[8];

			if (    (!edge_list[0]) && (edge_list[1]) && (!edge_list[2]) && (edge_list[3]) &&
					(!edge_list[4]) && (!edge_list[5]) && (!edge_list[6]) && (!edge_list[7]) &&
					(!edge_list[8]) && (edge_list[9]) && (!edge_list[10]) && (edge_list[11]) ) {

				cut_edge[0] = 3;
				cut_edge[1] = 1;
				cut_edge[2] = 11;
				cut_edge[3] = 9;

				connec_order_in[0] = 0;
				connec_order_in[1] = 1;
				connec_order_in[2] = 4;
				connec_order_in[3] = 5;
				connec_order_in[4] = 2;
				connec_order_in[5] = 3;
				connec_order_in[6] = 6;
				connec_order_in[7] = 7;

				connec_order_out1[0] = 0;
				connec_order_out1[1] = 1;
				connec_order_out1[2] = 4;
				connec_order_out1[3] = 5;
				connec_order_out1[4] = 3;
				connec_order_out1[5] = 2;
				connec_order_out1[6] = 7;
				connec_order_out1[7] = 6;

				connec_order_out2[0] = 0;
				connec_order_out2[1] = 1;
				connec_order_out2[2] = 4;
				connec_order_out2[3] = 5;
				connec_order_out2[4] = 2;
				connec_order_out2[5] = 3;
				connec_order_out2[6] = 6;
				connec_order_out2[7] = 7;

				elem->cut = 1;

			} else if (     (edge_list[0]) && (!edge_list[1]) && (edge_list[2]) && (!edge_list[3]) &&
					(!edge_list[4]) && (!edge_list[5]) && (!edge_list[6]) && (!edge_list[7]) &&
					(edge_list[8]) && (!edge_list[9]) && (edge_list[10]) && (!edge_list[11]) ) {

				cut_edge[0] = 0;
				cut_edge[1] = 2;
				cut_edge[2] = 8;
				cut_edge[3] = 10;

				connec_order_in[0] = 0;
				connec_order_in[1] = 3;
				connec_order_in[2] = 4;
				connec_order_in[3] = 7;
				connec_order_in[4] = 1;
				connec_order_in[5] = 2;
				connec_order_in[6] = 5;
				connec_order_in[7] = 6;

				connec_order_out1[0] = 0;
				connec_order_out1[1] = 3;
				connec_order_out1[2] = 4;
				connec_order_out1[3] = 7;
				connec_order_out1[4] = 1;
				connec_order_out1[5] = 2;
				connec_order_out1[6] = 5;
				connec_order_out1[7] = 6;

				connec_order_out2[0] = 0;
				connec_order_out2[1] = 3;
				connec_order_out2[2] = 4;
				connec_order_out2[3] = 7;
				connec_order_out2[4] = 1;
				connec_order_out2[5] = 2;
				connec_order_out2[6] = 5;
				connec_order_out2[7] = 6;

				elem->cut = 1;

			} else if (     (!edge_list[0]) && (!edge_list[1]) && (!edge_list[2]) && (!edge_list[3]) &&
					(edge_list[4]) && (edge_list[5]) && (edge_list[6]) && (edge_list[7]) &&
					(!edge_list[8]) && (!edge_list[9]) && (!edge_list[10]) && (!edge_list[11]) ) {

				cut_edge[0] = 4;
				cut_edge[1] = 5;
				cut_edge[2] = 6;
				cut_edge[3] = 7;

				connec_order_in[0] = 0;
				connec_order_in[1] = 1;
				connec_order_in[2] = 2;
				connec_order_in[3] = 3;
				connec_order_in[4] = 4;
				connec_order_in[5] = 5;
				connec_order_in[6] = 6;
				connec_order_in[7] = 7;

				connec_order_out1[0] = 0;
				connec_order_out1[1] = 1;
				connec_order_out1[2] = 2;
				connec_order_out1[3] = 3;
				connec_order_out1[4] = 4;
				connec_order_out1[5] = 5;
				connec_order_out1[6] = 6;
				connec_order_out1[7] = 7;

				connec_order_out2[0] = 0;
				connec_order_out2[1] = 1;
				connec_order_out2[2] = 2;
				connec_order_out2[3] = 3;
				connec_order_out2[4] = 4;
				connec_order_out2[5] = 5;
				connec_order_out2[6] = 6;
				connec_order_out2[7] = 7;

				elem->cut = 1;

			} else {
				elem->cut = -1;
				printf("Warning: Bad Element: %d, case 1\n",elements_ids[iel]);
				printf("Edges:  Elem. %d: %d %d %d %d %d %d %d %d %d %d %d %d %d\n", elements_ids[iel], edge_list[0],
						edge_list[1], edge_list[2], edge_list[3],
						edge_list[4], edge_list[5], edge_list[6],
						edge_list[7], edge_list[8], edge_list[9],
						edge_list[10], edge_list[11], ed_cont);

				continue;
			}

			GtsPoint *p0 = point[cut_edge[0]];
			GtsPoint *p1 = point[cut_edge[1]];
			GtsPoint *p2 = point[cut_edge[2]];
			GtsPoint *p3 = point[cut_edge[3]];

			//g_assert(p0 != NULL);
			//g_assert(p1 != NULL);
			//g_assert(p2 != NULL);
			//g_assert(p3 != NULL);

			if(p0==NULL) exit(1);
			if(p1==NULL) exit(1);
			if(p2==NULL) exit(1);
			if(p3==NULL) exit(1);


			conn_p[0] = AddPointOnEdge(Edge2GNode[cut_edge[0]], hash_nodes, mesh->local_n_nodes, p0, coords);
			conn_p[1] = AddPointOnEdge(Edge2GNode[cut_edge[1]], hash_nodes, mesh->local_n_nodes, p1, coords);
			conn_p[2] = AddPointOnEdge(Edge2GNode[cut_edge[2]], hash_nodes, mesh->local_n_nodes, p2, coords);
			conn_p[3] = AddPointOnEdge(Edge2GNode[cut_edge[3]], hash_nodes, mesh->local_n_nodes, p3, coords);

			octant_t *elem1 = (octant_t*) sc_array_index(&mesh->elements, elements_ids[iel]);

			elem1->nodes[connec_order_out1[0]].id = original_conn[connec_order_in[0]];
			elem1->nodes[connec_order_out1[1]].id = original_conn[connec_order_in[1]];
			elem1->nodes[connec_order_out1[2]].id = original_conn[connec_order_in[2]];
			elem1->nodes[connec_order_out1[3]].id = original_conn[connec_order_in[3]];

			elem1->nodes[connec_order_out1[4]].id = conn_p[0];
			elem1->nodes[connec_order_out1[5]].id = conn_p[1];
			elem1->nodes[connec_order_out1[6]].id = conn_p[2];
			elem1->nodes[connec_order_out1[7]].id = conn_p[3];

			octant_t* elem2 = (octant_t*) sc_array_push(&mesh->elements);

			elem2->nodes[connec_order_out2[0]].id = conn_p[0];
			elem2->nodes[connec_order_out2[1]].id = conn_p[1];
			elem2->nodes[connec_order_out2[2]].id = conn_p[2];
			elem2->nodes[connec_order_out2[3]].id = conn_p[3];

			elem2->nodes[connec_order_out2[4]].id = original_conn[connec_order_in[4]];
			elem2->nodes[connec_order_out2[5]].id = original_conn[connec_order_in[5]];
			elem2->nodes[connec_order_out2[6]].id = original_conn[connec_order_in[6]];
			elem2->nodes[connec_order_out2[7]].id = original_conn[connec_order_in[7]];

			elem2->cut = elem->cut;
			elem2->level = elem->level;



		}else if(elem->cut==2){

			// Check template 2.
			int conn_t2[6];
			int cut_edge[4];
			int cut_edge_s[2];
			int connec_order_in1[8];
			int connec_order_out1[8];
			int connec_order_in2[8];
			int connec_order_out2[8];
			int connec_order_in3[8];
			int connec_order_out3[8];

			//Edge 0
			if ( (!edge_list[0]) && (edge_list[1]) && (!edge_list[2]) && (edge_list[3]) &&
					(edge_list[4]) && (edge_list[5]) && (!edge_list[6]) && (!edge_list[7]) &&
					(!edge_list[8]) && (!edge_list[9]) && (!edge_list[10]) && (!edge_list[11]) &&
					(edge_list_s[4]) && (edge_list_s[6]) ) {
				//printf("Entrou na 0!\n");

				cut_edge[0] = 3;
				cut_edge[1] = 1;
				cut_edge[2] = 4;
				cut_edge[3] = 5;

				cut_edge_s[0] = 4;
				cut_edge_s[1] = 6;

				connec_order_in1[0] = 2;
				connec_order_in1[1] = 3;
				connec_order_in1[2] = 6;
				connec_order_in1[3] = 7;
				connec_order_in1[4] = 0;
				connec_order_in1[5] = 1;
				connec_order_in1[6] = 4;
				connec_order_in1[7] = 5;

				connec_order_out1[0] = 2;
				connec_order_out1[1] = 3;
				connec_order_out1[2] = 6;
				connec_order_out1[3] = 7;
				connec_order_out1[4] = 0;
				connec_order_out1[5] = 1;
				connec_order_out1[6] = 4;
				connec_order_out1[7] = 5;

				connec_order_in2[0] = 0;
				connec_order_in2[1] = 1;
				connec_order_in2[2] = 2;
				connec_order_in2[3] = 3;
				connec_order_in2[4] = 4;
				connec_order_in2[5] = 5;
				connec_order_in2[6] = 6;
				connec_order_in2[7] = 7;

				connec_order_out2[0] = 2;
				connec_order_out2[1] = 3;
				connec_order_out2[2] = 5;
				connec_order_out2[3] = 4;
				connec_order_out2[4] = 4;
				connec_order_out2[5] = 5;
				connec_order_out2[6] = 6;
				connec_order_out2[7] = 7;

				connec_order_in3[0] = 0;
				connec_order_in3[1] = 1;
				connec_order_in3[2] = 2;
				connec_order_in3[3] = 3;
				connec_order_in3[4] = 4;
				connec_order_in3[5] = 5;
				connec_order_in3[6] = 6;
				connec_order_in3[7] = 7;

				connec_order_out3[0] = 0;
				connec_order_out3[1] = 1;
				connec_order_out3[2] = 1;
				connec_order_out3[3] = 0;
				connec_order_out3[4] = 2;
				connec_order_out3[5] = 3;
				connec_order_out3[6] = 5;
				connec_order_out3[7] = 4;

				elem->cut = 2;

			}//Edge 1
			else if ((edge_list[0]) && (!edge_list[1]) && (edge_list[2]) && (!edge_list[3]) &&
					(!edge_list[4]) && (edge_list[5]) && (edge_list[6]) && (!edge_list[7]) &&
					(!edge_list[8]) && (!edge_list[9]) && (!edge_list[10]) && (!edge_list[11])&&
					(edge_list_s[1]) && (edge_list_s[3]) ) {
				//printf("Entrou na 1!\n");

				cut_edge[0] = 0;
				cut_edge[1] = 2;
				cut_edge[2] = 5;
				cut_edge[3] = 6;

				cut_edge_s[0] = 1;
				cut_edge_s[1] = 3;

				connec_order_in1[0] = 0;
				connec_order_in1[1] = 3;
				connec_order_in1[2] = 4;
				connec_order_in1[3] = 7;
				connec_order_in1[4] = 1;
				connec_order_in1[5] = 2;
				connec_order_in1[6] = 5;
				connec_order_in1[7] = 6;

				connec_order_out1[0] = 0;
				connec_order_out1[1] = 3;
				connec_order_out1[2] = 4;
				connec_order_out1[3] = 7;
				connec_order_out1[4] = 0;
				connec_order_out1[5] = 1;
				connec_order_out1[6] = 4;
				connec_order_out1[7] = 5;

				connec_order_in2[0] = 0;
				connec_order_in2[1] = 1;
				connec_order_in2[2] = 2;
				connec_order_in2[3] = 3;
				connec_order_in2[4] = 4;
				connec_order_in2[5] = 5;
				connec_order_in2[6] = 6;
				connec_order_in2[7] = 7;

				connec_order_out2[0] = 4;
				connec_order_out2[1] = 2;
				connec_order_out2[2] = 3;
				connec_order_out2[3] = 5;
				connec_order_out2[4] = 4;
				connec_order_out2[5] = 5;
				connec_order_out2[6] = 6;
				connec_order_out2[7] = 7;

				connec_order_in3[0] = 1;
				connec_order_in3[1] = 2;
				connec_order_in3[2] = 0;
				connec_order_in3[3] = 3;
				connec_order_in3[4] = 4;
				connec_order_in3[5] = 5;
				connec_order_in3[6] = 6;
				connec_order_in3[7] = 7;

				connec_order_out3[0] = 1;
				connec_order_out3[1] = 2;
				connec_order_out3[2] = 0;
				connec_order_out3[3] = 1;
				connec_order_out3[4] = 4;
				connec_order_out3[5] = 2;
				connec_order_out3[6] = 3;
				connec_order_out3[7] = 5;

				elem->cut = 2;

			}//Edge 2
			else if ((!edge_list[0]) && (edge_list[1]) && (!edge_list[2]) && (edge_list[3]) &&
					(!edge_list[4]) && (!edge_list[5]) && (edge_list[6]) && (edge_list[7]) &&
					(!edge_list[8]) && (!edge_list[9]) && (!edge_list[10]) && (!edge_list[11])&&
					(edge_list_s[5]) && (edge_list_s[7]) ) {
				//printf("Entrou na 2!\n");

				cut_edge[0] = 3;
				cut_edge[1] = 1;
				cut_edge[2] = 7;
				cut_edge[3] = 6;

				cut_edge_s[0] = 5;
				cut_edge_s[1] = 7;

				connec_order_in1[0] = 0;
				connec_order_in1[1] = 1;
				connec_order_in1[2] = 4;
				connec_order_in1[3] = 5;
				connec_order_in1[4] = 2;
				connec_order_in1[5] = 3;
				connec_order_in1[6] = 6;
				connec_order_in1[7] = 7;

				connec_order_out1[0] = 0;
				connec_order_out1[1] = 1;
				connec_order_out1[2] = 4;
				connec_order_out1[3] = 5;
				connec_order_out1[4] = 1;
				connec_order_out1[5] = 0;
				connec_order_out1[6] = 5;
				connec_order_out1[7] = 4;

				connec_order_in2[0] = 0;
				connec_order_in2[1] = 1;
				connec_order_in2[2] = 2;
				connec_order_in2[3] = 3;
				connec_order_in2[4] = 4;
				connec_order_in2[5] = 5;
				connec_order_in2[6] = 6;
				connec_order_in2[7] = 7;

				connec_order_out2[0] = 4;
				connec_order_out2[1] = 5;
				connec_order_out2[2] = 3;
				connec_order_out2[3] = 2;
				connec_order_out2[4] = 4;
				connec_order_out2[5] = 5;
				connec_order_out2[6] = 6;
				connec_order_out2[7] = 7;

				connec_order_in3[0] = 2;
				connec_order_in3[1] = 3;
				connec_order_in3[2] = 0;
				connec_order_in3[3] = 1;
				connec_order_in3[4] = 4;
				connec_order_in3[5] = 5;
				connec_order_in3[6] = 6;
				connec_order_in3[7] = 7;

				connec_order_out3[0] = 2;
				connec_order_out3[1] = 3;
				connec_order_out3[2] = 0;
				connec_order_out3[3] = 1;
				connec_order_out3[4] = 4;
				connec_order_out3[5] = 5;
				connec_order_out3[6] = 3;
				connec_order_out3[7] = 2;

				elem->cut = 2;

			}//Edge 3
			else if ((edge_list[0]) && (!edge_list[1]) && (edge_list[2]) && (!edge_list[3]) &&
					(edge_list[4]) && (!edge_list[5]) && (!edge_list[6]) && (edge_list[7]) &&
					(!edge_list[8]) && (!edge_list[9]) && (!edge_list[10]) && (!edge_list[11])&&
					(edge_list_s[0]) && (edge_list_s[2])) {
				//printf("Entrou na 3!\n");

				cut_edge[0] = 0;
				cut_edge[1] = 2;
				cut_edge[2] = 4;
				cut_edge[3] = 7;

				cut_edge_s[0] = 0;
				cut_edge_s[1] = 2;

				connec_order_in1[0] = 1;
				connec_order_in1[1] = 2;
				connec_order_in1[2] = 5;
				connec_order_in1[3] = 6;
				connec_order_in1[4] = 0;
				connec_order_in1[5] = 3;
				connec_order_in1[6] = 4;
				connec_order_in1[7] = 7;

				connec_order_out1[0] = 1;
				connec_order_out1[1] = 2;
				connec_order_out1[2] = 5;
				connec_order_out1[3] = 6;
				connec_order_out1[4] = 0;
				connec_order_out1[5] = 1;
				connec_order_out1[6] = 4;
				connec_order_out1[7] = 5;

				connec_order_in2[0] = 0;
				connec_order_in2[1] = 1;
				connec_order_in2[2] = 2;
				connec_order_in2[3] = 3;
				connec_order_in2[4] = 4;
				connec_order_in2[5] = 5;
				connec_order_in2[6] = 6;
				connec_order_in2[7] = 7;

				connec_order_out2[0] = 2;
				connec_order_out2[1] = 4;
				connec_order_out2[2] = 5;
				connec_order_out2[3] = 3;
				connec_order_out2[4] = 4;
				connec_order_out2[5] = 5;
				connec_order_out2[6] = 6;
				connec_order_out2[7] = 7;

				connec_order_in3[0] = 0;
				connec_order_in3[1] = 3;
				connec_order_in3[2] = 1;
				connec_order_in3[3] = 2;
				connec_order_in3[4] = 4;
				connec_order_in3[5] = 5;
				connec_order_in3[6] = 6;
				connec_order_in3[7] = 7;

				connec_order_out3[0] = 0;
				connec_order_out3[1] = 3;
				connec_order_out3[2] = 0;
				connec_order_out3[3] = 1;
				connec_order_out3[4] = 2;
				connec_order_out3[5] = 4;
				connec_order_out3[6] = 5;
				connec_order_out3[7] = 3;

				elem->cut = 2;

			}//Edge 4
			else if ((edge_list[0]) && (!edge_list[1]) && (!edge_list[2]) && (edge_list[3]) &&
					(!edge_list[4]) && (!edge_list[5]) && (!edge_list[6]) && (!edge_list[7]) &&
					(edge_list[8]) && (!edge_list[9]) && (!edge_list[10]) && (edge_list[11])&&
					(edge_list_s[8]) && (edge_list_s[10])) {
				//printf("Entrou na 4!\n");

				cut_edge[0] = 0;
				cut_edge[1] = 3;
				cut_edge[2] = 11;
				cut_edge[3] = 8;

				cut_edge_s[0] = 10;
				cut_edge_s[1] = 8;

				connec_order_in1[0] = 1;
				connec_order_in1[1] = 2;
				connec_order_in1[2] = 5;
				connec_order_in1[3] = 6;
				connec_order_in1[4] = 0;
				connec_order_in1[5] = 3;
				connec_order_in1[6] = 4;
				connec_order_in1[7] = 7;

				connec_order_out1[0] = 1;
				connec_order_out1[1] = 2;
				connec_order_out1[2] = 5;
				connec_order_out1[3] = 6;
				connec_order_out1[4] = 0;
				connec_order_out1[5] = 4;
				connec_order_out1[6] = 3;
				connec_order_out1[7] = 5;

				connec_order_in2[0] = 0;
				connec_order_in2[1] = 1;
				connec_order_in2[2] = 4;
				connec_order_in2[3] = 5;
				connec_order_in2[4] = 2;
				connec_order_in2[5] = 3;
				connec_order_in2[6] = 6;
				connec_order_in2[7] = 7;

				connec_order_out2[0] = 1;
				connec_order_out2[1] = 4;
				connec_order_out2[2] = 2;
				connec_order_out2[3] = 5;
				connec_order_out2[4] = 2;
				connec_order_out2[5] = 3;
				connec_order_out2[6] = 6;
				connec_order_out2[7] = 7;

				connec_order_in3[0] = 0;
				connec_order_in3[1] = 4;
				connec_order_in3[2] = 1;
				connec_order_in3[3] = 2;
				connec_order_in3[4] = 3;
				connec_order_in3[5] = 5;
				connec_order_in3[6] = 6;
				connec_order_in3[7] = 7;

				connec_order_out3[0] = 0;
				connec_order_out3[1] = 4;
				connec_order_out3[2] = 0;
				connec_order_out3[3] = 4;
				connec_order_out3[4] = 1;
				connec_order_out3[5] = 3;
				connec_order_out3[6] = 5;
				connec_order_out3[7] = 2;

				elem->cut = 2;

			}//Edge 5
			else if ((edge_list[0]) && (edge_list[1]) && (!edge_list[2]) && (!edge_list[3]) &&
					(!edge_list[4]) && (!edge_list[5]) && (!edge_list[6]) && (!edge_list[7]) &&
					(edge_list[8]) && (edge_list[9]) && (!edge_list[10]) && (!edge_list[11])&&
					(edge_list_s[9]) && (edge_list_s[11])) {
				//printf("Entrou na 5!\n");

				cut_edge[0] = 1;
				cut_edge[1] = 0;
				cut_edge[2] = 8;
				cut_edge[3] = 9;

				cut_edge_s[0] = 11;
				cut_edge_s[1] = 9;

				connec_order_in1[0] = 0;
				connec_order_in1[1] = 3;
				connec_order_in1[2] = 4;
				connec_order_in1[3] = 7;
				connec_order_in1[4] = 1;
				connec_order_in1[5] = 2;
				connec_order_in1[6] = 5;
				connec_order_in1[7] = 6;

				connec_order_out1[0] = 0;
				connec_order_out1[1] = 3;
				connec_order_out1[2] = 4;
				connec_order_out1[3] = 7;
				connec_order_out1[4] = 1;
				connec_order_out1[5] = 4;
				connec_order_out1[6] = 2;
				connec_order_out1[7] = 5;

				connec_order_in2[0] = 0;
				connec_order_in2[1] = 1;
				connec_order_in2[2] = 4;
				connec_order_in2[3] = 5;
				connec_order_in2[4] = 2;
				connec_order_in2[5] = 3;
				connec_order_in2[6] = 6;
				connec_order_in2[7] = 7;

				connec_order_out2[0] = 4;
				connec_order_out2[1] = 0;
				connec_order_out2[2] = 5;
				connec_order_out2[3] = 3;
				connec_order_out2[4] = 2;
				connec_order_out2[5] = 3;
				connec_order_out2[6] = 6;
				connec_order_out2[7] = 7;

				connec_order_in3[0] = 1;
				connec_order_in3[1] = 5;
				connec_order_in3[2] = 0;
				connec_order_in3[3] = 2;
				connec_order_in3[4] = 3;
				connec_order_in3[5] = 4;
				connec_order_in3[6] = 6;
				connec_order_in3[7] = 7;

				connec_order_out3[0] = 1;
				connec_order_out3[1] = 5;
				connec_order_out3[2] = 1;
				connec_order_out3[3] = 0;
				connec_order_out3[4] = 4;
				connec_order_out3[5] = 2;
				connec_order_out3[6] = 3;
				connec_order_out3[7] = 5;

				elem->cut = 2;

			}//Edge 6
			else if ((!edge_list[0]) && (edge_list[1]) && (edge_list[2]) && (!edge_list[3]) &&
					(!edge_list[4]) && (!edge_list[5]) && (!edge_list[6]) && (!edge_list[7]) &&
					(!edge_list[8]) && (edge_list[9]) && (edge_list[10]) && (!edge_list[11])&&
					(edge_list_s[8]) && (edge_list_s[10])) {
				//printf("Entrou na 6!\n");

				cut_edge[0] = 2;
				cut_edge[1] = 10;
				cut_edge[2] = 1;
				cut_edge[3] = 9;

				cut_edge_s[0] = 10;
				cut_edge_s[1] = 8;

				connec_order_in1[0] = 0;
				connec_order_in1[1] = 1;
				connec_order_in1[2] = 4;
				connec_order_in1[3] = 5;
				connec_order_in1[4] = 2;
				connec_order_in1[5] = 3;
				connec_order_in1[6] = 6;
				connec_order_in1[7] = 7;

				connec_order_out1[0] = 0;
				connec_order_out1[1] = 1;
				connec_order_out1[2] = 4;
				connec_order_out1[3] = 5;
				connec_order_out1[4] = 2;
				connec_order_out1[5] = 4;
				connec_order_out1[6] = 3;
				connec_order_out1[7] = 5;

				connec_order_in2[0] = 1;
				connec_order_in2[1] = 2;
				connec_order_in2[2] = 5;
				connec_order_in2[3] = 6;
				connec_order_in2[4] = 0;
				connec_order_in2[5] = 3;
				connec_order_in2[6] = 4;
				connec_order_in2[7] = 7;

				connec_order_out2[0] = 4;
				connec_order_out2[1] = 0;
				connec_order_out2[2] = 5;
				connec_order_out2[3] = 1;
				connec_order_out2[4] = 0;
				connec_order_out2[5] = 3;
				connec_order_out2[6] = 4;
				connec_order_out2[7] = 7;

				connec_order_in3[0] = 2;
				connec_order_in3[1] = 6;
				connec_order_in3[2] = 0;
				connec_order_in3[3] = 1;
				connec_order_in3[4] = 3;
				connec_order_in3[5] = 4;
				connec_order_in3[6] = 5;
				connec_order_in3[7] = 7;

				connec_order_out3[0] = 2;
				connec_order_out3[1] = 6;
				connec_order_out3[2] = 4;
				connec_order_out3[3] = 2;
				connec_order_out3[4] = 0;
				connec_order_out3[5] = 5;
				connec_order_out3[6] = 3;
				connec_order_out3[7] = 1;

				elem->cut = 2;

			}//Edge 7
			else if ((!edge_list[0]) && (!edge_list[1]) && (edge_list[2]) && (edge_list[3]) &&
					(!edge_list[4]) && (!edge_list[5]) && (!edge_list[6]) && (!edge_list[7]) &&
					(!edge_list[8]) && (!edge_list[9]) && (edge_list[10]) && (edge_list[11])&&
					(edge_list_s[9]) && (edge_list_s[11])) {
				//printf("Entrou na 7!\n");

				cut_edge[0] = 3;
				cut_edge[1] = 11;
				cut_edge[2] = 2;
				cut_edge[3] = 10;

				cut_edge_s[0] = 11;
				cut_edge_s[1] = 9;

				connec_order_in1[0] = 0;
				connec_order_in1[1] = 1;
				connec_order_in1[2] = 4;
				connec_order_in1[3] = 5;
				connec_order_in1[4] = 2;
				connec_order_in1[5] = 3;
				connec_order_in1[6] = 6;
				connec_order_in1[7] = 7;

				connec_order_out1[0] = 0;
				connec_order_out1[1] = 1;
				connec_order_out1[2] = 4;
				connec_order_out1[3] = 5;
				connec_order_out1[4] = 4;
				connec_order_out1[5] = 0;
				connec_order_out1[6] = 5;
				connec_order_out1[7] = 1;

				connec_order_in2[0] = 0;
				connec_order_in2[1] = 3;
				connec_order_in2[2] = 4;
				connec_order_in2[3] = 7;
				connec_order_in2[4] = 1;
				connec_order_in2[5] = 2;
				connec_order_in2[6] = 5;
				connec_order_in2[7] = 6;

				connec_order_out2[0] = 4;
				connec_order_out2[1] = 2;
				connec_order_out2[2] = 5;
				connec_order_out2[3] = 3;
				connec_order_out2[4] = 1;
				connec_order_out2[5] = 2;
				connec_order_out2[6] = 5;
				connec_order_out2[7] = 6;

				connec_order_in3[0] = 3;
				connec_order_in3[1] = 7;
				connec_order_in3[2] = 0;
				connec_order_in3[3] = 1;
				connec_order_in3[4] = 2;
				connec_order_in3[5] = 4;
				connec_order_in3[6] = 5;
				connec_order_in3[7] = 6;

				connec_order_out3[0] = 3;
				connec_order_out3[1] = 7;
				connec_order_out3[2] = 0;
				connec_order_out3[3] = 4;
				connec_order_out3[4] = 2;
				connec_order_out3[5] = 1;
				connec_order_out3[6] = 5;
				connec_order_out3[7] = 3;

				elem->cut = 2;

			}//Edge 8
			else if ((!edge_list[0]) && (!edge_list[1]) && (!edge_list[2]) && (!edge_list[3]) &&
					(edge_list[4]) && (edge_list[5]) && (!edge_list[6]) && (!edge_list[7]) &&
					(!edge_list[8]) && (edge_list[9]) && (!edge_list[10]) && (edge_list[11])&&
					(edge_list_s[5]) && (edge_list_s[7])) {
				//printf("Entrou na 8!\n");

				cut_edge[0] = 11;
				cut_edge[1] = 9;
				cut_edge[2] = 4;
				cut_edge[3] = 5;

				cut_edge_s[0] = 5;
				cut_edge_s[1] = 7;

				connec_order_in1[0] = 2;
				connec_order_in1[1] = 3;
				connec_order_in1[2] = 6;
				connec_order_in1[3] = 7;
				connec_order_in1[4] = 0;
				connec_order_in1[5] = 1;
				connec_order_in1[6] = 4;
				connec_order_in1[7] = 5;

				connec_order_out1[0] = 2;
				connec_order_out1[1] = 3;
				connec_order_out1[2] = 6;
				connec_order_out1[3] = 7;
				connec_order_out1[4] = 4;
				connec_order_out1[5] = 5;
				connec_order_out1[6] = 0;
				connec_order_out1[7] = 1;

				connec_order_in2[0] = 4;
				connec_order_in2[1] = 5;
				connec_order_in2[2] = 6;
				connec_order_in2[3] = 7;
				connec_order_in2[4] = 0;
				connec_order_in2[5] = 1;
				connec_order_in2[6] = 2;
				connec_order_in2[7] = 3;

				connec_order_out2[0] = 2;
				connec_order_out2[1] = 3;
				connec_order_out2[2] = 5;
				connec_order_out2[3] = 4;
				connec_order_out2[4] = 0;
				connec_order_out2[5] = 1;
				connec_order_out2[6] = 2;
				connec_order_out2[7] = 3;

				connec_order_in3[0] = 4;
				connec_order_in3[1] = 5;
				connec_order_in3[2] = 0;
				connec_order_in3[3] = 1;
				connec_order_in3[4] = 2;
				connec_order_in3[5] = 3;
				connec_order_in3[6] = 6;
				connec_order_in3[7] = 7;

				connec_order_out3[0] = 4;
				connec_order_out3[1] = 5;
				connec_order_out3[2] = 2;
				connec_order_out3[3] = 3;
				connec_order_out3[4] = 5;
				connec_order_out3[5] = 4;
				connec_order_out3[6] = 1;
				connec_order_out3[7] = 0;

				elem->cut = 2;

			}//Edge 9
			else if ((!edge_list[0]) && (!edge_list[1]) && (!edge_list[2]) && (!edge_list[3]) &&
					(!edge_list[4]) && (edge_list[5]) && (edge_list[6]) && (!edge_list[7]) &&
					(edge_list[8]) && (!edge_list[9]) && (edge_list[10]) && (!edge_list[11])&&
					(edge_list_s[0]) && (edge_list_s[2])) {
				//printf("Entrou na 9!\n");

				cut_edge[0] = 8;
				cut_edge[1] = 10;
				cut_edge[2] = 5;
				cut_edge[3] = 6;

				cut_edge_s[0] = 0;
				cut_edge_s[1] = 2;

				connec_order_in1[0] = 0;
				connec_order_in1[1] = 3;
				connec_order_in1[2] = 4;
				connec_order_in1[3] = 7;
				connec_order_in1[4] = 1;
				connec_order_in1[5] = 2;
				connec_order_in1[6] = 5;
				connec_order_in1[7] = 6;

				connec_order_out1[0] = 0;
				connec_order_out1[1] = 3;
				connec_order_out1[2] = 4;
				connec_order_out1[3] = 7;
				connec_order_out1[4] = 4;
				connec_order_out1[5] = 5;
				connec_order_out1[6] = 0;
				connec_order_out1[7] = 1;

				connec_order_in2[0] = 4;
				connec_order_in2[1] = 5;
				connec_order_in2[2] = 6;
				connec_order_in2[3] = 7;
				connec_order_in2[4] = 0;
				connec_order_in2[5] = 1;
				connec_order_in2[6] = 2;
				connec_order_in2[7] = 3;

				connec_order_out2[0] = 4;
				connec_order_out2[1] = 2;
				connec_order_out2[2] = 3;
				connec_order_out2[3] = 5;
				connec_order_out2[4] = 0;
				connec_order_out2[5] = 1;
				connec_order_out2[6] = 2;
				connec_order_out2[7] = 3;

				connec_order_in3[0] = 5;
				connec_order_in3[1] = 6;
				connec_order_in3[2] = 0;
				connec_order_in3[3] = 1;
				connec_order_in3[4] = 2;
				connec_order_in3[5] = 3;
				connec_order_in3[6] = 4;
				connec_order_in3[7] = 7;

				connec_order_out3[0] = 5;
				connec_order_out3[1] = 6;
				connec_order_out3[2] = 4;
				connec_order_out3[3] = 2;
				connec_order_out3[4] = 3;
				connec_order_out3[5] = 5;
				connec_order_out3[6] = 0;
				connec_order_out3[7] = 1;

				elem->cut = 2;

			}//Edge 10
			else if ((!edge_list[0]) && (!edge_list[1]) && (!edge_list[2]) && (!edge_list[3]) &&
					(!edge_list[4]) && (!edge_list[5]) && (edge_list[6]) && (edge_list[7]) &&
					(!edge_list[8]) && (edge_list[9]) && (!edge_list[10]) && (edge_list[11])&&
					(edge_list_s[4]) && (edge_list_s[6])) {
				//printf("Entrou na 10!\n");

				cut_edge[0] = 9;
				cut_edge[1] = 11;
				cut_edge[2] = 6;
				cut_edge[3] = 7;

				cut_edge_s[0] = 6;
				cut_edge_s[1] = 4;

				connec_order_in1[0] = 0;
				connec_order_in1[1] = 1;
				connec_order_in1[2] = 4;
				connec_order_in1[3] = 5;
				connec_order_in1[4] = 2;
				connec_order_in1[5] = 3;
				connec_order_in1[6] = 6;
				connec_order_in1[7] = 7;

				connec_order_out1[0] = 0;
				connec_order_out1[1] = 1;
				connec_order_out1[2] = 4;
				connec_order_out1[3] = 5;
				connec_order_out1[4] = 4;
				connec_order_out1[5] = 5;
				connec_order_out1[6] = 0;
				connec_order_out1[7] = 1;

				connec_order_in2[0] = 4;
				connec_order_in2[1] = 5;
				connec_order_in2[2] = 6;
				connec_order_in2[3] = 7;
				connec_order_in2[4] = 0;
				connec_order_in2[5] = 1;
				connec_order_in2[6] = 2;
				connec_order_in2[7] = 3;

				connec_order_out2[0] = 5;
				connec_order_out2[1] = 4;
				connec_order_out2[2] = 2;
				connec_order_out2[3] = 3;
				connec_order_out2[4] = 0;
				connec_order_out2[5] = 1;
				connec_order_out2[6] = 2;
				connec_order_out2[7] = 3;

				connec_order_in3[0] = 6;
				connec_order_in3[1] = 7;
				connec_order_in3[2] = 0;
				connec_order_in3[3] = 1;
				connec_order_in3[4] = 2;
				connec_order_in3[5] = 3;
				connec_order_in3[6] = 4;
				connec_order_in3[7] = 5;

				connec_order_out3[0] = 6;
				connec_order_out3[1] = 7;
				connec_order_out3[2] = 5;
				connec_order_out3[3] = 4;
				connec_order_out3[4] = 2;
				connec_order_out3[5] = 3;
				connec_order_out3[6] = 1;
				connec_order_out3[7] = 0;

				elem->cut = 2;

			}//Edge 11
			else if ((!edge_list[0]) && (!edge_list[1]) && (!edge_list[2]) && (!edge_list[3]) &&
					(edge_list[4]) && (!edge_list[5]) && (!edge_list[6]) && (edge_list[7]) &&
					(edge_list[8]) && (!edge_list[9]) && (edge_list[10]) && (!edge_list[11])&&
					(edge_list_s[1]) && (edge_list_s[3])) {
				//printf("Entrou na 11!\n");

				cut_edge[0] = 8;
				cut_edge[1] = 10;
				cut_edge[2] = 4;
				cut_edge[3] = 7;

				cut_edge_s[0] = 1;
				cut_edge_s[1] = 3;

				connec_order_in1[0] = 1;
				connec_order_in1[1] = 2;
				connec_order_in1[2] = 5;
				connec_order_in1[3] = 6;
				connec_order_in1[4] = 0;
				connec_order_in1[5] = 3;
				connec_order_in1[6] = 4;
				connec_order_in1[7] = 7;

				connec_order_out1[0] = 1;
				connec_order_out1[1] = 2;
				connec_order_out1[2] = 5;
				connec_order_out1[3] = 6;
				connec_order_out1[4] = 4;
				connec_order_out1[5] = 5;
				connec_order_out1[6] = 0;
				connec_order_out1[7] = 1;

				connec_order_in2[0] = 4;
				connec_order_in2[1] = 5;
				connec_order_in2[2] = 6;
				connec_order_in2[3] = 7;
				connec_order_in2[4] = 0;
				connec_order_in2[5] = 1;
				connec_order_in2[6] = 2;
				connec_order_in2[7] = 3;

				connec_order_out2[0] = 2;
				connec_order_out2[1] = 4;
				connec_order_out2[2] = 5;
				connec_order_out2[3] = 3;
				connec_order_out2[4] = 0;
				connec_order_out2[5] = 1;
				connec_order_out2[6] = 2;
				connec_order_out2[7] = 3;

				connec_order_in3[0] = 4;
				connec_order_in3[1] = 7;
				connec_order_in3[2] = 0;
				connec_order_in3[3] = 1;
				connec_order_in3[4] = 2;
				connec_order_in3[5] = 3;
				connec_order_in3[6] = 5;
				connec_order_in3[7] = 6;

				connec_order_out3[0] = 4;
				connec_order_out3[1] = 7;
				connec_order_out3[2] = 2;
				connec_order_out3[3] = 4;
				connec_order_out3[4] = 5;
				connec_order_out3[5] = 3;
				connec_order_out3[6] = 0;
				connec_order_out3[7] = 1;

				elem->cut = 2;

			} else {
				elem->cut = -1;
				printf("Warning: Bad Element: %d, case 2\n",elements_ids[iel]);
				printf("Edges:  Elem. %d: %d %d %d %d %d %d %d %d %d %d %d %d %d\n", elements_ids[iel], edge_list[0],
						edge_list[1], edge_list[2], edge_list[3],
						edge_list[4], edge_list[5], edge_list[6],
						edge_list[7], edge_list[8], edge_list[9],
						edge_list[10], edge_list[11], ed_cont);
				continue;
			}

			GtsPoint *p0 = point[cut_edge[0]];
			GtsPoint *p1 = point[cut_edge[1]];
			GtsPoint *p2 = point[cut_edge[2]];
			GtsPoint *p3 = point[cut_edge[3]];

			GtsPoint *p4 = point_s[cut_edge_s[0]];
			GtsPoint *p5 = point_s[cut_edge_s[1]];

			//g_assert(p0 != NULL);
			//g_assert(p1 != NULL);
			//g_assert(p2 != NULL);
			//g_assert(p3 != NULL);
			//g_assert(p4 != NULL);
			//g_assert(p5 != NULL);

			if(p0==NULL) exit(1);
			if(p1==NULL) exit(2);
			if(p2==NULL) exit(3);
			if(p3==NULL) exit(4);
			if(p4==NULL) exit(5);
			if(p5==NULL) exit(6);

			conn_t2[0] = AddPointOnEdge(Edge2GNode[cut_edge[0]], hash_nodes, mesh->local_n_nodes, p0, coords);
			conn_t2[1] = AddPointOnEdge(Edge2GNode[cut_edge[1]], hash_nodes, mesh->local_n_nodes, p1, coords);
			conn_t2[2] = AddPointOnEdge(Edge2GNode[cut_edge[2]], hash_nodes, mesh->local_n_nodes, p2, coords);
			conn_t2[3] = AddPointOnEdge(Edge2GNode[cut_edge[3]], hash_nodes, mesh->local_n_nodes, p3, coords);

			// add 2 extra points in the surface
			conn_t2[4] = AddPointOnEdge(Edge2GNode_s[cut_edge_s[0]], hash_nodes, mesh->local_n_nodes, p4, coords);
			conn_t2[5] = AddPointOnEdge(Edge2GNode_s[cut_edge_s[1]], hash_nodes, mesh->local_n_nodes, p5, coords);

			octant_t *elem1 = (octant_t*) sc_array_index(&mesh->elements, elements_ids[iel]);

			elem1->nodes[connec_order_in1[0]].id = original_conn[connec_order_out1[0]];
			elem1->nodes[connec_order_in1[1]].id = original_conn[connec_order_out1[1]];
			elem1->nodes[connec_order_in1[2]].id = original_conn[connec_order_out1[2]];
			elem1->nodes[connec_order_in1[3]].id = original_conn[connec_order_out1[3]];

			elem1->nodes[connec_order_in1[4]].id = conn_t2[connec_order_out1[4]];
			elem1->nodes[connec_order_in1[5]].id = conn_t2[connec_order_out1[5]];
			elem1->nodes[connec_order_in1[6]].id = conn_t2[connec_order_out1[6]];
			elem1->nodes[connec_order_in1[7]].id = conn_t2[connec_order_out1[7]];

			octant_t* elem2 = (octant_t*) sc_array_push(&mesh->elements);

			elem2->nodes[connec_order_in2[0]].id = conn_t2[connec_order_out2[0]];
			elem2->nodes[connec_order_in2[1]].id = conn_t2[connec_order_out2[1]];
			elem2->nodes[connec_order_in2[2]].id = conn_t2[connec_order_out2[2]];
			elem2->nodes[connec_order_in2[3]].id = conn_t2[connec_order_out2[3]];

			elem2->nodes[connec_order_in2[4]].id = original_conn[connec_order_out2[4]];
			elem2->nodes[connec_order_in2[5]].id = original_conn[connec_order_out2[5]];
			elem2->nodes[connec_order_in2[6]].id = original_conn[connec_order_out2[6]];
			elem2->nodes[connec_order_in2[7]].id = original_conn[connec_order_out2[7]];

			elem2->cut = elem->cut;
			elem2->level = elem->level;

			octant_t* elem3 = (octant_t*) sc_array_push(&mesh->elements);

			elem3->nodes[connec_order_in3[0]].id = original_conn[connec_order_out3[0]];
			elem3->nodes[connec_order_in3[1]].id = original_conn[connec_order_out3[1]];

			elem3->nodes[connec_order_in3[2]].id = conn_t2[connec_order_out3[2]];
			elem3->nodes[connec_order_in3[3]].id = conn_t2[connec_order_out3[3]];
			elem3->nodes[connec_order_in3[4]].id = conn_t2[connec_order_out3[4]];
			elem3->nodes[connec_order_in3[5]].id = conn_t2[connec_order_out3[5]];
			elem3->nodes[connec_order_in3[6]].id = conn_t2[connec_order_out3[6]];
			elem3->nodes[connec_order_in3[7]].id = conn_t2[connec_order_out3[7]];

			elem3->cut = elem->cut;
			elem3->level = elem->level;

		}else if(elem->cut==3){

			// Check template 3.
			int conn_t2[7];
			int cut_edge[3];
			int cut_edge_s[3];
			int cut_edge_v[1];
			int connec_order_in1[8];
			int connec_order_out1[8];
			int connec_order_in2[8];
			int connec_order_out2[8];
			int connec_order_in3[8];
			int connec_order_out3[8];
			int connec_order_in4[8];
			int connec_order_out4[8];


			//Corner 0
			if ((edge_list[0]) && (!edge_list[1]) && (!edge_list[2]) && (edge_list[3]) &&
					(edge_list[4]) && (!edge_list[5]) && (!edge_list[6]) && (!edge_list[7]) &&
					(!edge_list[8]) && (!edge_list[9]) && (!edge_list[10]) && (!edge_list[11]) &&
					(edge_list_s[4]) && (edge_list_s[0]) && (edge_list_s[10]) &&
					(edge_list_v[0]) ) {
				//printf("Entrou no 0!\n");

				cut_edge[0] = 4;
				cut_edge[1] = 0;
				cut_edge[2] = 3;

				cut_edge_s[0] = 0;
				cut_edge_s[1] = 4;
				cut_edge_s[2] = 10;

				cut_edge_v[0] = 0;

				connec_order_in1[0] = 4;
				connec_order_in1[1] = 5;
				connec_order_in1[2] = 6;
				connec_order_in1[3] = 7;
				connec_order_in1[4] = 0;
				connec_order_in1[5] = 1;
				connec_order_in1[6] = 2;
				connec_order_in1[7] = 3;

				connec_order_out1[0] = 4;
				connec_order_out1[1] = 5;
				connec_order_out1[2] = 6;
				connec_order_out1[3] = 7;
				connec_order_out1[4] = 0;
				connec_order_out1[5] = 3;
				connec_order_out1[6] = 6;
				connec_order_out1[7] = 4;

				connec_order_in2[0] = 0;
				connec_order_in2[1] = 3;
				connec_order_in2[2] = 4;
				connec_order_in2[3] = 7;
				connec_order_in2[4] = 1;
				connec_order_in2[5] = 2;
				connec_order_in2[6] = 5;
				connec_order_in2[7] = 6;

				connec_order_out2[0] = 1;
				connec_order_out2[1] = 5;
				connec_order_out2[2] = 3;
				connec_order_out2[3] = 6;
				connec_order_out2[4] = 1;
				connec_order_out2[5] = 2;
				connec_order_out2[6] = 5;
				connec_order_out2[7] = 6;

				connec_order_in3[0] = 0;
				connec_order_in3[1] = 1;
				connec_order_in3[2] = 4;
				connec_order_in3[3] = 5;
				connec_order_in3[4] = 2;
				connec_order_in3[5] = 3;
				connec_order_in3[6] = 6;
				connec_order_in3[7] = 7;

				connec_order_out3[0] = 2;
				connec_order_out3[1] = 5;
				connec_order_out3[2] = 4;
				connec_order_out3[3] = 6;
				connec_order_out3[4] = 2;
				connec_order_out3[5] = 3;
				connec_order_out3[6] = 6;
				connec_order_out3[7] = 7;

				connec_order_in4[0] = 0;
				connec_order_in4[1] = 1;
				connec_order_in4[2] = 2;
				connec_order_in4[3] = 3;
				connec_order_in4[4] = 4;
				connec_order_in4[5] = 5;
				connec_order_in4[6] = 6;
				connec_order_in4[7] = 7;

				connec_order_out4[0] = 0;
				connec_order_out4[1] = 1;
				connec_order_out4[2] = 5;
				connec_order_out4[3] = 2;
				connec_order_out4[4] = 0;
				connec_order_out4[5] = 3;
				connec_order_out4[6] = 6;
				connec_order_out4[7] = 4;

				elem->cut = 3;

			}

			//Corner 1
			else if ((edge_list[0]) && (edge_list[1]) && (!edge_list[2]) && (!edge_list[3]) &&
					(!edge_list[4]) && (edge_list[5]) && (!edge_list[6]) && (!edge_list[7]) &&
					(!edge_list[8]) && (!edge_list[9]) && (!edge_list[10]) && (!edge_list[11])&&
					(edge_list_s[1]) && (edge_list_s[6]) && (edge_list_s[11])  &&
					(edge_list_v[1]) ) {
				//printf("Entrou no 1!\n");

				cut_edge[0] = 5;
				cut_edge[1] = 1;
				cut_edge[2] = 0;

				cut_edge_s[0] = 6;
				cut_edge_s[1] = 1;
				cut_edge_s[2] = 11;

				cut_edge_v[0] = 1;

				connec_order_in1[0] = 4;
				connec_order_in1[1] = 5;
				connec_order_in1[2] = 6;
				connec_order_in1[3] = 7;
				connec_order_in1[4] = 0;
				connec_order_in1[5] = 1;
				connec_order_in1[6] = 2;
				connec_order_in1[7] = 3;

				connec_order_out1[0] = 4;
				connec_order_out1[1] = 5;
				connec_order_out1[2] = 6;
				connec_order_out1[3] = 7;
				connec_order_out1[4] = 4;
				connec_order_out1[5] = 0;
				connec_order_out1[6] = 3;
				connec_order_out1[7] = 6;

				connec_order_in2[0] = 0;
				connec_order_in2[1] = 1;
				connec_order_in2[2] = 4;
				connec_order_in2[3] = 5;
				connec_order_in2[4] = 2;
				connec_order_in2[5] = 3;
				connec_order_in2[6] = 6;
				connec_order_in2[7] = 7;

				connec_order_out2[0] = 5;
				connec_order_out2[1] = 1;
				connec_order_out2[2] = 6;
				connec_order_out2[3] = 3;
				connec_order_out2[4] = 2;
				connec_order_out2[5] = 3;
				connec_order_out2[6] = 6;
				connec_order_out2[7] = 7;

				connec_order_in3[0] = 1;
				connec_order_in3[1] = 2;
				connec_order_in3[2] = 5;
				connec_order_in3[3] = 6;
				connec_order_in3[4] = 0;
				connec_order_in3[5] = 3;
				connec_order_in3[6] = 4;
				connec_order_in3[7] = 7;

				connec_order_out3[0] = 2;
				connec_order_out3[1] = 5;
				connec_order_out3[2] = 4;
				connec_order_out3[3] = 6;
				connec_order_out3[4] = 0;
				connec_order_out3[5] = 3;
				connec_order_out3[6] = 4;
				connec_order_out3[7] = 7;

				connec_order_in4[0] = 1;
				connec_order_in4[1] = 0;
				connec_order_in4[2] = 2;
				connec_order_in4[3] = 3;
				connec_order_in4[4] = 4;
				connec_order_in4[5] = 5;
				connec_order_in4[6] = 6;
				connec_order_in4[7] = 7;

				connec_order_out4[0] = 1;
				connec_order_out4[1] = 2;
				connec_order_out4[2] = 1;
				connec_order_out4[3] = 5;
				connec_order_out4[4] = 4;
				connec_order_out4[5] = 0;
				connec_order_out4[6] = 3;
				connec_order_out4[7] = 6;

				elem->cut = 3;

			}

			//Corner 2
			else if ((!edge_list[0]) && (edge_list[1]) && (edge_list[2]) && (!edge_list[3]) &&
					(!edge_list[4]) && (!edge_list[5]) && (edge_list[6]) && (!edge_list[7]) &&
					(!edge_list[8]) && (!edge_list[9]) && (!edge_list[10]) && (!edge_list[11]) &&
					(edge_list_s[3]) && (edge_list_s[7]) && (edge_list_s[10])  &&
					(edge_list_v[2]) ) {
				//printf("Entrou no 2!\n");

				cut_edge[0] = 6;
				cut_edge[1] = 2;
				cut_edge[2] = 1;

				cut_edge_s[0] = 3;
				cut_edge_s[1] = 7;
				cut_edge_s[2] = 10;

				cut_edge_v[0] = 2;

				connec_order_in1[0] = 4;
				connec_order_in1[1] = 5;
				connec_order_in1[2] = 6;
				connec_order_in1[3] = 7;
				connec_order_in1[4] = 1;
				connec_order_in1[5] = 2;
				connec_order_in1[6] = 3;
				connec_order_in1[7] = 4;

				connec_order_out1[0] = 4;
				connec_order_out1[1] = 5;
				connec_order_out1[2] = 6;
				connec_order_out1[3] = 7;
				connec_order_out1[4] = 4;
				connec_order_out1[5] = 0;
				connec_order_out1[6] = 3;
				connec_order_out1[7] = 6;

				connec_order_in2[0] = 1;
				connec_order_in2[1] = 2;
				connec_order_in2[2] = 5;
				connec_order_in2[3] = 6;
				connec_order_in2[4] = 0;
				connec_order_in2[5] = 3;
				connec_order_in2[6] = 4;
				connec_order_in2[7] = 7;

				connec_order_out2[0] = 5;
				connec_order_out2[1] = 1;
				connec_order_out2[2] = 6;
				connec_order_out2[3] = 3;
				connec_order_out2[4] = 0;
				connec_order_out2[5] = 3;
				connec_order_out2[6] = 4;
				connec_order_out2[7] = 7;

				connec_order_in3[0] = 2;
				connec_order_in3[1] = 3;
				connec_order_in3[2] = 6;
				connec_order_in3[3] = 7;
				connec_order_in3[4] = 0;
				connec_order_in3[5] = 1;
				connec_order_in3[6] = 4;
				connec_order_in3[7] = 5;

				connec_order_out3[0] = 2;
				connec_order_out3[1] = 5;
				connec_order_out3[2] = 4;
				connec_order_out3[3] = 6;
				connec_order_out3[4] = 0;
				connec_order_out3[5] = 1;
				connec_order_out3[6] = 4;
				connec_order_out3[7] = 5;

				connec_order_in4[0] = 2;
				connec_order_in4[1] = 0;
				connec_order_in4[2] = 1;
				connec_order_in4[3] = 3;
				connec_order_in4[4] = 4;
				connec_order_in4[5] = 5;
				connec_order_in4[6] = 6;
				connec_order_in4[7] = 7;

				connec_order_out4[0] = 2;
				connec_order_out4[1] = 5;
				connec_order_out4[2] = 2;
				connec_order_out4[3] = 1;
				connec_order_out4[4] = 6;
				connec_order_out4[5] = 4;
				connec_order_out4[6] = 0;
				connec_order_out4[7] = 3;

				elem->cut = 3;

			}

			//Corner 3
			else if ((!edge_list[0]) && (!edge_list[1]) && (edge_list[2]) && (edge_list[3]) &&
					(!edge_list[4]) && (!edge_list[5]) && (!edge_list[6]) && (edge_list[7]) &&
					(!edge_list[8]) && (!edge_list[9]) && (!edge_list[10]) && (!edge_list[11]) &&
					(edge_list_s[5]) && (edge_list_s[2]) && (edge_list_s[11])  &&
					(edge_list_v[3]) ) {
				//printf("Entrou no 3!\n");

				cut_edge[0] = 7;
				cut_edge[1] = 3;
				cut_edge[2] = 2;

				cut_edge_s[0] = 5;
				cut_edge_s[1] = 2;
				cut_edge_s[2] = 11;

				cut_edge_v[0] = 3;

				connec_order_in1[0] = 4;
				connec_order_in1[1] = 5;
				connec_order_in1[2] = 6;
				connec_order_in1[3] = 7;
				connec_order_in1[4] = 0;
				connec_order_in1[5] = 1;
				connec_order_in1[6] = 2;
				connec_order_in1[7] = 3;

				connec_order_out1[0] = 4;
				connec_order_out1[1] = 5;
				connec_order_out1[2] = 6;
				connec_order_out1[3] = 7;
				connec_order_out1[4] = 3;
				connec_order_out1[5] = 6;
				connec_order_out1[6] = 4;
				connec_order_out1[7] = 0;

				connec_order_in2[0] = 2;
				connec_order_in2[1] = 3;
				connec_order_in2[2] = 6;
				connec_order_in2[3] = 7;
				connec_order_in2[4] = 0;
				connec_order_in2[5] = 1;
				connec_order_in2[6] = 4;
				connec_order_in2[7] = 5;

				connec_order_out2[0] = 5;
				connec_order_out2[1] = 1;
				connec_order_out2[2] = 6;
				connec_order_out2[3] = 3;
				connec_order_out2[4] = 0;
				connec_order_out2[5] = 1;
				connec_order_out2[6] = 4;
				connec_order_out2[7] = 5;

				connec_order_in3[0] = 0;
				connec_order_in3[1] = 3;
				connec_order_in3[2] = 4;
				connec_order_in3[3] = 7;
				connec_order_in3[4] = 1;
				connec_order_in3[5] = 2;
				connec_order_in3[6] = 5;
				connec_order_in3[7] = 6;

				connec_order_out3[0] = 5;
				connec_order_out3[1] = 2;
				connec_order_out3[2] = 6;
				connec_order_out3[3] = 4;
				connec_order_out3[4] = 1;
				connec_order_out3[5] = 2;
				connec_order_out3[6] = 5;
				connec_order_out3[7] = 6;

				connec_order_in4[0] = 3;
				connec_order_in4[1] = 0;
				connec_order_in4[2] = 1;
				connec_order_in4[3] = 2;
				connec_order_in4[4] = 4;
				connec_order_in4[5] = 5;
				connec_order_in4[6] = 6;
				connec_order_in4[7] = 7;

				connec_order_out4[0] = 3;
				connec_order_out4[1] = 1;
				connec_order_out4[2] = 5;
				connec_order_out4[3] = 2;
				connec_order_out4[4] = 3;
				connec_order_out4[5] = 6;
				connec_order_out4[6] = 4;
				connec_order_out4[7] = 0;

				elem->cut = 3;

			}


			//Corner 4
			else if ((!edge_list[0]) && (!edge_list[1]) && (!edge_list[2]) && (!edge_list[3]) &&
					(edge_list[4]) && (!edge_list[5]) && (!edge_list[6]) && (!edge_list[7]) &&
					(edge_list[8]) && (!edge_list[9]) && (!edge_list[10]) && (edge_list[11]) &&
					(edge_list_s[1]) && (edge_list_s[5]) && (edge_list_s[8])  &&
					(edge_list_v[2]) ) {
				//printf("Entrou no 4!\n");

				cut_edge[0] = 4;
				cut_edge[1] = 8;
				cut_edge[2] = 11;

				cut_edge_s[0] = 1;
				cut_edge_s[1] = 5;
				cut_edge_s[2] = 8;

				cut_edge_v[0] = 2;

				connec_order_in1[0] = 0;
				connec_order_in1[1] = 1;
				connec_order_in1[2] = 2;
				connec_order_in1[3] = 3;
				connec_order_in1[4] = 4;
				connec_order_in1[5] = 5;
				connec_order_in1[6] = 6;
				connec_order_in1[7] = 7;

				connec_order_out1[0] = 0;
				connec_order_out1[1] = 1;
				connec_order_out1[2] = 2;
				connec_order_out1[3] = 3;
				connec_order_out1[4] = 0;
				connec_order_out1[5] = 3;
				connec_order_out1[6] = 6;
				connec_order_out1[7] = 4;

				connec_order_in2[0] = 0;
				connec_order_in2[1] = 3;
				connec_order_in2[2] = 4;
				connec_order_in2[3] = 7;
				connec_order_in2[4] = 1;
				connec_order_in2[5] = 2;
				connec_order_in2[6] = 5;
				connec_order_in2[7] = 6;

				connec_order_out2[0] = 3;
				connec_order_out2[1] = 6;
				connec_order_out2[2] = 1;
				connec_order_out2[3] = 5;
				connec_order_out2[4] = 1;
				connec_order_out2[5] = 2;
				connec_order_out2[6] = 5;
				connec_order_out2[7] = 6;

				connec_order_in3[0] = 0;
				connec_order_in3[1] = 1;
				connec_order_in3[2] = 4;
				connec_order_in3[3] = 5;
				connec_order_in3[4] = 2;
				connec_order_in3[5] = 3;
				connec_order_in3[6] = 6;
				connec_order_in3[7] = 7;

				connec_order_out3[0] = 4;
				connec_order_out3[1] = 6;
				connec_order_out3[2] = 2;
				connec_order_out3[3] = 5;
				connec_order_out3[4] = 2;
				connec_order_out3[5] = 3;
				connec_order_out3[6] = 6;
				connec_order_out3[7] = 7;

				connec_order_in4[0] = 4;
				connec_order_in4[1] = 0;
				connec_order_in4[2] = 1;
				connec_order_in4[3] = 2;
				connec_order_in4[4] = 3;
				connec_order_in4[5] = 5;
				connec_order_in4[6] = 6;
				connec_order_in4[7] = 7;

				connec_order_out4[0] = 4;
				connec_order_out4[1] = 0;
				connec_order_out4[2] = 3;
				connec_order_out4[3] = 6;
				connec_order_out4[4] = 4;
				connec_order_out4[5] = 1;
				connec_order_out4[6] = 5;
				connec_order_out4[7] = 2;

				elem->cut = 3;

			}

			//Corner 5
			else if ((!edge_list[0]) && (!edge_list[1]) && (!edge_list[2]) && (!edge_list[3]) &&
					(!edge_list[4]) && (edge_list[5]) && (!edge_list[6]) && (!edge_list[7]) &&
					(edge_list[8]) && (edge_list[9]) && (!edge_list[10]) && (!edge_list[11]) &&
					(edge_list_s[7]) && (edge_list_s[0]) && (edge_list_s[9])  &&
					(edge_list_v[3]) ) {
				//printf("Entrou no 5!\n");

				cut_edge[0] = 5;
				cut_edge[1] = 9;
				cut_edge[2] = 8;

				cut_edge_s[0] = 7;
				cut_edge_s[1] = 0;
				cut_edge_s[2] = 9;

				cut_edge_v[0] = 3;

				connec_order_in1[0] = 0;
				connec_order_in1[1] = 1;
				connec_order_in1[2] = 2;
				connec_order_in1[3] = 3;
				connec_order_in1[4] = 4;
				connec_order_in1[5] = 5;
				connec_order_in1[6] = 6;
				connec_order_in1[7] = 7;

				connec_order_out1[0] = 0;
				connec_order_out1[1] = 1;
				connec_order_out1[2] = 2;
				connec_order_out1[3] = 3;
				connec_order_out1[4] = 4;
				connec_order_out1[5] = 0;
				connec_order_out1[6] = 3;
				connec_order_out1[7] = 6;

				connec_order_in2[0] = 0;
				connec_order_in2[1] = 1;
				connec_order_in2[2] = 4;
				connec_order_in2[3] = 5;
				connec_order_in2[4] = 2;
				connec_order_in2[5] = 3;
				connec_order_in2[6] = 6;
				connec_order_in2[7] = 7;

				connec_order_out2[0] = 6;
				connec_order_out2[1] = 3;
				connec_order_out2[2] = 5;
				connec_order_out2[3] = 1;
				connec_order_out2[4] = 2;
				connec_order_out2[5] = 3;
				connec_order_out2[6] = 6;
				connec_order_out2[7] = 7;

				connec_order_in3[0] = 1;
				connec_order_in3[1] = 2;
				connec_order_in3[2] = 5;
				connec_order_in3[3] = 6;
				connec_order_in3[4] = 0;
				connec_order_in3[5] = 3;
				connec_order_in3[6] = 4;
				connec_order_in3[7] = 7;

				connec_order_out3[0] = 4;
				connec_order_out3[1] = 6;
				connec_order_out3[2] = 2;
				connec_order_out3[3] = 5;
				connec_order_out3[4] = 0;
				connec_order_out3[5] = 3;
				connec_order_out3[6] = 4;
				connec_order_out3[7] = 7;

				connec_order_in4[0] = 5;
				connec_order_in4[1] = 0;
				connec_order_in4[2] = 1;
				connec_order_in4[3] = 2;
				connec_order_in4[4] = 3;
				connec_order_in4[5] = 4;
				connec_order_in4[6] = 6;
				connec_order_in4[7] = 7;

				connec_order_out4[0] = 5;
				connec_order_out4[1] = 4;
				connec_order_out4[2] = 0;
				connec_order_out4[3] = 3;
				connec_order_out4[4] = 6;
				connec_order_out4[5] = 2;
				connec_order_out4[6] = 1;
				connec_order_out4[7] = 5;

				elem->cut = 3;

			}

			//Corner 6
			else if ((!edge_list[0]) && (!edge_list[1]) && (!edge_list[2]) && (!edge_list[3]) &&
					(!edge_list[4]) && (!edge_list[5]) && (edge_list[6]) && (!edge_list[7]) &&
					(!edge_list[8]) && (edge_list[9]) && (edge_list[10]) && (!edge_list[11]) &&
					(edge_list_s[2]) && (edge_list_s[6]) && (edge_list_s[8])  &&
					(edge_list_v[0]) ) {
				//printf("Entrou no 6!\n");


				cut_edge[0] = 6;
				cut_edge[1] = 10;
				cut_edge[2] = 9;

				cut_edge_s[0] = 2;
				cut_edge_s[1] = 6;
				cut_edge_s[2] = 8;

				cut_edge_v[0] = 0;

				connec_order_in1[0] = 0;
				connec_order_in1[1] = 1;
				connec_order_in1[2] = 2;
				connec_order_in1[3] = 3;
				connec_order_in1[4] = 4;
				connec_order_in1[5] = 5;
				connec_order_in1[6] = 6;
				connec_order_in1[7] = 7;

				connec_order_out1[0] = 0;
				connec_order_out1[1] = 1;
				connec_order_out1[2] = 2;
				connec_order_out1[3] = 3;
				connec_order_out1[4] = 6;
				connec_order_out1[5] = 4;
				connec_order_out1[6] = 0;
				connec_order_out1[7] = 3;

				connec_order_in2[0] = 1;
				connec_order_in2[1] = 2;
				connec_order_in2[2] = 5;
				connec_order_in2[3] = 6;
				connec_order_in2[4] = 0;
				connec_order_in2[5] = 3;
				connec_order_in2[6] = 4;
				connec_order_in2[7] = 7;

				connec_order_out2[0] = 6;
				connec_order_out2[1] = 3;
				connec_order_out2[2] = 5;
				connec_order_out2[3] = 1;
				connec_order_out2[4] = 0;
				connec_order_out2[5] = 3;
				connec_order_out2[6] = 4;
				connec_order_out2[7] = 7;

				connec_order_in3[0] = 2;
				connec_order_in3[1] = 3;
				connec_order_in3[2] = 6;
				connec_order_in3[3] = 7;
				connec_order_in3[4] = 0;
				connec_order_in3[5] = 1;
				connec_order_in3[6] = 4;
				connec_order_in3[7] = 5;

				connec_order_out3[0] = 4;
				connec_order_out3[1] = 6;
				connec_order_out3[2] = 2;
				connec_order_out3[3] = 5;
				connec_order_out3[4] = 0;
				connec_order_out3[5] = 1;
				connec_order_out3[6] = 4;
				connec_order_out3[7] = 5;

				connec_order_in4[0] = 6;
				connec_order_in4[1] = 0;
				connec_order_in4[2] = 1;
				connec_order_in4[3] = 2;
				connec_order_in4[4] = 3;
				connec_order_in4[5] = 4;
				connec_order_in4[6] = 5;
				connec_order_in4[7] = 7;

				connec_order_out4[0] = 6;
				connec_order_out4[1] = 6;
				connec_order_out4[2] = 4;
				connec_order_out4[3] = 0;
				connec_order_out4[4] = 3;
				connec_order_out4[5] = 5;
				connec_order_out4[6] = 2;
				connec_order_out4[7] = 1;

				elem->cut = 3;

			}

			//Corner 7
			else if ((!edge_list[0]) && (!edge_list[1]) && (!edge_list[2]) && (!edge_list[3]) &&
					(!edge_list[4]) && (!edge_list[5]) && (!edge_list[6]) && (edge_list[7]) &&
					(!edge_list[8]) && (!edge_list[9]) && (edge_list[10]) && (edge_list[11]) &&
					(edge_list_s[4]) && (edge_list_s[3]) && (edge_list_s[9])  &&
					(edge_list_v[1]) ) {
				//printf("Entrou no 7!\n");


				cut_edge[0] = 7;
				cut_edge[1] = 11;
				cut_edge[2] = 10;

				cut_edge_s[0] = 4;
				cut_edge_s[1] = 3;
				cut_edge_s[2] = 9;

				cut_edge_v[0] = 1;

				connec_order_in1[0] = 0;
				connec_order_in1[1] = 1;
				connec_order_in1[2] = 2;
				connec_order_in1[3] = 3;
				connec_order_in1[4] = 4;
				connec_order_in1[5] = 5;
				connec_order_in1[6] = 6;
				connec_order_in1[7] = 7;

				connec_order_out1[0] = 0;
				connec_order_out1[1] = 1;
				connec_order_out1[2] = 2;
				connec_order_out1[3] = 3;
				connec_order_out1[4] = 3;
				connec_order_out1[5] = 6;
				connec_order_out1[6] = 4;
				connec_order_out1[7] = 0;

				connec_order_in2[0] = 2;
				connec_order_in2[1] = 3;
				connec_order_in2[2] = 6;
				connec_order_in2[3] = 7;
				connec_order_in2[4] = 0;
				connec_order_in2[5] = 1;
				connec_order_in2[6] = 4;
				connec_order_in2[7] = 5;

				connec_order_out2[0] = 6;
				connec_order_out2[1] = 3;
				connec_order_out2[2] = 5;
				connec_order_out2[3] = 1;
				connec_order_out2[4] = 0;
				connec_order_out2[5] = 1;
				connec_order_out2[6] = 4;
				connec_order_out2[7] = 5;

				connec_order_in3[0] = 0;
				connec_order_in3[1] = 3;
				connec_order_in3[2] = 4;
				connec_order_in3[3] = 7;
				connec_order_in3[4] = 1;
				connec_order_in3[5] = 2;
				connec_order_in3[6] = 5;
				connec_order_in3[7] = 6;

				connec_order_out3[0] = 6;
				connec_order_out3[1] = 4;
				connec_order_out3[2] = 5;
				connec_order_out3[3] = 2;
				connec_order_out3[4] = 1;
				connec_order_out3[5] = 2;
				connec_order_out3[6] = 5;
				connec_order_out3[7] = 6;

				connec_order_in4[0] = 7;
				connec_order_in4[1] = 0;
				connec_order_in4[2] = 1;
				connec_order_in4[3] = 2;
				connec_order_in4[4] = 3;
				connec_order_in4[5] = 4;
				connec_order_in4[6] = 5;
				connec_order_in4[7] = 6;

				connec_order_out4[0] = 7;
				connec_order_out4[1] = 3;
				connec_order_out4[2] = 6;
				connec_order_out4[3] = 4;
				connec_order_out4[4] = 0;
				connec_order_out4[5] = 1;
				connec_order_out4[6] = 5;
				connec_order_out4[7] = 2;

				elem->cut = 3;

			} else {
				elem->cut = -1;
				printf("Warning: Bad Element: %d, case 3\n",elements_ids[iel]);
				printf("Edges:  Elem. %d: %d %d %d %d %d %d %d %d %d %d %d %d %d\n", elements_ids[iel], edge_list[0],
						edge_list[1], edge_list[2], edge_list[3],
						edge_list[4], edge_list[5], edge_list[6],
						edge_list[7], edge_list[8], edge_list[9],
						edge_list[10], edge_list[11], ed_cont);
				continue;
			}

			GtsPoint *p0 = point[cut_edge[0]];
			GtsPoint *p1 = point[cut_edge[1]];
			GtsPoint *p2 = point[cut_edge[2]];

			GtsPoint *p3 = point_s[cut_edge_s[0]];
			GtsPoint *p4 = point_s[cut_edge_s[1]];
			GtsPoint *p5 = point_s[cut_edge_s[2]];

			GtsPoint *p6 = point_v[cut_edge_v[0]];

			//g_assert(p0 != NULL);
			//g_assert(p1 != NULL);
			//g_assert(p2 != NULL);
			//g_assert(p3 != NULL);
			//g_assert(p4 != NULL);
			//g_assert(p5 != NULL);
			//g_assert(p6 != NULL);

			if(p0==NULL) exit(1);
			if(p1==NULL) exit(2);
			if(p2==NULL) exit(3);
			if(p3==NULL) exit(4);
			if(p4==NULL) exit(5);
			if(p5==NULL) exit(6);
			if(p6==NULL) exit(7);

			conn_t2[0] = AddPointOnEdge(Edge2GNode[cut_edge[0]], hash_nodes, mesh->local_n_nodes, p0, coords);
			conn_t2[1] = AddPointOnEdge(Edge2GNode[cut_edge[1]], hash_nodes, mesh->local_n_nodes, p1, coords);
			conn_t2[2] = AddPointOnEdge(Edge2GNode[cut_edge[2]], hash_nodes, mesh->local_n_nodes, p2, coords);

			// add 3 extra points in the surface
			conn_t2[3] = AddPointOnEdge(Edge2GNode_s[cut_edge_s[0]], hash_nodes, mesh->local_n_nodes, p3, coords);
			conn_t2[4] = AddPointOnEdge(Edge2GNode_s[cut_edge_s[1]], hash_nodes, mesh->local_n_nodes, p4, coords);
			conn_t2[5] = AddPointOnEdge(Edge2GNode_s[cut_edge_s[2]], hash_nodes, mesh->local_n_nodes, p5, coords);

			// add 1 extra points in the volume
			conn_t2[6] = AddPointOnEdge(Edge2GNode_s[cut_edge_v[0]], hash_nodes, mesh->local_n_nodes, p6, coords);

			octant_t *elem1 = (octant_t*) sc_array_index(&mesh->elements, elements_ids[iel]);

			elem1->nodes[connec_order_in1[0]].id = original_conn[connec_order_out1[0]];
			elem1->nodes[connec_order_in1[1]].id = original_conn[connec_order_out1[1]];
			elem1->nodes[connec_order_in1[2]].id = original_conn[connec_order_out1[2]];
			elem1->nodes[connec_order_in1[3]].id = original_conn[connec_order_out1[3]];

			elem1->nodes[connec_order_in1[4]].id = conn_t2[connec_order_out1[4]];
			elem1->nodes[connec_order_in1[5]].id = conn_t2[connec_order_out1[5]];
			elem1->nodes[connec_order_in1[6]].id = conn_t2[connec_order_out1[6]];
			elem1->nodes[connec_order_in1[7]].id = conn_t2[connec_order_out1[7]];

			octant_t* elem2 = (octant_t*) sc_array_push(&mesh->elements);

			elem2->nodes[connec_order_in2[0]].id = conn_t2[connec_order_out2[0]];
			elem2->nodes[connec_order_in2[1]].id = conn_t2[connec_order_out2[1]];
			elem2->nodes[connec_order_in2[2]].id = conn_t2[connec_order_out2[2]];
			elem2->nodes[connec_order_in2[3]].id = conn_t2[connec_order_out2[3]];

			elem2->nodes[connec_order_in2[4]].id = original_conn[connec_order_out2[4]];
			elem2->nodes[connec_order_in2[5]].id = original_conn[connec_order_out2[5]];
			elem2->nodes[connec_order_in2[6]].id = original_conn[connec_order_out2[6]];
			elem2->nodes[connec_order_in2[7]].id = original_conn[connec_order_out2[7]];

			elem2->cut = elem->cut;
			elem2->level = elem->level;

			octant_t* elem3 = (octant_t*) sc_array_push(&mesh->elements);

			elem3->nodes[connec_order_in3[0]].id = conn_t2[connec_order_out3[0]];
			elem3->nodes[connec_order_in3[1]].id = conn_t2[connec_order_out3[1]];
			elem3->nodes[connec_order_in3[2]].id = conn_t2[connec_order_out3[2]];
			elem3->nodes[connec_order_in3[3]].id = conn_t2[connec_order_out3[3]];

			elem3->nodes[connec_order_in3[4]].id = original_conn[connec_order_out3[4]];
			elem3->nodes[connec_order_in3[5]].id = original_conn[connec_order_out3[5]];
			elem3->nodes[connec_order_in3[6]].id = original_conn[connec_order_out3[6]];
			elem3->nodes[connec_order_in3[7]].id = original_conn[connec_order_out3[7]];

			elem3->cut = elem->cut;
			elem3->level = elem->level;

			octant_t* elem4 = (octant_t*) sc_array_push(&mesh->elements);

			elem4->nodes[connec_order_in4[0]].id = original_conn[connec_order_out4[0]];

			elem4->nodes[connec_order_in4[1]].id = conn_t2[connec_order_out4[1]];
			elem4->nodes[connec_order_in4[2]].id = conn_t2[connec_order_out4[2]];
			elem4->nodes[connec_order_in4[3]].id = conn_t2[connec_order_out4[3]];
			elem4->nodes[connec_order_in4[4]].id = conn_t2[connec_order_out4[4]];
			elem4->nodes[connec_order_in4[5]].id = conn_t2[connec_order_out4[5]];
			elem4->nodes[connec_order_in4[6]].id = conn_t2[connec_order_out4[6]];
			elem4->nodes[connec_order_in4[7]].id = conn_t2[connec_order_out4[7]];

			elem4->cut = elem->cut;
			elem4->level = elem->level;


		}

		for (int edge = 0; edge < 4; edge++) {
			if (point_v[edge]) gts_object_destroy(GTS_OBJECT(point_v[edge]));
			point_v[edge] = NULL;
		}

		for (int edge = 0; edge < 12; edge++) {
			if (point_s[edge]) gts_object_destroy(GTS_OBJECT(point_s[edge]));
			point_s[edge] = NULL;
		}

		for (int edge = 0; edge < 12; edge++) {
			if (point[edge]) gts_object_destroy(GTS_OBJECT(point[edge]));
			point[edge] = NULL;
		}

	}

	mesh->local_n_elements = mesh->elements.elem_count;
	//mesh->local_n_nodes = coords.size()/3;

	//MPI_Allreduce(&mesh->local_n_elements, &mesh->total_n_elements, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
	//MPI_Allreduce(&mesh->local_n_nodes, &mesh->total_n_nodes, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);

	//	if (mesh->mpi_rank == 0) {
	//		printf("Total number of elements: %lld\n", mesh->local_n_elements);
	//		printf("Total number of nodes: %lld\n", mesh->local_n_nodes);
	//	}

	fclose(fdbg);
	sc_hash_array_destroy(hash_nodes);
}
