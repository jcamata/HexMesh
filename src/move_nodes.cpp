#include <gts.h>
#include <glib.h>
#include <vector>
#include <iostream>
using namespace std;
#include <set>
#include <algorithm>
#include <sc.h>
#include <sc_io.h>
#include <sc_containers.h>
#include <mpi.h>

#include "hexa.h"
#include "hilbert.h"
#include "refinement.h"

#include <ctime>


typedef struct {
	GPtrArray * array;
} ListOfPoints;

static void InsertPoint(GtsPoint * p, ListOfPoints * lp) {
	g_ptr_array_add(lp->array, p);
}

bool is_point_over_surface(GtsPoint * p, GNode * tree);

unsigned vertex_hash_id(const void *v, const void *u) {
	const octant_vertex_t *q = (const octant_vertex_t*) v;
	uint32_t a, b, c;

	a = (uint32_t) q->id;
	b = (uint32_t) 0;
	c = (uint32_t) 1;
	sc_hash_mix(a, b, c);
	sc_hash_final(a, b, c);
	return (unsigned) c;
}

int vertex_equal_id(const void *v, const void *u, const void *w) {
	const octant_vertex_t *e1 = (const octant_vertex_t*) v;
	const octant_vertex_t *e2 = (const octant_vertex_t*) u;

	return (unsigned) (e1->id==e2->id);
}

unsigned el_hash_id(const void *v, const void *u) {
	const octant_t *q = (const octant_t*) v;
	uint64_t a, b, c;

	a = (uint32_t) q->id;
	b = (uint32_t) 0;
	c = (uint32_t) 0;
	sc_hash_mix(a, b, c);
	sc_hash_final(a, b, c);
	return (unsigned) c;
}

int el_equal_id(const void *v, const void *u, const void *w) {
	const octant_t *e1 = (const octant_t*) v;
	const octant_t *e2 = (const octant_t*) u;

	return (unsigned) (e1->id == e2->id);

}

typedef struct {
	bitmask_t coord[3];
	int id;
} node_in_edge_t;

unsigned no_hash_fn1(const void *v, const void *u) {
	const node_in_edge_t *q = (const node_in_edge_t*) v;
	uint32_t a, b, c;

	a = (uint32_t) q->id;
	b = (uint32_t) 0;
	c = (uint32_t) 0;
	sc_hash_mix(a, b, c);
	sc_hash_final(a, b, c);
	return (unsigned) c;
}

int no_equal_fn1(const void *v, const void *u, const void *w) {
	const node_in_edge_t *e1 = (const node_in_edge_t*) v;
	const node_in_edge_t *e2 = (const node_in_edge_t*) u;

	return (unsigned) ((e1->id == e2->id));

}

void ProjectFreeNodes(hexa_tree_t* mesh,std::vector<double>& coords, std::vector<int>& nodes_b_mat){

	int Edge2GNode[12][2]={0};
	int Edge2GNode_s[12][2]={0};
	GtsSegment * segments[12]={0};
	GtsSegment * segments_s[12]={0};
	GtsPoint * point[12]={NULL};
	GtsPoint * point_s[12]={NULL};
	bool clamped = true;
	std::vector<double> aux;
	//int count = 0;

	//moving the nodes in the surface...
	if(true){
		//achando os pontos de onde a superficie corta o octree nas 6 superficies
		for (int ioc = 0; ioc < mesh->oct.elem_count; ++ioc) {
			octree_t* oct = (octree_t*)sc_array_index(&mesh->oct,ioc);

			int oc_count=0;
			for(int i =0;i<8;i++){
				if(oct->id[i]!=-1) {
					oc_count++;
				}
			}

			if(oc_count==8){
				for(int iel = 0; iel<8; iel++){
					octant_t* elem = (octant_t*)sc_array_index(&mesh->elements,oct->id[iel]);

					//verifica se as arestas foram cortadas
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
						if (list == NULL) continue;
						while (list) {
							GtsBBox *b = GTS_BBOX(list->data);
							point[edge] = SegmentTriangleIntersection(segments[edge], GTS_TRIANGLE(b->bounded));
							if (point[edge]) {
								break;
							}
							list = list->next;
						}
					}
					//verifica se as diagonais das faces foram cortadas
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
						if (list == NULL) continue;
						while (list) {
							GtsBBox *b = GTS_BBOX(list->data);
							point_s[edge] = SegmentTriangleIntersection(segments_s[edge], GTS_TRIANGLE(b->bounded));
							if (point_s[edge]) {
								break;
							}
							list = list->next;
						}
					}

					//teoricamente essa parte move os nos das faces :)
					//face x-
					if(oct->face[0]){
						bool edge_cut = false;
						if(iel==0){
							if(point[7]!=NULL){
								aux.push_back(elem->nodes[7].id);
								aux.push_back(point[7]->x);
								aux.push_back(point[7]->y);
								aux.push_back(point[7]->z);
								edge_cut = true;
							}else if(point[11]!=NULL){
								aux.push_back(elem->nodes[7].id);
								aux.push_back(point[11]->x);
								aux.push_back(point[11]->y);
								aux.push_back(point[11]->z);
								edge_cut = true;
							}
						}

						if(iel==3){
							if(point[4]!=NULL){
								aux.push_back(elem->nodes[4].id);
								aux.push_back(point[4]->x);
								aux.push_back(point[4]->y);
								aux.push_back(point[4]->z);
								edge_cut = true;
							}
							if(point[11]!=NULL){
								aux.push_back(elem->nodes[7].id);
								aux.push_back(point[11]->x);
								aux.push_back(point[11]->y);
								aux.push_back(point[11]->z);
								edge_cut = true;
							}
						}

						if(iel==4){
							if(point[7]!=NULL){
								aux.push_back(elem->nodes[3].id);
								aux.push_back(point[7]->x);
								aux.push_back(point[7]->y);
								aux.push_back(point[7]->z);
								edge_cut = true;
							}
							if(point[3]!=NULL){
								aux.push_back(elem->nodes[3].id);
								aux.push_back(point[3]->x);
								aux.push_back(point[3]->y);
								aux.push_back(point[3]->z);
								edge_cut = true;
							}
						}

						if(iel==7){
							if(point[4]!=NULL){
								aux.push_back(elem->nodes[0].id);
								aux.push_back(point[4]->x);
								aux.push_back(point[4]->y);
								aux.push_back(point[4]->z);
								edge_cut = true;
							}
							if(point[3]!=NULL){
								aux.push_back(elem->nodes[0].id);
								aux.push_back(point[3]->x);
								aux.push_back(point[3]->y);
								aux.push_back(point[3]->z);
								edge_cut = true;
							}

							if(!edge_cut){
								if(point_s[4]!=NULL){
									if(iel==0){
										aux.push_back(elem->nodes[7].id);
										aux.push_back(point_s[4]->x);
										aux.push_back(point_s[4]->y);
										aux.push_back(point_s[4]->z);
									}
									if(iel==7){
										aux.push_back(elem->nodes[0].id);
										aux.push_back(point_s[4]->x);
										aux.push_back(point_s[4]->y);
										aux.push_back(point_s[4]->z);
									}
								}
								if(point_s[5]!=NULL){
									if(iel==3){
										aux.push_back(elem->nodes[4].id);
										aux.push_back(point_s[5]->x);
										aux.push_back(point_s[5]->y);
										aux.push_back(point_s[5]->z);
									}
									if(iel==4){
										aux.push_back(elem->nodes[3].id);
										aux.push_back(point_s[5]->x);
										aux.push_back(point_s[5]->y);
										aux.push_back(point_s[5]->z);
									}
								}
							}
						}

					}

					//face x+
					if(oct->face[1]){
						bool edge_cut = false;

						if(iel==1){
							if(point[6]!=NULL){
								aux.push_back(elem->nodes[6].id);
								aux.push_back(point[6]->x);
								aux.push_back(point[6]->y);
								aux.push_back(point[6]->z);
								edge_cut = true;

							}else if(point[9]!=NULL){
								aux.push_back(elem->nodes[6].id);
								aux.push_back(point[9]->x);
								aux.push_back(point[9]->y);
								aux.push_back(point[9]->z);
								edge_cut = true;

							}
						}


						if(iel==2){
							if(point[5]!=NULL){
								aux.push_back(elem->nodes[5].id);
								aux.push_back(point[5]->x);
								aux.push_back(point[5]->y);
								aux.push_back(point[5]->z);
								edge_cut = true;

							}
							if(point[9]!=NULL){
								aux.push_back(elem->nodes[5].id);
								aux.push_back(point[9]->x);
								aux.push_back(point[9]->y);
								aux.push_back(point[9]->z);
								edge_cut = true;

							}
						}

						if(iel==5){
							if(point[1]!=NULL){
								aux.push_back(elem->nodes[2].id);
								aux.push_back(point[1]->x);
								aux.push_back(point[1]->y);
								aux.push_back(point[1]->z);
								edge_cut = true;

							}
							if(point[6]!=NULL){
								aux.push_back(elem->nodes[2].id);
								aux.push_back(point[6]->x);
								aux.push_back(point[6]->y);
								aux.push_back(point[6]->z);
								edge_cut = true;

							}
						}

						if(iel==6){
							if(point[1]!=NULL){
								aux.push_back(elem->nodes[1].id);
								aux.push_back(point[1]->x);
								aux.push_back(point[1]->y);
								aux.push_back(point[1]->z);
								edge_cut = true;

							}
							if(point[5]!=NULL){
								aux.push_back(elem->nodes[1].id);
								aux.push_back(point[5]->x);
								aux.push_back(point[5]->y);
								aux.push_back(point[5]->z);
								edge_cut = true;

							}
						}
						if(!edge_cut){
							if(point_s[6]!=NULL){
								if(iel==1){
									aux.push_back(elem->nodes[6].id);
									aux.push_back(point_s[6]->x);
									aux.push_back(point_s[6]->y);
									aux.push_back(point_s[6]->z);
								}
								if(iel==6){
									aux.push_back(elem->nodes[1].id);
									aux.push_back(point_s[6]->x);
									aux.push_back(point_s[6]->y);
									aux.push_back(point_s[6]->z);
								}
							}
							if(point_s[7]!=NULL){
								if(iel==2){
									aux.push_back(elem->nodes[5].id);
									aux.push_back(point_s[7]->x);
									aux.push_back(point_s[7]->y);
									aux.push_back(point_s[7]->z);

									if(iel==5){
										aux.push_back(elem->nodes[2].id);
										aux.push_back(point_s[7]->x);
										aux.push_back(point_s[7]->y);
										aux.push_back(point_s[7]->z);
									}
								}
							}
						}
					}

					//face y-
					if(oct->face[2]){
						bool edge_cut = false;

						if(iel==0){
							if(point[5]!=NULL){
								aux.push_back(elem->nodes[5].id);
								aux.push_back(point[5]->x);
								aux.push_back(point[5]->y);
								aux.push_back(point[5]->z);
								edge_cut = true;

							}else if(point[8]!=NULL){
								aux.push_back(elem->nodes[5].id);
								aux.push_back(point[8]->x);
								aux.push_back(point[8]->y);
								aux.push_back(point[8]->z);
								edge_cut = true;

							}
						}

						if(iel==1){
							if(point[4]!=NULL){
								aux.push_back(elem->nodes[4].id);
								aux.push_back(point[4]->x);
								aux.push_back(point[4]->y);
								aux.push_back(point[4]->z);
								edge_cut = true;

							}
							if(point[8]!=NULL){
								aux.push_back(elem->nodes[4].id);
								aux.push_back(point[8]->x);
								aux.push_back(point[8]->y);
								aux.push_back(point[8]->z);
								edge_cut = true;

							}
						}

						if(iel==4){
							if(point[5]!=NULL){
								aux.push_back(elem->nodes[1].id);
								aux.push_back(point[5]->x);
								aux.push_back(point[5]->y);
								aux.push_back(point[5]->z);
								edge_cut = true;

							}
							if(point[0]!=NULL){
								aux.push_back(elem->nodes[1].id);
								aux.push_back(point[0]->x);
								aux.push_back(point[0]->y);
								aux.push_back(point[0]->z);
								edge_cut = true;

							}
						}

						if(iel==5){
							if(point[4]!=NULL){
								aux.push_back(elem->nodes[0].id);
								aux.push_back(point[4]->x);
								aux.push_back(point[4]->y);
								aux.push_back(point[4]->z);
								edge_cut = true;

							}
							if(point[0]!=NULL){
								aux.push_back(elem->nodes[0].id);
								aux.push_back(point[0]->x);
								aux.push_back(point[0]->y);
								aux.push_back(point[0]->z);
								edge_cut = true;

							}
						}
						if(!edge_cut){
							if(point_s[0]!=NULL){
								if(iel==0){
									aux.push_back(elem->nodes[5].id);
									aux.push_back(point_s[0]->x);
									aux.push_back(point_s[0]->y);
									aux.push_back(point_s[0]->z);
								}
								if(iel==5){
									aux.push_back(elem->nodes[0].id);
									aux.push_back(point_s[0]->x);
									aux.push_back(point_s[0]->y);
									aux.push_back(point_s[0]->z);
								}
							}
							if(point_s[1]!=NULL){
								if(iel==1){
									aux.push_back(elem->nodes[4].id);
									aux.push_back(point_s[1]->x);
									aux.push_back(point_s[1]->y);
									aux.push_back(point_s[1]->z);
								}
								if(iel==4){
									aux.push_back(elem->nodes[1].id);
									aux.push_back(point_s[1]->x);
									aux.push_back(point_s[1]->y);
									aux.push_back(point_s[1]->z);
								}
							}
						}
					}

					//face y+
					if(oct->face[3]){
						bool edge_cut = false;

						if(iel==2){
							if(point[7]!=NULL){
								aux.push_back(elem->nodes[7].id);
								aux.push_back(point[7]->x);
								aux.push_back(point[7]->y);
								aux.push_back(point[7]->z);
								edge_cut = true;

							}else if(point[10]!=NULL){
								aux.push_back(elem->nodes[7].id);
								aux.push_back(point[10]->x);
								aux.push_back(point[10]->y);
								aux.push_back(point[10]->z);
								edge_cut = true;

							}
						}


						if(iel==3){
							if(point[6]!=NULL){
								aux.push_back(elem->nodes[6].id);
								aux.push_back(point[6]->x);
								aux.push_back(point[6]->y);
								aux.push_back(point[6]->z);
								edge_cut = true;

							}
							if(point[10]!=NULL){
								aux.push_back(elem->nodes[6].id);
								aux.push_back(point[10]->x);
								aux.push_back(point[10]->y);
								aux.push_back(point[10]->z);
								edge_cut = true;

							}
						}

						if(iel==6){
							if(point[7]!=NULL){
								aux.push_back(elem->nodes[3].id);
								aux.push_back(point[7]->x);
								aux.push_back(point[7]->y);
								aux.push_back(point[7]->z);
								edge_cut = true;

							}
							if(point[3]!=NULL){
								aux.push_back(elem->nodes[3].id);
								aux.push_back(point[3]->x);
								aux.push_back(point[3]->y);
								aux.push_back(point[3]->z);
								edge_cut = true;

							}
						}

						if(iel==7){
							if(point[6]!=NULL){
								aux.push_back(elem->nodes[2].id);
								aux.push_back(point[6]->x);
								aux.push_back(point[6]->y);
								aux.push_back(point[6]->z);
								edge_cut = true;

							}
							if(point[3]!=NULL){
								aux.push_back(elem->nodes[2].id);
								aux.push_back(point[3]->x);
								aux.push_back(point[3]->y);
								aux.push_back(point[3]->z);
								edge_cut = true;

							}
						}
						if(!edge_cut){
							if(point_s[3]!=NULL){
								if(iel==2){
									aux.push_back(elem->nodes[7].id);
									aux.push_back(point_s[3]->x);
									aux.push_back(point_s[3]->y);
									aux.push_back(point_s[3]->z);
								}
								if(iel==7){
									aux.push_back(elem->nodes[2].id);
									aux.push_back(point_s[3]->x);
									aux.push_back(point_s[3]->y);
									aux.push_back(point_s[3]->z);
								}
							}
							if(point_s[2]!=NULL){
								if(iel==3){
									aux.push_back(elem->nodes[6].id);
									aux.push_back(point_s[2]->x);
									aux.push_back(point_s[2]->y);
									aux.push_back(point_s[2]->z);
								}
								if(iel==6){
									aux.push_back(elem->nodes[3].id);
									aux.push_back(point_s[2]->x);
									aux.push_back(point_s[2]->y);
									aux.push_back(point_s[2]->z);
								}
							}
						}
					}

					//face z-
					if(oct->face[4]){
						bool edge_cut = false;

						if(iel==4){
							if(point[9]!=NULL){
								aux.push_back(elem->nodes[6].id);
								aux.push_back(point[9]->x);
								aux.push_back(point[9]->y);
								aux.push_back(point[9]->z);
								edge_cut = true;

							}else if(point[10]!=NULL){
								aux.push_back(elem->nodes[6].id);
								aux.push_back(point[10]->x);
								aux.push_back(point[10]->y);
								aux.push_back(point[10]->z);
								edge_cut = true;

							}
						}

						if(iel==5){
							if(point[11]!=NULL){
								aux.push_back(elem->nodes[7].id);
								aux.push_back(point[11]->x);
								aux.push_back(point[11]->y);
								aux.push_back(point[11]->z);
								edge_cut = true;

							}
							if(point[10]!=NULL){
								aux.push_back(elem->nodes[7].id);
								aux.push_back(point[10]->x);
								aux.push_back(point[10]->y);
								aux.push_back(point[10]->z);
								edge_cut = true;

							}
						}

						if(iel==6){
							if(point[8]!=NULL){
								aux.push_back(elem->nodes[4].id);
								aux.push_back(point[8]->x);
								aux.push_back(point[8]->y);
								aux.push_back(point[8]->z);
								edge_cut = true;

							}
							if(point[11]!=NULL){
								aux.push_back(elem->nodes[4].id);
								aux.push_back(point[11]->x);
								aux.push_back(point[11]->y);
								aux.push_back(point[11]->z);
								edge_cut = true;

							}
						}

						if(iel==7){
							if(point[8]!=NULL){
								aux.push_back(elem->nodes[5].id);
								aux.push_back(point[8]->x);
								aux.push_back(point[8]->y);
								aux.push_back(point[8]->z);
								edge_cut = true;

							}
							if(point[9]!=NULL){
								aux.push_back(elem->nodes[5].id);
								aux.push_back(point[9]->x);
								aux.push_back(point[9]->y);
								aux.push_back(point[9]->z);
								edge_cut = true;

							}
						}
						if(!edge_cut){
							if(point_s[8]!=NULL){
								if(iel==4){
									aux.push_back(elem->nodes[6].id);
									aux.push_back(point_s[8]->x);
									aux.push_back(point_s[8]->y);
									aux.push_back(point_s[8]->z);
								}
								if(iel==6){
									aux.push_back(elem->nodes[4].id);
									aux.push_back(point_s[8]->x);
									aux.push_back(point_s[8]->y);
									aux.push_back(point_s[8]->z);
								}
							}
							if(point_s[9]!=NULL){
								if(iel==5){
									aux.push_back(elem->nodes[7].id);
									aux.push_back(point_s[9]->x);
									aux.push_back(point_s[9]->y);
									aux.push_back(point_s[9]->z);
								}
								if(iel==7){
									aux.push_back(elem->nodes[5].id);
									aux.push_back(point_s[9]->x);
									aux.push_back(point_s[9]->y);
									aux.push_back(point_s[9]->z);
								}
							}
						}
					}

					//face z+
					if(oct->face[5]){
						bool edge_cut = false;

						if(iel==0){
							if(point[1]!=NULL){
								aux.push_back(elem->nodes[2].id);
								aux.push_back(point[1]->x);
								aux.push_back(point[1]->y);
								aux.push_back(point[1]->z);
								edge_cut = true;

							}else if(point[2]!=NULL){
								aux.push_back(elem->nodes[2].id);
								aux.push_back(point[2]->x);
								aux.push_back(point[2]->y);
								aux.push_back(point[2]->z);
								edge_cut = true;

							}
						}

						if(iel==1){
							if(point[3]!=NULL){
								aux.push_back(elem->nodes[3].id);
								aux.push_back(point[3]->x);
								aux.push_back(point[3]->y);
								aux.push_back(point[3]->z);
								edge_cut = true;

							}
							if(point[2]!=NULL){
								aux.push_back(elem->nodes[3].id);
								aux.push_back(point[2]->x);
								aux.push_back(point[2]->y);
								aux.push_back(point[2]->z);
								edge_cut = true;

							}
						}

						if(iel==2){
							if(point[0]!=NULL){
								aux.push_back(elem->nodes[0].id);
								aux.push_back(point[0]->x);
								aux.push_back(point[0]->y);
								aux.push_back(point[0]->z);
								edge_cut = true;

							}
							if(point[3]!=NULL){
								aux.push_back(elem->nodes[0].id);
								aux.push_back(point[3]->x);
								aux.push_back(point[3]->y);
								aux.push_back(point[3]->z);
								edge_cut = true;

							}
						}

						if(iel==3){
							if(point[0]!=NULL){
								aux.push_back(elem->nodes[1].id);
								aux.push_back(point[0]->x);
								aux.push_back(point[0]->y);
								aux.push_back(point[0]->z);
								edge_cut = true;

							}
							if(point[1]!=NULL){
								aux.push_back(elem->nodes[1].id);
								aux.push_back(point[1]->x);
								aux.push_back(point[1]->y);
								aux.push_back(point[1]->z);
								edge_cut = true;

							}
						}
						if(!edge_cut){

							if(point_s[10]!=NULL){
								if(iel==0){
									aux.push_back(elem->nodes[2].id);
									aux.push_back(point_s[10]->x);
									aux.push_back(point_s[10]->y);
									aux.push_back(point_s[10]->z);
								}
								if(iel==2){
									aux.push_back(elem->nodes[0].id);
									aux.push_back(point_s[10]->x);
									aux.push_back(point_s[10]->y);
									aux.push_back(point_s[10]->z);
								}
							}
							if(point_s[11]!=NULL){
								if(iel==1){
									aux.push_back(elem->nodes[3].id);
									aux.push_back(point_s[11]->x);
									aux.push_back(point_s[11]->y);
									aux.push_back(point_s[11]->z);
								}
								if(iel==3){
									aux.push_back(elem->nodes[1].id);
									aux.push_back(point_s[11]->x);
									aux.push_back(point_s[11]->y);
									aux.push_back(point_s[11]->z);
								}
							}
						}
					}
				}
			}

			if(oc_count==4){
				int node1, node2;
				GtsVertex *v1;
				GtsVertex *v2;
				int edge = 0;
				GtsBBox *sb;
				GSList* list;
				GtsBBox *b;

				if(oct->id[3]!=-1 && oct->id[7]!=-1){
					//essa é para x...
					octant_t* elem0 = (octant_t*)sc_array_index(&mesh->elements,oct->id[0]);
					octant_t* elem3 = (octant_t*)sc_array_index(&mesh->elements,oct->id[3]);
					octant_t* elem4 = (octant_t*)sc_array_index(&mesh->elements,oct->id[4]);
					octant_t* elem7 = (octant_t*)sc_array_index(&mesh->elements,oct->id[7]);

					//face x-
					if(true){
						node1 = elem0->nodes[3].id;
						node2 = elem4->nodes[7].id;

						v1 = gts_vertex_new(gts_vertex_class(), coords[node1 * 3 + 0], coords[node1 * 3 + 1], coords[node1 * 3 + 2]);
						v2 = gts_vertex_new(gts_vertex_class(), coords[node2 * 3 + 0], coords[node2 * 3 + 1], coords[node2 * 3 + 2]);

						point[edge] = NULL;
						segments[edge] = gts_segment_new(gts_segment_class(), v1, v2);
						sb = gts_bbox_segment(gts_bbox_class(), segments[edge]);
						list = gts_bb_tree_overlap(mesh->gdata.bbt, sb);
						//if (list == NULL) continue;
						while (list) {
							b = GTS_BBOX(list->data);
							point[edge] = SegmentTriangleIntersection(segments[edge], GTS_TRIANGLE(b->bounded));
							if (point[edge]) {
								break;
							}
							list = list->next;
						}

						if(point[edge]!=NULL){
							aux.push_back(elem0->nodes[7].id);
							aux.push_back(point[edge]->x);
							aux.push_back(point[edge]->y);
							aux.push_back(point[edge]->z);
						}
					}

					//face x+
					if(true){
						node1 = elem0->nodes[2].id;
						node2 = elem4->nodes[6].id;

						v1 = gts_vertex_new(gts_vertex_class(), coords[node1 * 3 + 0], coords[node1 * 3 + 1], coords[node1 * 3 + 2]);
						v2 = gts_vertex_new(gts_vertex_class(), coords[node2 * 3 + 0], coords[node2 * 3 + 1], coords[node2 * 3 + 2]);

						point[edge] = NULL;
						segments[edge] = gts_segment_new(gts_segment_class(), v1, v2);
						sb = gts_bbox_segment(gts_bbox_class(), segments[edge]);
						list = gts_bb_tree_overlap(mesh->gdata.bbt, sb);
						//if (list == NULL) continue;
						while (list) {
							b = GTS_BBOX(list->data);
							point[edge] = SegmentTriangleIntersection(segments[edge], GTS_TRIANGLE(b->bounded));
							if (point[edge]) {
								break;
							}
							list = list->next;
						}

						if(point[edge]!=NULL){
							aux.push_back(elem0->nodes[6].id);
							aux.push_back(point[edge]->x);
							aux.push_back(point[edge]->y);
							aux.push_back(point[edge]->z);
						}
					}

					//face y-
					if(true){
						node1 = elem0->nodes[1].id;
						node2 = elem4->nodes[5].id;

						v1 = gts_vertex_new(gts_vertex_class(), coords[node1 * 3 + 0], coords[node1 * 3 + 1], coords[node1 * 3 + 2]);
						v2 = gts_vertex_new(gts_vertex_class(), coords[node2 * 3 + 0], coords[node2 * 3 + 1], coords[node2 * 3 + 2]);

						point[edge] = NULL;
						segments[edge] = gts_segment_new(gts_segment_class(), v1, v2);
						sb = gts_bbox_segment(gts_bbox_class(), segments[edge]);
						list = gts_bb_tree_overlap(mesh->gdata.bbt, sb);
						//if (list == NULL) continue;
						while (list) {
							b = GTS_BBOX(list->data);
							point[edge] = SegmentTriangleIntersection(segments[edge], GTS_TRIANGLE(b->bounded));
							if (point[edge]) {
								break;
							}
							list = list->next;
						}

						if(point[edge]!=NULL){
							aux.push_back(elem0->nodes[5].id);
							aux.push_back(point[edge]->x);
							aux.push_back(point[edge]->y);
							aux.push_back(point[edge]->z);
						}
					}

					//face y+
					if(true){
						node1 = elem3->nodes[2].id;
						node2 = elem7->nodes[6].id;

						v1 = gts_vertex_new(gts_vertex_class(), coords[node1 * 3 + 0], coords[node1 * 3 + 1], coords[node1 * 3 + 2]);
						v2 = gts_vertex_new(gts_vertex_class(), coords[node2 * 3 + 0], coords[node2 * 3 + 1], coords[node2 * 3 + 2]);

						point[edge] = NULL;
						segments[edge] = gts_segment_new(gts_segment_class(), v1, v2);
						sb = gts_bbox_segment(gts_bbox_class(), segments[edge]);
						list = gts_bb_tree_overlap(mesh->gdata.bbt, sb);
						if (list == NULL) continue;
						while (list) {
							b = GTS_BBOX(list->data);
							point[edge] = SegmentTriangleIntersection(segments[edge], GTS_TRIANGLE(b->bounded));
							if (point[edge]) {
								break;
							}
							list = list->next;
						}

						if(point[edge]!=NULL){
							aux.push_back(elem3->nodes[6].id);
							aux.push_back(point[edge]->x);
							aux.push_back(point[edge]->y);
							aux.push_back(point[edge]->z);
						}
					}

					//face z-
					if(true){
						node1 = elem4->nodes[5].id;
						node2 = elem7->nodes[6].id;

						v1 = gts_vertex_new(gts_vertex_class(), coords[node1 * 3 + 0], coords[node1 * 3 + 1], coords[node1 * 3 + 2]);
						v2 = gts_vertex_new(gts_vertex_class(), coords[node2 * 3 + 0], coords[node2 * 3 + 1], coords[node2 * 3 + 2]);

						point[edge] = NULL;
						segments[edge] = gts_segment_new(gts_segment_class(), v1, v2);
						sb = gts_bbox_segment(gts_bbox_class(), segments[edge]);
						list = gts_bb_tree_overlap(mesh->gdata.bbt, sb);
						//if (list == NULL) continue;
						while (list) {
							b = GTS_BBOX(list->data);
							point[edge] = SegmentTriangleIntersection(segments[edge], GTS_TRIANGLE(b->bounded));
							if (point[edge]) {
								break;
							}
							list = list->next;
						}

						if(point[edge]!=NULL){
							aux.push_back(elem4->nodes[6].id);
							aux.push_back(point[edge]->x);
							aux.push_back(point[edge]->y);
							aux.push_back(point[edge]->z);
						}
					}

					//face z+
					if(true){
						node1 = elem0->nodes[1].id;
						node2 = elem3->nodes[2].id;

						v1 = gts_vertex_new(gts_vertex_class(), coords[node1 * 3 + 0], coords[node1 * 3 + 1], coords[node1 * 3 + 2]);
						v2 = gts_vertex_new(gts_vertex_class(), coords[node2 * 3 + 0], coords[node2 * 3 + 1], coords[node2 * 3 + 2]);

						point[edge] = NULL;
						segments[edge] = gts_segment_new(gts_segment_class(), v1, v2);
						sb = gts_bbox_segment(gts_bbox_class(), segments[edge]);
						list = gts_bb_tree_overlap(mesh->gdata.bbt, sb);
						//if (list == NULL) continue;
						while (list) {
							b = GTS_BBOX(list->data);
							point[edge] = SegmentTriangleIntersection(segments[edge], GTS_TRIANGLE(b->bounded));
							if (point[edge]) {
								break;
							}
							list = list->next;
						}

						if(point[edge]!=NULL){
							aux.push_back(elem0->nodes[2].id);
							aux.push_back(point[edge]->x);
							aux.push_back(point[edge]->y);
							aux.push_back(point[edge]->z);
						}
					}

				}

				if(oct->id[1]!=-1 && oct->id[5]!=-1){
					//essa é para y...
					octant_t* elem0 = (octant_t*)sc_array_index(&mesh->elements,oct->id[0]);
					octant_t* elem1 = (octant_t*)sc_array_index(&mesh->elements,oct->id[1]);
					octant_t* elem4 = (octant_t*)sc_array_index(&mesh->elements,oct->id[4]);
					octant_t* elem5 = (octant_t*)sc_array_index(&mesh->elements,oct->id[5]);

					//face x-
					if(true){
						node1 = elem0->nodes[3].id;
						node2 = elem4->nodes[7].id;

						v1 = gts_vertex_new(gts_vertex_class(), coords[node1 * 3 + 0], coords[node1 * 3 + 1], coords[node1 * 3 + 2]);
						v2 = gts_vertex_new(gts_vertex_class(), coords[node2 * 3 + 0], coords[node2 * 3 + 1], coords[node2 * 3 + 2]);

						point[edge] = NULL;
						segments[edge] = gts_segment_new(gts_segment_class(), v1, v2);
						sb = gts_bbox_segment(gts_bbox_class(), segments[edge]);
						list = gts_bb_tree_overlap(mesh->gdata.bbt, sb);
						//if (list == NULL) continue;
						while (list) {
							b = GTS_BBOX(list->data);
							point[edge] = SegmentTriangleIntersection(segments[edge], GTS_TRIANGLE(b->bounded));
							if (point[edge]) {
								break;
							}
							list = list->next;
						}

						if(point[edge]!=NULL){
							aux.push_back(elem0->nodes[7].id);
							aux.push_back(point[edge]->x);
							aux.push_back(point[edge]->y);
							aux.push_back(point[edge]->z);
						}
					}

					//face x+
					if(true){
						node1 = elem1->nodes[2].id;
						node2 = elem5->nodes[6].id;

						v1 = gts_vertex_new(gts_vertex_class(), coords[node1 * 3 + 0], coords[node1 * 3 + 1], coords[node1 * 3 + 2]);
						v2 = gts_vertex_new(gts_vertex_class(), coords[node2 * 3 + 0], coords[node2 * 3 + 1], coords[node2 * 3 + 2]);

						point[edge] = NULL;
						segments[edge] = gts_segment_new(gts_segment_class(), v1, v2);
						sb = gts_bbox_segment(gts_bbox_class(), segments[edge]);
						list = gts_bb_tree_overlap(mesh->gdata.bbt, sb);
						//if (list == NULL) continue;
						while (list) {
							b = GTS_BBOX(list->data);
							point[edge] = SegmentTriangleIntersection(segments[edge], GTS_TRIANGLE(b->bounded));
							if (point[edge]) {
								break;
							}
							list = list->next;
						}

						if(point[edge]!=NULL){
							aux.push_back(elem1->nodes[6].id);
							aux.push_back(point[edge]->x);
							aux.push_back(point[edge]->y);
							aux.push_back(point[edge]->z);
						}
					}

					//face y-
					if(true){
						node1 = elem0->nodes[4].id;
						node2 = elem1->nodes[5].id;

						v1 = gts_vertex_new(gts_vertex_class(), coords[node1 * 3 + 0], coords[node1 * 3 + 1], coords[node1 * 3 + 2]);
						v2 = gts_vertex_new(gts_vertex_class(), coords[node2 * 3 + 0], coords[node2 * 3 + 1], coords[node2 * 3 + 2]);

						point[edge] = NULL;
						segments[edge] = gts_segment_new(gts_segment_class(), v1, v2);
						sb = gts_bbox_segment(gts_bbox_class(), segments[edge]);
						list = gts_bb_tree_overlap(mesh->gdata.bbt, sb);
						//if (list == NULL) continue;
						while (list) {
							b = GTS_BBOX(list->data);
							point[edge] = SegmentTriangleIntersection(segments[edge], GTS_TRIANGLE(b->bounded));
							if (point[edge]) {
								break;
							}
							list = list->next;
						}

						if(point[edge]!=NULL){
							aux.push_back(elem0->nodes[5].id);
							aux.push_back(point[edge]->x);
							aux.push_back(point[edge]->y);
							aux.push_back(point[edge]->z);
						}
					}

					//face y+
					/*
				if(false){

					node1 = elem3->nodes[2].id;
					node2 = elem7->nodes[6].id;

					v1 = gts_vertex_new(gts_vertex_class(), coords[node1 * 3 + 0], coords[node1 * 3 + 1], coords[node1 * 3 + 2]);
					v2 = gts_vertex_new(gts_vertex_class(), coords[node2 * 3 + 0], coords[node2 * 3 + 1], coords[node2 * 3 + 2]);

					point[edge] = NULL;
					segments[edge] = gts_segment_new(gts_segment_class(), v1, v2);
					sb = gts_bbox_segment(gts_bbox_class(), segments[edge]);
					list = gts_bb_tree_overlap(mesh->gdata.bbt, sb);
					if (list == NULL) continue;
					while (list) {
						b = GTS_BBOX(list->data);
						point[edge] = SegmentTriangleIntersection(segments[edge], GTS_TRIANGLE(b->bounded));
						if (point[edge]) {
							break;
						}
						list = list->next;
					}

					if(point[edge]!=NULL){
						aux.push_back(elem3->nodes[6].id);
						aux.push_back(point[edge]->x);
						aux.push_back(point[edge]->y);
						aux.push_back(point[edge]->z);
					}
				}
					 */

					//face z-
					if(true){
						node1 = elem4->nodes[7].id;
						node2 = elem5->nodes[6].id;

						v1 = gts_vertex_new(gts_vertex_class(), coords[node1 * 3 + 0], coords[node1 * 3 + 1], coords[node1 * 3 + 2]);
						v2 = gts_vertex_new(gts_vertex_class(), coords[node2 * 3 + 0], coords[node2 * 3 + 1], coords[node2 * 3 + 2]);

						point[edge] = NULL;
						segments[edge] = gts_segment_new(gts_segment_class(), v1, v2);
						sb = gts_bbox_segment(gts_bbox_class(), segments[edge]);
						list = gts_bb_tree_overlap(mesh->gdata.bbt, sb);
						//if (list == NULL) continue;
						while (list) {
							b = GTS_BBOX(list->data);
							point[edge] = SegmentTriangleIntersection(segments[edge], GTS_TRIANGLE(b->bounded));
							if (point[edge]) {
								break;
							}
							list = list->next;
						}

						if(point[edge]!=NULL){
							aux.push_back(elem4->nodes[6].id);
							aux.push_back(point[edge]->x);
							aux.push_back(point[edge]->y);
							aux.push_back(point[edge]->z);
						}
					}

					//face z+
					if(true){
						node1 = elem0->nodes[3].id;
						node2 = elem1->nodes[2].id;

						v1 = gts_vertex_new(gts_vertex_class(), coords[node1 * 3 + 0], coords[node1 * 3 + 1], coords[node1 * 3 + 2]);
						v2 = gts_vertex_new(gts_vertex_class(), coords[node2 * 3 + 0], coords[node2 * 3 + 1], coords[node2 * 3 + 2]);

						point[edge] = NULL;
						segments[edge] = gts_segment_new(gts_segment_class(), v1, v2);
						sb = gts_bbox_segment(gts_bbox_class(), segments[edge]);
						list = gts_bb_tree_overlap(mesh->gdata.bbt, sb);
						//if (list == NULL) continue;
						while (list) {
							b = GTS_BBOX(list->data);
							point[edge] = SegmentTriangleIntersection(segments[edge], GTS_TRIANGLE(b->bounded));
							if (point[edge]) {
								break;
							}
							list = list->next;
						}

						if(point[edge]!=NULL){
							aux.push_back(elem0->nodes[2].id);
							aux.push_back(point[edge]->x);
							aux.push_back(point[edge]->y);
							aux.push_back(point[edge]->z);
						}
					}

				}
			}

			//TODO
			if(oc_count==2){
				printf("Trabalha o vagabundo\n");
			}

		}
	}

	//moving the nodes in the edges...
	if(true){
		//achando os pontos de onde a superficie corta o octree nas 12 arestas
		for (int ioc = 0; ioc < mesh->oct.elem_count; ++ioc) {
			octree_t* oct = (octree_t*)sc_array_index(&mesh->oct,ioc);

			int oc_count=0;
			for(int i =0;i<8;i++){
				if(oct->id[i]!=-1) {
					oc_count++;
				}
			}

			if(oc_count==8){
				for(int iel = 0; iel<8; iel++){
					octant_t* elem = (octant_t*)sc_array_index(&mesh->elements,oct->id[iel]);
					//verifica se as arestas foram cortadas
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
						if (list == NULL) continue;
						while (list) {
							GtsBBox *b = GTS_BBOX(list->data);
							point[edge] = SegmentTriangleIntersection(segments[edge], GTS_TRIANGLE(b->bounded));
							if (point[edge]) {
								break;
							}
							list = list->next;
						}
					}

					//debug
					if(false){
						if(ioc == 0 && (iel==0 || iel ==4)){
							for(int ide = 0; ide <12; ide++){
								printf("Aresta %d cortada? %d\n",ide,point[ide]);
								printf("Aresta do octree %d\n", oct->edge[ide]);
							}
						}
					}

					if(true){
						//teoricamente move os pontos na arestas da face z-...
						if(oct->edge[8]){
							if(point[8]!=NULL){
								if(iel==4){
									aux.push_back(elem->nodes[5].id);
									aux.push_back(point[8]->x);
									aux.push_back(point[8]->y);
									aux.push_back(point[8]->z);

								}
								if(iel==5){
									aux.push_back(elem->nodes[4].id);
									aux.push_back(point[8]->x);
									aux.push_back(point[8]->y);
									aux.push_back(point[8]->z);

								}
							}
						}

						if(oct->edge[9]){
							if(point[9]!=NULL){
								if(iel==5){
									aux.push_back(elem->nodes[6].id);
									aux.push_back(point[9]->x);
									aux.push_back(point[9]->y);
									aux.push_back(point[9]->z);

								}
								if(iel==6){
									aux.push_back(elem->nodes[5].id);
									aux.push_back(point[9]->x);
									aux.push_back(point[9]->y);
									aux.push_back(point[9]->z);

								}
							}
						}

						if(oct->edge[10]){
							if(point[10]!=NULL){
								if(iel==6){
									aux.push_back(elem->nodes[7].id);
									aux.push_back(point[10]->x);
									aux.push_back(point[10]->y);
									aux.push_back(point[10]->z);

								}
								if(iel==7){
									aux.push_back(elem->nodes[6].id);
									aux.push_back(point[10]->x);
									aux.push_back(point[10]->y);
									aux.push_back(point[10]->z);

								}
							}
						}

						if(oct->edge[11]){
							if(point[11]!=NULL){
								if(iel==7){
									aux.push_back(elem->nodes[4].id);
									aux.push_back(point[11]->x);
									aux.push_back(point[11]->y);
									aux.push_back(point[11]->z);

								}
								if(iel==4){
									aux.push_back(elem->nodes[7].id);
									aux.push_back(point[11]->x);
									aux.push_back(point[11]->y);
									aux.push_back(point[11]->z);

								}
							}
						}

						//teoricamente move os pontos na arestas verticais...
						if(oct->edge[4]){
							if(point[4]!=NULL){
								if(iel==0){
									aux.push_back(elem->nodes[4].id);
									aux.push_back(point[4]->x);
									aux.push_back(point[4]->y);
									aux.push_back(point[4]->z);

								}
								if(iel==4){
									aux.push_back(elem->nodes[0].id);
									aux.push_back(point[4]->x);
									aux.push_back(point[4]->y);
									aux.push_back(point[4]->z);

								}
							}
						}

						if(oct->edge[5]){
							if(point[5]!=NULL){
								if(iel==1){
									aux.push_back(elem->nodes[5].id);
									aux.push_back(point[5]->x);
									aux.push_back(point[5]->y);
									aux.push_back(point[5]->z);

								}
								if(iel==5){
									aux.push_back(elem->nodes[1].id);
									aux.push_back(point[5]->x);
									aux.push_back(point[5]->y);
									aux.push_back(point[5]->z);

								}
							}
						}

						if(oct->edge[6]){
							if(point[6]!=NULL){
								if(iel==2){
									aux.push_back(elem->nodes[6].id);
									aux.push_back(point[6]->x);
									aux.push_back(point[6]->y);
									aux.push_back(point[6]->z);

								}
								if(iel==6){
									aux.push_back(elem->nodes[2].id);
									aux.push_back(point[6]->x);
									aux.push_back(point[6]->y);
									aux.push_back(point[6]->z);

								}
							}
						}

						if(oct->edge[7]){
							if(point[7]!=NULL){
								if(iel==3){
									aux.push_back(elem->nodes[7].id);
									aux.push_back(point[7]->x);
									aux.push_back(point[7]->y);
									aux.push_back(point[7]->z);

								}
								if(iel==7){
									aux.push_back(elem->nodes[3].id);
									aux.push_back(point[7]->x);
									aux.push_back(point[7]->y);
									aux.push_back(point[7]->z);

								}
							}
						}

						//teoricamente move os pontos na arestas da face z+...
						if(oct->edge[0]){
							if(point[0]!=NULL){
								if(iel==0){
									aux.push_back(elem->nodes[1].id);
									aux.push_back(point[0]->x);
									aux.push_back(point[0]->y);
									aux.push_back(point[0]->z);

								}
								if(iel==1){
									aux.push_back(elem->nodes[0].id);
									aux.push_back(point[0]->x);
									aux.push_back(point[0]->y);
									aux.push_back(point[0]->z);

								}
							}
						}

						if(oct->edge[1]){
							if(point[1]!=NULL){
								if(iel==1){
									aux.push_back(elem->nodes[2].id);
									aux.push_back(point[1]->x);
									aux.push_back(point[1]->y);
									aux.push_back(point[1]->z);

								}
								if(iel==2){
									aux.push_back(elem->nodes[1].id);
									aux.push_back(point[1]->x);
									aux.push_back(point[1]->y);
									aux.push_back(point[1]->z);

								}
							}
						}

						if(oct->edge[2]){
							if(point[2]!=NULL){
								if(iel==2){
									aux.push_back(elem->nodes[3].id);
									aux.push_back(point[2]->x);
									aux.push_back(point[2]->y);
									aux.push_back(point[2]->z);

								}
								if(iel==3){
									aux.push_back(elem->nodes[2].id);
									aux.push_back(point[2]->x);
									aux.push_back(point[2]->y);
									aux.push_back(point[2]->z);

								}
							}
						}

						if(oct->edge[3]){
							if(point[3]!=NULL){
								if(iel==3){
									aux.push_back(elem->nodes[0].id);
									aux.push_back(point[3]->x);
									aux.push_back(point[3]->y);
									aux.push_back(point[3]->z);

								}
								if(iel==0){
									aux.push_back(elem->nodes[3].id);
									aux.push_back(point[3]->x);
									aux.push_back(point[3]->y);
									aux.push_back(point[3]->z);

								}
							}
						}
					}
				}

			}

			if(oc_count==4){
				int node1, node2;
				GtsVertex *v1;
				GtsVertex *v2;
				int edge = 0;
				GtsBBox *sb;
				GSList* list;
				GtsBBox *b;

				if(oct->id[3]!=-1 && oct->id[7]!=-1){
					//esse e para o x
					octant_t* elem0 = (octant_t*)sc_array_index(&mesh->elements,oct->id[0]);
					octant_t* elem3 = (octant_t*)sc_array_index(&mesh->elements,oct->id[3]);
					octant_t* elem4 = (octant_t*)sc_array_index(&mesh->elements,oct->id[4]);
					octant_t* elem7 = (octant_t*)sc_array_index(&mesh->elements,oct->id[7]);

					//edge 0
					if(true){
						node1 = elem0->nodes[0].id;
						node2 = elem0->nodes[1].id;

						v1 = gts_vertex_new(gts_vertex_class(), coords[node1 * 3 + 0], coords[node1 * 3 + 1], coords[node1 * 3 + 2]);
						v2 = gts_vertex_new(gts_vertex_class(), coords[node2 * 3 + 0], coords[node2 * 3 + 1], coords[node2 * 3 + 2]);

						point[edge] = NULL;
						segments[edge] = gts_segment_new(gts_segment_class(), v1, v2);
						sb = gts_bbox_segment(gts_bbox_class(), segments[edge]);
						list = gts_bb_tree_overlap(mesh->gdata.bbt, sb);
						//if (list == NULL) continue;
						while (list) {
							b = GTS_BBOX(list->data);
							point[edge] = SegmentTriangleIntersection(segments[edge], GTS_TRIANGLE(b->bounded));
							if (point[edge]) {
								break;
							}
							list = list->next;
						}

						if(point[edge]!=NULL){
							aux.push_back(elem0->nodes[0].id);
							aux.push_back(point[edge]->x);
							aux.push_back(point[edge]->y);
							aux.push_back(point[edge]->z);
						}

					}

					//edge 1 nope
					if(false){
						node1 = elem0->nodes[0].id;
						node2 = elem0->nodes[1].id;

						v1 = gts_vertex_new(gts_vertex_class(), coords[node1 * 3 + 0], coords[node1 * 3 + 1], coords[node1 * 3 + 2]);
						v2 = gts_vertex_new(gts_vertex_class(), coords[node2 * 3 + 0], coords[node2 * 3 + 1], coords[node2 * 3 + 2]);

						point[edge] = NULL;
						segments[edge] = gts_segment_new(gts_segment_class(), v1, v2);
						sb = gts_bbox_segment(gts_bbox_class(), segments[edge]);
						list = gts_bb_tree_overlap(mesh->gdata.bbt, sb);
						//if (list == NULL) continue;
						while (list) {
							b = GTS_BBOX(list->data);
							point[edge] = SegmentTriangleIntersection(segments[edge], GTS_TRIANGLE(b->bounded));
							if (point[edge]) {
								break;
							}
							list = list->next;
						}

						if(point[edge]!=NULL){
							aux.push_back(elem0->nodes[0].id);
							aux.push_back(point[edge]->x);
							aux.push_back(point[edge]->y);
							aux.push_back(point[edge]->z);
						}

					}

					//edge 2
					if(true){
						node1 = elem3->nodes[2].id;
						node2 = elem3->nodes[3].id;

						v1 = gts_vertex_new(gts_vertex_class(), coords[node1 * 3 + 0], coords[node1 * 3 + 1], coords[node1 * 3 + 2]);
						v2 = gts_vertex_new(gts_vertex_class(), coords[node2 * 3 + 0], coords[node2 * 3 + 1], coords[node2 * 3 + 2]);

						point[edge] = NULL;
						segments[edge] = gts_segment_new(gts_segment_class(), v1, v2);
						sb = gts_bbox_segment(gts_bbox_class(), segments[edge]);
						list = gts_bb_tree_overlap(mesh->gdata.bbt, sb);
						//if (list == NULL) continue;
						while (list) {
							b = GTS_BBOX(list->data);
							point[edge] = SegmentTriangleIntersection(segments[edge], GTS_TRIANGLE(b->bounded));
							if (point[edge]) {
								break;
							}
							list = list->next;
						}

						if(point[edge]!=NULL){
							aux.push_back(elem3->nodes[3].id);
							aux.push_back(point[edge]->x);
							aux.push_back(point[edge]->y);
							aux.push_back(point[edge]->z);
						}

					}

					//edge 3
					if(true){
						node1 = elem0->nodes[0].id;
						node2 = elem3->nodes[3].id;

						v1 = gts_vertex_new(gts_vertex_class(), coords[node1 * 3 + 0], coords[node1 * 3 + 1], coords[node1 * 3 + 2]);
						v2 = gts_vertex_new(gts_vertex_class(), coords[node2 * 3 + 0], coords[node2 * 3 + 1], coords[node2 * 3 + 2]);

						point[edge] = NULL;
						segments[edge] = gts_segment_new(gts_segment_class(), v1, v2);
						sb = gts_bbox_segment(gts_bbox_class(), segments[edge]);
						list = gts_bb_tree_overlap(mesh->gdata.bbt, sb);
						//if (list == NULL) continue;
						while (list) {
							b = GTS_BBOX(list->data);
							point[edge] = SegmentTriangleIntersection(segments[edge], GTS_TRIANGLE(b->bounded));
							if (point[edge]) {
								break;
							}
							list = list->next;
						}

						if(point[edge]!=NULL){
							aux.push_back(elem0->nodes[3].id);
							aux.push_back(point[edge]->x);
							aux.push_back(point[edge]->y);
							aux.push_back(point[edge]->z);
						}

					}

					//edge 4
					if(true){
						node1 = elem0->nodes[0].id;
						node2 = elem4->nodes[4].id;

						v1 = gts_vertex_new(gts_vertex_class(), coords[node1 * 3 + 0], coords[node1 * 3 + 1], coords[node1 * 3 + 2]);
						v2 = gts_vertex_new(gts_vertex_class(), coords[node2 * 3 + 0], coords[node2 * 3 + 1], coords[node2 * 3 + 2]);

						point[edge] = NULL;
						segments[edge] = gts_segment_new(gts_segment_class(), v1, v2);
						sb = gts_bbox_segment(gts_bbox_class(), segments[edge]);
						list = gts_bb_tree_overlap(mesh->gdata.bbt, sb);
						//if (list == NULL) continue;
						while (list) {
							b = GTS_BBOX(list->data);
							point[edge] = SegmentTriangleIntersection(segments[edge], GTS_TRIANGLE(b->bounded));
							if (point[edge]) {
								break;
							}
							list = list->next;
						}

						if(point[edge]!=NULL){
							aux.push_back(elem0->nodes[4].id);
							aux.push_back(point[edge]->x);
							aux.push_back(point[edge]->y);
							aux.push_back(point[edge]->z);
						}

					}

					//edge 5 nope
					if(false){
						node1 = elem0->nodes[0].id;
						node2 = elem0->nodes[1].id;

						v1 = gts_vertex_new(gts_vertex_class(), coords[node1 * 3 + 0], coords[node1 * 3 + 1], coords[node1 * 3 + 2]);
						v2 = gts_vertex_new(gts_vertex_class(), coords[node2 * 3 + 0], coords[node2 * 3 + 1], coords[node2 * 3 + 2]);

						point[edge] = NULL;
						segments[edge] = gts_segment_new(gts_segment_class(), v1, v2);
						sb = gts_bbox_segment(gts_bbox_class(), segments[edge]);
						list = gts_bb_tree_overlap(mesh->gdata.bbt, sb);
						//if (list == NULL) continue;
						while (list) {
							b = GTS_BBOX(list->data);
							point[edge] = SegmentTriangleIntersection(segments[edge], GTS_TRIANGLE(b->bounded));
							if (point[edge]) {
								break;
							}
							list = list->next;
						}

						if(point[edge]!=NULL){
							aux.push_back(elem0->nodes[0].id);
							aux.push_back(point[edge]->x);
							aux.push_back(point[edge]->y);
							aux.push_back(point[edge]->z);
						}

					}

					//edge 6 nope
					if(false){
						node1 = elem7->nodes[6].id;
						node2 = elem7->nodes[7].id;

						v1 = gts_vertex_new(gts_vertex_class(), coords[node1 * 3 + 0], coords[node1 * 3 + 1], coords[node1 * 3 + 2]);
						v2 = gts_vertex_new(gts_vertex_class(), coords[node2 * 3 + 0], coords[node2 * 3 + 1], coords[node2 * 3 + 2]);

						point[edge] = NULL;
						segments[edge] = gts_segment_new(gts_segment_class(), v1, v2);
						sb = gts_bbox_segment(gts_bbox_class(), segments[edge]);
						list = gts_bb_tree_overlap(mesh->gdata.bbt, sb);
						//if (list == NULL) continue;
						while (list) {
							b = GTS_BBOX(list->data);
							point[edge] = SegmentTriangleIntersection(segments[edge], GTS_TRIANGLE(b->bounded));
							if (point[edge]) {
								break;
							}
							list = list->next;
						}

						if(point[edge]!=NULL){
							aux.push_back(elem7->nodes[7].id);
							aux.push_back(point[edge]->x);
							aux.push_back(point[edge]->y);
							aux.push_back(point[edge]->z);
						}

					}

					//edge 7
					if(true){
						node1 = elem4->nodes[3].id;
						node2 = elem7->nodes[7].id;

						v1 = gts_vertex_new(gts_vertex_class(), coords[node1 * 3 + 0], coords[node1 * 3 + 1], coords[node1 * 3 + 2]);
						v2 = gts_vertex_new(gts_vertex_class(), coords[node2 * 3 + 0], coords[node2 * 3 + 1], coords[node2 * 3 + 2]);

						point[edge] = NULL;
						segments[edge] = gts_segment_new(gts_segment_class(), v1, v2);
						sb = gts_bbox_segment(gts_bbox_class(), segments[edge]);
						list = gts_bb_tree_overlap(mesh->gdata.bbt, sb);
						//if (list == NULL) continue;
						while (list) {
							b = GTS_BBOX(list->data);
							point[edge] = SegmentTriangleIntersection(segments[edge], GTS_TRIANGLE(b->bounded));
							if (point[edge]) {
								break;
							}
							list = list->next;
						}

						if(point[edge]!=NULL){
							aux.push_back(elem4->nodes[7].id);
							aux.push_back(point[edge]->x);
							aux.push_back(point[edge]->y);
							aux.push_back(point[edge]->z);
						}

					}

					//edge 8
					if(true){
						node1 = elem4->nodes[4].id;
						node2 = elem4->nodes[5].id;

						v1 = gts_vertex_new(gts_vertex_class(), coords[node1 * 3 + 0], coords[node1 * 3 + 1], coords[node1 * 3 + 2]);
						v2 = gts_vertex_new(gts_vertex_class(), coords[node2 * 3 + 0], coords[node2 * 3 + 1], coords[node2 * 3 + 2]);

						point[edge] = NULL;
						segments[edge] = gts_segment_new(gts_segment_class(), v1, v2);
						sb = gts_bbox_segment(gts_bbox_class(), segments[edge]);
						list = gts_bb_tree_overlap(mesh->gdata.bbt, sb);
						//if (list == NULL) continue;
						while (list) {
							b = GTS_BBOX(list->data);
							point[edge] = SegmentTriangleIntersection(segments[edge], GTS_TRIANGLE(b->bounded));
							if (point[edge]) {
								break;
							}
							list = list->next;
						}

						if(point[edge]!=NULL){
							aux.push_back(elem4->nodes[4].id);
							aux.push_back(point[edge]->x);
							aux.push_back(point[edge]->y);
							aux.push_back(point[edge]->z);
						}

					}

					//edge 9 nope
					if(false){
						node1 = elem0->nodes[0].id;
						node2 = elem0->nodes[1].id;

						v1 = gts_vertex_new(gts_vertex_class(), coords[node1 * 3 + 0], coords[node1 * 3 + 1], coords[node1 * 3 + 2]);
						v2 = gts_vertex_new(gts_vertex_class(), coords[node2 * 3 + 0], coords[node2 * 3 + 1], coords[node2 * 3 + 2]);

						point[edge] = NULL;
						segments[edge] = gts_segment_new(gts_segment_class(), v1, v2);
						sb = gts_bbox_segment(gts_bbox_class(), segments[edge]);
						list = gts_bb_tree_overlap(mesh->gdata.bbt, sb);
						//if (list == NULL) continue;
						while (list) {
							b = GTS_BBOX(list->data);
							point[edge] = SegmentTriangleIntersection(segments[edge], GTS_TRIANGLE(b->bounded));
							if (point[edge]) {
								break;
							}
							list = list->next;
						}

						if(point[edge]!=NULL){
							aux.push_back(elem0->nodes[0].id);
							aux.push_back(point[edge]->x);
							aux.push_back(point[edge]->y);
							aux.push_back(point[edge]->z);
						}

					}

					//edge 10
					if(true){
						node1 = elem7->nodes[6].id;
						node2 = elem7->nodes[7].id;

						v1 = gts_vertex_new(gts_vertex_class(), coords[node1 * 3 + 0], coords[node1 * 3 + 1], coords[node1 * 3 + 2]);
						v2 = gts_vertex_new(gts_vertex_class(), coords[node2 * 3 + 0], coords[node2 * 3 + 1], coords[node2 * 3 + 2]);

						point[edge] = NULL;
						segments[edge] = gts_segment_new(gts_segment_class(), v1, v2);
						sb = gts_bbox_segment(gts_bbox_class(), segments[edge]);
						list = gts_bb_tree_overlap(mesh->gdata.bbt, sb);
						//if (list == NULL) continue;
						while (list) {
							b = GTS_BBOX(list->data);
							point[edge] = SegmentTriangleIntersection(segments[edge], GTS_TRIANGLE(b->bounded));
							if (point[edge]) {
								break;
							}
							list = list->next;
						}

						if(point[edge]!=NULL){
							aux.push_back(elem7->nodes[7].id);
							aux.push_back(point[edge]->x);
							aux.push_back(point[edge]->y);
							aux.push_back(point[edge]->z);
						}

					}

					//edge 11
					if(true){
						node1 = elem4->nodes[4].id;
						node2 = elem7->nodes[7].id;

						v1 = gts_vertex_new(gts_vertex_class(), coords[node1 * 3 + 0], coords[node1 * 3 + 1], coords[node1 * 3 + 2]);
						v2 = gts_vertex_new(gts_vertex_class(), coords[node2 * 3 + 0], coords[node2 * 3 + 1], coords[node2 * 3 + 2]);

						point[edge] = NULL;
						segments[edge] = gts_segment_new(gts_segment_class(), v1, v2);
						sb = gts_bbox_segment(gts_bbox_class(), segments[edge]);
						list = gts_bb_tree_overlap(mesh->gdata.bbt, sb);
						//if (list == NULL) continue;
						while (list) {
							b = GTS_BBOX(list->data);
							point[edge] = SegmentTriangleIntersection(segments[edge], GTS_TRIANGLE(b->bounded));
							if (point[edge]) {
								break;
							}
							list = list->next;
						}

						if(point[edge]!=NULL){
							aux.push_back(elem4->nodes[7].id);
							aux.push_back(point[edge]->x);
							aux.push_back(point[edge]->y);
							aux.push_back(point[edge]->z);
						}

					}

				}

				if(oct->id[1]!=-1 && oct->id[5]!=-1){
					//esse e para o y
					octant_t* elem0 = (octant_t*)sc_array_index(&mesh->elements,oct->id[0]);
					octant_t* elem1 = (octant_t*)sc_array_index(&mesh->elements,oct->id[1]);
					octant_t* elem4 = (octant_t*)sc_array_index(&mesh->elements,oct->id[4]);
					octant_t* elem5 = (octant_t*)sc_array_index(&mesh->elements,oct->id[5]);

					//edge 0
					if(true){
						node1 = elem0->nodes[0].id;
						node2 = elem1->nodes[1].id;

						v1 = gts_vertex_new(gts_vertex_class(), coords[node1 * 3 + 0], coords[node1 * 3 + 1], coords[node1 * 3 + 2]);
						v2 = gts_vertex_new(gts_vertex_class(), coords[node2 * 3 + 0], coords[node2 * 3 + 1], coords[node2 * 3 + 2]);

						point[edge] = NULL;
						segments[edge] = gts_segment_new(gts_segment_class(), v1, v2);
						sb = gts_bbox_segment(gts_bbox_class(), segments[edge]);
						list = gts_bb_tree_overlap(mesh->gdata.bbt, sb);
						//if (list == NULL) continue;
						while (list) {
							b = GTS_BBOX(list->data);
							point[edge] = SegmentTriangleIntersection(segments[edge], GTS_TRIANGLE(b->bounded));
							if (point[edge]) {
								break;
							}
							list = list->next;
						}

						if(point[edge]!=NULL){
							aux.push_back(elem0->nodes[1].id);
							aux.push_back(point[edge]->x);
							aux.push_back(point[edge]->y);
							aux.push_back(point[edge]->z);
						}

					}

					//edge 1
					if(true){
						node1 = elem1->nodes[1].id;
						node2 = elem1->nodes[2].id;

						v1 = gts_vertex_new(gts_vertex_class(), coords[node1 * 3 + 0], coords[node1 * 3 + 1], coords[node1 * 3 + 2]);
						v2 = gts_vertex_new(gts_vertex_class(), coords[node2 * 3 + 0], coords[node2 * 3 + 1], coords[node2 * 3 + 2]);

						point[edge] = NULL;
						segments[edge] = gts_segment_new(gts_segment_class(), v1, v2);
						sb = gts_bbox_segment(gts_bbox_class(), segments[edge]);
						list = gts_bb_tree_overlap(mesh->gdata.bbt, sb);
						//if (list == NULL) continue;
						while (list) {
							b = GTS_BBOX(list->data);
							point[edge] = SegmentTriangleIntersection(segments[edge], GTS_TRIANGLE(b->bounded));
							if (point[edge]) {
								break;
							}
							list = list->next;
						}

						if(point[edge]!=NULL){
							aux.push_back(elem1->nodes[1].id);
							aux.push_back(point[edge]->x);
							aux.push_back(point[edge]->y);
							aux.push_back(point[edge]->z);
						}

					}

					//edge 2 nope
					if(false){
						node1 = elem0->nodes[2].id;
						node2 = elem0->nodes[3].id;

						v1 = gts_vertex_new(gts_vertex_class(), coords[node1 * 3 + 0], coords[node1 * 3 + 1], coords[node1 * 3 + 2]);
						v2 = gts_vertex_new(gts_vertex_class(), coords[node2 * 3 + 0], coords[node2 * 3 + 1], coords[node2 * 3 + 2]);

						point[edge] = NULL;
						segments[edge] = gts_segment_new(gts_segment_class(), v1, v2);
						sb = gts_bbox_segment(gts_bbox_class(), segments[edge]);
						list = gts_bb_tree_overlap(mesh->gdata.bbt, sb);
						//if (list == NULL) continue;
						while (list) {
							b = GTS_BBOX(list->data);
							point[edge] = SegmentTriangleIntersection(segments[edge], GTS_TRIANGLE(b->bounded));
							if (point[edge]) {
								break;
							}
							list = list->next;
						}

						if(point[edge]!=NULL){
							aux.push_back(elem0->nodes[3].id);
							aux.push_back(point[edge]->x);
							aux.push_back(point[edge]->y);
							aux.push_back(point[edge]->z);
						}

					}

					//edge 3
					if(true){
						node1 = elem0->nodes[0].id;
						node2 = elem0->nodes[3].id;

						v1 = gts_vertex_new(gts_vertex_class(), coords[node1 * 3 + 0], coords[node1 * 3 + 1], coords[node1 * 3 + 2]);
						v2 = gts_vertex_new(gts_vertex_class(), coords[node2 * 3 + 0], coords[node2 * 3 + 1], coords[node2 * 3 + 2]);

						point[edge] = NULL;
						segments[edge] = gts_segment_new(gts_segment_class(), v1, v2);
						sb = gts_bbox_segment(gts_bbox_class(), segments[edge]);
						list = gts_bb_tree_overlap(mesh->gdata.bbt, sb);
						//if (list == NULL) continue;
						while (list) {
							b = GTS_BBOX(list->data);
							point[edge] = SegmentTriangleIntersection(segments[edge], GTS_TRIANGLE(b->bounded));
							if (point[edge]) {
								break;
							}
							list = list->next;
						}

						if(point[edge]!=NULL){
							aux.push_back(elem0->nodes[0].id);
							aux.push_back(point[edge]->x);
							aux.push_back(point[edge]->y);
							aux.push_back(point[edge]->z);
						}

					}

					//edge 4
					if(true){
						node1 = elem0->nodes[0].id;
						node2 = elem4->nodes[4].id;

						v1 = gts_vertex_new(gts_vertex_class(), coords[node1 * 3 + 0], coords[node1 * 3 + 1], coords[node1 * 3 + 2]);
						v2 = gts_vertex_new(gts_vertex_class(), coords[node2 * 3 + 0], coords[node2 * 3 + 1], coords[node2 * 3 + 2]);

						point[edge] = NULL;
						segments[edge] = gts_segment_new(gts_segment_class(), v1, v2);
						sb = gts_bbox_segment(gts_bbox_class(), segments[edge]);
						list = gts_bb_tree_overlap(mesh->gdata.bbt, sb);
						//if (list == NULL) continue;
						while (list) {
							b = GTS_BBOX(list->data);
							point[edge] = SegmentTriangleIntersection(segments[edge], GTS_TRIANGLE(b->bounded));
							if (point[edge]) {
								break;
							}
							list = list->next;
						}

						if(point[edge]!=NULL){
							aux.push_back(elem0->nodes[4].id);
							aux.push_back(point[edge]->x);
							aux.push_back(point[edge]->y);
							aux.push_back(point[edge]->z);
						}

					}

					//edge 5
					if(true){
						node1 = elem1->nodes[1].id;
						node2 = elem5->nodes[5].id;

						v1 = gts_vertex_new(gts_vertex_class(), coords[node1 * 3 + 0], coords[node1 * 3 + 1], coords[node1 * 3 + 2]);
						v2 = gts_vertex_new(gts_vertex_class(), coords[node2 * 3 + 0], coords[node2 * 3 + 1], coords[node2 * 3 + 2]);

						point[edge] = NULL;
						segments[edge] = gts_segment_new(gts_segment_class(), v1, v2);
						sb = gts_bbox_segment(gts_bbox_class(), segments[edge]);
						list = gts_bb_tree_overlap(mesh->gdata.bbt, sb);
						//if (list == NULL) continue;
						while (list) {
							b = GTS_BBOX(list->data);
							point[edge] = SegmentTriangleIntersection(segments[edge], GTS_TRIANGLE(b->bounded));
							if (point[edge]) {
								break;
							}
							list = list->next;
						}

						if(point[edge]!=NULL){
							aux.push_back(elem1->nodes[5].id);
							aux.push_back(point[edge]->x);
							aux.push_back(point[edge]->y);
							aux.push_back(point[edge]->z);
						}

					}

					//edge 6 nope
					if(false){
						node1 = elem0->nodes[6].id;
						node2 = elem0->nodes[7].id;

						v1 = gts_vertex_new(gts_vertex_class(), coords[node1 * 3 + 0], coords[node1 * 3 + 1], coords[node1 * 3 + 2]);
						v2 = gts_vertex_new(gts_vertex_class(), coords[node2 * 3 + 0], coords[node2 * 3 + 1], coords[node2 * 3 + 2]);

						point[edge] = NULL;
						segments[edge] = gts_segment_new(gts_segment_class(), v1, v2);
						sb = gts_bbox_segment(gts_bbox_class(), segments[edge]);
						list = gts_bb_tree_overlap(mesh->gdata.bbt, sb);
						//if (list == NULL) continue;
						while (list) {
							b = GTS_BBOX(list->data);
							point[edge] = SegmentTriangleIntersection(segments[edge], GTS_TRIANGLE(b->bounded));
							if (point[edge]) {
								break;
							}
							list = list->next;
						}

						if(point[edge]!=NULL){
							aux.push_back(elem0->nodes[7].id);
							aux.push_back(point[edge]->x);
							aux.push_back(point[edge]->y);
							aux.push_back(point[edge]->z);
						}

					}

					//edge 7 nope
					if(false){
						node1 = elem0->nodes[3].id;
						node2 = elem0->nodes[7].id;

						v1 = gts_vertex_new(gts_vertex_class(), coords[node1 * 3 + 0], coords[node1 * 3 + 1], coords[node1 * 3 + 2]);
						v2 = gts_vertex_new(gts_vertex_class(), coords[node2 * 3 + 0], coords[node2 * 3 + 1], coords[node2 * 3 + 2]);

						point[edge] = NULL;
						segments[edge] = gts_segment_new(gts_segment_class(), v1, v2);
						sb = gts_bbox_segment(gts_bbox_class(), segments[edge]);
						list = gts_bb_tree_overlap(mesh->gdata.bbt, sb);
						//if (list == NULL) continue;
						while (list) {
							b = GTS_BBOX(list->data);
							point[edge] = SegmentTriangleIntersection(segments[edge], GTS_TRIANGLE(b->bounded));
							if (point[edge]) {
								break;
							}
							list = list->next;
						}

						if(point[edge]!=NULL){
							aux.push_back(elem0->nodes[7].id);
							aux.push_back(point[edge]->x);
							aux.push_back(point[edge]->y);
							aux.push_back(point[edge]->z);
						}

					}

					//edge 8
					if(true){
						node1 = elem4->nodes[4].id;
						node2 = elem5->nodes[5].id;

						v1 = gts_vertex_new(gts_vertex_class(), coords[node1 * 3 + 0], coords[node1 * 3 + 1], coords[node1 * 3 + 2]);
						v2 = gts_vertex_new(gts_vertex_class(), coords[node2 * 3 + 0], coords[node2 * 3 + 1], coords[node2 * 3 + 2]);

						point[edge] = NULL;
						segments[edge] = gts_segment_new(gts_segment_class(), v1, v2);
						sb = gts_bbox_segment(gts_bbox_class(), segments[edge]);
						list = gts_bb_tree_overlap(mesh->gdata.bbt, sb);
						//if (list == NULL) continue;
						while (list) {
							b = GTS_BBOX(list->data);
							point[edge] = SegmentTriangleIntersection(segments[edge], GTS_TRIANGLE(b->bounded));
							if (point[edge]) {
								break;
							}
							list = list->next;
						}

						if(point[edge]!=NULL){
							aux.push_back(elem4->nodes[5].id);
							aux.push_back(point[edge]->x);
							aux.push_back(point[edge]->y);
							aux.push_back(point[edge]->z);
						}

					}

					//edge 9
					if(true){
						node1 = elem5->nodes[4].id;
						node2 = elem5->nodes[6].id;

						v1 = gts_vertex_new(gts_vertex_class(), coords[node1 * 3 + 0], coords[node1 * 3 + 1], coords[node1 * 3 + 2]);
						v2 = gts_vertex_new(gts_vertex_class(), coords[node2 * 3 + 0], coords[node2 * 3 + 1], coords[node2 * 3 + 2]);

						point[edge] = NULL;
						segments[edge] = gts_segment_new(gts_segment_class(), v1, v2);
						sb = gts_bbox_segment(gts_bbox_class(), segments[edge]);
						list = gts_bb_tree_overlap(mesh->gdata.bbt, sb);
						//if (list == NULL) continue;
						while (list) {
							b = GTS_BBOX(list->data);
							point[edge] = SegmentTriangleIntersection(segments[edge], GTS_TRIANGLE(b->bounded));
							if (point[edge]) {
								break;
							}
							list = list->next;
						}

						if(point[edge]!=NULL){
							aux.push_back(elem5->nodes[5].id);
							aux.push_back(point[edge]->x);
							aux.push_back(point[edge]->y);
							aux.push_back(point[edge]->z);
						}

					}

					//edge 10 nope
					if(false){
						node1 = elem0->nodes[2].id;
						node2 = elem0->nodes[3].id;

						v1 = gts_vertex_new(gts_vertex_class(), coords[node1 * 3 + 0], coords[node1 * 3 + 1], coords[node1 * 3 + 2]);
						v2 = gts_vertex_new(gts_vertex_class(), coords[node2 * 3 + 0], coords[node2 * 3 + 1], coords[node2 * 3 + 2]);

						point[edge] = NULL;
						segments[edge] = gts_segment_new(gts_segment_class(), v1, v2);
						sb = gts_bbox_segment(gts_bbox_class(), segments[edge]);
						list = gts_bb_tree_overlap(mesh->gdata.bbt, sb);
						//if (list == NULL) continue;
						while (list) {
							b = GTS_BBOX(list->data);
							point[edge] = SegmentTriangleIntersection(segments[edge], GTS_TRIANGLE(b->bounded));
							if (point[edge]) {
								break;
							}
							list = list->next;
						}

						if(point[edge]!=NULL){
							aux.push_back(elem0->nodes[3].id);
							aux.push_back(point[edge]->x);
							aux.push_back(point[edge]->y);
							aux.push_back(point[edge]->z);
						}

					}

					//edge 11
					if(true){
						node1 = elem4->nodes[4].id;
						node2 = elem4->nodes[7].id;

						v1 = gts_vertex_new(gts_vertex_class(), coords[node1 * 3 + 0], coords[node1 * 3 + 1], coords[node1 * 3 + 2]);
						v2 = gts_vertex_new(gts_vertex_class(), coords[node2 * 3 + 0], coords[node2 * 3 + 1], coords[node2 * 3 + 2]);

						point[edge] = NULL;
						segments[edge] = gts_segment_new(gts_segment_class(), v1, v2);
						sb = gts_bbox_segment(gts_bbox_class(), segments[edge]);
						list = gts_bb_tree_overlap(mesh->gdata.bbt, sb);
						//if (list == NULL) continue;
						while (list) {
							b = GTS_BBOX(list->data);
							point[edge] = SegmentTriangleIntersection(segments[edge], GTS_TRIANGLE(b->bounded));
							if (point[edge]) {
								break;
							}
							list = list->next;
						}

						if(point[edge]!=NULL){
							aux.push_back(elem4->nodes[4].id);
							aux.push_back(point[edge]->x);
							aux.push_back(point[edge]->y);
							aux.push_back(point[edge]->z);
						}

					}

				}

			}

			//TODO
			if(oc_count==2){
				printf("Trabalha o vagabundo\n");
			}
		}
	}

	//moving the central node of the octree...
	if(true){
		//achando os pontos de onde a superficie corta o octree no central
		for (int ioc = 0; ioc < mesh->oct.elem_count; ++ioc) {
			octree_t* oct = (octree_t*)sc_array_index(&mesh->oct,ioc);

			int oc_count=0;
			for(int i =0;i<8;i++){
				if(oct->id[i]!=-1) {
					oc_count++;
				}
			}

			if(oc_count==8){
				octant_t* elem0 = (octant_t*)sc_array_index(&mesh->elements,oct->id[0]);
				octant_t* elem1 = (octant_t*)sc_array_index(&mesh->elements,oct->id[1]);
				octant_t* elem2 = (octant_t*)sc_array_index(&mesh->elements,oct->id[2]);
				octant_t* elem3 = (octant_t*)sc_array_index(&mesh->elements,oct->id[3]);

				octant_t* elem4 = (octant_t*)sc_array_index(&mesh->elements,oct->id[4]);
				octant_t* elem5 = (octant_t*)sc_array_index(&mesh->elements,oct->id[5]);
				octant_t* elem6 = (octant_t*)sc_array_index(&mesh->elements,oct->id[6]);
				octant_t* elem7 = (octant_t*)sc_array_index(&mesh->elements,oct->id[7]);

				int node1, node2;
				GtsVertex *v1;
				GtsVertex *v2;
				int edge = 0;
				GtsBBox *sb;
				GSList* list;
				GtsBBox *b;

				if(true){
					node1 = elem0->nodes[2].id;
					node2 = elem4->nodes[6].id;

					v1 = gts_vertex_new(gts_vertex_class(), coords[node1 * 3 + 0], coords[node1 * 3 + 1], coords[node1 * 3 + 2]);
					v2 = gts_vertex_new(gts_vertex_class(), coords[node2 * 3 + 0], coords[node2 * 3 + 1], coords[node2 * 3 + 2]);

					point[edge] = NULL;
					segments[edge] = gts_segment_new(gts_segment_class(), v1, v2);
					sb = gts_bbox_segment(gts_bbox_class(), segments[edge]);
					list = gts_bb_tree_overlap(mesh->gdata.bbt, sb);
					//if (list == NULL) continue;
					while (list) {
						b = GTS_BBOX(list->data);
						point[edge] = SegmentTriangleIntersection(segments[edge], GTS_TRIANGLE(b->bounded));
						if (point[edge]) {
							break;
						}
						list = list->next;
					}

					if(point[edge]!=NULL){
						aux.push_back(elem0->nodes[6].id);
						aux.push_back(point[edge]->x);
						aux.push_back(point[edge]->y);
						aux.push_back(point[edge]->z);
					}else{

						edge = 0; // 0 e 6
						if(true){
							node1 = elem0->nodes[0].id;
							node2 = elem6->nodes[6].id;

							v1 = gts_vertex_new(gts_vertex_class(), coords[node1 * 3 + 0], coords[node1 * 3 + 1], coords[node1 * 3 + 2]);
							v2 = gts_vertex_new(gts_vertex_class(), coords[node2 * 3 + 0], coords[node2 * 3 + 1], coords[node2 * 3 + 2]);

							point[edge] = NULL;
							segments[edge] = gts_segment_new(gts_segment_class(), v1, v2);
							sb = gts_bbox_segment(gts_bbox_class(), segments[edge]);
							list = gts_bb_tree_overlap(mesh->gdata.bbt, sb);
							//if (list == NULL) continue;
							while (list) {
								b = GTS_BBOX(list->data);
								point[edge] = SegmentTriangleIntersection(segments[edge], GTS_TRIANGLE(b->bounded));
								if (point[edge]) {
									break;
								}
								list = list->next;
							}
						}

						edge = 1; // 1 e 7
						if(true){
							node1 = elem1->nodes[1].id;
							node2 = elem7->nodes[7].id;

							v1 = gts_vertex_new(gts_vertex_class(), coords[node1 * 3 + 0], coords[node1 * 3 + 1], coords[node1 * 3 + 2]);
							v2 = gts_vertex_new(gts_vertex_class(), coords[node2 * 3 + 0], coords[node2 * 3 + 1], coords[node2 * 3 + 2]);

							point[edge] = NULL;
							segments[edge] = gts_segment_new(gts_segment_class(), v1, v2);
							sb = gts_bbox_segment(gts_bbox_class(), segments[edge]);
							list = gts_bb_tree_overlap(mesh->gdata.bbt, sb);
							//if (list == NULL) continue;
							while (list) {
								b = GTS_BBOX(list->data);
								point[edge] = SegmentTriangleIntersection(segments[edge], GTS_TRIANGLE(b->bounded));
								if (point[edge]) {
									break;
								}
								list = list->next;
							}
						}

						edge = 2; // 2 e 4
						if(true){
							node1 = elem2->nodes[2].id;
							node2 = elem4->nodes[4].id;

							v1 = gts_vertex_new(gts_vertex_class(), coords[node1 * 3 + 0], coords[node1 * 3 + 1], coords[node1 * 3 + 2]);
							v2 = gts_vertex_new(gts_vertex_class(), coords[node2 * 3 + 0], coords[node2 * 3 + 1], coords[node2 * 3 + 2]);

							point[edge] = NULL;
							segments[edge] = gts_segment_new(gts_segment_class(), v1, v2);
							sb = gts_bbox_segment(gts_bbox_class(), segments[edge]);
							list = gts_bb_tree_overlap(mesh->gdata.bbt, sb);
							//if (list == NULL) continue;
							while (list) {
								b = GTS_BBOX(list->data);
								point[edge] = SegmentTriangleIntersection(segments[edge], GTS_TRIANGLE(b->bounded));
								if (point[edge]) {
									break;
								}
								list = list->next;
							}
						}

						edge = 3; // 3 e 5
						if(true){
							node1 = elem3->nodes[3].id;
							node2 = elem5->nodes[5].id;

							v1 = gts_vertex_new(gts_vertex_class(), coords[node1 * 3 + 0], coords[node1 * 3 + 1], coords[node1 * 3 + 2]);
							v2 = gts_vertex_new(gts_vertex_class(), coords[node2 * 3 + 0], coords[node2 * 3 + 1], coords[node2 * 3 + 2]);

							point[edge] = NULL;
							segments[edge] = gts_segment_new(gts_segment_class(), v1, v2);
							sb = gts_bbox_segment(gts_bbox_class(), segments[edge]);
							list = gts_bb_tree_overlap(mesh->gdata.bbt, sb);
							//if (list == NULL) continue;
							while (list) {
								b = GTS_BBOX(list->data);
								point[edge] = SegmentTriangleIntersection(segments[edge], GTS_TRIANGLE(b->bounded));
								if (point[edge]) {
									break;
								}
								list = list->next;
							}
						}

						int np[4];
						int count_p = 0;
						for(int ip = 0; ip<4;ip++){
							if(point[ip]!=NULL){
								np[count_p] = ip;
								count_p++;
							}
						}

						if(count_p!=0){
							double xx, yy,zz;
							xx = 0;
							yy = 0;
							zz = 0;
							for(int ip = 0; ip<count_p;ip++){
								xx = point[np[ip]]->x+xx;
								yy = point[np[ip]]->y+yy;
								zz = point[np[ip]]->z+zz;
							}
							xx = xx/count_p;
							yy = yy/count_p;
							zz = zz/count_p;

							aux.push_back(elem0->nodes[6].id);
							aux.push_back(xx);
							aux.push_back(yy);
							aux.push_back(zz);

						}else{
							printf("Tu é mais burro que um tijolo 6 furos\n");
						}
					}
				}

			}

			if(oc_count==4){

				int node1, node2;
				GtsVertex *v1;
				GtsVertex *v2;
				int edge = 0;
				GtsBBox *sb;
				GSList* list;
				GtsBBox *b;

				if(oct->id[3]!=-1 && oct->id[7]!=-1){
					//essa é para x...
					octant_t* elem0 = (octant_t*)sc_array_index(&mesh->elements,oct->id[0]);
					octant_t* elem3 = (octant_t*)sc_array_index(&mesh->elements,oct->id[3]);
					octant_t* elem4 = (octant_t*)sc_array_index(&mesh->elements,oct->id[4]);
					octant_t* elem7 = (octant_t*)sc_array_index(&mesh->elements,oct->id[7]);

					edge = 0; // 0 e 7
					if(true){
						node1 = elem0->nodes[1].id;
						node2 = elem7->nodes[6].id;

						v1 = gts_vertex_new(gts_vertex_class(), coords[node1 * 3 + 0], coords[node1 * 3 + 1], coords[node1 * 3 + 2]);
						v2 = gts_vertex_new(gts_vertex_class(), coords[node2 * 3 + 0], coords[node2 * 3 + 1], coords[node2 * 3 + 2]);

						point[edge] = NULL;
						segments[edge] = gts_segment_new(gts_segment_class(), v1, v2);
						sb = gts_bbox_segment(gts_bbox_class(), segments[edge]);
						list = gts_bb_tree_overlap(mesh->gdata.bbt, sb);
						//if (list == NULL) continue;
						while (list) {
							b = GTS_BBOX(list->data);
							point[edge] = SegmentTriangleIntersection(segments[edge], GTS_TRIANGLE(b->bounded));
							if (point[edge]) {
								break;
							}
							list = list->next;
						}
					}

					edge = 1; // 3 e 4
					if(true){
						node1 = elem3->nodes[2].id;
						node2 = elem4->nodes[5].id;

						v1 = gts_vertex_new(gts_vertex_class(), coords[node1 * 3 + 0], coords[node1 * 3 + 1], coords[node1 * 3 + 2]);
						v2 = gts_vertex_new(gts_vertex_class(), coords[node2 * 3 + 0], coords[node2 * 3 + 1], coords[node2 * 3 + 2]);

						point[edge] = NULL;
						segments[edge] = gts_segment_new(gts_segment_class(), v1, v2);
						sb = gts_bbox_segment(gts_bbox_class(), segments[edge]);
						list = gts_bb_tree_overlap(mesh->gdata.bbt, sb);
						//if (list == NULL) continue;
						while (list) {
							b = GTS_BBOX(list->data);
							point[edge] = SegmentTriangleIntersection(segments[edge], GTS_TRIANGLE(b->bounded));
							if (point[edge]) {
								break;
							}
							list = list->next;
						}
					}

					int np[2];
					int count_p = 0;
					for(int ip = 0; ip<2;ip++){
						if(point[ip]!=NULL){
							np[count_p] = ip;
							count_p++;
						}
					}

					if(count_p!=0){
						double xx, yy,zz;
						xx = 0;
						yy = 0;
						zz = 0;
						for(int ip = 0; ip<count_p;ip++){
							xx = point[np[ip]]->x+xx;
							yy = point[np[ip]]->y+yy;
							zz = point[np[ip]]->z+zz;
						}
						xx = xx/count_p;
						yy = yy/count_p;
						zz = zz/count_p;

						aux.push_back(elem0->nodes[6].id);
						aux.push_back(xx);
						aux.push_back(yy);
						aux.push_back(zz);

					}else{
						printf("Tu é mais burro que um tijolo 6 furos\n");
					}
				}

				if(oct->id[1]!=-1 && oct->id[5]!=-1){
					//essa é para y...
					octant_t* elem0 = (octant_t*)sc_array_index(&mesh->elements,oct->id[0]);
					octant_t* elem1 = (octant_t*)sc_array_index(&mesh->elements,oct->id[1]);
					octant_t* elem4 = (octant_t*)sc_array_index(&mesh->elements,oct->id[4]);
					octant_t* elem5 = (octant_t*)sc_array_index(&mesh->elements,oct->id[5]);

					edge = 0; // 0 e 5
					if(true){
						node1 = elem0->nodes[3].id;
						node2 = elem5->nodes[6].id;

						v1 = gts_vertex_new(gts_vertex_class(), coords[node1 * 3 + 0], coords[node1 * 3 + 1], coords[node1 * 3 + 2]);
						v2 = gts_vertex_new(gts_vertex_class(), coords[node2 * 3 + 0], coords[node2 * 3 + 1], coords[node2 * 3 + 2]);

						point[edge] = NULL;
						segments[edge] = gts_segment_new(gts_segment_class(), v1, v2);
						sb = gts_bbox_segment(gts_bbox_class(), segments[edge]);
						list = gts_bb_tree_overlap(mesh->gdata.bbt, sb);
						//if (list == NULL) continue;
						while (list) {
							b = GTS_BBOX(list->data);
							point[edge] = SegmentTriangleIntersection(segments[edge], GTS_TRIANGLE(b->bounded));
							if (point[edge]) {
								break;
							}
							list = list->next;
						}
					}

					edge = 1; // 1 e 4
					if(true){
						node1 = elem1->nodes[2].id;
						node2 = elem4->nodes[7].id;

						v1 = gts_vertex_new(gts_vertex_class(), coords[node1 * 3 + 0], coords[node1 * 3 + 1], coords[node1 * 3 + 2]);
						v2 = gts_vertex_new(gts_vertex_class(), coords[node2 * 3 + 0], coords[node2 * 3 + 1], coords[node2 * 3 + 2]);

						point[edge] = NULL;
						segments[edge] = gts_segment_new(gts_segment_class(), v1, v2);
						sb = gts_bbox_segment(gts_bbox_class(), segments[edge]);
						list = gts_bb_tree_overlap(mesh->gdata.bbt, sb);
						//if (list == NULL) continue;
						while (list) {
							b = GTS_BBOX(list->data);
							point[edge] = SegmentTriangleIntersection(segments[edge], GTS_TRIANGLE(b->bounded));
							if (point[edge]) {
								break;
							}
							list = list->next;
						}
					}

					int np[2];
					int count_p = 0;
					for(int ip = 0; ip<2;ip++){
						if(point[ip]!=NULL){
							np[count_p] = ip;
							count_p++;
						}
					}

					if(count_p!=0){
						double xx, yy,zz;
						xx = 0;
						yy = 0;
						zz = 0;
						for(int ip = 0; ip<count_p;ip++){
							xx = point[np[ip]]->x+xx;
							yy = point[np[ip]]->y+yy;
							zz = point[np[ip]]->z+zz;
						}
						xx = xx/count_p;
						yy = yy/count_p;
						zz = zz/count_p;

						aux.push_back(elem0->nodes[6].id);
						aux.push_back(xx);
						aux.push_back(yy);
						aux.push_back(zz);

					}else{
						printf("Tu é mais burro que um tijolo 6 furos\n");
					}
				}

			}

			if(oc_count==2){
				printf("Trabalha o vagabundo\n");
			}
			/*
			point_v = NULL;
			int node1 = elem0->nodes[2].id;
			int node2 = elem4->nodes[6].id;

			GtsVertex *v1 = gts_vertex_new(gts_vertex_class(), coords[node1 * 3], coords[node1 * 3 + 1], coords[node1 * 3 + 2]+0.05*coords[node1 * 3 + 2]);
			GtsVertex *v2 = gts_vertex_new(gts_vertex_class(), coords[node2 * 3], coords[node2 * 3 + 1], coords[node2 * 3 + 2]-0.05*coords[node2 * 3 + 2]);

			segments_v = gts_segment_new(gts_segment_class(), v1, v2);
			GtsBBox *sb = gts_bbox_segment(gts_bbox_class(), segments_v);
			GSList* list = gts_bb_tree_overlap(mesh->gdata.bbt, sb);
			//if (list == NULL) continue;
			while (list) {
				GtsBBox *b = GTS_BBOX(list->data);
				point_v = SegmentTriangleIntersection(segments_v, GTS_TRIANGLE(b->bounded));
				if (point_v) {
					break;
				}
				list = list->next;
			}

			if(point_v!=NULL){
				int cnode = elem0->nodes[6].id;
				nodes_b_mat.push_back(cnode);
				coords[3*cnode+0] = point_v->x;
				coords[3*cnode+1] = point_v->y;
				coords[3*cnode+2] = point_v->z;
			}else{
				printf("Sou um incompetente que não sei como tratar um problema basico\n");
			}
			 */
		}
	}

	for(int iel = 0 ; iel<(aux.size()/4); iel++){
		int node = aux[4*iel+0];
		nodes_b_mat.push_back(node);
		coords[3*node+0] = aux[4*iel+1];
		coords[3*node+1] = aux[4*iel+2];
		coords[3*node+2] = aux[4*iel+3];
		//count++;
	}

	//creating a hash to remove duplicated nodes
	sc_hash_array_t* hash_FixedNodes = sc_hash_array_new(sizeof (node_t), edge_hash_fn, edge_equal_fn, &clamped);
	size_t position;
	node_t *r;
	node_t key;

	for(int ii = 0;ii < nodes_b_mat.size();ii++){
		//printf("No do caralho da porra:%d\n",nodes_b_mat[ii]);
		key.coord[0] = coords[3*nodes_b_mat[ii]+0];
		key.coord[1] = coords[3*nodes_b_mat[ii]+1];
		key.coord[2] = coords[3*nodes_b_mat[ii]+2];
		key.node_id = nodes_b_mat[ii];

		//printf("no fixo, id dele é:%d\n",nodes_b_mat[ii]);

		r = (node_t*) sc_hash_array_insert_unique(hash_FixedNodes, &key, &position);
		if (r != NULL) {
			r->coord[0] = key.coord[0];
			r->coord[1] = key.coord[1];
			r->coord[2] = key.coord[2];
			r->node_id = nodes_b_mat[ii];
		} else {
			//printf("Ta me achando com cara de palhaco?");
		}
	}

	nodes_b_mat.clear();
	for(int ii = 0;ii < hash_FixedNodes->a.elem_count ;ii++){
		node_t* node = (node_t*) sc_array_index (&hash_FixedNodes->a, ii);
		nodes_b_mat.push_back(node->node_id);
		//printf("nos que estao na hash:%d\n",nodes_b_mat[ii] );
	}
	printf("Total de %d nos fixos...\n",nodes_b_mat.size());


}

void FreeMovableNodes(hexa_tree_t* mesh){

	for(int iel = 0; iel < mesh->oct.elem_count; iel++){
		octree_t *oct = (octree_t*) sc_array_index(&mesh->oct, iel);

		for(int i=0; i<8; i++){
			if(oct->id[i]!=-1){
				octant_t *elem = (octant_t*) sc_array_index(&mesh->elements, oct->id[i]);

				//liberando as arestas agora
				//face superior
				if(oct->edge[0]){
					if(i==0){
						elem->nodes[1].fixed=2;
					}
					if(i==1){
						elem->nodes[0].fixed=2;
					}
				}
				if(oct->edge[1]){
					if(i==1){
						elem->nodes[2].fixed=2;
					}
					if(i==2){
						elem->nodes[1].fixed=2;
					}
				}
				if(oct->edge[2]){
					if(i==2){
						elem->nodes[3].fixed=2;
					}
					if(i==3){
						elem->nodes[2].fixed=2;
					}
				}
				if(oct->edge[3]){
					if(i==3){
						elem->nodes[0].fixed=2;
					}
					if(i==0){
						elem->nodes[3].fixed=2;
					}
				}

				//face inferior
				if(oct->edge[8]){
					if(i==4){
						elem->nodes[5].fixed=2;
					}
					if(i==5){
						elem->nodes[4].fixed=2;
					}
				}
				if(oct->edge[9]){
					if(i==5){
						elem->nodes[6].fixed=2;
					}
					if(i==6){
						elem->nodes[5].fixed=2;
					}
				}
				if(oct->edge[10]){
					if(i==6){
						elem->nodes[7].fixed=2;
					}
					if(i==7){
						elem->nodes[6].fixed=2;
					}
				}
				if(oct->edge[11]){
					if(i==7){
						elem->nodes[4].fixed=2;
					}
					if(i==4){
						elem->nodes[7].fixed=2;
					}
				}


				if(oct->edge[4]){
					if(i==0){
						elem->nodes[4].fixed=2;
					}
					if(i==4){
						elem->nodes[0].fixed=2;
					}
				}
				if(oct->edge[5]){
					if(i==1){
						elem->nodes[5].fixed=2;
					}
					if(i==5){
						elem->nodes[1].fixed=2;
					}
				}
				if(oct->edge[6]){
					if(i==2){
						elem->nodes[6].fixed=2;
					}
					if(i==6){
						elem->nodes[2].fixed=2;
					}
				}
				if(oct->edge[7]){
					if(i==3){
						elem->nodes[7].fixed=2;
					}
					if(i==7){
						elem->nodes[3].fixed=2;
					}
				}

				//liberando as faces agora
				if(oct->face[0]){
					if(i==0){
						elem->nodes[7].fixed=2;
					}
					if(i==4){
						elem->nodes[3].fixed=2;
					}
					if(i==3){
						elem->nodes[4].fixed=2;
					}
					if(i==7){
						elem->nodes[0].fixed=2;
					}
				}

				if(oct->face[1]){
					if(i==1){
						elem->nodes[6].fixed=2;
					}
					if(i==5){
						elem->nodes[2].fixed=2;
					}
					if(i==2){
						elem->nodes[5].fixed=2;
					}
					if(i==6){
						elem->nodes[1].fixed=2;
					}
				}

				if(oct->face[2]){
					if(i==0){
						elem->nodes[6].fixed=2;
					}
					if(i==1){
						elem->nodes[5].fixed=2;
					}
					if(i==4){
						elem->nodes[1].fixed=2;
					}
					if(i==5){
						elem->nodes[0].fixed=2;
					}
				}

				if(oct->face[3]){
					if(i==2){
						elem->nodes[7].fixed=2;
					}
					if(i==3){
						elem->nodes[6].fixed=2;
					}
					if(i==6){
						elem->nodes[3].fixed=2;
					}
					if(i==7){
						elem->nodes[2].fixed=2;
					}
				}

				if(oct->face[4]){
					if(i==4){
						elem->nodes[6].fixed=2;
					}
					if(i==5){
						elem->nodes[7].fixed=2;
					}
					if(i==6){
						elem->nodes[4].fixed=2;
					}
					if(i==7){
						elem->nodes[5].fixed=2;
					}
				}

				if(oct->face[5]){
					if(i==0){
						elem->nodes[2].fixed=2;
					}
					if(i==1){
						elem->nodes[3].fixed=2;
					}
					if(i==2){
						elem->nodes[0].fixed=2;
					}
					if(i==3){
						elem->nodes[1].fixed=2;
					}
				}


				if(oct->face[0] || oct->face[1] || oct->face[2] || oct->face[3] || oct->face[4] || oct->face[5]){
					//no central
					if(i==0){
						elem->nodes[6].fixed=2;
					}
					if(i==1){
						elem->nodes[7].fixed=2;
					}
					if(i==2){
						elem->nodes[4].fixed=2;
					}
					if(i==3){
						elem->nodes[5].fixed=2;
					}

					if(i==4){
						elem->nodes[2].fixed=2;
					}
					if(i==5){
						elem->nodes[3].fixed=2;
					}
					if(i==6){
						elem->nodes[0].fixed=2;
					}
					if(i==7){
						elem->nodes[1].fixed=2;
					}
				}
			}
		}
	}
}

void IdentifyMovableNodes(hexa_tree_t* mesh){

	for(int iel = 0; iel < mesh->oct.elem_count; iel++){
		octree_t *oct = (octree_t*) sc_array_index(&mesh->oct, iel);

		for(int i=0; i<8; i++){
			if(oct->id[i]!=-1){
				octant_t *elem = (octant_t*) sc_array_index(&mesh->elements, oct->id[i]);

				//identificando as arestas interceptadas
				//face superior
				if(i==0 || i ==1){
					if(elem->edge[0].ref){
						oct->edge[0]=true;
					}
				}
				if(i==1 || i ==2){
					if(elem->edge[1].ref){
						oct->edge[1]=true;
					}
				}
				if(i==2 || i ==3){
					if(elem->edge[2].ref){
						oct->edge[2]=true;
					}
				}
				if(i==3 || i ==0){
					if(elem->edge[3].ref){
						oct->edge[3]=true;
					}
				}

				//face inferior
				if(i==4 || i ==5){
					if(elem->edge[8].ref){
						oct->edge[8]=true;
					}
				}
				if(i==5 || i ==6){
					if(elem->edge[9].ref){
						oct->edge[9]=true;
					}
				}
				if(i==6 || i ==7){
					if(elem->edge[10].ref){
						oct->edge[10]=true;
					}
				}
				if(i==7 || i ==4){
					if(elem->edge[11].ref){
						oct->edge[12]=true;
					}
				}

				//arestas verticais
				if(i==0 || i ==4){
					if(elem->edge[4].ref){
						oct->edge[4]=true;
					}
				}
				if(i==1 || i ==5){
					if(elem->edge[5].ref){
						oct->edge[5]=true;
					}
				}
				if(i==2 || i ==6){
					if(elem->edge[6].ref){
						oct->edge[6]=true;
					}
				}
				if(i==3 || i ==7){
					if(elem->edge[7].ref){
						oct->edge[7]=true;
					}
				}
			}
		}

		//identificando as faces
		//x-
		if(oct->edge[4] || oct->edge[11] || oct->edge[7] || oct->edge[3]){
			oct->face[0]=true;
		}
		//x+
		if(oct->edge[5] || oct->edge[1] || oct->edge[6] || oct->edge[9]){
			oct->face[1]=true;
		}
		//y-
		if(oct->edge[0] || oct->edge[5] || oct->edge[8] || oct->edge[4]){
			oct->face[2]=true;
		}
		//y+
		if(oct->edge[2] || oct->edge[6] || oct->edge[10] || oct->edge[7]){
			oct->face[3]=true;
		}
		//z-
		if(oct->edge[8] || oct->edge[9] || oct->edge[10] || oct->edge[11]){
			oct->face[4]=true;
		}
		//z+
		if(oct->edge[0] || oct->edge[1] || oct->edge[2] || oct->edge[3]){
			oct->face[5]=true;
		}
	}
}

void DoOctree(hexa_tree_t* mesh){

	bool               clamped = true;
	size_t              position;
	sc_hash_array_t    *indep_vertex;
	indep_vertex    = (sc_hash_array_t *)sc_hash_array_new(sizeof (octant_vertex_t), vertex_hash_id, vertex_equal_id, &clamped);

	/////////////////
	// create the vertex structure
	for (int iel = 0; iel < mesh->elements.elem_count; ++iel) {

		octant_t *elem = (octant_t*) sc_array_index(&mesh->elements, iel);

		for (int ino = 0; ino < 8; ino++){

			octant_vertex_t key;
			key.id = elem->nodes[ino].id;

			octant_vertex_t* vert = (octant_vertex_t*) sc_hash_array_insert_unique (indep_vertex, &key, &position);
			if(vert != NULL){
				vert->id = elem->nodes[ino].id;
				vert->list_elem = 1;
				vert->elem[vert->list_elem-1] = elem->id;
			}else{
				vert = (octant_vertex_t*) sc_array_index(&indep_vertex->a, position);
				vert->elem[vert->list_elem] = elem->id;
				vert->list_elem++;
			}
		}
	}

	//creating the octree structure
	sc_array_init(&mesh->oct, sizeof(octree_t));
	sc_array_reset(&mesh->oct);

	//hash for the elements
	sc_hash_array_t * hash = (sc_hash_array_t *)sc_hash_array_new(sizeof (octant_t), el_hash_id, el_equal_id, &clamped);
	octant_t* r;
	bool elem_lookup;
	bool vert_lookup;

	//find the elements...
	for (int iel = 0; iel < mesh->elements.elem_count; ++iel) {
		octant_t *cut = (octant_t*) sc_array_index(&mesh->elements, iel);

		size_t position;
		octant_t key;
		key.id = cut->id;

		elem_lookup = sc_hash_array_lookup(hash, &key, &position);

		int temp_id[8];
		bool temp_cut = false;;

		for(int ii = 0; ii<8 ; ii++){
			temp_id[ii] = -1;
		}

		if(elem_lookup){

		}else{

			octant_vertex_t key1;
			size_t position1;
			key1.id = cut->nodes[6].id;

			vert_lookup = sc_hash_array_lookup(indep_vertex, &key1, &position1);

			if(vert_lookup){
				octant_vertex_t* vert = (octant_vertex_t*) sc_array_index(&indep_vertex->a, position1);

				//printf("entrei caraio\n");
				for(int iele = 0; iele < vert->list_elem; iele++){

					int test = vert->elem[iele];

					octant_t *elem = (octant_t*) sc_array_index(&mesh->elements, test);
					key.id = elem->id;

					if(cut->id==elem->id){
						r = (octant_t*) sc_hash_array_insert_unique(hash,&key,&position);
						temp_id[0] = elem->id;
						if(r!=NULL){
							r->id= elem->id;
						}else{
							r =  (octant_t*) sc_array_index(&hash->a,position);
						}
						if(elem->pad==-1 || cut->pad == -1){
							temp_cut = true;
						}
					}//x+

					if(cut->nodes[1].id == elem->nodes[0].id && cut->nodes[2].id == elem->nodes[3].id &&
							cut->nodes[5].id==elem->nodes[4].id && cut->nodes[6].id==elem->nodes[7].id){
						temp_id[1] = elem->id;
						r = (octant_t*) sc_hash_array_insert_unique(hash,&key,&position);
						if(r!=NULL){
							r->id= elem->id;
						}else{
							r =  (octant_t*) sc_array_index(&hash->a,position);
						}
						if(elem->pad==-1 || cut->pad == -1){
							temp_cut = true;
						}
					}//x+y+

					if(cut->nodes[2].id==elem->nodes[0].id && cut->nodes[6].id==elem->nodes[4].id){
						temp_id[2]= elem->id;
						r = (octant_t*) sc_hash_array_insert_unique(hash,&key,&position);
						if(r!=NULL){
							r->id= elem->id;
						}else{
							r =  (octant_t*) sc_array_index(&hash->a,position);
						}
						if(elem->pad==-1 || cut->pad == -1){
							temp_cut = true;
						}
					}//y+

					if(cut->nodes[3].id==elem->nodes[0].id && cut->nodes[2].id==elem->nodes[1].id &&
							cut->nodes[6].id==elem->nodes[5].id && cut->nodes[7].id==elem->nodes[4].id){
						temp_id[3]= elem->id;
						r = (octant_t*) sc_hash_array_insert_unique(hash,&key,&position);
						if(r!=NULL){
							r->id= elem->id;
						}else{
							r =  (octant_t*) sc_array_index(&hash->a,position);
						}
						if(elem->pad==-1 || cut->pad == -1){
							temp_cut = true;
						}
					}//z-

					if(cut->nodes[4].id==elem->nodes[0].id && cut->nodes[5].id==elem->nodes[1].id &&
							cut->nodes[6].id==elem->nodes[2].id && cut->nodes[7].id==elem->nodes[3].id){
						temp_id[4]= elem->id;
						r = (octant_t*) sc_hash_array_insert_unique(hash,&key,&position);
						if(r!=NULL){
							r->id= elem->id;
						}else{
							r =  (octant_t*) sc_array_index(&hash->a,position);
						}
						if(elem->pad==-1 || cut->pad == -1){
							temp_cut = true;
						}
					}//z- x+

					if(cut->nodes[5].id==elem->nodes[0].id && cut->nodes[6].id==elem->nodes[3].id){
						temp_id[5]= elem->id;
						r = (octant_t*) sc_hash_array_insert_unique(hash,&key,&position);
						if(r!=NULL){
							r->id= elem->id;
						}else{
							r =  (octant_t*) sc_array_index(&hash->a,position);
						}
						if(elem->pad==-1 || cut->pad == -1){
							temp_cut = true;
						}
					}//z- x+ y+

					if(cut->nodes[6].id==elem->nodes[0].id){
						temp_id[6]= elem->id;
						r = (octant_t*) sc_hash_array_insert_unique(hash,&key,&position);
						if(r!=NULL){
							r->id= elem->id;
						}else{
							r =  (octant_t*) sc_array_index(&hash->a,position);
						}
						if(elem->pad==-1 || cut->pad == -1){
							temp_cut = true;
						}
					} //z- y+

					if(cut->nodes[6].id==elem->nodes[1].id && cut->nodes[7].id==elem->nodes[0].id){
						temp_id[7]= elem->id;
						r = (octant_t*) sc_hash_array_insert_unique(hash,&key,&position);
						if(r!=NULL){
							r->id= elem->id;
						}else{
							r =  (octant_t*) sc_array_index(&hash->a,position);
						}
						if(elem->pad==-1 || cut->pad == -1){
							temp_cut = true;
						}
					}
				}
			}
		}

		if(temp_cut){
			//printf("Octree: ");
			//create the octree
			octree_t * oc = (octree_t*) sc_array_push(&mesh->oct);
			// fill with -1
			for(int i = 0; i<8; i++){
				oc->id[i] = temp_id[i];
				//printf("El: %d\n", temp_id[i]);
				oc->mat[i] = -1;
			}
		}
	}

	if(true){
		int count = 0;
		//seleciona apenas os octantes cortados
		for(int iel = 0; iel< mesh->oct.elem_count; iel++){
			octree_t *oc = (octree_t*) sc_array_index(&mesh->oct, iel);
			//inicia com arestas e faces como nao cortadas
			for(int i = 0; i<12;i++){
				oc->edge[i] = false;
			}
			for(int i = 0; i<6;i++){
				oc->face[i] = false;
			}

			for(int i = 0; i<8; i++){
				//printf("Octree numero:%d, elemento numero:%d\n",iel,oc->id[i]);
				if(oc->id[i]!=-1){
					octant_t *elem = (octant_t*) sc_array_index(&mesh->elements, oc->id[i]);
					elem->pad = oc->id[0]+1;
					//fix nodes
					if(i==0){
						//printf("Entrei para fixar o 0 para o elemento numero:%d\n",elem->id);
						elem->nodes[0].fixed=1;
					}else if(i==1){
						//printf("Entrei para fixar o 1 para o elemento numero:%d\n",elem->id);
						elem->nodes[1].fixed=1;
					}else if(i==2){
						//printf("Entrei para fixar o 2 para o elemento numero:%d\n",elem->id);
						elem->nodes[2].fixed=1;
					}else if(i==3){
						//printf("Entrei para fixar o 3 para o elemento numero:%d\n",elem->id);
						elem->nodes[3].fixed=1;
					}else if(i==4){
						//printf("Entrei para fixar o 4 para o elemento numero:%d\n",elem->id);
						elem->nodes[4].fixed=1;
					}else if(i==5){
						//printf("Entrei para fixar o 5 para o elemento numero:%d\n",elem->id);
						elem->nodes[5].fixed=1;
					}else if(i==6){
						//printf("Entrei para fixar o 6 para o elemento numero:%d\n",elem->id);
						elem->nodes[6].fixed=1;
					}else if(i==7){
						//printf("Entrei para fixar o 7 para o elemento numero:%d\n",elem->id);
						elem->nodes[7].fixed=1;
					}
				}
			}
		}
	}

	if(false){
		for(int iel = 0; iel< mesh->oct.elem_count; iel++){
			octree_t *oc = (octree_t*) sc_array_index(&mesh->oct, iel);
			printf("Octree numero:%d, foi cortado?:%d\n",iel,oc->cut);
			printf("Ids dos elementos: ");
			for(int i = 0; i<8; i++){
				printf("%d ",oc->id[i]);
			}
			printf("\n");
			printf("Edges do octree: ");

			for(int i = 0; i<12; i++){
				printf("%d ",oc->edge[i]);
			}
			printf("\n");
			printf("Faces do octree: ");

			for(int i = 0; i<6; i++){
				printf("%d ",oc->face[i]);
			}
			printf("\n");
		}
	}
}

void MovingNodes(hexa_tree_t* mesh, std::vector<double>& coords, std::vector<int>& nodes_b_mat, const char* surface) {

	int8_t* flag_nodes = (int8_t*) malloc(sizeof (int8_t) * mesh->local_n_nodes);

	time_t tstart, tend;

	memset(flag_nodes, 0, sizeof (int8_t) * mesh->local_n_nodes);

	tstart = time(0);
	printf(" Building the octree structure...\n");
	DoOctree(mesh);
	tend = time(0);
//	cout << "Time in DoOctree "<< difftime(tend, tstart) <<" second(s)."<< endl;


	tstart = time(0);
	printf(" Identifying the movable nodes...\n");
	IdentifyMovableNodes(mesh);
	tend = time(0);
	//cout << "Time in IdentifyMovableNodes "<< difftime(tend, tstart) <<" second(s)."<< endl;

	tstart = time(0);
	printf(" Free the movable nodes...\n");
	FreeMovableNodes(mesh);
	tend = time(0);
	//cout << "Time in FreeMovableNodes "<< difftime(tend, tstart) <<" second(s)."<< endl;

	tstart = time(0);
	printf(" Make the projection of the nodes into the surface...\n");
	nodes_b_mat.clear();
	ProjectFreeNodes(mesh,coords,nodes_b_mat);
	tend = time(0);
	//cout << "Time in ProjectFreeNodes "<< difftime(tend, tstart) <<" second(s)."<< endl;

	for (int i = 0; i < mesh->local_n_nodes; i++) {
		flag_nodes[i]=0;
	}

	//verifica a criacao os nos livres e fixos na criacao dos octrees
	if(true){
		//std::cout <<  mesh->elements.elem_count << " coisinhas para mover" << std::endl;
		for (int iel = 0; iel < mesh->elements.elem_count; ++iel) {
			octant_t *elem = (octant_t*) sc_array_index(&mesh->elements, iel);
			for(int ino = 0; ino<8;ino++){
				if(flag_nodes[elem->nodes[ino].id]==2){

				}else{
					flag_nodes[elem->nodes[ino].id] = elem->nodes[ino].fixed;
				}
			}
		}

		for (int i = 0; i < mesh->local_n_nodes; i++) {
			mesh->part_nodes[i] = flag_nodes[i];
		}
	}

	//verifica se os nos movidos estao fixos...
	for(int i = 0; i<nodes_b_mat.size();i++){
		mesh->part_nodes[nodes_b_mat[i]] = 1;
	}

	//gts_bb_tree_destroy(bbt_bathymetry, TRUE);

	//gts_object_destroy(GTS_OBJECT(bathymetry));

	free(flag_nodes);
}
