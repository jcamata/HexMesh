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

#include <ctime>

/*
typedef struct {
	GPtrArray * array;
} ListOfPoints;

static void InsertPoint(GtsPoint * p, ListOfPoints * lp) {
	g_ptr_array_add(lp->array, p);
}

bool is_point_over_surface(GtsPoint * p, GNode * tree);
 */
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

GtsPoint* FoundInterception(hexa_tree_t* mesh,std::vector<double>& coords,int node1, int node2){
	GtsPoint *point = NULL;
	GtsVertex *v1 = gts_vertex_new(gts_vertex_class(), coords[node1 * 3], coords[node1 * 3 + 1], coords[node1 * 3 + 2]);
	GtsVertex *v2 = gts_vertex_new(gts_vertex_class(), coords[node2 * 3], coords[node2 * 3 + 1], coords[node2 * 3 + 2]);

	GtsSegment *segments = gts_segment_new(gts_segment_class(), v1, v2);
	GtsBBox *sb = gts_bbox_segment(gts_bbox_class(), segments);
	GSList* list = gts_bb_tree_overlap(mesh->gdata.bbt, sb);
	//if (list == NULL) continue;
	while (list) {
		GtsBBox *b = GTS_BBOX(list->data);
		point = SegmentTriangleIntersection(segments, GTS_TRIANGLE(b->bounded));
		if (point) {
			break;
		}
		list = list->next;
	}
	return point;
}

void ProjectFreeNodes(hexa_tree_t* mesh,std::vector<double>& coords, std::vector<int>& nodes_b_mat){

	bool deb = false;
	int Edge2GNode[12][2]={0};
	GtsSegment * segments[12]={0};
	GtsPoint * point[12]={NULL};
	bool clamped = true;
	std::vector<double> aux;

	int coord_count = 0;

	//full octree
	if(true){
		//only for the complet octrees
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
								point[edge] = SegmentTriangleIntersectionCgal(segments[edge], GTS_TRIANGLE(b->bounded));
								if (point[edge]) {
									//printf("%f %f %f\n",point[edge]->x,point[edge]->y,point[edge]->z);
									//GtsPoint* pp = SegmentTriangleIntersection(segments[edge], GTS_TRIANGLE(b->bounded));
									//printf("%f %f %f\n",pp->x,pp->y,pp->z);
									//printf("break no %d\n",edge);
									break;
								}
								list = list->next;
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
			}
		}

		//TODO add a element counter to avoid redo the coords change
		//Changing the coords vector
		for(int iel = 0 ; iel<(aux.size()/4); iel++){
			int node = aux[4*iel+0];
			nodes_b_mat.push_back(node);
			coords[3*node+0] = aux[4*iel+1];
			coords[3*node+1] = aux[4*iel+2];
			coords[3*node+2] = aux[4*iel+3];
			coord_count = coord_count+4;
		}

		//perform the mean between the lines in order to obtain
		//the value of point that lies in the surface
		if(true){
			for (int ioc = 0; ioc < mesh->oct.elem_count; ++ioc) {
				octree_t* oct = (octree_t*)sc_array_index(&mesh->oct,ioc);

				int oc_count=0;
				for(int i =0;i<8;i++){
					if(oct->id[i]!=-1) {
						oc_count++;
					}
				}

				if(oc_count == 8){
					octant_t* elem0 = (octant_t*)sc_array_index(&mesh->elements,oct->id[0]);
					octant_t* elem2 = (octant_t*)sc_array_index(&mesh->elements,oct->id[2]);
					octant_t* elem5 = (octant_t*)sc_array_index(&mesh->elements,oct->id[5]);
					octant_t* elem7 = (octant_t*)sc_array_index(&mesh->elements,oct->id[7]);

					//face x-
					if(oct->face[0]){

						double xx = 0;
						double yy = 0;
						double zz = 0;
						int count = 0;
						int nodeId;

						if(oct->edge[3]){
							nodeId = elem0->nodes[3].id;
							xx += coords[3*nodeId+0];
							yy += coords[3*nodeId+1];
							zz += coords[3*nodeId+2];
							count++;
						}
						if(oct->edge[4]){
							nodeId = elem0->nodes[4].id;
							xx += coords[3*nodeId+0];
							yy += coords[3*nodeId+1];
							zz += coords[3*nodeId+2];
							count++;
						}
						if(oct->edge[7]){
							nodeId = elem7->nodes[3].id;
							xx += coords[3*nodeId+0];
							yy += coords[3*nodeId+1];
							zz += coords[3*nodeId+2];
							count++;
						}
						if(oct->edge[11]){
							nodeId = elem7->nodes[4].id;
							xx += coords[3*nodeId+0];
							yy += coords[3*nodeId+1];
							zz += coords[3*nodeId+2];
							count++;
						}

						//printf("Contador:%d\n",count);
						aux.push_back(elem0->nodes[7].id);
						aux.push_back(double(xx/count));
						aux.push_back(double(yy/count));
						aux.push_back(double(zz/count));

					}

					//face x+
					if(oct->face[1]){

						double xx = 0;
						double yy = 0;
						double zz = 0;
						int count = 0;
						int nodeId;

						if(oct->edge[1]){
							//printf("edge1\n");
							nodeId = elem2->nodes[1].id;
							xx += coords[3*nodeId+0];
							yy += coords[3*nodeId+1];
							zz += coords[3*nodeId+2];
							count++;
						}
						if(oct->edge[5]){
							//printf("edge5\n");
							nodeId = elem5->nodes[1].id;
							xx += coords[3*nodeId+0];
							yy += coords[3*nodeId+1];
							zz += coords[3*nodeId+2];
							count++;
						}
						if(oct->edge[6]){
							//printf("edge6\n");
							nodeId = elem2->nodes[6].id;
							xx += coords[3*nodeId+0];
							yy += coords[3*nodeId+1];
							zz += coords[3*nodeId+2];
							count++;
						}
						if(oct->edge[9]){
							//printf("edge9\n");
							nodeId = elem5->nodes[6].id;
							xx += coords[3*nodeId+0];
							yy += coords[3*nodeId+1];
							zz += coords[3*nodeId+2];
							count++;
						}

						//printf("Contador:%d\n",count);
						aux.push_back(elem5->nodes[2].id);
						aux.push_back(double(xx/count));
						aux.push_back(double(yy/count));
						aux.push_back(double(zz/count));

					}

					//face y-
					if(oct->face[2]){

						double xx = 0;
						double yy = 0;
						double zz = 0;
						int count = 0;
						int nodeId;

						if(oct->edge[0]){
							nodeId = elem0->nodes[1].id;
							xx += coords[3*nodeId+0];
							yy += coords[3*nodeId+1];
							zz += coords[3*nodeId+2];
							count++;
						}
						if(oct->edge[4]){
							nodeId = elem0->nodes[4].id;
							xx += coords[3*nodeId+0];
							yy += coords[3*nodeId+1];
							zz += coords[3*nodeId+2];
							count++;
						}
						if(oct->edge[5]){
							nodeId = elem5->nodes[1].id;
							xx += coords[3*nodeId+0];
							yy += coords[3*nodeId+1];
							zz += coords[3*nodeId+2];
							count++;
						}
						if(oct->edge[8]){
							nodeId = elem5->nodes[4].id;
							xx += coords[3*nodeId+0];
							yy += coords[3*nodeId+1];
							zz += coords[3*nodeId+2];
							count++;
						}

						//printf("Contador:%d\n",count);
						aux.push_back(elem0->nodes[5].id);
						aux.push_back(double(xx/count));
						aux.push_back(double(yy/count));
						aux.push_back(double(zz/count));

					}

					//face y+
					if(oct->face[3]){

						double xx = 0;
						double yy = 0;
						double zz = 0;
						int count = 0;
						int nodeId;

						if(oct->edge[2]){
							nodeId = elem2->nodes[3].id;
							xx += coords[3*nodeId+0];
							yy += coords[3*nodeId+1];
							zz += coords[3*nodeId+2];
							count++;
						}
						if(oct->edge[6]){
							nodeId = elem2->nodes[6].id;
							xx += coords[3*nodeId+0];
							yy += coords[3*nodeId+1];
							zz += coords[3*nodeId+2];
							count++;
						}
						if(oct->edge[7]){
							nodeId = elem7->nodes[3].id;
							xx += coords[3*nodeId+0];
							yy += coords[3*nodeId+1];
							zz += coords[3*nodeId+2];
							count++;
						}
						if(oct->edge[10]){
							nodeId = elem7->nodes[6].id;
							xx += coords[3*nodeId+0];
							yy += coords[3*nodeId+1];
							zz += coords[3*nodeId+2];
							count++;
						}

						//printf("Contador:%d\n",count);
						aux.push_back(elem2->nodes[7].id);
						aux.push_back(double(xx/count));
						aux.push_back(double(yy/count));
						aux.push_back(double(zz/count));

					}

					//face z-
					if(oct->face[4]){

						double xx = 0;
						double yy = 0;
						double zz = 0;
						int count = 0;
						int nodeId;

						if(oct->edge[8]){
							nodeId = elem5->nodes[4].id;
							xx += coords[3*nodeId+0];
							yy += coords[3*nodeId+1];
							zz += coords[3*nodeId+2];
							count++;
						}
						if(oct->edge[9]){
							nodeId = elem5->nodes[6].id;
							xx += coords[3*nodeId+0];
							yy += coords[3*nodeId+1];
							zz += coords[3*nodeId+2];
							count++;
						}
						if(oct->edge[10]){
							nodeId = elem7->nodes[6].id;
							xx += coords[3*nodeId+0];
							yy += coords[3*nodeId+1];
							zz += coords[3*nodeId+2];
							count++;
						}
						if(oct->edge[11]){
							nodeId = elem7->nodes[4].id;
							xx += coords[3*nodeId+0];
							yy += coords[3*nodeId+1];
							zz += coords[3*nodeId+2];
							count++;
						}

						//printf("Contador:%d\n",count);
						aux.push_back(elem5->nodes[7].id);
						aux.push_back(double(xx/count));
						aux.push_back(double(yy/count));
						aux.push_back(double(zz/count));

					}

					//face z+
					if(oct->face[5]){

						double xx = 0;
						double yy = 0;
						double zz = 0;
						int count = 0;
						int nodeId;

						if(oct->edge[0]){
							nodeId = elem0->nodes[1].id;
							xx += coords[3*nodeId+0];
							yy += coords[3*nodeId+1];
							zz += coords[3*nodeId+2];
							count++;
						}
						if(oct->edge[1]){
							nodeId = elem2->nodes[1].id;
							xx += coords[3*nodeId+0];
							yy += coords[3*nodeId+1];
							zz += coords[3*nodeId+2];
							count++;
						}
						if(oct->edge[2]){
							nodeId = elem2->nodes[3].id;
							xx += coords[3*nodeId+0];
							yy += coords[3*nodeId+1];
							zz += coords[3*nodeId+2];
							count++;
						}
						if(oct->edge[3]){
							nodeId = elem0->nodes[3].id;
							xx += coords[3*nodeId+0];
							yy += coords[3*nodeId+1];
							zz += coords[3*nodeId+2];
							count++;
						}

						//printf("Contador:%d\n",count);
						aux.push_back(elem0->nodes[2].id);
						aux.push_back(double(xx/count));
						aux.push_back(double(yy/count));
						aux.push_back(double(zz/count));

					}

				}else if(oc_count==4 && false){

				}else if(oc_count ==2 && false) {

				}else{
					printf("Error in move_node.cpp\n some error while moving one surface node\n");
					printf("Please verify the generation of the octree\n");
					//printf("Node %d\n",elem0->nodes[6].id);
					exit (EXIT_FAILURE);
				}

			}
		}

		//Changing the coords vector
		for(int iel = 0 ; iel<(aux.size()/4); iel++){
			int node = aux[4*iel+0];
			nodes_b_mat.push_back(node);
			coords[3*node+0] = aux[4*iel+1];
			coords[3*node+1] = aux[4*iel+2];
			coords[3*node+2] = aux[4*iel+3];
			coord_count = coord_count+4;
		}

		//perform the mean between the faces in order to obtain
		//the value of the central point
		if(true){
			for (int ioc = 0; ioc < mesh->oct.elem_count; ++ioc) {
				octree_t* oct = (octree_t*)sc_array_index(&mesh->oct,ioc);

				int oc_count=0;
				for(int i =0;i<8;i++){
					if(oct->id[i]!=-1) {
						oc_count++;
					}
				}

				if(oc_count == 8){
					octant_t* elem0 = (octant_t*)sc_array_index(&mesh->elements,oct->id[0]);
					octant_t* elem2 = (octant_t*)sc_array_index(&mesh->elements,oct->id[2]);
					octant_t* elem5 = (octant_t*)sc_array_index(&mesh->elements,oct->id[5]);

					double xx = 0;
					double yy = 0;
					double zz = 0;
					int count = 0;
					int nodeId;

					//face x-
					if(oct->face[0]){
						//printf("face0\n");
						nodeId = elem0->nodes[7].id;
						xx += coords[3*nodeId+0];
						yy += coords[3*nodeId+1];
						zz += coords[3*nodeId+2];
						count++;
					}

					//face x+
					if(oct->face[1]){
						//printf("face1\n");
						nodeId = elem2->nodes[5].id;
						xx += coords[3*nodeId+0];
						yy += coords[3*nodeId+1];
						zz += coords[3*nodeId+2];
						count++;
					}

					//face y-
					if(oct->face[2]){
						//printf("face2\n");
						nodeId = elem0->nodes[5].id;
						xx += coords[3*nodeId+0];
						yy += coords[3*nodeId+1];
						zz += coords[3*nodeId+2];
						count++;
					}

					//face y+
					if(oct->face[3]){
						//printf("face3\n");
						nodeId = elem2->nodes[7].id;
						xx += coords[3*nodeId+0];
						yy += coords[3*nodeId+1];
						zz += coords[3*nodeId+2];
						count++;
					}

					//face z-
					if(oct->face[4]){
						//printf("face4\n");
						nodeId = elem5->nodes[7].id;
						xx += coords[3*nodeId+0];
						yy += coords[3*nodeId+1];
						zz += coords[3*nodeId+2];
						count++;
					}

					//face z+
					if(oct->face[5]){
						//printf("face5\n");
						nodeId = elem0->nodes[2].id;
						xx += coords[3*nodeId+0];
						yy += coords[3*nodeId+1];
						zz += coords[3*nodeId+2];
						count++;
					}

					//printf("Contador:%d\n",count);
					aux.push_back(elem0->nodes[6].id);
					aux.push_back(double(xx/count));
					aux.push_back(double(yy/count));
					aux.push_back(double(zz/count));
				}else if(oc_count==4 && false){

				}else if(oc_count ==2 && false) {

				}else{
					printf("Error in move_node.cpp\n some error while moving the central node\n");
					//printf("Please verify the generation of the octree\n");
					//printf("Node %d\n",elem0->nodes[6].id);
					exit (EXIT_FAILURE);
				}

			}
		}

		//Changing the coords vector
		for(int iel = 0 ; iel<(aux.size()/4); iel++){
			int node = aux[4*iel+0];
			nodes_b_mat.push_back(node);
			coords[3*node+0] = aux[4*iel+1];
			coords[3*node+1] = aux[4*iel+2];
			coords[3*node+2] = aux[4*iel+3];
			coord_count=coord_count+4;
		}
	}

	//creating a hash to remove duplicated nodes
	sc_hash_array_t* hash_FixedNodes = sc_hash_array_new(sizeof (node_t), edge_hash_fn, edge_equal_fn, &clamped);
	size_t position;
	node_t *r;
	node_t key;

	for(int ii = 0;ii < nodes_b_mat.size();ii++){
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

		}
	}

	//clean and fix the nodes in the element structure
	for (int ioc = 0; ioc < mesh->oct.elem_count; ++ioc) {
		octree_t* oct = (octree_t*)sc_array_index(&mesh->oct,ioc);
		for(int iel = 0; iel<8; iel++){
			octant_t* elem = (octant_t*)sc_array_index(&mesh->elements,oct->id[iel]);
			for(int ino = 0; ino<8; ino++){
				elem->nodes[ino].fixed = 2;
				//elem->nodes[ino].color = 0;
				int node = elem->nodes[ino].id;
				key.coord[0] = coords[3*node+0];
				key.coord[1] = coords[3*node+1];
				key.coord[2] = coords[3*node+2];
				key.node_id = node;
				bool tre = sc_hash_array_lookup(hash_FixedNodes, &key, &position);
				if(tre){
					elem->nodes[ino].fixed = 1;
				}
			}
		}
	}

	nodes_b_mat.clear();
	for(int ii = 0;ii < hash_FixedNodes->a.elem_count ;ii++){
		node_t* node = (node_t*) sc_array_index (&hash_FixedNodes->a, ii);
		nodes_b_mat.push_back(node->node_id);
	}
	printf(" Total of %d fixed nodes\n",nodes_b_mat.size());

	//add the color of the node (help in the material selection)
	for (int ioc = 0; ioc < mesh->oct.elem_count; ++ioc){
		octree_t* oct = (octree_t*)sc_array_index(&mesh->oct,ioc);
		octant_t *elem[8];
		for(int iel = 0; iel < 8; iel++){
			elem[iel] = (octant_t *)sc_array_index(&mesh->elements,oct->id[iel]);
		}
		if(deb){
			for(int iedge = 0; iedge < 12; iedge++){
				if(oct->edge[iedge]) printf("Para o octree %d a aresta %d foi interceptada\n",ioc,iedge);
			}
		}
		for(int ino = 0; ino < 8; ino++){
			for(int iedge = 0; iedge < 3; iedge++){
				if(!oct->edge[VertexEdgeMap[ino][iedge]]){
					if(elem[ino]->nodes[ino].color == -1) elem[ino]->nodes[ino].color = 1;
					elem[OctNeighbourMap[ino][iedge]]->nodes[OctNeighbourMap[ino][iedge]].color = elem[ino]->nodes[ino].color;
				}else{
					if(deb)printf("No %d, iedge %d (aresta %d)\n",ino,iedge,VertexEdgeMap[ino][iedge]);
					if(elem[ino]->nodes[ino].color == 2){

					}else{
						if(elem[OctNeighbourMap[ino][iedge]]->nodes[OctNeighbourMap[ino][iedge]].color == -1) elem[OctNeighbourMap[ino][iedge]]->nodes[OctNeighbourMap[ino][iedge]].color = 2;
					}
				}
			}
		}
		if(deb){
			for(int iel = 0; iel < 8; iel++){
				printf("Node %d color %d\n",iel, elem[iel]->nodes[iel].color);
			}
		}
	}
}

void FreeMovableNodes(hexa_tree_t* mesh){
	//mesh->oct.elem_count
	for(int iel = 0; iel < 0; iel++){
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

	bool deb = false;

	for(int ino = 0; ino < mesh->nodes.elem_count; ino++){
		octant_node_t * node = (octant_node_t *) sc_array_index(&mesh->nodes, ino);
		node->fixed = 2;
	}

	for(int ioc = 0; ioc < mesh->oct.elem_count; ioc++){
		octree_t *oct = (octree_t*) sc_array_index(&mesh->oct, ioc);

		for(int iel=0; iel<8; iel++){
			octant_t *elem = (octant_t*) sc_array_index(&mesh->elements, oct->id[iel]);

			elem->pad = oct->id[0]+1;

			//fix all the nodes
			for(int ino = 0; ino <8; ino++){
				elem->nodes[ino].fixed = 2;
				elem->nodes[ino].color = -1;
			}


			//printf("Sou o elemento %d, o numero %d do octree\n",elem->id,iel);
			//identificando as arestas interceptadas and liberando os nos a serem movidos

			//face superior
			if(iel==0 || iel ==1){
				if(elem->edge[0].ref && iel == 0){
					oct->edge[0]=true;
					elem->nodes[1].fixed = 0;
				}else if(elem->edge[0].ref && iel == 1){
					oct->edge[0]=true;
					elem->nodes[0].fixed = 0;
				}
			}
			if(iel==1 || iel ==2){
				if(elem->edge[1].ref && iel == 1){
					oct->edge[1]=true;
					elem->nodes[2].fixed = 0;
				}else if(elem->edge[1].ref  && iel == 2){
					oct->edge[1]=true;
					elem->nodes[1].fixed = 0;
				}
			}
			if(iel==2 || iel ==3){
				if(elem->edge[2].ref && iel == 2){
					oct->edge[2]=true;
					elem->nodes[3].fixed = 0;
				}else if(elem->edge[2].ref && iel == 3){
					oct->edge[2]=true;
					elem->nodes[2].fixed = 0;
				}
			}
			if(iel==3 || iel ==0){
				if(elem->edge[3].ref && iel == 3){
					oct->edge[3]=true;
					elem->nodes[0].fixed = 0;
				}else if(elem->edge[3].ref && iel == 0){
					oct->edge[3]=true;
					elem->nodes[3].fixed = 0;
				}
			}

			//face inferior
			if(iel==4 || iel ==5){
				if(elem->edge[8].ref && iel == 4){
					oct->edge[8]=true;
					elem->nodes[5].fixed = 0;
				}else if(elem->edge[8].ref && iel == 5){
					oct->edge[8]=true;
					elem->nodes[4].fixed = 0;
				}
			}
			if(iel==5 || iel ==6){
				if(elem->edge[9].ref && iel == 5){
					oct->edge[9]=true;
					elem->nodes[6].fixed = 0;
				}else if(elem->edge[9].ref && iel == 6){
					oct->edge[9]=true;
					elem->nodes[5].fixed = 0;
				}
			}
			if(iel==6 || iel ==7){
				if(elem->edge[10].ref && iel == 6){
					oct->edge[10]=true;
					elem->nodes[7].fixed = 0;
				}else if(elem->edge[10].ref && iel == 7){
					oct->edge[10]=true;
					elem->nodes[6].fixed = 0;
				}
			}
			if(iel==7 || iel ==4){
				if(elem->edge[11].ref && iel == 7){
					oct->edge[11]=true;
					elem->nodes[4].fixed = 0;
				}else if(elem->edge[11].ref && iel == 4){
					oct->edge[11]=true;
					elem->nodes[7].fixed = 0;
				}
			}

			//arestas verticais
			if(iel==0 || iel ==4){
				if(elem->edge[4].ref && iel == 0){
					oct->edge[4]=true;
					elem->nodes[4].fixed = 0;
				}else if(elem->edge[4].ref && iel == 4){
					oct->edge[4]=true;
					elem->nodes[0].fixed = 0;
				}
			}
			if(iel==1 || iel ==5){
				if(elem->edge[5].ref && iel == 1){
					oct->edge[5]=true;
					elem->nodes[5].fixed = 0;
				}else if(elem->edge[5].ref && iel == 5){
					oct->edge[5]=true;
					elem->nodes[1].fixed = 0;
				}
			}
			if(iel==2 || iel ==6){
				if(elem->edge[6].ref && iel == 2){
					oct->edge[6]=true;
					elem->nodes[6].fixed = 0;
				}else if(elem->edge[6].ref && iel == 6){
					oct->edge[6]=true;
					elem->nodes[2].fixed = 0;
				}
			}
			if(iel==3 || iel ==7){
				if(elem->edge[7].ref && iel == 3){
					oct->edge[7]=true;
					elem->nodes[7].fixed = 0;
				}else if(elem->edge[7].ref && iel == 7){
					oct->edge[7]=true;
					elem->nodes[3].fixed = 0;
				}
			}

			//printf("%s\n", oct->edge[0] ? "true" : "false");
			//identificando as faces
			//x-
			if(oct->edge[4] || oct->edge[11] || oct->edge[7] || oct->edge[3]){
				oct->face[0]=true;
				if(iel == 0){
					elem->nodes[7].fixed = 0;
				}else if(iel == 3){
					elem->nodes[4].fixed = 0;
				}else if(iel == 4){
					elem->nodes[3].fixed = 0;
				}else if(iel == 7){
					elem->nodes[0].fixed = 0;
				}
			}
			//x+
			if(oct->edge[5] || oct->edge[1] || oct->edge[6] || oct->edge[9]){
				oct->face[1]=true;
				if(iel == 1){
					elem->nodes[6].fixed = 0;
				}else if(iel == 2){
					elem->nodes[5].fixed = 0;
				}else if(iel == 5){
					elem->nodes[2].fixed = 0;
				}else if(iel == 6){
					elem->nodes[1].fixed = 0;
				}
			}
			//y-
			if(oct->edge[0] || oct->edge[5] || oct->edge[8] || oct->edge[4]){
				oct->face[2]=true;
				if(iel == 0){
					elem->nodes[5].fixed = 0;
				}else if(iel == 1){
					elem->nodes[4].fixed = 0;
				}else if(iel == 4){
					elem->nodes[1].fixed = 0;
				}else if(iel == 5){
					elem->nodes[0].fixed = 0;
				}
			}
			//y+
			if(oct->edge[2] || oct->edge[6] || oct->edge[10] || oct->edge[7]){
				oct->face[3]=true;
				if(iel == 2){
					elem->nodes[7].fixed = 0;
				}else if(iel == 3){
					elem->nodes[6].fixed = 0;
				}else if(iel == 6){
					elem->nodes[3].fixed = 0;
				}else if(iel == 7){
					elem->nodes[2].fixed = 0;
				}
			}
			//z-
			if(oct->edge[8] || oct->edge[9] || oct->edge[10] || oct->edge[11]){
				oct->face[4]=true;
				if(iel == 4){
					elem->nodes[6].fixed = 0;
				}else if(iel == 5){
					elem->nodes[7].fixed = 0;
				}else if(iel == 6){
					elem->nodes[4].fixed = 0;
				}else if(iel == 7){
					elem->nodes[5].fixed = 0;
				}
			}
			//z+
			if(oct->edge[0] || oct->edge[1] || oct->edge[2] || oct->edge[3]){
				oct->face[5]=true;
				if(iel == 0){
					elem->nodes[2].fixed = 0;
				}else if(iel == 1){
					elem->nodes[3].fixed = 0;
				}else if(iel == 2){
					elem->nodes[0].fixed = 0;
				}else if(iel == 3){
					elem->nodes[1].fixed = 0;
				}
			}

			//free the central node
			if(iel == 0) elem->nodes[6].fixed = 0;

			for(int ino = 0; ino < 8; ino++){
				octant_node_t * node = (octant_node_t *) sc_array_index(&mesh->nodes, elem->nodes[ino].id);
				if(node->fixed == 2){
					node->fixed = elem->nodes[ino].fixed;
				}
			}
		}

	}

	if(deb){
		int count =0;
		for(int ino = 0; ino < mesh->nodes.elem_count; ino++){
			octant_node_t * node = (octant_node_t *) sc_array_index(&mesh->nodes, ino);
			mesh->part_nodes[ino] = node->fixed;
			if(node->fixed != 1){
				printf("no:%d \n",node->id);
				count++;
			}
		}
		printf("achei %d nos fixos\n",count);
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
	sc_hash_array_t * hashOctree = (sc_hash_array_t *)sc_hash_array_new(sizeof (octant_t), el_hash_id, el_equal_id, &clamped);
	octant_t* r;
	bool elem_lookup;
	bool vert_lookup;

	//find the elements...
	for (int iel = 0; iel < mesh->elements.elem_count; ++iel) {
		octant_t *cut = (octant_t*) sc_array_index(&mesh->elements, iel);

		size_t position;
		octant_t key;
		key.id = cut->id;

		elem_lookup = sc_hash_array_lookup(hashOctree, &key, &position);

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

				for(int iele = 0; iele < vert->list_elem; iele++){

					int test = vert->elem[iele];

					octant_t *elem = (octant_t*) sc_array_index(&mesh->elements, test);
					key.id = elem->id;

					if(cut->id==elem->id){
						r = (octant_t*) sc_hash_array_insert_unique(hashOctree,&key,&position);
						temp_id[0] = elem->id;
						if(r!=NULL){
							r->id= elem->id;
						}else{
							r =  (octant_t*) sc_array_index(&hashOctree->a,position);
						}
						if(elem->pad==-1 || cut->pad == -1){
							temp_cut = true;
						}
					}//x+

					if(cut->nodes[1].id == elem->nodes[0].id && cut->nodes[2].id == elem->nodes[3].id &&
							cut->nodes[5].id==elem->nodes[4].id && cut->nodes[6].id==elem->nodes[7].id){
						temp_id[1] = elem->id;
						r = (octant_t*) sc_hash_array_insert_unique(hashOctree,&key,&position);
						if(r!=NULL){
							r->id= elem->id;
						}else{
							r =  (octant_t*) sc_array_index(&hashOctree->a,position);
						}
						if(elem->pad==-1 || cut->pad == -1){
							temp_cut = true;
						}
					}//x+y+

					if(cut->nodes[2].id==elem->nodes[0].id && cut->nodes[6].id==elem->nodes[4].id){
						temp_id[2]= elem->id;
						r = (octant_t*) sc_hash_array_insert_unique(hashOctree,&key,&position);
						if(r!=NULL){
							r->id= elem->id;
						}else{
							r =  (octant_t*) sc_array_index(&hashOctree->a,position);
						}
						if(elem->pad==-1 || cut->pad == -1){
							temp_cut = true;
						}
					}//y+

					if(cut->nodes[3].id==elem->nodes[0].id && cut->nodes[2].id==elem->nodes[1].id &&
							cut->nodes[6].id==elem->nodes[5].id && cut->nodes[7].id==elem->nodes[4].id){
						temp_id[3]= elem->id;
						r = (octant_t*) sc_hash_array_insert_unique(hashOctree,&key,&position);
						if(r!=NULL){
							r->id= elem->id;
						}else{
							r =  (octant_t*) sc_array_index(&hashOctree->a,position);
						}
						if(elem->pad==-1 || cut->pad == -1){
							temp_cut = true;
						}
					}//z-

					if(cut->nodes[4].id==elem->nodes[0].id && cut->nodes[5].id==elem->nodes[1].id &&
							cut->nodes[6].id==elem->nodes[2].id && cut->nodes[7].id==elem->nodes[3].id){
						temp_id[4]= elem->id;
						r = (octant_t*) sc_hash_array_insert_unique(hashOctree,&key,&position);
						if(r!=NULL){
							r->id= elem->id;
						}else{
							r =  (octant_t*) sc_array_index(&hashOctree->a,position);
						}
						if(elem->pad==-1 || cut->pad == -1){
							temp_cut = true;
						}
					}//z- x+

					if(cut->nodes[5].id==elem->nodes[0].id && cut->nodes[6].id==elem->nodes[3].id){
						temp_id[5]= elem->id;
						r = (octant_t*) sc_hash_array_insert_unique(hashOctree,&key,&position);
						if(r!=NULL){
							r->id= elem->id;
						}else{
							r =  (octant_t*) sc_array_index(&hashOctree->a,position);
						}
						if(elem->pad==-1 || cut->pad == -1){
							temp_cut = true;
						}
					}//z- x+ y+

					if(cut->nodes[6].id==elem->nodes[0].id){
						temp_id[6]= elem->id;
						r = (octant_t*) sc_hash_array_insert_unique(hashOctree,&key,&position);
						if(r!=NULL){
							r->id= elem->id;
						}else{
							r =  (octant_t*) sc_array_index(&hashOctree->a,position);
						}
						if(elem->pad==-1 || cut->pad == -1){
							temp_cut = true;
						}
					} //z- y+

					if(cut->nodes[6].id==elem->nodes[1].id && cut->nodes[7].id==elem->nodes[0].id){
						temp_id[7]= elem->id;
						r = (octant_t*) sc_hash_array_insert_unique(hashOctree,&key,&position);
						if(r!=NULL){
							r->id= elem->id;
						}else{
							r =  (octant_t*) sc_array_index(&hashOctree->a,position);
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
				//oc->mat[i] = -1;
			}

			//inicia com arestas e faces como nao cortadas
			for(int iedge = 0; iedge<12;iedge++){
				oc->edge[iedge] = false;
			}
			for(int isurf = 0; isurf<6;isurf++){
				oc->face[isurf] = false;
			}
		}
	}

	//inicializacao e fixando nos dos vertices do octree mesh->oct.elem_count
	for(int ioc = 0; ioc< 0; ioc++){
		octree_t *oc = (octree_t*) sc_array_index(&mesh->oct, ioc);
		//inicia com arestas e faces como nao cortadas
		for(int iedge = 0; iedge<12;iedge++){
			oc->edge[iedge] = false;
		}
		for(int isurf = 0; isurf<6;isurf++){
			oc->face[isurf] = false;
		}
		/*
		for(int iel = 0; iel<8; iel++){
			//printf("Octree numero:%d, elemento numero:%d\n",iel,oc->id[i]);
			if(oc->id[iel]!=-1){
				octant_t *elem = (octant_t*) sc_array_index(&mesh->elements, oc->id[iel]);
				elem->pad = oc->id[0]+1;

				//fix nodes
				if(iel==0){
					//printf("Entrei para fixar o 0 para o elemento numero:%d\n",elem->id);
					elem->nodes[0].fixed=1;
				}else if(iel==1){
					//printf("Entrei para fixar o 1 para o elemento numero:%d\n",elem->id);
					elem->nodes[1].fixed=1;
				}else if(iel==2){
					//printf("Entrei para fixar o 2 para o elemento numero:%d\n",elem->id);
					elem->nodes[2].fixed=1;
				}else if(iel==3){
					//printf("Entrei para fixar o 3 para o elemento numero:%d\n",elem->id);
					elem->nodes[3].fixed=1;
				}else if(iel==4){
					//printf("Entrei para fixar o 4 para o elemento numero:%d\n",elem->id);
					elem->nodes[4].fixed=1;
				}else if(iel==5){
					//printf("Entrei para fixar o 5 para o elemento numero:%d\n",elem->id);
					elem->nodes[5].fixed=1;
				}else if(iel==6){
					//printf("Entrei para fixar o 6 para o elemento numero:%d\n",elem->id);
					elem->nodes[6].fixed=1;
				}else if(iel==7){
					//printf("Entrei para fixar o 7 para o elemento numero:%d\n",elem->id);
					elem->nodes[7].fixed=1;
				}
			}
		}
		 */
	}


	//extract the vertex structure
	//sc_hash_array_rip (indep_vertex,  &mesh->vertex);

	//for debug
	if(false){
		for(int ioc = 0; ioc< mesh->oct.elem_count; ioc++){
			octree_t *oc = (octree_t*) sc_array_index(&mesh->oct, ioc);
			printf("Octree numero:%d, foi cortado?:%d\n",ioc,oc->cut);
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

	bool deb = false;

	time_t tstart, tend;

	tstart = time(0);
	printf("    Building the octree structure...\n");
	DoOctree(mesh);
	tend = time(0);
	//	cout << "Time in DoOctree "<< difftime(tend, tstart) <<" second(s)."<< endl;

	tstart = time(0);
	printf("    Identifying the movable nodes...\n");
	IdentifyMovableNodes(mesh);
	tend = time(0);
	//cout << "Time in IdentifyMovableNodes "<< difftime(tend, tstart) <<" second(s)."<< endl;

	if(deb){
		char fdname[80];
		sprintf(fdname,"nos_livres_%04d.txt", mesh->mpi_rank);
		FILE* treta = fopen(fdname,"w");
		for(int ino = 0; ino < mesh->nodes.elem_count; ino ++){
			octant_node_t* node = (octant_node_t*) sc_array_index (&mesh->nodes, ino);
			if(node->fixed==0){
				int nnode = 3*node->id;
				fprintf(treta,"%f %f %f\n",coords[nnode+0],coords[nnode+1],coords[nnode+2]);
			}
		}
		fclose(treta);
	}

	if(deb){
		char fdname[80];
		sprintf(fdname,"debug_edges_%04d.txt", mesh->mpi_rank);
		FILE* treta = fopen(fdname,"w");
		for(int ioc = 0; ioc < mesh->oct.elem_count; ioc ++){
			octree_t* oct = (octree_t*) sc_array_index (&mesh->oct, ioc);
			for(int iel = 0; iel < 8; iel++){
				if(oct->id[iel] !=-1){
					octant_t* elem = (octant_t*) sc_array_index (&mesh->nodes, oct->id[iel]);
					//printf("%s\n", oct->edge[0] ? "true" : "false")
					fprintf(treta,"El %d\n", oct->id[iel]);
					for(int iedge = 0; iedge < 12;iedge++){
						fprintf(treta,"%s ", elem->edge[iedge].ref ? "T" : "F");
					}
					fprintf(treta,"\n");

				}
			}

		}
		fclose(treta);
	}
	//tstart = time(0);
	//printf(" Free the movable nodes...\n");
	//FreeMovableNodes(mesh);
	//tend = time(0);
	//cout << "Time in FreeMovableNodes "<< difftime(tend, tstart) <<" second(s)."<< endl;

	tstart = time(0);
	printf("    Make the projection of the nodes into the surface...\n");
	nodes_b_mat.clear();
	ProjectFreeNodes(mesh,coords,nodes_b_mat);
	tend = time(0);
	//cout << "Time in ProjectFreeNodes "<< difftime(tend, tstart) <<" second(s)."<< endl;

	if(deb){
		int8_t* flag_nodes = (int8_t*) malloc(sizeof (int8_t) * mesh->local_n_nodes);
		memset(flag_nodes, 0, sizeof (int8_t) * mesh->local_n_nodes);

		for (int i = 0; i < mesh->local_n_nodes; i++) {
			flag_nodes[i]=0;
		}

		//verifica a criacao os nos livres e fixos na criacao dos octrees

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
		free(flag_nodes);
		for(int i = 0; i<nodes_b_mat.size();i++){
			mesh->part_nodes[nodes_b_mat[i]] = 1;
		}
	}
}
