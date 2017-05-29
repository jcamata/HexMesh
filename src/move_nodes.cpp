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

#if 0
void Move_nodes(hexa_tree_t* mesh, std::vector<double>& coords, std::vector<int>& element_ids) {
	bool clamped = true;
	GtsSegment * segments[12];
	GtsPoint * point[12];
	int Edge2GNode[12][2];
	int original_conn[8];
	FILE * fdbg;

	int id_offset = 0;

	fdbg = fopen("move_nodes.txt", "w");
	int temp_1=0;

	for (int iel = 0; iel < element_ids.size(); ++iel) {
		octant_t *elem = (octant_t*) sc_array_index(&mesh->elements, element_ids[iel]);
		elem->pad = -1;
		//printf("Element: %d\n", element_ids[iel]);

		bool edge_list[12];
		int ed_cont = 0;

		for (int edge = 0; edge < 12; ++edge) {
			point[edge] = NULL;
			int node1 = elem->nodes[EdgeVerticesMap[edge][0]].id;
			int node2 = elem->nodes[EdgeVerticesMap[edge][1]].id;

			//fprintf(fdbg, " Nodes : %d, %d\n", node1, node2);
			//printf("Element: %d, Edges: %d\n", element_ids[iel],edge);

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
					fprintf(fdbg, "El: %d\n", element_ids[iel]);
					fprintf(fdbg, "edge, %d\n ", edge);
					fprintf(fdbg, "nodes  : %d, %d\n", node1, node2);
					fprintf(fdbg, "point : %f, %f, %f\n", point[edge]->x, point[edge]->y, point[edge]->z);
					edge_list[edge] = true;
					ed_cont++;
					break;
				}
				list = list->next;
			}
		}

		if(edge_list[4] && edge_list[5] && edge_list[6] && edge_list[7]){
			temp_1++;
			for(int i_node=4; i_node<8; i_node++){
				int node1 = elem->nodes[i_node].id;
				coords[node1 * 3 + 2] = point[i_node]->z;
			}
		}



		//
		//octant_node_t *gnode1 = (octant_node_t*) sc_array_index(&mesh->nodes, node1);
		//octant_node_t *gnode2 = (octant_node_t*) sc_array_index(&mesh->nodes, node2);
		/*
			if (edge_list[edge] && gnode1->color != 1 && gnode2->color != 1) {

				//TODO SEGFAULT aqui no gts point
				GtsPoint *p0 = gts_point_new(gts_point_class(),
						coords[node1 * 3],
						coords[node1 * 3 + 1],
						coords[node1 * 3 + 2]);

				GtsPoint *p1 = gts_point_new(gts_point_class(),
						coords[node2 * 3],
						coords[node2 * 3 + 1],
						coords[node2 * 3 + 2]);


				double d_c1 = gts_point_distance(point[edge], p0);
				double d_c2 = gts_point_distance(point[edge], p1);


				fprintf(fdbg, "El: %d\n", element_ids[iel]);
				fprintf(fdbg, "nodes  : %d, %d\n", node1, node2);
				fprintf(fdbg, "point : %f, %f, %f\n", point[edge]->x, point[edge]->y, point[edge]->z);
				fprintf(fdbg, "no 0 : %f, %f, %f\n", coords[node1 * 3], coords[node1 * 3 + 1], coords[node1 * 3 + 2]);
				fprintf(fdbg, "no 1 : %f, %f, %f\n", coords[node2 * 3], coords[node2 * 3 + 1], coords[node2 * 3 + 2]);
				fprintf(fdbg, "d_c1: %f, d_c2: %f\n", d_c1, d_c2);

				if(0){
					if (edge == 4 || edge == 5 || edge == 6 || edge == 7) {

						if (coords[node2 * 3 + 2] == 0) {
							coords[node1 * 3 + 2] = point[edge]->z;
						} else if (coords[node1 * 3 + 2] == 0) {
							coords[node2 * 3 + 2] = point[edge]->z;
						} else {

							if (d_c1 >= d_c2) {
								coords[node2 * 3 + 2] = point[edge]->z;
								gnode2->pad = 1;

							} else {
								coords[node1 * 3 + 2] = point[edge]->z;
								gnode1->pad = 1;
							}

						}

					} else {
						if (d_c1 >= d_c2) {
							coords[node2 * 3] = point[edge]->x;
							coords[node2 * 3 + 1] = point[edge]->y;
							gnode2->pad = 1;

						} else {
							coords[node1 * 3] = point[edge]->x;
							coords[node1 * 3 + 1] = point[edge]->y;
							gnode1->pad = 1;
						}
					}

				}
			}
		 */

	}
	printf("templates1: %d\n",temp_1);
	fclose(fdbg);

}
#endif

// Change the node positions to fit the surface.
void Move_nodes(hexa_tree_t* mesh, const char* surface_bathy, vector<double>& coords, std::vector<int>& element_ids) {

	GtsPoint *p;
	double dz;
	double d, d0;
	double zmax;

	mesh->gdata.s = SurfaceRead(surface_bathy);

	FILE *fout = fopen("surface.dat", "w");
	gts_surface_print_stats(mesh->gdata.s, fout);
	fclose(fout);

	// Get the surface bounding box
	mesh->gdata.bbox = gts_bbox_surface(gts_bbox_class(), mesh->gdata.s);

	// Change the box size to cut the external elements
	double factor = 0.05;
	double x_factor = (mesh->gdata.bbox->x2 - mesh->gdata.bbox->x1)*factor;
	double y_factor = (mesh->gdata.bbox->y2 - mesh->gdata.bbox->y1)*factor;

	mesh->gdata.bbox->x1 += x_factor;
	mesh->gdata.bbox->y1 += y_factor;

	mesh->gdata.bbox->x2 -= x_factor;
	mesh->gdata.bbox->y2 -= y_factor;

	double Lx = (mesh->gdata.bbox->x2 - mesh->gdata.bbox->x1);
	double Ly = (mesh->gdata.bbox->y2 - mesh->gdata.bbox->y1);
	double zmin = ((Lx < Ly) ? -Lx : -Ly);

	double el_size = Lx/mesh->ncellx;

	// Build the bounding box tree
	mesh->gdata.bbt = gts_bb_tree_surface(mesh->gdata.s);

	p = gts_point_new(gts_point_class(), 0.0, 0.0, mesh->gdata.bbox->z2);

	for (int iel = 0; iel < element_ids.size(); ++iel) {
		octant_t *elem = (octant_t*) sc_array_index(&mesh->elements, element_ids[iel]);

		for(int ino=0; ino<8; ino++){
			int node = elem->nodes[ino].id;

			if(elem->nodes[ino].z==0){
				coords[node * 3 + 2] = 0;
			}else{
				p->x = coords[3*node];
				p->y = coords[3*node+1];

				//d = gts_bb_tree_point_distance(mesh->gdata.bbt, p, distance, NULL);

				p->z = 0;
				d0 = gts_bb_tree_point_distance(mesh->gdata.bbt, p, distance, NULL);

				int n_el = ceil(d0/el_size/pow(3,elem->level-1));

				dz = d0/ (double) n_el;

				int z_el = ceil(coords[3*node+2]/el_size/pow(3,elem->level-1));

				z_el=abs(z_el)+1;

				double z = - dz * z_el;

				//printf("El: %d, d0: %f, n_el: %d, z_el: %d, z: %f\n",element_ids[iel],d0,n_el,z_el,z);

				coords[node * 3 + 2] = z;

				/*
				zmax = mesh->gdata.bbox->z2 - d;
				//printf("El: %d, d: %f, zmax: %f, bbz2: %f\n",element_ids[iel],d,zmax,mesh->gdata.bbox->z2);

				dz = (zmax - zmin) / (double) mesh->ncellz;
				//double z = zmax +  dz* elem->nodes[ino].z;
				//printf("El: %d, node: %d, nz %d, z:%f\n",element_ids[iel],node,elem->nodes[ino].z,z);
				//coords[node * 3 + 2] = z;
				 *
				 */

			}
		}
	}

}

