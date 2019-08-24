
#include <gts.h>
#include <glib.h>
//#include <gmodule.h>

using namespace std;
#include <sc.h>
#include <sc_io.h>
#include <sc_containers.h>

#include "hexa.h"

int bbox_is_stabbed_zinfty(GtsBBox * bb, GtsPoint * p) {

	g_return_val_if_fail((bb != NULL), 0);
	g_return_val_if_fail((p != NULL), 0);

	if (p->x < bb->x1 || p->x > bb->x2 ||
			p->y < bb->y1 || p->y > bb->y2 ||
			p->z < bb->z1)
		return 0;
	return 1;
}

GSList * bb_tree_stabbed_zinfty(GNode * tree, GtsPoint * p) {
	GSList * list = NULL;
	GtsBBox * bb;
	GNode * i;

	g_return_val_if_fail(tree != NULL, NULL);
	g_return_val_if_fail(p != NULL, NULL);

	bb = (GtsBBox *) tree->data;

	if (!bbox_is_stabbed_zinfty(bb, p))
		return NULL;
	if (tree->children == NULL) /* leaf node */
		return g_slist_prepend(NULL, bb);
	i = tree->children;
	while (i) {
		list = g_slist_concat(list, bb_tree_stabbed_zinfty(i, p));
		i = i->next;
	}
	return list;
}

bool is_point_over_surface(GtsPoint * p, GNode * tree) {

	g_return_val_if_fail((p != NULL), false);
	g_return_val_if_fail((tree != NULL), false);

	GSList* list = bb_tree_stabbed_zinfty(tree, p);
	return (list != NULL);

}

void Adjust_material(hexa_tree_t *mesh) {

	int min_n_mat = 100;
	int tot_n_mat = 0;
	int cdoub = 0;
	for (int iel = 0; iel < mesh->elements.elem_count; ++iel) {
		octant_t *elem = (octant_t*) sc_array_index(&mesh->elements, iel);
		if (elem->n_mat>tot_n_mat){
			tot_n_mat = elem->n_mat;
		}
		if (elem->n_mat<min_n_mat){
			min_n_mat = elem->n_mat;
		}
	}

	tot_n_mat=tot_n_mat+1;
	int idover[tot_n_mat];
	fill_n(idover, tot_n_mat, -1); //all elements are under the bathymetry
	for (int iel = 0; iel < mesh->elements.elem_count; ++iel) {
		octant_t *elem = (octant_t*) sc_array_index(&mesh->elements, iel);
		if (idover[elem->n_mat]==-1){
			idover[elem->n_mat]=0;
		}
	}
	int coff = 0;
	int offsetdoub[tot_n_mat];
	fill_n(offsetdoub, tot_n_mat,-1); //all elements are under the bathymetry
	for (int iel = 0; iel < tot_n_mat; ++iel) {
		if (idover[iel]==0){
			offsetdoub[iel]=coff;
			coff = coff+1;
		}
	}
	for (int iel = 0; iel < mesh->elements.elem_count; ++iel) {
		octant_t *elem = (octant_t*) sc_array_index(&mesh->elements, iel);
		elem->n_mat=offsetdoub[elem->n_mat];
	}
}

void Apply_material(hexa_tree_t *mesh, std::vector<double>& coords, std::vector<int>& element_ids, const char* surface_bathy) {

	bool over, over1;
	bool deb = false;

	GtsBBox* bbox = mesh->gdata.bbox;
	GNode *bbt_bathymetry = mesh->gdata.bbt;

	//for all the mesh
	for (int iel = 0; iel < mesh->elements.elem_count; ++iel) {
		octant_t *elem = (octant_t*) sc_array_index(&mesh->elements, iel);
		elem->n_mat = -1;
		GtsPoint * point;

		//getting the baricenter of the upper surface;
		//find the centroid of the upper surface
		double cord_in_x[8], cord_in_y[8], cord_in_z[8];
		//add the nodes in the coord vector
		for (int ii = 0; ii < 8; ii++) {
			cord_in_x[ii] = coords[3 * elem->nodes[ii].id + 0];
			cord_in_y[ii] = coords[3 * elem->nodes[ii].id + 1];
			cord_in_z[ii] = coords[3 * elem->nodes[ii].id + 2];
		}

		//superficie superior
		double cord_in_ref[3];
		cord_in_ref[0] = 0;
		cord_in_ref[1] = 0;
		cord_in_ref[2] = -1;
		point = LinearMapHex(cord_in_ref, cord_in_x, cord_in_y, cord_in_z);
		over = is_point_over_surface(point, bbt_bathymetry);

		//superficie inferior
		cord_in_ref[0] = 0;
		cord_in_ref[1] = 0;
		cord_in_ref[2] = 0.5;
		point = LinearMapHex(cord_in_ref, cord_in_x, cord_in_y, cord_in_z);
		over1 = is_point_over_surface(point, bbt_bathymetry);

		if(over && over1){
			elem->n_mat = 1;
		} else {
			elem->n_mat = 0;
		}
		gts_object_destroy(GTS_OBJECT(point));
	}

	//now we check only the elements in the interface region aka mesh->octree
	if(true){
		for (int ioc = 0; ioc <  mesh->oct.elem_count; ++ioc) {
			octree_t * oct = (octree_t*) sc_array_index(&mesh->oct, ioc);
			octant_t *elem[8];
			GtsPoint * point;

			int mat1 = 0;
			int mat2 = 0;
			int color1 = 0;
			int color2 = 0;
			int a1 = 0;
			int a2 = 0;
			int a3 = 0;
			int a4 = 0;
			for(int iel = 0; iel < 8; iel++){
				elem[iel] = (octant_t *)sc_array_index(&mesh->elements,oct->id[iel]);

				if(elem[iel]->n_mat == 0) mat1++;
				if(elem[iel]->n_mat == 1) mat2++;

				if(elem[iel]->nodes[iel].color == 1) color1++;
				if(elem[iel]->nodes[iel].color == 2) color2++;

				if(elem[iel]->n_mat == 0 && elem[iel]->nodes[iel].color == 1) a1++;
				if(elem[iel]->n_mat == 0 && elem[iel]->nodes[iel].color == 2) a2++;
				if(elem[iel]->n_mat == 1 && elem[iel]->nodes[iel].color == 1) a3++;
				if(elem[iel]->n_mat == 1 && elem[iel]->nodes[iel].color == 2) a4++;
			}

			if((mat1 == color1 || mat1 == color2) && (mat2 == color1 || mat2 == color2)){

			}else{
				//caso tipico, divisao em z...
				if((elem[0]->nodes[0].color == 1 && elem[1]->nodes[1].color == 1 &&
						elem[2]->nodes[2].color == 1 && elem[3]->nodes[3].color == 1 ) &&
						(elem[4]->nodes[4].color == 2 && elem[5]->nodes[5].color == 2 &&
								elem[6]->nodes[6].color == 2 && elem[7]->nodes[7].color == 2)){
					elem[0]->n_mat = 1;
					elem[1]->n_mat = 1;
					elem[2]->n_mat = 1;
					elem[3]->n_mat = 1;

					elem[4]->n_mat = 0;
					elem[5]->n_mat = 0;
					elem[6]->n_mat = 0;
					elem[7]->n_mat = 0;
				}else{
					for(int iel = 0; iel < 8; iel++){
						int node = elem[iel]->nodes[iel].id;
						double xx = coords[3*node+0];
						double yy = coords[3*node+1];
						double zz = coords[3*node+2];
						//if(iel<4) zz += 10;
						//if(iel>3) zz -= 10;

						/*
						gts_point_set(point,xx,yy,zz);
						over = is_point_over_surface(point, bbt_bathymetry);

						if(over){
							elem[iel]->n_mat = 1;
						}else{
							elem[iel]->n_mat = 0;
						}
						 */

						GtsPoint* p = NULL;
						GtsVertex *v1 = gts_vertex_new(gts_vertex_class(), xx, yy, zz);
						GtsVertex *v2 = gts_vertex_new(gts_vertex_class(), xx, yy, 2*bbox->z2);
						GtsSegment * segments = gts_segment_new(gts_segment_class(), v1, v2);
						GtsBBox *sb = gts_bbox_segment(gts_bbox_class(), segments);
						GSList* list = gts_bb_tree_overlap(mesh->gdata.bbt, sb);
						while (list) {
							GtsBBox *b = GTS_BBOX(list->data);
							p = SegmentTriangleIntersection(segments, GTS_TRIANGLE(b->bounded));
							if (p) {
								break;
							}
							list = list->next;
						}

						if(p){
							elem[iel]->n_mat = 0;
						}else{
							elem[iel]->n_mat = 1;
						}
						//printf("edge:%d, %s ",edge, elem->edge[edge].ref ? "T" : "F");

					}
					//printf("Octree %d El:%d, mat1:%d mat2:%d color1:%d color2:%d",ioc,oct->id[0],mat1,mat2,color1,color2);
					//printf(" a: %d %d %d %d\n",a1,a2,a3,a4);
				}

			}
		}

	}

	if(false){
		for (int ioc = 0; ioc <  mesh->oct.elem_count; ++ioc) {
			octree_t * oct = (octree_t*) sc_array_index(&mesh->oct, ioc);
			GtsPoint * point[2]={NULL};
			GtsSegment * segments={0};

			for(int ioc = 0; ioc < 8; ioc++){

				octant_t *elem = (octant_t*) sc_array_index(&mesh->elements, oct->id[ioc]);
				elem->n_mat = -1;

				GtsPoint * pp[5];
				//getting the baricenter of the upper surface;
				//find the centroid of the upper surface
				double cord_in_x[8], cord_in_y[8], cord_in_z[8];
				//add the nodes in the coord vector
				for (int ino = 0; ino < 8; ino++) {
					cord_in_x[ino] = coords[3 * elem->nodes[ino].id];
					cord_in_y[ino] = coords[3 * elem->nodes[ino].id + 1];
					cord_in_z[ino] = coords[3 * elem->nodes[ino].id + 2];
				}
				//superficie superior
				double cord_in_ref[3];
				cord_in_ref[0] = 1;
				cord_in_ref[1] = 1;
				cord_in_ref[2] = 0;
				pp[0] = LinearMapHex(cord_in_ref, cord_in_x, cord_in_y, cord_in_z);

				cord_in_ref[0] = 1;
				cord_in_ref[1] = -1;
				cord_in_ref[2] = 0;
				pp[1] = LinearMapHex(cord_in_ref, cord_in_x, cord_in_y, cord_in_z);

				cord_in_ref[0] = -1;
				cord_in_ref[1] = 1;
				cord_in_ref[2] = 0;
				pp[2] = LinearMapHex(cord_in_ref, cord_in_x, cord_in_y, cord_in_z);

				cord_in_ref[0] = -1;
				cord_in_ref[1] = -1;
				cord_in_ref[2] = 0;
				pp[3] = LinearMapHex(cord_in_ref, cord_in_x, cord_in_y, cord_in_z);

				cord_in_ref[0] = 0;
				cord_in_ref[1] = 0;
				cord_in_ref[2] = -0.9;
				pp[4] = LinearMapHex(cord_in_ref, cord_in_x, cord_in_y, cord_in_z);

				for(int ip = 0; ip < 5; ip++){
					GtsVertex *v1 = gts_vertex_new(gts_vertex_class(),pp[ip]->x ,pp[ip]->y ,pp[ip]->z);
					GtsVertex *v2 = gts_vertex_new(gts_vertex_class(),pp[ip]->x ,pp[ip]->y ,2*bbox->z2);
					segments = gts_segment_new(gts_segment_class(), v1, v2);
					point[ip] = NULL;
					GtsBBox *sb = gts_bbox_segment(gts_bbox_class(), segments);
					GSList* list = gts_bb_tree_overlap(mesh->gdata.bbt, sb);
					if (list == NULL) continue;
					while (list) {
						GtsBBox *b = GTS_BBOX(list->data);
						point[ip] = SegmentTriangleIntersection(segments, GTS_TRIANGLE(b->bounded));
						if (point[ip]) {
							break;
						}
						list = list->next;
					}
				}

				//if((point[0] || point[1] || point[2] || point[3]) && point[4]){
				if(point[0] && point[1] && point[2] && point[3] && point[4]){
					elem->n_mat = 0;
				}else{
					elem->n_mat = 1;
				}
			}
		}
	}

	for (int ino = 0; ino < mesh->nodes.elem_count; ++ino){
		mesh->part_nodes[ino] = 0;
	}
	for (int ioc = 0; ioc <  mesh->oct.elem_count; ++ioc) {
		octree_t * oct = (octree_t*) sc_array_index(&mesh->oct, ioc);
		for (int iel = 0; iel < 8; ++iel){
			octant_t * elem = (octant_t*) sc_array_index(&mesh->elements, oct->id[iel]);
			for(int ino = 0; ino < 8; ino++){
				mesh->part_nodes[elem->nodes[ino].id] += elem->nodes[ino].color;
			}
		}
	}
}
