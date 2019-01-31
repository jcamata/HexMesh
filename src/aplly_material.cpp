
#include <gts.h>
#include <glib.h>

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
	sc_array_t *elements = &mesh->elements;

	int min_n_mat = 100;
	int tot_n_mat = 0;
	int cdoub = 0;
	for (int iel = 0; iel < elements->elem_count; ++iel) {
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
	for (int iel = 0; iel < elements->elem_count; ++iel) {
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
	for (int iel = 0; iel < elements->elem_count; ++iel) {
		octant_t *elem = (octant_t*) sc_array_index(&mesh->elements, iel);
		elem->n_mat=offsetdoub[elem->n_mat];
	}
}
void Apply_material(hexa_tree_t *mesh, std::vector<double>& coords, std::vector<int>& element_ids, const char* surface_bathy) {

	sc_array_t *elements = &mesh->elements;

	bool over, over1;

	// Build the bounding box tree
	//mesh->gdata.s = SurfaceRead(surface_bathy);
	GtsSurface *bathymetry = SurfaceRead(surface_bathy);
	GNode *bbt_bathymetry = gts_bb_tree_surface(bathymetry);

	GtsBBox* bbox = gts_bbox_surface(gts_bbox_class(), bathymetry);

	if (mesh->mpi_rank == 0) {
		printf("Bounding box: \n");
		printf(" x ranges from %f to %f\n", bbox->x1, bbox->x2);
		printf(" y ranges from %f to %f\n", bbox->y1, bbox->y2);
		printf(" z ranges from %f to %f\n", bbox->z1, bbox->z2);
	}
	int mat1 = 0;
	int mat2 = 0;

	//for all the mesh
	for (int iel = 0; iel < elements->elem_count; ++iel) {
		octant_t *elem = (octant_t*) sc_array_index(&mesh->elements, iel);
		elem->n_mat = -1;
		GtsPoint * point;

		//getting the baricenter of the upper surface;
		//find the centroid of the upper surface
		double cord_in_x[8], cord_in_y[8], cord_in_z[8];
		//add the nodes in the coord vector
		for (int ii = 0; ii < 8; ii++) {
			cord_in_x[ii] = coords[3 * elem->nodes[ii].id];
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
			mat2++;
		} else {
			elem->n_mat = 0;
			mat1++;
		}
		gts_object_destroy(GTS_OBJECT(point));
	}

	//now we check only the elements in the interface region
	if(true){
		for (int ioc = 0; ioc <  mesh->oct.elem_count; ++ioc) {
			octree_t * oct = (octree_t*) sc_array_index(&mesh->oct, ioc);
			GtsPoint * point[2]={NULL};
			GtsSegment * segments={0};

			for(int i = 0; i<8; i++){
				if(oct->id[i]!=-1){
					octant_t *elem = (octant_t*) sc_array_index(&mesh->elements, oct->id[i]);
					elem->n_mat = -2;

					GtsPoint * pp[5];
					//getting the baricenter of the upper surface;
					//find the centroid of the upper surface
					double cord_in_x[8], cord_in_y[8], cord_in_z[8];
					//add the nodes in the coord vector
					for (int ii = 0; ii < 8; ii++) {
						cord_in_x[ii] = coords[3 * elem->nodes[ii].id];
						cord_in_y[ii] = coords[3 * elem->nodes[ii].id + 1];
						cord_in_z[ii] = coords[3 * elem->nodes[ii].id + 2];
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

					for(int ip = 0; ip<5;ip++){
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

					if((point[0] || point[1] || point[2] || point[3]) && point[4]){
						elem->n_mat = 0;
					}else{
						elem->n_mat = 1;
					}

				}
			}
		}
	}

	gts_object_destroy(GTS_OBJECT(bathymetry));

}
