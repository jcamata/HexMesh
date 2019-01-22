
#include <gts.h>
#include <glib.h>
#//include <vector>
//#include <iostream>
using namespace std;
//#include <set>
//#include <algorithm>
#include <sc.h>
#include <sc_io.h>
#include <sc_containers.h>

#include "hexa.h"
#include "refinement.h"

/**
 * bbox_is_stabbed:
 * @bb: a #GtsBBox.
 * @p: a #GtsPoint.
 *
 * Returns: %TRUE if the ray starting at @p and ending at (@p->x
 * @p->y, +infty) intersects with @bb, %FALSE otherwise.
 */
int bbox_is_stabbed_zinfty(GtsBBox * bb, GtsPoint * p) {

	g_return_val_if_fail((bb != NULL), 0);
	g_return_val_if_fail((p != NULL), 0);

	if (p->x < bb->x1 || p->x > bb->x2 ||
			p->y < bb->y1 || p->y > bb->y2 ||
			p->z < bb->z1)
		return 0;
	return 1;
}

/**
 * bb_tree_stabbed_zinfty:
 * @tree: a bounding box tree.
 * @p: a #GtsPoint.
 *
 * Returns: a list of bounding boxes, leaves of @tree which are
 * stabbed by the ray defined by @p (see gts_bbox_is_stabbed()).
 */
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
	//double * orientation;
	return (list != NULL);

	/*
    while(list)
    {
            GtsBBox     * b = (GtsBBox*)(list->data);
            GtsTriangle * t = (GtsTriangle*)(((GtsBBox*)(list->data))->bounded);
            // return one of the vertices of t, one of the edges of t or t if any
            // of these are stabbed by the ray starting at p (included) and
            // ending at (p->x, p->y, +infty), NULL otherwise. If the ray is contained
            // in the plane of the triangle NULL is also returned
            if(gts_triangle_is_stabbed(t,p, orientation));
            return 1;

            list = list->next;

    }
	 * */

	//return 0;

}


//Aplly the material properties to the elements

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
	//GtsPoint * p = gts_point_new(gts_point_class(),coords[n_id], coords[n_id + 1],coords[n_id + 2]);
	//GtsPoint * p = gts_point_new(gts_point_class(), 0, 0, 0);

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
		//if (over) {
			elem->n_mat = 1;
			//element_ids.push_back(iel);
			mat2++;
		} else {
			elem->n_mat = 0;
			mat1++;
		}
		gts_object_destroy(GTS_OBJECT(point));
	}

	//now we check the elements in the interface region
	if(false){
		for (int ioc = 0; ioc <  mesh->oct.elem_count; ++ioc) {
			octree_t * oct = (octree_t*) sc_array_index(&mesh->oct, ioc);
			GtsPoint * point[2]={NULL};
			GtsSegment * segments={0};

			for(int i = 0; i<8; i++){
				if((oct->cut==true) && (oct->id[i]!=-1) ){
					octant_t *elem = (octant_t*) sc_array_index(&mesh->elements, oct->id[i]);
					elem->n_mat = 30;//oct->mat[i];
					//printf("Elemento: %d \n",oct->id[i]);

					GtsPoint * pp[2];
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
					cord_in_ref[2] = -0.9;
					pp[0] = LinearMapHex(cord_in_ref, cord_in_x, cord_in_y, cord_in_z);

					cord_in_ref[0] = 0;
					cord_in_ref[1] = 0;
					cord_in_ref[2] = 0;
					pp[1] = LinearMapHex(cord_in_ref, cord_in_x, cord_in_y, cord_in_z);

					for(int ip = 0; ip<2;ip++){
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
								//printf("Edge: %d , node: %d coords: %f %f %f\n",edge,ii,point[edge]->x,point[edge]->y,point[edge]->z);
								break;
							}
							list = list->next;
						}
					}

					if(point[0] && point[1]){
						elem->n_mat = 0;
					}else{
						elem->n_mat = 1;
					}
				}
			}
		}
	}

	gts_object_destroy(GTS_OBJECT(bathymetry));

	printf("Apply material has found: \n MAT1: %d elements\n MAT2: %d elements\n", mat1, mat2);

}
