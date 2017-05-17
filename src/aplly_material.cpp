
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
int bbox_is_stabbed_zinfty (GtsBBox * bb, GtsPoint * p){

	g_return_val_if_fail ((bb != NULL), 0);
	g_return_val_if_fail ((p != NULL), 0);

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
GSList * bb_tree_stabbed_zinfty (GNode * tree, GtsPoint * p)
{
	GSList  * list = NULL;
	GtsBBox * bb;
	GNode   * i;

	g_return_val_if_fail (tree != NULL, NULL);
	g_return_val_if_fail (p != NULL, NULL   );

	bb = (GtsBBox *) tree->data;
        
	if (!bbox_is_stabbed_zinfty (bb, p))
		return NULL;
	if (tree->children == NULL) /* leaf node */
		return g_slist_prepend (NULL, bb);
	i = tree->children;
	while (i) {
		list = g_slist_concat (list, bb_tree_stabbed_zinfty (i, p));
		i = i->next;
	}
	return list;
}

bool is_point_over_surface(GtsPoint * p, GNode    * tree){

	g_return_val_if_fail ((p    != NULL), false);
	g_return_val_if_fail ((tree != NULL), false);

	GSList* list = bb_tree_stabbed_zinfty(tree, p);
	double * orientation;
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

	return 0;

}


//Aplly the material properties to the elements
void Apply_material(hexa_tree_t *mesh, std::vector<double>& coords, std::vector<int>& element_ids, const char* surface_bathy){

	sc_array_t *elements = &mesh->elements;

	bool over;

	// Build the bounding box tree
	//mesh->gdata.s = SurfaceRead(surface_bathy);
        GtsSurface *bathymetry     = SurfaceRead(surface_bathy);
        GNode      *bbt_bathymetry = gts_bb_tree_surface(bathymetry);
        	
        GtsBBox* bbox = gts_bbox_surface(gts_bbox_class(), bathymetry);
	if (mesh->mpi_rank == 0) {
		printf("Bounding box: \n");
		printf(" x ranges from %f to %f\n", bbox->x1,bbox->x2);
		printf(" y ranges from %f to %f\n", bbox->y1,bbox->y2);
		printf(" z ranges from %f to %f\n", bbox->z1,bbox->z2);
	}
        int mat1 = 0;
        int mat2 = 0;
	//GtsPoint * p = gts_point_new(gts_point_class(),coords[n_id], coords[n_id + 1],coords[n_id + 2]);
	GtsPoint * p = gts_point_new(gts_point_class(),0,0,0);

	for(int iel = 0; iel < elements->elem_count; ++iel) 
        {
		octant_t *elem = (octant_t*) sc_array_index(&mesh->elements, iel);
                elem->n_mat = -1;
           
                // Getting baricenter
		int xyz_min = elem->nodes[0].id;
                int xyz_max = elem->nodes[7].id;
                double dx    = 0.5*(coords[xyz_max]   - coords[xyz_min]  );
                double dy    = 0.5*(coords[xyz_max+1] - coords[xyz_min+1]);
                double dz    = 0.5*(coords[xyz_max+2] - coords[xyz_min+2]);
                
		gts_point_set(p, coords[xyz_min]+dx, coords[xyz_min]+dy,coords[xyz_min]+ dz);
                

		over = is_point_over_surface( p, bbt_bathymetry);

		if(over){
			elem->n_mat = 0;
                        mat1++;
		}else{
			elem->n_mat = 1;
                        mat2++;
		}
	}
        
        printf("Apply material has found: \n MAT1: %d elements\n MAT2: %d elements\n", mat1, mat2);
        
}
