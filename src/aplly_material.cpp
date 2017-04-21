
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

//Raytrace to found one position point
bool Point_is_under_surface (GtsPoint * p, GNode    * tree)  {

	bool is_under_surface = true;
	GtsBBox   * bb;
	GtsSegment *s;

	//g_return_val_if_fail ((p != NULL), false);
	if(p==NULL) exit(1);
	//g_return_val_if_fail ((tree != NULL), false);
	if(tree==NULL) exit(2);

	bb = (GtsBBox*) tree->data;
	double d;
	double dx = (bb->x2 - bb->x1);
	d = dx;
	double dy = (bb->y2 - bb->y1);
	d = (dy > d)? dy: d;
	double dz = (bb->z2 - bb->z1);
	d = (dz > d)? dz: d;

	GtsVertex *v1 = gts_vertex_new(gts_vertex_class(), p->x, p->y, p->z);
	GtsVertex *v2 = gts_vertex_new(gts_vertex_class(), p->x, p->y, bb->z1 - 1.1*d);

	s = gts_segment_new(gts_segment_class(), v1, v2);
	GtsBBox *sb = gts_bbox_segment(gts_bbox_class(), s);
	GSList* list = gts_bb_tree_overlap(tree, sb);
	while (list) {
		GtsTriangle * t = (GtsTriangle*)(((GtsBBox*)(list->data))->bounded);
		if(SegmentTriangleIntersection(s,t) != NULL) {
			is_under_surface = false;
			break;
		}
		list = list->next;
	}
	//g_slist_free (list);

	return is_under_surface;
}

//Aplly the material properties to the elements
void Material_apply(hexa_tree_t *mesh, std::vector<double>& coords, std::vector<int>& element_ids, const char* surface_bathy){

	sc_array_t *elements = &mesh->elements;

	/*
	bool under;
	GtsPoint * p;



	// Build the bounding box tree
	mesh->gdata.s = SurfaceRead(surface_bathy);
	mesh->gdata.bbt = gts_bb_tree_surface(mesh->gdata.s);
	mesh->gdata.bbox = gts_bbox_surface(gts_bbox_class(), mesh->gdata.s);

	//GNode*      gts_bb_tree_surface(mesh->gdata.s)

	for(int iel = 0; iel < elements->elem_count; ++iel) {
		octant_t *elem = (octant_t*) sc_array_index(&mesh->elements, iel);

		int n_id = elem->nodes[4].id;
		p = gts_point_set(p, coords[n_id], coords[n_id + 1],coords[n_id + 2]);
		under = Point_is_under_surface( p, gts_bb_tree_surface(mesh->gdata.s));

		if(under){
			elem->n_mat = 0;
		}else{
			elem->n_mat = 100;
		}
	}
	 */

	/*
	 bool Point_is_under_surface (GtsPoint * p, GNode    * tree)  {

	bool is_under_surface = true;
	GtsBBox   * bb;
	GtsSegment *s;

	//g_return_val_if_fail ((p != NULL), false);
	if(p==NULL) exit(1);
	//g_return_val_if_fail ((tree != NULL), false);
	if(tree==NULL) exit(2);

	bb = (GtsBBox*) tree->data;
	double d;
	double dx = (bb->x2 - bb->x1);
	d = dx;
	double dy = (bb->y2 - bb->y1);
	d = (dy > d)? dy: d;
	double dz = (bb->z2 - bb->z1);
	d = (dz > d)? dz: d;

	GtsVertex *v1 = gts_vertex_new(gts_vertex_class(), p->x, p->y, p->z);
	GtsVertex *v2 = gts_vertex_new(gts_vertex_class(), p->x, p->y, bb->z1 - 1.1*d);

	s = gts_segment_new(gts_segment_class(), v1, v2);
	GtsBBox *sb = gts_bbox_segment(gts_bbox_class(), s);
	GSList* list = gts_bb_tree_overlap(tree, sb);
	while (list) {
		GtsTriangle * t = (GtsTriangle*)(((GtsBBox*)(list->data))->bounded);
		if(SegmentTriangleIntersection(s,t) != NULL) {
			is_under_surface = false;
			break;
		}
		list = list->next;
	}
	//g_slist_free (list);

	return is_under_surface;

	 */


	for (int iel = 0; iel < element_ids.size(); ++iel) {
		octant_t *h = (octant_t*) sc_array_index(&mesh->elements, element_ids[iel]);

		if(h->z==0){
			h->pad = 3;

			if(h->pml_id==1){
				h->n_mat = 19;
			}else if(h->pml_id==2){
				h->n_mat = 20;
			}else if(h->pml_id==3){
				h->n_mat = 21;
			}else if(h->pml_id==4){
				h->n_mat = 22;
			}else if(h->pml_id==5){
				h->n_mat = 23;
			}else if(h->pml_id==6){
				h->n_mat = 24;
			}else if(h->pml_id==7){
				h->n_mat = 25;
			}else if(h->pml_id==8){
				h->n_mat = 26;
			}else if(h->pml_id==9){
				h->n_mat = 27;
			}else if(h->pml_id==10){
				h->n_mat = 28;
			}else if(h->pml_id==11){
				h->n_mat = 29;
			}else if(h->pml_id==12){
				h->n_mat = 30;
			}else if(h->pml_id==13){
				h->n_mat = 31;
			}else if(h->pml_id==14){
				h->n_mat = 32;
			}else if(h->pml_id==15){
				h->n_mat = 33;
			}else if(h->pml_id==16){
				h->n_mat = 34;
			}else if(h->pml_id==17){
				h->n_mat = 35;
			}else if(h->pml_id==18){
				h->n_mat = 36;
			}else if(h->pml_id==19){
				h->n_mat = 37;
			}else if(h->pml_id==20){
				h->n_mat = 38;
			}else if(h->pml_id==21){
				h->n_mat = 39;
			}else if(h->pml_id==22){
				h->n_mat = 40;
			}else if(h->pml_id==23){
				h->n_mat = 41;
			}else if(h->pml_id==24){
				h->n_mat = 42;
			}else if(h->pml_id==25){
				h->n_mat = 43;
			}else if(h->pml_id==26){
				h->n_mat = 44;
			}else{
				//with PML
				h->n_mat = 19;
				//without PML
				//h->n_mat = 2;
			}

		}else{
			h->pad=3;

			for(int iiel = 0; iiel < elements->elem_count; ++iiel) {
				octant_t *elem = (octant_t*) sc_array_index(&mesh->elements, iiel);

				if(h->nodes->z>=elem->nodes->z && elem->x==h->x && elem->y==h->y){
					elem->pad = 3;

					if(elem->pml_id==1){
						elem->n_mat = 19;
					}else if(elem->pml_id==2){
						elem->n_mat = 20;
					}else if(elem->pml_id==3){
						elem->n_mat = 21;
					}else if(elem->pml_id==4){
						elem->n_mat = 22;
					}else if(elem->pml_id==5){
						elem->n_mat = 23;
					}else if(elem->pml_id==6){
						elem->n_mat = 24;
					}else if(elem->pml_id==7){
						elem->n_mat = 25;
					}else if(elem->pml_id==8){
						elem->n_mat = 26;
					}else if(elem->pml_id==9){
						elem->n_mat = 27;
					}else if(elem->pml_id==10){
						elem->n_mat = 28;
					}else if(elem->pml_id==11){
						elem->n_mat = 29;
					}else if(elem->pml_id==12){
						elem->n_mat = 30;
					}else if(elem->pml_id==13){
						elem->n_mat = 31;
					}else if(elem->pml_id==14){
						elem->n_mat = 32;
					}else if(elem->pml_id==15){
						elem->n_mat = 33;
					}else if(elem->pml_id==16){
						elem->n_mat = 34;
					}else if(elem->pml_id==17){
						elem->n_mat = 35;
					}else if(elem->pml_id==18){
						elem->n_mat = 36;
					}else if(elem->pml_id==19){
						elem->n_mat = 37;
					}else if(elem->pml_id==20){
						elem->n_mat = 38;
					}else if(elem->pml_id==21){
						elem->n_mat = 39;
					}else if(elem->pml_id==22){
						elem->n_mat = 40;
					}else if(elem->pml_id==23){
						elem->n_mat = 41;
					}else if(elem->pml_id==24){
						elem->n_mat = 42;
					}else if(elem->pml_id==25){
						elem->n_mat = 43;
					}else if(elem->pml_id==26){
						elem->n_mat = 44;
					}else{
						//with PML
						elem->n_mat = 19;
						//without PML
						//elem->n_mat = 2;
					}
				}
			}
		}
	}
}
