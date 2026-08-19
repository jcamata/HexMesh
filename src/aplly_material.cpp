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
	std::vector<int> idover(tot_n_mat, -1);
	for (int iel = 0; iel < mesh->elements.elem_count; ++iel) {
		octant_t *elem = (octant_t*) sc_array_index(&mesh->elements, iel);
		if (idover[elem->n_mat]==-1){
			idover[elem->n_mat]=0;
		}
	}
	int coff = 0;
	std::vector<int> offsetdoub(tot_n_mat, -1);
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

void ClassifyOctreeCorners(hexa_tree_t *mesh, const std::vector<double>& coords)
{
	double z_max_ray = 2.0 * (mesh->tdata.bbox ? mesh->tdata.bbox->z2 : 5000.0);
	if (mesh->gdata.bbox && 2.0 * mesh->gdata.bbox->z2 > z_max_ray)
		z_max_ray = 2.0 * mesh->gdata.bbox->z2;

	for (int ioc = 0; ioc < mesh->oct.elem_count; ioc++) {
		octree_t *oct = (octree_t*) sc_array_index(&mesh->oct, ioc);

		bool complete = true;
		for (int i = 0; i < 8; i++)
			if (oct->id[i] < 0) { complete = false; break; }
		if (!complete) continue;

		for (int i = 0; i < 8; i++) {
			octant_t *elem = (octant_t*) sc_array_index(&mesh->elements, oct->id[i]);

			int nid = elem->nodes[i].id;
			double xx = coords[3*nid+0];
			double yy = coords[3*nid+1];
			double zz = coords[3*nid+2];

			GtsVertex *v1 = gts_vertex_new(gts_vertex_class(), xx, yy, zz);
			GtsVertex *v2 = gts_vertex_new(gts_vertex_class(), xx, yy, z_max_ray);
			GtsSegment *seg = gts_segment_new(gts_segment_class(), v1, v2);
			GtsBBox *sb = gts_bbox_segment(gts_bbox_class(), seg);

			int cuts_above = 0;
			for (size_t k = 0; k < mesh->gdata_vec.size(); k++) {
				if (!mesh->gdata_vec[k].bbt) continue;
				GSList *list = gts_bb_tree_overlap(mesh->gdata_vec[k].bbt, sb);
				bool hit = false;
				for (GSList *l = list; l; l = l->next) {
					GtsBBox *b = GTS_BBOX(l->data);
					GtsPoint *pt = SegmentTriangleIntersection(seg, GTS_TRIANGLE(b->bounded));
					if (pt) {
						hit = true;
						gts_object_destroy(GTS_OBJECT(pt));
						break;
					}
				}
				if (list) g_slist_free(list);
				if (hit) cuts_above++;
			}

			if (mesh->gdata_vec.size() == 1) {
				elem->n_mat = (cuts_above > 0) ? 0 : 1;
			} else {
				elem->n_mat = cuts_above;
			}

			gts_object_destroy(GTS_OBJECT(sb));
			gts_object_destroy(GTS_OBJECT(seg));
		}
	}
}

void Apply_material(hexa_tree_t *mesh, std::vector<double>& coords) {

	bool deb = false;
	if (mesh->input.interfaceNumber == 0 || mesh->gdata_vec.empty()) {
		for (int iel = 0; iel < mesh->elements.elem_count; ++iel) {
			octant_t *elem = (octant_t*) sc_array_index(&mesh->elements, iel);
			elem->n_mat = 0;
		}
	} else {
		double z_max_ray = 2.0 * (mesh->tdata.bbox ? mesh->tdata.bbox->z2 : 5000.0);
		if (mesh->gdata.bbox && 2.0 * mesh->gdata.bbox->z2 > z_max_ray)
			z_max_ray = 2.0 * mesh->gdata.bbox->z2;

		for (int iel = 0; iel < mesh->elements.elem_count; ++iel) {
			octant_t *elem = (octant_t*) sc_array_index(&mesh->elements, iel);

			double cord_in_x[8], cord_in_y[8], cord_in_z[8];
			for (int ii = 0; ii < 8; ii++) {
				cord_in_x[ii] = coords[3 * elem->nodes[ii].id + 0];
				cord_in_y[ii] = coords[3 * elem->nodes[ii].id + 1];
				cord_in_z[ii] = coords[3 * elem->nodes[ii].id + 2];
			}

			double cord_in_ref[3] = {0.0, 0.0, -1.0};
			GtsPoint *point = LinearMapHex(cord_in_ref, cord_in_x, cord_in_y, cord_in_z);

			GtsVertex *v1 = gts_vertex_new(gts_vertex_class(), point->x, point->y, point->z);
			GtsVertex *v2 = gts_vertex_new(gts_vertex_class(), point->x, point->y, z_max_ray);

			GtsSegment *segments = gts_segment_new(gts_segment_class(), v1, v2);
			GtsBBox *sb = gts_bbox_segment(gts_bbox_class(), segments);

			int cuts_above = 0;
			for (size_t k = 0; k < mesh->gdata_vec.size(); k++) {
				if (!mesh->gdata_vec[k].bbt) continue;
				GSList *list = gts_bb_tree_overlap(mesh->gdata_vec[k].bbt, sb);
				bool hit = false;
				for (GSList *l = list; l; l = l->next) {
					GtsBBox *b = GTS_BBOX(l->data);
					GtsPoint *pt = SegmentTriangleIntersection(segments, GTS_TRIANGLE(b->bounded));
					if (pt) {
						hit = true;
						gts_object_destroy(GTS_OBJECT(pt));
						break;
					}
				}
				if (list) g_slist_free(list);
				if (hit) cuts_above++;
			}

			if (mesh->gdata_vec.size() == 1) {
				elem->n_mat = (cuts_above > 0) ? 0 : 1;
			} else {
				elem->n_mat = cuts_above;
				if (mesh->input.nmat > 0 && elem->n_mat >= mesh->input.nmat) {
					elem->n_mat = mesh->input.nmat - 1;
				}
			}

			gts_object_destroy(GTS_OBJECT(sb));
			gts_object_destroy(GTS_OBJECT(segments));
			gts_object_destroy(GTS_OBJECT(point));
		}

		// The octree neighborhood map is only required by node-moving workflows.
		// When node moving is disabled, this pass may contain incomplete connectivity.
		if (mesh->input.movingNodes == 0 || mesh->oct.elem_count == 0) {
			return;
		}
		//now we check only the elements in the interface region aka mesh->octree mesh->oct.elem_count
		for (int ioc = 0; ioc < mesh->oct.elem_count; ++ioc)
		{
			octree_t * oct = (octree_t*) sc_array_index(&mesh->oct, ioc);
			octant_t *elem[8];
			bool valid_oct = true;

			for(int iel = 0; iel < 8; iel++){
				if (oct->id[iel] < 0 || oct->id[iel] >= mesh->elements.elem_count) {
					valid_oct = false;
					break;
				}
				elem[iel] = (octant_t *)sc_array_index(&mesh->elements,oct->id[iel]);
			}
			if (!valid_oct) continue;

			// Positional check: n_mat and color must agree on WHICH corner is which side, not
			// just how many of each (an aggregate-count match doesn't rule out a scrambled
			// element-to-color correspondence). color's 1/2 labels are arbitrary relative to
			// n_mat's 0/1 (GetOctreeBipartition always starts its BFS labeling at corner 0, see
			// moving_nodes.cpp), so try both label mappings before declaring a mismatch.
			// Elements with n_mat outside {0,1} (multi-material octree corners) never match this
			// binary color scheme and always fall through to the per-corner recompute below.
			bool consistent = false;
			for (int swap = 0; swap < 2 && !consistent; swap++) {
				bool ok = true;
				for (int iel = 0; iel < 8; iel++) {
					if (elem[iel]->n_mat != 0 && elem[iel]->n_mat != 1) { ok = false; break; }
					int expect_color = (elem[iel]->n_mat == 0) == (swap == 0) ? 1 : 2;
					if (elem[iel]->nodes[iel].color != expect_color) { ok = false; break; }
				}
				consistent = ok;
			}

			if(consistent){

			}else{
				// Previously: a "clean top/bottom split" color pattern short-circuited straight to a
				// hardcoded n_mat assignment (corners 0-3 -> 1, corners 4-7 -> 0), skipping the
				// geometric ray-cast below entirely. That assumed a fixed top/bottom <-> corner-index
				// correspondence, but (a) z orientation is inverted in this mesh (see h5 node reorder
				// {4,5,6,7,0,1,2,3} in the output writer) so "corners 0-3" isn't reliably one physical
				// side, and (b) `.color` comes from ClassifyOctreeCorners, which runs on PRE-warp
				// coordinates (moving_nodes.cpp's MovingNodes calls it right after DoOctree, before
				// WarpLatticeToCoastline/ProjectFreeNodes move nodes toward the real coastline) --
				// exactly the octants near the real, non-flat interface can look like a clean split in
				// that stale color topology while the POST-warp geometry (what n_mat should reflect)
				// no longer agrees. Always verify against current geometry instead of trusting the
				// color pattern's shape.
				{
					for(int iel = 0; iel < 8; iel++){
						int node = elem[iel]->nodes[iel].id;
						if (node < 0 || (3 * node + 2) >= (int) coords.size()) {
							continue;
						}
						double xx = coords[3*node+0];
						double yy = coords[3*node+1];
						double zz = coords[3*node+2];
						GtsVertex *v1 = gts_vertex_new(gts_vertex_class(), xx, yy, zz);
						GtsVertex *v2 = gts_vertex_new(gts_vertex_class(), xx, yy, z_max_ray);
						GtsSegment * segments = gts_segment_new(gts_segment_class(), v1, v2);
						GtsBBox *sb = gts_bbox_segment(gts_bbox_class(), segments);

						int cuts_above = 0;
						for (size_t k = 0; k < mesh->gdata_vec.size(); k++) {
							if (!mesh->gdata_vec[k].bbt) continue;
							GSList* list = gts_bb_tree_overlap(mesh->gdata_vec[k].bbt, sb);
							bool hit = false;
							for (GSList *l = list; l; l = l->next) {
								GtsBBox *b = GTS_BBOX(l->data);
								GtsPoint *pt = SegmentTriangleIntersection(segments, GTS_TRIANGLE(b->bounded));
								if (pt) {
									hit = true;
									gts_object_destroy(GTS_OBJECT(pt));
									break;
								}
							}
							if (list) g_slist_free(list);
							if (hit) cuts_above++;
						}
						if (mesh->gdata_vec.size() == 1) {
							elem[iel]->n_mat = (cuts_above > 0) ? 0 : 1;
						} else {
							elem[iel]->n_mat = cuts_above;
							if (mesh->input.nmat > 0 && elem[iel]->n_mat >= mesh->input.nmat) {
								elem[iel]->n_mat = mesh->input.nmat - 1;
							}
						}

						gts_object_destroy(GTS_OBJECT(sb));
						gts_object_destroy(GTS_OBJECT(segments));
					}
				}
			}
		}

		if(deb){
			for (int ino = 0; ino < mesh->nodes.elem_count; ++ino){
				mesh->part_nodes[ino] = 0;
			}
			for (int ioc = 0; ioc < mesh->oct.elem_count; ++ioc) {
				octree_t * oct = (octree_t*) sc_array_index(&mesh->oct, ioc);
				for (int iel = 0; iel < 8; ++iel){
					octant_t * elem = (octant_t*) sc_array_index(&mesh->elements, oct->id[iel]);
					for(int ino = 0; ino < 8; ino++){
						mesh->part_nodes[elem->nodes[ino].id] += elem->nodes[ino].color;
					}
				}
			}
		}
	}
}
