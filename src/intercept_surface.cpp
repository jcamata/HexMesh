
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

#include "hexa.h"
#include "hilbert.h"

// Read the gts file format and create a gts surface.
GtsSurface* SurfaceRead(const char* fname) {
	FILE *gts_file;
	GtsSurface *s;
	//GtsPoint *p;
	GtsFile *fp;

	gts_file = fopen(fname, "r");
	if (!gts_file) return NULL;
	fp = gts_file_new(gts_file);
	s = gts_surface_new(gts_surface_class(),
			gts_face_class(),
			gts_edge_class(),
			gts_vertex_class());
	if (gts_surface_read(s, fp)) {
		fputs("file on standard input is not a valid GTS file\n", stderr);
		fprintf(stderr, "stdin:%d:%d: %s\n", fp->line, fp->pos, fp->error);
		return NULL; /* failure */
	}

	gts_file_destroy(fp);
	return s;

}

// Compute the distance between a point and a triangle.
gdouble distance(GtsPoint *p, gpointer bounded) {
	GtsTriangle *t = (GtsTriangle*) bounded;
	return gts_point_triangle_distance(p, t);
}

// Change the node positions to fit the surface.
void GetMeshFromSurface(hexa_tree_t* mesh, const char* surface_topo, vector<double>& coords) {

	GtsPoint *p;
	double dx, dy, dz;
	double d;
	double zmax;

	sc_array_t *nodes = &mesh->nodes;

	mesh->gdata.s = SurfaceRead(surface_topo);

	FILE *fout = fopen("surface.dat", "w");
	gts_surface_print_stats(mesh->gdata.s, fout);
	fclose(fout);

	// Get the surface bounding box
	mesh->gdata.bbox = gts_bbox_surface(gts_bbox_class(), mesh->gdata.s);
	if (mesh->mpi_rank == 0) {
		printf("Bounding box: \n");
		printf(" x ranges from %f to %f\n", mesh->gdata.bbox->x1, mesh->gdata.bbox->x2);
		printf(" y ranges from %f to %f\n", mesh->gdata.bbox->y1, mesh->gdata.bbox->y2);
		printf(" z ranges from %f to %f\n", mesh->gdata.bbox->z1, mesh->gdata.bbox->z2);
	}

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

	// Get grid-spacing at x and y direction
	dx = (mesh->gdata.bbox->x2 - mesh->gdata.bbox->x1) / (double) mesh->ncellx;
	dy = (mesh->gdata.bbox->y2 - mesh->gdata.bbox->y1) / (double) mesh->ncelly;

	coords.resize(nodes->elem_count * 3);

	// Build the bounding box tree
	mesh->gdata.bbt = gts_bb_tree_surface(mesh->gdata.s);

	p = gts_point_new(gts_point_class(), 0.0, 0.0, mesh->gdata.bbox->z2);

	for (int i = 0; i < nodes->elem_count; ++i) {
		octant_node_t* n = (octant_node_t*) sc_array_index(nodes, i);
		p->x = mesh->gdata.bbox->x1 + n->x*dx;
		p->y = mesh->gdata.bbox->y1 + n->y*dy;

		d = gts_bb_tree_point_distance(mesh->gdata.bbt, p, distance, NULL);
		zmax = mesh->gdata.bbox->z2 - d;

		dz = (zmax - zmin) / (double) mesh->ncellz;
		double z = zmax - (n->z) * dz;

		coords[i * 3 + 0] = p->x;
		coords[i * 3 + 1] = p->y;
		coords[i * 3 + 2] = z;
	}
}

// Found the intercepted elements
void GetInterceptedElements(hexa_tree_t* mesh, std::vector<double>& coords, std::vector<int>& elements_ids, const char* surface_bathy) {
	sc_array_t *elements = &mesh->elements;
	GtsBBox *box;

	mesh->gdata.s = SurfaceRead(surface_bathy);
	mesh->gdata.bbox = gts_bbox_surface(gts_bbox_class(), mesh->gdata.s);
	mesh->gdata.bbt = gts_bb_tree_surface(mesh->gdata.s);

	box = gts_bbox_new(gts_bbox_class(), 0, 0, 0, 0, 1, 1, 1);

	for (int iel = 0; iel < elements->elem_count; ++iel) {

		octant_t *elem = (octant_t*) sc_array_index(&mesh->elements, iel);
		elem->pad = 0;

		box->x1 = box->y1 = box->z1 = 1.0E10;
		box->x2 = box->y2 = box->z2 = -1.0E10;

		for (int i = 0; i < 8; ++i) {
			octant_node_t* node = &elem->nodes[i];
			int id = node->id;
			double x = coords[id * 3];
			double y = coords[id * 3 + 1];
			double z = coords[id * 3 + 2];
			box->x1 = (x < box->x1) ? x : box->x1;
			box->y1 = (y < box->y1) ? y : box->y1;
			box->z1 = (z < box->z1) ? z : box->z1;
			box->x2 = (x > box->x2) ? x : box->x2;
			box->y2 = (y > box->y2) ? y : box->y2;
			box->z2 = (z > box->z2) ? z : box->z2;
		}

		if (gts_bb_tree_is_overlapping(mesh->gdata.bbt, box)) {
			elements_ids.push_back(iel);
			elem->pad = -1;
		}

		///////////
		GtsSegment * segments[12]={0};
		GtsPoint * point[12]={NULL};
		int ed_cont = 0;

		for (int edge = 0; edge < 12; ++edge) {
			point[edge] = NULL;
			int node1 = elem->nodes[EdgeVerticesMap[edge][0]].id;
			int node2 = elem->nodes[EdgeVerticesMap[edge][1]].id;
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
					elem->edge[edge].ref = true;
					elem->pad = -1;
					ed_cont++;
					break;
				}
				list = list->next;
			}
		}

		//Bounding box intercepted
		if(elem->pad == -1 && ed_cont == 0){
			elem->pad = 0;
			for (int edge = 0; edge < 12; ++edge) {
				elem->edge[edge].ref = false;
			}
		}

	}
}

//Found the intersection between a line and a triangle
GtsPoint* SegmentTriangleIntersection(GtsSegment * s, GtsTriangle * t) {
	GtsPoint * A, * B, * C, * D, * E;
	gint ABCE, ABCD, ADCE, ABDE, BCDE;
	GtsEdge * AB, * BC, * CA;
	gdouble a, b, c;

	//g_return_val_if_fail(s != NULL, NULL);
	//g_return_val_if_fail(t != NULL, NULL);
	if( s == NULL) return NULL;
	if (t == NULL) return NULL;

	gts_triangle_vertices_edges(t, NULL,
			(GtsVertex **) & A,
			(GtsVertex **) & B,
			(GtsVertex **) & C,
			&AB, &BC, &CA);
	D = GTS_POINT(s->v1);
	E = GTS_POINT(s->v2);

	ABCE = gts_point_orientation_3d_sos(A, B, C, E);
	ABCD = gts_point_orientation_3d_sos(A, B, C, D);
	if (ABCE < 0 || ABCD > 0) {
		GtsPoint * tmpp;
		gint tmp;

		tmpp = E;
		E = D;
		D = tmpp;
		tmp = ABCE;
		ABCE = ABCD;
		ABCD = tmp;
	}
	if (ABCE < 0 || ABCD > 0)
		return NULL;
	ADCE = gts_point_orientation_3d_sos(A, D, C, E);
	if (ADCE < 0)
		return NULL;
	ABDE = gts_point_orientation_3d_sos(A, B, D, E);
	if (ABDE < 0)
		return NULL;
	BCDE = gts_point_orientation_3d_sos(B, C, D, E);
	if (BCDE < 0)
		return NULL;
	a = gts_point_orientation_3d(A, B, C, E);
	b = gts_point_orientation_3d(A, B, C, D);
	if (a != b) {
		c = a / (a - b);
		return gts_point_new(gts_point_class(),
				E->x + c * (D->x - E->x),
				E->y + c * (D->y - E->y),
				E->z + c * (D->z - E->z));
	}
	/* D and E are contained within ABC */
#ifdef DEBUG
	fprintf(stderr,
			"segment: %p:%s triangle: %p:%s intersection\n"
			"D and E contained in ABC\n",
			s, GTS_NEDGE(s)->name, t, GTS_NFACE(t)->name);
#endif /* DEBUG */  
	//g_assert(a == 0.0);
	return gts_point_new(gts_point_class(),
			(E->x + D->x) / 2.,
			(E->y + D->y) / 2.,
			(E->z + D->z) / 2.);
}
