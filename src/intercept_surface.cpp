
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

//static unsigned nDims = 2;
//static unsigned nBits = 32;

typedef struct {
	bitmask_t coord[2];
	int node_id;
} node_in_edge_t;

unsigned edge_hash_fn(const void *v, const void *u) {
	const node_in_edge_t *q = (const node_in_edge_t*) v;
	uint32_t a, b, c;

	a = (uint32_t) q->coord[0];
	b = (uint32_t) q->coord[1];
	c = (uint32_t) 0;
	sc_hash_mix(a, b, c);
	sc_hash_final(a, b, c);
	return (unsigned) c;
}

int edge_equal_fn(const void *v, const void *u, const void *w) {
	const node_in_edge_t *e1 = (const node_in_edge_t*) v;
	const node_in_edge_t *e2 = (const node_in_edge_t*) u;

	return (unsigned) ((e1->coord[0] == e2->coord[0]) && (e2->coord[1] == e2->coord[1]));

}

void shift_n_position_from_index(sc_array_t *elements, uint32_t n, uint32_t index);

int EdgeVerticesMap[12][2] = {
		{0, 1}, // Edge 0
		{1, 2}, // Edge 1
		{2, 3}, //      2
		{3, 0},
		{0, 4},
		{1, 5},
		{2, 6},
		{3, 7},
		{4, 5},
		{5, 6},
		{6, 7},
		{7, 4}
};

int EdgeVerticesMap_surf_diagonal[12][2] = {
		{0, 5},
		{1, 4},
		{3, 6},
		{2, 7},
		{0, 7},
		{3, 4},
		{1, 6},
		{2, 5},
		{4, 6},
		{5, 7},
		{0, 2},
		{1, 3}
};

int EdgeVerticesMap_vol_diagonal[4][2] = {
		{0, 6},
		{1, 7},
		{2, 4},
		{3, 5}
};

int FaceEdgesMap[6][4] = {
		{4, 11, 7, 3},
		{5, 1, 6, 9},
		{0, 5, 8, 4},
		{2, 6, 10, 7},
		{8, 9, 10, 11},
		{0, 1, 2, 3}
};

int FaceNodesMap[6][4] = {
		{0,4,7,3},
		{1,5,6,2},
		{0,4,5,1},
		{2,6,7,3},
		{4,5,6,7},
		{0,1,2,3}
};

int EdgeEdgeMap[12][4] = {
		{3, 4, 1, 5}, // Edge 0
		{0, 5, 2, 6}, // Edge 1
		{1, 6, 3, 7}, //      2
		{2, 7, 0, 4,}, //3
		{0, 3, 8, 11}, //4
		{1, 0, 9, 8}, //5
		{2, 1, 10, 9}, //6
		{3, 2, 11, 10}, //7
		{11, 4, 9, 5}, //8
		{8, 5, 10, 6}, //9
		{9, 6, 11, 7}, //10
		{10, 7, 8, 9}//11
};

//unsigned long int edge_hash_function(int a, int b){
int edge_hash_function(int c, int d){
	//unsigned long int c;
	int e;
	int a;
	int b;

	if(c>d){
		a = c;
		b = d;
	}else{
		a=d;
		b=c;
	}

	if(a>=b){
		return e = a*a+a+b;
	}else{
		return e = a+b*b;
	}
}

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

//int AddPointOnEdge(int* nodes, sc_hash_array_t* hash, int &npoints, GtsPoint *p, std::vector<double> &coords);

void GetMeshFromSurface(hexa_tree_t* tree, const char* surface_topo, vector<double>& coords) {

	//GtsSurface *s;
	//GtsBBox *box;
	//GNode *t;
	GtsPoint *p;
	//int nx, ny, nz;
	double dx, dy, dz;
	double d;
	double zmax;

	sc_array_t *nodes = &tree->nodes;
	//sc_array_t *elements = &tree->elements;

	// Note that here we use a gts file.
	// There is a tool called stl2gts that convert STL files to GTS.
	// It is installed together with the gts library.
	tree->gdata.s = SurfaceRead(surface_topo);

	FILE *fout = fopen("surface.dat", "w");
	gts_surface_print_stats(tree->gdata.s, fout);
	fclose(fout);

	// Get the surface bounding box
	tree->gdata.bbox = gts_bbox_surface(gts_bbox_class(), tree->gdata.s);
	if (tree->mpi_rank == 0) {
		printf("Bounding box: \n");
		printf(" x ranges from %f to %f\n", tree->gdata.bbox->x1, tree->gdata.bbox->x2);
		printf(" y ranges from %f to %f\n", tree->gdata.bbox->y1, tree->gdata.bbox->y2);
		printf(" z ranges from %f to %f\n", tree->gdata.bbox->z1, tree->gdata.bbox->z2);
	}

	double factor = 0.02;
	double x_factor = (tree->gdata.bbox->x2 - tree->gdata.bbox->x1)*factor;
	double y_factor = (tree->gdata.bbox->y2 - tree->gdata.bbox->y1)*factor;

	tree->gdata.bbox->x1 += x_factor;
	tree->gdata.bbox->y1 += y_factor;

	tree->gdata.bbox->x2 -= x_factor;
	tree->gdata.bbox->y2 -= y_factor;

	double Lx = (tree->gdata.bbox->x2 - tree->gdata.bbox->x1);
	double Ly = (tree->gdata.bbox->y2 - tree->gdata.bbox->y1);
	double zmin = ((Lx < Ly) ? -Lx : -Ly);

	// Get grid-spacing at x and y direction
	dx = (tree->gdata.bbox->x2 - tree->gdata.bbox->x1) / (double) tree->ncellx;
	dy = (tree->gdata.bbox->y2 - tree->gdata.bbox->y1) / (double) tree->ncelly;

	coords.resize(nodes->elem_count * 3);

	// Build the bounding box tree
	tree->gdata.bbt = gts_bb_tree_surface(tree->gdata.s);

	p = gts_point_new(gts_point_class(), 0.0, 0.0, tree->gdata.bbox->z2);

	for (int i = 0; i < nodes->elem_count; ++i) {
		octant_node_t* n = (octant_node_t*) sc_array_index(nodes, i);
		p->x = tree->gdata.bbox->x1 + n->x*dx;
		p->y = tree->gdata.bbox->y1 + n->y*dy;

		d = gts_bb_tree_point_distance(tree->gdata.bbt, p, distance, NULL);
		zmax = tree->gdata.bbox->z2 - d;

		dz = (zmax - zmin) / (double) tree->ncellz;
		double z = zmax - (n->z) * dz;

		coords[i * 3 + 0] = p->x;
		coords[i * 3 + 1] = p->y;
		coords[i * 3 + 2] = z;
	}
}

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
		elem->tem = 0;

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

		//elem->n_mat = 1;
		if (gts_bb_tree_is_overlapping(mesh->gdata.bbt, box)) {
			elements_ids.push_back(iel);
			elem->pad = -1;
		}
	}
}

/*
 1. Para cada octante, 

 1.a vefificar se duas faces paralelas nao sao interceptadas pela
 superficie.
   => Aplica-se o template 1: divide elemento ao meio gerando dois novos elementos

 1.b: quatro faces inteceptadas e as restantes nao sao paralelas
   => Aplica-se o template 2:
           - Adicionar ponto na intersecao seguindo a projecao diagonal do ponto da esquina do elemento, 
           -  Tres novos elementos são criados. 
            Obs.: Possibilidade de formacao de triangulos.  Solucao: mover pontos.

 1.c: tres faces interceptadas: Veja template3 hexMesh/doc.
 */

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

void IdentifyTemplate(hexa_tree_t* mesh, const std::vector<double>& coords, std::vector<int>& elements_ids){

	int el_0 = 0;
	int el_1 = 0;
	int el_2 = 0;
	int el_3 = 0;
	int el_4 = 0;
	int el_5 = 0;
	int el_6 = 0;
	int el_7 = 0;
	int el_8 = 0;
	int el_9 = 0;
	int el_10 = 0;
	int el_11 = 0;

	int ed_cont;

	for (int iel = 0; iel < elements_ids.size(); ++iel) {

		octant_t *elem = (octant_t*) sc_array_index(&mesh->elements, elements_ids[iel]);
		ed_cont = 0;

		for (int edge = 0; edge < 12; ++edge) {
			if(elem->edge_ref[edge]){
				ed_cont++;
			}
		}

		//template 11
		elem->pad = 140;
		elem->tem = 11;

		// element verification
		// template 1
		if(ed_cont == 1){
			if ( (elem->edge_ref[0]) && (!elem->edge_ref[1]) && (!elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (!elem->edge_ref[5]) && (!elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (!elem->edge_ref[9]) && (!elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 10;
				elem->tem = 1;

			}
			else if ( (!elem->edge_ref[0]) && (elem->edge_ref[1]) && (!elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (!elem->edge_ref[5]) && (!elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (!elem->edge_ref[9]) && (!elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 11;
				elem->tem = 1;

			}
			else if ( (!elem->edge_ref[0]) && (!elem->edge_ref[1]) && (elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (!elem->edge_ref[5]) && (!elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (!elem->edge_ref[9]) && (!elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 12;
				elem->tem = 1;

			}
			else if ( (!elem->edge_ref[0]) && (!elem->edge_ref[1]) && (!elem->edge_ref[2]) && (elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (!elem->edge_ref[5]) && (!elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (!elem->edge_ref[9]) && (!elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 13;
				elem->tem = 1;

			}
			else if ( (!elem->edge_ref[0]) && (!elem->edge_ref[1]) && (!elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(elem->edge_ref[4]) && (!elem->edge_ref[5]) && (!elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (!elem->edge_ref[9]) && (!elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 14;
				elem->tem = 1;

			}
			else if ( (!elem->edge_ref[0]) && (!elem->edge_ref[1]) && (!elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (elem->edge_ref[5]) && (!elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (!elem->edge_ref[9]) && (!elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 15;
				elem->tem = 1;

			}
			else if ( (!elem->edge_ref[0]) && (!elem->edge_ref[1]) && (!elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (!elem->edge_ref[5]) && (elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (!elem->edge_ref[9]) && (!elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 16;
				elem->tem = 1;

			}
			else if ( (!elem->edge_ref[0]) && (!elem->edge_ref[1]) && (!elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (!elem->edge_ref[5]) && (!elem->edge_ref[6]) && (elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (!elem->edge_ref[9]) && (!elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 17;
				elem->tem = 1;

			}
			else if ( (!elem->edge_ref[0]) && (!elem->edge_ref[1]) && (!elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (!elem->edge_ref[5]) && (!elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(elem->edge_ref[8]) && (!elem->edge_ref[9]) && (!elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 18;
				elem->tem = 1;

			}
			else if ( (!elem->edge_ref[0]) && (!elem->edge_ref[1]) && (!elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (!elem->edge_ref[5]) && (!elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (elem->edge_ref[9]) && (!elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 19;
				elem->tem = 1;

			}
			else if ( (!elem->edge_ref[0]) && (!elem->edge_ref[1]) && (!elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (!elem->edge_ref[5]) && (!elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (!elem->edge_ref[9]) && (elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 20;
				elem->tem = 1;

			}
			else if ( (!elem->edge_ref[0]) && (!elem->edge_ref[1]) && (!elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (!elem->edge_ref[5]) && (!elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (!elem->edge_ref[9]) && (!elem->edge_ref[10]) && (elem->edge_ref[11])) {

				elem->pad = 21;
				elem->tem = 1;

			}
		}
		//template 2
		if(ed_cont == 4){
			if ( (elem->edge_ref[0]) && (!elem->edge_ref[1]) && (elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (!elem->edge_ref[5]) && (!elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(elem->edge_ref[8]) && (!elem->edge_ref[9]) && (elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 22;
				elem->tem = 2;

			}
			else if ( (!elem->edge_ref[0]) && (elem->edge_ref[1]) && (!elem->edge_ref[2]) && (elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (!elem->edge_ref[5]) && (!elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (elem->edge_ref[9]) && (!elem->edge_ref[10]) && (elem->edge_ref[11])) {

				elem->pad = 23;
				elem->tem = 2;

			}
			else if ( (!elem->edge_ref[0]) && (!elem->edge_ref[1]) && (!elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(elem->edge_ref[4]) && (elem->edge_ref[5]) && (elem->edge_ref[6]) && (elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (!elem->edge_ref[9]) && (!elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 24;
				elem->tem = 2;

			}
		}
		//template 3
		if(ed_cont == 2){
			if ( (elem->edge_ref[0]) && (!elem->edge_ref[1]) && (elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (!elem->edge_ref[5]) && (!elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (!elem->edge_ref[9]) && (!elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 25;
				elem->tem = 3;

			}
			else if ( (!elem->edge_ref[0]) && (elem->edge_ref[1]) && (!elem->edge_ref[2]) && (elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (!elem->edge_ref[5]) && (!elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (!elem->edge_ref[9]) && (!elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 26;
				elem->tem = 3;

			}
			else if ( (!elem->edge_ref[0]) && (!elem->edge_ref[1]) && (!elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (!elem->edge_ref[5]) && (!elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(elem->edge_ref[8]) && (!elem->edge_ref[9]) && (elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 27;
				elem->tem = 3;

			}
			else if ( (!elem->edge_ref[0]) && (!elem->edge_ref[1]) && (!elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (!elem->edge_ref[5]) && (!elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (elem->edge_ref[9]) && (!elem->edge_ref[10]) && (elem->edge_ref[11])) {

				elem->pad = 28;
				elem->tem = 3;

			}
			else if ( (!elem->edge_ref[0]) && (!elem->edge_ref[1]) && (!elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(elem->edge_ref[4]) && (elem->edge_ref[5]) && (!elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (!elem->edge_ref[9]) && (!elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 29;
				elem->tem = 3;

			}
			else if ( (!elem->edge_ref[0]) && (!elem->edge_ref[1]) && (!elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (elem->edge_ref[5]) && (elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (!elem->edge_ref[9]) && (!elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 30;
				elem->tem = 3;

			}
			else if ( (!elem->edge_ref[0]) && (!elem->edge_ref[1]) && (!elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (!elem->edge_ref[5]) && (elem->edge_ref[6]) && (elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (!elem->edge_ref[9]) && (!elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 31;
				elem->tem = 3;

			}
			else if ( (!elem->edge_ref[0]) && (!elem->edge_ref[1]) && (!elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(elem->edge_ref[4]) && (!elem->edge_ref[5]) && (!elem->edge_ref[6]) && (elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (!elem->edge_ref[9]) && (!elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 32;
				elem->tem = 3;

			}
			else if ( (elem->edge_ref[0]) && (!elem->edge_ref[1]) && (!elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (!elem->edge_ref[5]) && (!elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(elem->edge_ref[8]) && (!elem->edge_ref[9]) && (!elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 33;
				elem->tem = 3;

			}
			else if ( (!elem->edge_ref[0]) && (elem->edge_ref[1]) && (!elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (!elem->edge_ref[5]) && (!elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (elem->edge_ref[9]) && (!elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 34;
				elem->tem = 3;

			}
			else if ( (!elem->edge_ref[0]) && (!elem->edge_ref[1]) && (elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (!elem->edge_ref[5]) && (!elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (!elem->edge_ref[9]) && (elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 35;
				elem->tem = 3;

			}
			else if ( (!elem->edge_ref[0]) && (!elem->edge_ref[1]) && (!elem->edge_ref[2]) && (elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (!elem->edge_ref[5]) && (!elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (!elem->edge_ref[9]) && (!elem->edge_ref[10]) && (elem->edge_ref[11])) {

				elem->pad = 36;
				elem->tem = 3;

			}
		}
		//template 4
		if(ed_cont == 4){
			if ( (elem->edge_ref[0]) && (elem->edge_ref[1]) && (!elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (!elem->edge_ref[5]) && (!elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(elem->edge_ref[8]) && (elem->edge_ref[9]) && (!elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 37;
				elem->tem = 4;

			}
			else if ( (!elem->edge_ref[0]) && (elem->edge_ref[1]) && (elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (!elem->edge_ref[5]) && (!elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (elem->edge_ref[9]) && (elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 38;
				elem->tem = 4;

			}
			else if ( (!elem->edge_ref[0]) && (!elem->edge_ref[1]) && (elem->edge_ref[2]) && (elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (!elem->edge_ref[5]) && (!elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (!elem->edge_ref[9]) && (elem->edge_ref[10]) && (elem->edge_ref[11])) {

				elem->pad = 39;
				elem->tem = 4;

			}
			else if ( (elem->edge_ref[0]) && (!elem->edge_ref[1]) && (!elem->edge_ref[2]) && (elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (!elem->edge_ref[5]) && (!elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(elem->edge_ref[8]) && (!elem->edge_ref[9]) && (!elem->edge_ref[10]) && (elem->edge_ref[11])) {

				elem->pad = 40;
				elem->tem = 4;

			}
			else if ( (elem->edge_ref[0]) && (!elem->edge_ref[1]) && (elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(elem->edge_ref[4]) && (!elem->edge_ref[5]) && (!elem->edge_ref[6]) && (elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (!elem->edge_ref[9]) && (!elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 41;
				elem->tem = 4;

			}
			else if ( (elem->edge_ref[0]) && (!elem->edge_ref[1]) && (elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (elem->edge_ref[5]) && (elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (!elem->edge_ref[9]) && (!elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 42;
				elem->tem = 4;

			}
			else if ( (!elem->edge_ref[0]) && (!elem->edge_ref[1]) && (!elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (elem->edge_ref[5]) && (elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(elem->edge_ref[8]) && (!elem->edge_ref[9]) && (elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 43;
				elem->tem = 4;

			}
			else if ( (!elem->edge_ref[0]) && (!elem->edge_ref[1]) && (!elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(elem->edge_ref[4]) && (!elem->edge_ref[5]) && (!elem->edge_ref[6]) && (elem->edge_ref[7]) &&
					(elem->edge_ref[8]) && (!elem->edge_ref[9]) && (elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 44;
				elem->tem = 4;

			}
			else if ( (!elem->edge_ref[0]) && (elem->edge_ref[1]) && (!elem->edge_ref[2]) && (elem->edge_ref[3]) &&
					(elem->edge_ref[4]) && (elem->edge_ref[5]) && (!elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (!elem->edge_ref[9]) && (!elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 45;
				elem->tem = 4;

			}
			else if ( (!elem->edge_ref[0]) && (elem->edge_ref[1]) && (!elem->edge_ref[2]) && (elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (!elem->edge_ref[5]) && (elem->edge_ref[6]) && (elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (!elem->edge_ref[9]) && (!elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 46;
				elem->tem = 4;

			}
			else if ( (!elem->edge_ref[0]) && (!elem->edge_ref[1]) && (!elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (!elem->edge_ref[5]) && (elem->edge_ref[6]) && (elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (elem->edge_ref[9]) && (!elem->edge_ref[10]) && (elem->edge_ref[11])) {

				elem->pad = 47;
				elem->tem = 4;

			}
			else if ( (!elem->edge_ref[0]) && (!elem->edge_ref[1]) && (!elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(elem->edge_ref[4]) && (elem->edge_ref[5]) && (!elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (elem->edge_ref[9]) && (!elem->edge_ref[10]) && (elem->edge_ref[11])) {

				elem->pad = 48;
				elem->tem = 4;

			}
		}
		//template 5
		if(ed_cont == 4 || ed_cont == 3 || ed_cont == 2){
			if ( (elem->edge_ref[0]) && (elem->edge_ref[1]) && (elem->edge_ref[2]) && (elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (!elem->edge_ref[5]) && (!elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (!elem->edge_ref[9]) && (!elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 49;
				elem->tem = 5;

			}
			else if ( (!elem->edge_ref[0]) && (!elem->edge_ref[1]) && (!elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (!elem->edge_ref[5]) && (!elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(elem->edge_ref[8]) && (elem->edge_ref[9]) && (elem->edge_ref[10]) && (elem->edge_ref[11])) {

				elem->pad = 50;
				elem->tem = 5;

			}
			else if ( (elem->edge_ref[0]) && (!elem->edge_ref[1]) && (!elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(elem->edge_ref[4]) && (elem->edge_ref[5]) && (!elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(elem->edge_ref[8]) && (!elem->edge_ref[9]) && (!elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 51;
				elem->tem = 5;

			}
			else if ( (!elem->edge_ref[0]) && (elem->edge_ref[1]) && (!elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (elem->edge_ref[5]) && (elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (elem->edge_ref[9]) && (!elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 52;
				elem->tem = 5;

			}
			else if ( (!elem->edge_ref[0]) && (!elem->edge_ref[1]) && (elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (!elem->edge_ref[5]) && (elem->edge_ref[6]) && (elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (!elem->edge_ref[9]) && (elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 53;
				elem->tem = 5;

			}
			else if ( (!elem->edge_ref[0]) && (!elem->edge_ref[1]) && (!elem->edge_ref[2]) && (elem->edge_ref[3]) &&
					(elem->edge_ref[4]) && (!elem->edge_ref[5]) && (!elem->edge_ref[6]) && (elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (!elem->edge_ref[9]) && (!elem->edge_ref[10]) && (elem->edge_ref[11])) {

				elem->pad = 54;
				elem->tem = 5;

			}


			//template 5 com 3
			else if ( (elem->edge_ref[0]) && (elem->edge_ref[1]) && (!elem->edge_ref[2]) && (elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (!elem->edge_ref[5]) && (!elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (!elem->edge_ref[9]) && (!elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 55;
				elem->tem = 5;

			}
			else if ( (elem->edge_ref[0]) && (elem->edge_ref[1]) && (elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (!elem->edge_ref[5]) && (!elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (!elem->edge_ref[9]) && (!elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 56;
				elem->tem = 5;

			}
			else if ( (!elem->edge_ref[0]) && (elem->edge_ref[1]) && (elem->edge_ref[2]) && (elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (!elem->edge_ref[5]) && (!elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (!elem->edge_ref[9]) && (!elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 57;
				elem->tem = 5;

			}
			else if ( (elem->edge_ref[0]) && (!elem->edge_ref[1]) && (elem->edge_ref[2]) && (elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (!elem->edge_ref[5]) && (!elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (!elem->edge_ref[9]) && (!elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 58;
				elem->tem = 5;

			}
			else if ( (!elem->edge_ref[0]) && (!elem->edge_ref[1]) && (!elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (!elem->edge_ref[5]) && (!elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(elem->edge_ref[8]) && (elem->edge_ref[9]) && (elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 59;
				elem->tem = 5;

			}
			else if ( (!elem->edge_ref[0]) && (!elem->edge_ref[1]) && (!elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (!elem->edge_ref[5]) && (!elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(elem->edge_ref[8]) && (elem->edge_ref[9]) && (!elem->edge_ref[10]) && (elem->edge_ref[11])) {

				elem->pad = 60;
				elem->tem = 5;

			}
			else if ( (!elem->edge_ref[0]) && (!elem->edge_ref[1]) && (!elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (!elem->edge_ref[5]) && (!elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (elem->edge_ref[9]) && (elem->edge_ref[10]) && (elem->edge_ref[11])) {

				elem->pad = 61;
				elem->tem = 5;

			}
			else if ( (!elem->edge_ref[0]) && (!elem->edge_ref[1]) && (!elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (!elem->edge_ref[5]) && (!elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(elem->edge_ref[8]) && (!elem->edge_ref[9]) && (elem->edge_ref[10]) && (elem->edge_ref[11])) {

				elem->pad = 62;
				elem->tem = 5;

			}
			else if ( (elem->edge_ref[0]) && (!elem->edge_ref[1]) && (!elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(elem->edge_ref[4]) && (elem->edge_ref[5]) && (!elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (!elem->edge_ref[9]) && (!elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 63;
				elem->tem = 5;

			}
			else if ( (elem->edge_ref[0]) && (!elem->edge_ref[1]) && (!elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (elem->edge_ref[5]) && (!elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(elem->edge_ref[8]) && (!elem->edge_ref[9]) && (!elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 64;
				elem->tem = 5;

			}
			else if ( (!elem->edge_ref[0]) && (!elem->edge_ref[1]) && (!elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(elem->edge_ref[4]) && (elem->edge_ref[5]) && (!elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(elem->edge_ref[8]) && (!elem->edge_ref[9]) && (!elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 65;
				elem->tem = 5;

			}
			else if ( (elem->edge_ref[0]) && (!elem->edge_ref[1]) && (!elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(elem->edge_ref[4]) && (!elem->edge_ref[5]) && (!elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(elem->edge_ref[8]) && (!elem->edge_ref[9]) && (!elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 66;
				elem->tem = 5;

			}
			else if ( (!elem->edge_ref[0]) && (elem->edge_ref[1]) && (!elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (elem->edge_ref[5]) && (elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (!elem->edge_ref[9]) && (!elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 67;
				elem->tem = 5;

			}
			else if ( (!elem->edge_ref[0]) && (elem->edge_ref[1]) && (!elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (elem->edge_ref[5]) && (!elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (elem->edge_ref[9]) && (!elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 68;
				elem->tem = 5;

			}
			else if ( (!elem->edge_ref[0]) && (!elem->edge_ref[1]) && (!elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (elem->edge_ref[5]) && (elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (elem->edge_ref[9]) && (!elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 69;
				elem->tem = 5;

			}
			else if ( (!elem->edge_ref[0]) && (elem->edge_ref[1]) && (!elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (!elem->edge_ref[5]) && (elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (elem->edge_ref[9]) && (!elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 70;
				elem->tem = 5;

			}
			else if ( (!elem->edge_ref[0]) && (!elem->edge_ref[1]) && (elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (!elem->edge_ref[5]) && (elem->edge_ref[6]) && (elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (!elem->edge_ref[9]) && (!elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 71;
				elem->tem = 5;

			}
			else if ( (!elem->edge_ref[0]) && (!elem->edge_ref[1]) && (elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (!elem->edge_ref[5]) && (elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (!elem->edge_ref[9]) && (elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 72;
				elem->tem = 5;

			}
			else if ( (!elem->edge_ref[0]) && (!elem->edge_ref[1]) && (!elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (!elem->edge_ref[5]) && (elem->edge_ref[6]) && (elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (!elem->edge_ref[9]) && (elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 73;
				elem->tem = 5;

			}
			else if ( (!elem->edge_ref[0]) && (!elem->edge_ref[1]) && (elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (!elem->edge_ref[5]) && (!elem->edge_ref[6]) && (elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (!elem->edge_ref[9]) && (elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 74;
				elem->tem = 5;

			}
			else if ( (!elem->edge_ref[0]) && (!elem->edge_ref[1]) && (!elem->edge_ref[2]) && (elem->edge_ref[3]) &&
					(elem->edge_ref[4]) && (!elem->edge_ref[5]) && (!elem->edge_ref[6]) && (elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (!elem->edge_ref[9]) && (!elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 75;
				elem->tem = 5;

			}
			else if ( (!elem->edge_ref[0]) && (!elem->edge_ref[1]) && (!elem->edge_ref[2]) && (elem->edge_ref[3]) &&
					(elem->edge_ref[4]) && (!elem->edge_ref[5]) && (!elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (!elem->edge_ref[9]) && (!elem->edge_ref[10]) && (elem->edge_ref[11])) {

				elem->pad = 76;
				elem->tem = 5;

			}
			else if ( (!elem->edge_ref[0]) && (!elem->edge_ref[1]) && (!elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(elem->edge_ref[4]) && (!elem->edge_ref[5]) && (!elem->edge_ref[6]) && (elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (!elem->edge_ref[9]) && (!elem->edge_ref[10]) && (elem->edge_ref[11])) {

				elem->pad = 77;
				elem->tem = 5;

			}
			else if ( (!elem->edge_ref[0]) && (!elem->edge_ref[1]) && (!elem->edge_ref[2]) && (elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (!elem->edge_ref[5]) && (!elem->edge_ref[6]) && (elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (!elem->edge_ref[9]) && (!elem->edge_ref[10]) && (elem->edge_ref[11])) {

				elem->pad = 78;
				elem->tem = 5;

			}

			// templete 5 com 2
			else if ( (elem->edge_ref[0]) && (elem->edge_ref[1]) && (!elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (!elem->edge_ref[5]) && (!elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (!elem->edge_ref[9]) && (!elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 79;
				elem->tem = 5;

			}
			else if ( (!elem->edge_ref[0]) && (elem->edge_ref[1]) && (elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (!elem->edge_ref[5]) && (!elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (!elem->edge_ref[9]) && (!elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 80;
				elem->tem = 5;

			}
			else if ( (!elem->edge_ref[0]) && (!elem->edge_ref[1]) && (elem->edge_ref[2]) && (elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (!elem->edge_ref[5]) && (!elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (!elem->edge_ref[9]) && (!elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 81;
				elem->tem = 5;

			}
			else if ( (elem->edge_ref[0]) && (!elem->edge_ref[1]) && (!elem->edge_ref[2]) && (elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (!elem->edge_ref[5]) && (!elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (!elem->edge_ref[9]) && (!elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 82;
				elem->tem = 5;

			}
			else if ( (!elem->edge_ref[0]) && (!elem->edge_ref[1]) && (!elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (!elem->edge_ref[5]) && (!elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(elem->edge_ref[8]) && (elem->edge_ref[9]) && (!elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 83;
				elem->tem = 5;

			}
			else if ( (!elem->edge_ref[0]) && (!elem->edge_ref[1]) && (!elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (!elem->edge_ref[5]) && (!elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (elem->edge_ref[9]) && (elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 84;
				elem->tem = 5;

			}
			else if ( (!elem->edge_ref[0]) && (!elem->edge_ref[1]) && (!elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (!elem->edge_ref[5]) && (!elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (!elem->edge_ref[9]) && (elem->edge_ref[10]) && (elem->edge_ref[11])) {

				elem->pad = 85;
				elem->tem = 5;

			}
			else if ( (!elem->edge_ref[0]) && (!elem->edge_ref[1]) && (!elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (!elem->edge_ref[5]) && (!elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(elem->edge_ref[8]) && (!elem->edge_ref[9]) && (!elem->edge_ref[10]) && (elem->edge_ref[11])) {

				elem->pad = 86;
				elem->tem = 5;

			}
			else if ( (!elem->edge_ref[0]) && (!elem->edge_ref[1]) && (!elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (elem->edge_ref[5]) && (!elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (elem->edge_ref[9]) && (!elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 87;
				elem->tem = 5;

			}
			else if ( (!elem->edge_ref[0]) && (!elem->edge_ref[1]) && (!elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (!elem->edge_ref[5]) && (elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (elem->edge_ref[9]) && (!elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 88;
				elem->tem = 5;

			}
			else if ( (!elem->edge_ref[0]) && (elem->edge_ref[1]) && (!elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (!elem->edge_ref[5]) && (elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (!elem->edge_ref[9]) && (!elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 89;
				elem->tem = 5;

			}
			else if ( (!elem->edge_ref[0]) && (elem->edge_ref[1]) && (!elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (elem->edge_ref[5]) && (!elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (!elem->edge_ref[9]) && (!elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 90;
				elem->tem = 5;

			}
			else if ( (!elem->edge_ref[0]) && (!elem->edge_ref[1]) && (elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (!elem->edge_ref[5]) && (elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (!elem->edge_ref[9]) && (!elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 91;
				elem->tem = 5;

			}
			else if ( (!elem->edge_ref[0]) && (!elem->edge_ref[1]) && (!elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (!elem->edge_ref[5]) && (elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (!elem->edge_ref[9]) && (elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 92;
				elem->tem = 5;

			}
			else if ( (!elem->edge_ref[0]) && (!elem->edge_ref[1]) && (!elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (!elem->edge_ref[5]) && (!elem->edge_ref[6]) && (elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (!elem->edge_ref[9]) && (elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 93;
				elem->tem = 5;

			}
			else if ( (!elem->edge_ref[0]) && (!elem->edge_ref[1]) && (elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (!elem->edge_ref[5]) && (!elem->edge_ref[6]) && (elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (!elem->edge_ref[9]) && (!elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 94;
				elem->tem = 5;

			}
			else if ( (!elem->edge_ref[0]) && (!elem->edge_ref[1]) && (!elem->edge_ref[2]) && (elem->edge_ref[3]) &&
					(elem->edge_ref[4]) && (!elem->edge_ref[5]) && (!elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (!elem->edge_ref[9]) && (!elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 95;
				elem->tem = 5;

			}
			else if ( (!elem->edge_ref[0]) && (!elem->edge_ref[1]) && (!elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(elem->edge_ref[4]) && (!elem->edge_ref[5]) && (elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (!elem->edge_ref[9]) && (!elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 96;
				elem->tem = 5;

			}
			else if ( (!elem->edge_ref[0]) && (!elem->edge_ref[1]) && (!elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (!elem->edge_ref[5]) && (elem->edge_ref[6]) && (elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (!elem->edge_ref[9]) && (!elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 97;
				elem->tem = 5;

			}
			else if ( (!elem->edge_ref[0]) && (!elem->edge_ref[1]) && (!elem->edge_ref[2]) && (elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (!elem->edge_ref[5]) && (!elem->edge_ref[6]) && (elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (!elem->edge_ref[9]) && (!elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 98;
				elem->tem = 5;

			}
		}
		//template 6
		if(ed_cont == 6){
			if ( (elem->edge_ref[0]) && (elem->edge_ref[1]) && (!elem->edge_ref[2]) && (elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (!elem->edge_ref[5]) && (!elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(elem->edge_ref[8]) && (elem->edge_ref[9]) && (!elem->edge_ref[10]) && (elem->edge_ref[11])) {

				elem->pad = 99;
				elem->tem = 6;

			}
			else if ( (elem->edge_ref[0]) && (elem->edge_ref[1]) && (elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (!elem->edge_ref[5]) && (!elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(elem->edge_ref[8]) && (elem->edge_ref[9]) && (elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 100;
				elem->tem = 6;

			}
			else if ( (!elem->edge_ref[0]) && (elem->edge_ref[1]) && (elem->edge_ref[2]) && (elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (!elem->edge_ref[5]) && (!elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (elem->edge_ref[9]) && (elem->edge_ref[10]) && (elem->edge_ref[11])) {

				elem->pad = 101;
				elem->tem = 6;

			}
			else if ( (elem->edge_ref[0]) && (!elem->edge_ref[1]) && (elem->edge_ref[2]) && (elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (!elem->edge_ref[5]) && (!elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(elem->edge_ref[8]) && (!elem->edge_ref[9]) && (elem->edge_ref[10]) && (elem->edge_ref[11])) {

				elem->pad = 102;
				elem->tem = 6;

			}
			else if ( (elem->edge_ref[0]) && (!elem->edge_ref[1]) && (elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (elem->edge_ref[5]) && (elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(elem->edge_ref[8]) && (!elem->edge_ref[9]) && (elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 103;
				elem->tem = 6;

			}
			else if ( (!elem->edge_ref[0]) && (!elem->edge_ref[1]) && (!elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(elem->edge_ref[4]) && (elem->edge_ref[5]) && (elem->edge_ref[6]) && (elem->edge_ref[7]) &&
					(elem->edge_ref[8]) && (!elem->edge_ref[9]) && (elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 104;
				elem->tem = 6;

			}
			else if ( (elem->edge_ref[0]) && (!elem->edge_ref[1]) && (elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(elem->edge_ref[4]) && (!elem->edge_ref[5]) && (!elem->edge_ref[6]) && (elem->edge_ref[7]) &&
					(elem->edge_ref[8]) && (!elem->edge_ref[9]) && (elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 105;
				elem->tem = 6;

			}
			else if ( (elem->edge_ref[0]) && (!elem->edge_ref[1]) && (elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(elem->edge_ref[4]) && (elem->edge_ref[5]) && (elem->edge_ref[6]) && (elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (!elem->edge_ref[9]) && (!elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 106;
				elem->tem = 6;

			}
			else if ( (!elem->edge_ref[0]) && (elem->edge_ref[1]) && (!elem->edge_ref[2]) && (elem->edge_ref[3]) &&
					(elem->edge_ref[4]) && (elem->edge_ref[5]) && (elem->edge_ref[6]) && (elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (!elem->edge_ref[9]) && (!elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 107;
				elem->tem = 6;

			}
			else if ( (!elem->edge_ref[0]) && (elem->edge_ref[1]) && (!elem->edge_ref[2]) && (elem->edge_ref[3]) &&
					(elem->edge_ref[4]) && (elem->edge_ref[5]) && (!elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (elem->edge_ref[9]) && (!elem->edge_ref[10]) && (elem->edge_ref[11])) {

				elem->pad = 108;
				elem->tem = 6;

			}
			else if ( (!elem->edge_ref[0]) && (!elem->edge_ref[1]) && (!elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(elem->edge_ref[4]) && (elem->edge_ref[5]) && (elem->edge_ref[6]) && (elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (elem->edge_ref[9]) && (!elem->edge_ref[10]) && (elem->edge_ref[11])) {

				elem->pad = 109;
				elem->tem = 6;

			}
			else if ( (!elem->edge_ref[0]) && (elem->edge_ref[1]) && (!elem->edge_ref[2]) && (elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (!elem->edge_ref[5]) && (elem->edge_ref[6]) && (elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (elem->edge_ref[9]) && (!elem->edge_ref[10]) && (elem->edge_ref[11])) {

				elem->pad = 110;
				elem->tem = 6;

			}
		}
		//template 7
		if(ed_cont == 3){
			if ( (elem->edge_ref[0]) && (!elem->edge_ref[1]) && (!elem->edge_ref[2]) && (elem->edge_ref[3]) &&
					(elem->edge_ref[4]) && (!elem->edge_ref[5]) && (!elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (!elem->edge_ref[9]) && (!elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 111;
				elem->tem = 7;

			}
			else if ( (elem->edge_ref[0]) && (elem->edge_ref[1]) && (!elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (elem->edge_ref[5]) && (!elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (!elem->edge_ref[9]) && (!elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 112;
				elem->tem = 7;

			}
			else if ( (!elem->edge_ref[0]) && (elem->edge_ref[1]) && (elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (!elem->edge_ref[5]) && (elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (!elem->edge_ref[9]) && (!elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 113;
				elem->tem = 7;

			}
			else if ( (!elem->edge_ref[0]) && (!elem->edge_ref[1]) && (elem->edge_ref[2]) && (elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (!elem->edge_ref[5]) && (!elem->edge_ref[6]) && (elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (!elem->edge_ref[9]) && (!elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 114;
				elem->tem = 7;

			}
			else if ( (!elem->edge_ref[0]) && (!elem->edge_ref[1]) && (!elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(elem->edge_ref[4]) && (!elem->edge_ref[5]) && (!elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(elem->edge_ref[8]) && (!elem->edge_ref[9]) && (!elem->edge_ref[10]) && (elem->edge_ref[11])) {

				elem->pad = 115;
				elem->tem = 7;

			}
			else if ( (!elem->edge_ref[0]) && (!elem->edge_ref[1]) && (!elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (elem->edge_ref[5]) && (!elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(elem->edge_ref[8]) && (elem->edge_ref[9]) && (!elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 116;
				elem->tem = 7;

			}
			else if ( (!elem->edge_ref[0]) && (!elem->edge_ref[1]) && (!elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (!elem->edge_ref[5]) && (elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (elem->edge_ref[9]) && (elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 117;
				elem->tem = 7;

			}
			else if ( (!elem->edge_ref[0]) && (!elem->edge_ref[1]) && (!elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (!elem->edge_ref[5]) && (!elem->edge_ref[6]) && (elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (!elem->edge_ref[9]) && (elem->edge_ref[10]) && (elem->edge_ref[11])) {

				elem->pad = 118;
				elem->tem = 7;

			}
		}
		//template 8
		if(ed_cont == 8){
			if ( (elem->edge_ref[0]) && (elem->edge_ref[1]) && (elem->edge_ref[2]) && (elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (!elem->edge_ref[5]) && (!elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(elem->edge_ref[8]) && (elem->edge_ref[9]) && (elem->edge_ref[10]) && (elem->edge_ref[11])) {

				elem->pad = 119;
				elem->tem = 8;

			}
			else if ( (elem->edge_ref[0]) && (!elem->edge_ref[1]) && (elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(elem->edge_ref[4]) && (elem->edge_ref[5]) && (elem->edge_ref[6]) && (elem->edge_ref[7]) &&
					(elem->edge_ref[8]) && (!elem->edge_ref[9]) && (elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 120;
				elem->tem = 8;

			}
			else if ( (!elem->edge_ref[0]) && (elem->edge_ref[1]) && (!elem->edge_ref[2]) && (elem->edge_ref[3]) &&
					(elem->edge_ref[4]) && (elem->edge_ref[5]) && (elem->edge_ref[6]) && (elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (elem->edge_ref[9]) && (!elem->edge_ref[10]) && (elem->edge_ref[11])) {

				elem->pad = 121;
				elem->tem = 8;

			}
		}
		//template 9
		if(ed_cont == 8){
			if ( (elem->edge_ref[0]) && (elem->edge_ref[1]) && (elem->edge_ref[2]) && (elem->edge_ref[3]) &&
					(elem->edge_ref[4]) && (elem->edge_ref[5]) && (elem->edge_ref[6]) && (elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (!elem->edge_ref[9]) && (!elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 122;
				elem->tem = 9;

			}
			else if ( (!elem->edge_ref[0]) && (!elem->edge_ref[1]) && (!elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(elem->edge_ref[4]) && (elem->edge_ref[5]) && (elem->edge_ref[6]) && (elem->edge_ref[7]) &&
					(elem->edge_ref[8]) && (elem->edge_ref[9]) && (elem->edge_ref[10]) && (elem->edge_ref[11])) {

				elem->pad = 123;
				elem->tem = 9;

			}
			else if ( (elem->edge_ref[0]) && (elem->edge_ref[1]) && (!elem->edge_ref[2]) && (elem->edge_ref[3]) &&
					(elem->edge_ref[4]) && (elem->edge_ref[5]) && (!elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(elem->edge_ref[8]) && (elem->edge_ref[9]) && (!elem->edge_ref[10]) && (elem->edge_ref[11])) {

				elem->pad = 124;
				elem->tem = 9;

			}
			else if ( (!elem->edge_ref[0]) && (elem->edge_ref[1]) && (elem->edge_ref[2]) && (elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (!elem->edge_ref[5]) && (elem->edge_ref[6]) && (elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (elem->edge_ref[9]) && (elem->edge_ref[10]) && (elem->edge_ref[11])) {

				elem->pad = 125;
				elem->tem = 9;

			}
			else if ( (elem->edge_ref[0]) && (elem->edge_ref[1]) && (elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (elem->edge_ref[5]) && (elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(elem->edge_ref[8]) && (elem->edge_ref[9]) && (elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 126;
				elem->tem = 9;

			}
			else if ( (elem->edge_ref[0]) && (!elem->edge_ref[1]) && (elem->edge_ref[2]) && (elem->edge_ref[3]) &&
					(elem->edge_ref[4]) && (!elem->edge_ref[5]) && (!elem->edge_ref[6]) && (elem->edge_ref[7]) &&
					(elem->edge_ref[8]) && (!elem->edge_ref[9]) && (elem->edge_ref[10]) && (elem->edge_ref[11])) {

				elem->pad = 127;
				elem->tem = 9;

			}
		}
		//template 10
		if(ed_cont == 5){
			if ( (elem->edge_ref[0]) && (elem->edge_ref[1]) && (!elem->edge_ref[2]) && (elem->edge_ref[3]) &&
					(elem->edge_ref[4]) && (elem->edge_ref[5]) && (!elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (!elem->edge_ref[9]) && (!elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 128;
				elem->tem = 10;

			}
			else if ( (elem->edge_ref[0]) && (elem->edge_ref[1]) && (elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (elem->edge_ref[5]) && (elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (!elem->edge_ref[9]) && (!elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 129;
				elem->tem = 10;

			}
			else if ( (!elem->edge_ref[0]) && (elem->edge_ref[1]) && (elem->edge_ref[2]) && (elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (!elem->edge_ref[5]) && (elem->edge_ref[6]) && (elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (!elem->edge_ref[9]) && (!elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 130;
				elem->tem = 10;

			}
			else if ( (elem->edge_ref[0]) && (!elem->edge_ref[1]) && (elem->edge_ref[2]) && (elem->edge_ref[3]) &&
					(elem->edge_ref[4]) && (!elem->edge_ref[5]) && (!elem->edge_ref[6]) && (elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (!elem->edge_ref[9]) && (!elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 131;
				elem->tem = 10;

			}
			else if ( (elem->edge_ref[0]) && (!elem->edge_ref[1]) && (!elem->edge_ref[2]) && (elem->edge_ref[3]) &&
					(elem->edge_ref[4]) && (!elem->edge_ref[5]) && (!elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(elem->edge_ref[8]) && (!elem->edge_ref[9]) && (!elem->edge_ref[10]) && (elem->edge_ref[11])) {

				elem->pad = 132;
				elem->tem = 10;

			}
			else if ( (elem->edge_ref[0]) && (elem->edge_ref[1]) && (!elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (elem->edge_ref[5]) && (!elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(elem->edge_ref[8]) && (elem->edge_ref[9]) && (!elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 133;
				elem->tem = 10;

			}
			else if ( (!elem->edge_ref[0]) && (elem->edge_ref[1]) && (elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (!elem->edge_ref[5]) && (elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (elem->edge_ref[9]) && (elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 134;
				elem->tem = 10;

			}
			else if ( (!elem->edge_ref[0]) && (!elem->edge_ref[1]) && (elem->edge_ref[2]) && (elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (!elem->edge_ref[5]) && (!elem->edge_ref[6]) && (elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (!elem->edge_ref[9]) && (elem->edge_ref[10]) && (elem->edge_ref[11])) {

				elem->pad = 135;
				elem->tem = 10;

			}
			else if ( (!elem->edge_ref[0]) && (!elem->edge_ref[1]) && (!elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(elem->edge_ref[4]) && (elem->edge_ref[5]) && (!elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(elem->edge_ref[8]) && (elem->edge_ref[9]) && (!elem->edge_ref[10]) && (elem->edge_ref[11])) {

				elem->pad = 136;
				elem->tem = 10;

			}
			else if ( (!elem->edge_ref[0]) && (!elem->edge_ref[1]) && (!elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (elem->edge_ref[5]) && (elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(elem->edge_ref[8]) && (elem->edge_ref[9]) && (elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 137;
				elem->tem = 10;

			}
			else if ( (!elem->edge_ref[0]) && (!elem->edge_ref[1]) && (!elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (!elem->edge_ref[5]) && (elem->edge_ref[6]) && (elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (elem->edge_ref[9]) && (elem->edge_ref[10]) && (elem->edge_ref[11])) {

				elem->pad = 138;
				elem->tem = 10;

			}
			else if ( (!elem->edge_ref[0]) && (!elem->edge_ref[1]) && (!elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(elem->edge_ref[4]) && (!elem->edge_ref[5]) && (!elem->edge_ref[6]) && (elem->edge_ref[7]) &&
					(elem->edge_ref[8]) && (!elem->edge_ref[9]) && (elem->edge_ref[10]) && (elem->edge_ref[11])) {

				elem->pad = 139;
				elem->tem = 10;

			}
		}


		//elem->pad=-1;
		if(elem->pad>=10 & elem->pad<=21){el_1++;}
		if(elem->pad>=22 & elem->pad<=24){el_2++;}
		if(elem->pad>=25 & elem->pad<=36){el_3++;}
		if(elem->pad>=37 & elem->pad<=48){el_4++;}
		if(elem->pad>=49 & elem->pad<=98){el_5++;}
		if(elem->pad>=99 & elem->pad<=110){el_6++;}
		if(elem->pad>=111 & elem->pad<=118){el_7++;}
		if(elem->pad>=119 & elem->pad<=121){el_8++;}
		if(elem->pad>=122 & elem->pad<=127){el_9++;}
		if(elem->pad>=128 & elem->pad<=139){el_10++;}
		if(elem->pad==140 ){el_11++;}
	}



	int su = 0;
	su = el_1+el_2+el_3+el_4+el_5+el_6+el_7+el_8+el_9+el_10+el_11;

	if(true){
		printf("case_1: %d\n",el_1);
		printf("case_2: %d\n",el_2);
		printf("case_3: %d\n",el_3);
		printf("case_4: %d\n",el_4);
		printf("case_5: %d\n",el_5);
		printf("case_6: %d\n",el_6);
		printf("case_7: %d\n",el_7);
		printf("case_8: %d\n",el_8);
		printf("case_9: %d\n",el_9);
		printf("case_10: %d\n",el_10);
		printf("case_11: %d\n",el_11);
		printf("sum: %d\n",su);
		printf("case_0: %d\n",el_0);
	}
}



int AddPointOnEdge(int* nodes, sc_hash_array_t* hash, int &npoints, GtsPoint *p, std::vector<double> &coords) {
	size_t position;
	node_in_edge_t *r;
	node_in_edge_t key;
	key.coord[0] = nodes[0];
	key.coord[1] = nodes[1];

	r = (node_in_edge_t*) sc_hash_array_insert_unique(hash, &key, &position);
	if (r != NULL) {
		r->coord[0] = key.coord[0];
		r->coord[1] = key.coord[1];
		r->node_id = npoints;
		npoints++;
		coords.push_back(p->x);
		coords.push_back(p->y);
		coords.push_back(p->z);
		return r->node_id;
	} else {
		r = (node_in_edge_t*) sc_array_index(&hash->a, position);
		return r->node_id;
	}

}

void Edge_identification(hexa_tree_t* mesh, std::vector<int>& elements_ids, std::vector<int>& edges_ids) {

	for (int iel = 0; iel < elements_ids.size(); ++iel) {

		octant_t *elem = (octant_t*) sc_array_index(&mesh->elements, elements_ids[iel]);

		//template 11
		if(elem->pad==140){
			for (int edge = 0; edge < 12; ++edge) {
				edges_ids.push_back(elem->edge_id[edge]);
			}
		}
		//template 10
		else if(elem->pad==139){
			edges_ids.push_back(elem->edge_id[4]);
			edges_ids.push_back(elem->edge_id[7]);
			edges_ids.push_back(elem->edge_id[8]);
			edges_ids.push_back(elem->edge_id[10]);
			edges_ids.push_back(elem->edge_id[11]);
		}else if(elem->pad==138){
			edges_ids.push_back(elem->edge_id[6]);
			edges_ids.push_back(elem->edge_id[7]);
			edges_ids.push_back(elem->edge_id[9]);
			edges_ids.push_back(elem->edge_id[10]);
			edges_ids.push_back(elem->edge_id[11]);
		}else if(elem->pad==137){
			edges_ids.push_back(elem->edge_id[5]);
			edges_ids.push_back(elem->edge_id[6]);
			edges_ids.push_back(elem->edge_id[8]);
			edges_ids.push_back(elem->edge_id[9]);
			edges_ids.push_back(elem->edge_id[10]);
		}else if(elem->pad==136){
			edges_ids.push_back(elem->edge_id[4]);
			edges_ids.push_back(elem->edge_id[5]);
			edges_ids.push_back(elem->edge_id[8]);
			edges_ids.push_back(elem->edge_id[9]);
			edges_ids.push_back(elem->edge_id[11]);
		}else if(elem->pad==135){
			edges_ids.push_back(elem->edge_id[2]);
			edges_ids.push_back(elem->edge_id[3]);
			edges_ids.push_back(elem->edge_id[7]);
			edges_ids.push_back(elem->edge_id[10]);
			edges_ids.push_back(elem->edge_id[11]);
		}else if(elem->pad==134){
			edges_ids.push_back(elem->edge_id[1]);
			edges_ids.push_back(elem->edge_id[2]);
			edges_ids.push_back(elem->edge_id[6]);
			edges_ids.push_back(elem->edge_id[9]);
			edges_ids.push_back(elem->edge_id[10]);
		}else if(elem->pad==133){
			edges_ids.push_back(elem->edge_id[0]);
			edges_ids.push_back(elem->edge_id[1]);
			edges_ids.push_back(elem->edge_id[5]);
			edges_ids.push_back(elem->edge_id[8]);
			edges_ids.push_back(elem->edge_id[9]);
		}else if(elem->pad==132){
			edges_ids.push_back(elem->edge_id[0]);
			edges_ids.push_back(elem->edge_id[3]);
			edges_ids.push_back(elem->edge_id[4]);
			edges_ids.push_back(elem->edge_id[8]);
			edges_ids.push_back(elem->edge_id[11]);
		}else if(elem->pad==131){
			edges_ids.push_back(elem->edge_id[0]);
			edges_ids.push_back(elem->edge_id[2]);
			edges_ids.push_back(elem->edge_id[3]);
			edges_ids.push_back(elem->edge_id[4]);
			edges_ids.push_back(elem->edge_id[7]);
		}else if(elem->pad==130){
			edges_ids.push_back(elem->edge_id[1]);
			edges_ids.push_back(elem->edge_id[2]);
			edges_ids.push_back(elem->edge_id[3]);
			edges_ids.push_back(elem->edge_id[6]);
			edges_ids.push_back(elem->edge_id[7]);
		}else if(elem->pad==129){
			edges_ids.push_back(elem->edge_id[0]);
			edges_ids.push_back(elem->edge_id[1]);
			edges_ids.push_back(elem->edge_id[2]);
			edges_ids.push_back(elem->edge_id[5]);
			edges_ids.push_back(elem->edge_id[6]);
		}else if(elem->pad==128){
			edges_ids.push_back(elem->edge_id[0]);
			edges_ids.push_back(elem->edge_id[1]);
			edges_ids.push_back(elem->edge_id[3]);
			edges_ids.push_back(elem->edge_id[4]);
			edges_ids.push_back(elem->edge_id[5]);
		}
		//template 9
		else if(elem->pad==127){
			edges_ids.push_back(elem->edge_id[0]);
			edges_ids.push_back(elem->edge_id[2]);
			edges_ids.push_back(elem->edge_id[3]);
			edges_ids.push_back(elem->edge_id[4]);
			edges_ids.push_back(elem->edge_id[7]);
			edges_ids.push_back(elem->edge_id[8]);
			edges_ids.push_back(elem->edge_id[10]);
			edges_ids.push_back(elem->edge_id[11]);
		}else if(elem->pad==126){
			edges_ids.push_back(elem->edge_id[0]);
			edges_ids.push_back(elem->edge_id[1]);
			edges_ids.push_back(elem->edge_id[2]);
			edges_ids.push_back(elem->edge_id[5]);
			edges_ids.push_back(elem->edge_id[6]);
			edges_ids.push_back(elem->edge_id[8]);
			edges_ids.push_back(elem->edge_id[9]);
			edges_ids.push_back(elem->edge_id[10]);
		}else if(elem->pad==125){
			edges_ids.push_back(elem->edge_id[1]);
			edges_ids.push_back(elem->edge_id[2]);
			edges_ids.push_back(elem->edge_id[3]);
			edges_ids.push_back(elem->edge_id[6]);
			edges_ids.push_back(elem->edge_id[7]);
			edges_ids.push_back(elem->edge_id[9]);
			edges_ids.push_back(elem->edge_id[10]);
			edges_ids.push_back(elem->edge_id[11]);
		}else if(elem->pad==124){
			edges_ids.push_back(elem->edge_id[0]);
			edges_ids.push_back(elem->edge_id[1]);
			edges_ids.push_back(elem->edge_id[3]);
			edges_ids.push_back(elem->edge_id[4]);
			edges_ids.push_back(elem->edge_id[5]);
			edges_ids.push_back(elem->edge_id[8]);
			edges_ids.push_back(elem->edge_id[9]);
			edges_ids.push_back(elem->edge_id[11]);
		}else if(elem->pad==123){
			edges_ids.push_back(elem->edge_id[4]);
			edges_ids.push_back(elem->edge_id[5]);
			edges_ids.push_back(elem->edge_id[6]);
			edges_ids.push_back(elem->edge_id[7]);
			edges_ids.push_back(elem->edge_id[8]);
			edges_ids.push_back(elem->edge_id[9]);
			edges_ids.push_back(elem->edge_id[10]);
			edges_ids.push_back(elem->edge_id[11]);
		}else if(elem->pad==122){
			edges_ids.push_back(elem->edge_id[0]);
			edges_ids.push_back(elem->edge_id[1]);
			edges_ids.push_back(elem->edge_id[2]);
			edges_ids.push_back(elem->edge_id[3]);
			edges_ids.push_back(elem->edge_id[4]);
			edges_ids.push_back(elem->edge_id[5]);
			edges_ids.push_back(elem->edge_id[6]);
			edges_ids.push_back(elem->edge_id[7]);
		}
		//template 8
		else if(elem->pad==121){
			edges_ids.push_back(elem->edge_id[1]);
			edges_ids.push_back(elem->edge_id[3]);
			edges_ids.push_back(elem->edge_id[4]);
			edges_ids.push_back(elem->edge_id[5]);
			edges_ids.push_back(elem->edge_id[6]);
			edges_ids.push_back(elem->edge_id[7]);
			edges_ids.push_back(elem->edge_id[9]);
			edges_ids.push_back(elem->edge_id[11]);
		}else if(elem->pad==120){
			edges_ids.push_back(elem->edge_id[0]);
			edges_ids.push_back(elem->edge_id[2]);
			edges_ids.push_back(elem->edge_id[4]);
			edges_ids.push_back(elem->edge_id[5]);
			edges_ids.push_back(elem->edge_id[6]);
			edges_ids.push_back(elem->edge_id[7]);
			edges_ids.push_back(elem->edge_id[8]);
			edges_ids.push_back(elem->edge_id[10]);
		}else if(elem->pad==119){
			edges_ids.push_back(elem->edge_id[0]);
			edges_ids.push_back(elem->edge_id[1]);
			edges_ids.push_back(elem->edge_id[2]);
			edges_ids.push_back(elem->edge_id[3]);
			edges_ids.push_back(elem->edge_id[8]);
			edges_ids.push_back(elem->edge_id[9]);
			edges_ids.push_back(elem->edge_id[10]);
			edges_ids.push_back(elem->edge_id[11]);
		}
		//template 7
		else if(elem->pad==118){
			edges_ids.push_back(elem->edge_id[7]);
			edges_ids.push_back(elem->edge_id[10]);
			edges_ids.push_back(elem->edge_id[11]);
		}else if(elem->pad==117){
			edges_ids.push_back(elem->edge_id[6]);
			edges_ids.push_back(elem->edge_id[9]);
			edges_ids.push_back(elem->edge_id[10]);
		}else if(elem->pad==116){
			edges_ids.push_back(elem->edge_id[5]);
			edges_ids.push_back(elem->edge_id[8]);
			edges_ids.push_back(elem->edge_id[9]);
		}else if(elem->pad==115){
			edges_ids.push_back(elem->edge_id[4]);
			edges_ids.push_back(elem->edge_id[8]);
			edges_ids.push_back(elem->edge_id[11]);
		}else if(elem->pad==114){
			edges_ids.push_back(elem->edge_id[2]);
			edges_ids.push_back(elem->edge_id[3]);
			edges_ids.push_back(elem->edge_id[7]);
		}else if(elem->pad==113){
			edges_ids.push_back(elem->edge_id[1]);
			edges_ids.push_back(elem->edge_id[2]);
			edges_ids.push_back(elem->edge_id[6]);
		}else if(elem->pad==112){
			edges_ids.push_back(elem->edge_id[0]);
			edges_ids.push_back(elem->edge_id[1]);
			edges_ids.push_back(elem->edge_id[5]);
		}else if(elem->pad==111){
			edges_ids.push_back(elem->edge_id[0]);
			edges_ids.push_back(elem->edge_id[3]);
			edges_ids.push_back(elem->edge_id[4]);
		}
		//template 6
		else if(elem->pad==110){
			edges_ids.push_back(elem->edge_id[1]);
			edges_ids.push_back(elem->edge_id[3]);
			edges_ids.push_back(elem->edge_id[6]);
			edges_ids.push_back(elem->edge_id[7]);
			edges_ids.push_back(elem->edge_id[9]);
			edges_ids.push_back(elem->edge_id[11]);
		}else if(elem->pad==109){
			edges_ids.push_back(elem->edge_id[4]);
			edges_ids.push_back(elem->edge_id[5]);
			edges_ids.push_back(elem->edge_id[6]);
			edges_ids.push_back(elem->edge_id[7]);
			edges_ids.push_back(elem->edge_id[9]);
			edges_ids.push_back(elem->edge_id[11]);
		}else if(elem->pad==108){
			edges_ids.push_back(elem->edge_id[1]);
			edges_ids.push_back(elem->edge_id[3]);
			edges_ids.push_back(elem->edge_id[4]);
			edges_ids.push_back(elem->edge_id[5]);
			edges_ids.push_back(elem->edge_id[8]);
			edges_ids.push_back(elem->edge_id[11]);
		}else if(elem->pad==107){
			edges_ids.push_back(elem->edge_id[1]);
			edges_ids.push_back(elem->edge_id[3]);
			edges_ids.push_back(elem->edge_id[4]);
			edges_ids.push_back(elem->edge_id[5]);
			edges_ids.push_back(elem->edge_id[6]);
			edges_ids.push_back(elem->edge_id[7]);
		}else if(elem->pad==106){
			edges_ids.push_back(elem->edge_id[0]);
			edges_ids.push_back(elem->edge_id[2]);
			edges_ids.push_back(elem->edge_id[4]);
			edges_ids.push_back(elem->edge_id[5]);
			edges_ids.push_back(elem->edge_id[6]);
			edges_ids.push_back(elem->edge_id[7]);
		}else if(elem->pad==105){
			edges_ids.push_back(elem->edge_id[0]);
			edges_ids.push_back(elem->edge_id[2]);
			edges_ids.push_back(elem->edge_id[4]);
			edges_ids.push_back(elem->edge_id[7]);
			edges_ids.push_back(elem->edge_id[8]);
			edges_ids.push_back(elem->edge_id[10]);
		}else if(elem->pad==104){
			edges_ids.push_back(elem->edge_id[4]);
			edges_ids.push_back(elem->edge_id[5]);
			edges_ids.push_back(elem->edge_id[6]);
			edges_ids.push_back(elem->edge_id[7]);
			edges_ids.push_back(elem->edge_id[8]);
			edges_ids.push_back(elem->edge_id[10]);
		}else if(elem->pad==103){
			edges_ids.push_back(elem->edge_id[0]);
			edges_ids.push_back(elem->edge_id[2]);
			edges_ids.push_back(elem->edge_id[5]);
			edges_ids.push_back(elem->edge_id[6]);
			edges_ids.push_back(elem->edge_id[8]);
			edges_ids.push_back(elem->edge_id[10]);
		}else if(elem->pad==102){
			edges_ids.push_back(elem->edge_id[0]);
			edges_ids.push_back(elem->edge_id[2]);
			edges_ids.push_back(elem->edge_id[3]);
			edges_ids.push_back(elem->edge_id[8]);
			edges_ids.push_back(elem->edge_id[10]);
			edges_ids.push_back(elem->edge_id[11]);
		}else if(elem->pad==101){
			edges_ids.push_back(elem->edge_id[1]);
			edges_ids.push_back(elem->edge_id[2]);
			edges_ids.push_back(elem->edge_id[3]);
			edges_ids.push_back(elem->edge_id[9]);
			edges_ids.push_back(elem->edge_id[10]);
			edges_ids.push_back(elem->edge_id[11]);
		}else if(elem->pad==100){
			edges_ids.push_back(elem->edge_id[0]);
			edges_ids.push_back(elem->edge_id[1]);
			edges_ids.push_back(elem->edge_id[2]);
			edges_ids.push_back(elem->edge_id[8]);
			edges_ids.push_back(elem->edge_id[9]);
			edges_ids.push_back(elem->edge_id[10]);
		}else if(elem->pad==99){
			edges_ids.push_back(elem->edge_id[0]);
			edges_ids.push_back(elem->edge_id[1]);
			edges_ids.push_back(elem->edge_id[3]);
			edges_ids.push_back(elem->edge_id[8]);
			edges_ids.push_back(elem->edge_id[9]);
			edges_ids.push_back(elem->edge_id[11]);
		}
		//template 5
		else if(elem->pad==98){
			edges_ids.push_back(elem->edge_id[3]);
			edges_ids.push_back(elem->edge_id[7]);
		}else if(elem->pad==97){
			edges_ids.push_back(elem->edge_id[6]);
			edges_ids.push_back(elem->edge_id[7]);
		}else if(elem->pad==96){
			edges_ids.push_back(elem->edge_id[4]);
			edges_ids.push_back(elem->edge_id[6]);
		}else if(elem->pad==95){
			edges_ids.push_back(elem->edge_id[3]);
			edges_ids.push_back(elem->edge_id[4]);
		}else if(elem->pad==94){
			edges_ids.push_back(elem->edge_id[2]);
			edges_ids.push_back(elem->edge_id[7]);
		}else if(elem->pad==93){
			edges_ids.push_back(elem->edge_id[7]);
			edges_ids.push_back(elem->edge_id[10]);
		}else if(elem->pad==92){
			edges_ids.push_back(elem->edge_id[6]);
			edges_ids.push_back(elem->edge_id[10]);
		}else if(elem->pad==91){
			edges_ids.push_back(elem->edge_id[2]);
			edges_ids.push_back(elem->edge_id[6]);
		}else if(elem->pad==90){
			edges_ids.push_back(elem->edge_id[1]);
			edges_ids.push_back(elem->edge_id[5]);
		}else if(elem->pad==89){
			edges_ids.push_back(elem->edge_id[1]);
			edges_ids.push_back(elem->edge_id[6]);
		}else if(elem->pad==88){
			edges_ids.push_back(elem->edge_id[6]);
			edges_ids.push_back(elem->edge_id[9]);
		}else if(elem->pad==87){
			edges_ids.push_back(elem->edge_id[5]);
			edges_ids.push_back(elem->edge_id[9]);
		}else if(elem->pad==86){
			edges_ids.push_back(elem->edge_id[8]);
			edges_ids.push_back(elem->edge_id[11]);
		}else if(elem->pad==85){
			edges_ids.push_back(elem->edge_id[10]);
			edges_ids.push_back(elem->edge_id[11]);
		}else if(elem->pad==84){
			edges_ids.push_back(elem->edge_id[9]);
			edges_ids.push_back(elem->edge_id[10]);
		}else if(elem->pad==83){
			edges_ids.push_back(elem->edge_id[8]);
			edges_ids.push_back(elem->edge_id[9]);
		}else if(elem->pad==82){
			edges_ids.push_back(elem->edge_id[0]);
			edges_ids.push_back(elem->edge_id[3]);
		}else if(elem->pad==81){
			edges_ids.push_back(elem->edge_id[2]);
			edges_ids.push_back(elem->edge_id[3]);
		}else if(elem->pad==80){
			edges_ids.push_back(elem->edge_id[1]);
			edges_ids.push_back(elem->edge_id[2]);
		}else if(elem->pad==79){
			edges_ids.push_back(elem->edge_id[0]);
			edges_ids.push_back(elem->edge_id[1]);
		}else if(elem->pad==78){
			edges_ids.push_back(elem->edge_id[3]);
			edges_ids.push_back(elem->edge_id[7]);
			edges_ids.push_back(elem->edge_id[11]);
		}else if(elem->pad==77){
			edges_ids.push_back(elem->edge_id[4]);
			edges_ids.push_back(elem->edge_id[7]);
			edges_ids.push_back(elem->edge_id[11]);
		}else if(elem->pad==76){
			edges_ids.push_back(elem->edge_id[3]);
			edges_ids.push_back(elem->edge_id[4]);
			edges_ids.push_back(elem->edge_id[11]);
		}else if(elem->pad==75){
			edges_ids.push_back(elem->edge_id[3]);
			edges_ids.push_back(elem->edge_id[4]);
			edges_ids.push_back(elem->edge_id[7]);
		}else if(elem->pad==74){
			edges_ids.push_back(elem->edge_id[2]);
			edges_ids.push_back(elem->edge_id[7]);
			edges_ids.push_back(elem->edge_id[10]);
		}else if(elem->pad==73){
			edges_ids.push_back(elem->edge_id[6]);
			edges_ids.push_back(elem->edge_id[7]);
			edges_ids.push_back(elem->edge_id[10]);
		}else if(elem->pad==72){
			edges_ids.push_back(elem->edge_id[2]);
			edges_ids.push_back(elem->edge_id[6]);
			edges_ids.push_back(elem->edge_id[10]);
		}else if(elem->pad==71){
			edges_ids.push_back(elem->edge_id[2]);
			edges_ids.push_back(elem->edge_id[6]);
			edges_ids.push_back(elem->edge_id[7]);
		}else if(elem->pad==70){
			edges_ids.push_back(elem->edge_id[1]);
			edges_ids.push_back(elem->edge_id[6]);
			edges_ids.push_back(elem->edge_id[9]);
		}else if(elem->pad==69){
			edges_ids.push_back(elem->edge_id[5]);
			edges_ids.push_back(elem->edge_id[6]);
			edges_ids.push_back(elem->edge_id[9]);
		}else if(elem->pad==68){
			edges_ids.push_back(elem->edge_id[1]);
			edges_ids.push_back(elem->edge_id[5]);
			edges_ids.push_back(elem->edge_id[9]);
		}else if(elem->pad==67){
			edges_ids.push_back(elem->edge_id[1]);
			edges_ids.push_back(elem->edge_id[5]);
			edges_ids.push_back(elem->edge_id[6]);
		}else if(elem->pad==66){
			edges_ids.push_back(elem->edge_id[0]);
			edges_ids.push_back(elem->edge_id[4]);
			edges_ids.push_back(elem->edge_id[8]);
		}else if(elem->pad==65){
			edges_ids.push_back(elem->edge_id[4]);
			edges_ids.push_back(elem->edge_id[5]);
			edges_ids.push_back(elem->edge_id[8]);
		}else if(elem->pad==64){
			edges_ids.push_back(elem->edge_id[0]);
			edges_ids.push_back(elem->edge_id[5]);
			edges_ids.push_back(elem->edge_id[8]);
		}else if(elem->pad==63){
			edges_ids.push_back(elem->edge_id[0]);
			edges_ids.push_back(elem->edge_id[4]);
			edges_ids.push_back(elem->edge_id[5]);
		}else if(elem->pad==62){
			edges_ids.push_back(elem->edge_id[8]);
			edges_ids.push_back(elem->edge_id[10]);
			edges_ids.push_back(elem->edge_id[11]);
		}else if(elem->pad==61){
			edges_ids.push_back(elem->edge_id[9]);
			edges_ids.push_back(elem->edge_id[10]);
			edges_ids.push_back(elem->edge_id[11]);
		}else if(elem->pad==60){
			edges_ids.push_back(elem->edge_id[8]);
			edges_ids.push_back(elem->edge_id[9]);
			edges_ids.push_back(elem->edge_id[11]);
		}else if(elem->pad==59){
			edges_ids.push_back(elem->edge_id[8]);
			edges_ids.push_back(elem->edge_id[9]);
			edges_ids.push_back(elem->edge_id[10]);
		}else if(elem->pad==58){
			edges_ids.push_back(elem->edge_id[0]);
			edges_ids.push_back(elem->edge_id[2]);
			edges_ids.push_back(elem->edge_id[3]);
		}else if(elem->pad==57){
			edges_ids.push_back(elem->edge_id[1]);
			edges_ids.push_back(elem->edge_id[2]);
			edges_ids.push_back(elem->edge_id[3]);
		}else if(elem->pad==56){
			edges_ids.push_back(elem->edge_id[0]);
			edges_ids.push_back(elem->edge_id[1]);
			edges_ids.push_back(elem->edge_id[2]);
		}else if(elem->pad==55){
			edges_ids.push_back(elem->edge_id[0]);
			edges_ids.push_back(elem->edge_id[1]);
			edges_ids.push_back(elem->edge_id[3]);
		}else if(elem->pad==54){
			edges_ids.push_back(elem->edge_id[3]);
			edges_ids.push_back(elem->edge_id[4]);
			edges_ids.push_back(elem->edge_id[7]);
			edges_ids.push_back(elem->edge_id[11]);
		}else if(elem->pad==53){
			edges_ids.push_back(elem->edge_id[2]);
			edges_ids.push_back(elem->edge_id[6]);
			edges_ids.push_back(elem->edge_id[7]);
			edges_ids.push_back(elem->edge_id[10]);
		}else if(elem->pad==52){
			edges_ids.push_back(elem->edge_id[1]);
			edges_ids.push_back(elem->edge_id[5]);
			edges_ids.push_back(elem->edge_id[5]);
			edges_ids.push_back(elem->edge_id[9]);
		}else if(elem->pad==51){
			edges_ids.push_back(elem->edge_id[0]);
			edges_ids.push_back(elem->edge_id[4]);
			edges_ids.push_back(elem->edge_id[5]);
			edges_ids.push_back(elem->edge_id[8]);
		}else if(elem->pad==50){
			edges_ids.push_back(elem->edge_id[8]);
			edges_ids.push_back(elem->edge_id[9]);
			edges_ids.push_back(elem->edge_id[10]);
			edges_ids.push_back(elem->edge_id[11]);
		}else if(elem->pad==49){
			edges_ids.push_back(elem->edge_id[0]);
			edges_ids.push_back(elem->edge_id[1]);
			edges_ids.push_back(elem->edge_id[2]);
			edges_ids.push_back(elem->edge_id[3]);
		}
		//template 4
		else if(elem->pad==48){
			edges_ids.push_back(elem->edge_id[4]);
			edges_ids.push_back(elem->edge_id[5]);
			edges_ids.push_back(elem->edge_id[9]);
			edges_ids.push_back(elem->edge_id[11]);
		}else if(elem->pad==47){
			edges_ids.push_back(elem->edge_id[6]);
			edges_ids.push_back(elem->edge_id[7]);
			edges_ids.push_back(elem->edge_id[9]);
			edges_ids.push_back(elem->edge_id[11]);
		}else if(elem->pad==46){
			edges_ids.push_back(elem->edge_id[1]);
			edges_ids.push_back(elem->edge_id[3]);
			edges_ids.push_back(elem->edge_id[6]);
			edges_ids.push_back(elem->edge_id[7]);
		}else if(elem->pad==45){
			edges_ids.push_back(elem->edge_id[1]);
			edges_ids.push_back(elem->edge_id[3]);
			edges_ids.push_back(elem->edge_id[4]);
			edges_ids.push_back(elem->edge_id[5]);
		}else if(elem->pad==44){
			edges_ids.push_back(elem->edge_id[4]);
			edges_ids.push_back(elem->edge_id[7]);
			edges_ids.push_back(elem->edge_id[8]);
			edges_ids.push_back(elem->edge_id[10]);
		}else if(elem->pad==43){
			edges_ids.push_back(elem->edge_id[5]);
			edges_ids.push_back(elem->edge_id[6]);
			edges_ids.push_back(elem->edge_id[8]);
			edges_ids.push_back(elem->edge_id[10]);
		}else if(elem->pad==42){
			edges_ids.push_back(elem->edge_id[0]);
			edges_ids.push_back(elem->edge_id[2]);
			edges_ids.push_back(elem->edge_id[5]);
			edges_ids.push_back(elem->edge_id[6]);
		}else if(elem->pad==41){
			edges_ids.push_back(elem->edge_id[0]);
			edges_ids.push_back(elem->edge_id[2]);
			edges_ids.push_back(elem->edge_id[4]);
			edges_ids.push_back(elem->edge_id[7]);
		}else if(elem->pad==40){
			edges_ids.push_back(elem->edge_id[0]);
			edges_ids.push_back(elem->edge_id[3]);
			edges_ids.push_back(elem->edge_id[8]);
			edges_ids.push_back(elem->edge_id[11]);
		}else if(elem->pad==39){
			edges_ids.push_back(elem->edge_id[2]);
			edges_ids.push_back(elem->edge_id[3]);
			edges_ids.push_back(elem->edge_id[10]);
			edges_ids.push_back(elem->edge_id[11]);
		}else if(elem->pad==38){
			edges_ids.push_back(elem->edge_id[1]);
			edges_ids.push_back(elem->edge_id[2]);
			edges_ids.push_back(elem->edge_id[9]);
			edges_ids.push_back(elem->edge_id[10]);
		}else if(elem->pad==37){
			edges_ids.push_back(elem->edge_id[0]);
			edges_ids.push_back(elem->edge_id[1]);
			edges_ids.push_back(elem->edge_id[8]);
			edges_ids.push_back(elem->edge_id[9]);
		}
		//template 3
		else if(elem->pad==36){
			edges_ids.push_back(elem->edge_id[3]);
			edges_ids.push_back(elem->edge_id[11]);
		}else if(elem->pad==35){
			edges_ids.push_back(elem->edge_id[1]);
			edges_ids.push_back(elem->edge_id[10]);
		}else if(elem->pad==34){
			edges_ids.push_back(elem->edge_id[1]);
			edges_ids.push_back(elem->edge_id[9]);
		}else if(elem->pad==33){
			edges_ids.push_back(elem->edge_id[0]);
			edges_ids.push_back(elem->edge_id[8]);
		}else if(elem->pad==32){
			edges_ids.push_back(elem->edge_id[4]);
			edges_ids.push_back(elem->edge_id[7]);
		}else if(elem->pad==31){
			edges_ids.push_back(elem->edge_id[6]);
			edges_ids.push_back(elem->edge_id[7]);
		}else if(elem->pad==30){
			edges_ids.push_back(elem->edge_id[5]);
			edges_ids.push_back(elem->edge_id[6]);
		}else if(elem->pad==29){
			edges_ids.push_back(elem->edge_id[4]);
			edges_ids.push_back(elem->edge_id[5]);
		}else if(elem->pad==28){
			edges_ids.push_back(elem->edge_id[9]);
			edges_ids.push_back(elem->edge_id[11]);
		}else if(elem->pad==27){
			edges_ids.push_back(elem->edge_id[8]);
			edges_ids.push_back(elem->edge_id[10]);
		}else if(elem->pad==26){
			edges_ids.push_back(elem->edge_id[1]);
			edges_ids.push_back(elem->edge_id[3]);
		}else if(elem->pad==25){
			edges_ids.push_back(elem->edge_id[0]);
			edges_ids.push_back(elem->edge_id[2]);
		}
		//template 2
		else if(elem->pad==24){
			edges_ids.push_back(elem->edge_id[4]);
			edges_ids.push_back(elem->edge_id[5]);
			edges_ids.push_back(elem->edge_id[6]);
			edges_ids.push_back(elem->edge_id[7]);
		}else if(elem->pad==23){
			edges_ids.push_back(elem->edge_id[1]);
			edges_ids.push_back(elem->edge_id[3]);
			edges_ids.push_back(elem->edge_id[9]);
			edges_ids.push_back(elem->edge_id[11]);
		}else if(elem->pad==22){
			edges_ids.push_back(elem->edge_id[0]);
			edges_ids.push_back(elem->edge_id[2]);
			edges_ids.push_back(elem->edge_id[8]);
			edges_ids.push_back(elem->edge_id[10]);
		}
		//template 1
		else if(elem->pad==21){
			edges_ids.push_back(elem->edge_id[11]);
		}else if(elem->pad==20){
			edges_ids.push_back(elem->edge_id[10]);
		}else if(elem->pad==19){
			edges_ids.push_back(elem->edge_id[9]);
		}else if(elem->pad==18){
			edges_ids.push_back(elem->edge_id[8]);
		}else if(elem->pad==17){
			edges_ids.push_back(elem->edge_id[7]);
		}else if(elem->pad==16){
			edges_ids.push_back(elem->edge_id[6]);
		}else if(elem->pad==15){
			edges_ids.push_back(elem->edge_id[5]);
		}else if(elem->pad==14){
			edges_ids.push_back(elem->edge_id[4]);
		}else if(elem->pad==13){
			edges_ids.push_back(elem->edge_id[3]);
		}else if(elem->pad==12){
			edges_ids.push_back(elem->edge_id[2]);
		}else if(elem->pad==11){
			edges_ids.push_back(elem->edge_id[1]);
		}else if(elem->pad==10){
			edges_ids.push_back(elem->edge_id[0]);
		}
	}

	//TODO check the time of the two implementations
	printf(" Edges intercepted: %d\n", edges_ids.size());

	//remove duplicates edges
	//std::set<int> s;
	//unsigned size = edges_ids.size();
	//for( unsigned i = 0; i < size; ++i ) s.insert( edges_ids[i] );
	//edges_ids.assign( s.begin(), s.end() );
	//printf("Cleaned Edges intercepted: %d\n", edges_ids.size());

	std::sort( edges_ids.begin(), edges_ids.end() );
	edges_ids.erase( std::unique( edges_ids.begin(), edges_ids.end() ), edges_ids.end() );
	printf("Cleaned Edges intercepted: %d\n", edges_ids.size());

}

void Edge_propagation(hexa_tree_t* mesh, std::vector<int>& elements_ids, std::vector<int>& edges_ids) {

	for(int iel= 0; iel < mesh->total_n_elements; iel++ ){
		octant_t *elem = (octant_t*) sc_array_index(&mesh->elements, iel);

		for (int edge = 0; edge < 12; ++edge) {
			int node1 = elem->nodes[EdgeVerticesMap[edge][0]].id;
			int node2 = elem->nodes[EdgeVerticesMap[edge][1]].id;
			elem->edge_id[edge] = edge_hash_function(node1, node2);
			//elem->edge_ref[edge]=false;
			//printf("edge numero: %d, elemento numero: %d id da edge: %d\n",edge, iel,elem->edge_id[edge]);


			if(binary_search(edges_ids.begin(), edges_ids.end(), elem->edge_id[edge])){
				//printf("porra\n");
				//if(binary_search(elements_ids.begin(), elements_ids.end(), iel)){

				//}else{
					//printf("edge numero: %d, elemento numero: %d id da edge: %d\n",edge, iel,elem->edge_id[edge]);
					elem->edge_ref[edge]=true;
					elem->pad = -1;
					elements_ids.push_back(iel);
				//}
			}

			//for(int i1 = 0; i1 < edges_ids.size(); i1++){
				//for(int i2 = 0; i2 < elements_ids.size(); i2++){
					//if(iel==elements_ids[i2]){
						//break;
				//	}else{
					//	if(elem->edge_id[edge]==edges_ids[i1]){
							//printf("\n");
						//	printf("edge numero: %d, elemento numero: %d id da edge: %d\n",edge, iel,elem->edge_id[edge]);
							//printf("\n");
							//elements_ids.push_back(iel);
							//elem->edge_ref[edge]=true;
							//elem->pad = -1;
						//}
					//}
				//}
			//}

			// if(binary_search(vector.begin(), vector.end(), item)){
			// Found the item
		}

		/*
		if(std::find(edges_ids.begin(), edges_ids.end(), elem->edge_id[edge])!=edges_ids.end()){
			if(std::find(elements_ids.begin(), elements_ids.end(), iel)!=elements_ids.end()){

			}else{
				printf("\n");
				printf("edge numero: %d, elemento numero: %d id da edge: %d\n",edge, iel,elem->edge_id[edge]);
				printf("\n");

				elements_ids.push_back(iel);
				elem->edge_ref[edge]=true;
				elem->pad = -1;
			}
		}

	}
		 */
	}


	//for(int iel= 0; iel < mesh->total_n_elements; iel++ ){
	//	octant_t *elem = (octant_t*) sc_array_index(&mesh->elements, iel);

		//for (int edge = 0; edge < 12; ++edge) {

			//if (elem->edge_ref[edge] ){
				//elements_ids.push_back(iel);
			//}
		//}
	//}

	//cleaning the element vector
	std::sort( elements_ids.begin(), elements_ids.end() );
	elements_ids.erase( std::unique( elements_ids.begin(), elements_ids.end() ), elements_ids.end() );
}



void CheckTemplate(hexa_tree_t* mesh, const std::vector<double>& coords, std::vector<int>& elements_ids, bool flag) {

	std::vector<int> edges_ids;

	for (int iel = 0; iel < elements_ids.size(); ++iel) {

		octant_t *elem = (octant_t*) sc_array_index(&mesh->elements, elements_ids[iel]);

		int Edge2GNode[12][2]={0};
		GtsSegment * segments[12]={0};
		GtsPoint * point[12]={NULL};
		int ed_cont = 0;

		for (int edge = 0; edge < 12; ++edge) {
			point[edge] = NULL;
			elem->edge_ref[edge] = false;
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
					elem->edge_ref[edge] = true;
					ed_cont++;
					break;
				}
				list = list->next;
			}
		}

		//Bounding box intercepted
		if(elem->pad == -1 && ed_cont == 0){
			elements_ids.erase(elements_ids.begin() + iel);
			elem->pad = -2;
		}

		//clean points
		for (int edge = 0; edge < 12; edge++) {
			if (point[edge]) gts_object_destroy(GTS_OBJECT(point[edge]));
			point[edge] = NULL;
		}
	}

	for (int i = 0; i < 2; i++){
		printf("numero %d\n",i);
		printf(" Elements ref: %d\n", elements_ids.size());
		IdentifyTemplate(mesh, coords, elements_ids);
		Edge_identification( mesh, elements_ids, edges_ids);
		Edge_propagation( mesh, elements_ids, edges_ids);
	}
	printf(" Elements ref: %d\n", elements_ids.size());
	IdentifyTemplate(mesh, coords, elements_ids);

}

void Material_apply(hexa_tree_t *mesh, std::vector<double>& coords, std::vector<int>& element_ids, const char* surface_bathy){

	sc_array_t *elements = &mesh->elements;

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

void SideEL(hexa_tree_t *mesh,int x, int y, int z ,std::vector<int>& element_ids_local){

	sc_array_t *elements = &mesh->elements;

	for(int iel = 0; iel < elements->elem_count; ++iel) {

		octant_t *elem = (octant_t*) sc_array_index(&mesh->elements, iel);

		if(z==0){
			if((elem->x==x-1)&&(elem->y==y-1)&&(elem->z==z+1)){
				element_ids_local.push_back(iel);
			}else if((elem->x==x-1)&&(elem->y==y)&&(elem->z==z+1)) {
				element_ids_local.push_back(iel);
			}else if((elem->x==x-1)&&(elem->y==y+1)&&(elem->z==z+1)) {
				element_ids_local.push_back(iel);
			}else if((elem->x==x)&&(elem->y==y-1)&&(elem->z==z+1)) {
				element_ids_local.push_back(iel);
			}else if((elem->x==x)&&(elem->y==y)&&(elem->z==z+1)) {
				element_ids_local.push_back(iel);
			}else if((elem->x==x)&&(elem->y==y+1)&&(elem->z==z+1)) {
				element_ids_local.push_back(iel);
			}else if((elem->x==x+1)&&(elem->y==y-1)&&(elem->z==z+1)) {
				element_ids_local.push_back(iel);
			}else if((elem->x==x+1)&&(elem->y==y)&&(elem->z==z+1)) {
				element_ids_local.push_back(iel);
			}else if((elem->x==x+1)&&(elem->y==y+1)&&(elem->z==z+1)) {
				element_ids_local.push_back(iel);
			}else if((elem->x==x-1)&&(elem->y==y-1)&&(elem->z==z)){
				element_ids_local.push_back(iel);
			}else if((elem->x==x-1)&&(elem->y==y)&&(elem->z==z)) {
				element_ids_local.push_back(iel);
			}else if((elem->x==x-1)&&(elem->y==y+1)&&(elem->z==z)) {
				element_ids_local.push_back(iel);
			}else if((elem->x==x)&&(elem->y==y-1)&&(elem->z==z)) {
				element_ids_local.push_back(iel);
			}else if((elem->x==x)&&(elem->y==y)&&(elem->z==z)) {
				element_ids_local.push_back(iel);
			}else if((elem->x==x)&&(elem->y==y+1)&&(elem->z==z)) {
				element_ids_local.push_back(iel);
			}else if((elem->x==x+1)&&(elem->y==y-1)&&(elem->z==z)) {
				element_ids_local.push_back(iel);
			}else if((elem->x==x+1)&&(elem->y==y)&&(elem->z==z)) {
				element_ids_local.push_back(iel);
			}else if((elem->x==x+1)&&(elem->y==y+1)&&(elem->z==z)) {
				element_ids_local.push_back(iel);
			}
		}else{
			if((elem->x==x-1)&&(elem->y==y-1)&&(elem->z==z+1)){
				element_ids_local.push_back(iel);
			}else if((elem->x==x-1)&&(elem->y==y)&&(elem->z==z+1)) {
				element_ids_local.push_back(iel);
			}else if((elem->x==x-1)&&(elem->y==y+1)&&(elem->z==z+1)) {
				element_ids_local.push_back(iel);
			}else if((elem->x==x)&&(elem->y==y-1)&&(elem->z==z+1)) {
				element_ids_local.push_back(iel);
			}else if((elem->x==x)&&(elem->y==y)&&(elem->z==z+1)) {
				element_ids_local.push_back(iel);
			}else if((elem->x==x)&&(elem->y==y+1)&&(elem->z==z+1)) {
				element_ids_local.push_back(iel);
			}else if((elem->x==x+1)&&(elem->y==y-1)&&(elem->z==z+1)) {
				element_ids_local.push_back(iel);
			}else if((elem->x==x+1)&&(elem->y==y)&&(elem->z==z+1)) {
				element_ids_local.push_back(iel);
			}else if((elem->x==x+1)&&(elem->y==y+1)&&(elem->z==z+1)) {
				element_ids_local.push_back(iel);
			}else if((elem->x==x-1)&&(elem->y==y-1)&&(elem->z==z)){
				element_ids_local.push_back(iel);
			}else if((elem->x==x-1)&&(elem->y==y)&&(elem->z==z)) {
				element_ids_local.push_back(iel);
			}else if((elem->x==x-1)&&(elem->y==y+1)&&(elem->z==z)) {
				element_ids_local.push_back(iel);
			}else if((elem->x==x)&&(elem->y==y-1)&&(elem->z==z)) {
				element_ids_local.push_back(iel);
			}else if((elem->x==x)&&(elem->y==y)&&(elem->z==z)) {
				element_ids_local.push_back(iel);
			}else if((elem->x==x)&&(elem->y==y+1)&&(elem->z==z)) {
				element_ids_local.push_back(iel);
			}else if((elem->x==x+1)&&(elem->y==y-1)&&(elem->z==z)) {
				element_ids_local.push_back(iel);
			}else if((elem->x==x+1)&&(elem->y==y)&&(elem->z==z)) {
				element_ids_local.push_back(iel);
			}else if((elem->x==x+1)&&(elem->y==y+1)&&(elem->z==z)) {
				element_ids_local.push_back(iel);
			}else if((elem->x==x-1)&&(elem->y==y-1)&&(elem->z==z-1)){
				element_ids_local.push_back(iel);
			}else if((elem->x==x-1)&&(elem->y==y)&&(elem->z==z-1)) {
				element_ids_local.push_back(iel);
			}else if((elem->x==x-1)&&(elem->y==y+1)&&(elem->z==z-1)) {
				element_ids_local.push_back(iel);
			}else if((elem->x==x)&&(elem->y==y-1)&&(elem->z==z-1)) {
				element_ids_local.push_back(iel);
			}else if((elem->x==x)&&(elem->y==y)&&(elem->z==z-1)) {
				element_ids_local.push_back(iel);
			}else if((elem->x==x)&&(elem->y==y+1)&&(elem->z==z-1)) {
				element_ids_local.push_back(iel);
			}else if((elem->x==x+1)&&(elem->y==y-1)&&(elem->z==z-1)) {
				element_ids_local.push_back(iel);
			}else if((elem->x==x+1)&&(elem->y==y)&&(elem->z==z-1)) {
				element_ids_local.push_back(iel);
			}else if((elem->x==x+1)&&(elem->y==y+1)&&(elem->z==z-1)) {
				element_ids_local.push_back(iel);
			}
		}
	}
}

void template2Rand(hexa_tree_t *mesh, std::vector<double>& coords,octant_t * elem, double v){

	double dx = coords[3*elem->nodes[1].id  ] - coords[3*elem->nodes[0].id  ] ;
	double dy = coords[3*elem->nodes[2].id+1] - coords[3*elem->nodes[0].id+1] ;
	double dz = coords[3*elem->nodes[0].id+2] - coords[3*elem->nodes[4].id+2] ;

	dx = fabs(dx);
	dy = fabs(dy);
	dz = fabs(dz);

	for(int i = 0;i<8;i++){
		int v1 = (-50 + rand() % 100)*0.02;
		int v2 = (-50 + rand() % 100)*0.02;
		int v3 = (-50 + rand() % 100)*0.02;
		int node_change = 3*elem->nodes[i].id;
		coords[node_change  ] = coords[node_change  ] + v*v1*dx;
		coords[node_change+1] = coords[node_change+1] + v*v2*dy;
		if(elem->z>0){
			coords[node_change+2] = coords[node_change+2] + v*v3*dz;
		}else if(elem->z==0 && i>=4){
			coords[node_change+2] = coords[node_change+2] + v*v3*dz;
		}
	}
}

void ChangeTemplate(hexa_tree_t* mesh, std::vector<double>& coords, std::vector<int>& elements_ids) {
	std::vector<int> element_ids_local;
	double coords_orgi[24];
	double v_1 = 1;
	double v_2 = 1;
	int iter_max = 200;
	int ref=0;

#if 0
	double tt=1.02;


	FILE * fdbg;
	fdbg = fopen("Change_template.txt", "w");

	int nao_sei = 0;
	int n_iter=0;
	int var_aux[27];
	int var_aux_1[27];
	int el_4 = 0;

	srand (time(NULL));

	sc_array_t *elements = &mesh->elements;

	for (int iel = 0; iel < elements_ids.size(); ++iel) {
		octant_t *elem = (octant_t*) sc_array_index(&mesh->elements, elements_ids[iel]);
		v_1 = 1;
		//v_2 = 1;

		if(elem->pad==4){

			//printf("El:%d\n",elements_ids[iel]);
			element_ids_local.clear();
			SideEL(mesh, elem->x, elem->y, elem->z ,element_ids_local);

			fprintf(fdbg,"El: %d Case type 4\n",elements_ids[iel]);
			fprintf(fdbg,"El:");
			for(int co = 0; co<element_ids_local.size();co++){
				fprintf(fdbg,"%d ",element_ids_local[co]);
				if(element_ids_local[co]==elements_ids[iel]){ref=co;}
			}
			fprintf(fdbg,"\n");

			fprintf(fdbg,"REF:%d \n",ref);

			fprintf(fdbg,"Pad Original: ");
			for(int c = 0; c < element_ids_local.size(); c++){
				octant_t *h = (octant_t*) sc_array_index(&mesh->elements, element_ids_local[c]);
				var_aux[c]=h->pad;
				fprintf(fdbg,"%d ",var_aux[c]);
			}
			fprintf(fdbg,"\n");

			for(int n_nodes = 0; n_nodes<8;n_nodes++){
				octant_t *elem = (octant_t*) sc_array_index(&mesh->elements, elements_ids[iel]);
				int node = 3*elem->nodes[n_nodes].id;
				coords_orgi[3*n_nodes  ] = coords[node  ];
				coords_orgi[3*n_nodes+1] = coords[node+1];
				coords_orgi[3*n_nodes+2] = coords[node+2];
			}

			n_iter = 0;
			bool flag1 = true;
			bool flag2 = true;

			while(n_iter<iter_max && flag1){

				// move one node... try to change to template 1
#if 0
				template2Rand(mesh, coords, elem, v_1);
				CheckTemplate(mesh, coords, element_ids_local,false);
				octant_t *elem = (octant_t*) sc_array_index(&mesh->elements, elements_ids[iel]);

				fprintf(fdbg,"Pad local:    ");
				for(int co = 0; co < element_ids_local.size(); co++){
					octant_t *h = (octant_t*) sc_array_index(&mesh->elements, element_ids_local[co]);
					var_aux_1[co]=h->pad;
					fprintf(fdbg,"%d ",h->pad);
				}
				fprintf(fdbg,"\n");

				if(elem->pad!=4 && elem->pad!=-10){
					flag2=false;
					for(int co = 0; co < element_ids_local.size(); co++){
						if(element_ids_local[co]!=elements_ids[iel]){
							if(var_aux[co]==var_aux_1[co]){
								flag2=false;
							}else{
								if(var_aux_1[co]==4){
									flag2=true;
									break;
								}else if(var_aux_1[co]==-10){
									flag2=true;
									break;
								}
							}
						}
					}
				}

				if(!flag2){
					flag1=false;
					flag2=false;
					printf("To livreeeee rand El:%d, n_iter:%d\n",elements_ids[iel],n_iter);
					break;
				}

				//if(elem->pad==4){
				//	CuttedEdges (mesh, coords, elements_ids, iel, &point[12], &edge_list[12]);

				GtsSegment * segments[12]={0};
				int Edge2GNode[12][2];
				GtsPoint * point[12] = {NULL};
				bool edge_list[12] = {false};

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

					edge_list[edge] = false;

					if (list == NULL) continue;
					while (list) {
						GtsBBox *b = GTS_BBOX(list->data);
						point[edge] = SegmentTriangleIntersection(segments[edge], GTS_TRIANGLE(b->bounded));
						if (point[edge]) {
							edge_list[edge] = true;
							break;
						}
						list = list->next;
					}
				}

				double d_c1=0;
				double d_c2=0;

				if(edge_list[0]&&edge_list[1]){
					int node_change = 3*elem->nodes[1].id;

					GtsPoint *p0 = gts_point_new(gts_point_class(),
							coords[node_change],
							coords[node_change + 1],
							coords[node_change + 2]);

					d_c1 = gts_point_distance(point[0], p0);
					d_c2 = gts_point_distance(point[1], p0);

					coords[node_change  ] = coords[node_change  ] - d_c1*tt;
					coords[node_change+1] = coords[node_change+1] + d_c2*tt;
					//move 1

				}else if(edge_list[1]&&edge_list[2]){
					int node_change = 3*elem->nodes[2].id;

					GtsPoint *p0 = gts_point_new(gts_point_class(),
							coords[node_change],
							coords[node_change + 1],
							coords[node_change + 2]);

					d_c1 = gts_point_distance(point[2], p0);
					d_c2 = gts_point_distance(point[1], p0);

					coords[node_change  ] = coords[node_change  ] - d_c1*tt;
					coords[node_change+1] = coords[node_change+1] - d_c2*tt;
					//move 2

				}else if(edge_list[2]&&edge_list[3]){
					int node_change = 3*elem->nodes[3].id;

					GtsPoint *p0 = gts_point_new(gts_point_class(),
							coords[node_change],
							coords[node_change + 1],
							coords[node_change + 2]);

					d_c1 = gts_point_distance(point[2], p0);
					d_c2 = gts_point_distance(point[3], p0);

					coords[node_change  ] = coords[node_change  ] + d_c1*tt;
					coords[node_change+1] = coords[node_change+1] - d_c2*tt;
					//move 3

				}else if(edge_list[3]&&edge_list[0]){
					int node_change = 3*elem->nodes[0].id;

					GtsPoint *p0 = gts_point_new(gts_point_class(),
							coords[node_change],
							coords[node_change + 1],
							coords[node_change + 2]);

					d_c1 = gts_point_distance(point[0], p0);
					d_c2 = gts_point_distance(point[3], p0);

					coords[node_change  ] = coords[node_change  ] + d_c1*tt;
					coords[node_change+1] = coords[node_change+1] + d_c2*tt;
					//move 0

				}else if(edge_list[8]&&edge_list[9]){
					int node_change = 3*elem->nodes[5].id;

					GtsPoint *p0 = gts_point_new(gts_point_class(),
							coords[node_change],
							coords[node_change + 1],
							coords[node_change + 2]);

					d_c1 = gts_point_distance(point[8], p0);
					d_c2 = gts_point_distance(point[9], p0);

					coords[node_change  ] = coords[node_change  ] - d_c1*tt;
					coords[node_change+1] = coords[node_change+1] + d_c2*tt;
					//move 5

				}else if(edge_list[9]&&edge_list[10]){
					int node_change = 3*elem->nodes[6].id;

					GtsPoint *p0 = gts_point_new(gts_point_class(),
							coords[node_change],
							coords[node_change + 1],
							coords[node_change + 2]);

					d_c1 = gts_point_distance(point[10], p0);
					d_c2 = gts_point_distance(point[9], p0);

					coords[node_change  ] = coords[node_change  ] - d_c1*tt;
					coords[node_change+1] = coords[node_change+1] - d_c2*tt;
					//move 6

				}else if(edge_list[10]&&edge_list[11]){
					int node_change = 3*elem->nodes[7].id;

					GtsPoint *p0 = gts_point_new(gts_point_class(),
							coords[node_change],
							coords[node_change + 1],
							coords[node_change + 2]);

					d_c1 = gts_point_distance(point[10], p0);
					d_c2 = gts_point_distance(point[11], p0);

					coords[node_change  ] = coords[node_change  ] + d_c1*tt;
					coords[node_change+1] = coords[node_change+1] - d_c2*tt;
					//move 7

				}else if(edge_list[11]&&edge_list[8]){
					int node_change = 3*elem->nodes[4].id;

					GtsPoint *p0 = gts_point_new(gts_point_class(),
							coords[node_change],
							coords[node_change + 1],
							coords[node_change + 2]);

					d_c1 = gts_point_distance(point[8], p0);
					d_c2 = gts_point_distance(point[11], p0);

					coords[node_change  ] = coords[node_change  ] + d_c1*tt;
					coords[node_change+1] = coords[node_change+1] + d_c2*tt;
					//move 4

				}else if( edge_list[0]&&edge_list[4] ){
					int node_change = 3*elem->nodes[0].id;

					GtsPoint *p0 = gts_point_new(gts_point_class(),
							coords[node_change],
							coords[node_change + 1],
							coords[node_change + 2]);

					d_c1 = gts_point_distance(point[0], p0);
					d_c2 = gts_point_distance(point[4], p0);

					coords[node_change  ] = coords[node_change  ] + d_c1*tt;
					coords[node_change+2] = coords[node_change+2] - d_c2*tt;
					//move 0

				}else if( edge_list[3]&&edge_list[4] ){
					int node_change = 3*elem->nodes[0].id;

					GtsPoint *p0 = gts_point_new(gts_point_class(),
							coords[node_change],
							coords[node_change + 1],
							coords[node_change + 2]);

					d_c1 = gts_point_distance(point[3], p0);
					d_c2 = gts_point_distance(point[4], p0);

					coords[node_change+1] = coords[node_change+1] + d_c1*tt;
					coords[node_change+2] = coords[node_change+2] - d_c2*tt;
					//move 0

				}else if(edge_list[0]&&edge_list[5]){
					int node_change = 3*elem->nodes[1].id;

					GtsPoint *p0 = gts_point_new(gts_point_class(),
							coords[node_change],
							coords[node_change + 1],
							coords[node_change + 2]);

					d_c1 = gts_point_distance(point[0], p0);
					d_c2 = gts_point_distance(point[5], p0);

					coords[node_change  ] = coords[node_change  ] - d_c1*tt;
					coords[node_change+2] = coords[node_change+2] - d_c2*tt;
					//move 1

				}else if(edge_list[1]&&edge_list[5]){
					int node_change = 3*elem->nodes[1].id;

					GtsPoint *p0 = gts_point_new(gts_point_class(),
							coords[node_change],
							coords[node_change + 1],
							coords[node_change + 2]);

					d_c1 = gts_point_distance(point[1], p0);
					d_c2 = gts_point_distance(point[5], p0);

					coords[node_change+1] = coords[node_change+1] + d_c1*tt;
					coords[node_change+2] = coords[node_change+2] - d_c2*tt;
					//move 1

				}else if(edge_list[1]&&edge_list[6]){
					int node_change = 3*elem->nodes[2].id;

					GtsPoint *p0 = gts_point_new(gts_point_class(),
							coords[node_change],
							coords[node_change + 1],
							coords[node_change + 2]);

					d_c1 = gts_point_distance(point[1], p0);
					d_c2 = gts_point_distance(point[6], p0);

					coords[node_change+1] = coords[node_change+1] - d_c1*tt;
					coords[node_change+2] = coords[node_change+2] - d_c2*tt;
					//move 2

				}else if(edge_list[2]&&edge_list[6]){
					int node_change = 3*elem->nodes[2].id;

					GtsPoint *p0 = gts_point_new(gts_point_class(),
							coords[node_change],
							coords[node_change + 1],
							coords[node_change + 2]);

					d_c1 = gts_point_distance(point[2], p0);
					d_c2 = gts_point_distance(point[6], p0);

					coords[node_change  ] = coords[node_change  ] - d_c1*tt;
					coords[node_change+2] = coords[node_change+2] - d_c2*tt;
					//move 2

				}else if(edge_list[2]&&edge_list[7]){
					int node_change = 3*elem->nodes[3].id;

					GtsPoint *p0 = gts_point_new(gts_point_class(),
							coords[node_change],
							coords[node_change + 1],
							coords[node_change + 2]);

					d_c1 = gts_point_distance(point[2], p0);
					d_c2 = gts_point_distance(point[7], p0);

					coords[node_change  ] = coords[node_change  ] + d_c1*tt;
					coords[node_change+2] = coords[node_change+2] - d_c2*tt;
					//move 3

				}else if(edge_list[3]&&edge_list[7]){
					int node_change = 3*elem->nodes[3].id;

					GtsPoint *p0 = gts_point_new(gts_point_class(),
							coords[node_change],
							coords[node_change + 1],
							coords[node_change + 2]);

					d_c1 = gts_point_distance(point[3], p0);
					d_c2 = gts_point_distance(point[7], p0);

					coords[node_change+1] = coords[node_change+1] - d_c1*tt;
					coords[node_change+2] = coords[node_change+2] - d_c2*tt;
					//move 3

				}else if(edge_list[8]&&edge_list[4]){
					int node_change = 3*elem->nodes[4].id;

					GtsPoint *p0 = gts_point_new(gts_point_class(),
							coords[node_change],
							coords[node_change + 1],
							coords[node_change + 2]);

					d_c1 = gts_point_distance(point[8], p0);
					d_c2 = gts_point_distance(point[4], p0);

					coords[node_change  ] = coords[node_change  ] + d_c1*tt;
					coords[node_change+2] = coords[node_change+2] + d_c2*tt;
					//move 4

				}else if(edge_list[11]&&edge_list[4]){
					int node_change = 3*elem->nodes[4].id;

					GtsPoint *p0 = gts_point_new(gts_point_class(),
							coords[node_change],
							coords[node_change + 1],
							coords[node_change + 2]);

					d_c1 = gts_point_distance(point[11], p0);
					d_c2 = gts_point_distance(point[4], p0);

					coords[node_change+1] = coords[node_change+1] + d_c1*tt;
					coords[node_change+2] = coords[node_change+2] + d_c2*tt;
					//move 4

				}else if(edge_list[8]&&edge_list[5]){
					int node_change = 3*elem->nodes[5].id;

					GtsPoint *p0 = gts_point_new(gts_point_class(),
							coords[node_change],
							coords[node_change + 1],
							coords[node_change + 2]);

					d_c1 = gts_point_distance(point[8], p0);
					d_c2 = gts_point_distance(point[5], p0);

					coords[node_change  ] = coords[node_change  ] - d_c1*tt;
					coords[node_change+2] = coords[node_change+2] + d_c2*tt;
					//move 5

				}else if(edge_list[9]&&edge_list[5]){
					int node_change = 3*elem->nodes[5].id;

					GtsPoint *p0 = gts_point_new(gts_point_class(),
							coords[node_change],
							coords[node_change + 1],
							coords[node_change + 2]);

					d_c1 = gts_point_distance(point[9], p0);
					d_c2 = gts_point_distance(point[5], p0);

					coords[node_change+1] = coords[node_change+1] + d_c1*tt;
					coords[node_change+2] = coords[node_change+2] + d_c2*tt;
					//move 5

				}else if(edge_list[9]&&edge_list[6]){
					int node_change = 3*elem->nodes[6].id;

					GtsPoint *p0 = gts_point_new(gts_point_class(),
							coords[node_change],
							coords[node_change + 1],
							coords[node_change + 2]);

					d_c1 = gts_point_distance(point[9], p0);
					d_c2 = gts_point_distance(point[6], p0);

					coords[node_change+1] = coords[node_change+1] - d_c1*tt;
					coords[node_change+2] = coords[node_change+2] + d_c2*tt;
					//move 6

				}else if(edge_list[10]&&edge_list[6]){
					int node_change = 3*elem->nodes[6].id;

					GtsPoint *p0 = gts_point_new(gts_point_class(),
							coords[node_change],
							coords[node_change + 1],
							coords[node_change + 2]);

					d_c1 = gts_point_distance(point[10], p0);
					d_c2 = gts_point_distance(point[6], p0);

					coords[node_change  ] = coords[node_change  ] - d_c1*tt;
					coords[node_change+2] = coords[node_change+2] + d_c2*tt;
					//move 6

				}else if(edge_list[10]&&edge_list[7]){
					int node_change = 3*elem->nodes[7].id;

					GtsPoint *p0 = gts_point_new(gts_point_class(),
							coords[node_change],
							coords[node_change + 1],
							coords[node_change + 2]);

					d_c1 = gts_point_distance(point[10], p0);
					d_c2 = gts_point_distance(point[7], p0);

					coords[node_change  ] = coords[node_change  ] + d_c1*tt;
					coords[node_change+2] = coords[node_change+2] + d_c2*tt;
					//move 7

				}else if(edge_list[11]&&edge_list[7]){
					int node_change = 3*elem->nodes[7].id;

					GtsPoint *p0 = gts_point_new(gts_point_class(),
							coords[node_change],
							coords[node_change + 1],
							coords[node_change + 2]);

					d_c1 = gts_point_distance(point[11], p0);
					d_c2 = gts_point_distance(point[7], p0);

					coords[node_change+1] = coords[node_change+1] - d_c1*tt;
					coords[node_change+2] = coords[node_change+2] + d_c2*tt;
					//move 7
				}

				CheckTemplate(mesh, coords, element_ids_local,false);
				elem = (octant_t*) sc_array_index(&mesh->elements, elements_ids[iel]);

				fprintf(fdbg,"Pad local:    ");
				for(int co = 0; co < element_ids_local.size(); co++){
					octant_t *h = (octant_t*) sc_array_index(&mesh->elements, element_ids_local[co]);
					var_aux_1[co]=h->pad;
					fprintf(fdbg,"%d ",h->pad);
				}
				fprintf(fdbg,"\n");


				if(elem->pad!=4 && elem->pad!=-10){
					flag2=false;
					for(int co = 0; co < element_ids_local.size(); co++){
						if(element_ids_local[co]!=elements_ids[iel]){
							if(var_aux[co]==var_aux_1[co]){
								flag2=false;
							}else{
								if(var_aux_1[co]==4){
									flag2=true;
									break;
								}else if(var_aux_1[co]==-10){
									flag2=true;
									break;
								}
							}
						}
					}
				}
				//}

				if(!flag2){
					flag1=false;
					flag2=false;
					printf("To livreeeee mov 1 El:%d, n_iter:%d\n",elements_ids[iel],n_iter);
				}

				if( flag2 ){
					for(int n_nodes = 0; n_nodes<8;n_nodes++){
						int node = 3*elem->nodes[n_nodes].id;
						coords[node]   = coords_orgi[3*n_nodes];
						coords[node+1] = coords_orgi[3*n_nodes+1];
						coords[node+2] = coords_orgi[3*n_nodes+2];
					}
				}
#endif

#if 0
				template2Rand(mesh, coords, elem, v_1);
				CheckTemplate(mesh, coords, element_ids_local,false);
				elem = (octant_t*) sc_array_index(&mesh->elements, elements_ids[iel]);

				fprintf(fdbg,"Pad local:    ");
				for(int co = 0; co < element_ids_local.size(); co++){
					octant_t *h = (octant_t*) sc_array_index(&mesh->elements, element_ids_local[co]);
					var_aux_1[co]=h->pad;
					fprintf(fdbg,"%d ",h->pad);
				}
				fprintf(fdbg,"\n");

				if(elem->pad!=4 && elem->pad!=-10){
					flag2=false;
					for(int co = 0; co < element_ids_local.size(); co++){
						if(element_ids_local[co]!=elements_ids[iel]){
							if(var_aux[co]==var_aux_1[co]){
								flag2=false;
							}else{
								if(var_aux_1[co]==4){
									flag2=true;
									break;
								}else if(var_aux_1[co]==-10){
									flag2=true;
									break;
								}
							}
						}
					}
				}

				if(!flag2){
					flag1=false;
					flag2=false;
					printf("To livreeeee rand El:%d, n_iter:%d\n",elements_ids[iel],n_iter);
					break;
				}

				//template4to1(mesh, coords, elements_ids, iel, &point[12], &edge_list[12],1.1);
				// move two nodes "or" one edge
				//if(elem->pad==4){
				//	CuttedEdges (mesh, coords, elements_ids, iel, &point[12], &edge_list[12]);

				GtsSegment * segments[12]={0};
				int Edge2GNode[12][2];
				GtsPoint * point[12] = {NULL};
				bool edge_list[12] = {false};

				for (int edge = 0; edge < 12; ++edge) {
					point[edge] = NULL;

					int node1 = elem->nodes[EdgeVerticesMap[edge][0]].id;
					int node2 = elem->nodes[EdgeVerticesMap[edge][1]].id;
					Edge2GNode[edge][0] = node1 <= node2 ? node1 : node2;
					Edge2GNode[edge][1] = node1 >= node2 ? node1 : node2;

					GtsVertex *v1 = gts_vertex_new(gts_vertex_class(), coords[node1 * 3], coords[node1 * 3 + 1], coords[node1 * 3 + 2]);
					GtsVertex *v2 = gts_vertex_new(gts_vertex_class(), coords[node2 * 3], coords[node2 * 3 + 1], coords[node2 * 3 + 2]);

					segments[edge] = 0;
					segments[edge] = gts_segment_new(gts_segment_class(), v1, v2);

					GtsBBox *sb = gts_bbox_segment(gts_bbox_class(), segments[edge]);
					GSList* list = gts_bb_tree_overlap(mesh->gdata.bbt, sb);

					edge_list[edge] = false;

					if (list == NULL) continue;
					while (list) {
						GtsBBox *b = GTS_BBOX(list->data);
						point[edge] = SegmentTriangleIntersection(segments[edge], GTS_TRIANGLE(b->bounded));
						if (point[edge]) {
							edge_list[edge] = true;
							break;
						}
						list = list->next;
					}
				}

				d_c1=0;
				d_c2=0;

				if(edge_list[0]&&edge_list[1]){
					int node_change = 3*elem->nodes[1].id;

					GtsPoint *p0 = gts_point_new(gts_point_class(),
							coords[node_change],
							coords[node_change + 1],
							coords[node_change + 2]);

					d_c1 = gts_point_distance(point[0], p0);
					d_c2 = gts_point_distance(point[1], p0);

					coords[node_change  ] = coords[node_change  ] - d_c1*tt;
					coords[node_change+1] = coords[node_change+1] + d_c2*tt;

					node_change = 3*elem->nodes[5].id;
					coords[node_change  ] = coords[node_change  ] - d_c1*tt;
					coords[node_change+1] = coords[node_change+1] + d_c2*tt;

					//move 1

				}else if(edge_list[1]&&edge_list[2]){
					int node_change = 3*elem->nodes[2].id;

					GtsPoint *p0 = gts_point_new(gts_point_class(),
							coords[node_change],
							coords[node_change + 1],
							coords[node_change + 2]);

					d_c1 = gts_point_distance(point[2], p0);
					d_c2 = gts_point_distance(point[1], p0);

					coords[node_change  ] = coords[node_change  ] - d_c1*tt;
					coords[node_change+1] = coords[node_change+1] - d_c2*tt;

					node_change = 3*elem->nodes[6].id;
					coords[node_change  ] = coords[node_change  ] - d_c1*tt;
					coords[node_change+1] = coords[node_change+1] - d_c2*tt;
					//move 2

				}else if(edge_list[2]&&edge_list[3]){
					int node_change = 3*elem->nodes[3].id;

					GtsPoint *p0 = gts_point_new(gts_point_class(),
							coords[node_change],
							coords[node_change + 1],
							coords[node_change + 2]);

					d_c1 = gts_point_distance(point[2], p0);
					d_c2 = gts_point_distance(point[3], p0);

					coords[node_change  ] = coords[node_change  ] + d_c1*tt;
					coords[node_change+1] = coords[node_change+1] - d_c2*tt;

					node_change = 3*elem->nodes[7].id;
					coords[node_change  ] = coords[node_change  ] + d_c1*tt;
					coords[node_change+1] = coords[node_change+1] - d_c2*tt;
					//move 3

				}else if(edge_list[3]&&edge_list[0]){
					int node_change = 3*elem->nodes[0].id;

					GtsPoint *p0 = gts_point_new(gts_point_class(),
							coords[node_change],
							coords[node_change + 1],
							coords[node_change + 2]);

					d_c1 = gts_point_distance(point[0], p0);
					d_c2 = gts_point_distance(point[3], p0);

					coords[node_change  ] = coords[node_change  ] + d_c1*tt;
					coords[node_change+1] = coords[node_change+1] + d_c2*tt;

					node_change = 3*elem->nodes[4].id;
					coords[node_change  ] = coords[node_change  ] + d_c1*tt;
					coords[node_change+1] = coords[node_change+1] + d_c2*tt;
					//move 0

				}else if(edge_list[8]&&edge_list[9]){
					int node_change = 3*elem->nodes[5].id;

					GtsPoint *p0 = gts_point_new(gts_point_class(),
							coords[node_change],
							coords[node_change + 1],
							coords[node_change + 2]);

					d_c1 = gts_point_distance(point[8], p0);
					d_c2 = gts_point_distance(point[9], p0);

					coords[node_change  ] = coords[node_change  ] - d_c1*tt;
					coords[node_change+1] = coords[node_change+1] + d_c2*tt;

					node_change = 3*elem->nodes[1].id;
					coords[node_change  ] = coords[node_change  ] - d_c1*tt;
					coords[node_change+1] = coords[node_change+1] + d_c2*tt;
					//move 5

				}else if(edge_list[9]&&edge_list[10]){
					int node_change = 3*elem->nodes[6].id;

					GtsPoint *p0 = gts_point_new(gts_point_class(),
							coords[node_change],
							coords[node_change + 1],
							coords[node_change + 2]);

					d_c1 = gts_point_distance(point[10], p0);
					d_c2 = gts_point_distance(point[9], p0);

					coords[node_change  ] = coords[node_change  ] - d_c1*tt;
					coords[node_change+1] = coords[node_change+1] - d_c2*tt;

					node_change = 3*elem->nodes[2].id;
					coords[node_change  ] = coords[node_change  ] - d_c1*tt;
					coords[node_change+1] = coords[node_change+1] - d_c2*tt;
					//move 6

				}else if(edge_list[10]&&edge_list[11]){
					int node_change = 3*elem->nodes[7].id;

					GtsPoint *p0 = gts_point_new(gts_point_class(),
							coords[node_change],
							coords[node_change + 1],
							coords[node_change + 2]);

					d_c1 = gts_point_distance(point[10], p0);
					d_c2 = gts_point_distance(point[11], p0);

					coords[node_change  ] = coords[node_change  ] + d_c1*tt;
					coords[node_change+1] = coords[node_change+1] - d_c2*tt;

					node_change = 3*elem->nodes[3].id;
					coords[node_change  ] = coords[node_change  ] + d_c1*tt;
					coords[node_change+1] = coords[node_change+1] - d_c2*tt;
					//move 7

				}else if(edge_list[11]&&edge_list[8]){
					int node_change = 3*elem->nodes[4].id;

					GtsPoint *p0 = gts_point_new(gts_point_class(),
							coords[node_change],
							coords[node_change + 1],
							coords[node_change + 2]);

					d_c1 = gts_point_distance(point[8], p0);
					d_c2 = gts_point_distance(point[11], p0);

					coords[node_change  ] = coords[node_change  ] + d_c1*tt;
					coords[node_change+1] = coords[node_change+1] + d_c2*tt;

					node_change = 3*elem->nodes[0].id;
					coords[node_change  ] = coords[node_change  ] + d_c1*tt;
					coords[node_change+1] = coords[node_change+1] + d_c2*tt;
					//move 4

				}else if( edge_list[0]&&edge_list[4] ){
					int node_change = 3*elem->nodes[0].id;

					GtsPoint *p0 = gts_point_new(gts_point_class(),
							coords[node_change],
							coords[node_change + 1],
							coords[node_change + 2]);

					d_c1 = gts_point_distance(point[0], p0);
					d_c2 = gts_point_distance(point[4], p0);

					coords[node_change  ] = coords[node_change  ] + d_c1*tt;
					coords[node_change+2] = coords[node_change+2] - d_c2*tt;

					node_change = 3*elem->nodes[3].id;
					coords[node_change  ] = coords[node_change  ] + d_c1*tt;
					coords[node_change+2] = coords[node_change+2] - d_c2*tt;
					//move 0

				}else if( edge_list[3]&&edge_list[4] ){
					int node_change = 3*elem->nodes[0].id;

					GtsPoint *p0 = gts_point_new(gts_point_class(),
							coords[node_change],
							coords[node_change + 1],
							coords[node_change + 2]);

					d_c1 = gts_point_distance(point[3], p0);
					d_c2 = gts_point_distance(point[4], p0);

					coords[node_change+1] = coords[node_change+1] + d_c1*tt;
					coords[node_change+2] = coords[node_change+2] - d_c2*tt;

					node_change = 3*elem->nodes[1].id;
					coords[node_change+1] = coords[node_change+1] + d_c1*tt;
					coords[node_change+2] = coords[node_change+2] - d_c2*tt;
					//move 0

				}else if(edge_list[0]&&edge_list[5]){
					int node_change = 3*elem->nodes[1].id;

					GtsPoint *p0 = gts_point_new(gts_point_class(),
							coords[node_change],
							coords[node_change + 1],
							coords[node_change + 2]);

					d_c1 = gts_point_distance(point[0], p0);
					d_c2 = gts_point_distance(point[5], p0);

					coords[node_change  ] = coords[node_change  ] - d_c1*tt;
					coords[node_change+2] = coords[node_change+2] - d_c2*tt;

					node_change = 3*elem->nodes[2].id;
					coords[node_change  ] = coords[node_change  ] - d_c1*tt;
					coords[node_change+2] = coords[node_change+2] - d_c2*tt;
					//move 1

				}else if(edge_list[1]&&edge_list[5]){
					int node_change = 3*elem->nodes[1].id;

					GtsPoint *p0 = gts_point_new(gts_point_class(),
							coords[node_change],
							coords[node_change + 1],
							coords[node_change + 2]);

					d_c1 = gts_point_distance(point[1], p0);
					d_c2 = gts_point_distance(point[5], p0);

					coords[node_change+1] = coords[node_change+1] + d_c1*tt;
					coords[node_change+2] = coords[node_change+2] - d_c2*tt;

					node_change = 3*elem->nodes[0].id;
					coords[node_change+1] = coords[node_change+1] + d_c1*tt;
					coords[node_change+2] = coords[node_change+2] - d_c2*tt;
					//move 1

				}else if(edge_list[1]&&edge_list[6]){
					int node_change = 3*elem->nodes[2].id;

					GtsPoint *p0 = gts_point_new(gts_point_class(),
							coords[node_change],
							coords[node_change + 1],
							coords[node_change + 2]);

					d_c1 = gts_point_distance(point[1], p0);
					d_c2 = gts_point_distance(point[6], p0);

					coords[node_change+1] = coords[node_change+1] - d_c1*tt;
					coords[node_change+2] = coords[node_change+2] - d_c2*tt;

					node_change = 3*elem->nodes[3].id;
					coords[node_change+1] = coords[node_change+1] - d_c1*tt;
					coords[node_change+2] = coords[node_change+2] - d_c2*tt;
					//move 2

				}else if(edge_list[2]&&edge_list[6]){
					int node_change = 3*elem->nodes[2].id;

					GtsPoint *p0 = gts_point_new(gts_point_class(),
							coords[node_change],
							coords[node_change + 1],
							coords[node_change + 2]);

					d_c1 = gts_point_distance(point[2], p0);
					d_c2 = gts_point_distance(point[6], p0);

					coords[node_change  ] = coords[node_change  ] - d_c1*tt;
					coords[node_change+2] = coords[node_change+2] - d_c2*tt;

					node_change = 3*elem->nodes[1].id;
					coords[node_change  ] = coords[node_change  ] - d_c1*tt;
					coords[node_change+2] = coords[node_change+2] - d_c2*tt;
					//move 2

				}else if(edge_list[2]&&edge_list[7]){
					int node_change = 3*elem->nodes[3].id;

					GtsPoint *p0 = gts_point_new(gts_point_class(),
							coords[node_change],
							coords[node_change + 1],
							coords[node_change + 2]);

					d_c1 = gts_point_distance(point[2], p0);
					d_c2 = gts_point_distance(point[7], p0);

					coords[node_change  ] = coords[node_change  ] + d_c1*tt;
					coords[node_change+2] = coords[node_change+2] - d_c2*tt;

					node_change = 3*elem->nodes[0].id;
					coords[node_change  ] = coords[node_change  ] + d_c1*tt;
					coords[node_change+2] = coords[node_change+2] - d_c2*tt;
					//move 3

				}else if(edge_list[3]&&edge_list[7]){
					int node_change = 3*elem->nodes[3].id;

					GtsPoint *p0 = gts_point_new(gts_point_class(),
							coords[node_change],
							coords[node_change + 1],
							coords[node_change + 2]);

					d_c1 = gts_point_distance(point[3], p0);
					d_c2 = gts_point_distance(point[7], p0);

					coords[node_change+1] = coords[node_change+1] - d_c1*tt;
					coords[node_change+2] = coords[node_change+2] - d_c2*tt;

					node_change = 3*elem->nodes[2].id;
					coords[node_change+1] = coords[node_change+1] - d_c1*tt;
					coords[node_change+2] = coords[node_change+2] - d_c2*tt;
					//move 3

				}else if(edge_list[8]&&edge_list[4]){
					int node_change = 3*elem->nodes[4].id;

					GtsPoint *p0 = gts_point_new(gts_point_class(),
							coords[node_change],
							coords[node_change + 1],
							coords[node_change + 2]);

					d_c1 = gts_point_distance(point[8], p0);
					d_c2 = gts_point_distance(point[4], p0);

					coords[node_change  ] = coords[node_change  ] + d_c1*tt;
					coords[node_change+2] = coords[node_change+2] + d_c2*tt;

					node_change = 3*elem->nodes[7].id;
					coords[node_change  ] = coords[node_change  ] + d_c1*tt;
					coords[node_change+2] = coords[node_change+2] + d_c2*tt;
					//move 4

				}else if(edge_list[11]&&edge_list[4]){
					int node_change = 3*elem->nodes[4].id;

					GtsPoint *p0 = gts_point_new(gts_point_class(),
							coords[node_change],
							coords[node_change + 1],
							coords[node_change + 2]);

					d_c1 = gts_point_distance(point[11], p0);
					d_c2 = gts_point_distance(point[4], p0);

					coords[node_change+1] = coords[node_change+1] + d_c1*tt;
					coords[node_change+2] = coords[node_change+2] + d_c2*tt;

					node_change = 3*elem->nodes[5].id;
					coords[node_change+1] = coords[node_change+1] + d_c1*tt;
					coords[node_change+2] = coords[node_change+2] + d_c2*tt;
					//move 4

				}else if(edge_list[8]&&edge_list[5]){
					int node_change = 3*elem->nodes[5].id;

					GtsPoint *p0 = gts_point_new(gts_point_class(),
							coords[node_change],
							coords[node_change + 1],
							coords[node_change + 2]);

					d_c1 = gts_point_distance(point[8], p0);
					d_c2 = gts_point_distance(point[5], p0);

					coords[node_change  ] = coords[node_change  ] - d_c1*tt;
					coords[node_change+2] = coords[node_change+2] + d_c2*tt;

					node_change = 3*elem->nodes[6].id;
					coords[node_change  ] = coords[node_change  ] - d_c1*tt;
					coords[node_change+2] = coords[node_change+2] + d_c2*tt;
					//move 5

				}else if(edge_list[9]&&edge_list[5]){
					int node_change = 3*elem->nodes[5].id;

					GtsPoint *p0 = gts_point_new(gts_point_class(),
							coords[node_change],
							coords[node_change + 1],
							coords[node_change + 2]);

					d_c1 = gts_point_distance(point[9], p0);
					d_c2 = gts_point_distance(point[5], p0);

					coords[node_change+1] = coords[node_change+1] + d_c1*tt;
					coords[node_change+2] = coords[node_change+2] + d_c2*tt;

					node_change = 3*elem->nodes[4].id;
					coords[node_change+1] = coords[node_change+1] + d_c1*tt;
					coords[node_change+2] = coords[node_change+2] + d_c2*tt;
					//move 5

				}else if(edge_list[9]&&edge_list[6]){
					int node_change = 3*elem->nodes[6].id;

					GtsPoint *p0 = gts_point_new(gts_point_class(),
							coords[node_change],
							coords[node_change + 1],
							coords[node_change + 2]);

					d_c1 = gts_point_distance(point[9], p0);
					d_c2 = gts_point_distance(point[6], p0);

					coords[node_change+1] = coords[node_change+1] - d_c1*tt;
					coords[node_change+2] = coords[node_change+2] + d_c2*tt;

					node_change = 3*elem->nodes[7].id;
					coords[node_change+1] = coords[node_change+1] - d_c1*tt;
					coords[node_change+2] = coords[node_change+2] + d_c2*tt;
					//move 6

				}else if(edge_list[10]&&edge_list[6]){
					int node_change = 3*elem->nodes[6].id;

					GtsPoint *p0 = gts_point_new(gts_point_class(),
							coords[node_change],
							coords[node_change + 1],
							coords[node_change + 2]);

					d_c1 = gts_point_distance(point[10], p0);
					d_c2 = gts_point_distance(point[6], p0);

					coords[node_change  ] = coords[node_change  ] - d_c1*tt;
					coords[node_change+2] = coords[node_change+2] + d_c2*tt;

					node_change = 3*elem->nodes[5].id;
					coords[node_change  ] = coords[node_change  ] - d_c1*tt;
					coords[node_change+2] = coords[node_change+2] + d_c2*tt;
					//move 6

				}else if(edge_list[10]&&edge_list[7]){
					int node_change = 3*elem->nodes[7].id;

					GtsPoint *p0 = gts_point_new(gts_point_class(),
							coords[node_change],
							coords[node_change + 1],
							coords[node_change + 2]);

					d_c1 = gts_point_distance(point[10], p0);
					d_c2 = gts_point_distance(point[7], p0);

					coords[node_change  ] = coords[node_change  ] + d_c1*tt;
					coords[node_change+2] = coords[node_change+2] + d_c2*tt;

					node_change = 3*elem->nodes[4].id;
					coords[node_change  ] = coords[node_change  ] + d_c1*tt;
					coords[node_change+2] = coords[node_change+2] + d_c2*tt;
					//move 7

				}else if(edge_list[11]&&edge_list[7]){
					int node_change = 3*elem->nodes[7].id;

					GtsPoint *p0 = gts_point_new(gts_point_class(),
							coords[node_change],
							coords[node_change + 1],
							coords[node_change + 2]);

					d_c1 = gts_point_distance(point[11], p0);
					d_c2 = gts_point_distance(point[7], p0);

					coords[node_change+1] = coords[node_change+1] - d_c1*tt;
					coords[node_change+2] = coords[node_change+2] + d_c2*tt;

					node_change = 3*elem->nodes[6].id;
					coords[node_change+1] = coords[node_change+1] - d_c1*tt;
					coords[node_change+2] = coords[node_change+2] + d_c2*tt;
					//move 7
				}

				CheckTemplate(mesh, coords, element_ids_local,false);
				elem = (octant_t*) sc_array_index(&mesh->elements, elements_ids[iel]);

				fprintf(fdbg,"Pad local:    ");
				for(int co = 0; co < element_ids_local.size(); co++){
					octant_t *h = (octant_t*) sc_array_index(&mesh->elements, element_ids_local[co]);
					var_aux_1[co]=h->pad;
					fprintf(fdbg,"%d ",h->pad);
				}
				fprintf(fdbg,"\n");


				if(elem->pad!=4 && elem->pad!=-10){
					flag2=false;
					for(int co = 0; co < element_ids_local.size(); co++){
						if(element_ids_local[co]!=elements_ids[iel]){
							if(var_aux[co]==var_aux_1[co]){
								flag2=false;
							}else{
								if(var_aux_1[co]==4){
									flag2=true;
									break;
								}else if(var_aux_1[co]==-10){
									flag2=true;
									break;
								}
							}
						}
					}
				}
				//}

				if(!flag2){
					flag1=false;
					flag2=false;
					printf("To livreeeee mov 2 El:%d, n_iter:%d\n",elements_ids[iel],n_iter);
				}

				if( flag2 ){
					for(int n_nodes = 0; n_nodes<8;n_nodes++){
						int node = 3*elem->nodes[n_nodes].id;
						coords[node]   = coords_orgi[3*n_nodes];
						coords[node+1] = coords_orgi[3*n_nodes+1];
						coords[node+2] = coords_orgi[3*n_nodes+2];
					}
				}
#endif

#if 0
				template2Rand(mesh, coords, elem, v_1);
				CheckTemplate(mesh, coords, element_ids_local,false);
				elem = (octant_t*) sc_array_index(&mesh->elements, elements_ids[iel]);

				fprintf(fdbg,"Pad local:    ");
				for(int co = 0; co < element_ids_local.size(); co++){
					octant_t *h = (octant_t*) sc_array_index(&mesh->elements, element_ids_local[co]);
					var_aux_1[co]=h->pad;
					fprintf(fdbg,"%d ",h->pad);
				}
				fprintf(fdbg,"\n");

				if(elem->pad!=4 && elem->pad!=-10){
					flag2=false;
					for(int co = 0; co < element_ids_local.size(); co++){
						if(element_ids_local[co]!=elements_ids[iel]){
							if(var_aux[co]==var_aux_1[co]){
								flag2=false;
							}else{
								if(var_aux_1[co]==4){
									flag2=true;
									break;
								}else if(var_aux_1[co]==-10){
									flag2=true;
									break;
								}
							}
						}
					}
				}

				if(!flag2){
					flag1=false;
					flag2=false;
					printf("To livreeeee rand El:%d, n_iter:%d\n",elements_ids[iel],n_iter);
					break;
				}

				//template4to2(mesh, coords, elements_ids, iel, &point[12], &edge_list[12],1.1);
				if(elem->pad==4){

					GtsSegment * segments[12]={0};
					int Edge2GNode[12][2];
					GtsPoint * point[12] = {NULL};
					bool edge_list[12] = {false};

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

						edge_list[edge] = false;

						if (list == NULL) continue;
						while (list) {
							GtsBBox *b = GTS_BBOX(list->data);
							point[edge] = SegmentTriangleIntersection(segments[edge], GTS_TRIANGLE(b->bounded));
							if (point[edge]) {
								edge_list[edge] = true;
								break;
							}
							list = list->next;
						}
					}

					double d_c1;
					double d_c2;

					if(edge_list[0]&&edge_list[1]){
						int node_change = 3*elem->nodes[1].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[0], p0);
						d_c2 = gts_point_distance(point[1], p0);

						if(d_c1<=d_c2){
							coords[node_change  ] = coords[node_change  ] - d_c1*tt;
						}else{
							coords[node_change+1] = coords[node_change+1] + d_c2*tt;
						}
						//move 1

					}else if(edge_list[1]&&edge_list[2]){
						int node_change = 3*elem->nodes[2].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[2], p0);
						d_c2 = gts_point_distance(point[1], p0);

						if(d_c1<=d_c2){
							coords[node_change  ] = coords[node_change  ] - d_c1*tt;
						}else{
							coords[node_change+1] = coords[node_change+1] - d_c2*tt;
						}
						//move 2

					}else if(edge_list[2]&&edge_list[3]){
						int node_change = 3*elem->nodes[3].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[2], p0);
						d_c2 = gts_point_distance(point[3], p0);

						if(d_c1<=d_c2){
							coords[node_change  ] = coords[node_change  ] + d_c1*tt;
						}else{
							coords[node_change+1] = coords[node_change+1] - d_c2*tt;
						}
						//move 3

					}else if(edge_list[3]&&edge_list[0]){
						int node_change = 3*elem->nodes[0].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[0], p0);
						d_c2 = gts_point_distance(point[3], p0);

						if(d_c1<=d_c2){
							coords[node_change  ] = coords[node_change  ] + d_c1*tt;
						}else{
							coords[node_change+1] = coords[node_change+1] + d_c2*tt;
						}
						//move 0

					}else if(edge_list[8]&&edge_list[9]){
						int node_change = 3*elem->nodes[5].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[8], p0);
						d_c2 = gts_point_distance(point[9], p0);

						if(d_c1<=d_c2){
							coords[node_change  ] = coords[node_change  ] - d_c1*tt;
						}else{
							coords[node_change+1] = coords[node_change+1] + d_c2*tt;
						}
						//move 5

					}else if(edge_list[9]&&edge_list[10]){
						int node_change = 3*elem->nodes[6].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[10], p0);
						d_c2 = gts_point_distance(point[9], p0);

						if(d_c1<=d_c2){
							coords[node_change  ] = coords[node_change  ] - d_c1*tt;
						}else{
							coords[node_change+1] = coords[node_change+1] - d_c2*tt;
						}
						//move 6

					}else if(edge_list[10]&&edge_list[11]){
						int node_change = 3*elem->nodes[7].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[10], p0);
						d_c2 = gts_point_distance(point[11], p0);

						if(d_c1<=d_c2){
							coords[node_change  ] = coords[node_change  ] + d_c1*tt;
						}else{
							coords[node_change+1] = coords[node_change+1] - d_c2*tt;
						}
						//move 7

					}else if(edge_list[11]&&edge_list[8]){
						int node_change = 3*elem->nodes[4].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[8], p0);
						d_c2 = gts_point_distance(point[11], p0);

						if(d_c1<=d_c2){
							coords[node_change  ] = coords[node_change  ] + d_c1*tt;
						}else{
							coords[node_change+1] = coords[node_change+1] + d_c2*tt;
						}
						//move 4

					}else if( edge_list[0]&&edge_list[4] ){
						int node_change = 3*elem->nodes[0].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[0], p0);
						d_c2 = gts_point_distance(point[4], p0);

						if(d_c1<=d_c2){
							coords[node_change  ] = coords[node_change  ] + d_c1*tt;
						}else{
							coords[node_change+2] = coords[node_change+2] - d_c2*tt;
						}
						//move 0

					}else if( edge_list[3]&&edge_list[4] ){
						int node_change = 3*elem->nodes[0].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[3], p0);
						d_c2 = gts_point_distance(point[4], p0);

						if(d_c1<=d_c2){
							coords[node_change+1] = coords[node_change+1] + d_c1*tt;
						}else{
							coords[node_change+2] = coords[node_change+2] - d_c2*tt;
						}
						//move 0

					}else if(edge_list[0]&&edge_list[5]){
						int node_change = 3*elem->nodes[1].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[0], p0);
						d_c2 = gts_point_distance(point[5], p0);

						if(d_c1<=d_c2){
							coords[node_change  ] = coords[node_change  ] - d_c1*tt;
						}else{
							coords[node_change+2] = coords[node_change+2] - d_c2*tt;
						}
						//move 1

					}else if(edge_list[1]&&edge_list[5]){
						int node_change = 3*elem->nodes[1].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[1], p0);
						d_c2 = gts_point_distance(point[5], p0);

						if(d_c1<=d_c2){
							coords[node_change+1] = coords[node_change+1] + d_c1*tt;
						}else{
							coords[node_change+2] = coords[node_change+2] - d_c2*tt;
						}
						//move 1

					}else if(edge_list[1]&&edge_list[6]){
						int node_change = 3*elem->nodes[2].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[1], p0);
						d_c2 = gts_point_distance(point[6], p0);

						if(d_c1<=d_c2){
							coords[node_change+1] = coords[node_change+1] - d_c1*tt;
						}else{
							coords[node_change+2] = coords[node_change+2] - d_c2*tt;
						}
						//move 2

					}else if(edge_list[2]&&edge_list[6]){
						int node_change = 3*elem->nodes[2].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[2], p0);
						d_c2 = gts_point_distance(point[6], p0);

						if(d_c1<=d_c2){
							coords[node_change  ] = coords[node_change  ] - d_c1*tt;
						}else{
							coords[node_change+2] = coords[node_change+2] - d_c2*tt;
						}
						//move 2

					}else if(edge_list[2]&&edge_list[7]){
						int node_change = 3*elem->nodes[3].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[2], p0);
						d_c2 = gts_point_distance(point[7], p0);

						if(d_c1<=d_c2){
							coords[node_change  ] = coords[node_change  ] + d_c1*tt;
						}else{
							coords[node_change+2] = coords[node_change+2] - d_c2*tt;
						}
						//move 3

					}else if(edge_list[3]&&edge_list[7]){
						int node_change = 3*elem->nodes[3].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[3], p0);
						d_c2 = gts_point_distance(point[7], p0);

						if(d_c1<=d_c2){
							coords[node_change+1] = coords[node_change+1] - d_c1*tt;
						}else{
							coords[node_change+2] = coords[node_change+2] - d_c2*tt;
						}
						//move 3

					}else if(edge_list[8]&&edge_list[4]){
						int node_change = 3*elem->nodes[4].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[8], p0);
						d_c2 = gts_point_distance(point[4], p0);

						if(d_c1<=d_c2){
							coords[node_change  ] = coords[node_change  ] + d_c1*tt;
						}else{
							coords[node_change+2] = coords[node_change+2] + d_c2*tt;
						}
						//move 4

					}else if(edge_list[11]&&edge_list[4]){
						int node_change = 3*elem->nodes[4].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[11], p0);
						d_c2 = gts_point_distance(point[4], p0);

						if(d_c1<=d_c2){
							coords[node_change+1] = coords[node_change+1] + d_c1*tt;
						}else{
							coords[node_change+2] = coords[node_change+2] + d_c2*tt;
						}
						//move 4

					}else if(edge_list[8]&&edge_list[5]){
						int node_change = 3*elem->nodes[5].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[8], p0);
						d_c2 = gts_point_distance(point[5], p0);

						if(d_c1<=d_c2){
							coords[node_change  ] = coords[node_change  ] - d_c1*tt;
						}else{
							coords[node_change+2] = coords[node_change+2] + d_c2*tt;
						}
						//move 5

					}else if(edge_list[9]&&edge_list[5]){
						int node_change = 3*elem->nodes[5].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[9], p0);
						d_c2 = gts_point_distance(point[5], p0);

						if(d_c1<=d_c2){
							coords[node_change+1] = coords[node_change+1] + d_c1*tt;
						}else{
							coords[node_change+2] = coords[node_change+2] + d_c2*tt;
						}
						//move 5

					}else if(edge_list[9]&&edge_list[6]){
						int node_change = 3*elem->nodes[6].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[9], p0);
						d_c2 = gts_point_distance(point[6], p0);

						if(d_c1<=d_c2){
							coords[node_change+1] = coords[node_change+1] - d_c1*tt;
						}else{
							coords[node_change+2] = coords[node_change+2] + d_c2*tt;
						}
						//move 6

					}else if(edge_list[10]&&edge_list[6]){
						int node_change = 3*elem->nodes[6].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[10], p0);
						d_c2 = gts_point_distance(point[6], p0);

						if(d_c1<=d_c2){
							coords[node_change  ] = coords[node_change  ] - d_c1*tt;
						}else{
							coords[node_change+2] = coords[node_change+2] + d_c2*tt;
						}
						//move 6

					}else if(edge_list[10]&&edge_list[7]){
						int node_change = 3*elem->nodes[7].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[10], p0);
						d_c2 = gts_point_distance(point[7], p0);

						if(d_c1<=d_c2){
							coords[node_change  ] = coords[node_change  ] + d_c1*tt;
						}else{
							coords[node_change+2] = coords[node_change+2] + d_c2*tt;
						}
						//move 7

					}else if(edge_list[11]&&edge_list[7]){
						int node_change = 3*elem->nodes[7].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[11], p0);
						d_c2 = gts_point_distance(point[7], p0);

						if(d_c1<=d_c2){
							coords[node_change+1] = coords[node_change+1] - d_c1*tt;
						}else{
							coords[node_change+2] = coords[node_change+2] + d_c2*tt;
						}
						//move 7
					}

					CheckTemplate(mesh, coords, element_ids_local,false);
					elem = (octant_t*) sc_array_index(&mesh->elements, elements_ids[iel]);

					fprintf(fdbg,"Pad local:    ");
					for(int co = 0; co < element_ids_local.size(); co++){
						octant_t *h = (octant_t*) sc_array_index(&mesh->elements, element_ids_local[co]);
						var_aux_1[co]=h->pad;
						fprintf(fdbg,"%d ",h->pad);
					}
					fprintf(fdbg,"\n");

					if(elem->pad!=4 && elem->pad!=-10){
						flag2=false;
						for(int co = 0; co < element_ids_local.size(); co++){
							if(element_ids_local[co]!=elements_ids[iel]){
								if(var_aux[co]==var_aux_1[co]){
									flag2=false;
								}else{
									if(var_aux_1[co]==4){
										flag2=true;
										break;
									}else if(var_aux_1[co]==-10){
										flag2=true;
										break;
									}
								}
							}
						}
					}
				}

				if(!flag2){
					flag1=false;
					flag2=false;
					printf("To livreeeee mov 3 El:%d, n_iter:%d\n",elements_ids[iel],n_iter);
				}

				if( flag2 ){
					for(int n_nodes = 0; n_nodes<8;n_nodes++){
						int node = 3*elem->nodes[n_nodes].id;
						coords[node]   = coords_orgi[3*n_nodes];
						coords[node+1] = coords_orgi[3*n_nodes+1];
						coords[node+2] = coords_orgi[3*n_nodes+2];
					}
				}
#endif

#if 0 //isso aqui não funciona...

				if (elem->z!=0){

					if(n_iter==0){
						printf("entrou aqui EL: %d\n",elements_ids[iel]);
					}

					elem = (octant_t*) sc_array_index(&mesh->elements, elements_ids[iel]);

					double dz = (coords[3*elem->nodes[4].id+ 3]-coords[3*elem->nodes[0].id+ 3]);
					dz = abs(dz);

					double nn = (50 - rand() % 100)*0.1;

					for(int co = 0; co<8; co++){
						int node_change = 3*elem->nodes[co].id;
						coords[node_change + 3] = coords[node_change + 3] + dz*nn;
					}

					CheckTemplate(mesh, coords, element_ids_local,false);
					elem = (octant_t*) sc_array_index(&mesh->elements, elements_ids[iel]);

					fprintf(fdbg,"Pad local:    ");
					for(int co = 0; co < element_ids_local.size(); co++){
						octant_t *h = (octant_t*) sc_array_index(&mesh->elements, element_ids_local[co]);
						var_aux_1[co]=h->pad;
						fprintf(fdbg,"%d ",h->pad);
					}
					fprintf(fdbg,"\n");

					if(elem->pad!=4 && elem->pad!=-10){
						flag2=false;
						for(int co = 0; co < element_ids_local.size(); co++){
							if(element_ids_local[co]!=elements_ids[iel]){
								if(var_aux[co]==var_aux_1[co]){
									flag2=false;
								}else{
									if(var_aux_1[co]==4){
										flag2=true;
										break;
									}else if(var_aux_1[co]==-10){
										flag2=true;
										break;
									}
								}
							}
						}
					}

					if(!flag2){
						flag1=false;
						flag2=false;
						printf("To livreeeee todos juntinhos El:%d, n_iter:%d\n",elements_ids[iel],n_iter);
						break;
					}
				}

				if( flag2 ){
					for(int n_nodes = 0; n_nodes<8;n_nodes++){
						int node = 3*elem->nodes[n_nodes].id;
						coords[node]   = coords_orgi[3*n_nodes];
						coords[node+1] = coords_orgi[3*n_nodes+1];
						coords[node+2] = coords_orgi[3*n_nodes+2];
					}
				}
#endif

#if 0

				template2Rand(mesh, coords, elem, v_1);
				CheckTemplate(mesh, coords, element_ids_local,false);
				elem = (octant_t*) sc_array_index(&mesh->elements, elements_ids[iel]);

				fprintf(fdbg,"Pad local:    ");
				for(int co = 0; co < element_ids_local.size(); co++){
					octant_t *h = (octant_t*) sc_array_index(&mesh->elements, element_ids_local[co]);
					var_aux_1[co]=h->pad;
					fprintf(fdbg,"%d ",h->pad);
				}
				fprintf(fdbg,"\n");

				if(elem->pad!=4 && elem->pad!=-10){
					flag2=false;
					for(int co = 0; co < element_ids_local.size(); co++){
						if(element_ids_local[co]!=elements_ids[iel]){
							if(var_aux[co]==var_aux_1[co]){
								flag2=false;
							}else{
								if(var_aux_1[co]==4){
									flag2=true;
									break;
								}else if(var_aux_1[co]==-10){
									flag2=true;
									break;
								}
							}
						}
					}
				}

				if(!flag2){
					flag1=false;
					flag2=false;
					printf("To livreeeee rand El:%d, n_iter:%d\n",elements_ids[iel],n_iter);
					break;
				}

				if(elem->pad==4){
					//	CuttedEdges (mesh, coords, elements_ids, iel, &point[12], &edge_list[12]);

					GtsSegment * segments[12]={0};
					int Edge2GNode[12][2];
					GtsPoint * point[12] = {NULL};
					bool edge_list[12] = {false};

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

						edge_list[edge] = false;

						if (list == NULL) continue;
						while (list) {
							GtsBBox *b = GTS_BBOX(list->data);
							point[edge] = SegmentTriangleIntersection(segments[edge], GTS_TRIANGLE(b->bounded));
							if (point[edge]) {
								edge_list[edge] = true;
								break;
							}
							list = list->next;
						}
					}

					double d_c1;
					double d_c2;

					double ref1 = 	coords[3*elem->nodes[0].id ] ;
					double ref2 = 	coords[3*elem->nodes[0].id + 1 ] ;
					double ref3 = coords[3*elem->nodes[0].id + 3 ] ;
					double dx = coords[3*elem->nodes[1].id ] - coords[3*elem->nodes[0].id ];
					double dy = coords[3*elem->nodes[2].id +1] - coords[3*elem->nodes[1].id+1 ];
					double dz = coords[3*elem->nodes[0].id +2] - coords[3*elem->nodes[4].id+2 ];

					if(edge_list[0]&&edge_list[1]){
						int node_change = 3*elem->nodes[1].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[0], p0);
						d_c2 = gts_point_distance(point[1], p0);

						double alpha = atan(d_c1/d_c2);

						for(int co = 0; co<8;co++){
							node_change = 3*elem->nodes[co].id;
							coords[node_change  ] = coords[node_change   ] - ref1 - dx*0.5;
							coords[node_change +1  ] =coords[node_change +1  ] -ref2-dy/2;
							coords[node_change  ] = coords[node_change  ]*cos(alpha)  - coords[node_change+1]*sin(alpha);
							coords[node_change +1  ] = coords[node_change]*sin(alpha) +coords[node_change +1  ]*cos(alpha);
							coords[node_change  ] = coords[node_change   ] + ref1 + dx*0.5;
							coords[node_change +1  ] =coords[node_change +1  ] + ref2 + dy/2;
						}
						//move 1

					}else if(edge_list[1]&&edge_list[2]){
						int node_change = 3*elem->nodes[2].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[2], p0);
						d_c2 = gts_point_distance(point[1], p0);

						double alpha = atan(d_c1/d_c2);

						for(int co = 0; co<8;co++){
							node_change = 3*elem->nodes[co].id;
							coords[node_change  ] = coords[node_change   ] - ref1 - dx*0.5;
							coords[node_change +1  ] =coords[node_change +1  ] -ref2-dy/2;
							coords[node_change  ] = coords[node_change  ]*cos(alpha)  - coords[node_change+1]*sin(alpha);
							coords[node_change +1  ] = coords[node_change  ]*sin(alpha) +coords[node_change +1  ]*cos(alpha);
							coords[node_change  ] = coords[node_change   ] + ref1 + dx*0.5;
							coords[node_change +1  ] =coords[node_change +1  ] + ref2 + dy/2;
						}
						//move 2

					}else if(edge_list[2]&&edge_list[3]){
						int node_change = 3*elem->nodes[3].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[2], p0);
						d_c2 = gts_point_distance(point[3], p0);

						double alpha = atan(d_c1/d_c2);

						for(int co = 0; co<8;co++){
							node_change = 3*elem->nodes[co].id;
							coords[node_change  ] = coords[node_change   ] - ref1 - dx*0.5;
							coords[node_change +1  ] =coords[node_change +1  ] -ref2-dy/2;
							coords[node_change  ] = coords[node_change  ]*cos(alpha)  - coords[node_change+1]*sin(alpha);
							coords[node_change +1  ] = coords[node_change  ]*sin(alpha) +coords[node_change +1  ]*cos(alpha);
							coords[node_change  ] = coords[node_change   ] + ref1 + dx*0.5;
							coords[node_change +1  ] =coords[node_change +1  ] + ref2 + dy/2;
						}
						//move 3

					}else if(edge_list[3]&&edge_list[0]){
						int node_change = 3*elem->nodes[0].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[0], p0);
						d_c2 = gts_point_distance(point[3], p0);

						double alpha = atan(d_c1/d_c2);

						for(int co = 0; co<8;co++){
							node_change = 3*elem->nodes[co].id;
							coords[node_change  ] = coords[node_change   ] - ref1 - dx*0.5;
							coords[node_change +1  ] =coords[node_change +1  ] -ref2-dy/2;
							coords[node_change  ] = coords[node_change  ]*cos(alpha)  - coords[node_change+1]*sin(alpha);
							coords[node_change +1  ] = coords[node_change  ]*sin(alpha) +coords[node_change +1  ]*cos(alpha);
							coords[node_change  ] = coords[node_change   ] + ref1 + dx*0.5;
							coords[node_change +1  ] =coords[node_change +1  ] + ref2 + dy/2;
						}
						//move 0

					}else if(edge_list[8]&&edge_list[9]){
						int node_change = 3*elem->nodes[5].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[8], p0);
						d_c2 = gts_point_distance(point[9], p0);

						double alpha = atan(d_c1/d_c2);

						for(int co = 0; co<8;co++){
							node_change = 3*elem->nodes[co].id;
							coords[node_change  ] = coords[node_change   ] - ref1 - dx*0.5;
							coords[node_change +1  ] =coords[node_change +1  ] -ref2-dy/2;
							coords[node_change  ] = coords[node_change  ]*cos(alpha)  - coords[node_change+1]*sin(alpha);
							coords[node_change +1  ] = coords[node_change]*sin(alpha) +coords[node_change +1  ]*cos(alpha);
							coords[node_change  ] = coords[node_change   ] + ref1 + dx*0.5;
							coords[node_change +1  ] =coords[node_change +1  ] + ref2 + dy/2;
						}
						//move 5

					}else if(edge_list[9]&&edge_list[10]){
						int node_change = 3*elem->nodes[6].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[10], p0);
						d_c2 = gts_point_distance(point[9], p0);

						double alpha = atan(d_c1/d_c2);

						for(int co = 0; co<8;co++){
							node_change = 3*elem->nodes[co].id;
							coords[node_change  ] = coords[node_change   ] - ref1 - dx*0.5;
							coords[node_change +1  ] =coords[node_change +1  ] -ref2-dy/2;
							coords[node_change  ] = coords[node_change  ]*cos(alpha)  - coords[node_change+1]*sin(alpha);
							coords[node_change +1  ] = coords[node_change ]*sin(alpha) +coords[node_change +1  ]*cos(alpha);
							coords[node_change  ] = coords[node_change   ] + ref1 + dx*0.5;
							coords[node_change +1  ] =coords[node_change +1  ] + ref2 + dy/2;
						}
						//move 6

					}else if(edge_list[10]&&edge_list[11]){
						int node_change = 3*elem->nodes[7].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[10], p0);
						d_c2 = gts_point_distance(point[11], p0);

						double alpha = atan(d_c1/d_c2);

						for(int co = 0; co<8;co++){
							node_change = 3*elem->nodes[co].id;
							coords[node_change  ] = coords[node_change   ] - ref1 - dx*0.5;
							coords[node_change +1  ] =coords[node_change +1  ] -ref2-dy/2;
							coords[node_change  ] = coords[node_change  ]*cos(alpha)  - coords[node_change+1]*sin(alpha);
							coords[node_change +1  ] = coords[node_change]*sin(alpha) +coords[node_change +1  ]*cos(alpha);
							coords[node_change  ] = coords[node_change   ] + ref1 + dx*0.5;
							coords[node_change +1  ] =coords[node_change +1  ] + ref2 + dy/2;
						}
						//move 7

					}else if(edge_list[11]&&edge_list[8]){
						int node_change = 3*elem->nodes[4].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[8], p0);
						d_c2 = gts_point_distance(point[11], p0);

						double alpha = atan(d_c1/d_c2);

						for(int co = 0; co<8;co++){
							node_change = 3*elem->nodes[co].id;
							coords[node_change  ] = coords[node_change   ] - ref1 - dx*0.5;
							coords[node_change +1  ] =coords[node_change +1  ] -ref2-dy/2;
							coords[node_change  ] = coords[node_change  ]*cos(alpha)  - coords[node_change+1]*sin(alpha);
							coords[node_change +1  ] = coords[node_change ]*sin(alpha) +coords[node_change +1  ]*cos(alpha);
							coords[node_change  ] = coords[node_change   ] + ref1 + dx*0.5;
							coords[node_change +1  ] =coords[node_change +1  ] + ref2 + dy/2;
						}
						//move 4

					}else if( edge_list[0]&&edge_list[4] ){
						int node_change = 3*elem->nodes[0].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[0], p0);
						d_c2 = gts_point_distance(point[4], p0);

						double alpha = atan(d_c1/d_c2);

						for(int co = 0; co<8;co++){
							node_change = 3*elem->nodes[co].id;
							coords[node_change  ] = coords[node_change   ] - ref1 - dx*0.5;
							coords[node_change +2  ] =coords[node_change +2  ] -ref3-dz/2;
							coords[node_change  ] = coords[node_change  ]*cos(alpha)  - coords[node_change+2]*sin(alpha);
							coords[node_change +2  ] = coords[node_change]*sin(alpha) +coords[node_change +2  ]*cos(alpha);
							coords[node_change  ] = coords[node_change   ] + ref1 + dx*0.5;
							coords[node_change +2  ] =coords[node_change +1  ] + ref3 + dz/2;
						}
						//move 0

					}else if( edge_list[3]&&edge_list[4] ){
						int node_change = 3*elem->nodes[0].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[3], p0);
						d_c2 = gts_point_distance(point[4], p0);

						double alpha = atan(d_c1/d_c2);

						for(int co = 0; co<8;co++){
							node_change = 3*elem->nodes[co].id;
							coords[node_change + 1] = coords[node_change + 1] - ref2 - dy*0.5;
							coords[node_change + 2] = coords[node_change + 2] - ref3 - dz*0.5;
							coords[node_change + 1] = coords[node_change + 1]*cos(alpha) - coords[node_change + 2]*sin(alpha);
							coords[node_change + 2] = coords[node_change + 1]*sin(alpha) + coords[node_change + 2]*cos(alpha);
							coords[node_change + 1] = coords[node_change + 1] + ref2 + dy*0.5;
							coords[node_change + 2] = coords[node_change + 2] + ref3 + dz*0.5;
						}
						//move 0

					}else if(edge_list[0]&&edge_list[5]){
						int node_change = 3*elem->nodes[1].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[0], p0);
						d_c2 = gts_point_distance(point[5], p0);

						double alpha = atan(d_c1/d_c2);

						for(int co = 0; co<8;co++){
							node_change = 3*elem->nodes[co].id;
							coords[node_change    ] = coords[node_change    ] - ref1 - dx*0.5;
							coords[node_change + 2] = coords[node_change + 2] - ref3 - dz*0.5;
							coords[node_change    ] = coords[node_change    ]*cos(alpha) - coords[node_change + 2]*sin(alpha);
							coords[node_change + 2] = coords[node_change    ]*sin(alpha) + coords[node_change + 2]*cos(alpha);
							coords[node_change    ] = coords[node_change    ] + ref1 + dx*0.5;
							coords[node_change + 2] = coords[node_change + 2] + ref3 + dz*0.5;
						}
						//move 1

					}else if(edge_list[1]&&edge_list[5]){
						int node_change = 3*elem->nodes[1].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[1], p0);
						d_c2 = gts_point_distance(point[5], p0);

						double alpha = atan(d_c1/d_c2);

						for(int co = 0; co<8;co++){
							node_change = 3*elem->nodes[co].id;
							coords[node_change + 1] = coords[node_change + 1] - ref2 - dy*0.5;
							coords[node_change + 2] = coords[node_change + 2] - ref3 - dz*0.5;
							coords[node_change + 1] = coords[node_change + 1]*cos(alpha) - coords[node_change + 2]*sin(alpha);
							coords[node_change + 2] = coords[node_change + 1]*sin(alpha) + coords[node_change + 2]*cos(alpha);
							coords[node_change + 1] = coords[node_change + 1] + ref2 + dy*0.5;
							coords[node_change + 2] = coords[node_change + 2] + ref3 + dz*0.5;
						}
						//move 1

					}else if(edge_list[1]&&edge_list[6]){
						int node_change = 3*elem->nodes[2].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[1], p0);
						d_c2 = gts_point_distance(point[6], p0);

						double alpha = atan(d_c1/d_c2);

						for(int co = 0; co<8;co++){
							node_change = 3*elem->nodes[co].id;
							coords[node_change + 1] = coords[node_change + 1] - ref2 - dy*0.5;
							coords[node_change + 2] = coords[node_change + 2] - ref3 - dz*0.5;
							coords[node_change + 1] = coords[node_change + 1]*cos(alpha) - coords[node_change + 2]*sin(alpha);
							coords[node_change + 2] = coords[node_change + 1]*sin(alpha) + coords[node_change + 2]*cos(alpha);
							coords[node_change + 1] = coords[node_change + 1] + ref2 + dy*0.5;
							coords[node_change + 2] = coords[node_change + 2] + ref3 + dz*0.5;
						}
						//move 2

					}else if(edge_list[2]&&edge_list[6]){
						int node_change = 3*elem->nodes[2].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[2], p0);
						d_c2 = gts_point_distance(point[6], p0);

						double alpha = atan(d_c1/d_c2);

						for(int co = 0; co<8;co++){
							node_change = 3*elem->nodes[co].id;
							coords[node_change    ] = coords[node_change    ] - ref1 - dx*0.5;
							coords[node_change + 2] = coords[node_change + 2] - ref3 - dz*0.5;
							coords[node_change    ] = coords[node_change    ]*cos(alpha) - coords[node_change + 2]*sin(alpha);
							coords[node_change + 2] = coords[node_change    ]*sin(alpha) + coords[node_change + 2]*cos(alpha);
							coords[node_change    ] = coords[node_change    ] + ref1 + dx*0.5;
							coords[node_change + 2] = coords[node_change + 2] + ref3 + dz*0.5;
						}
						//move 2

					}else if(edge_list[2]&&edge_list[7]){
						int node_change = 3*elem->nodes[3].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[2], p0);
						d_c2 = gts_point_distance(point[7], p0);

						double alpha = atan(d_c1/d_c2);

						for(int co = 0; co<8;co++){
							node_change = 3*elem->nodes[co].id;
							coords[node_change    ] = coords[node_change    ] - ref1 - dx*0.5;
							coords[node_change + 2] = coords[node_change + 2] - ref3 - dz*0.5;
							coords[node_change    ] = coords[node_change    ]*cos(alpha) - coords[node_change + 2]*sin(alpha);
							coords[node_change + 2] = coords[node_change    ]*sin(alpha) + coords[node_change + 2]*cos(alpha);
							coords[node_change    ] = coords[node_change    ] + ref1 + dx*0.5;
							coords[node_change + 2] = coords[node_change + 2] + ref3 + dz*0.5;
						}
						//move 3

					}else if(edge_list[3]&&edge_list[7]){
						int node_change = 3*elem->nodes[3].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[3], p0);
						d_c2 = gts_point_distance(point[7], p0);

						double alpha = atan(d_c1/d_c2);

						for(int co = 0; co<8;co++){
							node_change = 3*elem->nodes[co].id;
							coords[node_change + 1] = coords[node_change + 1] - ref2 - dy*0.5;
							coords[node_change + 2] = coords[node_change + 2] - ref3 - dz*0.5;
							coords[node_change + 1] = coords[node_change + 1]*cos(alpha) - coords[node_change + 2]*sin(alpha);
							coords[node_change + 2] = coords[node_change + 1]*sin(alpha) + coords[node_change + 2]*cos(alpha);
							coords[node_change + 1] = coords[node_change + 1] + ref2 + dy*0.5;
							coords[node_change + 2] = coords[node_change + 2] + ref3 + dz*0.5;
						}
						//move 3

					}else if(edge_list[8]&&edge_list[4]){
						int node_change = 3*elem->nodes[4].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[8], p0);
						d_c2 = gts_point_distance(point[4], p0);

						double alpha = atan(d_c1/d_c2);

						for(int co = 0; co<8;co++){
							node_change = 3*elem->nodes[co].id;
							coords[node_change    ] = coords[node_change    ] - ref1 - dx*0.5;
							coords[node_change + 2] = coords[node_change + 2] - ref3 - dz*0.5;
							coords[node_change    ] = coords[node_change    ]*cos(alpha) - coords[node_change + 2]*sin(alpha);
							coords[node_change + 2] = coords[node_change    ]*sin(alpha) + coords[node_change + 2]*cos(alpha);
							coords[node_change    ] = coords[node_change    ] + ref1 + dx*0.5;
							coords[node_change + 2] = coords[node_change + 2] + ref3 + dz*0.5;
						}
						//move 4

					}else if(edge_list[11]&&edge_list[4]){
						int node_change = 3*elem->nodes[4].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[11], p0);
						d_c2 = gts_point_distance(point[4], p0);

						double alpha = atan(d_c1/d_c2);

						for(int co = 0; co<8;co++){
							node_change = 3*elem->nodes[co].id;
							coords[node_change + 1] = coords[node_change + 1] - ref2 - dy*0.5;
							coords[node_change + 2] = coords[node_change + 2] - ref3 - dz*0.5;
							coords[node_change + 1] = coords[node_change + 1]*cos(alpha) - coords[node_change + 2]*sin(alpha);
							coords[node_change + 2] = coords[node_change + 1]*sin(alpha) + coords[node_change + 2]*cos(alpha);
							coords[node_change + 1] = coords[node_change + 1] + ref2 + dy*0.5;
							coords[node_change + 2] = coords[node_change + 2] + ref3 + dz*0.5;
						}
						//move 4

					}else if(edge_list[8]&&edge_list[5]){
						int node_change = 3*elem->nodes[5].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[8], p0);
						d_c2 = gts_point_distance(point[5], p0);

						double alpha = atan(d_c1/d_c2);

						for(int co = 0; co<8;co++){
							node_change = 3*elem->nodes[co].id;
							coords[node_change    ] = coords[node_change    ] - ref1 - dx*0.5;
							coords[node_change + 2] = coords[node_change + 2] - ref3 - dz*0.5;
							coords[node_change    ] = coords[node_change    ]*cos(alpha) - coords[node_change + 2]*sin(alpha);
							coords[node_change + 2] = coords[node_change    ]*sin(alpha) + coords[node_change + 2]*cos(alpha);
							coords[node_change    ] = coords[node_change    ] + ref1 + dx*0.5;
							coords[node_change + 2] = coords[node_change + 2] + ref3 + dz*0.5;
						}
						//move 5

					}else if(edge_list[9]&&edge_list[5]){
						int node_change = 3*elem->nodes[5].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[9], p0);
						d_c2 = gts_point_distance(point[5], p0);

						double alpha = atan(d_c1/d_c2);

						for(int co = 0; co<8;co++){
							node_change = 3*elem->nodes[co].id;
							coords[node_change + 1] = coords[node_change + 1] - ref2 - dy*0.5;
							coords[node_change + 2] = coords[node_change + 2] - ref3 - dz*0.5;
							coords[node_change + 1] = coords[node_change + 1]*cos(alpha) - coords[node_change + 2]*sin(alpha);
							coords[node_change + 2] = coords[node_change + 1]*sin(alpha) + coords[node_change + 2]*cos(alpha);
							coords[node_change + 1] = coords[node_change + 1] + ref2 + dy*0.5;
							coords[node_change + 2] = coords[node_change + 2] + ref3 + dz*0.5;
						}
						//move 5

					}else if(edge_list[9]&&edge_list[6]){
						int node_change = 3*elem->nodes[6].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[9], p0);
						d_c2 = gts_point_distance(point[6], p0);

						double alpha = atan(d_c1/d_c2);

						for(int co = 0; co<8;co++){
							node_change = 3*elem->nodes[co].id;
							coords[node_change + 1] = coords[node_change + 1] - ref2 - dy*0.5;
							coords[node_change + 2] = coords[node_change + 2] - ref3 - dz*0.5;
							coords[node_change + 1] = coords[node_change + 1]*cos(alpha) - coords[node_change + 2]*sin(alpha);
							coords[node_change + 2] = coords[node_change + 1]*sin(alpha) + coords[node_change + 2]*cos(alpha);
							coords[node_change + 1] = coords[node_change + 1] + ref2 + dy*0.5;
							coords[node_change + 2] = coords[node_change + 2] + ref3 + dz*0.5;
						}
						//move 6

					}else if(edge_list[10]&&edge_list[6]){
						int node_change = 3*elem->nodes[6].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[10], p0);
						d_c2 = gts_point_distance(point[6], p0);

						double alpha = atan(d_c1/d_c2);

						for(int co = 0; co<8;co++){
							node_change = 3*elem->nodes[co].id;
							coords[node_change    ] = coords[node_change    ] - ref1 - dx*0.5;
							coords[node_change + 2] = coords[node_change + 2] - ref3 - dz*0.5;
							coords[node_change    ] = coords[node_change    ]*cos(alpha) - coords[node_change + 2]*sin(alpha);
							coords[node_change + 2] = coords[node_change    ]*sin(alpha) + coords[node_change + 2]*cos(alpha);
							coords[node_change    ] = coords[node_change    ] + ref1 + dx*0.5;
							coords[node_change + 2] = coords[node_change + 2] + ref3 + dz*0.5;
						}
						//move 6

					}else if(edge_list[10]&&edge_list[7]){
						int node_change = 3*elem->nodes[7].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[10], p0);
						d_c2 = gts_point_distance(point[7], p0);

						double alpha = atan(d_c1/d_c2);

						for(int co = 0; co<8;co++){
							node_change = 3*elem->nodes[co].id;
							coords[node_change    ] = coords[node_change    ] - ref1 - dx*0.5;
							coords[node_change + 2] = coords[node_change + 2] - ref3 - dz*0.5;
							coords[node_change    ] = coords[node_change    ]*cos(alpha) - coords[node_change + 2]*sin(alpha);
							coords[node_change + 2] = coords[node_change    ]*sin(alpha) + coords[node_change + 2]*cos(alpha);
							coords[node_change    ] = coords[node_change    ] + ref1 + dx*0.5;
							coords[node_change + 2] = coords[node_change + 2] + ref3 + dz*0.5;
						}
						//move 7

					}else if(edge_list[11]&&edge_list[7]){
						int node_change = 3*elem->nodes[7].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[11], p0);
						d_c2 = gts_point_distance(point[7], p0);

						double alpha = atan(d_c1/d_c2);

						for(int co = 0; co<8;co++){
							node_change = 3*elem->nodes[co].id;
							coords[node_change + 1] = coords[node_change + 1] - ref2 - dy*0.5;
							coords[node_change + 2] = coords[node_change + 2] - ref3 - dz*0.5;
							coords[node_change + 1] = coords[node_change + 1]*cos(alpha) - coords[node_change + 2]*sin(alpha);
							coords[node_change + 2] = coords[node_change + 1]*sin(alpha) + coords[node_change + 2]*cos(alpha);
							coords[node_change + 1] = coords[node_change + 1] + ref2 + dy*0.5;
							coords[node_change + 2] = coords[node_change + 2] + ref3 + dz*0.5;
						}
						//move 7
					}


					CheckTemplate(mesh, coords, element_ids_local,false);
					elem = (octant_t*) sc_array_index(&mesh->elements, elements_ids[iel]);

					fprintf(fdbg,"Pad local:    ");
					for(int co = 0; co < element_ids_local.size(); co++){
						octant_t *h = (octant_t*) sc_array_index(&mesh->elements, element_ids_local[co]);
						var_aux_1[co]=h->pad;
						fprintf(fdbg,"%d ",h->pad);
					}
					fprintf(fdbg,"\n");

					if(elem->pad!=4 && elem->pad!=-10){
						flag2=false;
						for(int co = 0; co < element_ids_local.size(); co++){
							if(element_ids_local[co]!=elements_ids[iel]){
								if(var_aux[co]==var_aux_1[co]){
									flag2=false;
								}else{
									if(var_aux_1[co]==4){
										flag2=true;
										break;
									}else if(var_aux_1[co]==-10){
										flag2=true;
										break;
									}
								}
							}
						}
					}
				}

				if(!flag2){
					flag1=false;
					flag2=false;
					printf("To livreeeee rand tipo 4 rodei a bahiana El:%d, n_iter:%d\n",elements_ids[iel],n_iter);
					break;
				}

				if( flag2 ){
					for(int n_nodes = 0; n_nodes<8;n_nodes++){
						int node = 3*elem->nodes[n_nodes].id;
						coords[node]   = coords_orgi[3*n_nodes];
						coords[node+1] = coords_orgi[3*n_nodes+1];
						coords[node+2] = coords_orgi[3*n_nodes+2];
					}
				}


#endif

#if 0
				template2Rand(mesh, coords, elem, v_1);
				CheckTemplate(mesh, coords, element_ids_local,false);
				elem = (octant_t*) sc_array_index(&mesh->elements, elements_ids[iel]);

				fprintf(fdbg,"Pad local:    ");
				for(int co = 0; co < element_ids_local.size(); co++){
					octant_t *h = (octant_t*) sc_array_index(&mesh->elements, element_ids_local[co]);
					var_aux_1[co]=h->pad;
					fprintf(fdbg,"%d ",h->pad);
				}
				fprintf(fdbg,"\n");

				if(elem->pad!=4 && elem->pad!=-10){
					flag2=false;
					for(int co = 0; co < element_ids_local.size(); co++){
						if(element_ids_local[co]!=elements_ids[iel]){
							if(var_aux[co]==var_aux_1[co]){
								flag2=false;
							}else{
								if(var_aux_1[co]==4){
									flag2=true;
									break;
								}else if(var_aux_1[co]==-10){
									flag2=true;
									break;
								}
							}
						}
					}
				}

				if(!flag2){
					flag1=false;
					flag2=false;
					printf("To livreeeee rand El:%d, n_iter:%d\n",elements_ids[iel],n_iter);
					break;
				}

				if(elem->pad==4){
					//	CuttedEdges (mesh, coords, elements_ids, iel, &point[12], &edge_list[12]);

					GtsSegment * segments[12]={0};
					int Edge2GNode[12][2];
					GtsPoint * point[12] = {NULL};
					bool edge_list[12] = {false};

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

						edge_list[edge] = false;

						if (list == NULL) continue;
						while (list) {
							GtsBBox *b = GTS_BBOX(list->data);
							point[edge] = SegmentTriangleIntersection(segments[edge], GTS_TRIANGLE(b->bounded));
							if (point[edge]) {
								edge_list[edge] = true;
								break;
							}
							list = list->next;
						}
					}

					double d_c1;
					double d_c2;

					double ref1 = 	coords[3*elem->nodes[0].id ] ;
					double ref2 = 	coords[3*elem->nodes[0].id + 1 ] ;
					double ref3 = coords[3*elem->nodes[0].id + 3 ] ;
					double dx = coords[3*elem->nodes[1].id ] - coords[3*elem->nodes[0].id ];
					double dy = coords[3*elem->nodes[2].id +1] - coords[3*elem->nodes[1].id+1 ];
					double dz = coords[3*elem->nodes[0].id +2] - coords[3*elem->nodes[4].id+2 ];

					if(edge_list[0]&&edge_list[1]){
						int node_change = 3*elem->nodes[1].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[0], p0);
						d_c2 = gts_point_distance(point[1], p0);

						double alpha = atan(d_c2/d_c1);

						for(int co = 0; co<8;co++){
							node_change = 3*elem->nodes[co].id;
							coords[node_change  ] = coords[node_change   ] - ref1 - dx*0.5;
							coords[node_change +1  ] =coords[node_change +1  ] -ref2-dy/2;
							coords[node_change  ] = coords[node_change  ]*cos(alpha)  - coords[node_change+1]*sin(alpha);
							coords[node_change +1  ] = coords[node_change]*sin(alpha) +coords[node_change +1  ]*cos(alpha);
							coords[node_change  ] = coords[node_change   ] + ref1 + dx*0.5;
							coords[node_change +1  ] =coords[node_change +1  ] + ref2 + dy/2;
						}
						//move 1

					}else if(edge_list[1]&&edge_list[2]){
						int node_change = 3*elem->nodes[2].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[2], p0);
						d_c2 = gts_point_distance(point[1], p0);

						double alpha = atan(d_c2/d_c1);

						for(int co = 0; co<8;co++){
							node_change = 3*elem->nodes[co].id;
							coords[node_change  ] = coords[node_change   ] - ref1 - dx*0.5;
							coords[node_change +1  ] =coords[node_change +1  ] -ref2-dy/2;
							coords[node_change  ] = coords[node_change  ]*cos(alpha)  - coords[node_change+1]*sin(alpha);
							coords[node_change +1  ] = coords[node_change  ]*sin(alpha) +coords[node_change +1  ]*cos(alpha);
							coords[node_change  ] = coords[node_change   ] + ref1 + dx*0.5;
							coords[node_change +1  ] =coords[node_change +1  ] + ref2 + dy/2;
						}
						//move 2

					}else if(edge_list[2]&&edge_list[3]){
						int node_change = 3*elem->nodes[3].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[2], p0);
						d_c2 = gts_point_distance(point[3], p0);

						double alpha = atan(d_c2/d_c1);

						for(int co = 0; co<8;co++){
							node_change = 3*elem->nodes[co].id;
							coords[node_change  ] = coords[node_change   ] - ref1 - dx*0.5;
							coords[node_change +1  ] =coords[node_change +1  ] -ref2-dy/2;
							coords[node_change  ] = coords[node_change  ]*cos(alpha)  - coords[node_change+1]*sin(alpha);
							coords[node_change +1  ] = coords[node_change  ]*sin(alpha) +coords[node_change +1  ]*cos(alpha);
							coords[node_change  ] = coords[node_change   ] + ref1 + dx*0.5;
							coords[node_change +1  ] =coords[node_change +1  ] + ref2 + dy/2;
						}
						//move 3

					}else if(edge_list[3]&&edge_list[0]){
						int node_change = 3*elem->nodes[0].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[0], p0);
						d_c2 = gts_point_distance(point[3], p0);

						double alpha = atan(d_c2/d_c1);

						for(int co = 0; co<8;co++){
							node_change = 3*elem->nodes[co].id;
							coords[node_change  ] = coords[node_change   ] - ref1 - dx*0.5;
							coords[node_change +1  ] =coords[node_change +1  ] -ref2-dy/2;
							coords[node_change  ] = coords[node_change  ]*cos(alpha)  - coords[node_change+1]*sin(alpha);
							coords[node_change +1  ] = coords[node_change  ]*sin(alpha) +coords[node_change +1  ]*cos(alpha);
							coords[node_change  ] = coords[node_change   ] + ref1 + dx*0.5;
							coords[node_change +1  ] =coords[node_change +1  ] + ref2 + dy/2;
						}
						//move 0

					}else if(edge_list[8]&&edge_list[9]){
						int node_change = 3*elem->nodes[5].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[8], p0);
						d_c2 = gts_point_distance(point[9], p0);

						double alpha = atan(d_c2/d_c1);

						for(int co = 0; co<8;co++){
							node_change = 3*elem->nodes[co].id;
							coords[node_change  ] = coords[node_change   ] - ref1 - dx*0.5;
							coords[node_change +1  ] =coords[node_change +1  ] -ref2-dy/2;
							coords[node_change  ] = coords[node_change  ]*cos(alpha)  - coords[node_change+1]*sin(alpha);
							coords[node_change +1  ] = coords[node_change]*sin(alpha) +coords[node_change +1  ]*cos(alpha);
							coords[node_change  ] = coords[node_change   ] + ref1 + dx*0.5;
							coords[node_change +1  ] =coords[node_change +1  ] + ref2 + dy/2;
						}
						//move 5

					}else if(edge_list[9]&&edge_list[10]){
						int node_change = 3*elem->nodes[6].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[10], p0);
						d_c2 = gts_point_distance(point[9], p0);

						double alpha = atan(d_c2/d_c1);

						for(int co = 0; co<8;co++){
							node_change = 3*elem->nodes[co].id;
							coords[node_change  ] = coords[node_change   ] - ref1 - dx*0.5;
							coords[node_change +1  ] =coords[node_change +1  ] -ref2-dy/2;
							coords[node_change  ] = coords[node_change  ]*cos(alpha)  - coords[node_change+1]*sin(alpha);
							coords[node_change +1  ] = coords[node_change ]*sin(alpha) +coords[node_change +1  ]*cos(alpha);
							coords[node_change  ] = coords[node_change   ] + ref1 + dx*0.5;
							coords[node_change +1  ] =coords[node_change +1  ] + ref2 + dy/2;
						}
						//move 6

					}else if(edge_list[10]&&edge_list[11]){
						int node_change = 3*elem->nodes[7].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[10], p0);
						d_c2 = gts_point_distance(point[11], p0);

						double alpha = atan(d_c2/d_c1);

						for(int co = 0; co<8;co++){
							node_change = 3*elem->nodes[co].id;
							coords[node_change  ] = coords[node_change   ] - ref1 - dx*0.5;
							coords[node_change +1  ] =coords[node_change +1  ] -ref2-dy/2;
							coords[node_change  ] = coords[node_change  ]*cos(alpha)  - coords[node_change+1]*sin(alpha);
							coords[node_change +1  ] = coords[node_change]*sin(alpha) +coords[node_change +1  ]*cos(alpha);
							coords[node_change  ] = coords[node_change   ] + ref1 + dx*0.5;
							coords[node_change +1  ] =coords[node_change +1  ] + ref2 + dy/2;
						}
						//move 7

					}else if(edge_list[11]&&edge_list[8]){
						int node_change = 3*elem->nodes[4].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[8], p0);
						d_c2 = gts_point_distance(point[11], p0);

						double alpha = atan(d_c2/d_c1);

						for(int co = 0; co<8;co++){
							node_change = 3*elem->nodes[co].id;
							coords[node_change  ] = coords[node_change   ] - ref1 - dx*0.5;
							coords[node_change +1  ] =coords[node_change +1  ] -ref2-dy/2;
							coords[node_change  ] = coords[node_change  ]*cos(alpha)  - coords[node_change+1]*sin(alpha);
							coords[node_change +1  ] = coords[node_change ]*sin(alpha) +coords[node_change +1  ]*cos(alpha);
							coords[node_change  ] = coords[node_change   ] + ref1 + dx*0.5;
							coords[node_change +1  ] =coords[node_change +1  ] + ref2 + dy/2;
						}
						//move 4

					}else if( edge_list[0]&&edge_list[4] ){
						int node_change = 3*elem->nodes[0].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[0], p0);
						d_c2 = gts_point_distance(point[4], p0);

						double alpha = atan(d_c2/d_c1);

						for(int co = 0; co<8;co++){
							node_change = 3*elem->nodes[co].id;
							coords[node_change  ] = coords[node_change   ] - ref1 - dx*0.5;
							coords[node_change +2  ] =coords[node_change +2  ] -ref3-dz/2;
							coords[node_change  ] = coords[node_change  ]*cos(alpha)  - coords[node_change+2]*sin(alpha);
							coords[node_change +2  ] = coords[node_change]*sin(alpha) +coords[node_change +2  ]*cos(alpha);
							coords[node_change  ] = coords[node_change   ] + ref1 + dx*0.5;
							coords[node_change +2  ] =coords[node_change +1  ] + ref3 + dz/2;
						}
						//move 0

					}else if( edge_list[3]&&edge_list[4] ){
						int node_change = 3*elem->nodes[0].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[3], p0);
						d_c2 = gts_point_distance(point[4], p0);

						double alpha = atan(d_c2/d_c1);

						for(int co = 0; co<8;co++){
							node_change = 3*elem->nodes[co].id;
							coords[node_change + 1] = coords[node_change + 1] - ref2 - dy*0.5;
							coords[node_change + 2] = coords[node_change + 2] - ref3 - dz*0.5;
							coords[node_change + 1] = coords[node_change + 1]*cos(alpha) - coords[node_change + 2]*sin(alpha);
							coords[node_change + 2] = coords[node_change + 1]*sin(alpha) + coords[node_change + 2]*cos(alpha);
							coords[node_change + 1] = coords[node_change + 1] + ref2 + dy*0.5;
							coords[node_change + 2] = coords[node_change + 2] + ref3 + dz*0.5;
						}
						//move 0

					}else if(edge_list[0]&&edge_list[5]){
						int node_change = 3*elem->nodes[1].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[0], p0);
						d_c2 = gts_point_distance(point[5], p0);

						double alpha = atan(d_c2/d_c1);

						for(int co = 0; co<8;co++){
							node_change = 3*elem->nodes[co].id;
							coords[node_change    ] = coords[node_change    ] - ref1 - dx*0.5;
							coords[node_change + 2] = coords[node_change + 2] - ref3 - dz*0.5;
							coords[node_change    ] = coords[node_change    ]*cos(alpha) - coords[node_change + 2]*sin(alpha);
							coords[node_change + 2] = coords[node_change    ]*sin(alpha) + coords[node_change + 2]*cos(alpha);
							coords[node_change    ] = coords[node_change    ] + ref1 + dx*0.5;
							coords[node_change + 2] = coords[node_change + 2] + ref3 + dz*0.5;
						}
						//move 1

					}else if(edge_list[1]&&edge_list[5]){
						int node_change = 3*elem->nodes[1].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[1], p0);
						d_c2 = gts_point_distance(point[5], p0);

						double alpha = atan(d_c2/d_c1);

						for(int co = 0; co<8;co++){
							node_change = 3*elem->nodes[co].id;
							coords[node_change + 1] = coords[node_change + 1] - ref2 - dy*0.5;
							coords[node_change + 2] = coords[node_change + 2] - ref3 - dz*0.5;
							coords[node_change + 1] = coords[node_change + 1]*cos(alpha) - coords[node_change + 2]*sin(alpha);
							coords[node_change + 2] = coords[node_change + 1]*sin(alpha) + coords[node_change + 2]*cos(alpha);
							coords[node_change + 1] = coords[node_change + 1] + ref2 + dy*0.5;
							coords[node_change + 2] = coords[node_change + 2] + ref3 + dz*0.5;
						}
						//move 1

					}else if(edge_list[1]&&edge_list[6]){
						int node_change = 3*elem->nodes[2].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[1], p0);
						d_c2 = gts_point_distance(point[6], p0);

						double alpha = atan(d_c2/d_c1);

						for(int co = 0; co<8;co++){
							node_change = 3*elem->nodes[co].id;
							coords[node_change + 1] = coords[node_change + 1] - ref2 - dy*0.5;
							coords[node_change + 2] = coords[node_change + 2] - ref3 - dz*0.5;
							coords[node_change + 1] = coords[node_change + 1]*cos(alpha) - coords[node_change + 2]*sin(alpha);
							coords[node_change + 2] = coords[node_change + 1]*sin(alpha) + coords[node_change + 2]*cos(alpha);
							coords[node_change + 1] = coords[node_change + 1] + ref2 + dy*0.5;
							coords[node_change + 2] = coords[node_change + 2] + ref3 + dz*0.5;
						}
						//move 2

					}else if(edge_list[2]&&edge_list[6]){
						int node_change = 3*elem->nodes[2].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[2], p0);
						d_c2 = gts_point_distance(point[6], p0);

						double alpha = atan(d_c2/d_c1);

						for(int co = 0; co<8;co++){
							node_change = 3*elem->nodes[co].id;
							coords[node_change    ] = coords[node_change    ] - ref1 - dx*0.5;
							coords[node_change + 2] = coords[node_change + 2] - ref3 - dz*0.5;
							coords[node_change    ] = coords[node_change    ]*cos(alpha) - coords[node_change + 2]*sin(alpha);
							coords[node_change + 2] = coords[node_change    ]*sin(alpha) + coords[node_change + 2]*cos(alpha);
							coords[node_change    ] = coords[node_change    ] + ref1 + dx*0.5;
							coords[node_change + 2] = coords[node_change + 2] + ref3 + dz*0.5;
						}
						//move 2

					}else if(edge_list[2]&&edge_list[7]){
						int node_change = 3*elem->nodes[3].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[2], p0);
						d_c2 = gts_point_distance(point[7], p0);

						double alpha = atan(d_c2/d_c1);

						for(int co = 0; co<8;co++){
							node_change = 3*elem->nodes[co].id;
							coords[node_change    ] = coords[node_change    ] - ref1 - dx*0.5;
							coords[node_change + 2] = coords[node_change + 2] - ref3 - dz*0.5;
							coords[node_change    ] = coords[node_change    ]*cos(alpha) - coords[node_change + 2]*sin(alpha);
							coords[node_change + 2] = coords[node_change    ]*sin(alpha) + coords[node_change + 2]*cos(alpha);
							coords[node_change    ] = coords[node_change    ] + ref1 + dx*0.5;
							coords[node_change + 2] = coords[node_change + 2] + ref3 + dz*0.5;
						}
						//move 3

					}else if(edge_list[3]&&edge_list[7]){
						int node_change = 3*elem->nodes[3].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[3], p0);
						d_c2 = gts_point_distance(point[7], p0);

						double alpha = atan(d_c2/d_c1);

						for(int co = 0; co<8;co++){
							node_change = 3*elem->nodes[co].id;
							coords[node_change + 1] = coords[node_change + 1] - ref2 - dy*0.5;
							coords[node_change + 2] = coords[node_change + 2] - ref3 - dz*0.5;
							coords[node_change + 1] = coords[node_change + 1]*cos(alpha) - coords[node_change + 2]*sin(alpha);
							coords[node_change + 2] = coords[node_change + 1]*sin(alpha) + coords[node_change + 2]*cos(alpha);
							coords[node_change + 1] = coords[node_change + 1] + ref2 + dy*0.5;
							coords[node_change + 2] = coords[node_change + 2] + ref3 + dz*0.5;
						}
						//move 3

					}else if(edge_list[8]&&edge_list[4]){
						int node_change = 3*elem->nodes[4].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[8], p0);
						d_c2 = gts_point_distance(point[4], p0);

						double alpha = atan(d_c2/d_c1);

						for(int co = 0; co<8;co++){
							node_change = 3*elem->nodes[co].id;
							coords[node_change    ] = coords[node_change    ] - ref1 - dx*0.5;
							coords[node_change + 2] = coords[node_change + 2] - ref3 - dz*0.5;
							coords[node_change    ] = coords[node_change    ]*cos(alpha) - coords[node_change + 2]*sin(alpha);
							coords[node_change + 2] = coords[node_change    ]*sin(alpha) + coords[node_change + 2]*cos(alpha);
							coords[node_change    ] = coords[node_change    ] + ref1 + dx*0.5;
							coords[node_change + 2] = coords[node_change + 2] + ref3 + dz*0.5;
						}
						//move 4

					}else if(edge_list[11]&&edge_list[4]){
						int node_change = 3*elem->nodes[4].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[11], p0);
						d_c2 = gts_point_distance(point[4], p0);

						double alpha = atan(d_c2/d_c1);

						for(int co = 0; co<8;co++){
							node_change = 3*elem->nodes[co].id;
							coords[node_change + 1] = coords[node_change + 1] - ref2 - dy*0.5;
							coords[node_change + 2] = coords[node_change + 2] - ref3 - dz*0.5;
							coords[node_change + 1] = coords[node_change + 1]*cos(alpha) - coords[node_change + 2]*sin(alpha);
							coords[node_change + 2] = coords[node_change + 1]*sin(alpha) + coords[node_change + 2]*cos(alpha);
							coords[node_change + 1] = coords[node_change + 1] + ref2 + dy*0.5;
							coords[node_change + 2] = coords[node_change + 2] + ref3 + dz*0.5;
						}
						//move 4

					}else if(edge_list[8]&&edge_list[5]){
						int node_change = 3*elem->nodes[5].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[8], p0);
						d_c2 = gts_point_distance(point[5], p0);

						double alpha = atan(d_c2/d_c1);

						for(int co = 0; co<8;co++){
							node_change = 3*elem->nodes[co].id;
							coords[node_change    ] = coords[node_change    ] - ref1 - dx*0.5;
							coords[node_change + 2] = coords[node_change + 2] - ref3 - dz*0.5;
							coords[node_change    ] = coords[node_change    ]*cos(alpha) - coords[node_change + 2]*sin(alpha);
							coords[node_change + 2] = coords[node_change    ]*sin(alpha) + coords[node_change + 2]*cos(alpha);
							coords[node_change    ] = coords[node_change    ] + ref1 + dx*0.5;
							coords[node_change + 2] = coords[node_change + 2] + ref3 + dz*0.5;
						}
						//move 5

					}else if(edge_list[9]&&edge_list[5]){
						int node_change = 3*elem->nodes[5].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[9], p0);
						d_c2 = gts_point_distance(point[5], p0);

						double alpha = atan(d_c2/d_c1);

						for(int co = 0; co<8;co++){
							node_change = 3*elem->nodes[co].id;
							coords[node_change + 1] = coords[node_change + 1] - ref2 - dy*0.5;
							coords[node_change + 2] = coords[node_change + 2] - ref3 - dz*0.5;
							coords[node_change + 1] = coords[node_change + 1]*cos(alpha) - coords[node_change + 2]*sin(alpha);
							coords[node_change + 2] = coords[node_change + 1]*sin(alpha) + coords[node_change + 2]*cos(alpha);
							coords[node_change + 1] = coords[node_change + 1] + ref2 + dy*0.5;
							coords[node_change + 2] = coords[node_change + 2] + ref3 + dz*0.5;
						}
						//move 5

					}else if(edge_list[9]&&edge_list[6]){
						int node_change = 3*elem->nodes[6].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[9], p0);
						d_c2 = gts_point_distance(point[6], p0);

						double alpha = atan(d_c2/d_c1);

						for(int co = 0; co<8;co++){
							node_change = 3*elem->nodes[co].id;
							coords[node_change + 1] = coords[node_change + 1] - ref2 - dy*0.5;
							coords[node_change + 2] = coords[node_change + 2] - ref3 - dz*0.5;
							coords[node_change + 1] = coords[node_change + 1]*cos(alpha) - coords[node_change + 2]*sin(alpha);
							coords[node_change + 2] = coords[node_change + 1]*sin(alpha) + coords[node_change + 2]*cos(alpha);
							coords[node_change + 1] = coords[node_change + 1] + ref2 + dy*0.5;
							coords[node_change + 2] = coords[node_change + 2] + ref3 + dz*0.5;
						}
						//move 6

					}else if(edge_list[10]&&edge_list[6]){
						int node_change = 3*elem->nodes[6].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[10], p0);
						d_c2 = gts_point_distance(point[6], p0);

						double alpha = atan(d_c2/d_c1);

						for(int co = 0; co<8;co++){
							node_change = 3*elem->nodes[co].id;
							coords[node_change    ] = coords[node_change    ] - ref1 - dx*0.5;
							coords[node_change + 2] = coords[node_change + 2] - ref3 - dz*0.5;
							coords[node_change    ] = coords[node_change    ]*cos(alpha) - coords[node_change + 2]*sin(alpha);
							coords[node_change + 2] = coords[node_change    ]*sin(alpha) + coords[node_change + 2]*cos(alpha);
							coords[node_change    ] = coords[node_change    ] + ref1 + dx*0.5;
							coords[node_change + 2] = coords[node_change + 2] + ref3 + dz*0.5;
						}
						//move 6

					}else if(edge_list[10]&&edge_list[7]){
						int node_change = 3*elem->nodes[7].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[10], p0);
						d_c2 = gts_point_distance(point[7], p0);

						double alpha = atan(d_c2/d_c1);

						for(int co = 0; co<8;co++){
							node_change = 3*elem->nodes[co].id;
							coords[node_change    ] = coords[node_change    ] - ref1 - dx*0.5;
							coords[node_change + 2] = coords[node_change + 2] - ref3 - dz*0.5;
							coords[node_change    ] = coords[node_change    ]*cos(alpha) - coords[node_change + 2]*sin(alpha);
							coords[node_change + 2] = coords[node_change    ]*sin(alpha) + coords[node_change + 2]*cos(alpha);
							coords[node_change    ] = coords[node_change    ] + ref1 + dx*0.5;
							coords[node_change + 2] = coords[node_change + 2] + ref3 + dz*0.5;
						}
						//move 7

					}else if(edge_list[11]&&edge_list[7]){
						int node_change = 3*elem->nodes[7].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[11], p0);
						d_c2 = gts_point_distance(point[7], p0);

						double alpha = atan(d_c2/d_c1);

						for(int co = 0; co<8;co++){
							node_change = 3*elem->nodes[co].id;
							coords[node_change + 1] = coords[node_change + 1] - ref2 - dy*0.5;
							coords[node_change + 2] = coords[node_change + 2] - ref3 - dz*0.5;
							coords[node_change + 1] = coords[node_change + 1]*cos(alpha) - coords[node_change + 2]*sin(alpha);
							coords[node_change + 2] = coords[node_change + 1]*sin(alpha) + coords[node_change + 2]*cos(alpha);
							coords[node_change + 1] = coords[node_change + 1] + ref2 + dy*0.5;
							coords[node_change + 2] = coords[node_change + 2] + ref3 + dz*0.5;
						}
						//move 7
					}


					CheckTemplate(mesh, coords, element_ids_local,false);
					elem = (octant_t*) sc_array_index(&mesh->elements, elements_ids[iel]);

					fprintf(fdbg,"Pad local:    ");
					for(int co = 0; co < element_ids_local.size(); co++){
						octant_t *h = (octant_t*) sc_array_index(&mesh->elements, element_ids_local[co]);
						var_aux_1[co]=h->pad;
						fprintf(fdbg,"%d ",h->pad);
					}
					fprintf(fdbg,"\n");

					if(elem->pad!=4 && elem->pad!=-10){
						flag2=false;
						for(int co = 0; co < element_ids_local.size(); co++){
							if(element_ids_local[co]!=elements_ids[iel]){
								if(var_aux[co]==var_aux_1[co]){
									flag2=false;
								}else{
									if(var_aux_1[co]==4){
										flag2=true;
										break;
									}else if(var_aux_1[co]==-10){
										flag2=true;
										break;
									}
								}
							}
						}
					}
				}

				if(!flag2){
					flag1=false;
					flag2=false;
					printf("To livreeeee rand tipo 4 rodei a bahiana El:%d, n_iter:%d\n",elements_ids[iel],n_iter);
					break;
				}

				if( flag2 ){
					for(int n_nodes = 0; n_nodes<8;n_nodes++){
						int node = 3*elem->nodes[n_nodes].id;
						coords[node]   = coords_orgi[3*n_nodes];
						coords[node+1] = coords_orgi[3*n_nodes+1];
						coords[node+2] = coords_orgi[3*n_nodes+2];
					}
				}


#endif

#if 0
				template2Rand(mesh, coords, elem, v_1);
				CheckTemplate(mesh, coords, element_ids_local,false);
				elem = (octant_t*) sc_array_index(&mesh->elements, elements_ids[iel]);

				fprintf(fdbg,"Pad local:    ");
				for(int co = 0; co < element_ids_local.size(); co++){
					octant_t *h = (octant_t*) sc_array_index(&mesh->elements, element_ids_local[co]);
					var_aux_1[co]=h->pad;
					fprintf(fdbg,"%d ",h->pad);
				}
				fprintf(fdbg,"\n");

				if(elem->pad!=4 && elem->pad!=-10){
					flag2=false;
					for(int co = 0; co < element_ids_local.size(); co++){
						if(element_ids_local[co]!=elements_ids[iel]){
							if(var_aux[co]==var_aux_1[co]){
								flag2=false;
							}else{
								if(var_aux_1[co]==4){
									flag2=true;
									break;
								}else if(var_aux_1[co]==-10){
									flag2=true;
									break;
								}
							}
						}
					}
				}

				if(!flag2){
					flag1=false;
					flag2=false;
					printf("To livreeeee rand El:%d, n_iter:%d\n",elements_ids[iel],n_iter);
					break;
				}

				if(elem->pad==4){
					//	CuttedEdges (mesh, coords, elements_ids, iel, &point[12], &edge_list[12]);

					GtsSegment * segments[12]={0};
					int Edge2GNode[12][2];
					GtsPoint * point[12] = {NULL};
					bool edge_list[12] = {false};

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

						edge_list[edge] = false;

						if (list == NULL) continue;
						while (list) {
							GtsBBox *b = GTS_BBOX(list->data);
							point[edge] = SegmentTriangleIntersection(segments[edge], GTS_TRIANGLE(b->bounded));
							if (point[edge]) {
								edge_list[edge] = true;
								break;
							}
							list = list->next;
						}
					}

					double d_c1;
					double d_c2;

					double dx = coords[3*elem->nodes[1].id    ] - coords[3*elem->nodes[0].id    ];
					double dy = coords[3*elem->nodes[2].id + 1] - coords[3*elem->nodes[1].id + 1];
					double dz = coords[3*elem->nodes[0].id + 2] - coords[3*elem->nodes[4].id + 2];

					if(edge_list[0]&&edge_list[1]){
						int node_change = 3*elem->nodes[1].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[0], p0);
						d_c2 = gts_point_distance(point[1], p0);

						if(d_c1>d_c2){
							for(int co = 0; co<4;co++){
								node_change = 3*elem->nodes[FaceNodesMap[5][co]].id;
								if(dx*0.5>d_c1){
									coords[node_change  ] = coords[node_change  ]-(d_c1-dx*1.1);
								}else{
									coords[node_change  ] = coords[node_change  ]+(d_c1-dx*1.1);
								}
							}
						}else{
							for(int co = 0; co<4;co++){
								node_change = 3*elem->nodes[FaceNodesMap[5][co]].id;
								if(dy*0.5>d_c2){
									coords[node_change +1  ] = coords[node_change +1  ] - (d_c2-dy*1.1);
								}else{
									coords[node_change +1  ] = coords[node_change +1  ] + (d_c2-dy*1.1);
								}
							}
						}
						//move 1

					}else if(edge_list[1]&&edge_list[2]){
						int node_change = 3*elem->nodes[2].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[2], p0);
						d_c2 = gts_point_distance(point[1], p0);

						if(d_c1>d_c2){
							for(int co = 0; co<4;co++){
								node_change = 3*elem->nodes[FaceNodesMap[5][co]].id;
								if(dx*0.5>d_c1){
									coords[node_change  ] = coords[node_change  ]-(d_c1-dx*1.1);
								}else{
									coords[node_change  ] = coords[node_change  ]+(d_c1-dx*1.1);
								}
							}
						}else{
							for(int co = 0; co<4;co++){
								node_change = 3*elem->nodes[FaceNodesMap[5][co]].id;
								if(dy*0.5>d_c2){
									coords[node_change +1  ] = coords[node_change +1  ] - (d_c2-dy*1.1);
								}else{
									coords[node_change +1  ] = coords[node_change +1  ] + (d_c2-dy*1.1);
								}
							}
						}
						//move 2

					}else if(edge_list[2]&&edge_list[3]){
						int node_change = 3*elem->nodes[3].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[2], p0);
						d_c2 = gts_point_distance(point[3], p0);

						if(d_c1>d_c2){
							for(int co = 0; co<4;co++){
								node_change = 3*elem->nodes[FaceNodesMap[5][co]].id;
								if(dx*0.5>d_c1){
									coords[node_change  ] = coords[node_change  ]-(d_c1-dx*1.1);
								}else{
									coords[node_change  ] = coords[node_change  ]+(d_c1-dx*1.1);
								}
							}
						}else{
							for(int co = 0; co<4;co++){
								node_change = 3*elem->nodes[FaceNodesMap[5][co]].id;
								if(dy*0.5>d_c2){
									coords[node_change +1  ] = coords[node_change +1  ] - (d_c2-dy*1.1);
								}else{
									coords[node_change +1  ] = coords[node_change +1  ] + (d_c2-dy*1.1);
								}
							}
						}
						//move 3

					}else if(edge_list[3]&&edge_list[0]){
						int node_change = 3*elem->nodes[0].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[0], p0);
						d_c2 = gts_point_distance(point[3], p0);

						if(d_c1>d_c2){
							for(int co = 0; co<4;co++){
								node_change = 3*elem->nodes[FaceNodesMap[5][co]].id;
								if(dx*0.5>d_c1){
									coords[node_change  ] = coords[node_change  ]-(d_c1-dx*1.1);
								}else{
									coords[node_change  ] = coords[node_change  ]+(d_c1-dx*1.1);
								}
							}
						}else{
							for(int co = 0; co<4;co++){
								node_change = 3*elem->nodes[FaceNodesMap[5][co]].id;
								if(dy*0.5>d_c2){
									coords[node_change +1  ] = coords[node_change +1  ] - (d_c2-dy*1.1);
								}else{
									coords[node_change +1  ] = coords[node_change +1  ] + (d_c2-dy*1.1);
								}
							}
						}
						//move 0

					}else if(edge_list[8]&&edge_list[9]){
						int node_change = 3*elem->nodes[5].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[8], p0);
						d_c2 = gts_point_distance(point[9], p0);

						if(d_c1>d_c2){
							for(int co = 0; co<4;co++){
								node_change = 3*elem->nodes[FaceNodesMap[4][co]].id;
								if(dx*0.5>d_c1){
									coords[node_change  ] = coords[node_change  ]-(d_c1-dx*1.1);
								}else{
									coords[node_change  ] = coords[node_change  ]+(d_c1-dx*1.1);
								}
							}
						}else{
							for(int co = 0; co<4;co++){
								node_change = 3*elem->nodes[FaceNodesMap[4][co]].id;
								if(dy*0.5>d_c2){
									coords[node_change +1  ] = coords[node_change +1  ] - (d_c2-dy*1.1);
								}else{
									coords[node_change +1  ] = coords[node_change +1  ] + (d_c2-dy*1.1);
								}
							}
						}
						//move 5

					}else if(edge_list[9]&&edge_list[10]){
						int node_change = 3*elem->nodes[6].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[10], p0);
						d_c2 = gts_point_distance(point[9], p0);

						if(d_c1>d_c2){
							for(int co = 0; co<4;co++){
								node_change = 3*elem->nodes[FaceNodesMap[4][co]].id;
								if(dx*0.5>d_c1){
									coords[node_change  ] = coords[node_change  ]-(d_c1-dx*1.1);
								}else{
									coords[node_change  ] = coords[node_change  ]+(d_c1-dx*1.1);
								}
							}
						}else{
							for(int co = 0; co<4;co++){
								node_change = 3*elem->nodes[FaceNodesMap[4][co]].id;
								if(dy*0.5>d_c2){
									coords[node_change +1  ] = coords[node_change +1  ] - (d_c2-dy*1.1);
								}else{
									coords[node_change +1  ] = coords[node_change +1  ] + (d_c2-dy*1.1);
								}
							}
						}
						//move 6

					}else if(edge_list[10]&&edge_list[11]){
						int node_change = 3*elem->nodes[7].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[10], p0);
						d_c2 = gts_point_distance(point[11], p0);

						if(d_c1>d_c2){
							for(int co = 0; co<4;co++){
								node_change = 3*elem->nodes[FaceNodesMap[4][co]].id;
								if(dx*0.5>d_c1){
									coords[node_change  ] = coords[node_change  ]-(d_c1-dx*1.1);
								}else{
									coords[node_change  ] = coords[node_change  ]+(d_c1-dx*1.1);
								}
							}
						}else{
							for(int co = 0; co<4;co++){
								node_change = 3*elem->nodes[FaceNodesMap[4][co]].id;
								if(dy*0.5>d_c2){
									coords[node_change +1  ] = coords[node_change +1  ] - (d_c2-dy*1.1);
								}else{
									coords[node_change +1  ] = coords[node_change +1  ] + (d_c2-dy*1.1);
								}
							}
						}
						//move 7

					}else if(edge_list[11]&&edge_list[8]){
						int node_change = 3*elem->nodes[4].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[8], p0);
						d_c2 = gts_point_distance(point[11], p0);

						if(d_c1>d_c2){
							for(int co = 0; co<4;co++){
								node_change = 3*elem->nodes[FaceNodesMap[4][co]].id;
								if(dx*0.5>d_c1){
									coords[node_change  ] = coords[node_change  ]-(d_c1-dx*1.1);
								}else{
									coords[node_change  ] = coords[node_change  ]+(d_c1-dx*1.1);
								}
							}
						}else{
							for(int co = 0; co<4;co++){
								node_change = 3*elem->nodes[FaceNodesMap[4][co]].id;
								if(dy*0.5>d_c2){
									coords[node_change + 1] = coords[node_change +1  ] - (d_c2-dy*1.1);
								}else{
									coords[node_change + 1] = coords[node_change +1  ] + (d_c2-dy*1.1);
								}
							}
						}
						//move 4

					}else if( edge_list[0]&&edge_list[4] ){
						int node_change = 3*elem->nodes[0].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[0], p0);
						d_c2 = gts_point_distance(point[4], p0);

						if(d_c1>d_c2){
							for(int co = 0; co<4;co++){
								node_change = 3*elem->nodes[FaceNodesMap[2][co]].id;
								if(dx*0.5>d_c1){
									coords[node_change  ] = coords[node_change  ]-(d_c1-dx*1.1);
								}else{
									coords[node_change  ] = coords[node_change  ]+(d_c1-dx*1.1);
								}
							}
						}else{
							for(int co = 0; co<4;co++){
								node_change = 3*elem->nodes[FaceNodesMap[2][co]].id;
								if(dz*0.5>d_c2){
									coords[node_change + 2] = coords[node_change +2  ] - (d_c2-dz*1.1);
								}else{
									coords[node_change + 2] = coords[node_change +2  ] + (d_c2-dz*1.1);
								}
							}
						}
						//move 0

					}else if( edge_list[3]&&edge_list[4] ){
						int node_change = 3*elem->nodes[0].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[3], p0);
						d_c2 = gts_point_distance(point[4], p0);

						if(d_c1>d_c2){
							for(int co = 0; co<4;co++){
								node_change = 3*elem->nodes[FaceNodesMap[0][co]].id;
								if(dy*0.5>d_c1){
									coords[node_change + 1] = coords[node_change + 1]-(d_c1-dy*1.1);
								}else{
									coords[node_change + 1] = coords[node_change + 1]+(d_c1-dy*1.1);
								}
							}
						}else{
							for(int co = 0; co<4;co++){
								node_change = 3*elem->nodes[FaceNodesMap[0][co]].id;
								if(dz*0.5>d_c2){
									coords[node_change + 2] = coords[node_change +2  ] - (d_c2-dz*1.1);
								}else{
									coords[node_change + 2] = coords[node_change +2  ] + (d_c2-dz*1.1);
								}
							}
						}

						//move 0

					}else if(edge_list[0]&&edge_list[5]){
						int node_change = 3*elem->nodes[1].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[0], p0);
						d_c2 = gts_point_distance(point[5], p0);

						if(d_c1>d_c2){
							for(int co = 0; co<4;co++){
								node_change = 3*elem->nodes[FaceNodesMap[2][co]].id;
								if(dx*0.5>d_c1){
									coords[node_change  ] = coords[node_change  ]-(d_c1-dx*1.1);
								}else{
									coords[node_change  ] = coords[node_change  ]+(d_c1-dx*1.1);
								}
							}
						}else{
							for(int co = 0; co<4;co++){
								node_change = 3*elem->nodes[FaceNodesMap[2][co]].id;
								if(dz*0.5>d_c2){
									coords[node_change + 2] = coords[node_change +2  ] - (d_c2-dz*1.1);
								}else{
									coords[node_change + 2] = coords[node_change +2  ] + (d_c2-dz*1.1);
								}
							}
						}
						//move 1

					}else if(edge_list[1]&&edge_list[5]){
						int node_change = 3*elem->nodes[1].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[1], p0);
						d_c2 = gts_point_distance(point[5], p0);

						if(d_c1>d_c2){
							for(int co = 0; co<4;co++){
								node_change = 3*elem->nodes[FaceNodesMap[1][co]].id;
								if(dy*0.5>d_c1){
									coords[node_change + 1] = coords[node_change + 1]-(d_c1-dy*1.1);
								}else{
									coords[node_change + 1] = coords[node_change + 1]+(d_c1-dy*1.1);
								}
							}
						}else{
							for(int co = 0; co<4;co++){
								node_change = 3*elem->nodes[FaceNodesMap[1][co]].id;
								if(dz*0.5>d_c2){
									coords[node_change + 2] = coords[node_change +2  ] - (d_c2-dz*1.1);
								}else{
									coords[node_change + 2] = coords[node_change +2  ] + (d_c2-dz*1.1);
								}
							}
						}
						//move 1

					}else if(edge_list[1]&&edge_list[6]){
						int node_change = 3*elem->nodes[2].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[1], p0);
						d_c2 = gts_point_distance(point[6], p0);

						if(d_c1>d_c2){
							for(int co = 0; co<4;co++){
								node_change = 3*elem->nodes[FaceNodesMap[1][co]].id;
								if(dy*0.5>d_c1){
									coords[node_change + 1] = coords[node_change + 1]-(d_c1-dy*1.1);
								}else{
									coords[node_change + 1] = coords[node_change + 1]+(d_c1-dy*1.1);
								}
							}
						}else{
							for(int co = 0; co<4;co++){
								node_change = 3*elem->nodes[FaceNodesMap[1][co]].id;
								if(dz*0.5>d_c2){
									coords[node_change + 2] = coords[node_change +2  ] - (d_c2-dz*1.1);
								}else{
									coords[node_change + 2] = coords[node_change +2  ] + (d_c2-dz*1.1);
								}
							}
						}
						//move 2

					}else if(edge_list[2]&&edge_list[6]){
						int node_change = 3*elem->nodes[2].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[2], p0);
						d_c2 = gts_point_distance(point[6], p0);

						if(d_c1>d_c2){
							for(int co = 0; co<4;co++){
								node_change = 3*elem->nodes[FaceNodesMap[3][co]].id;
								if(dx*0.5>d_c1){
									coords[node_change  ] = coords[node_change  ]-(d_c1-dx*1.1);
								}else{
									coords[node_change  ] = coords[node_change  ]+(d_c1-dx*1.1);
								}
							}
						}else{
							for(int co = 0; co<4;co++){
								node_change = 3*elem->nodes[FaceNodesMap[3][co]].id;
								if(dz*0.5>d_c2){
									coords[node_change + 2] = coords[node_change +2  ] - (d_c2-dz*1.1);
								}else{
									coords[node_change + 2] = coords[node_change +2  ] + (d_c2-dz*1.1);
								}
							}
						}
						//move 2

					}else if(edge_list[2]&&edge_list[7]){
						int node_change = 3*elem->nodes[3].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[2], p0);
						d_c2 = gts_point_distance(point[7], p0);

						if(d_c1>d_c2){
							for(int co = 0; co<4;co++){
								node_change = 3*elem->nodes[FaceNodesMap[3][co]].id;
								if(dx*0.5>d_c1){
									coords[node_change  ] = coords[node_change  ]-(d_c1-dx*1.1);
								}else{
									coords[node_change  ] = coords[node_change  ]+(d_c1-dx*1.1);
								}
							}
						}else{
							for(int co = 0; co<4;co++){
								node_change = 3*elem->nodes[FaceNodesMap[3][co]].id;
								if(dz*0.5>d_c2){
									coords[node_change + 2] = coords[node_change +2  ] - (d_c2-dz*1.1);
								}else{
									coords[node_change + 2] = coords[node_change +2  ] + (d_c2-dz*1.1);
								}
							}
						}
						//move 3

					}else if(edge_list[3]&&edge_list[7]){
						int node_change = 3*elem->nodes[3].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[3], p0);
						d_c2 = gts_point_distance(point[7], p0);

						if(d_c1>d_c2){
							for(int co = 0; co<4;co++){
								node_change = 3*elem->nodes[FaceNodesMap[0][co]].id;
								if(dy*0.5>d_c1){
									coords[node_change + 1] = coords[node_change + 1]-(d_c1-dy*1.1);
								}else{
									coords[node_change + 1] = coords[node_change + 1]+(d_c1-dy*1.1);
								}
							}
						}else{
							for(int co = 0; co<4;co++){
								node_change = 3*elem->nodes[FaceNodesMap[0][co]].id;
								if(dz*0.5>d_c2){
									coords[node_change + 2] = coords[node_change +2  ] - (d_c2-dz*1.1);
								}else{
									coords[node_change + 2] = coords[node_change +2  ] + (d_c2-dz*1.1);
								}
							}
						}
						//move 3

					}else if(edge_list[8]&&edge_list[4]){
						int node_change = 3*elem->nodes[4].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[8], p0);
						d_c2 = gts_point_distance(point[4], p0);

						if(d_c1>d_c2){
							for(int co = 0; co<4;co++){
								node_change = 3*elem->nodes[FaceNodesMap[2][co]].id;
								if(dx*0.5>d_c1){
									coords[node_change  ] = coords[node_change  ]-(d_c1-dx*1.1);
								}else{
									coords[node_change  ] = coords[node_change  ]+(d_c1-dx*1.1);
								}
							}
						}else{
							for(int co = 0; co<4;co++){
								node_change = 3*elem->nodes[FaceNodesMap[2][co]].id;
								if(dz*0.5>d_c2){
									coords[node_change + 2] = coords[node_change +2  ] - (d_c2-dz*1.1);
								}else{
									coords[node_change + 2] = coords[node_change +2  ] + (d_c2-dz*1.1);
								}
							}
						}
						//move 4

					}else if(edge_list[11]&&edge_list[4]){
						int node_change = 3*elem->nodes[4].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[11], p0);
						d_c2 = gts_point_distance(point[4], p0);

						if(d_c1>d_c2){
							for(int co = 0; co<4;co++){
								node_change = 3*elem->nodes[FaceNodesMap[0][co]].id;
								if(dy*0.5>d_c1){
									coords[node_change + 1] = coords[node_change + 1]-(d_c1-dy*1.1);
								}else{
									coords[node_change + 1] = coords[node_change + 1]+(d_c1-dy*1.1);
								}
							}
						}else{
							for(int co = 0; co<4;co++){
								node_change = 3*elem->nodes[FaceNodesMap[0][co]].id;
								if(dz*0.5>d_c2){
									coords[node_change + 2] = coords[node_change +2  ] - (d_c2-dz*1.1);
								}else{
									coords[node_change + 2] = coords[node_change +2  ] + (d_c2-dz*1.1);
								}
							}
						}
						//move 4

					}else if(edge_list[8]&&edge_list[5]){
						int node_change = 3*elem->nodes[5].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[8], p0);
						d_c2 = gts_point_distance(point[5], p0);

						if(d_c1>d_c2){
							for(int co = 0; co<4;co++){
								node_change = 3*elem->nodes[FaceNodesMap[2][co]].id;
								if(dx*0.5>d_c1){
									coords[node_change  ] = coords[node_change  ]-(d_c1-dx*1.1);
								}else{
									coords[node_change  ] = coords[node_change  ]+(d_c1-dx*1.1);
								}
							}
						}else{
							for(int co = 0; co<4;co++){
								node_change = 3*elem->nodes[FaceNodesMap[2][co]].id;
								if(dz*0.5>d_c2){
									coords[node_change + 2] = coords[node_change +2  ] - (d_c2-dz*1.1);
								}else{
									coords[node_change + 2] = coords[node_change +2  ] + (d_c2-dz*1.1);
								}
							}
						}
						//move 5

					}else if(edge_list[9]&&edge_list[5]){
						int node_change = 3*elem->nodes[5].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[9], p0);
						d_c2 = gts_point_distance(point[5], p0);

						if(d_c1>d_c2){
							for(int co = 0; co<4;co++){
								node_change = 3*elem->nodes[FaceNodesMap[1][co]].id;
								if(dy*0.5>d_c1){
									coords[node_change + 1] = coords[node_change + 1]-(d_c1-dy*1.1);
								}else{
									coords[node_change + 1] = coords[node_change + 1]+(d_c1-dy*1.1);
								}
							}
						}else{
							for(int co = 0; co<4;co++){
								node_change = 3*elem->nodes[FaceNodesMap[1][co]].id;
								if(dz*0.5>d_c2){
									coords[node_change + 2] = coords[node_change +2  ] - (d_c2-dz*1.1);
								}else{
									coords[node_change + 2] = coords[node_change +2  ] + (d_c2-dz*1.1);
								}
							}
						}
						//move 5

					}else if(edge_list[9]&&edge_list[6]){
						int node_change = 3*elem->nodes[6].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[9], p0);
						d_c2 = gts_point_distance(point[6], p0);

						if(d_c1>d_c2){
							for(int co = 0; co<4;co++){
								node_change = 3*elem->nodes[FaceNodesMap[1][co]].id;
								if(dy*0.5>d_c1){
									coords[node_change + 1] = coords[node_change + 1]-(d_c1-dy*1.1);
								}else{
									coords[node_change + 1] = coords[node_change + 1]+(d_c1-dy*1.1);
								}
							}
						}else{
							for(int co = 0; co<4;co++){
								node_change = 3*elem->nodes[FaceNodesMap[1][co]].id;
								if(dz*0.5>d_c2){
									coords[node_change + 2] = coords[node_change +2  ] - (d_c2-dz*1.1);
								}else{
									coords[node_change + 2] = coords[node_change +2  ] + (d_c2-dz*1.1);
								}
							}
						}
						//move 6

					}else if(edge_list[10]&&edge_list[6]){
						int node_change = 3*elem->nodes[6].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[10], p0);
						d_c2 = gts_point_distance(point[6], p0);

						if(d_c1>d_c2){
							for(int co = 0; co<4;co++){
								node_change = 3*elem->nodes[FaceNodesMap[3][co]].id;
								if(dx*0.5>d_c1){
									coords[node_change  ] = coords[node_change  ]-(d_c1-dx*1.1);
								}else{
									coords[node_change  ] = coords[node_change  ]+(d_c1-dx*1.1);
								}
							}
						}else{
							for(int co = 0; co<4;co++){
								node_change = 3*elem->nodes[FaceNodesMap[3][co]].id;
								if(dz*0.5>d_c2){
									coords[node_change + 2] = coords[node_change +2  ] - (d_c2-dz*1.1);
								}else{
									coords[node_change + 2] = coords[node_change +2  ] + (d_c2-dz*1.1);
								}
							}
						}
						//move 6

					}else if(edge_list[10]&&edge_list[7]){
						int node_change = 3*elem->nodes[7].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[10], p0);
						d_c2 = gts_point_distance(point[7], p0);

						if(d_c1>d_c2){
							for(int co = 0; co<4;co++){
								node_change = 3*elem->nodes[FaceNodesMap[3][co]].id;
								if(dx*0.5>d_c1){
									coords[node_change  ] = coords[node_change  ]-(d_c1-dx*1.1);
								}else{
									coords[node_change  ] = coords[node_change  ]+(d_c1-dx*1.1);
								}
							}
						}else{
							for(int co = 0; co<4;co++){
								node_change = 3*elem->nodes[FaceNodesMap[3][co]].id;
								if(dz*0.5>d_c2){
									coords[node_change + 2] = coords[node_change +2  ] - (d_c2-dz*1.1);
								}else{
									coords[node_change + 2] = coords[node_change +2  ] + (d_c2-dz*1.1);
								}
							}
						}
						//move 7

					}else if(edge_list[11]&&edge_list[7]){
						int node_change = 3*elem->nodes[7].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[11], p0);
						d_c2 = gts_point_distance(point[7], p0);

						if(d_c1>d_c2){
							for(int co = 0; co<4;co++){
								node_change = 3*elem->nodes[FaceNodesMap[0][co]].id;
								if(dy*0.5>d_c1){
									coords[node_change + 1] = coords[node_change + 1]-(d_c1-dy*1.1);
								}else{
									coords[node_change + 1] = coords[node_change + 1]+(d_c1-dy*1.1);
								}
							}
						}else{
							for(int co = 0; co<4;co++){
								node_change = 3*elem->nodes[FaceNodesMap[0][co]].id;
								if(dz*0.5>d_c2){
									coords[node_change + 2] = coords[node_change +2  ] - (d_c2-dz*1.1);
								}else{
									coords[node_change + 2] = coords[node_change +2  ] + (d_c2-dz*1.1);
								}
							}
						}
						//move 7
					}


					CheckTemplate(mesh, coords, element_ids_local,false);
					elem = (octant_t*) sc_array_index(&mesh->elements, elements_ids[iel]);

					fprintf(fdbg,"Pad local:    ");
					for(int co = 0; co < element_ids_local.size(); co++){
						octant_t *h = (octant_t*) sc_array_index(&mesh->elements, element_ids_local[co]);
						var_aux_1[co]=h->pad;
						fprintf(fdbg,"%d ",h->pad);
					}
					fprintf(fdbg,"\n");

					if(elem->pad!=4 && elem->pad!=-10){
						flag2=false;
						for(int co = 0; co < element_ids_local.size(); co++){
							if(element_ids_local[co]!=elements_ids[iel]){
								if(var_aux[co]==var_aux_1[co]){
									flag2=false;
								}else{
									if(var_aux_1[co]==4){
										flag2=true;
										break;
									}else if(var_aux_1[co]==-10){
										flag2=true;
										break;
									}
								}
							}
						}
					}
				}

				if(!flag2){
					flag1=false;
					flag2=false;
					printf("To livreeeee rand tipo 4 movi minha fesse todinha :) El:%d, n_iter:%d\n",elements_ids[iel],n_iter);
					break;
				}

				if( flag2 ){
					for(int n_nodes = 0; n_nodes<8;n_nodes++){
						int node = 3*elem->nodes[n_nodes].id;
						coords[node]   = coords_orgi[3*n_nodes];
						coords[node+1] = coords_orgi[3*n_nodes+1];
						coords[node+2] = coords_orgi[3*n_nodes+2];
					}
				}


#endif
				n_iter++;
				v_1= 1;
			}

			el_4++;

		}else if(elem->pad==-10){

			v_1 = 0.75;
			n_iter = 0;
			bool flag1 = false;
			bool flag2 = true;

			fprintf(fdbg,"El: %d not handle\n",elements_ids[iel]);
			//printf("El:%d\n",elements_ids[iel]);
			element_ids_local.clear();
			SideEL(mesh, elem->x, elem->y, elem->z ,element_ids_local);

			fprintf(fdbg,"El: %d Case type 10\n",elements_ids[iel]);
			fprintf(fdbg,"El:");
			for(int co = 0; co<element_ids_local.size();co++){
				fprintf(fdbg,"%d ",element_ids_local[co]);
				if(element_ids_local[co]==elements_ids[iel]){ref=co;}
			}
			fprintf(fdbg,"\n");

			fprintf(fdbg,"REF:%d \n",ref);

			fprintf(fdbg,"Pad Original: ");
			for(int c = 0; c < element_ids_local.size(); c++){
				octant_t *h = (octant_t*) sc_array_index(&mesh->elements, element_ids_local[c]);
				var_aux[c]=h->pad;
				fprintf(fdbg,"%d ",var_aux[c]);
			}
			fprintf(fdbg,"\n");

			for(int n_nodes = 0; n_nodes<8;n_nodes++){
				octant_t *elem = (octant_t*) sc_array_index(&mesh->elements, elements_ids[iel]);
				int node = 3*elem->nodes[n_nodes].id;
				coords_orgi[3*n_nodes  ] = coords[node  ];
				coords_orgi[3*n_nodes+1] = coords[node+1];
				coords_orgi[3*n_nodes+2] = coords[node+2];
			}

			while(n_iter<iter_max && flag1){
				template2Rand(mesh, coords, elem, v_1);
				CheckTemplate(mesh, coords, element_ids_local,false);
				elem = (octant_t*) sc_array_index(&mesh->elements, elements_ids[iel]);

				fprintf(fdbg,"Pad local:    ");
				for(int co = 0; co < element_ids_local.size(); co++){
					octant_t *h = (octant_t*) sc_array_index(&mesh->elements, element_ids_local[co]);
					var_aux_1[co]=h->pad;
					fprintf(fdbg,"%d ",h->pad);
				}
				fprintf(fdbg,"\n");

				if(elem->pad!=4 && elem->pad!=-10){
					flag2=false;
					for(int co = 0; co < element_ids_local.size(); co++){
						if(element_ids_local[co]!=elements_ids[iel]){
							if(var_aux[co]==var_aux_1[co]){
								flag2=false;
							}else{
								if(var_aux_1[co]==4){
									flag2=true;
									break;
								}else if(var_aux_1[co]==-10){
									flag2=true;
									break;
								}
							}
						}
					}
				}

				if(!flag2){
					flag1=false;
					flag2=false;
					printf("To livreeeee rand tipo 10 manuuu El:%d, n_iter:%d\n",elements_ids[iel],n_iter);
					break;
				}

				if( flag2 ){
					for(int n_nodes = 0; n_nodes<8;n_nodes++){
						int node = 3*elem->nodes[n_nodes].id;
						coords[node]   = coords_orgi[3*n_nodes];
						coords[node+1] = coords_orgi[3*n_nodes+1];
						coords[node+2] = coords_orgi[3*n_nodes+2];
					}
				}
				n_iter++;
			}
			nao_sei++;

		}

	}

	fclose(fdbg);
#endif
}

void ApplyTemplate(hexa_tree_t* mesh, std::vector<double>& coords, std::vector<int>& elements_ids) {
	bool clamped = true;
	bool face_intecepted[6];
	int conn_p[8];
	int original_conn[8];
	FILE * fdbg;

	int Edge2GNode[12][2]={0};
	int Edge2GNode_s[12][2]={0};
	int Edge2GNode_v[4][2]={0};

	GtsSegment * segments[12]={0};
	GtsSegment * segments_s[12]={0};
	GtsSegment * segments_v[4]={0};

	GtsPoint * point[12]={NULL};
	GtsPoint * point_s[12]={NULL};
	GtsPoint * point_v[4]={NULL};
	bool edge_list[12]={false};
	bool edge_list_s[12]={false};
	bool edge_list_v[4]={false};
	int ed_cont = 0;

	fdbg = fopen("intercepted_faces.dbg", "w");

	sc_hash_array_t* hash_nodes = sc_hash_array_new(sizeof (node_in_edge_t), edge_hash_fn, edge_equal_fn, &clamped);

	for (int iel = 0; iel < elements_ids.size(); ++iel) {

		octant_t *elem = (octant_t*) sc_array_index(&mesh->elements, elements_ids[iel]);

		for (int i = 0; i < 8; i++) original_conn[i] = elem->nodes[i].id;

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
			edge_list[edge] = false;
			if (list == NULL) continue;
			while (list) {
				GtsBBox *b = GTS_BBOX(list->data);
				point[edge] = SegmentTriangleIntersection(segments[edge], GTS_TRIANGLE(b->bounded));
				if (point[edge]) {
					edge_list[edge] = true;
					ed_cont++;
					break;
				}
				list = list->next;
			}
		}

		//check the diagonals in the surface and find the intersections
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
			edge_list_s[edge]=false;
			if (list == NULL) continue;
			while (list) {
				GtsBBox *b = GTS_BBOX(list->data);
				point_s[edge] = SegmentTriangleIntersection(segments_s[edge], GTS_TRIANGLE(b->bounded));
				if (point_s[edge]) {
					edge_list_s[edge]=true;
					break;
				}
				list = list->next;
			}
		}

		//check the diagonals in the volume and find the intersections
		for (int edge = 0; edge < 4; ++edge) {
			point_v[edge] = NULL;
			int node1 = elem->nodes[EdgeVerticesMap_vol_diagonal[edge][0]].id;
			int node2 = elem->nodes[EdgeVerticesMap_vol_diagonal[edge][1]].id;

			Edge2GNode_v[edge][0] = node1 <= node2 ? node1 : node2;
			Edge2GNode_v[edge][1] = node1 >= node2 ? node1 : node2;

			GtsVertex *v1 = gts_vertex_new(gts_vertex_class(), coords[node1 * 3], coords[node1 * 3 + 1], coords[node1 * 3 + 2]);
			GtsVertex *v2 = gts_vertex_new(gts_vertex_class(), coords[node2 * 3], coords[node2 * 3 + 1], coords[node2 * 3 + 2]);

			segments_v[edge] = gts_segment_new(gts_segment_class(), v1, v2);
			GtsBBox *sb = gts_bbox_segment(gts_bbox_class(), segments_v[edge]);
			GSList* list = gts_bb_tree_overlap(mesh->gdata.bbt, sb);
			edge_list_v[edge]=false;
			if (list == NULL) continue;
			while (list) {
				GtsBBox *b = GTS_BBOX(list->data);
				point_v[edge] = SegmentTriangleIntersection(segments_v[edge], GTS_TRIANGLE(b->bounded));
				if (point_v[edge]) {
					edge_list_v[edge]=true;
					break;
				}
				list = list->next;
			}
		}

		//check parallel faces
		for (int face = 0; face < 6; ++face) {
			face_intecepted[face] = false;
			for (int fe = 0; fe < 4; ++fe) {
				int edge = FaceEdgesMap[face][fe];
				if (point[edge] != NULL) {
					face_intecepted[face] = true;
					break;
				}
			}
			if (face_intecepted[face]) continue;
		}

		fprintf(fdbg, "Faces:  Elem. %d: %d %d %d %d %d %d\n", elements_ids[iel], face_intecepted[0],
				face_intecepted[1], face_intecepted[2],
				face_intecepted[3], face_intecepted[4], face_intecepted[5]);
		fprintf(fdbg, "Edges:  Elem. %d: %d %d %d %d %d %d %d %d %d %d %d %d %d\n", elements_ids[iel], edge_list[0],
				edge_list[1], edge_list[2], edge_list[3],
				edge_list[4], edge_list[5], edge_list[6],
				edge_list[7], edge_list[8], edge_list[9],
				edge_list[10], edge_list[11], ed_cont);

		int n_parallel_faces = (face_intecepted[0] && face_intecepted[1]) +
				(face_intecepted[2] && face_intecepted[3]) +
				(face_intecepted[4] && face_intecepted[5]);



		if(elem->pad == 22){
			//(edge_list[0]) && (!edge_list[1]) && (edge_list[2]) && (!edge_list[3]) &&
			//(!edge_list[4]) && (!edge_list[5]) && (!edge_list[6]) && (!edge_list[7]) &&
			//(edge_list[8]) && (!edge_list[9]) && (edge_list[10]) && (!edge_list[11])


		}else if(elem->pad == 23){
			//(!edge_list[0]) && (edge_list[1]) && (!edge_list[2]) && (edge_list[3]) &&
			//(!edge_list[4]) && (!edge_list[5]) && (!edge_list[6]) && (!edge_list[7]) &&
			//(!edge_list[8]) && (edge_list[9]) && (!edge_list[10]) && (edge_list[11])


		}else if(elem->pad == 24){
			//(!edge_list[0]) && (!edge_list[1]) && (!edge_list[2]) && (!edge_list[3]) &&
			//(edge_list[4]) && (edge_list[5]) && (edge_list[6]) && (edge_list[7]) &&
			//(!edge_list[8]) && (!edge_list[9]) && (!edge_list[10]) && (!edge_list[11])

			int cut_edge[4];

			cut_edge[0] = 4;
			cut_edge[1] = 5;
			cut_edge[2] = 6;
			cut_edge[3] = 7;

			GtsPoint *p0 = point[0];
			GtsPoint *p1 = point[0];
			GtsPoint *p2 = point[0];
			GtsPoint *p3 = point[0];
			GtsPoint *p4 = point[0];
			GtsPoint *p5 = point[0];
			GtsPoint *p6 = point[0];
			GtsPoint *p7 = point[0];

			if(p0==NULL) exit(1);
			if(p1==NULL) exit(1);
			if(p2==NULL) exit(1);
			if(p3==NULL) exit(1);
			if(p4==NULL) exit(1);
			if(p5==NULL) exit(1);
			if(p6==NULL) exit(1);
			if(p7==NULL) exit(1);

			conn_p[0] = AddPointOnEdge(Edge2GNode[cut_edge[0]], hash_nodes, mesh->local_n_nodes, p0, coords);
			conn_p[1] = AddPointOnEdge(Edge2GNode[cut_edge[0]], hash_nodes, mesh->local_n_nodes, p1, coords);
			conn_p[2] = AddPointOnEdge(Edge2GNode[cut_edge[1]], hash_nodes, mesh->local_n_nodes, p2, coords);
			conn_p[3] = AddPointOnEdge(Edge2GNode[cut_edge[1]], hash_nodes, mesh->local_n_nodes, p3, coords);
			conn_p[4] = AddPointOnEdge(Edge2GNode[cut_edge[2]], hash_nodes, mesh->local_n_nodes, p4, coords);
			conn_p[5] = AddPointOnEdge(Edge2GNode[cut_edge[2]], hash_nodes, mesh->local_n_nodes, p5, coords);
			conn_p[6] = AddPointOnEdge(Edge2GNode[cut_edge[3]], hash_nodes, mesh->local_n_nodes, p6, coords);
			conn_p[7] = AddPointOnEdge(Edge2GNode[cut_edge[3]], hash_nodes, mesh->local_n_nodes, p7, coords);

			octant_t *elem1 = (octant_t*) sc_array_index(&mesh->elements, elements_ids[iel]);

			elem1->nodes[0].id = original_conn[0];
			elem1->nodes[1].id = original_conn[1];
			elem1->nodes[2].id = original_conn[2];
			elem1->nodes[3].id = original_conn[3];

			elem1->nodes[4].id = conn_p[0];
			elem1->nodes[5].id = conn_p[1];
			elem1->nodes[6].id = conn_p[2];
			elem1->nodes[7].id = conn_p[3];

			elem1->level = elem->level;
			elem1->ref = 1;

			octant_t* elem2 = (octant_t*) sc_array_push(&mesh->elements);

			elem2->nodes[0].id = conn_p[0];
			elem2->nodes[1].id = conn_p[0];
			elem2->nodes[2].id = conn_p[0];
			elem2->nodes[3].id = conn_p[0];

			elem2->nodes[4].id = conn_p[0];
			elem2->nodes[5].id = conn_p[1];
			elem2->nodes[6].id = conn_p[2];
			elem2->nodes[7].id = conn_p[3];

			elem2->pad = elem->pad;
			elem2->level = elem->level;
			elem2->ref = 1;

			octant_t* elem3 = (octant_t*) sc_array_push(&mesh->elements);

			elem3->nodes[0].id = original_conn[0];
			elem3->nodes[1].id = original_conn[1];
			elem3->nodes[2].id = original_conn[2];
			elem3->nodes[3].id = original_conn[3];

			elem3->nodes[4].id = conn_p[0];
			elem3->nodes[5].id = conn_p[1];
			elem3->nodes[6].id = conn_p[2];
			elem3->nodes[7].id = conn_p[3];

			elem3->pad = elem->pad;
			elem3->level = elem->level;
			elem3->ref = 1;


		}


#if 0
		// verificação dos elementos
		if(elem->pad==1){
			if (n_parallel_faces == 2 && ed_cont == 4) {
				// Check template 1.
				int cut_edge[4];
				int connec_order_in[8];
				int connec_order_out1[8];
				int connec_order_out2[8];

				if (    (!edge_list[0]) && (edge_list[1]) && (!edge_list[2]) && (edge_list[3]) &&
						(!edge_list[4]) && (!edge_list[5]) && (!edge_list[6]) && (!edge_list[7]) &&
						(!edge_list[8]) && (edge_list[9]) && (!edge_list[10]) && (edge_list[11]) ) {

					cut_edge[0] = 3;
					cut_edge[1] = 1;
					cut_edge[2] = 11;
					cut_edge[3] = 9;

					connec_order_in[0] = 0;
					connec_order_in[1] = 1;
					connec_order_in[2] = 4;
					connec_order_in[3] = 5;
					connec_order_in[4] = 2;
					connec_order_in[5] = 3;
					connec_order_in[6] = 6;
					connec_order_in[7] = 7;

					connec_order_out1[0] = 0;
					connec_order_out1[1] = 1;
					connec_order_out1[2] = 4;
					connec_order_out1[3] = 5;
					connec_order_out1[4] = 3;
					connec_order_out1[5] = 2;
					connec_order_out1[6] = 7;
					connec_order_out1[7] = 6;

					connec_order_out2[0] = 0;
					connec_order_out2[1] = 1;
					connec_order_out2[2] = 4;
					connec_order_out2[3] = 5;
					connec_order_out2[4] = 2;
					connec_order_out2[5] = 3;
					connec_order_out2[6] = 6;
					connec_order_out2[7] = 7;

					elem->pad = 1;

				} else if (     (edge_list[0]) && (!edge_list[1]) && (edge_list[2]) && (!edge_list[3]) &&
						(!edge_list[4]) && (!edge_list[5]) && (!edge_list[6]) && (!edge_list[7]) &&
						(edge_list[8]) && (!edge_list[9]) && (edge_list[10]) && (!edge_list[11]) ) {

					cut_edge[0] = 0;
					cut_edge[1] = 2;
					cut_edge[2] = 8;
					cut_edge[3] = 10;

					connec_order_in[0] = 0;
					connec_order_in[1] = 3;
					connec_order_in[2] = 4;
					connec_order_in[3] = 7;
					connec_order_in[4] = 1;
					connec_order_in[5] = 2;
					connec_order_in[6] = 5;
					connec_order_in[7] = 6;

					connec_order_out1[0] = 0;
					connec_order_out1[1] = 3;
					connec_order_out1[2] = 4;
					connec_order_out1[3] = 7;
					connec_order_out1[4] = 1;
					connec_order_out1[5] = 2;
					connec_order_out1[6] = 5;
					connec_order_out1[7] = 6;

					connec_order_out2[0] = 0;
					connec_order_out2[1] = 3;
					connec_order_out2[2] = 4;
					connec_order_out2[3] = 7;
					connec_order_out2[4] = 1;
					connec_order_out2[5] = 2;
					connec_order_out2[6] = 5;
					connec_order_out2[7] = 6;

					elem->pad = 1;

				} else if (     (!edge_list[0]) && (!edge_list[1]) && (!edge_list[2]) && (!edge_list[3]) &&
						(edge_list[4]) && (edge_list[5]) && (edge_list[6]) && (edge_list[7]) &&
						(!edge_list[8]) && (!edge_list[9]) && (!edge_list[10]) && (!edge_list[11]) ) {

					cut_edge[0] = 4;
					cut_edge[1] = 5;
					cut_edge[2] = 6;
					cut_edge[3] = 7;

					connec_order_in[0] = 0;
					connec_order_in[1] = 1;
					connec_order_in[2] = 2;
					connec_order_in[3] = 3;
					connec_order_in[4] = 4;
					connec_order_in[5] = 5;
					connec_order_in[6] = 6;
					connec_order_in[7] = 7;

					connec_order_out1[0] = 0;
					connec_order_out1[1] = 1;
					connec_order_out1[2] = 2;
					connec_order_out1[3] = 3;
					connec_order_out1[4] = 4;
					connec_order_out1[5] = 5;
					connec_order_out1[6] = 6;
					connec_order_out1[7] = 7;

					connec_order_out2[0] = 0;
					connec_order_out2[1] = 1;
					connec_order_out2[2] = 2;
					connec_order_out2[3] = 3;
					connec_order_out2[4] = 4;
					connec_order_out2[5] = 5;
					connec_order_out2[6] = 6;
					connec_order_out2[7] = 7;

					elem->pad = 1;

				} else {
					elem->pad = -1;
					printf("Warning: Bad Element: %d, case 1\n",elements_ids[iel]);
					printf("Edges:  Elem. %d: %d %d %d %d %d %d %d %d %d %d %d %d %d\n", elements_ids[iel], edge_list[0],
							edge_list[1], edge_list[2], edge_list[3],
							edge_list[4], edge_list[5], edge_list[6],
							edge_list[7], edge_list[8], edge_list[9],
							edge_list[10], edge_list[11], ed_cont);

					continue;
				}

				GtsPoint *p0 = point[cut_edge[0]];
				GtsPoint *p1 = point[cut_edge[1]];
				GtsPoint *p2 = point[cut_edge[2]];
				GtsPoint *p3 = point[cut_edge[3]];

				//g_assert(p0 != NULL);
				//g_assert(p1 != NULL);
				//g_assert(p2 != NULL);
				//g_assert(p3 != NULL);

				if(p0==NULL) exit(1);
				if(p1==NULL) exit(1);
				if(p2==NULL) exit(1);
				if(p3==NULL) exit(1);


				conn_p[0] = AddPointOnEdge(Edge2GNode[cut_edge[0]], hash_nodes, mesh->local_n_nodes, p0, coords);
				conn_p[1] = AddPointOnEdge(Edge2GNode[cut_edge[1]], hash_nodes, mesh->local_n_nodes, p1, coords);
				conn_p[2] = AddPointOnEdge(Edge2GNode[cut_edge[2]], hash_nodes, mesh->local_n_nodes, p2, coords);
				conn_p[3] = AddPointOnEdge(Edge2GNode[cut_edge[3]], hash_nodes, mesh->local_n_nodes, p3, coords);

				octant_t *elem1 = (octant_t*) sc_array_index(&mesh->elements, elements_ids[iel]);

				elem1->nodes[connec_order_out1[0]].id = original_conn[connec_order_in[0]];
				elem1->nodes[connec_order_out1[1]].id = original_conn[connec_order_in[1]];
				elem1->nodes[connec_order_out1[2]].id = original_conn[connec_order_in[2]];
				elem1->nodes[connec_order_out1[3]].id = original_conn[connec_order_in[3]];

				elem1->nodes[connec_order_out1[4]].id = conn_p[0];
				elem1->nodes[connec_order_out1[5]].id = conn_p[1];
				elem1->nodes[connec_order_out1[6]].id = conn_p[2];
				elem1->nodes[connec_order_out1[7]].id = conn_p[3];

				octant_t* elem2 = (octant_t*) sc_array_push(&mesh->elements);

				elem2->nodes[connec_order_out2[0]].id = conn_p[0];
				elem2->nodes[connec_order_out2[1]].id = conn_p[1];
				elem2->nodes[connec_order_out2[2]].id = conn_p[2];
				elem2->nodes[connec_order_out2[3]].id = conn_p[3];

				elem2->nodes[connec_order_out2[4]].id = original_conn[connec_order_in[4]];
				elem2->nodes[connec_order_out2[5]].id = original_conn[connec_order_in[5]];
				elem2->nodes[connec_order_out2[6]].id = original_conn[connec_order_in[6]];
				elem2->nodes[connec_order_out2[7]].id = original_conn[connec_order_in[7]];

				elem2->pad = elem->pad;
				elem2->level = elem->level;


			}
		}else if(elem->pad==2){

			if (n_parallel_faces == 1 && ed_cont==4) {
				// Check template 2.

				int conn_t2[6];
				int cut_edge[4];
				int cut_edge_s[2];
				int connec_order_in1[8];
				int connec_order_out1[8];
				int connec_order_in2[8];
				int connec_order_out2[8];
				int connec_order_in3[8];
				int connec_order_out3[8];

				//Edge 0
				if ( (!edge_list[0]) && (edge_list[1]) && (!edge_list[2]) && (edge_list[3]) &&
						(edge_list[4]) && (edge_list[5]) && (!edge_list[6]) && (!edge_list[7]) &&
						(!edge_list[8]) && (!edge_list[9]) && (!edge_list[10]) && (!edge_list[11]) &&
						(edge_list_s[4]) && (edge_list_s[6]) ) {
					//printf("Entrou na 0!\n");

					cut_edge[0] = 3;
					cut_edge[1] = 1;
					cut_edge[2] = 4;
					cut_edge[3] = 5;

					cut_edge_s[0] = 4;
					cut_edge_s[1] = 6;

					connec_order_in1[0] = 2;
					connec_order_in1[1] = 3;
					connec_order_in1[2] = 6;
					connec_order_in1[3] = 7;
					connec_order_in1[4] = 0;
					connec_order_in1[5] = 1;
					connec_order_in1[6] = 4;
					connec_order_in1[7] = 5;

					connec_order_out1[0] = 2;
					connec_order_out1[1] = 3;
					connec_order_out1[2] = 6;
					connec_order_out1[3] = 7;
					connec_order_out1[4] = 0;
					connec_order_out1[5] = 1;
					connec_order_out1[6] = 4;
					connec_order_out1[7] = 5;

					connec_order_in2[0] = 0;
					connec_order_in2[1] = 1;
					connec_order_in2[2] = 2;
					connec_order_in2[3] = 3;
					connec_order_in2[4] = 4;
					connec_order_in2[5] = 5;
					connec_order_in2[6] = 6;
					connec_order_in2[7] = 7;

					connec_order_out2[0] = 2;
					connec_order_out2[1] = 3;
					connec_order_out2[2] = 5;
					connec_order_out2[3] = 4;
					connec_order_out2[4] = 4;
					connec_order_out2[5] = 5;
					connec_order_out2[6] = 6;
					connec_order_out2[7] = 7;

					connec_order_in3[0] = 0;
					connec_order_in3[1] = 1;
					connec_order_in3[2] = 2;
					connec_order_in3[3] = 3;
					connec_order_in3[4] = 4;
					connec_order_in3[5] = 5;
					connec_order_in3[6] = 6;
					connec_order_in3[7] = 7;

					connec_order_out3[0] = 0;
					connec_order_out3[1] = 1;
					connec_order_out3[2] = 1;
					connec_order_out3[3] = 0;
					connec_order_out3[4] = 2;
					connec_order_out3[5] = 3;
					connec_order_out3[6] = 5;
					connec_order_out3[7] = 4;

					elem->pad = 2;

				}//Edge 1
				else if ((edge_list[0]) && (!edge_list[1]) && (edge_list[2]) && (!edge_list[3]) &&
						(!edge_list[4]) && (edge_list[5]) && (edge_list[6]) && (!edge_list[7]) &&
						(!edge_list[8]) && (!edge_list[9]) && (!edge_list[10]) && (!edge_list[11])&&
						(edge_list_s[1]) && (edge_list_s[3]) ) {
					//printf("Entrou na 1!\n");

					cut_edge[0] = 0;
					cut_edge[1] = 2;
					cut_edge[2] = 5;
					cut_edge[3] = 6;

					cut_edge_s[0] = 1;
					cut_edge_s[1] = 3;

					connec_order_in1[0] = 0;
					connec_order_in1[1] = 3;
					connec_order_in1[2] = 4;
					connec_order_in1[3] = 7;
					connec_order_in1[4] = 1;
					connec_order_in1[5] = 2;
					connec_order_in1[6] = 5;
					connec_order_in1[7] = 6;

					connec_order_out1[0] = 0;
					connec_order_out1[1] = 3;
					connec_order_out1[2] = 4;
					connec_order_out1[3] = 7;
					connec_order_out1[4] = 0;
					connec_order_out1[5] = 1;
					connec_order_out1[6] = 4;
					connec_order_out1[7] = 5;

					connec_order_in2[0] = 0;
					connec_order_in2[1] = 1;
					connec_order_in2[2] = 2;
					connec_order_in2[3] = 3;
					connec_order_in2[4] = 4;
					connec_order_in2[5] = 5;
					connec_order_in2[6] = 6;
					connec_order_in2[7] = 7;

					connec_order_out2[0] = 4;
					connec_order_out2[1] = 2;
					connec_order_out2[2] = 3;
					connec_order_out2[3] = 5;
					connec_order_out2[4] = 4;
					connec_order_out2[5] = 5;
					connec_order_out2[6] = 6;
					connec_order_out2[7] = 7;

					connec_order_in3[0] = 1;
					connec_order_in3[1] = 2;
					connec_order_in3[2] = 0;
					connec_order_in3[3] = 3;
					connec_order_in3[4] = 4;
					connec_order_in3[5] = 5;
					connec_order_in3[6] = 6;
					connec_order_in3[7] = 7;

					connec_order_out3[0] = 1;
					connec_order_out3[1] = 2;
					connec_order_out3[2] = 0;
					connec_order_out3[3] = 1;
					connec_order_out3[4] = 4;
					connec_order_out3[5] = 2;
					connec_order_out3[6] = 3;
					connec_order_out3[7] = 5;

					elem->pad = 2;

				}//Edge 2
				else if ((!edge_list[0]) && (edge_list[1]) && (!edge_list[2]) && (edge_list[3]) &&
						(!edge_list[4]) && (!edge_list[5]) && (edge_list[6]) && (edge_list[7]) &&
						(!edge_list[8]) && (!edge_list[9]) && (!edge_list[10]) && (!edge_list[11])&&
						(edge_list_s[5]) && (edge_list_s[7]) ) {
					//printf("Entrou na 2!\n");

					cut_edge[0] = 3;
					cut_edge[1] = 1;
					cut_edge[2] = 7;
					cut_edge[3] = 6;

					cut_edge_s[0] = 5;
					cut_edge_s[1] = 7;

					connec_order_in1[0] = 0;
					connec_order_in1[1] = 1;
					connec_order_in1[2] = 4;
					connec_order_in1[3] = 5;
					connec_order_in1[4] = 2;
					connec_order_in1[5] = 3;
					connec_order_in1[6] = 6;
					connec_order_in1[7] = 7;

					connec_order_out1[0] = 0;
					connec_order_out1[1] = 1;
					connec_order_out1[2] = 4;
					connec_order_out1[3] = 5;
					connec_order_out1[4] = 1;
					connec_order_out1[5] = 0;
					connec_order_out1[6] = 5;
					connec_order_out1[7] = 4;

					connec_order_in2[0] = 0;
					connec_order_in2[1] = 1;
					connec_order_in2[2] = 2;
					connec_order_in2[3] = 3;
					connec_order_in2[4] = 4;
					connec_order_in2[5] = 5;
					connec_order_in2[6] = 6;
					connec_order_in2[7] = 7;

					connec_order_out2[0] = 4;
					connec_order_out2[1] = 5;
					connec_order_out2[2] = 3;
					connec_order_out2[3] = 2;
					connec_order_out2[4] = 4;
					connec_order_out2[5] = 5;
					connec_order_out2[6] = 6;
					connec_order_out2[7] = 7;

					connec_order_in3[0] = 2;
					connec_order_in3[1] = 3;
					connec_order_in3[2] = 0;
					connec_order_in3[3] = 1;
					connec_order_in3[4] = 4;
					connec_order_in3[5] = 5;
					connec_order_in3[6] = 6;
					connec_order_in3[7] = 7;

					connec_order_out3[0] = 2;
					connec_order_out3[1] = 3;
					connec_order_out3[2] = 0;
					connec_order_out3[3] = 1;
					connec_order_out3[4] = 4;
					connec_order_out3[5] = 5;
					connec_order_out3[6] = 3;
					connec_order_out3[7] = 2;

					elem->pad = 2;

				}//Edge 3
				else if ((edge_list[0]) && (!edge_list[1]) && (edge_list[2]) && (!edge_list[3]) &&
						(edge_list[4]) && (!edge_list[5]) && (!edge_list[6]) && (edge_list[7]) &&
						(!edge_list[8]) && (!edge_list[9]) && (!edge_list[10]) && (!edge_list[11])&&
						(edge_list_s[0]) && (edge_list_s[2])) {
					//printf("Entrou na 3!\n");

					cut_edge[0] = 0;
					cut_edge[1] = 2;
					cut_edge[2] = 4;
					cut_edge[3] = 7;

					cut_edge_s[0] = 0;
					cut_edge_s[1] = 2;

					connec_order_in1[0] = 1;
					connec_order_in1[1] = 2;
					connec_order_in1[2] = 5;
					connec_order_in1[3] = 6;
					connec_order_in1[4] = 0;
					connec_order_in1[5] = 3;
					connec_order_in1[6] = 4;
					connec_order_in1[7] = 7;

					connec_order_out1[0] = 1;
					connec_order_out1[1] = 2;
					connec_order_out1[2] = 5;
					connec_order_out1[3] = 6;
					connec_order_out1[4] = 0;
					connec_order_out1[5] = 1;
					connec_order_out1[6] = 4;
					connec_order_out1[7] = 5;

					connec_order_in2[0] = 0;
					connec_order_in2[1] = 1;
					connec_order_in2[2] = 2;
					connec_order_in2[3] = 3;
					connec_order_in2[4] = 4;
					connec_order_in2[5] = 5;
					connec_order_in2[6] = 6;
					connec_order_in2[7] = 7;

					connec_order_out2[0] = 2;
					connec_order_out2[1] = 4;
					connec_order_out2[2] = 5;
					connec_order_out2[3] = 3;
					connec_order_out2[4] = 4;
					connec_order_out2[5] = 5;
					connec_order_out2[6] = 6;
					connec_order_out2[7] = 7;

					connec_order_in3[0] = 0;
					connec_order_in3[1] = 3;
					connec_order_in3[2] = 1;
					connec_order_in3[3] = 2;
					connec_order_in3[4] = 4;
					connec_order_in3[5] = 5;
					connec_order_in3[6] = 6;
					connec_order_in3[7] = 7;

					connec_order_out3[0] = 0;
					connec_order_out3[1] = 3;
					connec_order_out3[2] = 0;
					connec_order_out3[3] = 1;
					connec_order_out3[4] = 2;
					connec_order_out3[5] = 4;
					connec_order_out3[6] = 5;
					connec_order_out3[7] = 3;

					elem->pad = 2;

				}//Edge 4
				else if ((edge_list[0]) && (!edge_list[1]) && (!edge_list[2]) && (edge_list[3]) &&
						(!edge_list[4]) && (!edge_list[5]) && (!edge_list[6]) && (!edge_list[7]) &&
						(edge_list[8]) && (!edge_list[9]) && (!edge_list[10]) && (edge_list[11])&&
						(edge_list_s[8]) && (edge_list_s[10])) {
					//printf("Entrou na 4!\n");

					cut_edge[0] = 0;
					cut_edge[1] = 3;
					cut_edge[2] = 11;
					cut_edge[3] = 8;

					cut_edge_s[0] = 10;
					cut_edge_s[1] = 8;

					connec_order_in1[0] = 1;
					connec_order_in1[1] = 2;
					connec_order_in1[2] = 5;
					connec_order_in1[3] = 6;
					connec_order_in1[4] = 0;
					connec_order_in1[5] = 3;
					connec_order_in1[6] = 4;
					connec_order_in1[7] = 7;

					connec_order_out1[0] = 1;
					connec_order_out1[1] = 2;
					connec_order_out1[2] = 5;
					connec_order_out1[3] = 6;
					connec_order_out1[4] = 0;
					connec_order_out1[5] = 4;
					connec_order_out1[6] = 3;
					connec_order_out1[7] = 5;

					connec_order_in2[0] = 0;
					connec_order_in2[1] = 1;
					connec_order_in2[2] = 4;
					connec_order_in2[3] = 5;
					connec_order_in2[4] = 2;
					connec_order_in2[5] = 3;
					connec_order_in2[6] = 6;
					connec_order_in2[7] = 7;

					connec_order_out2[0] = 1;
					connec_order_out2[1] = 4;
					connec_order_out2[2] = 2;
					connec_order_out2[3] = 5;
					connec_order_out2[4] = 2;
					connec_order_out2[5] = 3;
					connec_order_out2[6] = 6;
					connec_order_out2[7] = 7;

					connec_order_in3[0] = 0;
					connec_order_in3[1] = 4;
					connec_order_in3[2] = 1;
					connec_order_in3[3] = 2;
					connec_order_in3[4] = 3;
					connec_order_in3[5] = 5;
					connec_order_in3[6] = 6;
					connec_order_in3[7] = 7;

					connec_order_out3[0] = 0;
					connec_order_out3[1] = 4;
					connec_order_out3[2] = 0;
					connec_order_out3[3] = 4;
					connec_order_out3[4] = 1;
					connec_order_out3[5] = 3;
					connec_order_out3[6] = 5;
					connec_order_out3[7] = 2;

					elem->pad = 2;

				}//Edge 5
				else if ((edge_list[0]) && (edge_list[1]) && (!edge_list[2]) && (!edge_list[3]) &&
						(!edge_list[4]) && (!edge_list[5]) && (!edge_list[6]) && (!edge_list[7]) &&
						(edge_list[8]) && (edge_list[9]) && (!edge_list[10]) && (!edge_list[11])&&
						(edge_list_s[9]) && (edge_list_s[11])) {
					//printf("Entrou na 5!\n");

					cut_edge[0] = 1;
					cut_edge[1] = 0;
					cut_edge[2] = 8;
					cut_edge[3] = 9;

					cut_edge_s[0] = 11;
					cut_edge_s[1] = 9;

					connec_order_in1[0] = 0;
					connec_order_in1[1] = 3;
					connec_order_in1[2] = 4;
					connec_order_in1[3] = 7;
					connec_order_in1[4] = 1;
					connec_order_in1[5] = 2;
					connec_order_in1[6] = 5;
					connec_order_in1[7] = 6;

					connec_order_out1[0] = 0;
					connec_order_out1[1] = 3;
					connec_order_out1[2] = 4;
					connec_order_out1[3] = 7;
					connec_order_out1[4] = 1;
					connec_order_out1[5] = 4;
					connec_order_out1[6] = 2;
					connec_order_out1[7] = 5;

					connec_order_in2[0] = 0;
					connec_order_in2[1] = 1;
					connec_order_in2[2] = 4;
					connec_order_in2[3] = 5;
					connec_order_in2[4] = 2;
					connec_order_in2[5] = 3;
					connec_order_in2[6] = 6;
					connec_order_in2[7] = 7;

					connec_order_out2[0] = 4;
					connec_order_out2[1] = 0;
					connec_order_out2[2] = 5;
					connec_order_out2[3] = 3;
					connec_order_out2[4] = 2;
					connec_order_out2[5] = 3;
					connec_order_out2[6] = 6;
					connec_order_out2[7] = 7;

					connec_order_in3[0] = 1;
					connec_order_in3[1] = 5;
					connec_order_in3[2] = 0;
					connec_order_in3[3] = 2;
					connec_order_in3[4] = 3;
					connec_order_in3[5] = 4;
					connec_order_in3[6] = 6;
					connec_order_in3[7] = 7;

					connec_order_out3[0] = 1;
					connec_order_out3[1] = 5;
					connec_order_out3[2] = 1;
					connec_order_out3[3] = 0;
					connec_order_out3[4] = 4;
					connec_order_out3[5] = 2;
					connec_order_out3[6] = 3;
					connec_order_out3[7] = 5;

					elem->pad = 2;

				}//Edge 6
				else if ((!edge_list[0]) && (edge_list[1]) && (edge_list[2]) && (!edge_list[3]) &&
						(!edge_list[4]) && (!edge_list[5]) && (!edge_list[6]) && (!edge_list[7]) &&
						(!edge_list[8]) && (edge_list[9]) && (edge_list[10]) && (!edge_list[11])&&
						(edge_list_s[8]) && (edge_list_s[10])) {
					//printf("Entrou na 6!\n");

					cut_edge[0] = 2;
					cut_edge[1] = 10;
					cut_edge[2] = 1;
					cut_edge[3] = 9;

					cut_edge_s[0] = 10;
					cut_edge_s[1] = 8;

					connec_order_in1[0] = 0;
					connec_order_in1[1] = 1;
					connec_order_in1[2] = 4;
					connec_order_in1[3] = 5;
					connec_order_in1[4] = 2;
					connec_order_in1[5] = 3;
					connec_order_in1[6] = 6;
					connec_order_in1[7] = 7;

					connec_order_out1[0] = 0;
					connec_order_out1[1] = 1;
					connec_order_out1[2] = 4;
					connec_order_out1[3] = 5;
					connec_order_out1[4] = 2;
					connec_order_out1[5] = 4;
					connec_order_out1[6] = 3;
					connec_order_out1[7] = 5;

					connec_order_in2[0] = 1;
					connec_order_in2[1] = 2;
					connec_order_in2[2] = 5;
					connec_order_in2[3] = 6;
					connec_order_in2[4] = 0;
					connec_order_in2[5] = 3;
					connec_order_in2[6] = 4;
					connec_order_in2[7] = 7;

					connec_order_out2[0] = 4;
					connec_order_out2[1] = 0;
					connec_order_out2[2] = 5;
					connec_order_out2[3] = 1;
					connec_order_out2[4] = 0;
					connec_order_out2[5] = 3;
					connec_order_out2[6] = 4;
					connec_order_out2[7] = 7;

					connec_order_in3[0] = 2;
					connec_order_in3[1] = 6;
					connec_order_in3[2] = 0;
					connec_order_in3[3] = 1;
					connec_order_in3[4] = 3;
					connec_order_in3[5] = 4;
					connec_order_in3[6] = 5;
					connec_order_in3[7] = 7;

					connec_order_out3[0] = 2;
					connec_order_out3[1] = 6;
					connec_order_out3[2] = 4;
					connec_order_out3[3] = 2;
					connec_order_out3[4] = 0;
					connec_order_out3[5] = 5;
					connec_order_out3[6] = 3;
					connec_order_out3[7] = 1;

					elem->pad = 2;

				}//Edge 7
				else if ((!edge_list[0]) && (!edge_list[1]) && (edge_list[2]) && (edge_list[3]) &&
						(!edge_list[4]) && (!edge_list[5]) && (!edge_list[6]) && (!edge_list[7]) &&
						(!edge_list[8]) && (!edge_list[9]) && (edge_list[10]) && (edge_list[11])&&
						(edge_list_s[9]) && (edge_list_s[11])) {
					//printf("Entrou na 7!\n");

					cut_edge[0] = 3;
					cut_edge[1] = 11;
					cut_edge[2] = 2;
					cut_edge[3] = 10;

					cut_edge_s[0] = 11;
					cut_edge_s[1] = 9;

					connec_order_in1[0] = 0;
					connec_order_in1[1] = 1;
					connec_order_in1[2] = 4;
					connec_order_in1[3] = 5;
					connec_order_in1[4] = 2;
					connec_order_in1[5] = 3;
					connec_order_in1[6] = 6;
					connec_order_in1[7] = 7;

					connec_order_out1[0] = 0;
					connec_order_out1[1] = 1;
					connec_order_out1[2] = 4;
					connec_order_out1[3] = 5;
					connec_order_out1[4] = 4;
					connec_order_out1[5] = 0;
					connec_order_out1[6] = 5;
					connec_order_out1[7] = 1;

					connec_order_in2[0] = 0;
					connec_order_in2[1] = 3;
					connec_order_in2[2] = 4;
					connec_order_in2[3] = 7;
					connec_order_in2[4] = 1;
					connec_order_in2[5] = 2;
					connec_order_in2[6] = 5;
					connec_order_in2[7] = 6;

					connec_order_out2[0] = 4;
					connec_order_out2[1] = 2;
					connec_order_out2[2] = 5;
					connec_order_out2[3] = 3;
					connec_order_out2[4] = 1;
					connec_order_out2[5] = 2;
					connec_order_out2[6] = 5;
					connec_order_out2[7] = 6;

					connec_order_in3[0] = 3;
					connec_order_in3[1] = 7;
					connec_order_in3[2] = 0;
					connec_order_in3[3] = 1;
					connec_order_in3[4] = 2;
					connec_order_in3[5] = 4;
					connec_order_in3[6] = 5;
					connec_order_in3[7] = 6;

					connec_order_out3[0] = 3;
					connec_order_out3[1] = 7;
					connec_order_out3[2] = 0;
					connec_order_out3[3] = 4;
					connec_order_out3[4] = 2;
					connec_order_out3[5] = 1;
					connec_order_out3[6] = 5;
					connec_order_out3[7] = 3;

					elem->pad = 2;

				}//Edge 8
				else if ((!edge_list[0]) && (!edge_list[1]) && (!edge_list[2]) && (!edge_list[3]) &&
						(edge_list[4]) && (edge_list[5]) && (!edge_list[6]) && (!edge_list[7]) &&
						(!edge_list[8]) && (edge_list[9]) && (!edge_list[10]) && (edge_list[11])&&
						(edge_list_s[5]) && (edge_list_s[7])) {
					//printf("Entrou na 8!\n");

					cut_edge[0] = 11;
					cut_edge[1] = 9;
					cut_edge[2] = 4;
					cut_edge[3] = 5;

					cut_edge_s[0] = 5;
					cut_edge_s[1] = 7;

					connec_order_in1[0] = 2;
					connec_order_in1[1] = 3;
					connec_order_in1[2] = 6;
					connec_order_in1[3] = 7;
					connec_order_in1[4] = 0;
					connec_order_in1[5] = 1;
					connec_order_in1[6] = 4;
					connec_order_in1[7] = 5;

					connec_order_out1[0] = 2;
					connec_order_out1[1] = 3;
					connec_order_out1[2] = 6;
					connec_order_out1[3] = 7;
					connec_order_out1[4] = 4;
					connec_order_out1[5] = 5;
					connec_order_out1[6] = 0;
					connec_order_out1[7] = 1;

					connec_order_in2[0] = 4;
					connec_order_in2[1] = 5;
					connec_order_in2[2] = 6;
					connec_order_in2[3] = 7;
					connec_order_in2[4] = 0;
					connec_order_in2[5] = 1;
					connec_order_in2[6] = 2;
					connec_order_in2[7] = 3;

					connec_order_out2[0] = 2;
					connec_order_out2[1] = 3;
					connec_order_out2[2] = 5;
					connec_order_out2[3] = 4;
					connec_order_out2[4] = 0;
					connec_order_out2[5] = 1;
					connec_order_out2[6] = 2;
					connec_order_out2[7] = 3;

					connec_order_in3[0] = 4;
					connec_order_in3[1] = 5;
					connec_order_in3[2] = 0;
					connec_order_in3[3] = 1;
					connec_order_in3[4] = 2;
					connec_order_in3[5] = 3;
					connec_order_in3[6] = 6;
					connec_order_in3[7] = 7;

					connec_order_out3[0] = 4;
					connec_order_out3[1] = 5;
					connec_order_out3[2] = 2;
					connec_order_out3[3] = 3;
					connec_order_out3[4] = 5;
					connec_order_out3[5] = 4;
					connec_order_out3[6] = 1;
					connec_order_out3[7] = 0;

					elem->pad = 2;

				}//Edge 9
				else if ((!edge_list[0]) && (!edge_list[1]) && (!edge_list[2]) && (!edge_list[3]) &&
						(!edge_list[4]) && (edge_list[5]) && (edge_list[6]) && (!edge_list[7]) &&
						(edge_list[8]) && (!edge_list[9]) && (edge_list[10]) && (!edge_list[11])&&
						(edge_list_s[0]) && (edge_list_s[2])) {
					//printf("Entrou na 9!\n");

					cut_edge[0] = 8;
					cut_edge[1] = 10;
					cut_edge[2] = 5;
					cut_edge[3] = 6;

					cut_edge_s[0] = 0;
					cut_edge_s[1] = 2;

					connec_order_in1[0] = 0;
					connec_order_in1[1] = 3;
					connec_order_in1[2] = 4;
					connec_order_in1[3] = 7;
					connec_order_in1[4] = 1;
					connec_order_in1[5] = 2;
					connec_order_in1[6] = 5;
					connec_order_in1[7] = 6;

					connec_order_out1[0] = 0;
					connec_order_out1[1] = 3;
					connec_order_out1[2] = 4;
					connec_order_out1[3] = 7;
					connec_order_out1[4] = 4;
					connec_order_out1[5] = 5;
					connec_order_out1[6] = 0;
					connec_order_out1[7] = 1;

					connec_order_in2[0] = 4;
					connec_order_in2[1] = 5;
					connec_order_in2[2] = 6;
					connec_order_in2[3] = 7;
					connec_order_in2[4] = 0;
					connec_order_in2[5] = 1;
					connec_order_in2[6] = 2;
					connec_order_in2[7] = 3;

					connec_order_out2[0] = 4;
					connec_order_out2[1] = 2;
					connec_order_out2[2] = 3;
					connec_order_out2[3] = 5;
					connec_order_out2[4] = 0;
					connec_order_out2[5] = 1;
					connec_order_out2[6] = 2;
					connec_order_out2[7] = 3;

					connec_order_in3[0] = 5;
					connec_order_in3[1] = 6;
					connec_order_in3[2] = 0;
					connec_order_in3[3] = 1;
					connec_order_in3[4] = 2;
					connec_order_in3[5] = 3;
					connec_order_in3[6] = 4;
					connec_order_in3[7] = 7;

					connec_order_out3[0] = 5;
					connec_order_out3[1] = 6;
					connec_order_out3[2] = 4;
					connec_order_out3[3] = 2;
					connec_order_out3[4] = 3;
					connec_order_out3[5] = 5;
					connec_order_out3[6] = 0;
					connec_order_out3[7] = 1;

					elem->pad = 2;

				}//Edge 10
				else if ((!edge_list[0]) && (!edge_list[1]) && (!edge_list[2]) && (!edge_list[3]) &&
						(!edge_list[4]) && (!edge_list[5]) && (edge_list[6]) && (edge_list[7]) &&
						(!edge_list[8]) && (edge_list[9]) && (!edge_list[10]) && (edge_list[11])&&
						(edge_list_s[4]) && (edge_list_s[6])) {
					//printf("Entrou na 10!\n");

					cut_edge[0] = 9;
					cut_edge[1] = 11;
					cut_edge[2] = 6;
					cut_edge[3] = 7;

					cut_edge_s[0] = 6;
					cut_edge_s[1] = 4;

					connec_order_in1[0] = 0;
					connec_order_in1[1] = 1;
					connec_order_in1[2] = 4;
					connec_order_in1[3] = 5;
					connec_order_in1[4] = 2;
					connec_order_in1[5] = 3;
					connec_order_in1[6] = 6;
					connec_order_in1[7] = 7;

					connec_order_out1[0] = 0;
					connec_order_out1[1] = 1;
					connec_order_out1[2] = 4;
					connec_order_out1[3] = 5;
					connec_order_out1[4] = 4;
					connec_order_out1[5] = 5;
					connec_order_out1[6] = 0;
					connec_order_out1[7] = 1;

					connec_order_in2[0] = 4;
					connec_order_in2[1] = 5;
					connec_order_in2[2] = 6;
					connec_order_in2[3] = 7;
					connec_order_in2[4] = 0;
					connec_order_in2[5] = 1;
					connec_order_in2[6] = 2;
					connec_order_in2[7] = 3;

					connec_order_out2[0] = 5;
					connec_order_out2[1] = 4;
					connec_order_out2[2] = 2;
					connec_order_out2[3] = 3;
					connec_order_out2[4] = 0;
					connec_order_out2[5] = 1;
					connec_order_out2[6] = 2;
					connec_order_out2[7] = 3;

					connec_order_in3[0] = 6;
					connec_order_in3[1] = 7;
					connec_order_in3[2] = 0;
					connec_order_in3[3] = 1;
					connec_order_in3[4] = 2;
					connec_order_in3[5] = 3;
					connec_order_in3[6] = 4;
					connec_order_in3[7] = 5;

					connec_order_out3[0] = 6;
					connec_order_out3[1] = 7;
					connec_order_out3[2] = 5;
					connec_order_out3[3] = 4;
					connec_order_out3[4] = 2;
					connec_order_out3[5] = 3;
					connec_order_out3[6] = 1;
					connec_order_out3[7] = 0;

					elem->pad = 2;

				}//Edge 11
				else if ((!edge_list[0]) && (!edge_list[1]) && (!edge_list[2]) && (!edge_list[3]) &&
						(edge_list[4]) && (!edge_list[5]) && (!edge_list[6]) && (edge_list[7]) &&
						(edge_list[8]) && (!edge_list[9]) && (edge_list[10]) && (!edge_list[11])&&
						(edge_list_s[1]) && (edge_list_s[3])) {
					//printf("Entrou na 11!\n");

					cut_edge[0] = 8;
					cut_edge[1] = 10;
					cut_edge[2] = 4;
					cut_edge[3] = 7;

					cut_edge_s[0] = 1;
					cut_edge_s[1] = 3;

					connec_order_in1[0] = 1;
					connec_order_in1[1] = 2;
					connec_order_in1[2] = 5;
					connec_order_in1[3] = 6;
					connec_order_in1[4] = 0;
					connec_order_in1[5] = 3;
					connec_order_in1[6] = 4;
					connec_order_in1[7] = 7;

					connec_order_out1[0] = 1;
					connec_order_out1[1] = 2;
					connec_order_out1[2] = 5;
					connec_order_out1[3] = 6;
					connec_order_out1[4] = 4;
					connec_order_out1[5] = 5;
					connec_order_out1[6] = 0;
					connec_order_out1[7] = 1;

					connec_order_in2[0] = 4;
					connec_order_in2[1] = 5;
					connec_order_in2[2] = 6;
					connec_order_in2[3] = 7;
					connec_order_in2[4] = 0;
					connec_order_in2[5] = 1;
					connec_order_in2[6] = 2;
					connec_order_in2[7] = 3;

					connec_order_out2[0] = 2;
					connec_order_out2[1] = 4;
					connec_order_out2[2] = 5;
					connec_order_out2[3] = 3;
					connec_order_out2[4] = 0;
					connec_order_out2[5] = 1;
					connec_order_out2[6] = 2;
					connec_order_out2[7] = 3;

					connec_order_in3[0] = 4;
					connec_order_in3[1] = 7;
					connec_order_in3[2] = 0;
					connec_order_in3[3] = 1;
					connec_order_in3[4] = 2;
					connec_order_in3[5] = 3;
					connec_order_in3[6] = 5;
					connec_order_in3[7] = 6;

					connec_order_out3[0] = 4;
					connec_order_out3[1] = 7;
					connec_order_out3[2] = 2;
					connec_order_out3[3] = 4;
					connec_order_out3[4] = 5;
					connec_order_out3[5] = 3;
					connec_order_out3[6] = 0;
					connec_order_out3[7] = 1;

					elem->pad = 2;

				} else {
					elem->pad = -1;
					printf("Warning: Bad Element: %d, case 2\n",elements_ids[iel]);
					printf("Edges:  Elem. %d: %d %d %d %d %d %d %d %d %d %d %d %d %d\n", elements_ids[iel], edge_list[0],
							edge_list[1], edge_list[2], edge_list[3],
							edge_list[4], edge_list[5], edge_list[6],
							edge_list[7], edge_list[8], edge_list[9],
							edge_list[10], edge_list[11], ed_cont);
					continue;
				}

				GtsPoint *p0 = point[cut_edge[0]];
				GtsPoint *p1 = point[cut_edge[1]];
				GtsPoint *p2 = point[cut_edge[2]];
				GtsPoint *p3 = point[cut_edge[3]];

				GtsPoint *p4 = point_s[cut_edge_s[0]];
				GtsPoint *p5 = point_s[cut_edge_s[1]];

				//g_assert(p0 != NULL);
				//g_assert(p1 != NULL);
				//g_assert(p2 != NULL);
				//g_assert(p3 != NULL);
				//g_assert(p4 != NULL);
				//g_assert(p5 != NULL);

				if(p0==NULL) exit(1);
				if(p1==NULL) exit(2);
				if(p2==NULL) exit(3);
				if(p3==NULL) exit(4);
				if(p4==NULL) exit(5);
				if(p5==NULL) exit(6);

				conn_t2[0] = AddPointOnEdge(Edge2GNode[cut_edge[0]], hash_nodes, mesh->local_n_nodes, p0, coords);
				conn_t2[1] = AddPointOnEdge(Edge2GNode[cut_edge[1]], hash_nodes, mesh->local_n_nodes, p1, coords);
				conn_t2[2] = AddPointOnEdge(Edge2GNode[cut_edge[2]], hash_nodes, mesh->local_n_nodes, p2, coords);
				conn_t2[3] = AddPointOnEdge(Edge2GNode[cut_edge[3]], hash_nodes, mesh->local_n_nodes, p3, coords);

				// add 2 extra points in the surface
				conn_t2[4] = AddPointOnEdge(Edge2GNode_s[cut_edge_s[0]], hash_nodes, mesh->local_n_nodes, p4, coords);
				conn_t2[5] = AddPointOnEdge(Edge2GNode_s[cut_edge_s[1]], hash_nodes, mesh->local_n_nodes, p5, coords);

				octant_t *elem1 = (octant_t*) sc_array_index(&mesh->elements, elements_ids[iel]);

				elem1->nodes[connec_order_in1[0]].id = original_conn[connec_order_out1[0]];
				elem1->nodes[connec_order_in1[1]].id = original_conn[connec_order_out1[1]];
				elem1->nodes[connec_order_in1[2]].id = original_conn[connec_order_out1[2]];
				elem1->nodes[connec_order_in1[3]].id = original_conn[connec_order_out1[3]];

				elem1->nodes[connec_order_in1[4]].id = conn_t2[connec_order_out1[4]];
				elem1->nodes[connec_order_in1[5]].id = conn_t2[connec_order_out1[5]];
				elem1->nodes[connec_order_in1[6]].id = conn_t2[connec_order_out1[6]];
				elem1->nodes[connec_order_in1[7]].id = conn_t2[connec_order_out1[7]];

				octant_t* elem2 = (octant_t*) sc_array_push(&mesh->elements);

				elem2->nodes[connec_order_in2[0]].id = conn_t2[connec_order_out2[0]];
				elem2->nodes[connec_order_in2[1]].id = conn_t2[connec_order_out2[1]];
				elem2->nodes[connec_order_in2[2]].id = conn_t2[connec_order_out2[2]];
				elem2->nodes[connec_order_in2[3]].id = conn_t2[connec_order_out2[3]];

				elem2->nodes[connec_order_in2[4]].id = original_conn[connec_order_out2[4]];
				elem2->nodes[connec_order_in2[5]].id = original_conn[connec_order_out2[5]];
				elem2->nodes[connec_order_in2[6]].id = original_conn[connec_order_out2[6]];
				elem2->nodes[connec_order_in2[7]].id = original_conn[connec_order_out2[7]];

				elem2->pad = elem->pad;
				elem2->level = elem->level;

				octant_t* elem3 = (octant_t*) sc_array_push(&mesh->elements);

				elem3->nodes[connec_order_in3[0]].id = original_conn[connec_order_out3[0]];
				elem3->nodes[connec_order_in3[1]].id = original_conn[connec_order_out3[1]];

				elem3->nodes[connec_order_in3[2]].id = conn_t2[connec_order_out3[2]];
				elem3->nodes[connec_order_in3[3]].id = conn_t2[connec_order_out3[3]];
				elem3->nodes[connec_order_in3[4]].id = conn_t2[connec_order_out3[4]];
				elem3->nodes[connec_order_in3[5]].id = conn_t2[connec_order_out3[5]];
				elem3->nodes[connec_order_in3[6]].id = conn_t2[connec_order_out3[6]];
				elem3->nodes[connec_order_in3[7]].id = conn_t2[connec_order_out3[7]];

				elem3->pad = elem->pad;
				elem3->level = elem->level;
			}
		}else if(elem->pad==3){

			if (n_parallel_faces == 0 && ed_cont == 3) {
				// Check template 3.

				int conn_t2[7];
				int cut_edge[3];
				int cut_edge_s[3];
				int cut_edge_v[1];
				int connec_order_in1[8];
				int connec_order_out1[8];
				int connec_order_in2[8];
				int connec_order_out2[8];
				int connec_order_in3[8];
				int connec_order_out3[8];
				int connec_order_in4[8];
				int connec_order_out4[8];


				//Corner 0
				if ((edge_list[0]) && (!edge_list[1]) && (!edge_list[2]) && (edge_list[3]) &&
						(edge_list[4]) && (!edge_list[5]) && (!edge_list[6]) && (!edge_list[7]) &&
						(!edge_list[8]) && (!edge_list[9]) && (!edge_list[10]) && (!edge_list[11]) &&
						(edge_list_s[4]) && (edge_list_s[0]) && (edge_list_s[10]) &&
						(edge_list_v[0]) ) {
					//printf("Entrou no 0!\n");

					cut_edge[0] = 4;
					cut_edge[1] = 0;
					cut_edge[2] = 3;

					cut_edge_s[0] = 0;
					cut_edge_s[1] = 4;
					cut_edge_s[2] = 10;

					cut_edge_v[0] = 0;

					connec_order_in1[0] = 4;
					connec_order_in1[1] = 5;
					connec_order_in1[2] = 6;
					connec_order_in1[3] = 7;
					connec_order_in1[4] = 0;
					connec_order_in1[5] = 1;
					connec_order_in1[6] = 2;
					connec_order_in1[7] = 3;

					connec_order_out1[0] = 4;
					connec_order_out1[1] = 5;
					connec_order_out1[2] = 6;
					connec_order_out1[3] = 7;
					connec_order_out1[4] = 0;
					connec_order_out1[5] = 3;
					connec_order_out1[6] = 6;
					connec_order_out1[7] = 4;

					connec_order_in2[0] = 0;
					connec_order_in2[1] = 3;
					connec_order_in2[2] = 4;
					connec_order_in2[3] = 7;
					connec_order_in2[4] = 1;
					connec_order_in2[5] = 2;
					connec_order_in2[6] = 5;
					connec_order_in2[7] = 6;

					connec_order_out2[0] = 1;
					connec_order_out2[1] = 5;
					connec_order_out2[2] = 3;
					connec_order_out2[3] = 6;
					connec_order_out2[4] = 1;
					connec_order_out2[5] = 2;
					connec_order_out2[6] = 5;
					connec_order_out2[7] = 6;

					connec_order_in3[0] = 0;
					connec_order_in3[1] = 1;
					connec_order_in3[2] = 4;
					connec_order_in3[3] = 5;
					connec_order_in3[4] = 2;
					connec_order_in3[5] = 3;
					connec_order_in3[6] = 6;
					connec_order_in3[7] = 7;

					connec_order_out3[0] = 2;
					connec_order_out3[1] = 5;
					connec_order_out3[2] = 4;
					connec_order_out3[3] = 6;
					connec_order_out3[4] = 2;
					connec_order_out3[5] = 3;
					connec_order_out3[6] = 6;
					connec_order_out3[7] = 7;

					connec_order_in4[0] = 0;
					connec_order_in4[1] = 1;
					connec_order_in4[2] = 2;
					connec_order_in4[3] = 3;
					connec_order_in4[4] = 4;
					connec_order_in4[5] = 5;
					connec_order_in4[6] = 6;
					connec_order_in4[7] = 7;

					connec_order_out4[0] = 0;
					connec_order_out4[1] = 1;
					connec_order_out4[2] = 5;
					connec_order_out4[3] = 2;
					connec_order_out4[4] = 0;
					connec_order_out4[5] = 3;
					connec_order_out4[6] = 6;
					connec_order_out4[7] = 4;

					elem->pad = 3;

				}

				//Corner 1
				else if ((edge_list[0]) && (edge_list[1]) && (!edge_list[2]) && (!edge_list[3]) &&
						(!edge_list[4]) && (edge_list[5]) && (!edge_list[6]) && (!edge_list[7]) &&
						(!edge_list[8]) && (!edge_list[9]) && (!edge_list[10]) && (!edge_list[11])&&
						(edge_list_s[1]) && (edge_list_s[6]) && (edge_list_s[11])  &&
						(edge_list_v[1]) ) {
					//printf("Entrou no 1!\n");

					cut_edge[0] = 5;
					cut_edge[1] = 1;
					cut_edge[2] = 0;

					cut_edge_s[0] = 6;
					cut_edge_s[1] = 1;
					cut_edge_s[2] = 11;

					cut_edge_v[0] = 1;

					connec_order_in1[0] = 4;
					connec_order_in1[1] = 5;
					connec_order_in1[2] = 6;
					connec_order_in1[3] = 7;
					connec_order_in1[4] = 0;
					connec_order_in1[5] = 1;
					connec_order_in1[6] = 2;
					connec_order_in1[7] = 3;

					connec_order_out1[0] = 4;
					connec_order_out1[1] = 5;
					connec_order_out1[2] = 6;
					connec_order_out1[3] = 7;
					connec_order_out1[4] = 4;
					connec_order_out1[5] = 0;
					connec_order_out1[6] = 3;
					connec_order_out1[7] = 6;

					connec_order_in2[0] = 0;
					connec_order_in2[1] = 1;
					connec_order_in2[2] = 4;
					connec_order_in2[3] = 5;
					connec_order_in2[4] = 2;
					connec_order_in2[5] = 3;
					connec_order_in2[6] = 6;
					connec_order_in2[7] = 7;

					connec_order_out2[0] = 5;
					connec_order_out2[1] = 1;
					connec_order_out2[2] = 6;
					connec_order_out2[3] = 3;
					connec_order_out2[4] = 2;
					connec_order_out2[5] = 3;
					connec_order_out2[6] = 6;
					connec_order_out2[7] = 7;

					connec_order_in3[0] = 1;
					connec_order_in3[1] = 2;
					connec_order_in3[2] = 5;
					connec_order_in3[3] = 6;
					connec_order_in3[4] = 0;
					connec_order_in3[5] = 3;
					connec_order_in3[6] = 4;
					connec_order_in3[7] = 7;

					connec_order_out3[0] = 2;
					connec_order_out3[1] = 5;
					connec_order_out3[2] = 4;
					connec_order_out3[3] = 6;
					connec_order_out3[4] = 0;
					connec_order_out3[5] = 3;
					connec_order_out3[6] = 4;
					connec_order_out3[7] = 7;

					connec_order_in4[0] = 1;
					connec_order_in4[1] = 0;
					connec_order_in4[2] = 2;
					connec_order_in4[3] = 3;
					connec_order_in4[4] = 4;
					connec_order_in4[5] = 5;
					connec_order_in4[6] = 6;
					connec_order_in4[7] = 7;

					connec_order_out4[0] = 1;
					connec_order_out4[1] = 2;
					connec_order_out4[2] = 1;
					connec_order_out4[3] = 5;
					connec_order_out4[4] = 4;
					connec_order_out4[5] = 0;
					connec_order_out4[6] = 3;
					connec_order_out4[7] = 6;

					elem->pad = 3;

				}

				//Corner 2
				else if ((!edge_list[0]) && (edge_list[1]) && (edge_list[2]) && (!edge_list[3]) &&
						(!edge_list[4]) && (!edge_list[5]) && (edge_list[6]) && (!edge_list[7]) &&
						(!edge_list[8]) && (!edge_list[9]) && (!edge_list[10]) && (!edge_list[11]) &&
						(edge_list_s[3]) && (edge_list_s[7]) && (edge_list_s[10])  &&
						(edge_list_v[2]) ) {
					//printf("Entrou no 2!\n");

					cut_edge[0] = 6;
					cut_edge[1] = 2;
					cut_edge[2] = 1;

					cut_edge_s[0] = 3;
					cut_edge_s[1] = 7;
					cut_edge_s[2] = 10;

					cut_edge_v[0] = 2;

					connec_order_in1[0] = 4;
					connec_order_in1[1] = 5;
					connec_order_in1[2] = 6;
					connec_order_in1[3] = 7;
					connec_order_in1[4] = 1;
					connec_order_in1[5] = 2;
					connec_order_in1[6] = 3;
					connec_order_in1[7] = 4;

					connec_order_out1[0] = 4;
					connec_order_out1[1] = 5;
					connec_order_out1[2] = 6;
					connec_order_out1[3] = 7;
					connec_order_out1[4] = 4;
					connec_order_out1[5] = 0;
					connec_order_out1[6] = 3;
					connec_order_out1[7] = 6;

					connec_order_in2[0] = 1;
					connec_order_in2[1] = 2;
					connec_order_in2[2] = 5;
					connec_order_in2[3] = 6;
					connec_order_in2[4] = 0;
					connec_order_in2[5] = 3;
					connec_order_in2[6] = 4;
					connec_order_in2[7] = 7;

					connec_order_out2[0] = 5;
					connec_order_out2[1] = 1;
					connec_order_out2[2] = 6;
					connec_order_out2[3] = 3;
					connec_order_out2[4] = 0;
					connec_order_out2[5] = 3;
					connec_order_out2[6] = 4;
					connec_order_out2[7] = 7;

					connec_order_in3[0] = 2;
					connec_order_in3[1] = 3;
					connec_order_in3[2] = 6;
					connec_order_in3[3] = 7;
					connec_order_in3[4] = 0;
					connec_order_in3[5] = 1;
					connec_order_in3[6] = 4;
					connec_order_in3[7] = 5;

					connec_order_out3[0] = 2;
					connec_order_out3[1] = 5;
					connec_order_out3[2] = 4;
					connec_order_out3[3] = 6;
					connec_order_out3[4] = 0;
					connec_order_out3[5] = 1;
					connec_order_out3[6] = 4;
					connec_order_out3[7] = 5;

					connec_order_in4[0] = 2;
					connec_order_in4[1] = 0;
					connec_order_in4[2] = 1;
					connec_order_in4[3] = 3;
					connec_order_in4[4] = 4;
					connec_order_in4[5] = 5;
					connec_order_in4[6] = 6;
					connec_order_in4[7] = 7;

					connec_order_out4[0] = 2;
					connec_order_out4[1] = 5;
					connec_order_out4[2] = 2;
					connec_order_out4[3] = 1;
					connec_order_out4[4] = 6;
					connec_order_out4[5] = 4;
					connec_order_out4[6] = 0;
					connec_order_out4[7] = 3;

					elem->pad = 3;

				}

				//Corner 3
				else if ((!edge_list[0]) && (!edge_list[1]) && (edge_list[2]) && (edge_list[3]) &&
						(!edge_list[4]) && (!edge_list[5]) && (!edge_list[6]) && (edge_list[7]) &&
						(!edge_list[8]) && (!edge_list[9]) && (!edge_list[10]) && (!edge_list[11]) &&
						(edge_list_s[5]) && (edge_list_s[2]) && (edge_list_s[11])  &&
						(edge_list_v[3]) ) {
					//printf("Entrou no 3!\n");

					cut_edge[0] = 7;
					cut_edge[1] = 3;
					cut_edge[2] = 2;

					cut_edge_s[0] = 5;
					cut_edge_s[1] = 2;
					cut_edge_s[2] = 11;

					cut_edge_v[0] = 3;

					connec_order_in1[0] = 4;
					connec_order_in1[1] = 5;
					connec_order_in1[2] = 6;
					connec_order_in1[3] = 7;
					connec_order_in1[4] = 0;
					connec_order_in1[5] = 1;
					connec_order_in1[6] = 2;
					connec_order_in1[7] = 3;

					connec_order_out1[0] = 4;
					connec_order_out1[1] = 5;
					connec_order_out1[2] = 6;
					connec_order_out1[3] = 7;
					connec_order_out1[4] = 3;
					connec_order_out1[5] = 6;
					connec_order_out1[6] = 4;
					connec_order_out1[7] = 0;

					connec_order_in2[0] = 2;
					connec_order_in2[1] = 3;
					connec_order_in2[2] = 6;
					connec_order_in2[3] = 7;
					connec_order_in2[4] = 0;
					connec_order_in2[5] = 1;
					connec_order_in2[6] = 4;
					connec_order_in2[7] = 5;

					connec_order_out2[0] = 5;
					connec_order_out2[1] = 1;
					connec_order_out2[2] = 6;
					connec_order_out2[3] = 3;
					connec_order_out2[4] = 0;
					connec_order_out2[5] = 1;
					connec_order_out2[6] = 4;
					connec_order_out2[7] = 5;

					connec_order_in3[0] = 0;
					connec_order_in3[1] = 3;
					connec_order_in3[2] = 4;
					connec_order_in3[3] = 7;
					connec_order_in3[4] = 1;
					connec_order_in3[5] = 2;
					connec_order_in3[6] = 5;
					connec_order_in3[7] = 6;

					connec_order_out3[0] = 5;
					connec_order_out3[1] = 2;
					connec_order_out3[2] = 6;
					connec_order_out3[3] = 4;
					connec_order_out3[4] = 1;
					connec_order_out3[5] = 2;
					connec_order_out3[6] = 5;
					connec_order_out3[7] = 6;

					connec_order_in4[0] = 3;
					connec_order_in4[1] = 0;
					connec_order_in4[2] = 1;
					connec_order_in4[3] = 2;
					connec_order_in4[4] = 4;
					connec_order_in4[5] = 5;
					connec_order_in4[6] = 6;
					connec_order_in4[7] = 7;

					connec_order_out4[0] = 3;
					connec_order_out4[1] = 1;
					connec_order_out4[2] = 5;
					connec_order_out4[3] = 2;
					connec_order_out4[4] = 3;
					connec_order_out4[5] = 6;
					connec_order_out4[6] = 4;
					connec_order_out4[7] = 0;

					elem->pad = 3;

				}


				//Corner 4
				else if ((!edge_list[0]) && (!edge_list[1]) && (!edge_list[2]) && (!edge_list[3]) &&
						(edge_list[4]) && (!edge_list[5]) && (!edge_list[6]) && (!edge_list[7]) &&
						(edge_list[8]) && (!edge_list[9]) && (!edge_list[10]) && (edge_list[11]) &&
						(edge_list_s[1]) && (edge_list_s[5]) && (edge_list_s[8])  &&
						(edge_list_v[2]) ) {
					//printf("Entrou no 4!\n");

					cut_edge[0] = 4;
					cut_edge[1] = 8;
					cut_edge[2] = 11;

					cut_edge_s[0] = 1;
					cut_edge_s[1] = 5;
					cut_edge_s[2] = 8;

					cut_edge_v[0] = 2;

					connec_order_in1[0] = 0;
					connec_order_in1[1] = 1;
					connec_order_in1[2] = 2;
					connec_order_in1[3] = 3;
					connec_order_in1[4] = 4;
					connec_order_in1[5] = 5;
					connec_order_in1[6] = 6;
					connec_order_in1[7] = 7;

					connec_order_out1[0] = 0;
					connec_order_out1[1] = 1;
					connec_order_out1[2] = 2;
					connec_order_out1[3] = 3;
					connec_order_out1[4] = 0;
					connec_order_out1[5] = 3;
					connec_order_out1[6] = 6;
					connec_order_out1[7] = 4;

					connec_order_in2[0] = 0;
					connec_order_in2[1] = 3;
					connec_order_in2[2] = 4;
					connec_order_in2[3] = 7;
					connec_order_in2[4] = 1;
					connec_order_in2[5] = 2;
					connec_order_in2[6] = 5;
					connec_order_in2[7] = 6;

					connec_order_out2[0] = 3;
					connec_order_out2[1] = 6;
					connec_order_out2[2] = 1;
					connec_order_out2[3] = 5;
					connec_order_out2[4] = 1;
					connec_order_out2[5] = 2;
					connec_order_out2[6] = 5;
					connec_order_out2[7] = 6;

					connec_order_in3[0] = 0;
					connec_order_in3[1] = 1;
					connec_order_in3[2] = 4;
					connec_order_in3[3] = 5;
					connec_order_in3[4] = 2;
					connec_order_in3[5] = 3;
					connec_order_in3[6] = 6;
					connec_order_in3[7] = 7;

					connec_order_out3[0] = 4;
					connec_order_out3[1] = 6;
					connec_order_out3[2] = 2;
					connec_order_out3[3] = 5;
					connec_order_out3[4] = 2;
					connec_order_out3[5] = 3;
					connec_order_out3[6] = 6;
					connec_order_out3[7] = 7;

					connec_order_in4[0] = 4;
					connec_order_in4[1] = 0;
					connec_order_in4[2] = 1;
					connec_order_in4[3] = 2;
					connec_order_in4[4] = 3;
					connec_order_in4[5] = 5;
					connec_order_in4[6] = 6;
					connec_order_in4[7] = 7;

					connec_order_out4[0] = 4;
					connec_order_out4[1] = 0;
					connec_order_out4[2] = 3;
					connec_order_out4[3] = 6;
					connec_order_out4[4] = 4;
					connec_order_out4[5] = 1;
					connec_order_out4[6] = 5;
					connec_order_out4[7] = 2;

					elem->pad = 3;

				}

				//Corner 5
				else if ((!edge_list[0]) && (!edge_list[1]) && (!edge_list[2]) && (!edge_list[3]) &&
						(!edge_list[4]) && (edge_list[5]) && (!edge_list[6]) && (!edge_list[7]) &&
						(edge_list[8]) && (edge_list[9]) && (!edge_list[10]) && (!edge_list[11]) &&
						(edge_list_s[7]) && (edge_list_s[0]) && (edge_list_s[9])  &&
						(edge_list_v[3]) ) {
					//printf("Entrou no 5!\n");

					cut_edge[0] = 5;
					cut_edge[1] = 9;
					cut_edge[2] = 8;

					cut_edge_s[0] = 7;
					cut_edge_s[1] = 0;
					cut_edge_s[2] = 9;

					cut_edge_v[0] = 3;

					connec_order_in1[0] = 0;
					connec_order_in1[1] = 1;
					connec_order_in1[2] = 2;
					connec_order_in1[3] = 3;
					connec_order_in1[4] = 4;
					connec_order_in1[5] = 5;
					connec_order_in1[6] = 6;
					connec_order_in1[7] = 7;

					connec_order_out1[0] = 0;
					connec_order_out1[1] = 1;
					connec_order_out1[2] = 2;
					connec_order_out1[3] = 3;
					connec_order_out1[4] = 4;
					connec_order_out1[5] = 0;
					connec_order_out1[6] = 3;
					connec_order_out1[7] = 6;

					connec_order_in2[0] = 0;
					connec_order_in2[1] = 1;
					connec_order_in2[2] = 4;
					connec_order_in2[3] = 5;
					connec_order_in2[4] = 2;
					connec_order_in2[5] = 3;
					connec_order_in2[6] = 6;
					connec_order_in2[7] = 7;

					connec_order_out2[0] = 6;
					connec_order_out2[1] = 3;
					connec_order_out2[2] = 5;
					connec_order_out2[3] = 1;
					connec_order_out2[4] = 2;
					connec_order_out2[5] = 3;
					connec_order_out2[6] = 6;
					connec_order_out2[7] = 7;

					connec_order_in3[0] = 1;
					connec_order_in3[1] = 2;
					connec_order_in3[2] = 5;
					connec_order_in3[3] = 6;
					connec_order_in3[4] = 0;
					connec_order_in3[5] = 3;
					connec_order_in3[6] = 4;
					connec_order_in3[7] = 7;

					connec_order_out3[0] = 4;
					connec_order_out3[1] = 6;
					connec_order_out3[2] = 2;
					connec_order_out3[3] = 5;
					connec_order_out3[4] = 0;
					connec_order_out3[5] = 3;
					connec_order_out3[6] = 4;
					connec_order_out3[7] = 7;

					connec_order_in4[0] = 5;
					connec_order_in4[1] = 0;
					connec_order_in4[2] = 1;
					connec_order_in4[3] = 2;
					connec_order_in4[4] = 3;
					connec_order_in4[5] = 4;
					connec_order_in4[6] = 6;
					connec_order_in4[7] = 7;

					connec_order_out4[0] = 5;
					connec_order_out4[1] = 4;
					connec_order_out4[2] = 0;
					connec_order_out4[3] = 3;
					connec_order_out4[4] = 6;
					connec_order_out4[5] = 2;
					connec_order_out4[6] = 1;
					connec_order_out4[7] = 5;

					elem->pad = 3;

				}

				//Corner 6
				else if ((!edge_list[0]) && (!edge_list[1]) && (!edge_list[2]) && (!edge_list[3]) &&
						(!edge_list[4]) && (!edge_list[5]) && (edge_list[6]) && (!edge_list[7]) &&
						(!edge_list[8]) && (edge_list[9]) && (edge_list[10]) && (!edge_list[11]) &&
						(edge_list_s[2]) && (edge_list_s[6]) && (edge_list_s[8])  &&
						(edge_list_v[0]) ) {
					//printf("Entrou no 6!\n");


					cut_edge[0] = 6;
					cut_edge[1] = 10;
					cut_edge[2] = 9;

					cut_edge_s[0] = 2;
					cut_edge_s[1] = 6;
					cut_edge_s[2] = 8;

					cut_edge_v[0] = 0;

					connec_order_in1[0] = 0;
					connec_order_in1[1] = 1;
					connec_order_in1[2] = 2;
					connec_order_in1[3] = 3;
					connec_order_in1[4] = 4;
					connec_order_in1[5] = 5;
					connec_order_in1[6] = 6;
					connec_order_in1[7] = 7;

					connec_order_out1[0] = 0;
					connec_order_out1[1] = 1;
					connec_order_out1[2] = 2;
					connec_order_out1[3] = 3;
					connec_order_out1[4] = 6;
					connec_order_out1[5] = 4;
					connec_order_out1[6] = 0;
					connec_order_out1[7] = 3;

					connec_order_in2[0] = 1;
					connec_order_in2[1] = 2;
					connec_order_in2[2] = 5;
					connec_order_in2[3] = 6;
					connec_order_in2[4] = 0;
					connec_order_in2[5] = 3;
					connec_order_in2[6] = 4;
					connec_order_in2[7] = 7;

					connec_order_out2[0] = 6;
					connec_order_out2[1] = 3;
					connec_order_out2[2] = 5;
					connec_order_out2[3] = 1;
					connec_order_out2[4] = 0;
					connec_order_out2[5] = 3;
					connec_order_out2[6] = 4;
					connec_order_out2[7] = 7;

					connec_order_in3[0] = 2;
					connec_order_in3[1] = 3;
					connec_order_in3[2] = 6;
					connec_order_in3[3] = 7;
					connec_order_in3[4] = 0;
					connec_order_in3[5] = 1;
					connec_order_in3[6] = 4;
					connec_order_in3[7] = 5;

					connec_order_out3[0] = 4;
					connec_order_out3[1] = 6;
					connec_order_out3[2] = 2;
					connec_order_out3[3] = 5;
					connec_order_out3[4] = 0;
					connec_order_out3[5] = 1;
					connec_order_out3[6] = 4;
					connec_order_out3[7] = 5;

					connec_order_in4[0] = 6;
					connec_order_in4[1] = 0;
					connec_order_in4[2] = 1;
					connec_order_in4[3] = 2;
					connec_order_in4[4] = 3;
					connec_order_in4[5] = 4;
					connec_order_in4[6] = 5;
					connec_order_in4[7] = 7;

					connec_order_out4[0] = 6;
					connec_order_out4[1] = 6;
					connec_order_out4[2] = 4;
					connec_order_out4[3] = 0;
					connec_order_out4[4] = 3;
					connec_order_out4[5] = 5;
					connec_order_out4[6] = 2;
					connec_order_out4[7] = 1;

					elem->pad = 3;

				}

				//Corner 7
				else if ((!edge_list[0]) && (!edge_list[1]) && (!edge_list[2]) && (!edge_list[3]) &&
						(!edge_list[4]) && (!edge_list[5]) && (!edge_list[6]) && (edge_list[7]) &&
						(!edge_list[8]) && (!edge_list[9]) && (edge_list[10]) && (edge_list[11]) &&
						(edge_list_s[4]) && (edge_list_s[3]) && (edge_list_s[9])  &&
						(edge_list_v[1]) ) {
					//printf("Entrou no 7!\n");


					cut_edge[0] = 7;
					cut_edge[1] = 11;
					cut_edge[2] = 10;

					cut_edge_s[0] = 4;
					cut_edge_s[1] = 3;
					cut_edge_s[2] = 9;

					cut_edge_v[0] = 1;

					connec_order_in1[0] = 0;
					connec_order_in1[1] = 1;
					connec_order_in1[2] = 2;
					connec_order_in1[3] = 3;
					connec_order_in1[4] = 4;
					connec_order_in1[5] = 5;
					connec_order_in1[6] = 6;
					connec_order_in1[7] = 7;

					connec_order_out1[0] = 0;
					connec_order_out1[1] = 1;
					connec_order_out1[2] = 2;
					connec_order_out1[3] = 3;
					connec_order_out1[4] = 3;
					connec_order_out1[5] = 6;
					connec_order_out1[6] = 4;
					connec_order_out1[7] = 0;

					connec_order_in2[0] = 2;
					connec_order_in2[1] = 3;
					connec_order_in2[2] = 6;
					connec_order_in2[3] = 7;
					connec_order_in2[4] = 0;
					connec_order_in2[5] = 1;
					connec_order_in2[6] = 4;
					connec_order_in2[7] = 5;

					connec_order_out2[0] = 6;
					connec_order_out2[1] = 3;
					connec_order_out2[2] = 5;
					connec_order_out2[3] = 1;
					connec_order_out2[4] = 0;
					connec_order_out2[5] = 1;
					connec_order_out2[6] = 4;
					connec_order_out2[7] = 5;

					connec_order_in3[0] = 0;
					connec_order_in3[1] = 3;
					connec_order_in3[2] = 4;
					connec_order_in3[3] = 7;
					connec_order_in3[4] = 1;
					connec_order_in3[5] = 2;
					connec_order_in3[6] = 5;
					connec_order_in3[7] = 6;

					connec_order_out3[0] = 6;
					connec_order_out3[1] = 4;
					connec_order_out3[2] = 5;
					connec_order_out3[3] = 2;
					connec_order_out3[4] = 1;
					connec_order_out3[5] = 2;
					connec_order_out3[6] = 5;
					connec_order_out3[7] = 6;

					connec_order_in4[0] = 7;
					connec_order_in4[1] = 0;
					connec_order_in4[2] = 1;
					connec_order_in4[3] = 2;
					connec_order_in4[4] = 3;
					connec_order_in4[5] = 4;
					connec_order_in4[6] = 5;
					connec_order_in4[7] = 6;

					connec_order_out4[0] = 7;
					connec_order_out4[1] = 3;
					connec_order_out4[2] = 6;
					connec_order_out4[3] = 4;
					connec_order_out4[4] = 0;
					connec_order_out4[5] = 1;
					connec_order_out4[6] = 5;
					connec_order_out4[7] = 2;

					elem->pad = 3;

				} else {
					elem->pad = -1;
					printf("Warning: Bad Element: %d, case 3\n",elements_ids[iel]);
					printf("Edges:  Elem. %d: %d %d %d %d %d %d %d %d %d %d %d %d %d\n", elements_ids[iel], edge_list[0],
							edge_list[1], edge_list[2], edge_list[3],
							edge_list[4], edge_list[5], edge_list[6],
							edge_list[7], edge_list[8], edge_list[9],
							edge_list[10], edge_list[11], ed_cont);
					continue;
				}

				GtsPoint *p0 = point[cut_edge[0]];
				GtsPoint *p1 = point[cut_edge[1]];
				GtsPoint *p2 = point[cut_edge[2]];

				GtsPoint *p3 = point_s[cut_edge_s[0]];
				GtsPoint *p4 = point_s[cut_edge_s[1]];
				GtsPoint *p5 = point_s[cut_edge_s[2]];

				GtsPoint *p6 = point_v[cut_edge_v[0]];

				//g_assert(p0 != NULL);
				//g_assert(p1 != NULL);
				//g_assert(p2 != NULL);
				//g_assert(p3 != NULL);
				//g_assert(p4 != NULL);
				//g_assert(p5 != NULL);
				//g_assert(p6 != NULL);

				if(p0==NULL) exit(1);
				if(p1==NULL) exit(2);
				if(p2==NULL) exit(3);
				if(p3==NULL) exit(4);
				if(p4==NULL) exit(5);
				if(p5==NULL) exit(6);
				if(p6==NULL) exit(7);

				conn_t2[0] = AddPointOnEdge(Edge2GNode[cut_edge[0]], hash_nodes, mesh->local_n_nodes, p0, coords);
				conn_t2[1] = AddPointOnEdge(Edge2GNode[cut_edge[1]], hash_nodes, mesh->local_n_nodes, p1, coords);
				conn_t2[2] = AddPointOnEdge(Edge2GNode[cut_edge[2]], hash_nodes, mesh->local_n_nodes, p2, coords);

				// add 3 extra points in the surface
				conn_t2[3] = AddPointOnEdge(Edge2GNode_s[cut_edge_s[0]], hash_nodes, mesh->local_n_nodes, p3, coords);
				conn_t2[4] = AddPointOnEdge(Edge2GNode_s[cut_edge_s[1]], hash_nodes, mesh->local_n_nodes, p4, coords);
				conn_t2[5] = AddPointOnEdge(Edge2GNode_s[cut_edge_s[2]], hash_nodes, mesh->local_n_nodes, p5, coords);

				// add 1 extra points in the volume
				conn_t2[6] = AddPointOnEdge(Edge2GNode_s[cut_edge_v[0]], hash_nodes, mesh->local_n_nodes, p6, coords);

				octant_t *elem1 = (octant_t*) sc_array_index(&mesh->elements, elements_ids[iel]);

				elem1->nodes[connec_order_in1[0]].id = original_conn[connec_order_out1[0]];
				elem1->nodes[connec_order_in1[1]].id = original_conn[connec_order_out1[1]];
				elem1->nodes[connec_order_in1[2]].id = original_conn[connec_order_out1[2]];
				elem1->nodes[connec_order_in1[3]].id = original_conn[connec_order_out1[3]];

				elem1->nodes[connec_order_in1[4]].id = conn_t2[connec_order_out1[4]];
				elem1->nodes[connec_order_in1[5]].id = conn_t2[connec_order_out1[5]];
				elem1->nodes[connec_order_in1[6]].id = conn_t2[connec_order_out1[6]];
				elem1->nodes[connec_order_in1[7]].id = conn_t2[connec_order_out1[7]];

				octant_t* elem2 = (octant_t*) sc_array_push(&mesh->elements);

				elem2->nodes[connec_order_in2[0]].id = conn_t2[connec_order_out2[0]];
				elem2->nodes[connec_order_in2[1]].id = conn_t2[connec_order_out2[1]];
				elem2->nodes[connec_order_in2[2]].id = conn_t2[connec_order_out2[2]];
				elem2->nodes[connec_order_in2[3]].id = conn_t2[connec_order_out2[3]];

				elem2->nodes[connec_order_in2[4]].id = original_conn[connec_order_out2[4]];
				elem2->nodes[connec_order_in2[5]].id = original_conn[connec_order_out2[5]];
				elem2->nodes[connec_order_in2[6]].id = original_conn[connec_order_out2[6]];
				elem2->nodes[connec_order_in2[7]].id = original_conn[connec_order_out2[7]];

				elem2->pad = elem->pad;
				elem2->level = elem->level;

				octant_t* elem3 = (octant_t*) sc_array_push(&mesh->elements);

				elem3->nodes[connec_order_in3[0]].id = conn_t2[connec_order_out3[0]];
				elem3->nodes[connec_order_in3[1]].id = conn_t2[connec_order_out3[1]];
				elem3->nodes[connec_order_in3[2]].id = conn_t2[connec_order_out3[2]];
				elem3->nodes[connec_order_in3[3]].id = conn_t2[connec_order_out3[3]];

				elem3->nodes[connec_order_in3[4]].id = original_conn[connec_order_out3[4]];
				elem3->nodes[connec_order_in3[5]].id = original_conn[connec_order_out3[5]];
				elem3->nodes[connec_order_in3[6]].id = original_conn[connec_order_out3[6]];
				elem3->nodes[connec_order_in3[7]].id = original_conn[connec_order_out3[7]];

				elem3->pad = elem->pad;
				elem3->level = elem->level;

				octant_t* elem4 = (octant_t*) sc_array_push(&mesh->elements);

				elem4->nodes[connec_order_in4[0]].id = original_conn[connec_order_out4[0]];

				elem4->nodes[connec_order_in4[1]].id = conn_t2[connec_order_out4[1]];
				elem4->nodes[connec_order_in4[2]].id = conn_t2[connec_order_out4[2]];
				elem4->nodes[connec_order_in4[3]].id = conn_t2[connec_order_out4[3]];
				elem4->nodes[connec_order_in4[4]].id = conn_t2[connec_order_out4[4]];
				elem4->nodes[connec_order_in4[5]].id = conn_t2[connec_order_out4[5]];
				elem4->nodes[connec_order_in4[6]].id = conn_t2[connec_order_out4[6]];
				elem4->nodes[connec_order_in4[7]].id = conn_t2[connec_order_out4[7]];

				elem4->pad = elem->pad;
				elem4->level = elem->level;
			}

		}

#endif

		for (int edge = 0; edge < 4; edge++) {
			if (point_v[edge]) gts_object_destroy(GTS_OBJECT(point_v[edge]));
			point_v[edge] = NULL;
		}

		for (int edge = 0; edge < 12; edge++) {
			if (point_s[edge]) gts_object_destroy(GTS_OBJECT(point_s[edge]));
			point_s[edge] = NULL;
		}

		for (int edge = 0; edge < 12; edge++) {
			if (point[edge]) gts_object_destroy(GTS_OBJECT(point[edge]));
			point[edge] = NULL;
		}

	}

	mesh->local_n_elements = mesh->elements.elem_count;
	//mesh->local_n_nodes = coords.size()/3;

	//MPI_Allreduce(&mesh->local_n_elements, &mesh->total_n_elements, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
	//MPI_Allreduce(&mesh->local_n_nodes, &mesh->total_n_nodes, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);

	//	if (mesh->mpi_rank == 0) {
	//		printf("Total number of elements: %lld\n", mesh->local_n_elements);
	//		printf("Total number of nodes: %lld\n", mesh->local_n_nodes);
	//	}

	fclose(fdbg);
	sc_hash_array_destroy(hash_nodes);
}
