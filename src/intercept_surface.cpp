
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
#include "cgal_h.h"
#include "hilbert.h"

// 		Class used to save and process the intersections.
//		It inherits Boost's visitor class

class IntersectionPointsVisitor_3
		: public boost::static_visitor<void>
{
protected:
	// --- Protected members
	int nbOfIntersections;
public:
	// --- Public members
	std::vector<ExactPoint_3> intersection_vertices;

	// --- Constructors
	// Empty constructor.
	IntersectionPointsVisitor_3()
	{
		intersection_vertices.resize(1);
		nbOfIntersections = 0;
	};

	// --- Visitor operators
	// Intersection is a point
	void operator()(const ExactPoint_3& p)
	{
		intersection_vertices[nbOfIntersections] = p;
		++nbOfIntersections;
	};

	// Intersection is a segment
	void operator()(const ExactSegment_3& s)
	{
		intersection_vertices[nbOfIntersections] = s.source();
		++nbOfIntersections;
		intersection_vertices[nbOfIntersections] = s.target();
		++nbOfIntersections;
	};

	void clear()
	{
		nbOfIntersections = 0;
	};

	int size()
	{
		return nbOfIntersections;
	};

	// Intersection typedef (boost::variant)
	typedef CGAL::cpp11::result_of<ExactKernel::Intersect_3(ExactTriangle_3, ExactSegment_3)>::type
			Triangle_3_Intersection_Variant;
};

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

	mesh->tdata.s = SurfaceRead(surface_topo);

	FILE *fout = fopen("surface.dat", "w");
	gts_surface_print_stats(mesh->tdata.s, fout);
	fclose(fout);

	// Get the surface bounding box
	mesh->tdata.bbox = gts_bbox_surface(gts_bbox_class(), mesh->tdata.s);
	if (mesh->mpi_rank == 0) {
		printf("Bounding box: \n");
		printf(" x ranges from %f to %f\n", mesh->tdata.bbox->x1, mesh->tdata.bbox->x2);
		printf(" y ranges from %f to %f\n", mesh->tdata.bbox->y1, mesh->tdata.bbox->y2);
		printf(" z ranges from %f to %f\n", mesh->tdata.bbox->z1, mesh->tdata.bbox->z2);
	}

	// Change the box size to cut the external elements
	double factor = 0.05;
	double x_factor = (mesh->tdata.bbox->x2 - mesh->tdata.bbox->x1)*factor;
	double y_factor = (mesh->tdata.bbox->y2 - mesh->tdata.bbox->y1)*factor;

	mesh->tdata.bbox->x1 += x_factor;
	mesh->tdata.bbox->y1 += y_factor;

	mesh->tdata.bbox->x2 -= x_factor;
	mesh->tdata.bbox->y2 -= y_factor;

	double Lx = (mesh->tdata.bbox->x2 - mesh->tdata.bbox->x1);
	double Ly = (mesh->tdata.bbox->y2 - mesh->tdata.bbox->y1);
	double zmin = ((Lx < Ly) ? -Lx : -Ly);

	// Get grid-spacing at x and y direction
	dx = (mesh->tdata.bbox->x2 - mesh->tdata.bbox->x1) / (double) mesh->ncellx;
	dy = (mesh->tdata.bbox->y2 - mesh->tdata.bbox->y1) / (double) mesh->ncelly;

	coords.resize(nodes->elem_count * 3);

	// Build the bounding box tree
	mesh->tdata.bbt = gts_bb_tree_surface(mesh->tdata.s);

	p = gts_point_new(gts_point_class(), 0.0, 0.0, mesh->tdata.bbox->z2);

	for (int i = 0; i < nodes->elem_count; ++i) {
		octant_node_t* n = (octant_node_t*) sc_array_index(nodes, i);
		p->x = mesh->tdata.bbox->x1 + n->x*dx;
		p->y = mesh->tdata.bbox->y1 + n->y*dy;

		d = gts_bb_tree_point_distance(mesh->tdata.bbt, p, distance, NULL);
		zmax = mesh->tdata.bbox->z2 - d;

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
		//printf("Element:%d\n",elem->id);
		for (int edge = 0; edge < 12; ++edge) {
			point[edge] = NULL;
			elem->edge[edge].ref = false;
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
			//printf("edge:%d, %s ",edge, elem->edge[edge].ref ? "T" : "F");
		}
		//printf("\n");

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

GtsPoint* SegmentTriangleIntersectionCgal(GtsSegment * s, GtsTriangle * t){

	GtsPoint * A, * B, * C, * D, * E;
	GtsEdge * AB, * BC, * CA;
	GtsPoint * out ;

	//construction exact triangle
	gts_triangle_vertices_edges(t, NULL,
			(GtsVertex **) & A,
			(GtsVertex **) & B,
			(GtsVertex **) & C,
			&AB, &BC, &CA);

	//ExactTriangle_3 triangle;
	ExactTriangle_3 triangle = ExactTriangle_3(	ExactPoint_3( A->x, A->y, A->z),
			ExactPoint_3( B->x, B->y, B->z),
			ExactPoint_3( C->x, C->y, C->z));

	//construction exact segment
	D = GTS_POINT(s->v1);
	E = GTS_POINT(s->v2);

	//ExactSegment_3 segment;
	ExactSegment_3 segment = ExactSegment_3(ExactPoint_3(D->x,D->y,D->z), ExactPoint_3(E->x,E->y,E->z));

	//do the intersection
	//"Triangle_3_Intersection_Variant" is a pretty ugly typedef from CGAL,
	// used to save data tha can be either a bool, an ExactPoint_3 or an
	// ExactSegment_3, depending on the context. If you're using C++11
	// - and you should! - an simple "auto" does the same job.
	//Triangle_3_Intersection_Variant triangle_segment_intersect;

	// This "visitor" will contain the intersection data, and its operator()
	// behaves differently depending on whenever the input is an
	// ExactPoint_3 or an ExactSegment_3
	IntersectionPointsVisitor_3 intersection_data;

	Triangle_3_Intersection_Variant triangle_segment_intersect = CGAL::intersection(segment,triangle);
	std::vector<double> aux;

	bool DoIntersect = false;
	if( triangle_segment_intersect ) // Asked for a bool, got a bool
	{
		DoIntersect = true;
		boost::apply_visitor(intersection_data, *triangle_segment_intersect);
		const ExactPoint_3 p = ExactPoint_3(intersection_data.intersection_vertices[0].x(),
				intersection_data.intersection_vertices[0].y(),
				intersection_data.intersection_vertices[0].z());

		//const ExactPoint_3 p = boost::get<ExactPoint_3>(*triangle_segment_intersect);
		ExactKernel_to_Kernel toinexact;
		Point_3 pp = toinexact(p);
		//Point_3 pp = Point_3(0,0,0);
		//std::cout << "   " << p << std::endl;
		//std::cout << "   " << pp << std::endl;
		//std::cout << "   " << CGAL::to_double(pp.x()) << "   " << CGAL::to_double(pp.y())
		//<< "   " << CGAL::to_double(pp.z()) << std::endl;

		const double xx = CGAL::to_double(pp.x());
		const double yy = CGAL::to_double(pp.y());
		const double zz = CGAL::to_double(pp.z());
		//std::cout << "   " << xx << "   " << yy
		//<< "   " << zz << std::endl;
		//aux.push_back(xx);
		//aux.push_back(yy);
		//aux.push_back(zz);

		out = gts_point_new(gts_point_class(),xx,yy,zz);

		//printf("%f %f %f\n",out->x,out->y,out->z);
	}else{
		// No intersection was found, do nothing
	}

	if(DoIntersect){
		return out;
		//return gts_point_new(gts_point_class(),out->x,out->y,out->z);
	}else{
		return NULL;
	}

}
