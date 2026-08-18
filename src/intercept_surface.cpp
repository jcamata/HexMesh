#include <gts.h>
#include <glib.h>
#include <vector>
#include <iostream>
#include <unordered_map>
#include <cmath>
using namespace std;
#include <set>
#include <algorithm>
#include <sc.h>
#include <sc_io.h>
#include <sc_containers.h>

#include "hexa.h"
#include "cgal_h.h"
#include "hilbert.h"

#include <variant>

// CGAL::intersection returns boost::optional<boost::variant<...>> up to CGAL 5.x and
// std::optional<std::variant<...>> from CGAL 6. Both are visited below.
#if CGAL_VERSION_NR >= 1060000000
#define HEXA_VISIT(vis, var) std::visit((vis), (var))
#else
#define HEXA_VISIT(vis, var) boost::apply_visitor((vis), (var))
#endif

//         Class used to save and process the intersections.
//        Visitor for the variant returned by CGAL::intersection

class IntersectionPointsVisitor_3
{
protected:
    // --- Protected members
    int nbOfIntersections;

public:
    // boost::apply_visitor (CGAL 5.x) needs this; harmless for std::visit
    using result_type = void;

    // --- Public members
    std::vector<ExactPoint_3> intersection_vertices;

    // --- Constructors
    // Empty constructor.
    IntersectionPointsVisitor_3()
    {
        intersection_vertices.resize(2);
        nbOfIntersections = 0;
    };

    // --- Visitor operators
    // Intersection is a point
    void operator()(const ExactPoint_3 &p)
    {
        intersection_vertices[nbOfIntersections] = p;
        ++nbOfIntersections;
    };

    // Intersection is a segment
    void operator()(const ExactSegment_3 &s)
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
};

namespace
{

    struct VertexIndex
    {
        std::vector<GtsVertex *> verts;
        std::unordered_map<GtsVertex *, size_t> index;
    };

    static int CollectVertex(gpointer item, gpointer data)
    {
        VertexIndex *vtx = static_cast<VertexIndex *>(data);
        GtsVertex *v = static_cast<GtsVertex *>(item);
        if (vtx->index.emplace(v, vtx->verts.size()).second)
        {
            vtx->verts.push_back(v);
        }
        return FALSE;
    }

    struct EdgeAdjacency
    {
        VertexIndex *vtx;
        std::vector<std::vector<int>> *adj;
    };

    static int CollectEdgeAdjacency(gpointer item, gpointer data)
    {
        EdgeAdjacency *ctx = static_cast<EdgeAdjacency *>(data);
        GtsEdge *edge = static_cast<GtsEdge *>(item);
        GtsSegment *seg = GTS_SEGMENT(edge);
        auto it1 = ctx->vtx->index.find(seg->v1);
        auto it2 = ctx->vtx->index.find(seg->v2);
        if (it1 == ctx->vtx->index.end() || it2 == ctx->vtx->index.end())
            return FALSE;
        int i = static_cast<int>(it1->second);
        int j = static_cast<int>(it2->second);
        (*ctx->adj)[i].push_back(j);
        (*ctx->adj)[j].push_back(i);
        return FALSE;
    }

    struct EdgeLengthAccumulator
    {
        double sum = 0.0;
        uint64_t count = 0;
    };

    static int AccumulateEdgeLength(gpointer item, gpointer data)
    {
        EdgeLengthAccumulator *acc = static_cast<EdgeLengthAccumulator *>(data);
        GtsEdge *edge = static_cast<GtsEdge *>(item);
        GtsSegment *seg = GTS_SEGMENT(edge);
        GtsPoint *p1 = GTS_POINT(seg->v1);
        GtsPoint *p2 = GTS_POINT(seg->v2);
        double dx = p1->x - p2->x;
        double dy = p1->y - p2->y;
        double dz = p1->z - p2->z;
        acc->sum += std::sqrt(dx * dx + dy * dy + dz * dz);
        acc->count++;
        return FALSE;
    }

    static double AverageEdgeLength(GtsSurface *s)
    {
        EdgeLengthAccumulator acc;
        gts_surface_foreach_edge(s, AccumulateEdgeLength, &acc);
        if (acc.count == 0)
            return 0.0;
        return acc.sum / static_cast<double>(acc.count);
    }

    // Laplacian smoothing of the INPUT surfaces is off. It moves vertices in x, y and z, so
    // it rounds the vertical coastline wall out of vertical and low-passes the sea floor --
    // the opposite of what the surface is built for. Note the iteration count grows as the
    // surface gets finer (iters = ceil(target_h / mean_edge)), so a better bathymetry source
    // would be smoothed harder. Set to true to restore the old behaviour.
    static const bool kSmoothInputSurfaces = false;

    static void SmoothGtsSurfaceLaplacian(GtsSurface *s, double target_h, int max_iters = 50, double lambda = 0.5)
    {
        if (!s || target_h <= 0.0)
            return;

        VertexIndex vtx;
        gts_surface_foreach_vertex(s, CollectVertex, &vtx);
        const size_t n = vtx.verts.size();
        if (n == 0)
            return;

        std::vector<std::vector<int>> adj(n);
        EdgeAdjacency ctx{&vtx, &adj};
        gts_surface_foreach_edge(s, CollectEdgeAdjacency, &ctx);

        const double mean_edge = AverageEdgeLength(s);
        if (mean_edge <= 0.0)
            return;
        int iters = static_cast<int>(std::ceil(target_h / mean_edge));
        if (iters < 1)
            iters = 1;
        if (iters > max_iters)
            iters = max_iters;

        std::vector<double> nx(n), ny(n), nz(n);

        for (int iter = 0; iter < iters; ++iter)
        {
            for (size_t i = 0; i < n; ++i)
            {
                const auto &nei = adj[i];
                if (nei.empty())
                {
                    GtsPoint *p = GTS_POINT(vtx.verts[i]);
                    nx[i] = p->x;
                    ny[i] = p->y;
                    nz[i] = p->z;
                    continue;
                }
                double ax = 0.0, ay = 0.0, az = 0.0;
                for (int j : nei)
                {
                    GtsPoint *pj = GTS_POINT(vtx.verts[j]);
                    ax += pj->x;
                    ay += pj->y;
                    az += pj->z;
                }
                const double inv = 1.0 / static_cast<double>(nei.size());
                ax *= inv;
                ay *= inv;
                az *= inv;
                GtsPoint *p = GTS_POINT(vtx.verts[i]);
                nx[i] = p->x + lambda * (ax - p->x);
                ny[i] = p->y + lambda * (ay - p->y);
                nz[i] = p->z + lambda * (az - p->z);
            }

            for (size_t i = 0; i < n; ++i)
            {
                GtsPoint *p = GTS_POINT(vtx.verts[i]);
                p->x = nx[i];
                p->y = ny[i];
                p->z = nz[i];
            }
        }
    }

}

// Read the gts file format and create a gts surface.
GtsSurface *SurfaceRead(const char *fname)
{
    FILE *gts_file;
    GtsSurface *s;
    // GtsPoint *p;
    GtsFile *fp;

    gts_file = fopen(fname, "r");
    if (!gts_file)
        return NULL;
    fp = gts_file_new(gts_file);
    s = gts_surface_new(gts_surface_class(),
                        gts_face_class(),
                        gts_edge_class(),
                        gts_vertex_class());
    if (gts_surface_read(s, fp))
    {
        fputs("file on standard input is not a valid GTS file\n", stderr);
        fprintf(stderr, "stdin:%d:%d: %s\n", fp->line, fp->pos, fp->error);
        return NULL; /* failure */
    }

    gts_file_destroy(fp);
    return s;
}

// Compute the distance between a point and a triangle.
gdouble distance(GtsPoint *p, gpointer bounded)
{
    GtsTriangle *t = (GtsTriangle *)bounded;
    return gts_point_triangle_distance(p, t);
}

// Change the node positions to fit the surface.
void GetMeshFromSurface(hexa_tree_t *mesh, const char *surface_topo, vector<double> &coords)
{

    GtsPoint *p;
    double dx, dy, dz;
    double d;
    double zmax;
    sc_array_t *nodes = &mesh->nodes;

    mesh->tdata.s = SurfaceRead(surface_topo);

    if (mesh->gdata.s)
    {
        GtsBBox *tmp_bbox = gts_bbox_surface(gts_bbox_class(), mesh->gdata.s);
        if (tmp_bbox)
        {
            double hx = (tmp_bbox->x2 - tmp_bbox->x1) / static_cast<double>(mesh->ncellx);
            double hy = (tmp_bbox->y2 - tmp_bbox->y1) / static_cast<double>(mesh->ncelly);
            double hz = mesh->input.z / static_cast<double>(mesh->ncellz);
            double h = 2.1 * std::min(hx, std::min(hy, hz));
            if (kSmoothInputSurfaces) {
                printf("Smoothing bathymetry surface with target edge length hx: %f hy: %f hz: %f h: %f\n", hx, hy, hz, h);
                SmoothGtsSurfaceLaplacian(mesh->gdata.s, h);
            }
        }
        mesh->gdata.bbox = gts_bbox_surface(gts_bbox_class(), mesh->gdata.s);
    }

    FILE *fout = fopen("surfaceOut.dat", "w");
    gts_surface_print_stats(mesh->tdata.s, fout);
    fclose(fout);

    // Get the surface bounding box
    mesh->tdata.bbox = gts_bbox_surface(gts_bbox_class(), mesh->tdata.s);
    if (mesh->mpi_rank == 0)
    {
        printf("Bounding box: \n");
        printf(" x ranges from %f to %f\n", mesh->tdata.bbox->x1, mesh->tdata.bbox->x2);
        printf(" y ranges from %f to %f\n", mesh->tdata.bbox->y1, mesh->tdata.bbox->y2);
        printf(" z ranges from %f to %f\n", mesh->tdata.bbox->z1, mesh->tdata.bbox->z2);
    }

    // Change the box size to cut the external elements
    double factor = 0.02;
    double x_factor = (mesh->tdata.bbox->x2 - mesh->tdata.bbox->x1) * factor;
    double y_factor = (mesh->tdata.bbox->y2 - mesh->tdata.bbox->y1) * factor;

    mesh->tdata.bbox->x1 += x_factor;
    mesh->tdata.bbox->y1 += y_factor;

    mesh->tdata.bbox->x2 -= x_factor;
    mesh->tdata.bbox->y2 -= y_factor;

    double Lx = (mesh->tdata.bbox->x2 - mesh->tdata.bbox->x1);
    double Ly = (mesh->tdata.bbox->y2 - mesh->tdata.bbox->y1);
    // double zmin = ((Lx < Ly) ? -Lx : -Ly);
    double zmin = -mesh->input.z;

    // Get grid-spacing at x and y direction
    dx = (mesh->tdata.bbox->x2 - mesh->tdata.bbox->x1) / (double)mesh->ncellx;
    dy = (mesh->tdata.bbox->y2 - mesh->tdata.bbox->y1) / (double)mesh->ncelly;

    coords.resize(nodes->elem_count * 3);

    // Build the bounding box tree
    mesh->tdata.bbt = gts_bb_tree_surface(mesh->tdata.s);

    p = gts_point_new(gts_point_class(), 0.0, 0.0, mesh->tdata.bbox->z2);

    for (int i = 0; i < nodes->elem_count; ++i)
    {
        octant_node_t *n = (octant_node_t *)sc_array_index(nodes, i);
        p->x = mesh->tdata.bbox->x1 + n->x * dx;
        p->y = mesh->tdata.bbox->y1 + n->y * dy;

        d = gts_bb_tree_point_distance(mesh->tdata.bbt, p, distance, NULL);
        zmax = mesh->tdata.bbox->z2 - d;

        dz = (zmax - zmin) / (double)mesh->ncellz;
        double z = zmax - (n->z) * dz;

        coords[i * 3 + 0] = p->x;
        coords[i * 3 + 1] = p->y;
        coords[i * 3 + 2] = z;
    }
}

// Found the intercepted elements
void GetInterceptedElements(hexa_tree_t *mesh, std::vector<double> &coords, std::vector<int> &elements_ids, const char *surface_bathy)
{
    sc_array_t *elements = &mesh->elements;
    GtsBBox *box;

    if (mesh->input.inter_files.empty() && surface_bathy && strlen(surface_bathy) > 0) {
        mesh->input.inter_files.push_back(std::string(surface_bathy));
    }
    if (mesh->input.inter_files.empty() && !mesh->input.inter.empty()) {
        mesh->input.inter_files.push_back(mesh->input.inter);
    }

    mesh->gdata_vec.resize(mesh->input.inter_files.size());

    for (size_t k = 0; k < mesh->input.inter_files.size(); k++) {
        mesh->gdata_vec[k].s = SurfaceRead(mesh->input.inter_files[k].c_str());
        if (mesh->gdata_vec[k].s)
        {
            GtsBBox *tmp_bbox = gts_bbox_surface(gts_bbox_class(), mesh->gdata_vec[k].s);
            if (tmp_bbox)
            {
                double hx = (tmp_bbox->x2 - tmp_bbox->x1) / static_cast<double>(mesh->ncellx);
                double hy = (tmp_bbox->y2 - tmp_bbox->y1) / static_cast<double>(mesh->ncelly);
                double hz = mesh->input.z / static_cast<double>(mesh->ncellz);
                double h = 2.1 * std::min(hx, std::min(hy, hz));
                if (kSmoothInputSurfaces) {
                    printf("Smoothing interface %zu surface with target edge length hx: %f hy: %f hz: %f h: %f\n", k+1, hx, hy, hz, h);
                    SmoothGtsSurfaceLaplacian(mesh->gdata_vec[k].s, h);
                }
            }
            mesh->gdata_vec[k].bbox = gts_bbox_surface(gts_bbox_class(), mesh->gdata_vec[k].s);
            mesh->gdata_vec[k].bbt = gts_bb_tree_surface(mesh->gdata_vec[k].s);
        }
    }

    if (!mesh->gdata_vec.empty()) {
        mesh->gdata = mesh->gdata_vec[0];
    }

    box = gts_bbox_new(gts_bbox_class(), 0, 0, 0, 0, 1, 1, 1);

    for (int iel = 0; iel < elements->elem_count; ++iel)
    {

        octant_t *elem = (octant_t *)sc_array_index(&mesh->elements, iel);
        elem->pad = 0;

        box->x1 = box->y1 = box->z1 = 1.0E10;
        box->x2 = box->y2 = box->z2 = -1.0E10;

        for (int i = 0; i < 8; ++i)
        {
            octant_node_t *node = &elem->nodes[i];
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

        bool any_overlap = false;
        for (size_t k = 0; k < mesh->gdata_vec.size(); k++) {
            if (mesh->gdata_vec[k].bbt && gts_bb_tree_is_overlapping(mesh->gdata_vec[k].bbt, box)) {
                any_overlap = true;
                break;
            }
        }

        if (any_overlap)
        {
            elements_ids.push_back(iel);
            elem->pad = -1;
        }

        ///////////
        GtsSegment *segments[12] = {0};
        GtsPoint *point[12] = {NULL};
        int ed_cont = 0;
        // printf("Element:%d\n",elem->id);
        for (int edge = 0; edge < 12; ++edge)
        {
            point[edge] = NULL;
            elem->edge[edge].ref = false;
            int node1 = elem->nodes[EdgeVerticesMap[edge][0]].id;
            int node2 = elem->nodes[EdgeVerticesMap[edge][1]].id;
            double x1 = coords[node1 * 3], y1 = coords[node1 * 3 + 1], z1 = coords[node1 * 3 + 2];
            double x2 = coords[node2 * 3], y2 = coords[node2 * 3 + 1], z2 = coords[node2 * 3 + 2];
            // Extend segment 2% beyond both endpoints to catch surface intersections that
            // land exactly at the element boundary (endpoint-at-surface floating-point miss).
            const double ext = 0.02;
            double dx = x2 - x1, dy = y2 - y1, dz = z2 - z1;
            GtsVertex *v1 = gts_vertex_new(gts_vertex_class(), x1 - ext * dx, y1 - ext * dy, z1 - ext * dz);
            GtsVertex *v2 = gts_vertex_new(gts_vertex_class(), x2 + ext * dx, y2 + ext * dy, z2 + ext * dz);
            segments[edge] = gts_segment_new(gts_segment_class(), v1, v2);
            GtsBBox *sb = gts_bbox_segment(gts_bbox_class(), segments[edge]);

            for (size_t k = 0; k < mesh->gdata_vec.size(); k++) {
                if (!mesh->gdata_vec[k].bbt) continue;
                GSList *list = gts_bb_tree_overlap(mesh->gdata_vec[k].bbt, sb);
                if (list == NULL) continue;
                while (list)
                {
                    GtsBBox *b = GTS_BBOX(list->data);
                    if (mesh->input.CgalUse)
                    {
                        point[edge] = SegmentTriangleIntersectionCgal(segments[edge], GTS_TRIANGLE(b->bounded));
                    }
                    else
                    {
                        point[edge] = SegmentTriangleIntersection(segments[edge], GTS_TRIANGLE(b->bounded));
                    }
                    if (point[edge])
                    {
                        elem->edge[edge].ref = true;
                        elem->pad = -1;
                        ed_cont++;
                        break;
                    }
                    list = list->next;
                }
                if (point[edge]) break;
            }
            // printf("edge:%d, %s ",edge, elem->edge[edge].ref ? "T" : "F");
        }
        // printf("\n");

        // Bounding box intercepted
        if (elem->pad == -1 && ed_cont == 0)
        {
            elem->pad = 0;
            for (int edge = 0; edge < 12; ++edge)
            {
                elem->edge[edge].ref = false;
            }
        }
    }
}

// Found the intersection between a line and a triangle
GtsPoint *SegmentTriangleIntersection(GtsSegment *s, GtsTriangle *t)
{
    GtsPoint *A, *B, *C, *D, *E;
    gint ABCE, ABCD, ADCE, ABDE, BCDE;
    GtsEdge *AB, *BC, *CA;
    gdouble a, b, c;

    // g_return_val_if_fail(s != NULL, NULL);
    // g_return_val_if_fail(t != NULL, NULL);
    if (s == NULL)
        return NULL;
    if (t == NULL)
        return NULL;

    gts_triangle_vertices_edges(t, NULL,
                                (GtsVertex **)&A,
                                (GtsVertex **)&B,
                                (GtsVertex **)&C,
                                &AB, &BC, &CA);
    D = GTS_POINT(s->v1);
    E = GTS_POINT(s->v2);

    ABCE = gts_point_orientation_3d_sos(A, B, C, E);
    ABCD = gts_point_orientation_3d_sos(A, B, C, D);
    if (ABCE < 0 || ABCD > 0)
    {
        GtsPoint *tmpp;
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
    if (a != b)
    {
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
    // g_assert(a == 0.0);
    return gts_point_new(gts_point_class(),
                         (E->x + D->x) / 2.,
                         (E->y + D->y) / 2.,
                         (E->z + D->z) / 2.);
}

// Found the intersection between a line and a triangle
GtsPoint *SegmentTriangleIntersectionCgal(GtsSegment *s, GtsTriangle *t)
{

    GtsPoint *A, *B, *C, *D, *E;
    GtsEdge *AB, *BC, *CA;
    GtsPoint *out;

    // construction exact triangle
    gts_triangle_vertices_edges(t, NULL,
                                (GtsVertex **)&A,
                                (GtsVertex **)&B,
                                (GtsVertex **)&C,
                                &AB, &BC, &CA);

    // ExactTriangle_3 triangle;
    ExactTriangle_3 triangle = ExactTriangle_3(ExactPoint_3(A->x, A->y, A->z),
                                               ExactPoint_3(B->x, B->y, B->z),
                                               ExactPoint_3(C->x, C->y, C->z));

    // construction exact segment
    D = GTS_POINT(s->v1);
    E = GTS_POINT(s->v2);

    // ExactSegment_3 segment;
    ExactSegment_3 segment = ExactSegment_3(ExactPoint_3(D->x, D->y, D->z), ExactPoint_3(E->x, E->y, E->z));

    // do the intersection
    //"Triangle_3_Intersection_Variant" is a pretty ugly typedef from CGAL,
    //  used to save data tha can be either a bool, an ExactPoint_3 or an
    //  ExactSegment_3, depending on the context. If you're using C++11
    //  - and you should! - an simple "auto" does the same job.
    // Triangle_3_Intersection_Variant triangle_segment_intersect;

    // This "visitor" will contain the intersection data, and its operator()
    // behaves differently depending on whenever the input is an
    // ExactPoint_3 or an ExactSegment_3
    IntersectionPointsVisitor_3 intersection_data;

    auto triangle_segment_intersect = CGAL::intersection(segment, triangle);
    std::vector<double> aux;

    bool DoIntersect = false;
    if (triangle_segment_intersect) // Asked for a bool, got a bool
    {
        DoIntersect = true;
        HEXA_VISIT(intersection_data, *triangle_segment_intersect);
        const ExactPoint_3 p = ExactPoint_3(intersection_data.intersection_vertices[0].x(),
                                            intersection_data.intersection_vertices[0].y(),
                                            intersection_data.intersection_vertices[0].z());

        // const ExactPoint_3 p = boost::get<ExactPoint_3>(*triangle_segment_intersect);
        ExactKernel_to_Kernel toinexact;
        Point_3 pp = toinexact(p);
        // Point_3 pp = Point_3(0,0,0);
        // std::cout << "   " << p << std::endl;
        // std::cout << "   " << pp << std::endl;
        // std::cout << "   " << CGAL::to_double(pp.x()) << "   " << CGAL::to_double(pp.y())
        //<< "   " << CGAL::to_double(pp.z()) << std::endl;

        const double xx = CGAL::to_double(pp.x());
        const double yy = CGAL::to_double(pp.y());
        const double zz = CGAL::to_double(pp.z());
        // std::cout << "   " << xx << "   " << yy
        //<< "   " << zz << std::endl;
        // aux.push_back(xx);
        // aux.push_back(yy);
        // aux.push_back(zz);

        out = gts_point_new(gts_point_class(), xx, yy, zz);

        // printf("%f %f %f\n",out->x,out->y,out->z);
    }
    else
    {
        // No intersection was found, do nothing
    }

    if (DoIntersect)
    {
        return out;
        // return gts_point_new(gts_point_class(),out->x,out->y,out->z);
    }
    else
    {
        return NULL;
    }
}
