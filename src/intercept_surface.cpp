
#include <gts.h>
#include <glib.h>
#include <vector>
using namespace std;

#include "hexa.h"


// Read the gts file format and create a gts surface.
GtsSurface* SurfaceRead(const char* fname)
{
    FILE      *gts_file;
    GtsSurface *s;
    GtsPoint   *p;
    GtsFile    *fp;
    
    gts_file = fopen(fname, "r");
    if(!gts_file) return NULL;
    fp = gts_file_new(gts_file);
    s  = gts_surface_new(gts_surface_class (),
		         gts_face_class (),
		         gts_edge_class (),
		         gts_vertex_class ());
    if (gts_surface_read (s, fp)) {
        fputs ("file on standard input is not a valid GTS file\n", stderr);
        fprintf (stderr, "stdin:%d:%d: %s\n", fp->line, fp->pos, fp->error);
        return NULL; /* failure */
    }
    
    gts_file_destroy(fp);
    return s;
    
}

// Compute the distance between a point and a triangle.
gdouble  distance(GtsPoint *p, gpointer bounded)
{
    GtsTriangle *t = (GtsTriangle*) bounded;
    return gts_point_triangle_distance(p, t);
}

/*
void GetElementNodes(octant_t* n, octant_node_t* noel[])
{
    for(int i =0; i < 8; i++)
    {
        noel[i] = (octant_node_t*) sc_array_index(&n->nodes, i);
    }
  
}
*/

void GetMeshFromSurface(hexa_tree_t* tree, const char* surface, vector<double>& coords)
{

    GtsSurface *s;
    GtsBBox    *box;
    GNode      *t;
    GtsPoint   *p;
    int    nx,ny,nz;
    double dx, dy, dz;
    sc_array_t *nodes    = &tree->nodes;
    sc_array_t *elements = &tree->elements;
    
    
    // Note that here we use a gts file.
    // There is a tool called stl2gts that convert STL files to GTS.
    // It is installed together with the gts library.
    s   = SurfaceRead(surface);
    
    // Get the surface bounding box
    box = gts_bbox_surface (gts_bbox_class (), s);
    if(tree->mpi_rank == 0) {
        printf("Bounding box: \n");
        printf(" x ranges from %f to %f\n", box->x1, box->x2);
        printf(" y ranges from %f to %f\n", box->y1, box->y2);
        printf(" z ranges from %f to %f\n", box->z1, box->z2);
    }
    
    double Lx = (box->x2 - box->x1);
    double Ly = (box->y2 - box->y1);
    double zmin = ((Lx < Ly)?-Lx:-Ly);
    // Get grid-spacing at x and y direction
    dx = (box->x2 - box->x1)/(double)tree->ncellx;
    dy = (box->y2 - box->y1)/(double)tree->ncelly;

    coords.resize(nodes->elem_count*3);
    
    // Build the bounding box tree
    t = gts_bb_tree_surface(s);
    p = gts_point_new(gts_point_class(), 0.0 , 0.0 , box->z2);
    
    for(int i=0; i < nodes->elem_count; ++i)
    {
        octant_node_t* n = (octant_node_t*) sc_array_index(nodes,i);
        p->x = box->x1 + n->x*dx;
        p->y = box->y1 + n->y*dy;
        double d    = gts_bb_tree_point_distance(t,p,distance,NULL);
        double zmax = box->z2 - d;
        dz = (zmax-zmin)/(double)tree->ncellz;
        double z = zmax - n->z*dz;
        /*
        if(z < 0.0) {
            zmax = 0.0;
            dz = (zmax-zmin)/(double)tree->ncellz;
            z  = zmax - n->z*dz;
         }
         */
        
        coords[i*3 + 0] = p->x;
        coords[i*3 + 1] = p->y;
        coords[i*3 + 2] = z;
    }
    
    int fluid_elem = 0;
    int solid_elem = 0;
    for(int i=0; i < elements->elem_count; ++i)
    {
        
        octant_t *h     = (octant_t*) sc_array_index(elements, i);
        octant_node_t* noels = h->nodes;
        h->pad = 0;
        // Calculando octant baricenter
        
        int I   = h->nodes[0].id;
        int IJ  = h->nodes[2].id;
        int IJK = h->nodes[5].id;
        double dx =  0.5*(coords[IJ*3  ] - coords[I*3]);
        double dy =  0.5*(coords[IJ*3+1] - coords[I*3+1]);
        double dz =  0.5*(coords[IJK*3+2]- coords[I*3+2]);
        double bx = coords[i*3]   + dx;
        double by = coords[i*3+1] + dy;
        double bz = coords[i*3+2] + dz;
        if(bz < 0.0) {
            p->x = bx;
            p->y = by;
            double d    = gts_bb_tree_point_distance(t,p,distance,NULL);
            if(bz+d>0.0) { 
                h->pad = 1;
                fluid_elem++;
            } else
            {
                solid_elem++;
            }
        } else
            solid_elem++;
    }
    
    printf("Number of fluid elements: %d\n", fluid_elem);
    printf("Number of solid elements: %d\n", solid_elem);
    
}