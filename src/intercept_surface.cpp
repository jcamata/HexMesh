
#include <gts.h>
#include <glib.h>
#include <vector>
using namespace std;
#include <sc_io.h>
#include <sc_containers.h>

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


void GetMeshFromSurface(hexa_tree_t* tree, const char* surface, vector<double>& coords)
{

    GtsSurface *s;
    GtsBBox    *box;
    GNode      *t  ;
    GtsPoint   *p;
    int    nx,ny,nz;
    double dx, dy, dz;
    sc_array_t *nodes    = &tree->nodes;
    sc_array_t *elements = &tree->elements;
    
    
    // Note that here we use a gts file.
    // There is a tool called stl2gts that convert STL files to GTS.
    // It is installed together with the gts library.
    tree->gdata.s   = SurfaceRead(surface);
    
    // Get the surface bounding box
    tree->gdata.bbox = gts_bbox_surface (gts_bbox_class (), tree->gdata.s);
    if(tree->mpi_rank == 0) {
        printf("Bounding box: \n");
        printf(" x ranges from %f to %f\n", tree->gdata.bbox->x1, tree->gdata.bbox->x2);
        printf(" y ranges from %f to %f\n", tree->gdata.bbox->y1, tree->gdata.bbox->y2);
        printf(" z ranges from %f to %f\n", tree->gdata.bbox->z1, tree->gdata.bbox->z2);
    }
    
    double Lx = (tree->gdata.bbox->x2 - tree->gdata.bbox->x1);
    double Ly = (tree->gdata.bbox->y2 - tree->gdata.bbox->y1);
    double zmin = ((Lx < Ly)?-Lx:-Ly);

    // Get grid-spacing at x and y direction
    dx = (tree->gdata.bbox->x2 - tree->gdata.bbox->x1)/(double)tree->ncellx;
    dy = (tree->gdata.bbox->y2 - tree->gdata.bbox->y1)/(double)tree->ncelly;

    coords.resize(nodes->elem_count*3);
    
    // Build the bounding box tree
    tree->gdata.bbt = gts_bb_tree_surface(tree->gdata.s);
    p = gts_point_new(gts_point_class(), 0.0 , 0.0 , tree->gdata.bbox->z2);
    
    //GSList *list = gts_bb_tree_overlap(t, box);

    for(int i=0; i < nodes->elem_count; ++i)
    {
        octant_node_t* n = (octant_node_t*) sc_array_index(nodes,i);
        p->x = tree->gdata.bbox->x1 + n->x*dx;
        p->y = tree->gdata.bbox->y1 + n->y*dy;
        double d    = gts_bb_tree_point_distance(tree->gdata.bbt,p,distance,NULL);
        double zmax = tree->gdata.bbox->z2 - d;
        
        // implict water level (z=0).

        if(zmax < 0.0) zmax = 0.0;
        dz = (zmax-zmin)/(double)tree->ncellz;
        double z = zmax - (n->z)*dz;
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
    
    
}



void GetInterceptedElements(hexa_tree_t* mesh, std::vector<double>& coords, std::vector<int>& elements_ids)
{
    sc_array_t *elements = &mesh->elements;
    GtsBBox    *box;

    box = gts_bbox_new (gts_bbox_class (),0,0,0,0,1,1,1);

    for(int iel; iel < elements->elem_count; ++iel)
    {

        octant_t *elem    = (octant_t*) sc_array_index(&mesh->elements, iel);

        box->x1 = box->y1= box->z1  = 1.0E10;
        box->x2 = box->y2 = box->z2 = -1.0E10;
        for(int i = 0; i < 8; ++i)
        {
            octant_node_t* node = &elem->nodes[i];
            int id = node->id;
            double x = coords[id*3];
            double y = coords[id*3+1];
            double z = coords[id*3+2];
            box->x1 = (x < box->x1)?x:box->x1;
            box->y1 = (y < box->y1)?y:box->y1;
            box->z1 = (z < box->z1)?z:box->z1;
            box->x2 = (x > box->x2)?x:box->x2;
            box->y2 = (y > box->y2)?y:box->y2;
            box->z2 = (z > box->z2)?z:box->z2;
        }    
        elem->pad = 0;

        if(gts_bb_tree_is_overlapping(mesh->gdata.bbt, box)) 
        {
            elements_ids.push_back(iel);
            elem->pad = 1;
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

 1.c: tres faces interceptadas: Vega template3 hexMesh/doc.
















*/