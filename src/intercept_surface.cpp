
#include <gts.h>
#include <glib.h>
#include <vector>
using namespace std;
#include <sc_io.h>
#include <sc_containers.h>

#include "hexa.h"


int EdgeVerticesMap[12][2] = {
    {0,1}, // Edge 0
    {1,2}, // Edge 1
    {2,3}, //      2
    {3,0},
    {0,4},
    {1,5},
    {2,6},
    {3,7},
    {4,5},
    {5,6},
    {6,7},
    {7,4}
};

int FaceEdgesMap[6][4] = {
    {4,11,7,3 },
    {5,1 ,6,9 },
    {0,5,8,4  },
    {2,6,10,7 },
    {0,1,2,3  },
    {8,9,10,11}
};


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

 1.c: tres faces interceptadas: Veja template3 hexMesh/doc.


*/


GtsPoint* SegmentTriangleIntersection (GtsSegment * s, GtsTriangle * t)
{
  GtsPoint * A, * B, * C, * D, * E;
  gint ABCE, ABCD, ADCE, ABDE, BCDE;
  GtsEdge * AB, * BC, * CA;
  gdouble a, b, c;

  g_return_val_if_fail (s != NULL, NULL);
  g_return_val_if_fail (t != NULL, NULL);

  gts_triangle_vertices_edges (t, NULL, 
			       (GtsVertex **) &A, 
			       (GtsVertex **) &B, 
			       (GtsVertex **) &C, 
			       &AB, &BC, &CA);
  D = GTS_POINT (s->v1);
  E = GTS_POINT (s->v2);

  ABCE = gts_point_orientation_3d_sos (A, B, C, E);
  ABCD = gts_point_orientation_3d_sos (A, B, C, D);
  if (ABCE < 0 || ABCD > 0) {
    GtsPoint * tmpp;
    gint tmp;

    tmpp = E; E = D; D = tmpp;
    tmp = ABCE; ABCE = ABCD; ABCD = tmp;
  }
  if (ABCE < 0 || ABCD > 0)
    return NULL;
  ADCE = gts_point_orientation_3d_sos (A, D, C, E);
  if (ADCE < 0)
    return NULL;
  ABDE = gts_point_orientation_3d_sos (A, B, D, E);
  if (ABDE < 0)
    return NULL;
  BCDE = gts_point_orientation_3d_sos (B, C, D, E);
  if (BCDE < 0)
    return NULL;
  a = gts_point_orientation_3d (A, B, C, E);
  b = gts_point_orientation_3d (A, B, C, D);
  if (a != b) {
    c = a/(a - b);
    return gts_point_new (gts_point_class(),
			  E->x + c*(D->x - E->x),
			  E->y + c*(D->y - E->y),
			  E->z + c*(D->z - E->z));
  }
  /* D and E are contained within ABC */
#ifdef DEBUG
  fprintf (stderr, 
	   "segment: %p:%s triangle: %p:%s intersection\n"
	   "D and E contained in ABC\n",
	   s, GTS_NEDGE (s)->name, t, GTS_NFACE (t)->name);
#endif /* DEBUG */  
  g_assert (a == 0.0); 
  return gts_point_new (gts_point_class(),
			(E->x + D->x)/2.,
			(E->y + D->y)/2.,
			(E->z + D->z)/2.);
}


void CheckTemplate(hexa_tree_t* mesh, std::vector<double>& coords, std::vector<int>& elements_ids)
{
    GtsSegment* segments[12];
    bool        face_intecepted[6];
    GtsPoint*   point[12];
    
    for(int i = 0; i < elements_ids.size(); ++i)
    {
        octant_t *elem    = (octant_t*) sc_array_index(&mesh->elements, elements_ids[i]);
        for(int edge =0; edge < 12; ++edge)
        {
            int node1 = elem->nodes[EdgeVerticesMap[edge][0]].id;
            int node2 = elem->nodes[EdgeVerticesMap[edge][1]].id;
            
            GtsVertex *v1  = gts_vertex_new(gts_vertex_class(), coords[node1*3], coords[node1*3+1], coords[node1*3+2]);
            GtsVertex *v2  = gts_vertex_new(gts_vertex_class(), coords[node2*3], coords[node2*3+1], coords[node2*3+2]);
            segments[edge] = gts_segment_new(gts_segment_class(), v1, v2);
            GtsBBox   *sb  = gts_bbox_segment(gts_bbox_class(), segments[edge]);
            GSList* list   = gts_bb_tree_overlap(mesh->gdata.bbt, sb);
            if(list == NULL) continue;
            while(list)
            {
                GtsTriangle* tri = GTS_TRIANGLE(list->data);
                point[edge] = SegmentTriangleIntersection (segments[edge], tri);
                if(point[edge]) break;
                list = list->next;
            }
            
        }
        
        
        //check parallel faces
        for(int face=0; face < 6; ++face)
        {
            face_intecepted[face] = false;
            for(int fe=0; fe < 4; ++fe) {
               int edge = FaceEdgesMap[face][fe];
               if( point[edge] != NULL ) {
                   face_intecepted[face] = true;
                   break;
               }
            }
            if(face_intecepted[face]) continue;
        }
        
        int n_parallel_faces = (face_intecepted[0] && face_intecepted[1]) + 
                               (face_intecepted[2] && face_intecepted[3]) +
                               (face_intecepted[4] && face_intecepted[5]);
        
        if(n_parallel_faces == 4 ) 
        {
            // Apply template 1.
            elem->pad = elem->pad + 1;
            if((!face_intecepted[0]) && (!face_intecepted[1])) 
            {
                // New points are created splitting edges 0,2,8,10
                // Two new elements are created
                //  NP1 = Point splitting Edge1  
                GtsPoint *np1 = point[0];
                GtsPoint *np2 = point[2];
                GtsPoint *np3 = point[8];
                GtsPoint *np4 = point[10];
                
                g_assert(np1 != NULL);
                g_assert(np2 != NULL);
                g_assert(np3 != NULL);
                g_assert(np4 != NULL);
                // Remove the parent element form the list and
                // introduce two new children elements
                // Elem 1:
                // N0, NP1, NP2, N3, N4, NP3, NP4, N7
                // Elem 2:
                // NP1, N1, N2, NP2, NP3, N5, N6, NP4
                int conn[8];

                
                conn[0] = elem->nodes[0].id;
                conn[3] = elem->nodes[3].id;
                conn[4] = elem->nodes[4].id;
                conn[7] = elem->nodes[7].id;
                    
                coords.push_back(np1->x);
                coords.push_back(np1->y);
                coords.push_back(np1->z);      
                conn[1] = mesh->local_n_nodes;
                mesh->local_n_nodes++;
                
                coords.push_back(np2->x);
                coords.push_back(np2->y);
                coords.push_back(np2->z);      
                conn[2] = mesh->local_n_nodes;
                mesh->local_n_nodes++;
                
                coords.push_back(np3->x);
                coords.push_back(np3->y);
                coords.push_back(np3->z);      
                conn[5] = mesh->local_n_nodes;
                mesh->local_n_nodes++;
                
                coords.push_back(np4->x);
                coords.push_back(np4->y);
                coords.push_back(np4->z);      
                conn[6] = mesh->local_n_nodes;
                mesh->local_n_nodes++;
                
                
            } else if((!face_intecepted[2]) && (!face_intecepted[3]))
            {
                
                
            } else if((!face_intecepted[4]) && (!face_intecepted[5]))
            {
                
            } 
        }
        
    }
    
}
