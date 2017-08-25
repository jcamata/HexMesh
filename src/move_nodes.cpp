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
#include <mpi.h>

#include "hexa.h"
#include "hilbert.h"
#include "refinement.h"


typedef struct {
  GPtrArray * array;
} ListOfPoints;

static void InsertPoint(GtsPoint * p, ListOfPoints * lp) {
    g_ptr_array_add(lp->array, p);
}

bool is_point_over_surface(GtsPoint * p, GNode * tree);



#if 0

void Move_nodes(hexa_tree_t* mesh, std::vector<double>& coords, std::vector<int>& element_ids) {
    bool clamped = true;
    GtsSegment * segments[12];
    GtsPoint * point[12];
    int Edge2GNode[12][2];
    int original_conn[8];
    FILE * fdbg;

    int id_offset = 0;

    fdbg = fopen("move_nodes.txt", "w");
    int temp_1 = 0;

    for (int iel = 0; iel < element_ids.size(); ++iel) {
        octant_t *elem = (octant_t*) sc_array_index(&mesh->elements, element_ids[iel]);
        elem->pad = -1;
        //printf("Element: %d\n", element_ids[iel]);

        bool edge_list[12];
        int ed_cont = 0;

        for (int edge = 0; edge < 12; ++edge) {
            point[edge] = NULL;
            int node1 = elem->nodes[EdgeVerticesMap[edge][0]].id;
            int node2 = elem->nodes[EdgeVerticesMap[edge][1]].id;

            //fprintf(fdbg, " Nodes : %d, %d\n", node1, node2);
            //printf("Element: %d, Edges: %d\n", element_ids[iel],edge);

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
                    fprintf(fdbg, "El: %d\n", element_ids[iel]);
                    fprintf(fdbg, "edge, %d\n ", edge);
                    fprintf(fdbg, "nodes  : %d, %d\n", node1, node2);
                    fprintf(fdbg, "point : %f, %f, %f\n", point[edge]->x, point[edge]->y, point[edge]->z);
                    edge_list[edge] = true;
                    ed_cont++;
                    break;
                }
                list = list->next;
            }
        }

        if (edge_list[4] && edge_list[5] && edge_list[6] && edge_list[7]) {
            temp_1++;
            for (int i_node = 4; i_node < 8; i_node++) {
                int node1 = elem->nodes[i_node].id;
                coords[node1 * 3 + 2] = point[i_node]->z;
            }
        }



        //
        //octant_node_t *gnode1 = (octant_node_t*) sc_array_index(&mesh->nodes, node1);
        //octant_node_t *gnode2 = (octant_node_t*) sc_array_index(&mesh->nodes, node2);
        /*
                if (edge_list[edge] && gnode1->color != 1 && gnode2->color != 1) {

                        //TODO SEGFAULT aqui no gts point
                        GtsPoint *p0 = gts_point_new(gts_point_class(),
                                        coords[node1 * 3],
                                        coords[node1 * 3 + 1],
                                        coords[node1 * 3 + 2]);

                        GtsPoint *p1 = gts_point_new(gts_point_class(),
                                        coords[node2 * 3],
                                        coords[node2 * 3 + 1],
                                        coords[node2 * 3 + 2]);


                        double d_c1 = gts_point_distance(point[edge], p0);
                        double d_c2 = gts_point_distance(point[edge], p1);


                        fprintf(fdbg, "El: %d\n", element_ids[iel]);
                        fprintf(fdbg, "nodes  : %d, %d\n", node1, node2);
                        fprintf(fdbg, "point : %f, %f, %f\n", point[edge]->x, point[edge]->y, point[edge]->z);
                        fprintf(fdbg, "no 0 : %f, %f, %f\n", coords[node1 * 3], coords[node1 * 3 + 1], coords[node1 * 3 + 2]);
                        fprintf(fdbg, "no 1 : %f, %f, %f\n", coords[node2 * 3], coords[node2 * 3 + 1], coords[node2 * 3 + 2]);
                        fprintf(fdbg, "d_c1: %f, d_c2: %f\n", d_c1, d_c2);

                        if(0){
                                if (edge == 4 || edge == 5 || edge == 6 || edge == 7) {

                                        if (coords[node2 * 3 + 2] == 0) {
                                                coords[node1 * 3 + 2] = point[edge]->z;
                                        } else if (coords[node1 * 3 + 2] == 0) {
                                                coords[node2 * 3 + 2] = point[edge]->z;
                                        } else {

                                                if (d_c1 >= d_c2) {
                                                        coords[node2 * 3 + 2] = point[edge]->z;
                                                        gnode2->pad = 1;

                                                } else {
                                                        coords[node1 * 3 + 2] = point[edge]->z;
                                                        gnode1->pad = 1;
                                                }

                                        }

                                } else {
                                        if (d_c1 >= d_c2) {
                                                coords[node2 * 3] = point[edge]->x;
                                                coords[node2 * 3 + 1] = point[edge]->y;
                                                gnode2->pad = 1;

                                        } else {
                                                coords[node1 * 3] = point[edge]->x;
                                                coords[node1 * 3 + 1] = point[edge]->y;
                                                gnode1->pad = 1;
                                        }
                                }

                        }
                }
         */

    }
    printf("templates1: %d\n", temp_1);
    fclose(fdbg);

}


// Change the node positions to fit the surface.

void Move_nodes(hexa_tree_t* mesh, const char* surface_bathy, vector<double>& coords, std::vector<int>& element_ids) {

    GtsPoint *p;
    double dz;
    double d, d0;
    double zmax;

    mesh->gdata.s = SurfaceRead(surface_bathy);

    FILE *fout = fopen("surface.dat", "w");
    gts_surface_print_stats(mesh->gdata.s, fout);
    fclose(fout);

    // Get the surface bounding box
    mesh->gdata.bbox = gts_bbox_surface(gts_bbox_class(), mesh->gdata.s);

    // Change the box size to cut the external elements
    double factor = 0.05;
    double x_factor = (mesh->gdata.bbox->x2 - mesh->gdata.bbox->x1) * factor;
    double y_factor = (mesh->gdata.bbox->y2 - mesh->gdata.bbox->y1) * factor;

    mesh->gdata.bbox->x1 += x_factor;
    mesh->gdata.bbox->y1 += y_factor;

    mesh->gdata.bbox->x2 -= x_factor;
    mesh->gdata.bbox->y2 -= y_factor;

    double Lx = (mesh->gdata.bbox->x2 - mesh->gdata.bbox->x1);
    double Ly = (mesh->gdata.bbox->y2 - mesh->gdata.bbox->y1);
    double zmin = ((Lx < Ly) ? -Lx : -Ly);

    double el_size = Lx / mesh->ncellx;

    // Build the bounding box tree
    mesh->gdata.bbt = gts_bb_tree_surface(mesh->gdata.s);

    p = gts_point_new(gts_point_class(), 0.0, 0.0, mesh->gdata.bbox->z2);

    for (int iel = 0; iel < element_ids.size(); ++iel) {
        octant_t *elem = (octant_t*) sc_array_index(&mesh->elements, element_ids[iel]);

        for (int ino = 0; ino < 8; ino++) {
            int node = elem->nodes[ino].id;

            if (elem->nodes[ino].z == 0) {
                coords[node * 3 + 2] = 0;
            } else {
                p->x = coords[3 * node];
                p->y = coords[3 * node + 1];

                //d = gts_bb_tree_point_distance(mesh->gdata.bbt, p, distance, NULL);

                p->z = 0;
                d0 = gts_bb_tree_point_distance(mesh->gdata.bbt, p, distance, NULL);

                int n_el = ceil(d0 / el_size / pow(3, elem->level - 1));

                dz = d0 / (double) n_el;

                int z_el = ceil(coords[3 * node + 2] / el_size / pow(3, elem->level - 1));

                z_el = abs(z_el) + 1;

                double z = -dz * z_el;

                //printf("El: %d, d0: %f, n_el: %d, z_el: %d, z: %f\n",element_ids[iel],d0,n_el,z_el,z);

                coords[node * 3 + 2] = z;

                /*
                zmax = mesh->gdata.bbox->z2 - d;
                //printf("El: %d, d: %f, zmax: %f, bbz2: %f\n",element_ids[iel],d,zmax,mesh->gdata.bbox->z2);

                dz = (zmax - zmin) / (double) mesh->ncellz;
                //double z = zmax +  dz* elem->nodes[ino].z;
                //printf("El: %d, node: %d, nz %d, z:%f\n",element_ids[iel],node,elem->nodes[ino].z,z);
                //coords[node * 3 + 2] = z;
                 *
                 */

            }
        }
    }

}

#endif

//Aplly the material properties to the elements

void FindNodesBetweenMaterials(hexa_tree_t *mesh, std::vector<double>& coords, std::vector<int>& nodes_b_mat, int8_t* flag_nodes) {

    sc_array_t *elements = &mesh->elements;

    bool over;
 
    memset(flag_nodes, 0, sizeof (int8_t) * mesh->local_n_nodes);

    for (int iel = 0; iel < elements->elem_count; ++iel) {
        octant_t *elem = (octant_t*) sc_array_index(&mesh->elements, iel);
        if (elem->n_mat == 0)
            for (int ino = 0; ino < 8; ++ino) {
                int id = elem->nodes[ino].id;
                if (flag_nodes[id] == 0) flag_nodes[id] = 1;

            }

    }

    for (int iel = 0; iel < elements->elem_count; ++iel) {
        octant_t *elem = (octant_t*) sc_array_index(&mesh->elements, iel);
        if (elem->n_mat == 1)
            for (int ino = 0; ino < 8; ++ino) {
                int id = elem->nodes[ino].id;
                if (flag_nodes[id] == 1) {
                    flag_nodes[id] = 2;
                    //nodes_b_mat.push_back(id);
                }
            }

    }


    for (int iel = 0; iel < elements->elem_count; ++iel) {
        octant_t *elem = (octant_t*) sc_array_index(&mesh->elements, iel);
        if (elem->n_mat == 1) {

            //Vertical edges:
            int edge11 = elem->nodes[0].id;
            int edge12 = elem->nodes[4].id;
            if ((flag_nodes[edge11] == 2) && (flag_nodes[edge12] == 2)) flag_nodes[edge11] = 0;

            int edge21 = elem->nodes[1].id;
            int edge22 = elem->nodes[5].id;
            if ((flag_nodes[edge21] == 2) && (flag_nodes[edge22] == 2)) flag_nodes[edge21] = 0;

            int edge31 = elem->nodes[2].id;
            int edge32 = elem->nodes[6].id;
            if ((flag_nodes[edge31] == 2) && (flag_nodes[edge32] == 2)) flag_nodes[edge31] = 0;

            int edge41 = elem->nodes[3].id;
            int edge42 = elem->nodes[7].id;
            if ((flag_nodes[edge41] == 2) && (flag_nodes[edge42] == 2)) flag_nodes[edge41] = 0;
        }
    }


    for (int i = 0; i < mesh->local_n_nodes; i++) {

        if (flag_nodes[i] == 2)
            nodes_b_mat.push_back(i);

    }


    //free(flag_nodes);

}

// Change the node positions to fit the surface.

void ProjectNodes(hexa_tree_t* mesh, vector<double>& coords, std::vector<int>& nodes_b_mat, int8_t* flag_nodes) {

    // Change the box size to cut the external elements
    double factor = 0.05;
    double x_factor = (mesh->gdata.bbox->x2 - mesh->gdata.bbox->x1) * factor;
    double y_factor = (mesh->gdata.bbox->y2 - mesh->gdata.bbox->y1) * factor;

    mesh->gdata.bbox->x1 += x_factor;
    mesh->gdata.bbox->y1 += y_factor;

    mesh->gdata.bbox->x2 -= x_factor;
    mesh->gdata.bbox->y2 -= y_factor;

    double Lx = (mesh->gdata.bbox->x2 - mesh->gdata.bbox->x1);
    double Ly = (mesh->gdata.bbox->y2 - mesh->gdata.bbox->y1);
    double zmin = ((Lx < Ly) ? -Lx : -Ly);

    ListOfPoints lst;
    lst.array = g_ptr_array_new();

    gts_surface_foreach_vertex(mesh->gdata.s, (GtsFunc) InsertPoint, &lst);

    GNode* kdtree = gts_kdtree_new(lst.array, NULL);

    g_ptr_array_free(lst.array, TRUE);


    double el_size = Lx / mesh->ncellx;

    GtsPoint *p = gts_point_new(gts_point_class(), 0.0, 0.0, 0.0);

    for (int i = 0; i < nodes_b_mat.size(); i++) {
        int node = nodes_b_mat[i];
        p->x = coords[3 * node + 0];
        p->y = coords[3 * node + 1];
        p->z = coords[3 * node + 2];


        if (p->z < 0.0) {

            GtsBBox * bbox;
            GSList * selected, * j;

            double epsilon = gts_bb_tree_point_distance(mesh->gdata.bbt, p, distance, NULL);
            /* build bounding box */
            bbox = gts_bbox_new(gts_bbox_class(),
                    p,
                    p->x - epsilon,
                    p->y - epsilon,
                    p->z - epsilon,
                    p->x + epsilon,
                    p->y + epsilon,
                    p->z + epsilon);

            /* select vertices which are inside bbox using kdtree */
            j = gts_kdtree_range(kdtree, bbox, NULL);
            //std::cout << "projecting node ... " << node << std::endl;

            double min_dist = epsilon;
            GtsPoint *pmin = NULL;
            while (j) {
                GtsPoint * v = (GtsPoint *) j->data;
                double d = gts_point_distance(v, p);
                if (d < min_dist)
                    pmin = v;
                j = j->next;
            }

            if (pmin != NULL) {
                p->x = pmin->x;
                p->z = pmin->y;
                p->z = pmin->z;
            } else {
                bool over = is_point_over_surface(p, mesh->gdata.bbt);
                if (over)
                    p->z -= epsilon;
                else
                    p->z += epsilon;
            }


            if (p->z >= 0.0) p->z = 0.0;

 

        } else
            p->z = 0;

        coords[3 * node + 0] = p->x;
        coords[3 * node + 1] = p->y;
        coords[3 * node + 2] = p->z;

    }

    for (int iel = 0; iel < mesh->elements.elem_count; ++iel) {
        octant_t *elem = (octant_t*) sc_array_index(&mesh->elements, iel);
        if (elem->n_mat == 1) {

            //Vertical edges:

            int edge11 = elem->nodes[0].id;
            int edge12 = elem->nodes[4].id;
            {
                int z1 = coords[3 * edge11 + 2];
                int z2 = coords[3 * edge12 + 2];
                if (z2 > z1) {
                    double tmp = z1;
                    coords[3 * edge11 + 2] = z2;
                    coords[3 * edge12 + 2] = tmp;
                    int i = flag_nodes[edge11];
                    flag_nodes[edge11] = flag_nodes[edge12];
                    flag_nodes[edge12] = i;

                }

            }


            int edge21 = elem->nodes[1].id;
            int edge22 = elem->nodes[5].id;
            {
                int z1 = coords[3 * edge21 + 2];
                int z2 = coords[3 * edge22 + 2];
                if (z2 > z1) {
                    double tmp = z1;
                    coords[3 * edge21 + 2] = z2;
                    coords[3 * edge22 + 2] = tmp;
                    int i = flag_nodes[edge21];
                    flag_nodes[edge21] = flag_nodes[edge22];
                    flag_nodes[edge22] = i;

                }

            }

            int edge31 = elem->nodes[2].id;
            int edge32 = elem->nodes[6].id;
            {
                int z1 = coords[3 * edge31 + 2];
                int z2 = coords[3 * edge32 + 2];
                if (z2 > z1) {
                    double tmp = z1;
                    coords[3 * edge31 + 2] = z2;
                    coords[3 * edge32 + 2] = tmp;
                    int i = flag_nodes[edge31];
                    flag_nodes[edge31] = flag_nodes[edge32];
                    flag_nodes[edge32] = i;
                }

            }


            int edge41 = elem->nodes[3].id;
            int edge42 = elem->nodes[7].id;
            {
                int z1 = coords[3 * edge41 + 2];
                int z2 = coords[3 * edge42 + 2];
                if (z2 > z1) {
                    double tmp = z1;
                    coords[3 * edge41 + 2] = z2;
                    coords[3 * edge42 + 2] = tmp;
                    int i = flag_nodes[edge41];
                    flag_nodes[edge41] = flag_nodes[edge42];
                    flag_nodes[edge42] = i;
                }

            }

        }
    }

    gts_kdtree_destroy(kdtree);

    nodes_b_mat.clear();

    for (int i = 0; i < mesh->local_n_nodes; i++) {
        mesh->part_nodes[i] = flag_nodes[i];
        if (flag_nodes[i] == 2)
            nodes_b_mat.push_back(i);
    }


}

void MovingNodes(hexa_tree_t* mesh, std::vector<double>& coords, std::vector<int>& nodes_b_mat, const char* surface) {

    int8_t* flag_nodes = (int8_t*) malloc(sizeof (int8_t) * mesh->local_n_nodes);

    GtsSurface *bathymetry     = SurfaceRead(surface);
    GNode      *bbt_bathymetry = gts_bb_tree_surface(bathymetry);


    FindNodesBetweenMaterials(mesh, coords, nodes_b_mat, flag_nodes);
    ProjectNodes(mesh, coords, nodes_b_mat, flag_nodes);

    for (int iel = 0; iel < mesh->elements.elem_count; ++iel) {
        octant_t *elem = (octant_t*) sc_array_index(&mesh->elements, iel);

        if (elem->n_mat == 1) {
            
            GtsPoint * point;
            //getting the baricenter of the upper surface;
            //find the centroid of the upper surface
            double cord_in_x[8], cord_in_y[8], cord_in_z[8];
            //add the nodes in the coord vector
            for (int ii = 0; ii < 8; ii++) {
                cord_in_x[ii] = coords[3 * elem->nodes[ii].id];
                cord_in_y[ii] = coords[3 * elem->nodes[ii].id + 1];
                cord_in_z[ii] = coords[3 * elem->nodes[ii].id + 2];
            }

            double cord_in_ref[3];
            cord_in_ref[0] = 0;
            cord_in_ref[1] = 0;
            cord_in_ref[2] = 0;
            point = LinearMapHex(cord_in_ref, cord_in_x, cord_in_y, cord_in_z);

            bool over = is_point_over_surface(point, bbt_bathymetry);

            if (!over) elem->n_mat = 0;
            
            gts_object_destroy(GTS_OBJECT(point));

        }
    }

    gts_bb_tree_destroy(bbt_bathymetry, TRUE);
    
    gts_object_destroy(GTS_OBJECT(bathymetry));
    
    free(flag_nodes);

}

