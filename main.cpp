/* 
 * File:   main.cpp
 * Author: camata
 *
 * Created on March 19, 2015, 1:42 PM
 */

#include <cstdlib>
#include <vector>

#include <sc.h>
#include "hexa.h"

void GetMeshFromSurface(hexa_tree_t* tree, const char* surface, std::vector<double>& coords);
void GetInterceptedElements(hexa_tree_t* tree, std::vector<double>& coords, std::vector<int> &element_ids);
void CheckTemplate(hexa_tree_t* mesh, std::vector<double>& coords, std::vector<int>& elements_ids);
/*
 * 
 */

int main(int argc, char** argv) {

    hexa_tree_t mesh;

    std::vector<double> coords;
    std::vector<int>  element_ids;
    int l = atoi(argv[1]);

    hexa_init(argc, argv, &mesh);

    hexa_tree_init(&mesh,l);
    hexa_tree_cube(&mesh);
    //hexa_debug_face_hanging(&mesh);
    
    hexa_mesh(&mesh);
    GetMeshFromSurface(&mesh,"simple.gts", coords);
    GetInterceptedElements(&mesh,coords,element_ids);
    CheckTemplate(&mesh, coords, element_ids);
    // Add PML elements

    printf(" Elements intercepted: %d\n", element_ids.size() );
    
    hexa_mesh_write_vtk(&mesh,"mesh", &coords);
    //hexa_mesh_write_unv(&mesh,"teste", &coords);
    //hexa_mesh_write_vtk(&mesh,"template", NULL);
    //hexa_mesh_write_unv(&mesh,"teste", NULL);
    hexa_tree_destroy(&mesh);
    hexa_finalize(&mesh);
    return 0;
}

void add_pml_elements(hexa_tree_t *mesh)
{
    
}
