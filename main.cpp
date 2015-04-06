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


/*
 * 
 */
int main(int argc, char** argv) {

    hexa_tree_t mesh;
    hexa_init(argc, argv, &mesh);
    
    hexa_tree_init(&mesh);
    hexa_tree_cube(&mesh,1);
    hexa_mesh(&mesh);
    hexa_mesh_write_vtk(&mesh,"mesh");
    hexa_mesh_write_unv(&mesh,"teste");
    hexa_tree_destroy(&mesh);
    hexa_finalize(&mesh);
    return 0;
}

