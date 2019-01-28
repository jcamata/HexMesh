/* 
 * File:   main.cpp
 * Author: camata
 *
 * Created on March 19, 2015, 1:42 PM
 */

#include <cstdlib>
#include <vector>
#include <mpi.h>

#include <sc.h>
#include <sc_io.h>
#include <sc_containers.h>
#include "hexa.h"
#include "hilbert.h"
#include "refinement.h"


/*
 * 
 */

int main(int argc, char** argv) {

	hexa_tree_t mesh;

	std::vector<double> coords;
	std::vector<int> element_ids;
	std::vector<int> nodes_b_mat;

	int l = atoi(argv[1]);

	hexa_init(argc, argv, &mesh);

	hexa_tree_init(&mesh, l);
	hexa_tree_cube(&mesh);

	hexa_mesh(&mesh);

	const char * bathy = "./input/bathy_Pipo_small.gts";
	const char * topo =  "./input/topo_Pipo_small.gts";

	// Note that here we use a gts file.
	// There is a tool called stl2gts that convert STL files to GTS.
	// It is installed together with the gts library.
	GetMeshFromSurface(&mesh, topo, coords);

	GetInterceptedElements(&mesh, coords, element_ids, bathy);
	printf(" Elements intercepted: %lld\n\n", element_ids.size());

	printf(" Project nodes to the surface\n\n");
	MovingNodes(&mesh,coords, nodes_b_mat,bathy);

	printf(" Applying material \n\n");
	element_ids.clear();
	Apply_material(&mesh, coords, element_ids, bathy);

	printf(" Untangle meshes\n\n");
	UntagleMesh(&mesh, coords, nodes_b_mat);

	printf(" Optimization of the mesh\n\n");
	MeshOpt(&mesh,coords,nodes_b_mat);

	printf(" Extrude elements\n\n");
	ExtrudePMLElements(&mesh,coords);

	//clean vectors
	std::vector<int>().swap(element_ids);
	std::vector<int>().swap(nodes_b_mat);

	printf(" Writing output files \n\n");
	hexa_mesh_write_vtk(&mesh, "mesh", &coords);
	hexa_mesh_write_msh(&mesh, "mesh", &coords);
	hexa_mesh_write_h5(&mesh,"mesh", coords);

	hexa_mesh_write_vtk(&mesh, "test",NULL);

	printf(" Cleaning variables \n\n");

	hexa_tree_destroy(&mesh);
	hexa_finalize(&mesh);
	std::vector<double>().swap(coords);

	return 0;
}

