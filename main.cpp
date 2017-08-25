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

void GetMeshFromSurface(hexa_tree_t* tree, const char* surface_topo, std::vector<double>& coords);
void GetInterceptedElements(hexa_tree_t* mesh, std::vector<double>& coords, std::vector<int>& elements_ids, const char* surface_bathy);
void CheckOctreeTemplate(hexa_tree_t* mesh, const std::vector<double>& coords, std::vector<int>& elements_ids, bool flag);
void ApplyOctreeTemplate(hexa_tree_t* mesh, std::vector<double>& coords, std::vector<int>& elements_ids);

//void ChangeTemplate(hexa_tree_t* mesh, std::vector<double>& coords, std::vector<int>& elements_ids);
//void ApplyTemplate(hexa_tree_t* mesh, std::vector<double>& coords, std::vector<int>& elements_ids);

void Apply_material(hexa_tree_t* mesh, std::vector<double>& coords,std::vector<int>& element_ids, const char* surface_bathy);
void Move_nodes(hexa_tree_t* tree, const char* surface_bathy, std::vector<double>& coords, std::vector<int>& element_ids);

void AddPMLElements(hexa_tree_t* mesh);
void ExtrudePMLElements(hexa_tree_t* mesh, std::vector<double>& coords);

void IdentifyTemplate(hexa_tree_t* mesh, const std::vector<double>& coords, std::vector<int>& elements_ids);

void MovingNodes(hexa_tree_t* mesh, std::vector<double>& coords, std::vector<int>& nodes_b_mat, const char* surface);

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

	//hexa_debug_face_hanging(&mesh);
	//AddPMLElements(&mesh);

	hexa_mesh(&mesh);

	if(0){
		// Note that here we use a gts file.
		// There is a tool called stl2gts that convert STL files to GTS.
		// It is installed together with the gts library.
		GetMeshFromSurface(&mesh, "./input/topo_Arg_small.gts", coords);
		GetInterceptedElements(&mesh, coords, element_ids, "./input/bathy_Arg_small.gts");

		printf(" Elements intercepted: %lld\n\n", element_ids.size());

		//printf(" Check and propagate 27-tree templates\n\n");
		CheckOctreeTemplate(&mesh, coords, element_ids, true);

		//printf(" Apply 27-tree templates\n\n");
		ApplyOctreeTemplate(&mesh, coords, element_ids);

		printf(" Applying material \n\n");
		element_ids.clear();
		Apply_material(&mesh, coords, element_ids, "./input/bathy_Arg_small.gts");

		//printf(" Project nodes to the bathymetry\n\n");
		//Move_nodes(&mesh,"./input/bathy_Arg_small.gts", coords,element_ids);
	}

	if(1){
		// Note that here we use a gts file.
		// There is a tool called stl2gts that convert STL files to GTS.
		// It is installed together with the gts library.
		GetMeshFromSurface(&mesh, "./input/topo_Pipo_small.gts", coords);
		GetInterceptedElements(&mesh, coords, element_ids, "./input/bathy_Pipo_small.gts");

		printf(" Elements intercepted: %lld\n\n", element_ids.size());

		printf(" Check and propagate 27-tree templates\n\n");
		CheckOctreeTemplate(&mesh, coords, element_ids, true);

		//printf(" Apply 27-tree templates\n\n");
		//ApplyOctreeTemplate(&mesh, coords, element_ids);

		//printf(" Applying material \n\n");
		element_ids.clear();
		Apply_material(&mesh, coords, element_ids, "./input/bathy_Pipo_small.gts");

//		printf(" Project nodes to the bathymetry\n\n");
	//	Move_nodes(&mesh,"./input/bathy_Pipo_small.gts", coords,element_ids);

                MovingNodes(&mesh,coords, nodes_b_mat,"./input/bathy_Pipo_small.gts");
                        

	}

	printf(" Writing output files \n\n");
	hexa_mesh_write_vtk(&mesh, "mesh", &coords);
        
	//hexa_mesh_write_msh(&mesh, "mesh", &coords);
	//hexa_mesh_write_h5(&mesh,"mesh", coords);

	printf(" Cleaning variables \n\n");

	hexa_tree_destroy(&mesh);
	hexa_finalize(&mesh);

	return 0;
}

//void add_pml_elements(hexa_tree_t *mesh) {
//
//}
