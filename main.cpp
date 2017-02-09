/* 
 * File:   main.cpp
 * Author: camata
 *
 * Created on March 19, 2015, 1:42 PM
 */

#include <cstdlib>
#include <vector>

#include <sc.h>
#include <sc_io.h>
#include <sc_containers.h>
#include "hexa.h"
#include "hilbert.h"

void GetMeshFromSurface(hexa_tree_t* tree, const char* surface_topo, std::vector<double>& coords);
void GetInterceptedElements(hexa_tree_t* mesh, std::vector<double>& coords, std::vector<int>& elements_ids, const char* surface_bathy);
void CheckTemplate(hexa_tree_t* mesh, const std::vector<double>& coords, std::vector<int>& elements_ids, bool flag);

//void ChangeTemplate(hexa_tree_t* mesh, std::vector<double>& coords, std::vector<int>& elements_ids);
//void ApplyTemplate(hexa_tree_t* mesh, std::vector<double>& coords, std::vector<int>& elements_ids);

void Material_apply(hexa_tree_t *mesh, std::vector<double>& coords,std::vector<int>& element_ids, const char* surface_bathy);


void AddPMLElements(hexa_tree_t* mesh);
void ExtrudePMLElements(hexa_tree_t* mesh, std::vector<double>& coords);

/*
 * 
 */

int main(int argc, char** argv) {

	hexa_tree_t mesh;

	std::vector<double> coords;
	std::vector<int> element_ids;

	int l = atoi(argv[1]);

	hexa_init(argc, argv, &mesh);

	hexa_tree_init(&mesh, l);
	hexa_tree_cube(&mesh);

	//hexa_debug_face_hanging(&mesh);
	AddPMLElements(&mesh);

	hexa_mesh(&mesh);
	GetMeshFromSurface(&mesh, "./input/topo_Arg_ref.gts", coords);
	GetInterceptedElements(&mesh, coords, element_ids, "./input/bathy_Arg_ref.gts");
	//GetMeshFromSurface(&mesh, "./input/topo_Arg_small.gts", coords);
	//GetInterceptedElements(&mesh, coords, element_ids, "./input/bathy_Arg_small.gts");
	//GetMeshFromSurface(&mesh, "./input/topo_Pipo_ref.gts", coords);
	//GetInterceptedElements(&mesh, coords, element_ids, "./input/bathy_Pipo_ref.gts");
	printf(" Elements intercepted: %d\n", element_ids.size());

	Material_apply(&mesh, coords, element_ids, "./input/bathy_Arg_small.gts");
//	printf("Check Template \n");
	//CheckTemplate(&mesh, coords, element_ids, true);
//	printf(" Elements ref: %d\n", element_ids.size());

	//ApplyTemplate(&mesh, coords, element_ids);

	//ExtrudePMLElements(&mesh,coords);

	hexa_mesh_write_vtk(&mesh, "mesh", &coords);
	hexa_mesh_write_msh(&mesh, "mesh", &coords);

	//hexa_mesh_write_vtk(&mesh,"template", NULL);
	//hexa_mesh_write_msh(&mesh,"teste", NULL);


	//hexa_mesh_write_unv(&mesh, "mesh", &coords);
	//hexa_mesh_write_unv(&mesh,"teste", NULL);

	hexa_tree_destroy(&mesh);
	hexa_finalize(&mesh);
	return 0;
}

//void add_pml_elements(hexa_tree_t *mesh) {
//
//}
