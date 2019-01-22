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


int CheckTemplate(hexa_tree_t* mesh, const std::vector<double>& coords, std::vector<int>& elements_ids, bool flag);
void CutTemplate(hexa_tree_t* mesh, std::vector<double>& coords, std::vector<int>& elements_ids);

void Apply_material(hexa_tree_t* mesh, std::vector<double>& coords,std::vector<int>& element_ids, const char* surface_bathy);
void ExtrudePMLElements(hexa_tree_t* mesh, std::vector<double>& coords);

void IdentifyTemplate(hexa_tree_t* mesh, const std::vector<double>& coords, std::vector<int>& elements_ids);
void MovingNodes(hexa_tree_t* mesh, std::vector<double>& coords, std::vector<int>& nodes_b_mat, const char* surface);
void MeshOpt(hexa_tree_t* mesh, std::vector<double>& coords, std::vector<int> material_fixed_nodes);
void UntagleMesh(hexa_tree_t* mesh, std::vector<double>& coords, std::vector<int> material_fixed_nodes);
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

	//const char * bathy = "./input/bathy_Arg_small.gts";
	//const char * bathy = "./input/bathy_Arg_ref.gts";
	const char * bathy = "./input/bathy_Pipo_small.gts";
	//const char * bathy = "./input/bathy_pipo.gts";
	//const char * bathy = "./input/bathy_test.gts";

	//const char * topo = "./input/topo_Arg_small.gts";
	//const char * topo = "./input/topo_Arg_ref.gts";
	const char * topo =  "./input/topo_Pipo_small.gts";
	//const char * topo =  "./input/topo_pipo.gts";
	//const char * topo = "./input/topo_test.gts";

	// Note that here we use a gts file.
	// There is a tool called stl2gts that convert STL files to GTS.
	// It is installed together with the gts library.
	GetMeshFromSurface(&mesh, topo, coords);

	GetInterceptedElements(&mesh, coords, element_ids, bathy);
	printf(" Elements intercepted: %lld\n\n", element_ids.size());

	//printf(" Check the method\n");
	//el_not_handle = CheckTemplate(&mesh, coords, element_ids ,true);

	//if(el_not_handle == 0){
	//printf("Cut templates\n");
	//CutTemplate(&mesh, coords, element_ids);
	//}else{
/*
	printf(" Check and propagate 27-tree templates\n");
	CheckOctreeTemplate(&mesh, coords, element_ids, true);

	printf(" Apply 27-tree templates\n");
	ApplyOctreeTemplate(&mesh, coords, element_ids);
*/
    //redo the mesh conectivity
    //sem ela não posso aplicar o moving nodes
	element_ids.clear();
	GetInterceptedElements(&mesh, coords, element_ids, bathy);
	printf(" Elements intercepted: %lld\n\n", element_ids.size());

	printf(" Project nodes to the surface\n\n");
	MovingNodes(&mesh,coords, nodes_b_mat,bathy);

	printf(" Applying material \n\n");
	element_ids.clear();
	Apply_material(&mesh, coords, element_ids, bathy);

	//}

	//printf(" Untangle meshes\n\n");
	//UntagleMesh(&mesh, coords, nodes_b_mat);

	//printf(" Optimization of the mesh\n\n");
	//MeshOpt(&mesh,coords,nodes_b_mat);

	printf(" Extrude elements\n\n");
	ExtrudePMLElements(&mesh,coords);

	//clean vectors
	std::vector<int>().swap(element_ids);
	std::vector<int>().swap(nodes_b_mat);

	printf(" Writing output files \n\n");
	hexa_mesh_write_vtk(&mesh, "mesh", &coords);
	//hexa_mesh_write_msh(&mesh, "mesh", &coords);
	//hexa_mesh_write_h5(&mesh,"mesh", coords);

	//hexa_mesh_write_vtk(&mesh, "test",NULL);

	printf(" Cleaning variables \n\n");

	//hexa_mesh_destroy(&mesh);
	hexa_tree_destroy(&mesh);
	hexa_finalize(&mesh);
	std::vector<double>().swap(coords);

	return 0;
}

