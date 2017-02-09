#include <iostream>
#include <vector>
#include <fstream>
#include <iomanip>
using namespace std;

#include "hexa.h"

void hexa_mesh_write_msh(hexa_tree_t* mesh, const char* root_name, std::vector<double> *coords = NULL) {

	int element;
	int element_type = 5;
	int i;
	int i2;
	int node;
	int tag_num;
	int tag1;

	int m = 3;
	int element_order = 1;
	int element_node[3];

	char filename[80];
	sprintf(filename, "%s_%04d_%04d.msh",root_name , mesh->mpi_size, mesh->mpi_rank);

	std::ofstream out_file(filename);
	if (!out_file.good()) {
		cerr << "ERROR: output file not good" << endl;
		return;
	}

	//
	//  Write the data.
	//
	out_file << "$MeshFormat\n";
	out_file << "2.2 0 8\n";
	out_file << "$EndMeshFormat\n";

	int n_nodes = mesh->nodes.elem_count;
	sc_array_t* nodes = &mesh->nodes;

	out_file << "$Nodes\n";
	out_file << n_nodes << "\n";

	if (coords == NULL) {
		for (int node = 0; node < n_nodes; ++node) {
			octant_node_t* n = (octant_node_t*) sc_array_index(nodes, node);

			out_file << node + 1
					<< "  " << (double) n->x
					<< "  " << (double) n->y
					<< "  " << (double) n->z << "\n";

		}
	} else {
		for ( node = 0; node < n_nodes; node++ ){

			//out_file << node + 1
			//		<< "  " << (double) ((*coords)[node * 3 + 0]-(*coords)[0 * 3 + 0])
			//		<< "  " << (double) ((*coords)[node * 3 + 1]-(*coords)[0 * 3 + 1])
			//		<< "  " << (double) (*coords)[node * 3 + 2] << "\n";

			out_file << node + 1
					<< "  " << (double) (*coords)[node * 3 + 0]
					<< "  " << (double) (*coords)[node * 3 + 1]
					<< "  " << (double) (*coords)[node * 3 + 2] << "\n";

		}
	}

	out_file << "$EndNodes\n";

	int n_elements = mesh->elements.elem_count;
	sc_array_t* elements = &mesh->elements;

	unsigned int assign_elem_nodes[8];

	if (coords == NULL) {
		assign_elem_nodes[0] = 0;
		assign_elem_nodes[1] = 1;
		assign_elem_nodes[2] = 2;
		assign_elem_nodes[3] = 3;
		assign_elem_nodes[4] = 4;
		assign_elem_nodes[5] = 5;
		assign_elem_nodes[6] = 6;
		assign_elem_nodes[7] = 7;
	} else {
		assign_elem_nodes[0] = 4;
		assign_elem_nodes[1] = 5;
		assign_elem_nodes[2] = 6;
		assign_elem_nodes[3] = 7;
		assign_elem_nodes[4] = 0;
		assign_elem_nodes[5] = 1;
		assign_elem_nodes[6] = 2;
		assign_elem_nodes[7] = 3;
	}

	tag_num = 2;
	tag1 = 0;
	out_file << "$Elements\n";
	out_file << n_elements << "\n";
	for ( element = 0; element < n_elements; element++ ){
		octant_t* h = (octant_t*) sc_array_index(elements, element);

		out_file << element + 1
				<< "  " << element_type
				<< "  " << tag_num
				<< "  " << h->n_mat
				<< "  " << h->pad;

		for (int j = 0; j < 8; j++) {
			octant_node_t* n = &h->nodes[assign_elem_nodes[j]];
			out_file << "  " <<  (n->id + 1);
		}

		out_file << "\n";
	}

	out_file << "$EndElements\n";

	out_file.close ( );

	return;
}
