
#include <iostream>
#include <vector>
#include <fstream>
#include <iomanip>
using namespace std;

#include "hexa.h" 

void UNVIO_WriteNodes(std::ofstream &out_file, hexa_tree_t* mesh, vector<double> *coords = NULL) {
    int n_nodes = mesh->nodes.elem_count;
    sc_array_t* nodes = &mesh->nodes;

    int _nodes_dataset_label = 2411;
    // Write beginning of dataset
    out_file << "    -1\n"
            << "  "
            << _nodes_dataset_label
            << '\n';
    unsigned int
    exp_coord_sys_dummy = 1, // export coordinate sys. (not supported yet)
            disp_coord_sys_dummy = 1, // displacement coordinate sys. (not supp. yet)
            color_dummy = 11; // color(not supported yet)

    // Use scientific notation with captial E and 16 digits for printing out the coordinates
    out_file << std::scientific << std::setprecision(16) << std::uppercase;

    if (coords == NULL) {
        for (int inode = 0; inode < n_nodes; ++inode) {
            octant_node_t* n = (octant_node_t*) sc_array_index(nodes, inode);

            out_file << std::setw(10) << (inode + 1)
                    << std::setw(10) << exp_coord_sys_dummy
                    << std::setw(10) << disp_coord_sys_dummy
                    << std::setw(10) << color_dummy
                    << '\n';
            out_file << std::setw(25) << (double) n->x
                    << std::setw(25) << (double) n->y
                    << std::setw(25) << (double) n->z
                    << '\n';
        }
    } else {
        for (int inode = 0; inode < n_nodes; ++inode) {
            out_file << std::setw(10) << (inode + 1)
                    << std::setw(10) << exp_coord_sys_dummy
                    << std::setw(10) << disp_coord_sys_dummy
                    << std::setw(10) << color_dummy
                    << '\n';
            out_file << std::setw(25) << (double) (*coords)[inode * 3 + 0]
                    << std::setw(25) << (double) (*coords)[inode * 3 + 1]
                    << std::setw(25) << (double) (*coords)[inode * 3 + 2]
                    << '\n';
        }
    }
    // Write end of dataset
    out_file << "    -1\n";
}

void UNVIO_WriteElements(std::ofstream &out_file, hexa_tree_t* mesh) {
    int n_elements = mesh->elements.elem_count;
    sc_array_t* elements = &mesh->elements;

    int _elements_dataset_label = 2412;
    // Write beginning of dataset
    out_file << "    -1\n"
            << "  "
            << _elements_dataset_label
            << "\n";
    unsigned
    fe_descriptor_id = 0, // FE descriptor id
            phys_prop_tab_dummy = 1, // physical property (not supported yet)
            mat_prop_tab_dummy = 0, // material property 
            color_dummy = 7; // color (not supported yet)

    unsigned int assign_elem_nodes[8];
    fe_descriptor_id = 115; // Solid Linear Brick
    /*
    assign_elem_nodes[0] = 0;
    assign_elem_nodes[1] = 4;
    assign_elem_nodes[2] = 5;
    assign_elem_nodes[3] = 1;
    assign_elem_nodes[4] = 3;
    assign_elem_nodes[5] = 7;
    assign_elem_nodes[6] = 6;
    assign_elem_nodes[7] = 2;
 * 
     
    assign_elem_nodes[0] = 4; 
    assign_elem_nodes[1] = 5; 
    assign_elem_nodes[2] = 6;
    assign_elem_nodes[3] = 7;
    assign_elem_nodes[4] = 0;
    assign_elem_nodes[5] = 1;
    assign_elem_nodes[6] = 2;
    assign_elem_nodes[7] = 3;
    */
    assign_elem_nodes[0] = 0; 
    assign_elem_nodes[1] = 1; 
    assign_elem_nodes[2] = 2;
    assign_elem_nodes[3] = 3;
    assign_elem_nodes[4] = 4;
    assign_elem_nodes[5] = 5;
    assign_elem_nodes[6] = 6;
    assign_elem_nodes[7] = 7;
    
    
    for (int iel = 0; iel < n_elements; ++iel) {
        octant_t* h = (octant_t*) sc_array_index(elements, iel);

        out_file << std::setw(10) << iel + 1 // element ID (start in 1 not 0)
                << std::setw(10) << fe_descriptor_id // type of element
                << std::setw(10) << phys_prop_tab_dummy // not supported
                << std::setw(10) << h->n_mat//mat_prop_tab_dummy  // not supported
                << std::setw(10) << color_dummy // not supported
                << std::setw(10) << 8 // No. of nodes per element
                << '\n';

        for (int j = 0; j < 8; j++) {
            octant_node_t* n = &h->nodes[assign_elem_nodes[j]];
            out_file << std::setw(10) << (n->id + 1);
        }
        out_file << '\n';
    }
    // Write end of dataset
    out_file << "    -1\n";
}

void UNVIO_WritePhysicalVolumes(std::ofstream &out_file, hexa_tree_t* mesh) {
    
    sc_array_t* elements = &mesh->elements;

    int n_elements = mesh->elements.elem_count;
    int total_n_mat = mesh->total_n_mat;
    int _physicalVolume_dataset_label = 2477;
    // Write beginning of dataset
    out_file << "    -1\n"
            << "  "
            << _physicalVolume_dataset_label<< "\n"
            << "    -1\n"
            << "    -1\n"
            << "    -1\n";
    unsigned
    active_constraint_set = 0, // physical property (not supported yet)
            active_restraint_set = 0, // material property 
            active_load_set = 0, // color (not supported yet)
            active_dof_set = 0,
            active_temperature_set = 0,
            active_contact_set = 0,
            entity_type_code = 8,
            entity_node_leaf_id = 0,
            entity_component_id = 0;

    std::vector<int> count(total_n_mat);
    std::vector<int> el_tab(n_elements);
    int counter = 1;
    for (int n_mater = 1; n_mater <= total_n_mat; ++n_mater) {
        count[n_mater] = 0;
        for (int iel = 0; iel < n_elements; ++iel) {
            octant_t* h = (octant_t*) sc_array_index(elements, iel);
            if (n_mater == h->n_mat) {
                count[n_mater] = count[n_mater] + 1;
                el_tab[ counter ] = iel + 1;
                counter = counter + 1;
            }
        }
    }

    counter = 1;

    for (int n_mater = 1; n_mater <= total_n_mat; ++n_mater) {

        if (count[n_mater] != 0) {
            //describe the data set
            out_file << std::setw(10) << n_mater//  set number
                    << std::setw(10) << active_constraint_set // 
                    << std::setw(10) << active_restraint_set // 
                    << std::setw(10) << active_load_set // 
                    << std::setw(10) << active_dof_set // 
                    << std::setw(10) << active_temperature_set // 
                    << std::setw(10) << active_contact_set // 
                    << std::setw(10) << count[n_mater]// number of elements in the set
                    << '\n';
            //set name
            out_file << "PhysicalVolume"<<n_mater // set name
                    << '\n';
            //set components
            if (count[n_mater] % 2 == 0) {
                for (int els = counter; els < counter + count[n_mater]; els = els + 2) {
                    int tag_id_1 = el_tab[ 2 * els ];
                    int tag_id_2 = el_tab[ 2 * els + 1 ];
                    out_file << std::setw(10) << entity_type_code//  set number
                            << std::setw(10) << tag_id_1 // 
                            << std::setw(10) << entity_node_leaf_id // 
                            << std::setw(10) << entity_component_id //
                            << std::setw(10) << entity_type_code//  set number
                            << std::setw(10) << tag_id_2 // 
                            << std::setw(10) << entity_node_leaf_id // 
                            << std::setw(10) << entity_component_id //
                            << '\n';
                }
            }

            if (count[n_mater] % 2 != 0) {
                if (count[n_mater] > 1) {
                    for (int els = counter; els < counter + count[n_mater] - 1; els = els + 2) {
                        int tag_id_1 = el_tab[ els ];
                        int tag_id_2 = el_tab[ els + 1];
                        out_file << std::setw(10) << entity_type_code//  set number
                                << std::setw(10) << tag_id_1 // 
                                << std::setw(10) << entity_node_leaf_id // 
                                << std::setw(10) << entity_component_id //
                                << std::setw(10) << entity_type_code//  set number
                                << std::setw(10) << tag_id_2 // 
                                << std::setw(10) << entity_node_leaf_id // 
                                << std::setw(10) << entity_component_id //
                                << '\n';
                    }
                }
                int tag_id_1 = el_tab[ counter + count[n_mater] - 1 ];
                out_file << std::setw(10) << entity_type_code//  set number
                        << std::setw(10) << tag_id_1 // 
                        << std::setw(10) << entity_node_leaf_id // 
                        << std::setw(10) << entity_component_id //
                        << '\n';
            }
            counter = counter + count[n_mater];
        }
    }
    // Write end of dataset
    out_file << "    -1\n";
}

void hexa_mesh_write_unv(hexa_tree_t* mesh, const char* root_name, std::vector<double> *coords = NULL) {


    char filename[80];
    sprintf(filename, "%s_%04d_%04d.unv", root_name, mesh->mpi_size, mesh->mpi_rank);
    std::ofstream out_file(filename);
    if (!out_file.good()) {
        cerr << "ERROR: output file not good" << endl;
        return;
    }

    UNVIO_WriteNodes(out_file, mesh, coords);
    UNVIO_WriteElements(out_file, mesh);
    UNVIO_WritePhysicalVolumes(out_file, mesh);
    out_file.close();

}