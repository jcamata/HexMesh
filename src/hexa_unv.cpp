
#include <iostream>
#include <vector>
#include <fstream>
#include <iomanip>
using namespace std;

#include "hexa.h" 


void UNVIO_WriteNodes(std::ofstream &out_file, hexa_tree_t* mesh)
{
    int n_nodes       = mesh->nodes.elem_count;
    sc_array_t* nodes = &mesh->nodes;
    
    int _nodes_dataset_label = 2411;
    // Write beginning of dataset
    out_file << "    -1\n"
             << "  "
             << _nodes_dataset_label
             << '\n';
     unsigned int
        exp_coord_sys_dummy  = 1, // export coordinate sys. (not supported yet)
        disp_coord_sys_dummy = 1, // displacement coordinate sys. (not supp. yet)
        color_dummy          = 11; // color(not supported yet)
  
     // Use scientific notation with captial E and 16 digits for printing out the coordinates
    out_file << std::scientific << std::setprecision(16) << std::uppercase;
    
    for(int inode = 0; inode < n_nodes; ++inode)
    {
        octant_node_t* n = (octant_node_t*) sc_array_index(nodes,inode);
        
        out_file << std::setw(10) << (inode+1)
                 << std::setw(10) << exp_coord_sys_dummy
                 << std::setw(10) << disp_coord_sys_dummy
                 << std::setw(10) << color_dummy
                 << '\n';
        out_file << std::setw(25) << (double) n->x
                 << std::setw(25) << (double) n->y
                 << std::setw(25) << (double) n->z
                 << '\n';
    }
    // Write end of dataset
    out_file << "    -1\n";
}

void UNVIO_WriteElements(std::ofstream &out_file, hexa_tree_t* mesh)
{
    int n_elements       = mesh->elements.elem_count;
    sc_array_t* elements = &mesh->elements;
    
    int  _elements_dataset_label = 2412;
    // Write beginning of dataset
    out_file << "    -1\n"
             << "  "
             << _elements_dataset_label
             << "\n";
     unsigned
        fe_descriptor_id    = 0,    // FE descriptor id
        phys_prop_tab_dummy = 1, // physical property (not supported yet)
        mat_prop_tab_dummy  = 0,  // material property 
        color_dummy         = 7;         // color (not supported yet)
     
     unsigned int assign_elem_nodes[8];
     fe_descriptor_id = 115; // Solid Linear Brick
     assign_elem_nodes[0] = 0;
     assign_elem_nodes[1] = 4;
     assign_elem_nodes[2] = 5;
     assign_elem_nodes[3] = 1;
     assign_elem_nodes[4] = 3;
     assign_elem_nodes[5] = 7;
     assign_elem_nodes[6] = 6;
     assign_elem_nodes[7] = 2;
     
     for(int  iel = 0; iel < n_elements; ++iel)
     {
          octant_t* h = (octant_t*) sc_array_index(elements, iel);
          
          out_file << std::setw(10)   << iel            // element ID
                << std::setw(10) << fe_descriptor_id    // type of element
                << std::setw(10) << phys_prop_tab_dummy // not supported
                << std::setw(10) << mat_prop_tab_dummy  // not supported
                << std::setw(10) << color_dummy         // not supported
                << std::setw(10) << 8                   // No. of nodes per element
                << '\n';
          
          for(int j =0; j < 8; j++)
          {
              octant_node_t* n = (octant_node_t*) sc_array_index(&h->nodes, assign_elem_nodes[j]);
              out_file << std::setw(10) << (n->id+1);
          }
          out_file << '\n';
     }
     // Write end of dataset
    out_file << "    -1\n";
}

void hexa_mesh_write_unv(hexa_tree_t* mesh, const char* root_name)
{

    char filename[80];
    sprintf(filename,"%s_%04d_%04d.unv",root_name, mesh->mpi_size, mesh->mpi_rank);
    std::ofstream out_file (filename);
    if(! out_file.good())
    {
        cerr << "ERROR: output file not good" << endl;
        return;
    }
        
    UNVIO_WriteNodes(out_file, mesh);
    UNVIO_WriteElements(out_file,mesh);
    out_file.close();
       
}