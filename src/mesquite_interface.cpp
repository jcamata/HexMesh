/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include <sstream>
#include <vector>
#include <mpi.h>
#include <assert.h>

#include "MeshImpl.hpp"
#include "MeshUtil.hpp"
#include "MsqTimer.hpp"
#include "MsqDebug.hpp"
#include "Mesquite.hpp"
#include "MsqError.hpp"
#include "Vector3D.hpp"
#include "InstructionQueue.hpp"
#include "LaplaceWrapper.hpp"
#include "PatchData.hpp"
#include "TerminationCriterion.hpp"
#include "QualityAssessor.hpp"

/* Mesquite includes */
#include "ParallelMeshImpl.hpp"
#include "ParallelHelper.hpp"
#include "MeshImplTags.hpp"


// algorithms
#include "Randomize.hpp"
#include "ConditionNumberQualityMetric.hpp"
#include "UntangleBetaQualityMetric.hpp"
#include "LPtoPTemplate.hpp"
#include "LInfTemplate.hpp"
#include "SteepestDescent.hpp"
#include "ConjugateGradient.hpp"
#include "PlanarDomain.hpp"


#include "hexa.h"

using namespace Mesquite;


void MeshOpt(hexa_tree_t* mesh, std::vector<double> coords, std::vector<int> material_fixed_nodes)
{
    Mesquite::MsqPrintError err(std::cout);
    
    int nvertices = mesh->local_n_nodes;
    int nelem     = mesh->local_n_elements;
    
    std::vector<int> conn(8*nelem);
    bool   *fixed_nodes = (bool*)malloc(nvertices*sizeof(bool));
    size_t *gid         = (size_t*) malloc(nvertices*sizeof(size_t));
    int    *pid         = (int*) malloc (nvertices*sizeof(int));
    
    //std::vector<bool> fixed_nodes(nvertices);
    
    int c = 0;
    int tmp[8];
    for(int  i =0; i < mesh->elements.elem_count; i++)
    {
        octant_t * elem = (octant_t*) sc_array_index(&mesh->elements, i);
        for(int j = 0; j < 8; j++) 
        {
            tmp[j] = elem->nodes[j].id;
            //c++;
        }
        conn[c] = tmp[4]; c++;
        conn[c] = tmp[5]; c++;
        conn[c] = tmp[6]; c++;
        conn[c] = tmp[7]; c++;
        conn[c] = tmp[0]; c++;
        conn[c] = tmp[1]; c++;
        conn[c] = tmp[2]; c++;
        conn[c] = tmp[3]; c++;
    }
    
    assert(mesh->nodes.elem_count == mesh->local_n_nodes);
    

    int bnd_nodes = 0;
    for(int i = 0; i < mesh->nodes.elem_count; i++)
    {
        fixed_nodes[i] = false;
        gid[i] = (size_t) mesh->global_id[i];
        pid[i] = mesh->part_nodes[i];
        octant_node_t * node = (octant_node_t*) sc_array_index(&mesh->nodes, i);
        assert(node != NULL);
        
        
        if(node->x == 0 || node->y == 0 || node->z == 0 || 
           node->x == mesh->ncellx || node->y == mesh->ncelly || node->z == (mesh->ncellz-1) )
        //if(node->z == 0 || node->z == (mesh->ncellz-1) )
        {
            fixed_nodes[i] = true; 
            bnd_nodes++;
        }
        
        if(node->z == mesh->ncellz)
            std::cout<<"Nodes: " << node->id << std::endl;
    } 
    
    //std::cout<<"BND Nodes: " << bnd_nodes << std::endl;
    
    for(int i = 0; i < material_fixed_nodes.size(); i++)
    {
        int id = material_fixed_nodes[i];
        fixed_nodes[id] = true;
    }
    
    
    Mesquite::MeshImpl mesq_mesh(nvertices,nelem,Mesquite::HEXAHEDRON, &fixed_nodes[0], &coords[0], &conn[0]);
    
    int default_gid = -1;
    int default_pid = 0;
    mesq_mesh.tag_create("GLOBAL_ID"   ,Mesquite::Mesh::HANDLE,1,&default_gid, err);
    mesq_mesh.tag_create("PROCESSOR_ID",Mesquite::Mesh::INT   ,1,&default_pid, err);
    
    TagHandle tag_processor_id = mesq_mesh.tag_get("PROCESSOR_ID", err);
    TagHandle tag_global_id    = mesq_mesh.tag_get("GLOBAL_ID", err);
    std::vector<ParallelMeshImpl::VertexHandle> vertices;   
     
    mesq_mesh.get_all_vertices(vertices, err);
    mesq_mesh.tag_set_vertex_data(tag_global_id   ,vertices.size(),&vertices[0], gid, err);
    mesq_mesh.tag_set_vertex_data(tag_processor_id,vertices.size(),&vertices[0],pid,err);
   
    
    { 
        std::ostringstream out_name;
        out_name << "original_mesh." << mesh->mpi_size << "." << mesh->mpi_rank << ".vtk";
        mesq_mesh.write_vtk(out_name.str().c_str(), err);
    }
    
    
/*    
    Mesquite::ParallelMeshImpl parallel_mesh(&mesq_mesh,"GLOBAL_ID","PROCESSOR_ID");
    Mesquite::ParallelHelperImpl helper;
    helper.set_communicator(MPI_COMM_WORLD);
    helper.set_parallel_mesh(&parallel_mesh);
    parallel_mesh.set_parallel_helper(&helper);
*/
    // do Laplacian smooth 
    Mesquite::LaplaceWrapper optimizer;
    optimizer.set_vertex_movement_limit_factor(1.e-6);
    optimizer.set_iteration_limit(100);
    optimizer.enable_culling(false);
    optimizer.run_instructions(&mesq_mesh, err);
    if (err) {std::cerr << err << std::endl; }

    std::ostringstream out_name;

    out_name << "smoothed_mesh." << mesh->mpi_size << "." << mesh->mpi_rank << ".vtk";
    mesq_mesh.write_vtk(out_name.str().c_str(), err);
    
    free(fixed_nodes);
    free(pid);
    free(gid);
    
}