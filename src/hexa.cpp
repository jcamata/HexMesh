
#include <math.h>
#include "hexa.h"

// Compute node id based on cartesian coordinates.
inline int get_node_id(int nx, int ny, int i, int j, int k)
{
        return (k*((nx+1)*(ny+1))+j*(nx+1)+i+1);
}


inline int get_hexa_id(int nx, int ny, int i, int j, int k)
{
    return (k*((nx)*(ny))+j*(nx)+i+1);
}



void hexa_tree_init(hexa_tree_t* mesh)
{
    sc_array_init(&mesh->elements, sizeof(octant_t));    
    //mesh->hexa_pool      =  sc_mempool_new(sizeof(octant_t));
    //mesh->hexa_node_pool =  sc_mempool_new(sizeof(octant_node_t));
}

void hexa_element_init(octant_t *elem)
{
    elem->level = 0;
    elem->x     = elem->y = elem->z;
    sc_array_init(&elem->nodes, sizeof(octant_node_t));
    sc_array_resize(&elem->nodes, 8);
}

void hexa_element_conn(octant_t* h, int i, int j, int k, int step,  int level)
{
    octant_node_t *node;
    h->x = i;
    h->y = j;
    h->z = k;
    h->level = level;
    //node 1
    node = (octant_node_t*) sc_array_index(&h->nodes,0);
    node->x = i;
    node->y = j;
    node->z = k;
    
    //node 2
    node = (octant_node_t*) sc_array_index(&h->nodes,1);
    node->x = i+step;
    node->y = j;
    node->z = k;
    
    //node 3
    node = (octant_node_t*) sc_array_index(&h->nodes,2);
    node->x = i+step;
    node->y = j+step;
    node->z = k;
    
    //node 4
    node = (octant_node_t*) sc_array_index(&h->nodes,3);
    node->x = i      ;
    node->y = j+step;
    node->z = k;
    
    //node 5
    node = (octant_node_t*) sc_array_index(&h->nodes,4);
    node->x = i;
    node->y = j;
    node->z = k+step;
    
    //node 6
    node = (octant_node_t*) sc_array_index(&h->nodes,5);
    node->x = i+step;
    node->y = j;
    node->z = k+step;
    
    //node 7
    node = (octant_node_t*) sc_array_index(&h->nodes,6);
    node->x = i+step;
    node->y = j+step;
    node->z = k+step;
    
    //node 8
    node = (octant_node_t*) sc_array_index(&h->nodes,7);
    node->x = i     ;
    node->y = j+step;
    node->z = k+step;
   
}

void hexa_get_27tree(hexa_tree_t* mesh, int nx, int ny, int nz, int step, int level)
{
    //int x,y,z;
    for(int k=0, z=nz; k < 3; k++, z+=level)
        for(int j=0, y=ny; j < 3;j++, y+=level)
            for(int i=0,x=nx ; i < 3; i++, x+=level)
            {
                octant_t * elem = (octant_t*) sc_array_push(&mesh->elements);
                hexa_element_init(elem);
                hexa_element_conn(elem,x,y,z, step, level);
            }
}


void hexa_refinement_layer(hexa_tree_t* mesh, int nz, int coarse_step, int internal_step, int level)
{
    for(int ny=0; ny < mesh->ncelly; ny+=coarse_step) {
            for(int nx=0; nx < mesh->ncellx; nx+=coarse_step)
            {
                hexa_get_27tree(mesh,nx,ny,nz,internal_step,level);
            }
        }
    
}


void hexa_uniform_layer(hexa_tree_t* mesh, int nz, int coarse_step, int internal_step, int level)
{
     for(int ny=mesh->y_start; ny < mesh->y_end; ny +=coarse_step) {
            for(int nx=mesh->x_start; nx < mesh->x_end; nx+=coarse_step)
            {
                octant_t * elem = (octant_t*) sc_array_push(&mesh->elements);
                hexa_element_init(elem);
                hexa_element_conn(elem,nx,ny,nz, coarse_step, level);
            }
     }
}


void hexa_transient_layer(hexa_tree_t* mesh, int nz, int coarse_step, int internal_step, int level)
{
     for(int ny=mesh->y_start; ny < mesh->y_end; ny +=coarse_step) {
            for(int nx=mesh->x_start; nx < mesh->x_end; nx+=coarse_step)
            {
                hexa_transition_element(mesh,nx,ny,nz,internal_step,level);
            }
        }
    
}


void hexa_tree_cube(hexa_tree_t* mesh, int max_tree_levels)
{
    
    int32_t nx, ny, nz, coarse_step, internal_step;
    
    mesh->ncellx = (int32_t) pow(3, max_tree_levels);
    mesh->ncelly = mesh->ncellx;
    mesh->ncellz = mesh->ncellx;
    mesh->max_levels = max_tree_levels;
    
    internal_step = 1;
    coarse_step   = 1;
    int level     = 1;
    int layer     = 0;
    int nlayer    = 0; 
    
    hexa_processors_interval(mesh);
    
    nz = 0;
    int nz_test = nz+internal_step;
    while( (nz+internal_step) <= mesh->ncellz)
    {
        if((nlayer+1)%9 == 0) {
            coarse_step*=3;
            hexa_transient_layer(mesh,nz,coarse_step, internal_step, level);
            internal_step*=3;
        } else {
            hexa_uniform_layer(mesh,nz,coarse_step, internal_step,level);
        }
        
        nz+=internal_step;
        
        nlayer++;
    }
    
    //printf("Number of octants: %ld\n", mesh->elements.elem_count);
  
}



void hexa_tree_destroy(hexa_tree_t* mesh)
{
    sc_array_reset(&mesh->elements);
    sc_array_reset(&mesh->nodes);
    //sc_mempool_destroy (mesh->hexa_pool);
    //sc_mempool_destroy(mesh->hexa_node_pool);
}




