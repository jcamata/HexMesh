
#include <stdlib.h>
#include <stdio.h>
#include <sc_io.h>
#include <sc_containers.h>

#include "hexa.h"
#include "pml.h"

void hexa_element_init(octant_t *elem);


void hexa_transition_element(hexa_tree_t* mesh, int i, int j, int k, int step, int level)
{
     octant_node_t *node;
     octant_t * h;
     // Creating 13 news hexahedra
   
     //ELEment #1
     h = (octant_t*) sc_array_push(&mesh->elements);
     hexa_element_init(h);
     h->x = i;
     h->y = j;
     h->z = k;
     h->level = level;
             
     //----------------------------
     h->pad = 0;
     //SetElemPML(mesh,h,step);
    
     //node 1
     node = &h->nodes[0];
     node->x = i;
     node->y = j;
     node->z = k;
    
     //node 2
     node = &h->nodes[1];
     node->x = i+step;
     node->y = j;
     node->z = k;
    
     //node 3
     node = &h->nodes[2];
     node->x = i+step;
     node->y = j+step;
     node->z = k;
    
    //node 4
    node = &h->nodes[3];
    node->x = i      ;
    node->y = j+step;
    node->z = k;
    
    //node 5
    node = &h->nodes[4];
    node->x = i;
    node->y = j;
    node->z = k+3*step;
    
    //node 6
    node = &h->nodes[5];
    node->x = i+step;
    node->y = j;
    node->z = k+2*step;
    
    //node 7
    node = &h->nodes[6];
    node->x = i+step;
    node->y = j+step;
    node->z = k+step;
    
    //node 8
    node = &h->nodes[7];
    node->x = i       ;
    node->y = j+step  ;
    node->z = k+2*step;
    


    // Elemento 2 
     h = (octant_t*) sc_array_push(&mesh->elements);
     hexa_element_init(h);
     h->x = i+step;
     h->y = j;
     h->z = k;
     h->level = level;
     h->pad = 0;
     //SetElemPML(mesh,h,step);
     
          //node 1
     node = &h->nodes[0];
     node->x = i+step;
     node->y = j;
     node->z = k;
    
     //node 2
     node = &h->nodes[1];
     node->x = i+2*step;
     node->y = j;
     node->z = k;
    
     //node 3
     node = &h->nodes[2];
     node->x = i+2*step;
     node->y = j+step;
     node->z = k;
    
    //node 4
    node = &h->nodes[3];
    node->x = i+step      ;
    node->y = j+step;
    node->z = k;
    
    
    //node 5
    node = &h->nodes[4];
    node->x = i+step;
    node->y = j;
    node->z = k+2*step;
    
    //node 6
    node = &h->nodes[5];
    node->x = i+2*step;
    node->y = j;
    node->z = k+2*step;
    
    //node 7
    node = &h->nodes[6];
    node->x = i+2*step;
    node->y = j+step;
    node->z = k+step;
    
    //node 8
    node = &h->nodes[7];
    node->x = i+step  ;
    node->y = j+step  ;
    node->z = k+step;
    

    
    // Elemento 3
     h = (octant_t*) sc_array_push(&mesh->elements);
     hexa_element_init(h);
     h->x = i+2*step;
     h->y = j;
     h->z = k;
     h->level = level;
     h->pad = 0;
     //SetElemPML(mesh,h,step);
     
     //node 1
     node = &h->nodes[0];
     node->x = i+2*step;
     node->y = j;
     node->z = k;
    
     //node 2
     node = &h->nodes[1];
     node->x = i+3*step;
     node->y = j;
     node->z = k;
    
     //node 3
     node = &h->nodes[2];
     node->x = i+3*step;
     node->y = j+step;
     node->z = k;
    
    //node 4
    node = &h->nodes[3];
    node->x = i+2*step      ;
    node->y = j+step;
    node->z = k;
    
    
    //node 5
    node = &h->nodes[4];
    node->x = i+2*step;
    node->y = j;
    node->z = k+2*step;
    
    //node 6
    node = &h->nodes[5];
    node->x = i+3*step;
    node->y = j;
    node->z = k+3*step;
    
    //node 7
    node = &h->nodes[6];
    node->x = i+3*step;
    node->y = j+step;
    node->z = k+2*step;
    
    //node 8
    node = &h->nodes[7];
    node->x = i+2*step  ;
    node->y = j+step  ;
    node->z = k+step;
    

    
    // Elemento 4 
     h = (octant_t*) sc_array_push(&mesh->elements);
     hexa_element_init(h);
     h->x = i;
     h->y = j+step;
     h->z = k;
     h->level = level;
     h->pad = 0;
     //SetElemPML(mesh,h,step);
          //node 1
     node = &h->nodes[0];
     node->x = i;
     node->y = j+step;
     node->z = k;
    
     //node 2
     node = &h->nodes[1];
     node->x = i+step;
     node->y = j+step;
     node->z = k;
    
     //node 3
     node = &h->nodes[2];
     node->x = i+step;
     node->y = j+2*step;
     node->z = k;
    
    //node 4
    node = &h->nodes[3];
    node->x = i      ;
    node->y = j+2*step;
    node->z = k;
    
    //node 5
    node = &h->nodes[4];
    node->x = i;
    node->y = j+step;
    node->z = k+2*step;
    
    //node 6
    node = &h->nodes[5];
    node->x = i+step;
    node->y = j+step;
    node->z = k+step;
    
    //node 7
    node = &h->nodes[6];
    node->x = i+step;
    node->y = j+2*step;
    node->z = k+step;
    
    //node 8
    node = &h->nodes[7];
    node->x = i   ;
    node->y = j+2*step;
    node->z = k+2*step;
    

    
    // Elemento 5 
     h = (octant_t*) sc_array_push(&mesh->elements);
     hexa_element_init(h);
     h->x = i+step;
     h->y = j+step;
     h->z = k;
     h->level = level;
     h->pad = 0;
     //SetElemPML(mesh,h,step);
     
          //node 1
     node = &h->nodes[0];
     node->x = i+step;
     node->y = j+step;
     node->z = k;
    
     //node 2
     node = &h->nodes[1];
     node->x = i+2*step;
     node->y = j+step;
     node->z = k;
    
     //node 3
     node = &h->nodes[2];
     node->x = i+2*step;
     node->y = j+2*step;
     node->z = k;
    
    //node 4
    node = &h->nodes[3];
    node->x = i+step      ;
    node->y = j+2*step;
    node->z = k;
    
    
    //node 5
    node = &h->nodes[4];
    node->x = i+step;
    node->y = j+step;
    node->z = k+step;
    
    //node 6
    node = &h->nodes[5];
    node->x = i+2*step;
    node->y = j+step;
    node->z = k+step;
    
    //node 7
    node = &h->nodes[6];
    node->x = i+2*step;
    node->y = j+2*step;
    node->z = k+step;
    
    //node 8
    node = &h->nodes[7];
    node->x = i+step  ;
    node->y = j+2*step  ;
    node->z = k+step;



    
    // Elemento 6 
     h = (octant_t*) sc_array_push(&mesh->elements);
     hexa_element_init(h);
     h->x = i+2*step;
     h->y = j+step;
     h->z = k;
     h->level = level;
     h->pad = 0;
     //SetElemPML(mesh,h,step);
     
          //node 1
     node = &h->nodes[0];
     node->x = i+2*step;
     node->y = j+step;
     node->z = k;
    
     //node 2
     node = &h->nodes[1];
     node->x = i+3*step;
     node->y = j+step;
     node->z = k;
    
     //node 3
     node = &h->nodes[2];
     node->x = i+3*step;
     node->y = j+2*step;
     node->z = k;
    
    //node 4
    node = &h->nodes[3];
    node->x = i+2*step      ;
    node->y = j+2*step;
    node->z = k;
    
    
    //node 5
    node = &h->nodes[4];
    node->x = i+2*step;
    node->y = j+step;
    node->z = k+step;
    
    //node 6
    node = &h->nodes[5];
    node->x = i+3*step;
    node->y = j+step;
    node->z = k+2*step;
    
    //node 7
    node = &h->nodes[6];
    node->x = i+3*step;
    node->y = j+2*step;
    node->z = k+2*step;
    
    //node 8
    node = &h->nodes[7];
    node->x = i+2*step  ;
    node->y = j+2*step  ;
    node->z = k+step;
 

   
  
    // Elemento 7 
     h = (octant_t*) sc_array_push(&mesh->elements);
     hexa_element_init(h);
     h->x = i;
     h->y = j+2*step;
     h->z = k;
     h->level = level;
     h->pad = 0;
     //SetElemPML(mesh,h,step);
     
          //node 1
     node = &h->nodes[0];
     node->x = i;
     node->y = j+2*step;
     node->z = k;
    
     //node 2
     node = &h->nodes[1];
     node->x = i+step;
     node->y = j+2*step;
     node->z = k;
    
     //node 3
     node = &h->nodes[2];
     node->x = i+step;
     node->y = j+3*step;
     node->z = k;
    
    //node 4
    node = &h->nodes[3];
    node->x = i      ;
    node->y = j+3*step;
    node->z = k;
    
    
    //node 5
    node = &h->nodes[4];
    node->x = i;
    node->y = j+2*step;
    node->z = k+2*step;
    
    //node 6
    node = &h->nodes[5];
    node->x = i+step;
    node->y = j+2*step;
    node->z = k+step;
    
    //node 7
    node = &h->nodes[6];
    node->x = i+step;
    node->y = j+3*step;
    node->z = k+2*step;
    
    //node 8
    node = &h->nodes[7];
    node->x = i  ;
    node->y = j+3*step  ;
    node->z = k+3*step;
    

   
    
    // Elemento 8 
     h = (octant_t*) sc_array_push(&mesh->elements);
     hexa_element_init(h);
     h->x = i+step;
     h->y = j+2*step;
     h->z = k;
     h->level = level;
     h->pad = 0;
     //SetElemPML(mesh,h,step);
     
          //node 1
     node = &h->nodes[0];
     node->x = i+step;
     node->y = j+2*step;
     node->z = k;
    
     //node 2
     node = &h->nodes[1];
     node->x = i+2*step;
     node->y = j+2*step;
     node->z = k;
    
     //node 3
     node = &h->nodes[2];
     node->x = i+2*step;
     node->y = j+3*step;
     node->z = k;
    
    //node 4
    node = &h->nodes[3];
    node->x = i+step      ;
    node->y = j+3*step;
    node->z = k;
    
    
    //node 5
    node = &h->nodes[4];
    node->x = i+step;
    node->y = j+2*step;
    node->z = k+step;
    
    //node 6
    node = &h->nodes[5];
    node->x = i+2*step;
    node->y = j+2*step;
    node->z = k+step;
    
    //node 7
    node = &h->nodes[6];
    node->x = i+2*step;
    node->y = j+3*step;
    node->z = k+2*step;
    
    //node 8
    node = &h->nodes[7];
    node->x = i+step  ;
    node->y = j+3*step  ;
    node->z = k+2*step;
    
     
   
    
    // Elemento 9 
     h = (octant_t*) sc_array_push(&mesh->elements);
     hexa_element_init(h);
     h->x = i+step;
     h->y = j+2*step;
     h->z = k;
     h->level = level;
     h->pad = 0;
     //SetElemPML(mesh,h,step);
     
          //node 1
     node = &h->nodes[0];
     node->x = i+2*step;
     node->y = j+2*step;
     node->z = k;
    
     //node 2
     node = &h->nodes[1];
     node->x = i+3*step;
     node->y = j+2*step;
     node->z = k;
    
     //node 3
     node = &h->nodes[2];
     node->x = i+3*step;
     node->y = j+3*step;
     node->z = k;
    
    //node 4
    node = &h->nodes[3];
    node->x = i+2*step      ;
    node->y = j+3*step;
    node->z = k;
    
    
    //node 5
    node = &h->nodes[4];
    node->x = i+2*step;
    node->y = j+2*step;
    node->z = k+step;
    
    //node 6
    node = &h->nodes[5];
    node->x = i+3*step;
    node->y = j+2*step;
    node->z = k+2*step;
    
    //node 7
    node = &h->nodes[6];
    node->x = i+3*step;
    node->y = j+3*step;
    node->z = k+3*step;
    
    //node 8
    node = &h->nodes[7];
    node->x = i+2*step  ;
    node->y = j+3*step  ;
    node->z = k+2*step;
    

    
    // Elemento 10 
     h = (octant_t*) sc_array_push(&mesh->elements);
     hexa_element_init(h);
     h->x = i+step;
     h->y = j+step;
     h->z = k+step;
     h->level = level;
     h->pad = 0;
     //SetElemPML(mesh,h,step);
     
          //node 1
     node = &h->nodes[0];
     node->x = i+step;
     node->y = j+step;
     node->z = k+step;
    
     //node 2
     node = &h->nodes[1];
     node->x = i+2*step;
     node->y = j+step;
     node->z = k+step;
    
     //node 3
     node = &h->nodes[2];
     node->x = i+2*step;
     node->y = j+2*step;
     node->z = k+step;
    
    //node 4
    node = &h->nodes[3];
    node->x = i+step      ;
    node->y = j+2*step;
    node->z = k+step;
    
    
    //node 5
    node = &h->nodes[4];
    node->x = i+step;
    node->y = j;
    node->z = k+2*step;
    
    //node 6
    node = &h->nodes[5];
    node->x = i+2*step;
    node->y = j;
    node->z = k+2*step;
    
    //node 7
    node = &h->nodes[6];
    node->x = i+2*step;
    node->y = j+3*step;
    node->z = k+2*step;
    
    //node 8
    node = &h->nodes[7];
    node->x = i+step  ;
    node->y = j+3*step  ;
    node->z = k+2*step;
   
 

 
    
    // Elemento 11 
     h = (octant_t*) sc_array_push(&mesh->elements);
     hexa_element_init(h);
     h->x = i;
     h->y = j+step;
     h->z = k+2*step;
     h->level = level;
     h->pad = 0;
     //SetElemPML(mesh,h,step);
     
     //node 1
     node = &h->nodes[0];
     node->x = i;
     node->y = j+step;
     node->z = k+2*step;
    
     //node 2
     node = &h->nodes[1];
     node->x = i+step;
     node->y = j+step;
     node->z = k+step;
    
     //node 3
     node = &h->nodes[2];
     node->x = i+step;
     node->y = j+2*step;
     node->z = k+step;
    
    //node 4
    node = &h->nodes[3];
    node->x = i       ;
    node->y = j+2*step;
    node->z = k+2*step;
    
    
    //node 5
    node = &h->nodes[4];
    node->x = i;
    node->y = j;
    node->z = k+3*step;
    
    //node 6
    node = &h->nodes[5];
    node->x = i+step;
    node->y = j;
    node->z = k+2*step;
    
    //node 7
    node = &h->nodes[6];
    node->x = i+step;
    node->y = j+3*step;
    node->z = k+2*step;
    
    //node 8
    node = &h->nodes[7];
    node->x = i  ;
    node->y = j+3*step  ;
    node->z = k+3*step;

  
    
 

    // Elemento 12 
     h = (octant_t*) sc_array_push(&mesh->elements);
     hexa_element_init(h);
     h->x = i+2*step;
     h->y = j+step;
     h->z = k+step;
     h->level = level;
     h->pad = 0;
     //SetElemPML(mesh,h,step);
     
          //node 1
     node = &h->nodes[0];
     node->x = i+2*step;
     node->y = j+step;
     node->z = k+step;
    
     //node 2
     node = &h->nodes[1];
     node->x = i+3*step;
     node->y = j+step;
     node->z = k+2*step;
    
     //node 3
     node = &h->nodes[2];
     node->x = i+3*step;
     node->y = j+2*step;
     node->z = k+2*step;
    
    //node 4
    node = &h->nodes[3];
    node->x = i+2*step;
    node->y = j+2*step;
    node->z = k+ step;
    
    
    //node 5
    node = &h->nodes[4];
    node->x = i+2*step;
    node->y = j  ;
    node->z = k+2*step;
    
    //node 6
    node = &h->nodes[5];
    node->x = i+3*step;
    node->y = j    ;
    node->z = k+3*step;
    
    //node 7
    node = &h->nodes[6];
    node->x = i+3*step;
    node->y = j+3*step;
    node->z = k+3*step;
    
    //node 8
    node = &h->nodes[7];
    node->x = i+2*step  ;
    node->y = j+3*step  ;
    node->z = k+2*step;
    

    
    // Elemento 13 
     h = (octant_t*) sc_array_push(&mesh->elements);
     hexa_element_init(h);
     h->x = i+step;
     h->y = j+step;
     h->z = k+2*step;
     h->level = level;
     h->pad = 0;
     //SetElemPML(mesh,h,step);
     
          //node 1
     node = &h->nodes[0];
     node->x = i+step;
     node->y = j;
     node->z = k+2*step;
    
     //node 2
     node = &h->nodes[1];
     node->x = i+2*step;
     node->y = j;
     node->z = k+2*step;
    
     //node 3
     node = &h->nodes[2];
     node->x = i+2*step;
     node->y = j+3*step;
     node->z = k+2*step;
    
    //node 4
    node = &h->nodes[3];
    node->x = i+step  ;
    node->y = j+3*step;
    node->z = k+2*step;
    
    
    //node 5
    node = &h->nodes[4];
    node->x = i;
    node->y = j;
    node->z = k+3*step;
    
    //node 6
    node = &h->nodes[5];
    node->x = i+3*step;
    node->y = j;
    node->z = k+3*step;
    
    //node 7
    node = &h->nodes[6];
    node->x = i+3*step;
    node->y = j+3*step;
    node->z = k+3*step;
    
    //node 8
    node = &h->nodes[7];
    node->x = i;
    node->y = j+3*step;
    node->z = k+3*step;


    
}

void hexa_debug_face_hanging(hexa_tree_t* mesh)
{
    hexa_transition_element(mesh,0,0,0,1,0);
}