
#include "hexa.h"
#include "pml.h"

void SetPML(hexa_tree_t* tree, octant_t *elem, int step)
{
    if(elem->x == 0)                    elem->pad |= PML_X0;
    if(elem->x == (tree->ncellx-step))  elem->pad |= PML_X1;
    if(elem->y == 0)                    elem->pad |= PML_Y0;
    if(elem->y == (tree->ncelly-step))  elem->pad |= PML_Y1;
    if(elem->z == 0)                    elem->pad |= PML_Z0;
    if(elem->z == (tree->ncellz-step))  elem->pad |= PML_Z1;
}

inline void SetPMLMask(int8_t* mask, int8_t pad)
{
    
    mask[PML_CORNER_X0Y0Z0] = ((pad&PML_X0==PML_X0) && (pad&PML_Y0==PML_Y0) && (pad&PML_Z0==PML_Z0));
    mask[PML_CORNER_X1Y0Z0] = ((pad&PML_X1==PML_X1) && (pad&PML_Y0==PML_Y0) && (pad&PML_Z0==PML_Z0));
    mask[PML_CORNER_X0Y1Z0] = ((pad&PML_X0==PML_X0) && (pad&PML_Y1==PML_Y1) && (pad&PML_Z0==PML_Z0));
    mask[PML_CORNER_X1Y1Z0] = ((pad&PML_X1==PML_X1) && (pad&PML_Y1==PML_Y1) && (pad&PML_Z0==PML_Z0));
    mask[PML_CORNER_X0Y0Z1] = ((pad&PML_X0==PML_X0) && (pad&PML_Y0==PML_Y0) && (pad&PML_Z1==PML_Z1));
    mask[PML_CORNER_X1Y0Z1] = ((pad&PML_X1==PML_X1) && (pad&PML_Y0==PML_Y0) && (pad&PML_Z1==PML_Z1));
    mask[PML_CORNER_X0Y1Z1] = ((pad&PML_X0==PML_X0) && (pad&PML_Y1==PML_Y1) && (pad&PML_Z1==PML_Z1));
    mask[PML_CORNER_X1Y1Z1] = ((pad&PML_X1==PML_X1) && (pad&PML_Y1==PML_Y1) && (pad&PML_Z1==PML_Z1));
    
    mask[PML_EDGE_Z0_X0] = ((pad&PML_X0==PML_X0) && (pad&PML_Z0==PML_Z0));
    mask[PML_EDGE_Z0_X1] = ((pad&PML_X1==PML_X1) && (pad&PML_Z0==PML_Z0));
    mask[PML_EDGE_Z0_Y0] = ((pad&PML_Y0==PML_Y0) && (pad&PML_Z0==PML_Z0));
    mask[PML_EDGE_Z0_Y1] = ((pad&PML_Y1==PML_Y1) && (pad&PML_Z0==PML_Z0));
    
    mask[PML_EDGE_X0_Y0] = ((pad&PML_X0==PML_X0) && (pad&PML_Y0==PML_Y0));
    mask[PML_EDGE_X0_Y1] = ((pad&PML_X0==PML_X1) && (pad&PML_Y1==PML_Y1));
    mask[PML_EDGE_X1_Y0] = ((pad&PML_X1==PML_X1) && (pad&PML_Y0==PML_Y0));
    mask[PML_EDGE_X1_Y1] = ((pad&PML_X1==PML_X1) && (pad&PML_Y1==PML_Y1));
    
    mask[PML_EDGE_Z1_X0] = ((pad&PML_X0==PML_X0) && (pad&PML_Z1==PML_Z1));
    mask[PML_EDGE_Z1_X1] = ((pad&PML_X1==PML_X1) && (pad&PML_Z1==PML_Z1));
    mask[PML_EDGE_Z1_Y0] = ((pad&PML_Y0==PML_Y0) && (pad&PML_Z1==PML_Z1));
    mask[PML_EDGE_Z1_Y1] = ((pad&PML_Y1==PML_Y1) && (pad&PML_Z1==PML_Z1));
    
    mask[PML_FACE_X0] = (pad&PML_X0==PML_X0);
    mask[PML_FACE_X1] = (pad&PML_X1==PML_X1);
    mask[PML_FACE_Y0] = (pad&PML_Y0==PML_Y0);
    mask[PML_FACE_Y1] = (pad&PML_Y1==PML_Y1);
    mask[PML_FACE_Z0] = (pad&PML_Z0==PML_Z0);
    mask[PML_FACE_Z1] = (pad&PML_Z1==PML_Z1);
    
}

void AddPMLElements(hexa_tree_t* tree)
{
   
    int8_t mask[NPML];
    sc_array_t *elements = &tree->elements;
    int nelem = elements->elem_count;
    for(int i =0 ; i < nelem ; ++i)
    {
        octant_t* elem = (octant_t*) sc_array_index(elements, i);
        SetPMLMask(mask,elem->pad);
        
        if(mask[PML_CORNER_X0Y0Z0]) {
           octant_t* pml_e = (octant_t*) sc_array_push(elements);
           pml_e->pad = PML_CORNER_X0Y0Z0;
           continue;
        }
        
        if(mask[PML_CORNER_X1Y0Z0])
        {
           octant_t* pml_e = (octant_t*) sc_array_push(elements);
           pml_e->pad = PML_CORNER_X1Y0Z0;
           for(int i = 0; i < 8; i++)
           pml_e->nodes[0] = elem->nodes[elem_to_pml_map[PML_CORNER_X1Y0Z0][i]];
           continue;
        }
 
        if(mask[PML_CORNER_X0Y1Z0])
        {
           octant_t* pml_e = (octant_t*) sc_array_push(elements);
           pml_e->pad = PML_CORNER_X0Y1Z0;
           continue;
        }
        
        if(mask[PML_CORNER_X1Y1Z0])
        {
           octant_t* pml_e = (octant_t*) sc_array_push(elements);
           pml_e->pad = PML_CORNER_X1Y1Z0;
           continue;
        }
        
        if(mask[PML_CORNER_X0Y0Z1])
        {
           octant_t* pml_e = (octant_t*) sc_array_push(elements);
           pml_e->pad = PML_CORNER_X0Y0Z1;
           continue;
        }
        
        if(mask[PML_CORNER_X1Y0Z1])
        {
           octant_t* pml_e = (octant_t*) sc_array_push(elements);
           pml_e->pad = PML_CORNER_X1Y0Z1;
           continue;
        }
        
        if(mask[PML_CORNER_X0Y1Z1])
        {
           octant_t* pml_e = (octant_t*) sc_array_push(elements);
           pml_e->pad = PML_CORNER_X0Y1Z1;
           continue;
        }  
        
        if(mask[PML_CORNER_X1Y1Z1])
        {
           octant_t* pml_e = (octant_t*) sc_array_push(elements);
           pml_e->pad = PML_CORNER_X1Y1Z1;
           continue;
        }
        
        if(mask[PML_EDGE_Z0_X0])
        {
           octant_t* pml_e = (octant_t*) sc_array_push(elements);
           pml_e->pad = PML_EDGE_Z0_X0;
           continue;
        }
        
        if(mask[PML_EDGE_Z0_X1])
        {
           octant_t* pml_e = (octant_t*) sc_array_push(elements);
           pml_e->pad = PML_EDGE_Z0_X1;
           continue;
        }
        
        if(mask[PML_EDGE_Z0_Y0])
        {
           octant_t* pml_e = (octant_t*) sc_array_push(elements);
           pml_e->pad = PML_EDGE_Z0_Y0;
           continue;
        }
        
        if(mask[PML_EDGE_Z0_Y1])
        {
           octant_t* pml_e = (octant_t*) sc_array_push(elements);
           pml_e->pad = PML_EDGE_Z0_Y1;
           continue;
        }        
              
        if(mask[PML_EDGE_X0_Y0])
        {
           octant_t* pml_e = (octant_t*) sc_array_push(elements);
           pml_e->pad = PML_EDGE_X0_Y0;
           continue;
        }
        
        if(mask[PML_EDGE_X0_Y1])
        {
           octant_t* pml_e = (octant_t*) sc_array_push(elements);
           pml_e->pad = PML_EDGE_X0_Y1;
           continue;
        }
        
        if(mask[PML_EDGE_X1_Y0])
        {
           octant_t* pml_e = (octant_t*) sc_array_push(elements);
           pml_e->pad = PML_EDGE_X1_Y0;
           continue;
        }
        
        if(mask[PML_EDGE_X1_Y1])
        {
           octant_t* pml_e = (octant_t*) sc_array_push(elements);
           pml_e->pad = PML_EDGE_X1_Y0;
           continue;
        }        
        
        if(mask[PML_EDGE_Z1_X0])
        {
           octant_t* pml_e = (octant_t*) sc_array_push(elements);
           pml_e->pad = PML_EDGE_Z1_X0;
           continue;
        }  
        
        if(mask[PML_EDGE_Z1_X1])
        {
           octant_t* pml_e = (octant_t*) sc_array_push(elements);
           pml_e->pad = PML_EDGE_Z1_X1;
           continue;
        } 
        
        if(mask[PML_EDGE_Z1_Y0])
        {
           octant_t* pml_e = (octant_t*) sc_array_push(elements);
           pml_e->pad = PML_EDGE_Z1_Y0;
           continue;
        }
        
        if(mask[PML_EDGE_Z1_Y1])
        {
           octant_t* pml_e = (octant_t*) sc_array_push(elements);
           pml_e->pad = PML_EDGE_Z1_Y1;
           continue;
        }
        
        if(mask[PML_FACE_X0])
        {
           octant_t* pml_e = (octant_t*) sc_array_push(elements);
           pml_e->pad = PML_FACE_X0;
           continue;
        }
        
        if(mask[PML_FACE_X1])
        {
           octant_t* pml_e = (octant_t*) sc_array_push(elements);
           pml_e->pad = PML_FACE_X1;
           continue;
        }       
        
        if(mask[PML_FACE_Y0])
        {
           octant_t* pml_e = (octant_t*) sc_array_push(elements);
           pml_e->pad = PML_FACE_Y0;
           continue;
        } 
        
        if(mask[PML_FACE_Y1])
        {
           octant_t* pml_e = (octant_t*) sc_array_push(elements);
           pml_e->pad = PML_FACE_Y1;
           continue;
        }
        
        if(mask[PML_FACE_Z0])
        {
           octant_t* pml_e = (octant_t*) sc_array_push(elements);
           pml_e->pad = PML_FACE_Z0;
           continue;
        } 
        
        if(mask[PML_FACE_Z1])
        {
           octant_t* pml_e = (octant_t*) sc_array_push(elements);
           pml_e->pad = PML_FACE_Z1;
           continue;
        }  
        
    }
    
}




