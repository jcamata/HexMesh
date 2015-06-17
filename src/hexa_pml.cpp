
#include "hexa.h"
#include "pml.h"

int8_t SetNodePML(hexa_tree_t* tree, octant_node_t* node)
{
    int8_t pad = 0;
    if(node->x == 0)               pad |= PML_X0;
    if(node->x == tree->ncellx)    pad |= PML_X1;
    if(node->y == 0)               pad |= PML_Y0;
    if(node->y == tree->ncelly)    pad |= PML_Y1;
    //if(elem->z == 0)             pad |= PML_Z0;
    if(node->z == tree->ncellz)    pad |= PML_Z1;
    return pad;
}

inline int isX0(int8_t pad)
{
    return ((pad&PML_X0)==PML_X0);
}

inline int isX1(int8_t pad)
{
    return ((pad&PML_X1)==PML_X1);
}

inline int isY0(int8_t pad)
{
    return ((pad&PML_Y0)==PML_Y0);
}

inline int isY1(int8_t pad)
{
    return ((pad&PML_Y1)==PML_Y1);
}

inline int isZ0(int8_t pad)
{
    return ((pad&PML_Z0)==PML_Z0);
}

inline int isZ1(int8_t pad)
{
    return ((pad&PML_Z1)==PML_Z1);
}

void SetElemPML(hexa_tree_t* tree, octant_t *elem, int step)
{

    int8_t pad[8] = {0,0,0,0,0,0,0,0};
    for(int i = 0; i <8; ++i)
    {
        pad[i] = SetNodePML(tree,&elem->nodes[i]);
    }
    
    for(int face= 0; face < 6; ++face)
    {
        int no1 = face_map[face][0];
        int no2 = face_map[face][1];
        int no3 = face_map[face][2];
        int no4 = face_map[face][3];
        
        if(isX0(pad[no1]) && isX0(pad[no2]) && isX0(pad[no3]) && isX0(pad[no4]) )
            elem->pad |= PML_X0;
        
        if(isX1(pad[no1]) && isX1(pad[no2]) && isX1(pad[no3]) && isX1(pad[no4]) )
            elem->pad |= PML_X1;
                
        if(isY0(pad[no1]) && isY0(pad[no2]) && isY0(pad[no3]) && isY0(pad[no4]) )
            elem->pad |= PML_Y0;
        
        if(isY1(pad[no1]) && isY1(pad[no2]) && isY1(pad[no3]) && isY1(pad[no4]) )
            elem->pad |= PML_Y1;
        
        if(isZ0(pad[no1]) && isZ0(pad[no2]) && isZ0(pad[no3]) && isZ0(pad[no4]) )
            elem->pad |= PML_Z0;
        
        if(isZ1(pad[no1]) && isZ1(pad[no2]) && isZ1(pad[no3]) && isZ1(pad[no4]) )
            elem->pad |= PML_Z1;
        
    }
    
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
           for(int i = 0; i < 8; i++) {
               int id = elem_to_pml_map[PML_CORNER_X1Y0Z0][i];
               if(id >=0 )
                   pml_e->nodes[i] = elem->nodes[id];
           }
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

