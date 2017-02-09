/* 
 * File:   pml.h
 * Author: camata
 *
 * Created on April 13, 2015, 11:30 AM
 */

#ifndef PML_H
#define	PML_H

#define NPML 27

// Corner PML
#define PML_CORNER_X0Y0Z0 0
#define PML_CORNER_X1Y0Z0 1
#define PML_CORNER_X0Y1Z0 2
#define PML_CORNER_X1Y1Z0 3
#define PML_CORNER_X0Y0Z1 4
#define PML_CORNER_X1Y0Z1 5
#define PML_CORNER_X0Y1Z1 6
#define PML_CORNER_X1Y1Z1 7
// Edge PMLs
//
#define PML_EDGE_Z0_X0         8
#define PML_EDGE_Z0_X1         9
#define PML_EDGE_Z0_Y0         10
#define PML_EDGE_Z0_Y1         11
#define PML_EDGE_X0_Y0         12
#define PML_EDGE_X0_Y1         13
#define PML_EDGE_X1_Y0         14
#define PML_EDGE_X1_Y1         15
#define PML_EDGE_Z1_X0         16
#define PML_EDGE_Z1_X1         17
#define PML_EDGE_Z1_Y0         18
#define PML_EDGE_Z1_Y1         19
// Face PMLS
#define PML_FACE_X0  20
#define PML_FACE_X1  21
#define PML_FACE_Y0  22
#define PML_FACE_Y1  23
#define PML_FACE_Z0  24
#define PML_FACE_Z1  25

#define XM (mesh->ncellx+2)
#define XP (mesh->ncellx+3)

#define YM (mesh->ncelly+2)
#define YP (mesh->ncelly+3)

#define ZM (mesh->ncellz+2)
#define ZP (mesh->ncellz+3)

typedef enum {
    PML_NULL = 0,
    PML_X0   = 1 << 1,
    PML_X1   = 1 << 2,
    PML_Y0   = 1 << 3,
    PML_Y1   = 1 << 4,
    PML_Z0   = 1 << 5,
    PML_Z1   = 1 << 6
} pml_mask_t;

static int face_map[6][4] = 
{
    {0,4,7,3}, // 
    {1,2,6,5},
    {0,1,5,4},
    {3,2,6,7},
    {0,1,2,3},
    {4,5,6,7}
};


static int elem_to_pml_map[NPML][8] = {
                          /* 0    1   2   3   4   5   6   7   */
 /* PML_CORNER_X0Y0Z0 */   {-1 , -1, -1, -1, -1, -1, 0, -1},
 /* PML_CORNER_X1Y0Z0 */   {-1 , -1, -1, -1, -1, -1, -1, -1},
 /* PML_CORNER_X0Y1Z0 */   {-1 , -1, -1,  3, -1, -1, -1, -1},
 /* PML_CORNER_X1Y1Z0 */   {-1 , -1,  2, -1, -1, -1, -1, -1},
 /* PML_CORNER_X0Y0Z1 */   {-1 , -1, -1, -1,  4, -1, -1, -1},  
 /* PML_CORNER_X1Y0Z1 */   {-1 , -1, -1, -1, -1,  5, -1, -1},
 /* PML_CORNER_X0Y1Z1 */   {-1 , -1, -1, -1, -1, -1, -1,  7},
 /* PML_CORNER_X1Y1Z1 */   {-1 , -1, -1, -1, -1, -1,  6, -1},
                           
};

void SetElemPML(hexa_tree_t* tree, octant_t *elem, int step);

#endif	/* PML_H */

