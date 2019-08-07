
#include <stdlib.h>
#include <stdio.h>
#include <sc.h>
#include <sc_io.h>
#include <sc_containers.h>

#include "hexa.h"
#include "pml.h"

void hexa_element_init(octant_t *elem);


void hexa_transition_element(hexa_tree_t* mesh, int i, int j, int k, int step, int level, int ext)
{
	octant_node_t *node;
	octant_t * h;
	// Creating 13 news hexahedra

	//Elemento 1
	if(true){
		h = (octant_t*) sc_array_push(&mesh->elements);
		hexa_element_init(h);
		h->x = i;
		h->y = j;
		h->z = k;
		h->level = level;
		if(ext == 1) h->boundary = true;
		if(ext == 3) h->boundary = true;

		if(ext == 5) h->boundary = true;
		if(ext == 6) h->boundary = true;
		if(ext == 8) h->boundary = true;
		h->pad = 0;
		h->tem = 5;
		//SetElemPML(mesh,h,step);
		//----------------------------

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
	}

	//Elemento 2
	if(true){
		h = (octant_t*) sc_array_push(&mesh->elements);
		hexa_element_init(h);
		h->x = i+step;
		h->y = j;
		h->z = k;
		h->level = level;
		h->pad = 0;
		h->tem = 5;
		if(ext == 1) h->boundary = true;

		if(ext == 5) h->boundary = true;
		if(ext == 6) h->boundary = true;
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
		node->x = i+step;
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
	}

	//Elemento 3
	if(true){
		h = (octant_t*) sc_array_push(&mesh->elements);
		hexa_element_init(h);
		h->x = i+2*step;
		h->y = j;
		h->z = k;
		h->level = level;
		h->pad = 0;
		h->tem = 5;
		if(ext == 1) h->boundary = true;
		if(ext == 4) h->boundary = true;

		if(ext == 5) h->boundary = true;
		if(ext == 6) h->boundary = true;
		if(ext == 7) h->boundary = true;
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
		node->x = i+2*step;
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
	}

	//Elemento 4
	if(true){
		h = (octant_t*) sc_array_push(&mesh->elements);
		hexa_element_init(h);
		h->x = i;
		h->y = j+step;
		h->z = k;
		h->level = level;
		h->pad = 0;
		h->tem = 5;
		if(ext == 3) h->boundary = true;

		if(ext == 5) h->boundary = true;
		if(ext == 8) h->boundary = true;
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
		node->x = i;
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
	}

	//Elemento 5 //central
	if(true){
		h = (octant_t*) sc_array_push(&mesh->elements);
		hexa_element_init(h);
		h->x = i+step;
		h->y = j+step;
		h->z = k;
		h->level = level;
		h->pad = 0;
		h->tem = 5;
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
		node->x = i+step;
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
	}

	//Elemento 6
	if(true){
		h = (octant_t*) sc_array_push(&mesh->elements);
		hexa_element_init(h);
		h->x = i+2*step;
		h->y = j+step;
		h->z = k;
		h->level = level;
		h->pad = 0;
		h->tem = 5;
		if(ext == 4) h->boundary = true;

		if(ext == 6) h->boundary = true;
		if(ext == 7) h->boundary = true;
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
		node->x = i+2*step;
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
	}

	//Elemento 7
	if(true){
		h = (octant_t*) sc_array_push(&mesh->elements);
		hexa_element_init(h);
		h->x = i;
		h->y = j+2*step;
		h->z = k;
		h->level = level;
		h->pad = 0;
		h->tem = 5;
		if(ext == 3) h->boundary = true;
		if(ext == 2) h->boundary = true;

		if(ext == 5) h->boundary = true;
		if(ext == 7) h->boundary = true;
		if(ext == 8) h->boundary = true;
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
	}

	//Elemento 8
	if(true){
		h = (octant_t*) sc_array_push(&mesh->elements);
		hexa_element_init(h);
		h->x = i+step;
		h->y = j+2*step;
		h->z = k;
		h->level = level;
		h->pad = 0;
		h->tem = 5;
		if(ext == 2) h->boundary = true;

		if(ext == 7) h->boundary = true;
		if(ext == 8) h->boundary = true;
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
	}

	//Elemento 9
	if(true){
		h = (octant_t*) sc_array_push(&mesh->elements);
		hexa_element_init(h);
		h->x = i+2*step;
		h->y = j+2*step;
		h->z = k;
		h->level = level;
		h->pad = 0;
		h->tem = 5;
		if(ext == 2) h->boundary = true;
		if(ext == 4) h->boundary = true;

		if(ext == 6) h->boundary = true;
		if(ext == 7) h->boundary = true;
		if(ext == 8) h->boundary = true;
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
	}

	//Elemento 10 //bottom center
	if(true){
		h = (octant_t*) sc_array_push(&mesh->elements);
		hexa_element_init(h);
		h->x = i+step;
		h->y = j+step;
		h->z = k+step;
		h->level = level;
		h->pad = 0;
		h->tem = 5;
		//h->boundary = true;
		//if(ext == 1) h->boundary = true;
		//if(ext == 2) h->boundary = true;
		//if(ext == 3) h->boundary = true;
		//if(ext == 4) h->boundary = true;
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
		node->x = i+step;
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
	}

	//Elemento 11
	if(true){
		h = (octant_t*) sc_array_push(&mesh->elements);
		hexa_element_init(h);
		h->x = i;
		h->y = j+step;
		h->z = k+2*step;
		h->level = level;
		h->pad = 0;
		h->tem = 5;
		if(ext==3) h->boundary = true;
		if(ext==5) h->boundary = true;
		if(ext==8) h->boundary = true;
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
	}

	//Elemento 12
	if(true){
		h = (octant_t*) sc_array_push(&mesh->elements);
		hexa_element_init(h);
		h->x = i+2*step;
		h->y = j+step;
		h->z = k+step;
		h->level = level;
		h->pad = 0;
		h->tem = 5;
		if(ext==4) h->boundary = true;
		if(ext==6) h->boundary = true;
		if(ext==7) h->boundary = true;
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
	}

	//Elemento 13
	if(true){
		h = (octant_t*) sc_array_push(&mesh->elements);
		hexa_element_init(h);
		h->x = i+step;
		h->y = j;//+step
		h->z = k+2*step;
		h->level = level;
		h->pad = 0;
		h->tem = 5;
		if(ext==2 || ext == 1) h->boundary = true;
		if(ext==5) h->boundary = true;
		if(ext==6) h->boundary = true;
		if(ext==7) h->boundary = true;
		if(ext==8) h->boundary = true;
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

}
