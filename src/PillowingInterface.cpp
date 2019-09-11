#include <gts.h>
#include <glib.h>
#include <vector>
#include <iostream>
using namespace std;
#include <set>
#include <algorithm>
#include <sc.h>
#include <sc_io.h>
#include <sc_containers.h>
#include <mpi.h>
#include <cassert>

#include "hexa.h"
#include "hilbert.h"

#include <ctime>
#include <chrono>

typedef struct pillow
{
	uint64_t id;
	uint64_t a;
	uint64_t b;
	//	uint64_t elem_a;
	//	uint64_t elem_b;
	int32_t x,y,z;
	int8_t  list_elem;
	int8_t  list_face[8];
	int32_t  elem[8]; // we must check this values
	int8_t   face[8][3]; // we must check this values
} pillow_t;

unsigned pillow_hash_fn (const void *v, const void *u)
{
	const pillow_t *q = (const pillow_t*) v;
	uint32_t            a, b, c;

	a = (uint32_t) q->x;
	b = (uint32_t) q->y;
	c = (uint32_t) q->z;
	sc_hash_mix (a, b, c);
	sc_hash_final (a, b, c);
	return (unsigned) c;
}

int pillow_equal_fn (const void *v1, const void *v2, const void *u)
{
	const pillow_t *q1 = (const pillow_t*) v1;
	const pillow_t *q2 = (const pillow_t*) v2;
	return (q1->x == q2->x && q1->y == q2->y && q1->z == q2->z);
}

int AddPoint(hexa_tree_t* mesh, sc_hash_array_t* hash_nodes, GtsPoint *p, std::vector<double> &coords, int x, int y, int z) {
	size_t position;
	octant_node_t *r;
	octant_node_t key;

	key.x = x;
	key.y = y;
	key.z = z;

	r = (octant_node_t*) sc_hash_array_insert_unique(hash_nodes, &key, &position);

	if (r != NULL) {
		r->x = x;
		r->y = y;
		r->z = z;
		r->id = mesh->nodes.elem_count;
		octant_node_t* n = (octant_node_t*) sc_array_push(&mesh->nodes);
		n->id = r->id;
		n->x = x;
		n->y = y;
		n->z = z;
		n->color = -1;
		n->fixed = 0;

		const double xx = p->x;
		const double yy = p->y;
		const double zz = p->z;

		coords.push_back(xx);
		coords.push_back(yy);
		coords.push_back(zz);
		return r->id;
	} else {
		r = (octant_node_t*) sc_array_index(&hash_nodes->a, position);
		return r->id;
	}
}

GtsPoint* LinearMapHex(const double* cord_in_ref, const double* cord_in_x, const double* cord_in_y, const double* cord_in_z){

	double N[8];
	GtsPoint* point;
	double out[3];


	N[0] = (1-cord_in_ref[0])*(1-cord_in_ref[1])*(1-cord_in_ref[2])/double(8);
	N[1] = (1+cord_in_ref[0])*(1-cord_in_ref[1])*(1-cord_in_ref[2])/double(8);
	N[2] = (1+cord_in_ref[0])*(1+cord_in_ref[1])*(1-cord_in_ref[2])/double(8);
	N[3] = (1-cord_in_ref[0])*(1+cord_in_ref[1])*(1-cord_in_ref[2])/double(8);

	N[4] = (1-cord_in_ref[0])*(1-cord_in_ref[1])*(1+cord_in_ref[2])/double(8);
	N[5] = (1+cord_in_ref[0])*(1-cord_in_ref[1])*(1+cord_in_ref[2])/double(8);
	N[6] = (1+cord_in_ref[0])*(1+cord_in_ref[1])*(1+cord_in_ref[2])/double(8);
	N[7] = (1-cord_in_ref[0])*(1+cord_in_ref[1])*(1+cord_in_ref[2])/double(8);

	out[0] = 0;
	out[1] = 0;
	out[2] = 0;

	for(int i=0;i<8;i++){
		out[0] = N[i]*cord_in_x[i] + out[0] ;
		out[1] = N[i]*cord_in_y[i] + out[1];
		out[2] = N[i]*cord_in_z[i] + out[2];
	}

	point = gts_point_new(gts_point_class(),out[0],out[1],out[2]);

	return point;
}

void RedoNodeMapping(hexa_tree_t* mesh)
{

	int factor = 12;
	//just the multiplication
	// it allow us add int points in the mesh
	// keeping a structured mesh
	for(int iel = 0; iel <mesh->elements.elem_count; iel++)
	{
		octant_t * elem = (octant_t*) sc_array_index (&mesh->elements, iel);
		for(int ino = 0; ino < 8; ino++)
		{
			elem->nodes[ino].x = factor*elem->nodes[ino].x;
			elem->nodes[ino].y = factor*elem->nodes[ino].y;
			elem->nodes[ino].z = factor*elem->nodes[ino].z;
		}
		elem->x = 4*elem->x;
		elem->y = 4*elem->y;
		elem->z = 4*elem->z;
	}

	for(int ino = 0; ino <mesh->nodes.elem_count; ino++)
	{
		octant_node_t * node = (octant_node_t*) sc_array_index (&mesh->nodes, ino);
		node->x = factor*node->x;
		node->y = factor*node->y;
		node->z = factor*node->z;
	}
	mesh->x_start = factor*mesh->x_start;
	mesh->y_start = factor*mesh->y_start;
	mesh->x_end = factor*mesh->x_end;
	mesh->y_end = factor*mesh->y_end;

	mesh->ncellx = 4*mesh->ncellx;
	mesh->ncelly = 4*mesh->ncelly;
	mesh->max_z = 4*mesh->max_z;
}


void CopyPropEl(hexa_tree_t* mesh, int id, octant_t *elem1)
{

	octant_t *elem = (octant_t*) sc_array_index(&mesh->elements, id);

	elem1->level = -1;
	elem1->tem = elem->tem;
	elem1->pad = elem->pad;
	elem1->n_mat = elem->n_mat;
	elem1->pml_id = elem->pml_id;
	elem1->father = id;
	elem1->boundary = elem->boundary;

	for(int ino = 0; ino < 8; ino++)
	{
		elem1->nodes[ino].color = elem->nodes[ino].color;
		elem1->nodes[ino].fixed = elem->nodes[ino].fixed;
		elem1->nodes[ino].x = elem->nodes[ino].x;
		elem1->nodes[ino].y = elem->nodes[ino].y;
		elem1->nodes[ino].z = elem->nodes[ino].z;
	}

	for(int iedge = 0; iedge < 12; iedge++)
	{
		elem1->edge[iedge].coord[0] = elem->edge[iedge].coord[0];
		elem1->edge[iedge].coord[1] = elem->edge[iedge].coord[1];
		elem1->edge[iedge].id = elem->edge[iedge].id;
		elem1->edge[iedge].ref = elem->edge[iedge].ref;
	}

	for(int isurf = 0; isurf < 6; isurf++)
	{
		elem1->surf[isurf].ext = elem->surf[isurf].ext;
	}

	elem1->x=elem->x;
	elem1->y=elem->y;
	elem1->z=elem->z;
}

void SurfaceIdentification(hexa_tree_t* mesh, std::vector<double>& coords){

	//assign a color for the nodes...

	// free node
	// color= 0 && fixed = 0;

	//exterior surface fixed nodes
	//fixed = 1 && color = 1

	//exterior global nodes
	//color && fixed = -1;
	// x- = 2 x+ = 3
	// y- = 4 y+ = 5
	// z- = 6 z+ = 7

	//exterior local nodes
	//color && fixed = -1;
	// x- = -2 x+ = -3
	// y- = -4 y+ = -5
	// z- = -6 z+ = -7

	//exterior global edges
	//color && fixed = -1;
	// z-y- = 0+11
	// z-x+ = 1+11
	// z-y+ = 2+11
	// z-x- = 3+11

	// x-y- = 4+11
	// x+y- = 5+11
	// x+y+ = 6+11
	// x-y+ = 7+11

	// z+y- = 8+11
	// z+x+ = 9+11
	// z+y+ = 10+11
	// z+x- = 11+11

	//exterior local edges
	//color && fixed = -1;
	// z-y- = -0-11
	// z-x+ = -1-11
	// z-y+ = -2-11
	// z-x- = -3-11

	// x-y- = -4-11
	// x+y- = -5-11
	// x+y+ = -6-11
	// x-y+ = -7-11

	// z+y- = -8-11
	// z+x+ = -9-11
	// z+y+ = -10-11
	// z+x- = -11-11

	bool deb = false;
	bool clamped = true;
	//TODO remove vertex hash
	//vertex hash
	sc_hash_array_t*	  vertex_hash  = (sc_hash_array_t *)sc_hash_array_new(sizeof (octant_vertex_t), vertex_hash_id, vertex_equal_id, &clamped);

	//fazendo vertex hash & assign the color to the local nodes
	for(int iel = 0; iel < mesh->elements.elem_count ; iel++){
		size_t  position;
		octant_t *elem = (octant_t*) sc_array_index(&mesh->elements, iel);
		for (int ino = 0; ino < 8; ino++){
			//build the hash
			octant_vertex_t key;
			key.id = elem->nodes[ino].id;
			octant_vertex_t* vert = (octant_vertex_t*) sc_hash_array_insert_unique (vertex_hash, &key, &position);
			if(vert != NULL){
				vert->id = elem->nodes[ino].id;
				vert->list_elem = 1;
				vert->elem[vert->list_elem-1] = elem->id;
			}else{
				vert = (octant_vertex_t*) sc_array_index(&vertex_hash->a, position);
				vert->elem[vert->list_elem] = elem->id;
				vert->list_elem++;
			}

			//setting free the free nodes
			if(elem->nodes[ino].fixed == 0) elem->nodes[ino].color = 0;

			//assign the color for the local nodes...

			//surface and edges
			if(elem->nodes[ino].x == mesh->x_start)
			{
				elem->nodes[ino].fixed = -1;
				elem->nodes[ino].color = -2;
			}

			if(elem->nodes[ino].x == mesh->x_end)
			{
				elem->nodes[ino].fixed = -1;
				elem->nodes[ino].color = -3;
			}

			if(elem->nodes[ino].y == mesh->y_start)
			{
				elem->nodes[ino].fixed = -1;
				elem->nodes[ino].color = -4;
			}

			if(elem->nodes[ino].y == mesh->y_end)
			{
				elem->nodes[ino].fixed = -1;
				elem->nodes[ino].color = -5;
			}

			if(elem->nodes[ino].z == 0)
			{
				elem->nodes[ino].fixed = -1;
				elem->nodes[ino].color = -6;
			}

			if(elem->nodes[ino].z == 3*mesh->max_z)
			{
				elem->nodes[ino].fixed = -1;
				elem->nodes[ino].color = -7;
			}

			// z-y- = -0-11
			if(elem->nodes[ino].z == 0 && elem->nodes[ino].y == mesh->y_start)
			{
				elem->nodes[ino].fixed = -1;
				elem->nodes[ino].color = -0-11;
			}
			// z-x+ = -1-11
			if(elem->nodes[ino].z == 0 && elem->nodes[ino].x == mesh->x_end)
			{
				elem->nodes[ino].fixed = -1;
				elem->nodes[ino].color = -1-11;
			}
			// z-y+ = -2-11
			if(elem->nodes[ino].z == 0 && elem->nodes[ino].y == mesh->y_end)
			{
				elem->nodes[ino].fixed = -1;
				elem->nodes[ino].color = -2-11;
			}
			// z-x- = -3-11
			if(elem->nodes[ino].z == 0 && elem->nodes[ino].x == mesh->x_start)
			{
				elem->nodes[ino].fixed = -1;
				elem->nodes[ino].color = -3-11;
			}


			// x-y- = -4-11
			if(elem->nodes[ino].x == mesh->x_start && elem->nodes[ino].y == mesh->y_start)
			{
				elem->nodes[ino].fixed = -1;
				elem->nodes[ino].color = -4-11;
			}
			// x+y- = -5-11
			if(elem->nodes[ino].x == mesh->x_end && elem->nodes[ino].y == mesh->y_start)
			{
				elem->nodes[ino].fixed = -1;
				elem->nodes[ino].color = -5-11;
			}
			// x+y+ = -6-11
			if(elem->nodes[ino].x == mesh->x_end && elem->nodes[ino].y == mesh->y_end)
			{
				elem->nodes[ino].fixed = -1;
				elem->nodes[ino].color = -6-11;
			}
			// x-y+ = -7-11
			if(elem->nodes[ino].x == mesh->x_start && elem->nodes[ino].y == mesh->y_end)
			{
				elem->nodes[ino].fixed = -1;
				elem->nodes[ino].color = -7-11;
			}

			// z+y- = -8-11
			if(elem->nodes[ino].z == 3*mesh->max_z && elem->nodes[ino].y == mesh->y_start)
			{
				elem->nodes[ino].fixed = -1;
				elem->nodes[ino].color = -8-11;
			}
			// z+x+ = -9-11
			if(elem->nodes[ino].z == 3*mesh->max_z && elem->nodes[ino].x == mesh->x_end)
			{
				elem->nodes[ino].fixed = -1;
				elem->nodes[ino].color = -9-11;
			}
			// z+y+ = -10-11
			if(elem->nodes[ino].z == 3*mesh->max_z && elem->nodes[ino].y == mesh->y_end)
			{
				elem->nodes[ino].fixed = -1;
				elem->nodes[ino].color = -10-11;
			}
			// z+x- = -11-11
			if(elem->nodes[ino].z == 3*mesh->max_z && elem->nodes[ino].x == mesh->x_start)
			{
				elem->nodes[ino].fixed = -1;
				elem->nodes[ino].color = -11-11;
			}
		}
	}

	sc_array_init(&mesh->outsurf, sizeof(octant_t));

	sc_array_t toto;
	sc_array_init(&toto, sizeof(octant_t));

	//id global exterior surface
	for(int iel = 0; iel < mesh->elements.elem_count; iel++){
		octant_t* elemOrig = (octant_t*) sc_array_index(&mesh->elements, iel);
		octant_t * elem = (octant_t*) sc_array_push(&toto);
		hexa_element_copy(elemOrig,elem);

		if(elem->boundary){

			if(true){
				int isurf;
				isurf = 0;
				elem->surf[isurf].ext = false;
				if(elem->nodes[FaceNodesMap[isurf][0]].x == 0 && elem->nodes[FaceNodesMap[isurf][1]].x == 0 &&
						elem->nodes[FaceNodesMap[isurf][2]].x == 0 && elem->nodes[FaceNodesMap[isurf][3]].x == 0){
					elem->surf[isurf].ext = true;
					if(deb) for(int ino = 0; ino < 4; ino++) mesh->part_nodes[elem->nodes[FaceNodesMap[isurf][ino]].id] = isurf+20;
				}

				isurf = 1;
				elem->surf[isurf].ext = false;
				if(elem->nodes[FaceNodesMap[isurf][0]].x == 3*mesh->ncellx && elem->nodes[FaceNodesMap[isurf][1]].x == 3*mesh->ncellx &&
						elem->nodes[FaceNodesMap[isurf][2]].x == 3*mesh->ncellx && elem->nodes[FaceNodesMap[isurf][3]].x == 3*mesh->ncellx){
					elem->surf[isurf].ext = true;
					if(deb) for(int ino = 0; ino < 4; ino++) mesh->part_nodes[elem->nodes[FaceNodesMap[isurf][ino]].id] = isurf+20;
				}

				isurf = 2;
				elem->surf[isurf].ext = false;
				if(elem->nodes[FaceNodesMap[isurf][0]].y == 0 && elem->nodes[FaceNodesMap[isurf][1]].y == 0 &&
						elem->nodes[FaceNodesMap[isurf][2]].y == 0 && elem->nodes[FaceNodesMap[isurf][3]].y == 0){
					elem->surf[isurf].ext = true;
					if(deb) for(int ino = 0; ino < 4; ino++) mesh->part_nodes[elem->nodes[FaceNodesMap[isurf][ino]].id] = isurf+20;
				}

				isurf = 3;
				elem->surf[isurf].ext = false;
				if(elem->nodes[FaceNodesMap[isurf][0]].y == 3*mesh->ncelly && elem->nodes[FaceNodesMap[isurf][1]].y == 3*mesh->ncelly &&
						elem->nodes[FaceNodesMap[isurf][2]].y == 3*mesh->ncelly && elem->nodes[FaceNodesMap[isurf][3]].y == 3*mesh->ncelly){
					elem->surf[isurf].ext = true;
					if(deb) for(int ino = 0; ino < 4; ino++) mesh->part_nodes[elem->nodes[FaceNodesMap[isurf][ino]].id] = isurf+20;
				}

				isurf = 4;
				elem->surf[isurf].ext = false;
				if(elem->nodes[FaceNodesMap[isurf][0]].z == 3*mesh->max_z && elem->nodes[FaceNodesMap[isurf][1]].z == 3*mesh->max_z &&
						elem->nodes[FaceNodesMap[isurf][2]].z == 3*mesh->max_z && elem->nodes[FaceNodesMap[isurf][3]].z == 3*mesh->max_z){
					elem->surf[isurf].ext = true;
					if(deb) for(int ino = 0; ino < 4; ino++) mesh->part_nodes[elem->nodes[FaceNodesMap[isurf][ino]].id] = isurf+20;
				}

				isurf = 5;
				elem->surf[isurf].ext = false;
				if(elem->nodes[FaceNodesMap[isurf][0]].z == 0 && elem->nodes[FaceNodesMap[isurf][1]].z == 0 &&
						elem->nodes[FaceNodesMap[isurf][2]].z == 0 && elem->nodes[FaceNodesMap[isurf][3]].z == 0){
					elem->surf[isurf].ext = true;
					if(deb) for(int ino = 0; ino < 4; ino++) mesh->part_nodes[elem->nodes[FaceNodesMap[isurf][ino]].id] = isurf+20;
				}
			}

			bool aux = false;
			for(int isurf = 0; isurf < 6; isurf++){
				if(elem->surf[isurf].ext) aux = true;
			}
			if(aux){
				octant_t* elem1 = (octant_t*) sc_array_push(&mesh->outsurf);
				//hexa_element_init(elem1);

				elem1->level = -1;
				elem1->id = elem->id;
				elem1->tem = elem->tem;
				elem1->pad = elem->pad;
				elem1->n_mat = elem->n_mat;
				elem1->pml_id = elem->pml_id;
				elem1->father = elem->id;
				elem1->boundary = elem->boundary;
				elem1->x=elem->x;
				elem1->y=elem->y;
				elem1->z=elem->z;

				for(int ino = 0; ino < 8; ino++){
					elem1->nodes[ino].color = -elem->nodes[ino].color;
					elem1->nodes[ino].fixed = 0;
					if(elem->nodes[ino].fixed == 1)
					{
						elem1->nodes[ino].fixed = 1;
					}else
					{
						elem1->nodes[ino].fixed = -elem->nodes[ino].fixed;
					}
					elem1->nodes[ino].id = elem->nodes[ino].id;
					elem1->nodes[ino].x = elem->nodes[ino].x;
					elem1->nodes[ino].y = elem->nodes[ino].y;
					elem1->nodes[ino].z = elem->nodes[ino].z;
				}

				for(int iedge = 0; iedge < 12; iedge++){
					elem1->edge[iedge].coord[0] = elem->edge[iedge].coord[0];
					elem1->edge[iedge].coord[1] = elem->edge[iedge].coord[1];
					elem1->edge[iedge].id = elem->edge[iedge].id;
					elem1->edge[iedge].ref = false;
				}

				for(int isurf = 0; isurf < 6; isurf++){
					elem1->surf[isurf].ext = elem->surf[isurf].ext;
				}
			}
		}

		sc_array_reset(&toto);
	}

	//id global exterior edges
	for(int iel = 0; iel < mesh->outsurf.elem_count; iel++){
		octant_t* elem = (octant_t*) sc_array_index (&mesh->outsurf, iel);

		//edge 0
		int iedge;
		iedge = 0;
		if(elem->nodes[EdgeVerticesMap[iedge][0]].y == 0 && elem->nodes[EdgeVerticesMap[iedge][1]].y == 0){
			if(elem->nodes[EdgeVerticesMap[iedge][0]].z == 0 && elem->nodes[EdgeVerticesMap[iedge][1]].z == 0){
				elem->edge[iedge].ref = true;
				if(deb) mesh->part_nodes[elem->nodes[EdgeVerticesMap[iedge][0]].id] = iedge;
				if(deb) mesh->part_nodes[elem->nodes[EdgeVerticesMap[iedge][1]].id] = iedge;
			}
		}

		iedge = 1;
		if(elem->nodes[EdgeVerticesMap[iedge][0]].x == 3*mesh->ncellx && elem->nodes[EdgeVerticesMap[iedge][1]].x == 3*mesh->ncellx){
			if(elem->nodes[EdgeVerticesMap[iedge][0]].z == 0 && elem->nodes[EdgeVerticesMap[iedge][1]].z == 0){
				elem->edge[iedge].ref = true;
				if(deb) mesh->part_nodes[elem->nodes[EdgeVerticesMap[iedge][0]].id] = iedge;
				if(deb) mesh->part_nodes[elem->nodes[EdgeVerticesMap[iedge][1]].id] = iedge;
			}
		}

		iedge = 2;
		if(elem->nodes[EdgeVerticesMap[iedge][0]].y == 3*mesh->ncelly && elem->nodes[EdgeVerticesMap[iedge][1]].y == 3*mesh->ncelly){
			if(elem->nodes[EdgeVerticesMap[iedge][0]].z == 0 && elem->nodes[EdgeVerticesMap[iedge][1]].z == 0){
				elem->edge[iedge].ref = true;
				if(deb) mesh->part_nodes[elem->nodes[EdgeVerticesMap[iedge][0]].id] = iedge;
				if(deb) mesh->part_nodes[elem->nodes[EdgeVerticesMap[iedge][1]].id] = iedge;
			}
		}

		iedge = 3;
		if(elem->nodes[EdgeVerticesMap[iedge][0]].x == 0 && elem->nodes[EdgeVerticesMap[iedge][1]].x == 0){
			if(elem->nodes[EdgeVerticesMap[iedge][0]].z == 0 && elem->nodes[EdgeVerticesMap[iedge][1]].z == 0){
				elem->edge[iedge].ref = true;
				if(deb) mesh->part_nodes[elem->nodes[EdgeVerticesMap[iedge][0]].id] = iedge;
				if(deb) mesh->part_nodes[elem->nodes[EdgeVerticesMap[iedge][1]].id] = iedge;
			}
		}

		iedge = 4;
		if(elem->nodes[EdgeVerticesMap[iedge][0]].x == 0 && elem->nodes[EdgeVerticesMap[iedge][1]].x == 0){
			if(elem->nodes[EdgeVerticesMap[iedge][0]].y == 0 && elem->nodes[EdgeVerticesMap[iedge][1]].y == 0){
				elem->edge[iedge].ref = true;
				if(deb) mesh->part_nodes[elem->nodes[EdgeVerticesMap[iedge][0]].id] = iedge;
				if(deb) mesh->part_nodes[elem->nodes[EdgeVerticesMap[iedge][1]].id] = iedge;
			}
		}

		iedge = 5;
		if(elem->nodes[EdgeVerticesMap[iedge][0]].x == 3*mesh->ncellx && elem->nodes[EdgeVerticesMap[iedge][1]].x == 3*mesh->ncellx){
			if(elem->nodes[EdgeVerticesMap[iedge][0]].y == 0 && elem->nodes[EdgeVerticesMap[iedge][1]].y == 0){
				elem->edge[iedge].ref = true;
				if(deb) mesh->part_nodes[elem->nodes[EdgeVerticesMap[iedge][0]].id] = iedge;
				if(deb) mesh->part_nodes[elem->nodes[EdgeVerticesMap[iedge][1]].id] = iedge;
			}
		}

		iedge = 6;
		if(elem->nodes[EdgeVerticesMap[iedge][0]].y == 3*mesh->ncelly && elem->nodes[EdgeVerticesMap[iedge][1]].y == 3*mesh->ncelly){
			if(elem->nodes[EdgeVerticesMap[iedge][0]].x == 3*mesh->ncellx && elem->nodes[EdgeVerticesMap[iedge][1]].x == 3*mesh->ncellx){
				elem->edge[iedge].ref = true;
				if(deb) mesh->part_nodes[elem->nodes[EdgeVerticesMap[iedge][0]].id] = iedge;
				if(deb) mesh->part_nodes[elem->nodes[EdgeVerticesMap[iedge][1]].id] = iedge;
			}
		}

		iedge = 7;
		if(elem->nodes[EdgeVerticesMap[iedge][0]].x == 0 && elem->nodes[EdgeVerticesMap[iedge][1]].x == 0){
			if(elem->nodes[EdgeVerticesMap[iedge][0]].y == 3*mesh->ncelly && elem->nodes[EdgeVerticesMap[iedge][1]].y == 3*mesh->ncelly){
				elem->edge[iedge].ref = true;
				if(deb) mesh->part_nodes[elem->nodes[EdgeVerticesMap[iedge][0]].id] = iedge;
				if(deb) mesh->part_nodes[elem->nodes[EdgeVerticesMap[iedge][1]].id] = iedge;
			}
		}

		iedge = 8;
		if(elem->nodes[EdgeVerticesMap[iedge][0]].y == 0 && elem->nodes[EdgeVerticesMap[iedge][1]].y == 0){
			if(elem->nodes[EdgeVerticesMap[iedge][0]].z == 3*mesh->max_z && elem->nodes[EdgeVerticesMap[iedge][1]].z == 3*mesh->max_z){
				elem->edge[iedge].ref = true;
				if(deb) mesh->part_nodes[elem->nodes[EdgeVerticesMap[iedge][0]].id] = iedge;
				if(deb) mesh->part_nodes[elem->nodes[EdgeVerticesMap[iedge][1]].id] = iedge;
			}
		}

		iedge = 9;
		if(elem->nodes[EdgeVerticesMap[iedge][0]].x == 3*mesh->ncellx && elem->nodes[EdgeVerticesMap[iedge][1]].x == 3*mesh->ncellx){
			if(elem->nodes[EdgeVerticesMap[iedge][0]].z == 3*mesh->max_z && elem->nodes[EdgeVerticesMap[iedge][1]].z == 3*mesh->max_z){
				elem->edge[iedge].ref = true;
				if(deb) mesh->part_nodes[elem->nodes[EdgeVerticesMap[iedge][0]].id] = iedge;
				if(deb) mesh->part_nodes[elem->nodes[EdgeVerticesMap[iedge][1]].id] = iedge;
			}
		}

		iedge = 10;
		if(elem->nodes[EdgeVerticesMap[iedge][0]].y == 3*mesh->ncelly && elem->nodes[EdgeVerticesMap[iedge][1]].y == 3*mesh->ncelly){
			if(elem->nodes[EdgeVerticesMap[iedge][0]].z == 3*mesh->max_z && elem->nodes[EdgeVerticesMap[iedge][1]].z == 3*mesh->max_z){
				elem->edge[iedge].ref = true;
				if(deb) mesh->part_nodes[elem->nodes[EdgeVerticesMap[iedge][0]].id] = iedge;
				if(deb) mesh->part_nodes[elem->nodes[EdgeVerticesMap[iedge][1]].id] = iedge;
			}
		}

		iedge = 11;
		if(elem->nodes[EdgeVerticesMap[iedge][0]].x == 0 && elem->nodes[EdgeVerticesMap[iedge][1]].x == 0){
			if(elem->nodes[EdgeVerticesMap[iedge][0]].z == 3*mesh->max_z && elem->nodes[EdgeVerticesMap[iedge][1]].z == 3*mesh->max_z){
				elem->edge[iedge].ref = true;
				if(deb) mesh->part_nodes[elem->nodes[EdgeVerticesMap[iedge][0]].id] = iedge;
				if(deb) mesh->part_nodes[elem->nodes[EdgeVerticesMap[iedge][1]].id] = iedge;
			}
		}

	}

	//id global exterior nodes
	for(int iel = 0; iel < mesh->outsurf.elem_count; iel++){
		octant_t* elem = (octant_t*) sc_array_index (&mesh->outsurf, iel);
		if(elem->edge[4].ref){
			if(elem->nodes[0].z == 0) elem->nodes[0].fixed = -1;
			if(deb && elem->nodes[0].z == 0) mesh->part_nodes[elem->nodes[0].id] = -1;
			if(elem->nodes[4].z == 3*mesh->max_z) elem->nodes[4].fixed = -5;
			if(deb && elem->nodes[4].z == 3*mesh->max_z) mesh->part_nodes[elem->nodes[4].id] = -5;
		}
		if(elem->edge[5].ref){
			if(elem->nodes[1].z == 0) elem->nodes[1].fixed = -2;
			if(deb && elem->nodes[1].z == 0) mesh->part_nodes[elem->nodes[1].id] = -2;
			if(elem->nodes[5].z == 3*mesh->max_z) elem->nodes[5].fixed = -6;
			if(deb && elem->nodes[5].z == 3*mesh->max_z) mesh->part_nodes[elem->nodes[5].id] = -6;
		}
		if(elem->edge[6].ref){
			if(elem->nodes[2].z == 0) elem->nodes[2].fixed = -3;
			if(deb && elem->nodes[2].z == 0) mesh->part_nodes[elem->nodes[2].id] = -3;
			if(elem->nodes[6].z == 3*mesh->max_z) elem->nodes[6].fixed = -7;
			if(deb && elem->nodes[6].z == 3*mesh->max_z) mesh->part_nodes[elem->nodes[6].id] = -7;
		}
		if(elem->edge[7].ref){
			if(elem->nodes[3].z == 0) elem->nodes[3].fixed = -4;
			if(deb && elem->nodes[3].z == 0) mesh->part_nodes[elem->nodes[3].id] = -4;
			if(elem->nodes[7].z == 3*mesh->max_z) elem->nodes[7].fixed = -8;
			if(deb && elem->nodes[7].z == 3*mesh->max_z) mesh->part_nodes[elem->nodes[7].id] = -8;
		}
	}

	if(deb){
		//debug
		for(int iel = 0; iel < mesh->outsurf.elem_count; iel++){
			octant_t* elem = (octant_t*) sc_array_index (&mesh->outsurf, iel);
			for(int ino = 0; ino < 8; ino++){
				if(elem->nodes[ino].fixed == 1) mesh->part_nodes[elem->nodes[ino].id] = 1;
			}
		}

		for(int iel = 0; iel < mesh->outsurf.elem_count; iel++){
			octant_t* elem = (octant_t*) sc_array_index (&mesh->outsurf, iel);

			printf("Sou o elemento: %d\n",elem->id);
			printf("Surface: \n",elem->id);
			for(int isurf = 0; isurf < 6; isurf++) printf("%s ", elem->surf[isurf].ext ? "true" : "false");

			printf("\nEdge: \n",elem->id);
			for(int isurf = 0; isurf < 12; isurf++) printf("%d ",elem->edge[isurf].ref);

			printf("\nNode: \n",elem->id);
			for(int isurf = 0; isurf < 8; isurf++) printf("%d ",elem->nodes[isurf].fixed);
			printf("\n",elem->id);

		}
	}
	//printf("Numero de elementos de superficie %d\n", mesh->outsurf.elem_count);

}


void ShrinkSet(hexa_tree_t* mesh, std::vector<double>& coords, std::vector<int>& nodes_b_mat){

	bool deb = false;
	bool clamped = true;
	//criando a hash de nos para evitar nos duplicados no pillowing
	sc_hash_array_t*   hash_nodes  = (sc_hash_array_t *)sc_hash_array_new(sizeof(octant_node_t), node_hash_fn , node_equal_fn, &clamped);
	//hash nodes
	for(int ino = 0; ino < mesh->nodes.elem_count; ino++)
	{
		size_t position;
		octant_node_t *r;
		octant_node_t key;
		octant_node_t* node = (octant_node_t*) sc_array_index (&mesh->nodes, ino);
		key.x = node->x;
		key.y = node->y;
		key.z = node->z;
		r = (octant_node_t*) sc_hash_array_insert_unique(hash_nodes, &key, &position);
		if(r!=NULL){
			r->x = node->x;
			r->y = node->y;
			r->z = node->z;
			r->id = node->id;
		}else{
			printf("Verificar o no numero %d\n",node->id);
			octant_node_t* node_i = (octant_node_t*) sc_array_index (&hash_nodes->a, position);
			printf("Ele foi confundido com o no %d\n", node_i->id);
		}
	}

	assert(hash_nodes->a.elem_count == mesh->nodes.elem_count);
	//hash nodes nodes_b_mat
	sc_hash_array_t*   hash_b_mat  = (sc_hash_array_t *)sc_hash_array_new(sizeof(octant_node_t), node_hash_fn, node_equal_fn, &clamped);

	//fazendo hash nos nodes_b_mat
	for(int ino = 0; ino < nodes_b_mat.size(); ino++)
	{
		size_t position;
		octant_node_t *r;
		octant_node_t key;
		octant_node_t* node = (octant_node_t*) sc_array_index (&mesh->nodes, nodes_b_mat[ino]);
		key.x = node->x;
		key.y = node->y;
		key.z = node->z;

		r = (octant_node_t*) sc_hash_array_insert_unique(hash_b_mat, &key, &position);
		if(r!=NULL){
			r->x = key.x;
			r->y = key.y;
			r->z = key.z;
			r->id = node->id;
		}else{
			printf("Verificar o no numero %d\n",node);
		}
	}

	sc_hash_array_t*   pillow  = (sc_hash_array_t *)sc_hash_array_new(sizeof(pillow_t), pillow_hash_fn, pillow_equal_fn, &clamped);
	for(int ioc = 0; ioc < mesh->oct.elem_count; ioc++)
	{
		octree_t* oct = (octree_t*) sc_array_index(&mesh->oct,ioc);

		//en train de faire la hash du pillowing
		for(int iel = 0; iel < 8; iel++)
		{
			octant_t* elem = (octant_t*) sc_array_index(&mesh->elements, oct->id[iel]);
			for (int ino = 0; ino < 8; ino++)
			{
				//check si le noeud est sur la hahs
				//hash_b_mat
				size_t position;
				octant_node_t key;
				key.x = elem->nodes[ino].x;
				key.y = elem->nodes[ino].y;
				key.z = elem->nodes[ino].z;
				key.id = elem->nodes[ino].id;

				bool lnode = sc_hash_array_lookup(hash_b_mat, &key, &position);

				if(lnode)
				{
					if(deb) printf("Eu achei o no %d (%d) no el %d (%d) do octree %d\n",ino,key.id,iel,elem->id,ioc);

					//on commance a ajouter les noeuds dans la hahs de pillow
					pillow_t keyP;
					size_t positionP;
					keyP.x = elem->nodes[ino].x;
					keyP.y = elem->nodes[ino].y;
					keyP.z = elem->nodes[ino].z;
					pillow_t* p = (pillow_t*) sc_hash_array_insert_unique(pillow, &keyP, &positionP);
					if(p!=NULL)
					{
						p->id = elem->nodes[ino].id;
						p->a = -1;
						p->b = -1;
						p->x = keyP.x;
						p->y = keyP.y;
						p->z = keyP.z;
						p->elem[0] = elem->id;
						p->list_elem = 1;
						// on doit gerer les surfaces que sont fixes...
						p->list_face[0] = 0;
						for(int isurf = 0; isurf < 3; isurf++)
						{
							bool surf = true;
							if(deb) printf("Testando superficie %d do elemento %d nome do no %d",VertexSurfMap[ino][isurf], iel,elem->nodes[ino].id);
							for(int ive = 0; ive < 4; ive++)
							{
								size_t positionSurfVert;
								octant_node_t SurfVert;
								int refnode = FaceNodesMap[VertexSurfMap[ino][isurf]][ive];
								SurfVert.x = elem->nodes[refnode].x;
								SurfVert.y = elem->nodes[refnode].y;
								SurfVert.z = elem->nodes[refnode].z;
								bool lvert = sc_hash_array_lookup(hash_b_mat, &SurfVert, &positionSurfVert);
								if(!lvert) surf = false;
							}
							if(deb) printf("e ela é %d\n",surf);
							if(surf) p->face[0][p->list_face[0]] = VertexSurfMap[ino][isurf];
							if(surf) p->list_face[0]++;
						}
					}
					else
					{
						pillow_t* p1 = (pillow_t*) sc_array_index(&pillow->a,positionP);
						p1->elem[p1->list_elem] = elem->id;
						p1->list_face[p1->list_elem] = 0;

						for(int isurf = 0; isurf < 3; isurf++)
						{
							bool surf = true;
							if(deb) printf("Testando superficie %d do elemento %d nome do no %d ",VertexSurfMap[ino][isurf], iel,elem->nodes[ino].id);
							for(int ive = 0; ive < 4; ive++)
							{
								size_t positionSurfVert;
								octant_node_t SurfVert;
								int refnode = FaceNodesMap[VertexSurfMap[ino][isurf]][ive];
								SurfVert.x = elem->nodes[refnode].x;
								SurfVert.y = elem->nodes[refnode].y;
								SurfVert.z = elem->nodes[refnode].z;
								bool lvert = sc_hash_array_lookup(hash_b_mat, &SurfVert, &positionSurfVert);
								if(!lvert) surf = false;
							}
							if(deb) printf("e ela é %d\n",surf);
							if(surf) p1->face[p1->list_elem][p1->list_face[p1->list_elem]] = VertexSurfMap[ino][isurf];
							if(surf) p1->list_face[p1->list_elem]++;
						}
						p1->list_elem++;
					}
				}
			}
		}

		if(deb)
		{
			printf("Nombre de éléments dans Pillow:%d\n", pillow->a.elem_count);
			for(int ino = 0; ino < pillow->a.elem_count; ino ++)
			{
				pillow_t* p = (pillow_t*) sc_array_index(&pillow->a, ino);
				printf("    Noeud %d x:%d y:%d z:%d\n",p->id,p->x,p->y,p->z);
				printf("        les éléments que partagent le noeud sont: ");
				for(int i = 0; i < p->list_elem; i++)
				{
					printf("%d ",p->elem[i]);
				}
				printf("\n");
				printf("        les surfaces sur les éléments que partagent le noeud sont: \n");
				for(int i = 0; i < p->list_elem; i++)
				{
					printf("        pour le élément:%d Surface:",p->elem[i]);
					for(int j = 0; j < p->list_face[i]; j++) printf(" %d",p->face[i][j]);
					printf("\n");
				}
			}
		}

		//faire la géneration des noeuds
		for(int ino = 0; ino < pillow->a.elem_count; ino ++)
		{
			pillow_t* p = (pillow_t*) sc_array_index(&pillow->a, ino);
			if(deb)printf("    Noeud %d x:%d y:%d z:%d\n",p->id,p->x,p->y,p->z);

			//on cherche ce que on plus de surface
			int c = 0;
			std::vector<int>auxsurf;
			for(int iel = 0; iel < p->list_elem; iel++)
			{
				const int b = p->list_face[iel];
				c = std::max(c,b);
				if(deb)printf("        pour le élément:%d Surface:",p->elem[iel]);
				if(deb)for(int j = 0; j < p->list_face[iel]; j++) printf(" %d", p->face[iel][j]);
				if(deb)printf("\n");
			}

			//on cherche le élément avec plus de surface pour génerer les noeuds
			for(int iel = 0; iel < p->list_elem; iel++)
			{
				if(p->list_face[iel] == c && c!=0)
				{
					for(int isurf = 0; isurf < p->list_face[iel]; isurf++)
					{
						//printf("%d %d\n",p->id,p->face[iel][isurf]);
						auxsurf.push_back(p->face[iel][isurf]);
					}
					break;
				}
			}

			std::sort(auxsurf.begin(), auxsurf.end());
			if(auxsurf.size()!=0)
			{

				//le noeud base:
				int32_t x = p->x;
				int32_t y = p->y;
				int32_t z = p->z;
				//noeud a
				int32_t xa = p->x;
				int32_t ya = p->y;
				int32_t za = p->z;
				//noeud b
				int32_t xb = p->x;
				int32_t yb = p->y;
				int32_t zb = p->z;

				if(auxsurf.size() == 1)
				{
					if(auxsurf[0] == 0 || auxsurf[0] == 1 )
					{
						xa += 6;
						xb -= 6;
					}
					if(auxsurf[0] == 2 || auxsurf[0] == 3 )
					{
						ya += 6;
						yb -= 6;
					}
					if(auxsurf[0] == 4 || auxsurf[0] == 5 )
					{
						za += 6;
						zb -= 6;
					}
				}
				if(auxsurf.size() == 2)
				{
					if(auxsurf[0] == 0 || auxsurf[0] == 1)
					{
						xa += 6;
						xb -= 6;
						if(auxsurf[1] == 2 || auxsurf[1] == 3)
						{
							ya += 6;
							yb -= 6;
						}
						if(auxsurf[1] == 4 || auxsurf[1] == 5)
						{
							za += 6;
							zb -= 6;
						}
					}
					if(auxsurf[0] == 2 || auxsurf[0] == 3)
					{
						ya += 6;
						yb -= 6;
						if(auxsurf[1] == 4 || auxsurf[1] == 5)
						{
							za += 6;
							zb -= 6;
						}
					}
				}
				if(auxsurf.size() == 3)
				{
					xa += 6;
					xb -= 6;
					ya += 6;
					yb -= 6;
					za += 6;
					zb -= 6;
				}

				//if(deb)printf("%d %d\n",p->id,auxsurf.size());
				//if(deb)printf("%d %d %d\n",xa,ya,za);
				//if(deb)printf("%d %d %d\n",xb,yb,zb);

				//ajouter sur la hash de noeuds
				size_t position;
				octant_node_t key;
				key.x = xa;
				key.y = ya;
				key.z = za;
				octant_node_t* ra = (octant_node_t*) sc_hash_array_insert_unique(hash_nodes, &key, &position);
				if(ra!=NULL){
					ra->x = xa;
					ra->y = ya;
					ra->z = za;
					ra->id = hash_nodes->a.elem_count-1;
					p->a = ra->id;
					ra->fixed = 0;
					ra->color = 0;
					coords.push_back(0);
					coords.push_back(0);
					coords.push_back(0);
				}else{
					octant_node_t* ra = (octant_node_t*) sc_array_index(&hash_nodes->a,position);
					p->a = ra->id;
				}

				key.x = xb;
				key.y = yb;
				key.z = zb;
				octant_node_t* rb = (octant_node_t*) sc_hash_array_insert_unique(hash_nodes, &key, &position);
				if(rb!=NULL){
					rb->x = xb;
					rb->y = yb;
					rb->z = zb;
					rb->id = hash_nodes->a.elem_count-1;
					p->b = rb->id;
					rb->fixed = 0;
					rb->color = 0;
					coords.push_back(0);
					coords.push_back(0);
					coords.push_back(0);
				}else{
					octant_node_t* rb = (octant_node_t*) sc_array_index(&hash_nodes->a,position);
					p->b = rb->id;
				}
			}

			//trouver les coords
			octant_node_t * a = (octant_node_t*) sc_array_index(&hash_nodes->a,p->a);
			octant_node_t * b = (octant_node_t*) sc_array_index(&hash_nodes->a,p->b);
			bool ina[8];
			bool inb[8];
			for(int iel = 0; iel < 8; iel++)
			{
				octant_t* elem = (octant_t*) sc_array_index(&mesh->elements, oct->id[iel]);
				ina[iel] = false;
				inb[iel] = false;
				if(elem->nodes[0].x <= a->x && elem->nodes[6].x >= a->x &&
						elem->nodes[0].y <= a->y && elem->nodes[6].y >= a->y &&
						elem->nodes[0].z <= a->z && elem->nodes[6].z >= a->z) ina[iel] = true;
				if(elem->nodes[0].x <= b->x && elem->nodes[6].x >= b->x &&
						elem->nodes[0].y <= b->y && elem->nodes[6].y >= b->y &&
						elem->nodes[0].z <= b->z && elem->nodes[6].z >= b->z) inb[iel] = true;
				if(deb)printf("%d %d\n",elem->id,p->id);
				if(deb)printf("%d %d %d %d %s\n",a->id,a->x,a->y,a->z, ina[iel] ? "true" : "false");
				if(deb)printf("%d %d %d %d %s\n",b->id,b->x,b->y,b->z, inb[iel] ? "true" : "false");
				if(deb)printf("%d %d %d %d %d %d\n",elem->nodes[0].x,elem->nodes[6].x,elem->nodes[0].y,elem->nodes[6].y,
						elem->nodes[0].z,elem->nodes[6].z);


				if(ina[iel] || inb[iel])
				{
					double cord_in_ref[3];
					int nnode;
					double ref_in_x[8],ref_in_y[8],ref_in_z[8];
					int color = 0;
					int colorc = 0;
					for (int ino = 0; ino < 8; ino++)
					{
						ref_in_x[ino]=coords[3*elem->nodes[ino].id+0];
						ref_in_y[ino]=coords[3*elem->nodes[ino].id+1];
						ref_in_z[ino]=coords[3*elem->nodes[ino].id+2];
						if(elem->nodes[ino].color != 0)
						{
							color += elem->nodes[ino].color;
							colorc++;
						}
					}

					if(ina[iel])
					{
						cord_in_ref[0] = (a->x - elem->nodes[0].x)/6 -1;
						cord_in_ref[1] = (a->y - elem->nodes[0].y)/6 -1;
						cord_in_ref[2] = (a->z - elem->nodes[0].z)/6 -1;
						nnode = a->id;
						a->color = int(color)/int(colorc);
					}
					if(inb[iel])
					{
						cord_in_ref[0] = (b->x - elem->nodes[0].x)/6 -1;
						cord_in_ref[1] = (b->y - elem->nodes[0].y)/6 -1;
						cord_in_ref[2] = (b->z - elem->nodes[0].z)/6 -1;
						nnode = p->b;
						b->color = int(color)/int(colorc);
					}
					GtsPoint* ref_point = LinearMapHex(cord_in_ref, ref_in_x, ref_in_y, ref_in_z);

					coords[3*nnode + 0] = ref_point->x;
					coords[3*nnode + 1] = ref_point->y;
					coords[3*nnode + 2] = ref_point->z;
				}
			}
		}

		//faire la criation des éléments pour le pillowing
		sc_array_t toto;
		sc_array_init(&toto, sizeof(octant_t));
		for(int iel = 0; iel < 8; iel++)
		{
			octant_t* elemOrig = (octant_t*) sc_array_index(&mesh->elements, oct->id[iel]);
			octant_t * elem = (octant_t*) sc_array_push(&toto);

			hexa_element_copy(elemOrig,elem);

			//on doit voir si le élément demande de extrusion...
			//extrusion se fait à chaque surface du élément
			bool surf[6];
			std::vector<int> aux;
			int surf_count = 0;
			for(int isurf = 0; isurf < 6; isurf++)
			{
				surf[isurf] = true;
				for(int ino = 0; ino < 4; ino++)
				{
					size_t position;
					octant_node_t key;
					key.x = elem->nodes[FaceNodesMap[isurf][ino]].x;
					key.y = elem->nodes[FaceNodesMap[isurf][ino]].y;
					key.z = elem->nodes[FaceNodesMap[isurf][ino]].z;

					bool lnode = sc_hash_array_lookup(hash_b_mat, &key, &position);
					if(!lnode) surf[isurf] = false;
				}
				if(surf[isurf]) surf_count++;
				if(surf[isurf]) aux.push_back(isurf);
			}

			if(deb) printf("Para o elemento %d eu vou extrudar %d elementos\n",elem->id,aux.size());

			//on garde la connectivite original...
			int connec[8];
			for(int ino = 0; ino <8; ino++){
				connec[ino] = elem->nodes[ino].id;
			}

			for(int isurf = 0; isurf < aux.size(); isurf++)
			{
				int new_nodes[4];
				int ina,inb;
				for(int ino = 0; ino < 4; ino++)
				{
					ina = inb = 0;
					size_t position;
					pillow_t key;
					key.x = elem->nodes[FaceNodesMap[aux[isurf]][ino]].x;
					key.y = elem->nodes[FaceNodesMap[aux[isurf]][ino]].y;
					key.z = elem->nodes[FaceNodesMap[aux[isurf]][ino]].z;

					bool lnode = sc_hash_array_lookup(pillow, &key, &position);
					pillow_t* p = (pillow_t*) sc_array_index(&pillow->a,position);
					octant_node_t * a = (octant_node_t*) sc_array_index(&hash_nodes->a,p->a);
					octant_node_t * b = (octant_node_t*) sc_array_index(&hash_nodes->a,p->b);
					if(elem->nodes[0].x <= a->x && elem->nodes[1].x >= a->x) ina++;
					if(elem->nodes[0].y <= a->y && elem->nodes[3].y >= a->y) ina++;
					if(elem->nodes[0].z <= a->z && elem->nodes[4].z >= a->z) ina++;
					if(elem->nodes[0].x <= b->x && elem->nodes[1].x >= b->x) inb++;
					if(elem->nodes[0].y <= b->y && elem->nodes[3].y >= b->y) inb++;
					if(elem->nodes[0].z <= b->z && elem->nodes[4].z >= b->z) inb++;

					if(ina > inb)
					{
						new_nodes[ino] = a->id;
					}
					else
					{
						new_nodes[ino] = b->id;
					}
				}

				//créer l'élément
				octant_t* pelem;
				pelem = (octant_t*) sc_array_push(&mesh->elements);
				pelem->id = mesh->elements.elem_count+1;
				pelem->n_mat = elem->n_mat;
				pelem->pad = elem->pad;
				pelem->x = elem->x;
				pelem->y = elem->y;
				pelem->z = elem->z;

				//for(int ino = 0; ino<8; ino++){
				//	pelem->nodes[ino].fixed= 0;
				//	pelem->nodes[ino].x = elem->nodes[ino].x;
				//	pelem->nodes[ino].y = elem->nodes[ino].y;
				//	pelem->nodes[ino].z = elem->nodes[ino].z;
				//}

				for(int ino = 0; ino < 4; ino++)
				{
					pelem->nodes[FaceNodesMap[aux[isurf]][ino]].id = connec[FaceNodesMap[aux[isurf]][ino]];
					pelem->nodes[FaceNodesMap[aux[isurf]][ino]].fixed = elem->nodes[FaceNodesMap[aux[isurf]][ino]].fixed;
					pelem->nodes[FaceNodesMap_inv[aux[isurf]][ino]].id= new_nodes[ino];
					elem->nodes[FaceNodesMap[aux[isurf]][ino]].id = new_nodes[ino];
				}

			}

			//chercher le noeuds pour trouver la bonne connectivite dans les éléments
			for (int ino = 0; ino < 8; ino++)
			{
				int ina = 0;
				int inb = 0;
				pillow_t keyP;
				size_t positionP;
				keyP.x = elem->nodes[ino].x;
				keyP.y = elem->nodes[ino].y;
				keyP.z = elem->nodes[ino].z;
				bool lnode = sc_hash_array_lookup(pillow, &keyP, &positionP);
				if(lnode)
				{
					pillow_t * p = (pillow_t*) sc_array_index(&pillow->a,positionP);
					octant_node_t * a = (octant_node_t*) sc_array_index(&hash_nodes->a,p->a);
					octant_node_t * b = (octant_node_t*) sc_array_index(&hash_nodes->a,p->b);
					if(deb)printf("%d %d %d %d %d %d\n",elem->id,p->id,p->a, a->id,p->b,b->id);

					if(deb)printf("%d %d %d\n",a->x,a->y,a->z);
					if(deb)printf("%d %d %d\n",b->x,b->y,b->z);
					if(deb)printf("%d %d %d %d %d %d\n",elem->nodes[0].x,elem->nodes[1].x,elem->nodes[0].y,elem->nodes[3].y,elem->nodes[0].z,elem->nodes[4].z);
					if(elem->nodes[0].x <= a->x && elem->nodes[1].x >= a->x) ina++;
					if(elem->nodes[0].y <= a->y && elem->nodes[3].y >= a->y) ina++;
					if(elem->nodes[0].z <= a->z && elem->nodes[4].z >= a->z) ina++;
					if(elem->nodes[0].x <= b->x && elem->nodes[1].x >= b->x) inb++;
					if(elem->nodes[0].y <= b->y && elem->nodes[3].y >= b->y) inb++;
					if(elem->nodes[0].z <= b->z && elem->nodes[4].z >= b->z) inb++;
					if(deb)printf("%d %d\n",ina, inb);

					if(ina > inb)
					{
						elem->nodes[ino].id = p->a;
						elem->nodes[ino].x = a->x;
						elem->nodes[ino].y = a->y;
						elem->nodes[ino].z = a->z;
					}
					else
					{
						elem->nodes[ino].id = p->b;
						elem->nodes[ino].x = b->x;
						elem->nodes[ino].y = b->y;
						elem->nodes[ino].z = b->z;
					}

				}
			}
			sc_array_reset(&toto);
		}


	}

	sc_hash_array_destroy(pillow);

	//faire la atualisation de mon tableau de noeuds
	sc_array_reset(&mesh->nodes);
	sc_hash_array_rip(hash_nodes,&mesh->nodes);

}

void PillowingInterface(hexa_tree_t* mesh, std::vector<double>& coords, std::vector<int>& nodes_b_mat)
{
	int elem_old = mesh->elements.elem_count;
	int nodes_old = mesh->nodes.elem_count;
	bool clamped = true;
	fprintf(mesh->profile,"Time inside PillowingInterface\n");

	//redo the mapping in the nodes
	auto start = std::chrono::steady_clock::now( );
	printf("     Redo Node Mapping\n");
	RedoNodeMapping(mesh);
	auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>( std::chrono::steady_clock::now( ) - start );
	fprintf(mesh->profile,"    Redo the mapping in the nodes %lld millisecond(s).\n",elapsed.count());
	//std::cout << "Redo the mapping in the nodes "<< elapsed.count() <<" millisecond(s)."<< std::endl;

	//Identify the boundaries: global and local
	start = std::chrono::steady_clock::now( );
	printf("     Surface Identification\n");
	SurfaceIdentification(mesh, coords);
	fprintf(mesh->profile,"    Time in SurfaceIdentification %lld millisecond(s).\n",elapsed.count());
	//std::cout << "Time SurfaceIdentification "<< elapsed.count() <<" millisecond(s)."<< std::endl;

	//Identify the boundaries: global and local
	start = std::chrono::steady_clock::now( );
	printf("     Pillow Layer\n");
	ShrinkSet(mesh, coords,nodes_b_mat);
	fprintf(mesh->profile,"    Time in PillowLayer %lld millisecond(s).\n",elapsed.count());
	//std::cout << "Time SurfaceIdentification "<< elapsed.count() <<" millisecond(s)."<< std::endl;

	//update the vectors
	mesh->local_n_elements = mesh->elements.elem_count;
	mesh->local_n_nodes = mesh->nodes.elem_count;
	MPI_Allreduce(&mesh->local_n_elements, &mesh->total_n_elements, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(&mesh->local_n_nodes, &mesh->total_n_nodes, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
	printf("     %d nodes were created in the pillowing process\n",mesh->nodes.elem_count - nodes_old);
	printf("     %d elements were created in the pillowing process\n",mesh->elements.elem_count - elem_old);

	free(mesh->part_nodes);
	mesh->part_nodes = (int*) malloc (mesh->local_n_nodes*sizeof(int));
	for (int ino = 0; ino < mesh->local_n_nodes; ino++) {
		mesh->part_nodes[ino]=mesh->mpi_rank;
	}
}

/*
 * for(int ino = 0; ino < mesh->nodes.elem_count; ino++)
	{
		octant_node_t * node = (octant_node_t*) sc_array_index(&mesh->nodes,ino);
		printf("%d %d %d\n",node->x,node->y,node->z);
	}



	sc_array_t toto;
	sc_array_init(&toto, sizeof(octant_t));
	//faire la criation des éléments
	for(int iel = 0; iel < 8; iel++)
	{
		octant_t* elemOrig = (octant_t*) sc_array_index(&mesh->elements, oct->id[iel]);
		octant_t * elem = (octant_t*) sc_array_push(&toto);

		hexa_element_copy(elemOrig,elem);
		//en train de faire l'élément de base pour pouvoir faire la coupure
		double cord_in_x[8],cord_in_y[8],cord_in_z[8];
		double ref_in_x[8], ref_in_y[8], ref_in_z[8];
		GtsPoint* ref_point[8];
		GtsPoint* point[8];
		int conn[8];
		for (int ino = 0; ino < 8; ino++)
		{
			cord_in_x[ino]=coords[3*elem->nodes[ino].id+0];
			cord_in_y[ino]=coords[3*elem->nodes[ino].id+1];
			cord_in_z[ino]=coords[3*elem->nodes[ino].id+2];
			octant_node_t* node = (octant_node_t*) sc_array_index(&mesh->nodes, elem->nodes[ino].id);
			ref_in_x[ino] = node->x;
			ref_in_y[ino] = node->y;
			ref_in_z[ino] = node->z;
		}

		//extrusion se fait à chaque surface du élément
		bool surf[6];
		std::vector<int> aux;
		int surf_count = 0;
		for(int isurf = 0; isurf < 6; isurf++)
		{
			surf[isurf] = true;
			for(int ino = 0; ino < 4; ino++)
			{
				size_t position;
				octant_node_t key;
				key.x = elem->nodes[FaceNodesMap[isurf][ino]].x;
				key.y = elem->nodes[FaceNodesMap[isurf][ino]].y;
				key.z = elem->nodes[FaceNodesMap[isurf][ino]].z;

				bool lnode = sc_hash_array_lookup(hash_b_mat, &key, &position);
				if(!lnode) surf[isurf] = false;
			}
			if(surf[isurf]) surf_count++;
			if(surf[isurf]) aux.push_back(isurf);
		}

		if(deb) printf("O elemento %d(%d) do octree %d vai gerar %d novos elementos\n",iel,elem->id,ioc,surf_count);

		if(surf_count != 0)
		{
			if(deb)printf("Para o elemento %d\n", iel);
			for(int isurf = 0; isurf < aux.size(); isurf++)
			{

				//on cherche le pillow sur le 4 noeuds
				pillow_t* p[4];
				for(int ino = 0; ino < 4; ino ++)
				{
					std::vector<int> auxsurf;
					pillow_t key;
					size_t position;
					key.x = elem->nodes[FaceNodesMap[aux[isurf]][ino]].x;
					key.y = elem->nodes[FaceNodesMap[aux[isurf]][ino]].y;
					key.z = elem->nodes[FaceNodesMap[aux[isurf]][ino]].z;
					bool lnode = sc_hash_array_lookup(pillow, &key, &position);
					p[ino] = (pillow_t*) sc_array_index(&pillow->a, position);
					//faire la generation de 4 Noeud
					for(int i = 0; i < p[ino]->list_elem;i++)
					{
						if(p[ino]->elem[i] == elem->id)
						{
							for(int j = 0; j < p[ino]->list_face[i]; j++)
							{
								auxsurf.push_back(p[ino]->face[i][j]);
								if(deb)printf("%d \n",p[ino]->face[i][j]);
							}
						}
					}
					if(deb)printf("Vou ter que mover o no %d em %d direcoes\n",p[ino]->id,auxsurf.size());

					//fazendo o sort do vetor
					std::sort(auxsurf.begin(), auxsurf.end());
					double cord_in_ref[3];
					if(auxsurf.size() == 1)
					{
						if(auxsurf[0] == 0 || auxsurf[0] == 1)
						{
							cord_in_ref[0] = 0;
							cord_in_ref[1] = coord_ref[FaceNodesMap[auxsurf[0]][ino]][1];
							cord_in_ref[2] = coord_ref[FaceNodesMap[auxsurf[0]][ino]][2];
						}
						if(auxsurf[0] == 2 || auxsurf[0] == 3)
						{
							cord_in_ref[0] = coord_ref[FaceNodesMap[auxsurf[0]][ino]][0];
							cord_in_ref[1] = 0;
							cord_in_ref[2] = coord_ref[FaceNodesMap[auxsurf[0]][ino]][2];
						}
						if(auxsurf[0] == 4 || auxsurf[0] == 5)
						{
							cord_in_ref[0] = coord_ref[FaceNodesMap[auxsurf[0]][ino]][0];
							cord_in_ref[1] = coord_ref[FaceNodesMap[auxsurf[0]][ino]][1];
							cord_in_ref[2] = 0;
						}
					}
					if(auxsurf.size() == 2)
					{
						if(auxsurf[0] == 0 || auxsurf[0] == 1)
						{
							cord_in_ref[0] = 0;
							if(auxsurf[1] == 2 || auxsurf[1] == 3)
							{
								cord_in_ref[1] = 0;
								cord_in_ref[2] = coord_ref[FaceNodesMap[auxsurf[0]][ino]][2];
							}
							if(auxsurf[1] == 4 || auxsurf[1] == 5)
							{
								cord_in_ref[1] = coord_ref[FaceNodesMap[auxsurf[0]][ino]][1];
								cord_in_ref[2] = 0;
							}
						}
						if(auxsurf[0] == 2 || auxsurf[0] == 3)
						{
							cord_in_ref[1] = 0;
							if(auxsurf[1] == 4 || auxsurf[1] == 5)
							{
								cord_in_ref[0] = coord_ref[FaceNodesMap[auxsurf[0]][ino]][0];
								cord_in_ref[2] = 0;
							}
						}
					}
					if(auxsurf.size() == 3)
					{
						if(auxsurf[0] == 0 || auxsurf[0] == 1)
						{
							cord_in_ref[0] = 0;
							if(auxsurf[1] == 2 || auxsurf[1] == 3)
							{
								cord_in_ref[1] = 0;
								if(auxsurf[2] == 2 || auxsurf[2] == 3)
								{
									cord_in_ref[2] = 0;
								}
							}
						}
					}


					//////
					ref_point[ino] = LinearMapHex(cord_in_ref, ref_in_x, ref_in_y, ref_in_z);
					double auxx = ref_point[ino]->x;
					double auxy = ref_point[ino]->y;
					double auxz = ref_point[ino]->z;
					int x = round(auxx);
					int y = round(auxy);
					int z = round(auxz);
					printf("%d %d %d\n",x,y,z);

					point[FaceNodesMap[auxsurf[isurf]][ino]] = LinearMapHex(cord_in_ref, cord_in_x,cord_in_y,cord_in_z);
					conn[FaceNodesMap[auxsurf[isurf]][ino]] = -1;
					conn[FaceNodesMap[auxsurf[isurf]][ino]] =
							AddPoint( mesh, hash_nodes,point[FaceNodesMap[auxsurf[isurf]][ino]], coords, x, y, z);
					elem->nodes[FaceNodesMap_inv[auxsurf[isurf]][ino]].id = conn[FaceNodesMap[auxsurf[isurf]][ino]] ;


				}




			}
		}
	}


	//mesh->oct.elem_count
	for(int ioc = 0; ioc < 0; ioc++)
	{
		octree_t* oct = (octree_t*) sc_array_index(&mesh->oct,ioc);

		//en train de faire la hash du pillowing
		sc_hash_array_t*   pillow  = (sc_hash_array_t *)sc_hash_array_new(sizeof(pillow_t), pillow_hash_fn, pillow_equal_fn, &clamped);
		for(int iel = 0; iel < 8; iel++)
		{
			octant_t* elem = (octant_t*) sc_array_index(&mesh->elements, oct->id[iel]);

			for (int ino = 0; ino < 8; ino++)
			{
				size_t position;
				octant_node_t key;
				key.x = elem->nodes[ino].x;
				key.y = elem->nodes[ino].y;
				key.z = elem->nodes[ino].z;
				bool lnode = sc_hash_array_lookup(hash_b_mat, &key, &position);

				if(lnode)
				{
					pillow_t keyP;
					keyP.x = elem->nodes[ino].x;
					keyP.y = elem->nodes[ino].y;
					keyP.z = elem->nodes[ino].z;
					pillow_t* p = (pillow_t*) sc_hash_array_insert_unique(pillow, &keyP, &position);
					if(p!=NULL)
					{
						p->id = elem->nodes[ino].id;
						p->x = elem->nodes[ino].x;
						p->y = elem->nodes[ino].y;
						p->z = elem->nodes[ino].z;
						p->list_elem = 1;
						p->elem[p->list_elem-1] = elem->id;
						p->list_face[iel] = 0;

						for(int isurf = 0; isurf < 3; isurf++)
						{
							//printf("Node %d surface %d\n",ino, VertexSurfMap[ino][isurf]);
							bool surf = true;
							//printf("Node %d surface %d\n",ino, VertexSurfMap[ino][isurf]);
							//printf("Verificando os id ref: ");
							for(int ive = 0; ive < 4; ive++)
							{
								size_t position;
								octant_node_t key;
								//printf("%d ",FaceNodesMap[VertexSurfMap[ino][isurf]][ive]);

								key.x = elem->nodes[FaceNodesMap[VertexSurfMap[ino][isurf]][ive]].x;
								key.y = elem->nodes[FaceNodesMap[VertexSurfMap[ino][isurf]][ive]].y;
								key.z = elem->nodes[FaceNodesMap[VertexSurfMap[ino][isurf]][ive]].z;

								bool lnode = sc_hash_array_lookup(hash_b_mat, &key, &position);
								if(!lnode) surf = false;
							}
							//printf("\n");

							if(surf) p->face[iel][p->list_face[iel]] = VertexSurfMap[ino][isurf];
							if(surf) p->list_face[iel]++;
						}

					}
					else
					{
						p = (pillow_t*) sc_array_index(&pillow->a, position);
						p->elem[p->list_elem] = elem->id;
						p->list_face[iel] = 0;
						p->list_elem++;

						for(int isurf = 0; isurf < 3; isurf++)
						{
							bool surf = true;
							//printf("Node %d surface %d\n",ino, VertexSurfMap[ino][isurf]);
							//printf("Verificando os id ref: ");

							for(int ive = 0; ive < 4; ive++)
							{
								size_t position;
								octant_node_t key;
								//printf("%d ",FaceNodesMap[VertexSurfMap[ino][isurf]][ive]);
								key.x = elem->nodes[FaceNodesMap[VertexSurfMap[ino][isurf]][ive]].x;
								key.y = elem->nodes[FaceNodesMap[VertexSurfMap[ino][isurf]][ive]].y;
								key.z = elem->nodes[FaceNodesMap[VertexSurfMap[ino][isurf]][ive]].z;
								bool lnode = sc_hash_array_lookup(hash_b_mat, &key, &position);
								if(!lnode) surf = false;
							}
							//printf("\n");
							if(surf) p->face[iel][p->list_face[iel]] = VertexSurfMap[ino][isurf];
							if(surf) p->list_face[iel]++;
						}


					}
				}
			}
		}

		for(int ino = 0; ino < 0; ino ++)
		{
			pillow_t* p = (pillow_t*) sc_array_index(&pillow->a, ino);
			for(int i = 0; i < p->list_elem; i++){
				octant_t* elem = (octant_t*) sc_array_index(&mesh->elements, p->elem[i]);
				bool surf[6];
				for(int isurf = 0; isurf < 6; isurf++)
				{
					surf[isurf] = true;
					for(int ino = 0; ino < 4; ino ++)
					{
						octant_node_t key;
						size_t position;
						key.x = elem->nodes[FaceNodesMap[isurf][ino]].x;
						key.y = elem->nodes[FaceNodesMap[isurf][ino]].y;
						key.z = elem->nodes[FaceNodesMap[isurf][ino]].z;
						bool lnode = sc_hash_array_lookup(hash_b_mat, &key, &position);
						if(!lnode) surf[isurf] = false;
					}

					if(surf[isurf])
					{
						for(int ino = 0; ino < 4; ino ++)
						{
							if(elem->nodes[FaceNodesMap[isurf][ino]].id == p->id){
								p->face[i][p->list_face[i]] = isurf;
								p->list_face[i]++;
							}
						}
					}
				}
			}
		}

		//////////////////////////////
		//DEBUG
		//////////////////////////////
		if(true)
		{
			printf("Numero de elementos no pillow %d\n",pillow->a.elem_count);
			for(int ino = 0; ino < pillow->a.elem_count; ino ++)
			{
				pillow_t* p = (pillow_t*) sc_array_index(&pillow->a, ino);
				printf("    No %d x:%d y:%d z:%d\n",p->id,p->x,p->y,p->z);
				for(int i = 0; i < p->list_elem; i++){
					printf("    elemento:%d numero de faces:%d\n",p->elem[i],p->list_face[i]);
					for(int j = 0; j < p->list_face[i]; j++)
					{
						printf("    face:%d\n",p->face[i][j]);
					}
					if(p->list_face[i] == 0 && false)
					{
						octant_t* elem = (octant_t*) sc_array_index(&mesh->elements, p->elem[i]);
						bool surf[6];
						for(int isurf = 0; isurf < 6; isurf++)
						{
							surf[isurf] = true;
							for(int ino = 0; ino < 4; ino ++)
							{
								octant_node_t key;
								size_t position;
								key.x = elem->nodes[FaceNodesMap[isurf][ino]].x;
								key.y = elem->nodes[FaceNodesMap[isurf][ino]].y;
								key.z = elem->nodes[FaceNodesMap[isurf][ino]].z;
								bool lnode = sc_hash_array_lookup(hash_b_mat, &key, &position);
								if(!lnode) surf[isurf] = false;
							}
							if(surf[isurf]) printf("%d ", isurf);
						}
						printf("\n");
					}

				}

			}
		}
		//////////////////////////////
		//DEBUG
		//////////////////////////////

		sc_array_t toto;
		sc_array_init(&toto, sizeof(octant_t));
		//faire la criation des éléments
		for(int iel = 0; iel < 0; iel++)
		{
			octant_t* elemOrig = (octant_t*) sc_array_index(&mesh->elements, oct->id[iel]);
			octant_t * elem = (octant_t*) sc_array_push(&toto);

			hexa_element_copy(elemOrig,elem);

			//voir si on a le surface
			//pour faire la extrusion du élémént
			//ou des éléments
			bool surf[6];
			int surf_count =0;
			for(int isurf = 0; isurf < 6; isurf++)
			{
				surf[isurf] = true;
				for(int ino = 0; ino < 4; ino++)
				{
					size_t position;
					octant_node_t key;
					key.x = elem->nodes[FaceNodesMap[isurf][ino]].x;
					key.y = elem->nodes[FaceNodesMap[isurf][ino]].y;
					key.z = elem->nodes[FaceNodesMap[isurf][ino]].z;

					bool lnode = sc_hash_array_lookup(hash_b_mat, &key, &position);
					if(!lnode) surf[isurf] = false;
				}
				if(surf[isurf]) surf_count++;
			}

			//en train de faire l'élément de base pour pouvoir faire la coupure
			double cord_in_x[8],cord_in_y[8],cord_in_z[8];
			double ref_in_x[8], ref_in_y[8], ref_in_z[8];
			GtsPoint* ref_point[8];
			GtsPoint* point[8];
			int conn[8];
			for (int ino = 0; ino < 8; ino++)
			{
				cord_in_x[ino]=coords[3*elem->nodes[ino].id+0] ;
				cord_in_y[ino]=coords[3*elem->nodes[ino].id+1] ;
				cord_in_z[ino]=coords[3*elem->nodes[ino].id+2] ;
				octant_node_t* node = (octant_node_t*) sc_array_index(&mesh->nodes, elem->nodes[ino].id);
				ref_in_x[ino] = node->x;
				ref_in_y[ino] = node->y;
				ref_in_z[ino] = node->z;
			}

			if(surf_count!=0)
			{
				for(int isurf = 0; isurf < 6; isurf++)
				{
					if(surf[isurf])
					{
						for(int ino = 0; ino < 4; ino ++)
						{
							pillow_t key;
							size_t position;
							key.x = elem->nodes[FaceNodesMap[isurf][ino]].x;
							key.y = elem->nodes[FaceNodesMap[isurf][ino]].y;
							key.z = elem->nodes[FaceNodesMap[isurf][ino]].z;
							bool lnode = sc_hash_array_lookup(pillow, &key, &position);
							if(lnode)
							{
								pillow_t* p = (pillow_t*) sc_array_index(&pillow->a, position);


							}
						}
					}

				}
			}







			for(int isurf = 0; isurf < 6; isurf++)
			{
				int face = isurf;
				//fazendo o shrink em funcao das faces
				if(surf_oct[face])
				{
					for(int ino = 0; ino < 4; ino ++)
					{
						// 1D extrude
						if(count[FaceNodesMap[face][ino]] == 1)
						{
							double cord_in_ref[3];
							if(isurf == 0 ||isurf == 1)
							{
								cord_in_ref[0] = 0;
								cord_in_ref[1] = coord_ref[FaceNodesMap[face][ino]][1];
								cord_in_ref[2] = coord_ref[FaceNodesMap[face][ino]][2];
							}
							if(isurf == 2 ||isurf == 3)
							{
								cord_in_ref[1] = 0;
								cord_in_ref[0] = coord_ref[FaceNodesMap[face][ino]][0];
								cord_in_ref[2] = coord_ref[FaceNodesMap[face][ino]][2];
							}
							if(isurf == 4 ||isurf == 5)
							{
								cord_in_ref[2] = 0;
								cord_in_ref[1] = coord_ref[FaceNodesMap[face][ino]][1];
								cord_in_ref[0] = coord_ref[FaceNodesMap[face][ino]][0];
							}
							ref_point[ino] = LinearMapHex(cord_in_ref, ref_in_x, ref_in_y, ref_in_z);
							double auxx = ref_point[ino]->x;
							double auxy = ref_point[ino]->y;
							double auxz = ref_point[ino]->z;
							int x = round(auxx);
							int y = round(auxy);
							int z = round(auxz);
							point[FaceNodesMap[face][ino]] = LinearMapHex(cord_in_ref, cord_in_x,cord_in_y,cord_in_z);
							conn[FaceNodesMap[face][ino]] = -1;
							conn[FaceNodesMap[face][ino]] = AddPoint( mesh, hash_nodes,
									point[FaceNodesMap[face][ino]], coords, x, y, z);
							elem->nodes[FaceNodesMap[face][ino]].id = conn[FaceNodesMap[face][ino]] ;
						}
					}
				}


			}

			sc_array_reset(&toto);

		}
		sc_hash_array_destroy(pillow);
	}
 */
/*
unsigned edget_id_hash(const void *v, const void *u) {
	const edge_t *q = (const edge_t*) v;
	uint64_t a, b, c;

	a = (uint64_t) q->id;
	b = (uint64_t) 1;
	c = (uint64_t) 1;
	sc_hash_mix(a, b, c);
	sc_hash_final(a, b, c);
	return (unsigned) c;
}

int edget_id_equal(const void *v, const void *u, const void *w) {
	const edge_t *e1 = (const edge_t*) v;
	const edge_t *e2 = (const edge_t*) u;

	return (unsigned) ((e1->id == e2->id));

}

unsigned id_hash(const void *v, const void *u) {
	const octant_edge_t *q = (const octant_edge_t*) v;
	uint64_t a, b, c;

	a = (uint64_t) q->id;
	b = (uint64_t) 1;
	c = (uint64_t) 1;
	sc_hash_mix(a, b, c);
	sc_hash_final(a, b, c);
	return (unsigned) c;
}

int id_equal(const void *v, const void *u, const void *w) {
	const octant_edge_t *e1 = (const octant_edge_t*) v;
	const octant_edge_t *e2 = (const octant_edge_t*) u;

	return (unsigned) ((e1->id == e2->id));

}

void edge_add(uint64_t id, sc_hash_array_t* hash_edge_ref) {
	size_t position;
	octant_edge_t *r;
	octant_edge_t key;
	key.id = id;

	r = (octant_edge_t*) sc_hash_array_insert_unique(hash_edge_ref, &key, &position);
	if(r != NULL){
		r->id = key.id;
		r->ref = true;
	}
}
 */

/*
void BuildHash(hexa_tree_t* mesh, std::vector<double>& coords, std::vector<int>& nodes_b_mat,sc_hash_array_t*hash_edge_ref){


	//hash nodes
	for(int ino = 0; ino < mesh->nodes.elem_count; ino++){
		size_t position;
		octant_node_t *r;
		octant_node_t key;
		octant_node_t* node = (octant_node_t*) sc_array_index (&mesh->nodes, ino);
		key.x = node->x;
		key.y = node->y;
		key.z = node->z;
		r = (node_t*) sc_hash_array_insert_unique(hash_nodes, &key, &position);
		if(r!=NULL){
			r->x = node->x;
			r->y = node->y;
			r->z = node->z;
			r->id = node->id;
		}else{
			printf("Verificar o no numero %d\n",node->id);
			octant_node_t* node_i = (octant_node_t*) sc_array_index (&hash_nodes->a, position);
			printf("Ele foi confundido com o no %d\n", node_i->id);
		}
	}

	assert(hash_nodes->a.elem_count == mesh->nodes.elem_count);

	//fazendo vertex hash
	for(int ioc = 0 ; ioc < mesh->oct.elem_count; ioc++){
		octree_t * oc = (octree_t*) sc_array_index (&mesh->oct, ioc);
		for(int iel = 0; iel < 8 ; iel++){
			size_t  position;
			octant_t *elem = (octant_t*) sc_array_index(&mesh->elements, oc->id[iel]);

			for (int ino = 0; ino < 8; ino++){
				octant_vertex_t key;
				key.id = elem->nodes[ino].id;
				octant_vertex_t* vert = (octant_vertex_t*) sc_hash_array_insert_unique (vertex_hash, &key, &position);
				if(vert != NULL){
					vert->id = elem->nodes[ino].id;
					vert->list_elem = 1;
					vert->elem[vert->list_elem-1] = elem->id;
				}else{
					vert = (octant_vertex_t*) sc_array_index(&vertex_hash->a, position);
					vert->elem[vert->list_elem] = elem->id;
					vert->list_elem++;
				}
			}
		}
	}



	bool clamped = true;
	//hash nodes nodes_b_mat
	sc_hash_array_t*   hash_b_mat  = (sc_hash_array_t *)sc_hash_array_new(sizeof(node_t), edge_hash_fn, edge_equal_fn, &clamped);


	//fazendo hash nos nodes_b_mat
	for(int n = 0; n < nodes_b_mat.size(); n++){
		size_t position;
		node_t *r;
		node_t key;
		int node = nodes_b_mat[n];
		key.coord[0] = coords[3*node+0];
		key.coord[1] = coords[3*node+1];
		key.coord[2] = coords[3*node+2];
		key.node_id = node;

		r = (node_t*) sc_hash_array_insert_unique(hash_b_mat, &key, &position);
		if(r!=NULL){
			r->coord[0] = coords[3*node+0];
			r->coord[1] = coords[3*node+1];
			r->coord[2] = coords[3*node+2];
			r->node_id = node;
			r->flag = true;
		}else{
			printf("Verificar o no numero %d\n",node);
		}
	}
	//fazendo edge hash
	//mesh->oct.elem_count
	for(int ioc = 0 ; ioc < mesh->oct.elem_count; ioc++){
		octree_t * oc = (octree_t*) sc_array_index (&mesh->oct, ioc);
		for(int iel = 0; iel < 8 ; iel++){
			octant_t *elem = (octant_t*) sc_array_index(&mesh->elements, oc->id[iel]);
			size_t  position;
			int ncount = 0;
			int aux[8];
			for(int ino = 0; ino < 8; ino++){
				node_t key;
				int node = elem->nodes[ino].id;
				key.coord[0] = coords[3*node+0];
				key.coord[1] = coords[3*node+1];
				key.coord[2] = coords[3*node+2];
				key.node_id = node;

				bool nodel = sc_hash_array_lookup(hash_b_mat, &key, &position);
				if(nodel){
					aux[ncount] = ino;
					//VertexEdgeMap[ino][0];
					//VertexEdgeMap[ino][1];
					//VertexEdgeMap[ino][2];
					ncount++;
				}
			}

			int id;
			if(ncount == 1 && false){

				elem->edge[VertexEdgeMap[aux[0]][0]].ref = true;
				elem->edge[VertexEdgeMap[aux[0]][1]].ref = true;
				elem->edge[VertexEdgeMap[aux[0]][2]].ref = true;

				id = elem->edge[VertexEdgeMap[aux[0]][0]].id;
				edge_add(id, hash_edge_ref);

				id = elem->edge[VertexEdgeMap[aux[0]][1]].id;
				edge_add(id, hash_edge_ref);

				id = elem->edge[VertexEdgeMap[aux[0]][2]].id;
				edge_add(id, hash_edge_ref);
			}

			if(ncount != 1){
				std::vector<int> listaux;
				std::vector<int> list;
				//printf("vetor auxiliar: ");
				for(int ive = 0; ive < ncount; ive++){
					for(int iedge = 0; iedge < 3; iedge++){
						list.push_back(VertexEdgeMap[aux[ive]][iedge]);
						//printf("%d ", VertexEdgeMap[aux[ive]][iedge]);
					}
				}
				//printf("\n");

				for(int i = 0; i < list.size(); i++){
					for(int j = i + 1; j < list.size(); j++){
						if(list[i] == list[j]){
							listaux.push_back(list[i]);
						}
					}
				}

				//printf("vetor: ");
				//for(int i = 0; i < list.size(); i++){
				//	printf("%d ", list[i]);
				//}
				//printf("\n");

				for(int i = 0; i < list.size(); i++){
					bool aux = true;
					for(int j = 0; j < listaux.size(); j++){
						if(list[i]==listaux[j]){
							aux = false;
						}
					}

					if(aux){
						elem->edge[list[i]].ref = true;
						id = elem->edge[list[i]].id;
						edge_add(id, hash_edge_ref);
					}
				}
			}

			//printf("numero de nos do elemento: %d eh %d\n",elem->id,ncount);

		}

	}

}


void IdentifyTemplate(hexa_tree_t* mesh, std::vector<int>& elements_ids, sc_hash_array_t* hash_edge_ref){

	int el_1 = 0;
	int el_2 = 0;
	int el_3 = 0;
	int el_4 = 0;
	int el_5 = 0;
	int el_6 = 0;
	int el_7 = 0;
	int el_8 = 0;
	int el_9 = 0;
	int el_10 = 0;
	int el_11 = 0;

	//update the element information
	for (int iel = 0; iel < elements_ids.size(); ++iel)
	{
		size_t position;
		octant_t *elem = (octant_t*) sc_array_index(&mesh->elements, elements_ids[iel]);
		for(int iedge = 0; iedge < 12; iedge++)
		{
			octant_edge_t key;
			key.id = elem->edge[iedge].id;
			bool out = sc_hash_array_lookup(hash_edge_ref, &key, &position);
			elem->edge[iedge].ref = out;
		}
	}

	for (int iel = 0; iel < elements_ids.size(); ++iel) {

		octant_t *elem = (octant_t*) sc_array_index(&mesh->elements, elements_ids[iel]);
		int ed_cont = 0;

		for (int edge = 0; edge < 12; ++edge) {
			if(elem->edge[edge].ref){
				ed_cont++;
			}
		}

		//template 11
		elem->pad = 144;
		elem->tem = 11;

		// element verification
		// template 1
		if(ed_cont == 1){
			if ( (elem->edge[0].ref) && (!elem->edge[1].ref) && (!elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (!elem->edge[5].ref) && (!elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (!elem->edge[9].ref) && (!elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 10;
				elem->tem = 1;

			}
			else if ( (!elem->edge[0].ref) && (elem->edge[1].ref) && (!elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (!elem->edge[5].ref) && (!elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (!elem->edge[9].ref) && (!elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 11;
				elem->tem = 1;

			}
			else if ( (!elem->edge[0].ref) && (!elem->edge[1].ref) && (elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (!elem->edge[5].ref) && (!elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (!elem->edge[9].ref) && (!elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 12;
				elem->tem = 1;

			}
			else if ( (!elem->edge[0].ref) && (!elem->edge[1].ref) && (!elem->edge[2].ref) && (elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (!elem->edge[5].ref) && (!elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (!elem->edge[9].ref) && (!elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 13;
				elem->tem = 1;

			}
			else if ( (!elem->edge[0].ref) && (!elem->edge[1].ref) && (!elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(elem->edge[4].ref) && (!elem->edge[5].ref) && (!elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (!elem->edge[9].ref) && (!elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 14;
				elem->tem = 1;

			}
			else if ( (!elem->edge[0].ref) && (!elem->edge[1].ref) && (!elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (elem->edge[5].ref) && (!elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (!elem->edge[9].ref) && (!elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 15;
				elem->tem = 1;

			}
			else if ( (!elem->edge[0].ref) && (!elem->edge[1].ref) && (!elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (!elem->edge[5].ref) && (elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (!elem->edge[9].ref) && (!elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 16;
				elem->tem = 1;

			}
			else if ( (!elem->edge[0].ref) && (!elem->edge[1].ref) && (!elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (!elem->edge[5].ref) && (!elem->edge[6].ref) && (elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (!elem->edge[9].ref) && (!elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 17;
				elem->tem = 1;

			}
			else if ( (!elem->edge[0].ref) && (!elem->edge[1].ref) && (!elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (!elem->edge[5].ref) && (!elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(elem->edge[8].ref) && (!elem->edge[9].ref) && (!elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 18;
				elem->tem = 1;

			}
			else if ( (!elem->edge[0].ref) && (!elem->edge[1].ref) && (!elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (!elem->edge[5].ref) && (!elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (elem->edge[9].ref) && (!elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 19;
				elem->tem = 1;

			}
			else if ( (!elem->edge[0].ref) && (!elem->edge[1].ref) && (!elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (!elem->edge[5].ref) && (!elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (!elem->edge[9].ref) && (elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 20;
				elem->tem = 1;

			}
			else if ( (!elem->edge[0].ref) && (!elem->edge[1].ref) && (!elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (!elem->edge[5].ref) && (!elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (!elem->edge[9].ref) && (!elem->edge[10].ref) && (elem->edge[11].ref)) {

				elem->pad = 21;
				elem->tem = 1;

			}
		}
		//template 2
		if(ed_cont == 4){
			if ( (elem->edge[0].ref) && (!elem->edge[1].ref) && (elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (!elem->edge[5].ref) && (!elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(elem->edge[8].ref) && (!elem->edge[9].ref) && (elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 22;
				elem->tem = 2;

			}
			else if ( (!elem->edge[0].ref) && (elem->edge[1].ref) && (!elem->edge[2].ref) && (elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (!elem->edge[5].ref) && (!elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (elem->edge[9].ref) && (!elem->edge[10].ref) && (elem->edge[11].ref)) {

				elem->pad = 23;
				elem->tem = 2;

			}
			else if ( (!elem->edge[0].ref) && (!elem->edge[1].ref) && (!elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(elem->edge[4].ref) && (elem->edge[5].ref) && (elem->edge[6].ref) && (elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (!elem->edge[9].ref) && (!elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 24;
				elem->tem = 2;

			}
		}
		//template 3
		if(ed_cont == 2){
			if ( (elem->edge[0].ref) && (!elem->edge[1].ref) && (elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (!elem->edge[5].ref) && (!elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (!elem->edge[9].ref) && (!elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 25;
				elem->tem = 3;

			}
			else if ( (!elem->edge[0].ref) && (elem->edge[1].ref) && (!elem->edge[2].ref) && (elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (!elem->edge[5].ref) && (!elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (!elem->edge[9].ref) && (!elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 26;
				elem->tem = 3;

			}
			else if ( (!elem->edge[0].ref) && (!elem->edge[1].ref) && (!elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (!elem->edge[5].ref) && (!elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(elem->edge[8].ref) && (!elem->edge[9].ref) && (elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 27;
				elem->tem = 3;

			}
			else if ( (!elem->edge[0].ref) && (!elem->edge[1].ref) && (!elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (!elem->edge[5].ref) && (!elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (elem->edge[9].ref) && (!elem->edge[10].ref) && (elem->edge[11].ref)) {

				elem->pad = 28;
				elem->tem = 3;

			}
			else if ( (!elem->edge[0].ref) && (!elem->edge[1].ref) && (!elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(elem->edge[4].ref) && (elem->edge[5].ref) && (!elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (!elem->edge[9].ref) && (!elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 29;
				elem->tem = 3;

			}
			else if ( (!elem->edge[0].ref) && (!elem->edge[1].ref) && (!elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (elem->edge[5].ref) && (elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (!elem->edge[9].ref) && (!elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 30;
				elem->tem = 3;

			}
			else if ( (!elem->edge[0].ref) && (!elem->edge[1].ref) && (!elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (!elem->edge[5].ref) && (elem->edge[6].ref) && (elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (!elem->edge[9].ref) && (!elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 31;
				elem->tem = 3;

			}
			else if ( (!elem->edge[0].ref) && (!elem->edge[1].ref) && (!elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(elem->edge[4].ref) && (!elem->edge[5].ref) && (!elem->edge[6].ref) && (elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (!elem->edge[9].ref) && (!elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 32;
				elem->tem = 3;

			}
			else if ( (elem->edge[0].ref) && (!elem->edge[1].ref) && (!elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (!elem->edge[5].ref) && (!elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(elem->edge[8].ref) && (!elem->edge[9].ref) && (!elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 33;
				elem->tem = 3;

			}
			else if ( (!elem->edge[0].ref) && (elem->edge[1].ref) && (!elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (!elem->edge[5].ref) && (!elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (elem->edge[9].ref) && (!elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 34;
				elem->tem = 3;

			}
			else if ( (!elem->edge[0].ref) && (!elem->edge[1].ref) && (elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (!elem->edge[5].ref) && (!elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (!elem->edge[9].ref) && (elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 35;
				elem->tem = 3;

			}
			else if ( (!elem->edge[0].ref) && (!elem->edge[1].ref) && (!elem->edge[2].ref) && (elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (!elem->edge[5].ref) && (!elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (!elem->edge[9].ref) && (!elem->edge[10].ref) && (elem->edge[11].ref)) {

				elem->pad = 36;
				elem->tem = 3;

			}
		}
		//template 4
		if(ed_cont == 4){
			if ( (elem->edge[0].ref) && (elem->edge[1].ref) && (!elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (!elem->edge[5].ref) && (!elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(elem->edge[8].ref) && (elem->edge[9].ref) && (!elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 37;
				elem->tem = 4;

			}
			else if ( (!elem->edge[0].ref) && (elem->edge[1].ref) && (elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (!elem->edge[5].ref) && (!elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (elem->edge[9].ref) && (elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 38;
				elem->tem = 4;

			}
			else if ( (!elem->edge[0].ref) && (!elem->edge[1].ref) && (elem->edge[2].ref) && (elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (!elem->edge[5].ref) && (!elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (!elem->edge[9].ref) && (elem->edge[10].ref) && (elem->edge[11].ref)) {

				elem->pad = 39;
				elem->tem = 4;

			}
			else if ( (elem->edge[0].ref) && (!elem->edge[1].ref) && (!elem->edge[2].ref) && (elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (!elem->edge[5].ref) && (!elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(elem->edge[8].ref) && (!elem->edge[9].ref) && (!elem->edge[10].ref) && (elem->edge[11].ref)) {

				elem->pad = 40;
				elem->tem = 4;

			}
			else if ( (elem->edge[0].ref) && (!elem->edge[1].ref) && (elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(elem->edge[4].ref) && (!elem->edge[5].ref) && (!elem->edge[6].ref) && (elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (!elem->edge[9].ref) && (!elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 41;
				elem->tem = 4;

			}
			else if ( (elem->edge[0].ref) && (!elem->edge[1].ref) && (elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (elem->edge[5].ref) && (elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (!elem->edge[9].ref) && (!elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 42;
				elem->tem = 4;

			}
			else if ( (!elem->edge[0].ref) && (!elem->edge[1].ref) && (!elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (elem->edge[5].ref) && (elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(elem->edge[8].ref) && (!elem->edge[9].ref) && (elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 43;
				elem->tem = 4;

			}
			else if ( (!elem->edge[0].ref) && (!elem->edge[1].ref) && (!elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(elem->edge[4].ref) && (!elem->edge[5].ref) && (!elem->edge[6].ref) && (elem->edge[7].ref) &&
					(elem->edge[8].ref) && (!elem->edge[9].ref) && (elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 44;
				elem->tem = 4;

			}
			else if ( (!elem->edge[0].ref) && (elem->edge[1].ref) && (!elem->edge[2].ref) && (elem->edge[3].ref) &&
					(elem->edge[4].ref) && (elem->edge[5].ref) && (!elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (!elem->edge[9].ref) && (!elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 45;
				elem->tem = 4;

			}
			else if ( (!elem->edge[0].ref) && (elem->edge[1].ref) && (!elem->edge[2].ref) && (elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (!elem->edge[5].ref) && (elem->edge[6].ref) && (elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (!elem->edge[9].ref) && (!elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 46;
				elem->tem = 4;

			}
			else if ( (!elem->edge[0].ref) && (!elem->edge[1].ref) && (!elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (!elem->edge[5].ref) && (elem->edge[6].ref) && (elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (elem->edge[9].ref) && (!elem->edge[10].ref) && (elem->edge[11].ref)) {

				elem->pad = 47;
				elem->tem = 4;

			}
			else if ( (!elem->edge[0].ref) && (!elem->edge[1].ref) && (!elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(elem->edge[4].ref) && (elem->edge[5].ref) && (!elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (elem->edge[9].ref) && (!elem->edge[10].ref) && (elem->edge[11].ref)) {

				elem->pad = 48;
				elem->tem = 4;

			}
		}
		//template 5
		if(ed_cont == 4 || ed_cont == 3 || ed_cont == 2){
			//template 5 com 4
			if ( (elem->edge[0].ref) && (elem->edge[1].ref) && (elem->edge[2].ref) && (elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (!elem->edge[5].ref) && (!elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (!elem->edge[9].ref) && (!elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 49;
				elem->tem = 5;

			}
			else if ( (!elem->edge[0].ref) && (!elem->edge[1].ref) && (!elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (!elem->edge[5].ref) && (!elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(elem->edge[8].ref) && (elem->edge[9].ref) && (elem->edge[10].ref) && (elem->edge[11].ref)) {

				elem->pad = 50;
				elem->tem = 5;

			}
			else if ( (elem->edge[0].ref) && (!elem->edge[1].ref) && (!elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(elem->edge[4].ref) && (elem->edge[5].ref) && (!elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(elem->edge[8].ref) && (!elem->edge[9].ref) && (!elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 51;
				elem->tem = 5;

			}
			else if ( (!elem->edge[0].ref) && (elem->edge[1].ref) && (!elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (elem->edge[5].ref) && (elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (elem->edge[9].ref) && (!elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 52;
				elem->tem = 5;

			}
			else if ( (!elem->edge[0].ref) && (!elem->edge[1].ref) && (elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (!elem->edge[5].ref) && (elem->edge[6].ref) && (elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (!elem->edge[9].ref) && (elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 53;
				elem->tem = 5;

			}
			else if ( (!elem->edge[0].ref) && (!elem->edge[1].ref) && (!elem->edge[2].ref) && (elem->edge[3].ref) &&
					(elem->edge[4].ref) && (!elem->edge[5].ref) && (!elem->edge[6].ref) && (elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (!elem->edge[9].ref) && (!elem->edge[10].ref) && (elem->edge[11].ref)) {

				elem->pad = 54;
				elem->tem = 5;

			}


			//template 5 com 3
			else if ( (elem->edge[0].ref) && (elem->edge[1].ref) && (!elem->edge[2].ref) && (elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (!elem->edge[5].ref) && (!elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (!elem->edge[9].ref) && (!elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 55;
				elem->tem = 5;

			}
			else if ( (elem->edge[0].ref) && (elem->edge[1].ref) && (elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (!elem->edge[5].ref) && (!elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (!elem->edge[9].ref) && (!elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 56;
				elem->tem = 5;

			}
			else if ( (!elem->edge[0].ref) && (elem->edge[1].ref) && (elem->edge[2].ref) && (elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (!elem->edge[5].ref) && (!elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (!elem->edge[9].ref) && (!elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 57;
				elem->tem = 5;

			}
			else if ( (elem->edge[0].ref) && (!elem->edge[1].ref) && (elem->edge[2].ref) && (elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (!elem->edge[5].ref) && (!elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (!elem->edge[9].ref) && (!elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 58;
				elem->tem = 5;

			}
			else if ( (!elem->edge[0].ref) && (!elem->edge[1].ref) && (!elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (!elem->edge[5].ref) && (!elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(elem->edge[8].ref) && (elem->edge[9].ref) && (elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 59;
				elem->tem = 5;

			}
			else if ( (!elem->edge[0].ref) && (!elem->edge[1].ref) && (!elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (!elem->edge[5].ref) && (!elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(elem->edge[8].ref) && (elem->edge[9].ref) && (!elem->edge[10].ref) && (elem->edge[11].ref)) {

				elem->pad = 60;
				elem->tem = 5;

			}
			else if ( (!elem->edge[0].ref) && (!elem->edge[1].ref) && (!elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (!elem->edge[5].ref) && (!elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (elem->edge[9].ref) && (elem->edge[10].ref) && (elem->edge[11].ref)) {

				elem->pad = 61;
				elem->tem = 5;

			}
			else if ( (!elem->edge[0].ref) && (!elem->edge[1].ref) && (!elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (!elem->edge[5].ref) && (!elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(elem->edge[8].ref) && (!elem->edge[9].ref) && (elem->edge[10].ref) && (elem->edge[11].ref)) {

				elem->pad = 62;
				elem->tem = 5;

			}
			else if ( (elem->edge[0].ref) && (!elem->edge[1].ref) && (!elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(elem->edge[4].ref) && (elem->edge[5].ref) && (!elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (!elem->edge[9].ref) && (!elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 63;
				elem->tem = 5;

			}
			else if ( (elem->edge[0].ref) && (!elem->edge[1].ref) && (!elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (elem->edge[5].ref) && (!elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(elem->edge[8].ref) && (!elem->edge[9].ref) && (!elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 64;
				elem->tem = 5;

			}
			else if ( (!elem->edge[0].ref) && (!elem->edge[1].ref) && (!elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(elem->edge[4].ref) && (elem->edge[5].ref) && (!elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(elem->edge[8].ref) && (!elem->edge[9].ref) && (!elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 65;
				elem->tem = 5;

			}
			else if ( (elem->edge[0].ref) && (!elem->edge[1].ref) && (!elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(elem->edge[4].ref) && (!elem->edge[5].ref) && (!elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(elem->edge[8].ref) && (!elem->edge[9].ref) && (!elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 66;
				elem->tem = 5;

			}
			else if ( (!elem->edge[0].ref) && (elem->edge[1].ref) && (!elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (elem->edge[5].ref) && (elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (!elem->edge[9].ref) && (!elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 67;
				elem->tem = 5;

			}
			else if ( (!elem->edge[0].ref) && (elem->edge[1].ref) && (!elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (elem->edge[5].ref) && (!elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (elem->edge[9].ref) && (!elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 68;
				elem->tem = 5;

			}
			else if ( (!elem->edge[0].ref) && (!elem->edge[1].ref) && (!elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (elem->edge[5].ref) && (elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (elem->edge[9].ref) && (!elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 69;
				elem->tem = 5;

			}
			else if ( (!elem->edge[0].ref) && (elem->edge[1].ref) && (!elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (!elem->edge[5].ref) && (elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (elem->edge[9].ref) && (!elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 70;
				elem->tem = 5;

			}
			else if ( (!elem->edge[0].ref) && (!elem->edge[1].ref) && (elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (!elem->edge[5].ref) && (elem->edge[6].ref) && (elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (!elem->edge[9].ref) && (!elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 71;
				elem->tem = 5;

			}
			else if ( (!elem->edge[0].ref) && (!elem->edge[1].ref) && (elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (!elem->edge[5].ref) && (elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (!elem->edge[9].ref) && (elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 72;
				elem->tem = 5;

			}
			else if ( (!elem->edge[0].ref) && (!elem->edge[1].ref) && (!elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (!elem->edge[5].ref) && (elem->edge[6].ref) && (elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (!elem->edge[9].ref) && (elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 73;
				elem->tem = 5;

			}
			else if ( (!elem->edge[0].ref) && (!elem->edge[1].ref) && (elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (!elem->edge[5].ref) && (!elem->edge[6].ref) && (elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (!elem->edge[9].ref) && (elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 74;
				elem->tem = 5;

			}
			else if ( (!elem->edge[0].ref) && (!elem->edge[1].ref) && (!elem->edge[2].ref) && (elem->edge[3].ref) &&
					(elem->edge[4].ref) && (!elem->edge[5].ref) && (!elem->edge[6].ref) && (elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (!elem->edge[9].ref) && (!elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 75;
				elem->tem = 5;

			}
			else if ( (!elem->edge[0].ref) && (!elem->edge[1].ref) && (!elem->edge[2].ref) && (elem->edge[3].ref) &&
					(elem->edge[4].ref) && (!elem->edge[5].ref) && (!elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (!elem->edge[9].ref) && (!elem->edge[10].ref) && (elem->edge[11].ref)) {

				elem->pad = 76;
				elem->tem = 5;

			}
			else if ( (!elem->edge[0].ref) && (!elem->edge[1].ref) && (!elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(elem->edge[4].ref) && (!elem->edge[5].ref) && (!elem->edge[6].ref) && (elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (!elem->edge[9].ref) && (!elem->edge[10].ref) && (elem->edge[11].ref)) {

				elem->pad = 77;
				elem->tem = 5;

			}
			else if ( (!elem->edge[0].ref) && (!elem->edge[1].ref) && (!elem->edge[2].ref) && (elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (!elem->edge[5].ref) && (!elem->edge[6].ref) && (elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (!elem->edge[9].ref) && (!elem->edge[10].ref) && (elem->edge[11].ref)) {

				elem->pad = 78;
				elem->tem = 5;

			}

			// templete 5 com 2
			else if ( (elem->edge[0].ref) && (elem->edge[1].ref) && (!elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (!elem->edge[5].ref) && (!elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (!elem->edge[9].ref) && (!elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 79;
				elem->tem = 5;

			}
			else if ( (!elem->edge[0].ref) && (elem->edge[1].ref) && (elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (!elem->edge[5].ref) && (!elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (!elem->edge[9].ref) && (!elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 80;
				elem->tem = 5;

			}
			else if ( (!elem->edge[0].ref) && (!elem->edge[1].ref) && (elem->edge[2].ref) && (elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (!elem->edge[5].ref) && (!elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (!elem->edge[9].ref) && (!elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 81;
				elem->tem = 5;

			}
			else if ( (elem->edge[0].ref) && (!elem->edge[1].ref) && (!elem->edge[2].ref) && (elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (!elem->edge[5].ref) && (!elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (!elem->edge[9].ref) && (!elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 82;
				elem->tem = 5;

			}
			else if ( (!elem->edge[0].ref) && (!elem->edge[1].ref) && (!elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (!elem->edge[5].ref) && (!elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(elem->edge[8].ref) && (elem->edge[9].ref) && (!elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 83;
				elem->tem = 5;

			}
			else if ( (!elem->edge[0].ref) && (!elem->edge[1].ref) && (!elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (!elem->edge[5].ref) && (!elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (elem->edge[9].ref) && (elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 84;
				elem->tem = 5;

			}
			else if ( (!elem->edge[0].ref) && (!elem->edge[1].ref) && (!elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (!elem->edge[5].ref) && (!elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (!elem->edge[9].ref) && (elem->edge[10].ref) && (elem->edge[11].ref)) {

				elem->pad = 85;
				elem->tem = 5;

			}
			else if ( (!elem->edge[0].ref) && (!elem->edge[1].ref) && (!elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (!elem->edge[5].ref) && (!elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(elem->edge[8].ref) && (!elem->edge[9].ref) && (!elem->edge[10].ref) && (elem->edge[11].ref)) {

				elem->pad = 86;
				elem->tem = 5;

			}
			else if ( (elem->edge[0].ref) && (!elem->edge[1].ref) && (!elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(elem->edge[4].ref) && (!elem->edge[5].ref) && (!elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (!elem->edge[9].ref) && (!elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 87;
				elem->tem = 5;

			}
			else if ( (elem->edge[0].ref) && (!elem->edge[1].ref) && (!elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (elem->edge[5].ref) && (!elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (!elem->edge[9].ref) && (!elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 88;
				elem->tem = 5;

			}
			else if ( (!elem->edge[0].ref) && (!elem->edge[1].ref) && (!elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (elem->edge[5].ref) && (!elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(elem->edge[8].ref) && (!elem->edge[9].ref) && (!elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 89;
				elem->tem = 5;

			}
			else if ( (!elem->edge[0].ref) && (!elem->edge[1].ref) && (!elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(elem->edge[4].ref) && (!elem->edge[5].ref) && (!elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(elem->edge[8].ref) && (!elem->edge[9].ref) && (!elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 90;
				elem->tem = 5;

			}
			else if ( (!elem->edge[0].ref) && (!elem->edge[1].ref) && (!elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (elem->edge[5].ref) && (!elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (elem->edge[9].ref) && (!elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 91;
				elem->tem = 5;

			}
			else if ( (!elem->edge[0].ref) && (!elem->edge[1].ref) && (!elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (!elem->edge[5].ref) && (elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (elem->edge[9].ref) && (!elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 92;
				elem->tem = 5;

			}
			else if ( (!elem->edge[0].ref) && (elem->edge[1].ref) && (!elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (!elem->edge[5].ref) && (elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (!elem->edge[9].ref) && (!elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 93;
				elem->tem = 5;

			}
			else if ( (!elem->edge[0].ref) && (elem->edge[1].ref) && (!elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (elem->edge[5].ref) && (!elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (!elem->edge[9].ref) && (!elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 94;
				elem->tem = 5;

			}
			else if ( (!elem->edge[0].ref) && (!elem->edge[1].ref) && (elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (!elem->edge[5].ref) && (elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (!elem->edge[9].ref) && (!elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 95;
				elem->tem = 5;

			}
			else if ( (!elem->edge[0].ref) && (!elem->edge[1].ref) && (!elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (!elem->edge[5].ref) && (elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (!elem->edge[9].ref) && (elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 96;
				elem->tem = 5;

			}
			else if ( (!elem->edge[0].ref) && (!elem->edge[1].ref) && (!elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (!elem->edge[5].ref) && (!elem->edge[6].ref) && (elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (!elem->edge[9].ref) && (elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 97;
				elem->tem = 5;

			}
			else if ( (!elem->edge[0].ref) && (!elem->edge[1].ref) && (elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (!elem->edge[5].ref) && (!elem->edge[6].ref) && (elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (!elem->edge[9].ref) && (!elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 98;
				elem->tem = 5;

			}
			else if ( (!elem->edge[0].ref) && (!elem->edge[1].ref) && (!elem->edge[2].ref) && (elem->edge[3].ref) &&
					(elem->edge[4].ref) && (!elem->edge[5].ref) && (!elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (!elem->edge[9].ref) && (!elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 99;
				elem->tem = 5;

			}
			else if ( (!elem->edge[0].ref) && (!elem->edge[1].ref) && (!elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(elem->edge[4].ref) && (!elem->edge[5].ref) && (elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (!elem->edge[9].ref) && (!elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 100;
				elem->tem = 5;

			}
			else if ( (!elem->edge[0].ref) && (!elem->edge[1].ref) && (!elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (!elem->edge[5].ref) && (elem->edge[6].ref) && (elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (!elem->edge[9].ref) && (!elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 101;
				elem->tem = 5;

			}
			else if ( (!elem->edge[0].ref) && (!elem->edge[1].ref) && (!elem->edge[2].ref) && (elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (!elem->edge[5].ref) && (!elem->edge[6].ref) && (elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (!elem->edge[9].ref) && (!elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 102;
				elem->tem = 5;

			}
		}
		//template 6
		if(ed_cont == 6){
			if ( (elem->edge[0].ref) && (elem->edge[1].ref) && (!elem->edge[2].ref) && (elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (!elem->edge[5].ref) && (!elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(elem->edge[8].ref) && (elem->edge[9].ref) && (!elem->edge[10].ref) && (elem->edge[11].ref)) {

				elem->pad = 103;
				elem->tem = 6;

			}
			else if ( (elem->edge[0].ref) && (elem->edge[1].ref) && (elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (!elem->edge[5].ref) && (!elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(elem->edge[8].ref) && (elem->edge[9].ref) && (elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 104;
				elem->tem = 6;

			}
			else if ( (!elem->edge[0].ref) && (elem->edge[1].ref) && (elem->edge[2].ref) && (elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (!elem->edge[5].ref) && (!elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (elem->edge[9].ref) && (elem->edge[10].ref) && (elem->edge[11].ref)) {

				elem->pad = 105;
				elem->tem = 6;

			}
			else if ( (elem->edge[0].ref) && (!elem->edge[1].ref) && (elem->edge[2].ref) && (elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (!elem->edge[5].ref) && (!elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(elem->edge[8].ref) && (!elem->edge[9].ref) && (elem->edge[10].ref) && (elem->edge[11].ref)) {

				elem->pad = 106;
				elem->tem = 6;

			}
			else if ( (elem->edge[0].ref) && (!elem->edge[1].ref) && (elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (elem->edge[5].ref) && (elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(elem->edge[8].ref) && (!elem->edge[9].ref) && (elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 107;
				elem->tem = 6;

			}
			else if ( (!elem->edge[0].ref) && (!elem->edge[1].ref) && (!elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(elem->edge[4].ref) && (elem->edge[5].ref) && (elem->edge[6].ref) && (elem->edge[7].ref) &&
					(elem->edge[8].ref) && (!elem->edge[9].ref) && (elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 108;
				elem->tem = 6;

			}
			else if ( (elem->edge[0].ref) && (!elem->edge[1].ref) && (elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(elem->edge[4].ref) && (!elem->edge[5].ref) && (!elem->edge[6].ref) && (elem->edge[7].ref) &&
					(elem->edge[8].ref) && (!elem->edge[9].ref) && (elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 109;
				elem->tem = 6;

			}
			else if ( (elem->edge[0].ref) && (!elem->edge[1].ref) && (elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(elem->edge[4].ref) && (elem->edge[5].ref) && (elem->edge[6].ref) && (elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (!elem->edge[9].ref) && (!elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 110;
				elem->tem = 6;

			}
			else if ( (!elem->edge[0].ref) && (elem->edge[1].ref) && (!elem->edge[2].ref) && (elem->edge[3].ref) &&
					(elem->edge[4].ref) && (elem->edge[5].ref) && (elem->edge[6].ref) && (elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (!elem->edge[9].ref) && (!elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 111;
				elem->tem = 6;

			}
			else if ( (!elem->edge[0].ref) && (elem->edge[1].ref) && (!elem->edge[2].ref) && (elem->edge[3].ref) &&
					(elem->edge[4].ref) && (elem->edge[5].ref) && (!elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (elem->edge[9].ref) && (!elem->edge[10].ref) && (elem->edge[11].ref)) {

				elem->pad = 112;
				elem->tem = 6;

			}
			else if ( (!elem->edge[0].ref) && (!elem->edge[1].ref) && (!elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(elem->edge[4].ref) && (elem->edge[5].ref) && (elem->edge[6].ref) && (elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (elem->edge[9].ref) && (!elem->edge[10].ref) && (elem->edge[11].ref)) {

				elem->pad = 113;
				elem->tem = 6;

			}
			else if ( (!elem->edge[0].ref) && (elem->edge[1].ref) && (!elem->edge[2].ref) && (elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (!elem->edge[5].ref) && (elem->edge[6].ref) && (elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (elem->edge[9].ref) && (!elem->edge[10].ref) && (elem->edge[11].ref)) {

				elem->pad = 114;
				elem->tem = 6;

			}
		}
		//template 7
		if(ed_cont == 3){
			if ( (elem->edge[0].ref) && (!elem->edge[1].ref) && (!elem->edge[2].ref) && (elem->edge[3].ref) &&
					(elem->edge[4].ref) && (!elem->edge[5].ref) && (!elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (!elem->edge[9].ref) && (!elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 115;
				elem->tem = 7;

			}
			else if ( (elem->edge[0].ref) && (elem->edge[1].ref) && (!elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (elem->edge[5].ref) && (!elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (!elem->edge[9].ref) && (!elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 116;
				elem->tem = 7;

			}
			else if ( (!elem->edge[0].ref) && (elem->edge[1].ref) && (elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (!elem->edge[5].ref) && (elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (!elem->edge[9].ref) && (!elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 117;
				elem->tem = 7;

			}
			else if ( (!elem->edge[0].ref) && (!elem->edge[1].ref) && (elem->edge[2].ref) && (elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (!elem->edge[5].ref) && (!elem->edge[6].ref) && (elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (!elem->edge[9].ref) && (!elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 118;
				elem->tem = 7;

			}
			else if ( (!elem->edge[0].ref) && (!elem->edge[1].ref) && (!elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(elem->edge[4].ref) && (!elem->edge[5].ref) && (!elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(elem->edge[8].ref) && (!elem->edge[9].ref) && (!elem->edge[10].ref) && (elem->edge[11].ref)) {

				elem->pad = 119;
				elem->tem = 7;

			}
			else if ( (!elem->edge[0].ref) && (!elem->edge[1].ref) && (!elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (elem->edge[5].ref) && (!elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(elem->edge[8].ref) && (elem->edge[9].ref) && (!elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 120;
				elem->tem = 7;

			}
			else if ( (!elem->edge[0].ref) && (!elem->edge[1].ref) && (!elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (!elem->edge[5].ref) && (elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (elem->edge[9].ref) && (elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 121;
				elem->tem = 7;

			}
			else if ( (!elem->edge[0].ref) && (!elem->edge[1].ref) && (!elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (!elem->edge[5].ref) && (!elem->edge[6].ref) && (elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (!elem->edge[9].ref) && (elem->edge[10].ref) && (elem->edge[11].ref)) {

				elem->pad = 122;
				elem->tem = 7;

			}
		}
		//template 8
		if(ed_cont == 8){
			if ( (elem->edge[0].ref) && (elem->edge[1].ref) && (elem->edge[2].ref) && (elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (!elem->edge[5].ref) && (!elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(elem->edge[8].ref) && (elem->edge[9].ref) && (elem->edge[10].ref) && (elem->edge[11].ref)) {

				elem->pad = 123;
				elem->tem = 8;

			}
			else if ( (elem->edge[0].ref) && (!elem->edge[1].ref) && (elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(elem->edge[4].ref) && (elem->edge[5].ref) && (elem->edge[6].ref) && (elem->edge[7].ref) &&
					(elem->edge[8].ref) && (!elem->edge[9].ref) && (elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 124;
				elem->tem = 8;

			}
			else if ( (!elem->edge[0].ref) && (elem->edge[1].ref) && (!elem->edge[2].ref) && (elem->edge[3].ref) &&
					(elem->edge[4].ref) && (elem->edge[5].ref) && (elem->edge[6].ref) && (elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (elem->edge[9].ref) && (!elem->edge[10].ref) && (elem->edge[11].ref)) {

				elem->pad = 125;
				elem->tem = 8;

			}
		}
		//template 9
		if(ed_cont == 8){
			if ( (elem->edge[0].ref) && (elem->edge[1].ref) && (elem->edge[2].ref) && (elem->edge[3].ref) &&
					(elem->edge[4].ref) && (elem->edge[5].ref) && (elem->edge[6].ref) && (elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (!elem->edge[9].ref) && (!elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 126;
				elem->tem = 9;

			}
			else if ( (!elem->edge[0].ref) && (!elem->edge[1].ref) && (!elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(elem->edge[4].ref) && (elem->edge[5].ref) && (elem->edge[6].ref) && (elem->edge[7].ref) &&
					(elem->edge[8].ref) && (elem->edge[9].ref) && (elem->edge[10].ref) && (elem->edge[11].ref)) {

				elem->pad = 127;
				elem->tem = 9;

			}
			else if ( (elem->edge[0].ref) && (elem->edge[1].ref) && (!elem->edge[2].ref) && (elem->edge[3].ref) &&
					(elem->edge[4].ref) && (elem->edge[5].ref) && (!elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(elem->edge[8].ref) && (elem->edge[9].ref) && (!elem->edge[10].ref) && (elem->edge[11].ref)) {

				elem->pad = 128;
				elem->tem = 9;

			}
			else if ( (!elem->edge[0].ref) && (elem->edge[1].ref) && (elem->edge[2].ref) && (elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (!elem->edge[5].ref) && (elem->edge[6].ref) && (elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (elem->edge[9].ref) && (elem->edge[10].ref) && (elem->edge[11].ref)) {

				elem->pad = 129;
				elem->tem = 9;

			}
			else if ( (elem->edge[0].ref) && (elem->edge[1].ref) && (elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (elem->edge[5].ref) && (elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(elem->edge[8].ref) && (elem->edge[9].ref) && (elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 130;
				elem->tem = 9;

			}
			else if ( (elem->edge[0].ref) && (!elem->edge[1].ref) && (elem->edge[2].ref) && (elem->edge[3].ref) &&
					(elem->edge[4].ref) && (!elem->edge[5].ref) && (!elem->edge[6].ref) && (elem->edge[7].ref) &&
					(elem->edge[8].ref) && (!elem->edge[9].ref) && (elem->edge[10].ref) && (elem->edge[11].ref)) {

				elem->pad = 131;
				elem->tem = 9;

			}
		}
		//template 10
		if(ed_cont == 5){
			if ( (elem->edge[0].ref) && (elem->edge[1].ref) && (!elem->edge[2].ref) && (elem->edge[3].ref) &&
					(elem->edge[4].ref) && (elem->edge[5].ref) && (!elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (!elem->edge[9].ref) && (!elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 132;
				elem->tem = 10;

			}
			else if ( (elem->edge[0].ref) && (elem->edge[1].ref) && (elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (elem->edge[5].ref) && (elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (!elem->edge[9].ref) && (!elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 133;
				elem->tem = 10;

			}
			else if ( (!elem->edge[0].ref) && (elem->edge[1].ref) && (elem->edge[2].ref) && (elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (!elem->edge[5].ref) && (elem->edge[6].ref) && (elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (!elem->edge[9].ref) && (!elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 134;
				elem->tem = 10;

			}
			else if ( (elem->edge[0].ref) && (!elem->edge[1].ref) && (elem->edge[2].ref) && (elem->edge[3].ref) &&
					(elem->edge[4].ref) && (!elem->edge[5].ref) && (!elem->edge[6].ref) && (elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (!elem->edge[9].ref) && (!elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 135;
				elem->tem = 10;

			}
			else if ( (elem->edge[0].ref) && (!elem->edge[1].ref) && (!elem->edge[2].ref) && (elem->edge[3].ref) &&
					(elem->edge[4].ref) && (!elem->edge[5].ref) && (!elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(elem->edge[8].ref) && (!elem->edge[9].ref) && (!elem->edge[10].ref) && (elem->edge[11].ref)) {

				elem->pad = 136;
				elem->tem = 10;

			}
			else if ( (elem->edge[0].ref) && (elem->edge[1].ref) && (!elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (elem->edge[5].ref) && (!elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(elem->edge[8].ref) && (elem->edge[9].ref) && (!elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 137;
				elem->tem = 10;

			}
			else if ( (!elem->edge[0].ref) && (elem->edge[1].ref) && (elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (!elem->edge[5].ref) && (elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (elem->edge[9].ref) && (elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 138;
				elem->tem = 10;

			}
			else if ( (!elem->edge[0].ref) && (!elem->edge[1].ref) && (elem->edge[2].ref) && (elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (!elem->edge[5].ref) && (!elem->edge[6].ref) && (elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (!elem->edge[9].ref) && (elem->edge[10].ref) && (elem->edge[11].ref)) {

				elem->pad = 139;
				elem->tem = 10;

			}
			else if ( (!elem->edge[0].ref) && (!elem->edge[1].ref) && (!elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(elem->edge[4].ref) && (elem->edge[5].ref) && (!elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(elem->edge[8].ref) && (elem->edge[9].ref) && (!elem->edge[10].ref) && (elem->edge[11].ref)) {

				elem->pad = 140;
				elem->tem = 10;

			}
			else if ( (!elem->edge[0].ref) && (!elem->edge[1].ref) && (!elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (elem->edge[5].ref) && (elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(elem->edge[8].ref) && (elem->edge[9].ref) && (elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 141;
				elem->tem = 10;

			}
			else if ( (!elem->edge[0].ref) && (!elem->edge[1].ref) && (!elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (!elem->edge[5].ref) && (elem->edge[6].ref) && (elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (elem->edge[9].ref) && (elem->edge[10].ref) && (elem->edge[11].ref)) {

				elem->pad = 142;
				elem->tem = 10;

			}
			else if ( (!elem->edge[0].ref) && (!elem->edge[1].ref) && (!elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(elem->edge[4].ref) && (!elem->edge[5].ref) && (!elem->edge[6].ref) && (elem->edge[7].ref) &&
					(elem->edge[8].ref) && (!elem->edge[9].ref) && (elem->edge[10].ref) && (elem->edge[11].ref)) {

				elem->pad = 143;
				elem->tem = 10;

			}
		}
		if(elem->tem== 1){el_1++;}
		if(elem->tem== 2){el_2++;}
		if(elem->tem== 3){el_3++;}
		if(elem->tem== 4){el_4++;}
		if(elem->tem== 5){el_5++;}
		if(elem->tem== 6){el_6++;}
		if(elem->tem== 7){el_7++;}
		if(elem->tem== 8){el_8++;}
		if(elem->tem== 9){el_9++;}
		if(elem->tem==10){el_10++;}
		if(elem->tem==11){el_11++;}
	}

	int su = 0;
	su = el_1+el_2+el_3+el_4+el_5+el_6+el_7+el_8+el_9+el_10+el_11;

#ifdef HEXA_DEBUG_
	fprintf(mesh->fdbg,"case_1: %d\n",el_1);
	fprintf(mesh->fdbg,"case_2: %d\n",el_2);
	fprintf(mesh->fdbg,"case_3: %d\n",el_3);
	fprintf(mesh->fdbg,"case_4: %d\n",el_4);
	fprintf(mesh->fdbg,"case_5: %d\n",el_5);
	fprintf(mesh->fdbg,"case_6: %d\n",el_6);
	fprintf(mesh->fdbg,"case_7: %d\n",el_7);
	fprintf(mesh->fdbg,"case_8: %d\n",el_8);
	fprintf(mesh->fdbg,"case_9: %d\n",el_9);
	fprintf(mesh->fdbg,"case_10: %d\n",el_10);
	fprintf(mesh->fdbg,"case_11: %d\n",el_11);
	fprintf(mesh->fdbg,"sum: %d\n",su);
#endif
}

void Edge_identification(hexa_tree_t* mesh, std::vector<int>& elements_ids, sc_hash_array_t* hash_edge_ref)
{

	for (int iel = 0; iel < elements_ids.size(); ++iel)
	{

		octant_t *elem = (octant_t*) sc_array_index(&mesh->elements, elements_ids[iel]);
		//template 11
		if(elem->pad==144)
		{
			for (int edge = 0; edge < 12; ++edge) edge_add(elem->edge[edge].id, hash_edge_ref );
		}
		//template 10
		else if(elem->pad==143)
		{
			edge_add(elem->edge[4].id, hash_edge_ref );
			edge_add(elem->edge[7].id, hash_edge_ref );
			edge_add(elem->edge[8].id, hash_edge_ref );
			edge_add(elem->edge[10].id, hash_edge_ref );
			edge_add(elem->edge[11].id, hash_edge_ref );
		}
		else if(elem->pad==142)
		{
			edge_add(elem->edge[6].id, hash_edge_ref );
			edge_add(elem->edge[7].id, hash_edge_ref );
			edge_add(elem->edge[9].id, hash_edge_ref );
			edge_add(elem->edge[10].id, hash_edge_ref );
			edge_add(elem->edge[11].id, hash_edge_ref );
		}
		else if(elem->pad==141)
		{
			edge_add(elem->edge[5].id, hash_edge_ref );
			edge_add(elem->edge[6].id, hash_edge_ref );
			edge_add(elem->edge[8].id, hash_edge_ref );
			edge_add(elem->edge[9].id, hash_edge_ref );
			edge_add(elem->edge[10].id, hash_edge_ref );
		}
		else if(elem->pad==140)
		{
			edge_add(elem->edge[4].id, hash_edge_ref );
			edge_add(elem->edge[5].id, hash_edge_ref );
			edge_add(elem->edge[8].id, hash_edge_ref );
			edge_add(elem->edge[9].id, hash_edge_ref );
			edge_add(elem->edge[11].id, hash_edge_ref );
		}
		else if(elem->pad==139)
		{
			edge_add(elem->edge[2].id, hash_edge_ref );
			edge_add(elem->edge[3].id, hash_edge_ref );
			edge_add(elem->edge[7].id, hash_edge_ref );
			edge_add(elem->edge[10].id, hash_edge_ref );
			edge_add(elem->edge[11].id, hash_edge_ref );
		}
		else if(elem->pad==138)
		{
			edge_add(elem->edge[1].id, hash_edge_ref );
			edge_add(elem->edge[2].id, hash_edge_ref );
			edge_add(elem->edge[6].id, hash_edge_ref );
			edge_add(elem->edge[9].id, hash_edge_ref );
			edge_add(elem->edge[10].id, hash_edge_ref );
		}
		else if(elem->pad==137)
		{
			edge_add(elem->edge[0].id, hash_edge_ref );
			edge_add(elem->edge[1].id, hash_edge_ref );
			edge_add(elem->edge[5].id, hash_edge_ref );
			edge_add(elem->edge[8].id, hash_edge_ref );
			edge_add(elem->edge[9].id, hash_edge_ref );
		}
		else if(elem->pad==136)
		{
			edge_add(elem->edge[0].id, hash_edge_ref );
			edge_add(elem->edge[3].id, hash_edge_ref );
			edge_add(elem->edge[4].id, hash_edge_ref );
			edge_add(elem->edge[8].id, hash_edge_ref );
			edge_add(elem->edge[11].id, hash_edge_ref );
		}
		else if(elem->pad==135)
		{
			edge_add(elem->edge[0].id, hash_edge_ref );
			edge_add(elem->edge[2].id, hash_edge_ref );
			edge_add(elem->edge[3].id, hash_edge_ref );
			edge_add(elem->edge[4].id, hash_edge_ref );
			edge_add(elem->edge[7].id, hash_edge_ref );
		}
		else if(elem->pad==134)
		{
			edge_add(elem->edge[1].id, hash_edge_ref );
			edge_add(elem->edge[2].id, hash_edge_ref );
			edge_add(elem->edge[3].id, hash_edge_ref );
			edge_add(elem->edge[6].id, hash_edge_ref );
			edge_add(elem->edge[7].id, hash_edge_ref );
		}
		else if(elem->pad==133)
		{
			edge_add(elem->edge[0].id, hash_edge_ref );
			edge_add(elem->edge[1].id, hash_edge_ref );
			edge_add(elem->edge[2].id, hash_edge_ref );
			edge_add(elem->edge[5].id, hash_edge_ref );
			edge_add(elem->edge[6].id, hash_edge_ref );
		}
		else if(elem->pad==132)
		{
			edge_add(elem->edge[0].id, hash_edge_ref );
			edge_add(elem->edge[1].id, hash_edge_ref );
			edge_add(elem->edge[3].id, hash_edge_ref );
			edge_add(elem->edge[4].id, hash_edge_ref );
			edge_add(elem->edge[5].id, hash_edge_ref );
		}
		//template 9
		else if(elem->pad==131){
			edge_add(elem->edge[0].id, hash_edge_ref );
			edge_add(elem->edge[2].id, hash_edge_ref );
			edge_add(elem->edge[3].id, hash_edge_ref );
			edge_add(elem->edge[4].id, hash_edge_ref );
			edge_add(elem->edge[7].id, hash_edge_ref );
			edge_add(elem->edge[8].id, hash_edge_ref );
			edge_add(elem->edge[10].id, hash_edge_ref );
			edge_add(elem->edge[11].id, hash_edge_ref );
		}else if(elem->pad==130){
			edge_add(elem->edge[0].id, hash_edge_ref );
			edge_add(elem->edge[1].id, hash_edge_ref );
			edge_add(elem->edge[3].id, hash_edge_ref );
			edge_add(elem->edge[5].id, hash_edge_ref );
			edge_add(elem->edge[6].id, hash_edge_ref );
			edge_add(elem->edge[8].id, hash_edge_ref );
			edge_add(elem->edge[9].id, hash_edge_ref );
			edge_add(elem->edge[11].id, hash_edge_ref );
		}else if(elem->pad==129){
			edge_add(elem->edge[1].id, hash_edge_ref );
			edge_add(elem->edge[2].id, hash_edge_ref );
			edge_add(elem->edge[3].id, hash_edge_ref );
			edge_add(elem->edge[6].id, hash_edge_ref );
			edge_add(elem->edge[7].id, hash_edge_ref );
			edge_add(elem->edge[9].id, hash_edge_ref );
			edge_add(elem->edge[10].id, hash_edge_ref );
			edge_add(elem->edge[11].id, hash_edge_ref );
		}else if(elem->pad==128){
			edge_add(elem->edge[0].id, hash_edge_ref );
			edge_add(elem->edge[1].id, hash_edge_ref );
			edge_add(elem->edge[3].id, hash_edge_ref );
			edge_add(elem->edge[5].id, hash_edge_ref );
			edge_add(elem->edge[8].id, hash_edge_ref );
			edge_add(elem->edge[9].id, hash_edge_ref );
			edge_add(elem->edge[9].id, hash_edge_ref );
			edge_add(elem->edge[11].id, hash_edge_ref );
		}else if(elem->pad==127){
			edge_add(elem->edge[4].id, hash_edge_ref );
			edge_add(elem->edge[5].id, hash_edge_ref );
			edge_add(elem->edge[6].id, hash_edge_ref );
			edge_add(elem->edge[7].id, hash_edge_ref );
			edge_add(elem->edge[8].id, hash_edge_ref );
			edge_add(elem->edge[9].id, hash_edge_ref );
			edge_add(elem->edge[10].id, hash_edge_ref );
			edge_add(elem->edge[11].id, hash_edge_ref );
		}else if(elem->pad==126){
			edge_add(elem->edge[0].id, hash_edge_ref );
			edge_add(elem->edge[1].id, hash_edge_ref );
			edge_add(elem->edge[2].id, hash_edge_ref );
			edge_add(elem->edge[3].id, hash_edge_ref );
			edge_add(elem->edge[4].id, hash_edge_ref );
			edge_add(elem->edge[5].id, hash_edge_ref );
			edge_add(elem->edge[6].id, hash_edge_ref );
			edge_add(elem->edge[7].id, hash_edge_ref );
		}
		//template 8
		else if(elem->pad==125){
			edge_add(elem->edge[1].id, hash_edge_ref );
			edge_add(elem->edge[3].id, hash_edge_ref );
			edge_add(elem->edge[4].id, hash_edge_ref );
			edge_add(elem->edge[5].id, hash_edge_ref );
			edge_add(elem->edge[6].id, hash_edge_ref );
			edge_add(elem->edge[7].id, hash_edge_ref );
			edge_add(elem->edge[9].id, hash_edge_ref );
			edge_add(elem->edge[11].id, hash_edge_ref );
		}else if(elem->pad==124){
			edge_add(elem->edge[0].id, hash_edge_ref );
			edge_add(elem->edge[2].id, hash_edge_ref );
			edge_add(elem->edge[4].id, hash_edge_ref );
			edge_add(elem->edge[5].id, hash_edge_ref );
			edge_add(elem->edge[6].id, hash_edge_ref );
			edge_add(elem->edge[7].id, hash_edge_ref );
			edge_add(elem->edge[8].id, hash_edge_ref );
			edge_add(elem->edge[10].id, hash_edge_ref );;
		}else if(elem->pad==123){
			edge_add(elem->edge[0].id, hash_edge_ref );
			edge_add(elem->edge[1].id, hash_edge_ref );
			edge_add(elem->edge[2].id, hash_edge_ref );
			edge_add(elem->edge[3].id, hash_edge_ref );
			edge_add(elem->edge[8].id, hash_edge_ref );
			edge_add(elem->edge[9].id, hash_edge_ref );
			edge_add(elem->edge[10].id, hash_edge_ref );
			edge_add(elem->edge[11].id, hash_edge_ref );
		}
		//template 7
		else if(elem->pad==122){
			edge_add(elem->edge[7].id, hash_edge_ref );
			edge_add(elem->edge[10].id, hash_edge_ref );
			edge_add(elem->edge[11].id, hash_edge_ref );
		}else if(elem->pad==121){
			edge_add(elem->edge[6].id, hash_edge_ref );
			edge_add(elem->edge[9].id, hash_edge_ref );
			edge_add(elem->edge[10].id, hash_edge_ref );
		}else if(elem->pad==120){
			edge_add(elem->edge[5].id, hash_edge_ref );
			edge_add(elem->edge[8].id, hash_edge_ref );
			edge_add(elem->edge[9].id, hash_edge_ref );
		}else if(elem->pad==119){
			edge_add(elem->edge[4].id, hash_edge_ref );
			edge_add(elem->edge[8].id, hash_edge_ref );
			edge_add(elem->edge[11].id, hash_edge_ref );
		}else if(elem->pad==118){
			edge_add(elem->edge[2].id, hash_edge_ref );
			edge_add(elem->edge[3].id, hash_edge_ref );
			edge_add(elem->edge[7].id, hash_edge_ref );
		}else if(elem->pad==117){
			edge_add(elem->edge[1].id, hash_edge_ref );
			edge_add(elem->edge[2].id, hash_edge_ref );
			edge_add(elem->edge[6].id, hash_edge_ref );
		}else if(elem->pad==116){
			edge_add(elem->edge[0].id, hash_edge_ref );
			edge_add(elem->edge[1].id, hash_edge_ref );
			edge_add(elem->edge[5].id, hash_edge_ref );
		}else if(elem->pad==115){
			edge_add(elem->edge[0].id, hash_edge_ref );
			edge_add(elem->edge[3].id, hash_edge_ref );
			edge_add(elem->edge[4].id, hash_edge_ref );
		}
		//template 6
		else if(elem->pad==114){
			edge_add(elem->edge[1].id, hash_edge_ref );
			edge_add(elem->edge[3].id, hash_edge_ref );
			edge_add(elem->edge[6].id, hash_edge_ref );
			edge_add(elem->edge[7].id, hash_edge_ref );
			edge_add(elem->edge[9].id, hash_edge_ref );
			edge_add(elem->edge[11].id, hash_edge_ref );
		}else if(elem->pad==113){
			edge_add(elem->edge[4].id, hash_edge_ref );
			edge_add(elem->edge[5].id, hash_edge_ref );
			edge_add(elem->edge[6].id, hash_edge_ref );
			edge_add(elem->edge[7].id, hash_edge_ref );
			edge_add(elem->edge[9].id, hash_edge_ref );
			edge_add(elem->edge[11].id, hash_edge_ref );
		}else if(elem->pad==112){
			edge_add(elem->edge[1].id, hash_edge_ref );
			edge_add(elem->edge[3].id, hash_edge_ref );
			edge_add(elem->edge[4].id, hash_edge_ref );
			edge_add(elem->edge[5].id, hash_edge_ref );
			edge_add(elem->edge[8].id, hash_edge_ref );
			edge_add(elem->edge[11].id, hash_edge_ref );
		}else if(elem->pad==111){
			edge_add(elem->edge[1].id, hash_edge_ref );
			edge_add(elem->edge[3].id, hash_edge_ref );
			edge_add(elem->edge[4].id, hash_edge_ref );
			edge_add(elem->edge[5].id, hash_edge_ref );
			edge_add(elem->edge[6].id, hash_edge_ref );
			edge_add(elem->edge[7].id, hash_edge_ref );
		}else if(elem->pad==110){
			edge_add(elem->edge[0].id, hash_edge_ref );
			edge_add(elem->edge[2].id, hash_edge_ref );
			edge_add(elem->edge[4].id, hash_edge_ref );
			edge_add(elem->edge[5].id, hash_edge_ref );
			edge_add(elem->edge[6].id, hash_edge_ref );
			edge_add(elem->edge[7].id, hash_edge_ref );
		}else if(elem->pad==109){
			edge_add(elem->edge[0].id, hash_edge_ref );
			edge_add(elem->edge[2].id, hash_edge_ref );
			edge_add(elem->edge[4].id, hash_edge_ref );
			edge_add(elem->edge[7].id, hash_edge_ref );
			edge_add(elem->edge[8].id, hash_edge_ref );
			edge_add(elem->edge[10].id, hash_edge_ref );
		}else if(elem->pad==108){
			edge_add(elem->edge[4].id, hash_edge_ref );
			edge_add(elem->edge[5].id, hash_edge_ref );
			edge_add(elem->edge[6].id, hash_edge_ref );
			edge_add(elem->edge[7].id, hash_edge_ref );
			edge_add(elem->edge[8].id, hash_edge_ref );
			edge_add(elem->edge[10].id, hash_edge_ref );
		}else if(elem->pad==107){
			edge_add(elem->edge[0].id, hash_edge_ref );
			edge_add(elem->edge[2].id, hash_edge_ref );
			edge_add(elem->edge[5].id, hash_edge_ref );
			edge_add(elem->edge[6].id, hash_edge_ref );
			edge_add(elem->edge[8].id, hash_edge_ref );
			edge_add(elem->edge[10].id, hash_edge_ref );
		}else if(elem->pad==106){
			edge_add(elem->edge[0].id, hash_edge_ref );
			edge_add(elem->edge[2].id, hash_edge_ref );
			edge_add(elem->edge[3].id, hash_edge_ref );
			edge_add(elem->edge[8].id, hash_edge_ref );
			edge_add(elem->edge[10].id, hash_edge_ref );
			edge_add(elem->edge[11].id, hash_edge_ref );
		}else if(elem->pad==105){
			edge_add(elem->edge[1].id, hash_edge_ref );
			edge_add(elem->edge[2].id, hash_edge_ref );
			edge_add(elem->edge[3].id, hash_edge_ref );
			edge_add(elem->edge[9].id, hash_edge_ref );
			edge_add(elem->edge[10].id, hash_edge_ref );
			edge_add(elem->edge[11].id, hash_edge_ref );;
		}else if(elem->pad==104){
			edge_add(elem->edge[0].id, hash_edge_ref );
			edge_add(elem->edge[1].id, hash_edge_ref );
			edge_add(elem->edge[2].id, hash_edge_ref );
			edge_add(elem->edge[8].id, hash_edge_ref );
			edge_add(elem->edge[9].id, hash_edge_ref );
			edge_add(elem->edge[10].id, hash_edge_ref );
		}else if(elem->pad==103){
			edge_add(elem->edge[0].id, hash_edge_ref );
			edge_add(elem->edge[1].id, hash_edge_ref );
			edge_add(elem->edge[3].id, hash_edge_ref );
			edge_add(elem->edge[8].id, hash_edge_ref );
			edge_add(elem->edge[9].id, hash_edge_ref );
			edge_add(elem->edge[11].id, hash_edge_ref );
		}
		//template 5
		else if(elem->pad==102){
			edge_add(elem->edge[3].id, hash_edge_ref );
			edge_add(elem->edge[7].id, hash_edge_ref );
		}else if(elem->pad==101){
			edge_add(elem->edge[6].id, hash_edge_ref );
			edge_add(elem->edge[7].id, hash_edge_ref );
		}else if(elem->pad==100){
			edge_add(elem->edge[4].id, hash_edge_ref );
			edge_add(elem->edge[6].id, hash_edge_ref );
		}else if(elem->pad==99){
			edge_add(elem->edge[3].id, hash_edge_ref );
			edge_add(elem->edge[4].id, hash_edge_ref );
		}else if(elem->pad==98){
			edge_add(elem->edge[2].id, hash_edge_ref );
			edge_add(elem->edge[7].id, hash_edge_ref );
		}else if(elem->pad==97){
			edge_add(elem->edge[7].id, hash_edge_ref );
			edge_add(elem->edge[10].id, hash_edge_ref );
		}else if(elem->pad==96){
			edge_add(elem->edge[6].id, hash_edge_ref );
			edge_add(elem->edge[10].id, hash_edge_ref );
		}else if(elem->pad==95){
			edge_add(elem->edge[2].id, hash_edge_ref );
			edge_add(elem->edge[6].id, hash_edge_ref );;
		}else if(elem->pad==94){
			edge_add(elem->edge[1].id, hash_edge_ref );
			edge_add(elem->edge[5].id, hash_edge_ref );
		}else if(elem->pad==93){
			edge_add(elem->edge[1].id, hash_edge_ref );
			edge_add(elem->edge[6].id, hash_edge_ref );
		}else if(elem->pad==92){
			edge_add(elem->edge[6].id, hash_edge_ref );
			edge_add(elem->edge[9].id, hash_edge_ref );
		}else if(elem->pad==91){
			edge_add(elem->edge[5].id, hash_edge_ref );
			edge_add(elem->edge[9].id, hash_edge_ref );
		}else if(elem->pad==90){
			edge_add(elem->edge[0].id, hash_edge_ref );
			edge_add(elem->edge[4].id, hash_edge_ref );
		}else if(elem->pad==89){
			edge_add(elem->edge[0].id, hash_edge_ref );
			edge_add(elem->edge[5].id, hash_edge_ref );
		}else if(elem->pad==88){
			edge_add(elem->edge[5].id, hash_edge_ref );
			edge_add(elem->edge[8].id, hash_edge_ref );;
		}else if(elem->pad==87){
			edge_add(elem->edge[8].id, hash_edge_ref );
			edge_add(elem->edge[4].id, hash_edge_ref );
		}else if(elem->pad==86){
			edge_add(elem->edge[8].id, hash_edge_ref );
			edge_add(elem->edge[11].id, hash_edge_ref );
		}else if(elem->pad==85){
			edge_add(elem->edge[10].id, hash_edge_ref );
			edge_add(elem->edge[11].id, hash_edge_ref );
		}else if(elem->pad==84){
			edge_add(elem->edge[9].id, hash_edge_ref );
			edge_add(elem->edge[10].id, hash_edge_ref );
		}else if(elem->pad==83){
			edge_add(elem->edge[8].id, hash_edge_ref );
			edge_add(elem->edge[9].id, hash_edge_ref );
		}else if(elem->pad==82){
			edge_add(elem->edge[0].id, hash_edge_ref );
			edge_add(elem->edge[3].id, hash_edge_ref );
		}else if(elem->pad==81){
			edge_add(elem->edge[2].id, hash_edge_ref );
			edge_add(elem->edge[3].id, hash_edge_ref );
		}else if(elem->pad==80){
			edge_add(elem->edge[1].id, hash_edge_ref );
			edge_add(elem->edge[2].id, hash_edge_ref );
		}else if(elem->pad==79){
			edge_add(elem->edge[0].id, hash_edge_ref );
			edge_add(elem->edge[1].id, hash_edge_ref );
		}else if(elem->pad==78){
			edge_add(elem->edge[3].id, hash_edge_ref );
			edge_add(elem->edge[7].id, hash_edge_ref );
			edge_add(elem->edge[11].id, hash_edge_ref );
		}else if(elem->pad==77){
			edge_add(elem->edge[4].id, hash_edge_ref );
			edge_add(elem->edge[7].id, hash_edge_ref );
			edge_add(elem->edge[11].id, hash_edge_ref );
		}else if(elem->pad==76){
			edge_add(elem->edge[3].id, hash_edge_ref );
			edge_add(elem->edge[4].id, hash_edge_ref );
			edge_add(elem->edge[11].id, hash_edge_ref );
		}else if(elem->pad==75){
			edge_add(elem->edge[3].id, hash_edge_ref );
			edge_add(elem->edge[4].id, hash_edge_ref );
			edge_add(elem->edge[7].id, hash_edge_ref );
		}else if(elem->pad==74){
			edge_add(elem->edge[2].id, hash_edge_ref );
			edge_add(elem->edge[7].id, hash_edge_ref );
			edge_add(elem->edge[10].id, hash_edge_ref );
		}else if(elem->pad==73){
			edge_add(elem->edge[6].id, hash_edge_ref );
			edge_add(elem->edge[7].id, hash_edge_ref );
			edge_add(elem->edge[10].id, hash_edge_ref );
		}else if(elem->pad==72){
			edge_add(elem->edge[2].id, hash_edge_ref );
			edge_add(elem->edge[6].id, hash_edge_ref );
			edge_add(elem->edge[10].id, hash_edge_ref );
		}else if(elem->pad==71){
			edge_add(elem->edge[2].id, hash_edge_ref );
			edge_add(elem->edge[6].id, hash_edge_ref );
			edge_add(elem->edge[7].id, hash_edge_ref );
		}else if(elem->pad==70){
			edge_add(elem->edge[1].id, hash_edge_ref );
			edge_add(elem->edge[6].id, hash_edge_ref );
			edge_add(elem->edge[9].id, hash_edge_ref );
		}else if(elem->pad==69){
			edge_add(elem->edge[5].id, hash_edge_ref );
			edge_add(elem->edge[6].id, hash_edge_ref );
			edge_add(elem->edge[9].id, hash_edge_ref );
		}else if(elem->pad==68){
			edge_add(elem->edge[1].id, hash_edge_ref );
			edge_add(elem->edge[5].id, hash_edge_ref );
			edge_add(elem->edge[9].id, hash_edge_ref );
		}else if(elem->pad==67){
			edge_add(elem->edge[1].id, hash_edge_ref );
			edge_add(elem->edge[5].id, hash_edge_ref );
			edge_add(elem->edge[6].id, hash_edge_ref );
		}else if(elem->pad==66){
			edge_add(elem->edge[0].id, hash_edge_ref );
			edge_add(elem->edge[4].id, hash_edge_ref );
			edge_add(elem->edge[8].id, hash_edge_ref );
		}else if(elem->pad==65){
			edge_add(elem->edge[4].id, hash_edge_ref );
			edge_add(elem->edge[5].id, hash_edge_ref );
			edge_add(elem->edge[8].id, hash_edge_ref );
		}else if(elem->pad==64){
			edge_add(elem->edge[0].id, hash_edge_ref );
			edge_add(elem->edge[5].id, hash_edge_ref );
			edge_add(elem->edge[8].id, hash_edge_ref );
		}else if(elem->pad==63){
			edge_add(elem->edge[0].id, hash_edge_ref );
			edge_add(elem->edge[4].id, hash_edge_ref );
			edge_add(elem->edge[5].id, hash_edge_ref );
		}else if(elem->pad==62){
			edge_add(elem->edge[8].id, hash_edge_ref );
			edge_add(elem->edge[10].id, hash_edge_ref );
			edge_add(elem->edge[11].id, hash_edge_ref );
		}else if(elem->pad==61){
			edge_add(elem->edge[9].id, hash_edge_ref );
			edge_add(elem->edge[10].id, hash_edge_ref );
			edge_add(elem->edge[11].id, hash_edge_ref );
		}else if(elem->pad==60){
			edge_add(elem->edge[8].id, hash_edge_ref );
			edge_add(elem->edge[9].id, hash_edge_ref );
			edge_add(elem->edge[11].id, hash_edge_ref );
		}else if(elem->pad==59){
			edge_add(elem->edge[8].id, hash_edge_ref );
			edge_add(elem->edge[9].id, hash_edge_ref );
			edge_add(elem->edge[10].id, hash_edge_ref );
		}else if(elem->pad==58){
			edge_add(elem->edge[0].id, hash_edge_ref );
			edge_add(elem->edge[2].id, hash_edge_ref );
			edge_add(elem->edge[11].id, hash_edge_ref );
		}else if(elem->pad==57){
			edge_add(elem->edge[1].id, hash_edge_ref );
			edge_add(elem->edge[2].id, hash_edge_ref );
			edge_add(elem->edge[3].id, hash_edge_ref );
		}else if(elem->pad==56){
			edge_add(elem->edge[0].id, hash_edge_ref );
			edge_add(elem->edge[1].id, hash_edge_ref );
			edge_add(elem->edge[2].id, hash_edge_ref );
		}else if(elem->pad==55){
			edge_add(elem->edge[0].id, hash_edge_ref );
			edge_add(elem->edge[1].id, hash_edge_ref );
			edge_add(elem->edge[3].id, hash_edge_ref );
		}else if(elem->pad==54){
			edge_add(elem->edge[3].id, hash_edge_ref );
			edge_add(elem->edge[4].id, hash_edge_ref );
			edge_add(elem->edge[7].id, hash_edge_ref );
			edge_add(elem->edge[11].id, hash_edge_ref );
		}else if(elem->pad==53){
			edge_add(elem->edge[2].id, hash_edge_ref );
			edge_add(elem->edge[6].id, hash_edge_ref );
			edge_add(elem->edge[7].id, hash_edge_ref );
			edge_add(elem->edge[10].id, hash_edge_ref );
		}else if(elem->pad==52){
			edge_add(elem->edge[1].id, hash_edge_ref );
			edge_add(elem->edge[5].id, hash_edge_ref );
			edge_add(elem->edge[6].id, hash_edge_ref );
			edge_add(elem->edge[9].id, hash_edge_ref );
		}else if(elem->pad==51){
			edge_add(elem->edge[0].id, hash_edge_ref );
			edge_add(elem->edge[4].id, hash_edge_ref );
			edge_add(elem->edge[5].id, hash_edge_ref );
			edge_add(elem->edge[8].id, hash_edge_ref );
		}else if(elem->pad==50){
			edge_add(elem->edge[8].id, hash_edge_ref );
			edge_add(elem->edge[9].id, hash_edge_ref );
			edge_add(elem->edge[10].id, hash_edge_ref );
			edge_add(elem->edge[11].id, hash_edge_ref );
		}else if(elem->pad==49){
			edge_add(elem->edge[0].id, hash_edge_ref );
			edge_add(elem->edge[1].id, hash_edge_ref );
			edge_add(elem->edge[2].id, hash_edge_ref );
			edge_add(elem->edge[3].id, hash_edge_ref );
		}
		//template 4
		else if(elem->pad==48){
			edge_add(elem->edge[4].id, hash_edge_ref );
			edge_add(elem->edge[5].id, hash_edge_ref );
			edge_add(elem->edge[9].id, hash_edge_ref );
			edge_add(elem->edge[11].id, hash_edge_ref );
		}else if(elem->pad==47){
			edge_add(elem->edge[6].id, hash_edge_ref );
			edge_add(elem->edge[7].id, hash_edge_ref );
			edge_add(elem->edge[9].id, hash_edge_ref );
			edge_add(elem->edge[11].id, hash_edge_ref );
		}else if(elem->pad==46){
			edge_add(elem->edge[1].id, hash_edge_ref );
			edge_add(elem->edge[3].id, hash_edge_ref );
			edge_add(elem->edge[6].id, hash_edge_ref );
			edge_add(elem->edge[7].id, hash_edge_ref );
		}else if(elem->pad==45){
			edge_add(elem->edge[1].id, hash_edge_ref );
			edge_add(elem->edge[3].id, hash_edge_ref );
			edge_add(elem->edge[4].id, hash_edge_ref );
			edge_add(elem->edge[5].id, hash_edge_ref );
		}else if(elem->pad==44){
			edge_add(elem->edge[4].id, hash_edge_ref );
			edge_add(elem->edge[7].id, hash_edge_ref );
			edge_add(elem->edge[8].id, hash_edge_ref );
			edge_add(elem->edge[10].id, hash_edge_ref );
		}else if(elem->pad==43){
			edge_add(elem->edge[5].id, hash_edge_ref );
			edge_add(elem->edge[6].id, hash_edge_ref );
			edge_add(elem->edge[8].id, hash_edge_ref );
			edge_add(elem->edge[10].id, hash_edge_ref );
		}else if(elem->pad==42){
			edge_add(elem->edge[0].id, hash_edge_ref );
			edge_add(elem->edge[2].id, hash_edge_ref );
			edge_add(elem->edge[5].id, hash_edge_ref );
			edge_add(elem->edge[6].id, hash_edge_ref );
		}else if(elem->pad==41){
			edge_add(elem->edge[0].id, hash_edge_ref );
			edge_add(elem->edge[2].id, hash_edge_ref );
			edge_add(elem->edge[4].id, hash_edge_ref );
			edge_add(elem->edge[7].id, hash_edge_ref );
		}else if(elem->pad==40){
			edge_add(elem->edge[0].id, hash_edge_ref );
			edge_add(elem->edge[3].id, hash_edge_ref );
			edge_add(elem->edge[8].id, hash_edge_ref );
			edge_add(elem->edge[11].id, hash_edge_ref );
		}else if(elem->pad==39){
			edge_add(elem->edge[2].id, hash_edge_ref );
			edge_add(elem->edge[3].id, hash_edge_ref );
			edge_add(elem->edge[10].id, hash_edge_ref );
			edge_add(elem->edge[11].id, hash_edge_ref );
		}else if(elem->pad==38){
			edge_add(elem->edge[1].id, hash_edge_ref );
			edge_add(elem->edge[2].id, hash_edge_ref );
			edge_add(elem->edge[9].id, hash_edge_ref );
			edge_add(elem->edge[10].id, hash_edge_ref );
		}else if(elem->pad==37){
			edge_add(elem->edge[0].id, hash_edge_ref );
			edge_add(elem->edge[1].id, hash_edge_ref );
			edge_add(elem->edge[8].id, hash_edge_ref );
			edge_add(elem->edge[9].id, hash_edge_ref );
		}
		//template 3
		else if(elem->pad==36){
			edge_add(elem->edge[3].id, hash_edge_ref );
			edge_add(elem->edge[11].id, hash_edge_ref );
		}else if(elem->pad==35){
			edge_add(elem->edge[1].id, hash_edge_ref );
			edge_add(elem->edge[10].id, hash_edge_ref );
		}else if(elem->pad==34){
			edge_add(elem->edge[1].id, hash_edge_ref );
			edge_add(elem->edge[9].id, hash_edge_ref );
		}else if(elem->pad==33){
			edge_add(elem->edge[0].id, hash_edge_ref );
			edge_add(elem->edge[8].id, hash_edge_ref );
		}else if(elem->pad==32){
			edge_add(elem->edge[4].id, hash_edge_ref );
			edge_add(elem->edge[7].id, hash_edge_ref );
		}else if(elem->pad==31){
			edge_add(elem->edge[6].id, hash_edge_ref );
			edge_add(elem->edge[7].id, hash_edge_ref );
		}else if(elem->pad==30){
			edge_add(elem->edge[5].id, hash_edge_ref );
			edge_add(elem->edge[6].id, hash_edge_ref );
		}else if(elem->pad==29){
			edge_add(elem->edge[4].id, hash_edge_ref );
			edge_add(elem->edge[5].id, hash_edge_ref );
		}else if(elem->pad==28){
			edge_add(elem->edge[9].id, hash_edge_ref );
			edge_add(elem->edge[11].id, hash_edge_ref );
		}else if(elem->pad==27){
			edge_add(elem->edge[8].id, hash_edge_ref );
			edge_add(elem->edge[10].id, hash_edge_ref );
		}else if(elem->pad==26){
			edge_add(elem->edge[1].id, hash_edge_ref );
			edge_add(elem->edge[3].id, hash_edge_ref );
		}else if(elem->pad==25){
			edge_add(elem->edge[0].id, hash_edge_ref );
			edge_add(elem->edge[2].id, hash_edge_ref );
		}
		//template 2
		else if(elem->pad==24){
			edge_add(elem->edge[4].id, hash_edge_ref );
			edge_add(elem->edge[5].id, hash_edge_ref );
			edge_add(elem->edge[6].id, hash_edge_ref );
			edge_add(elem->edge[7].id, hash_edge_ref );
		}else if(elem->pad==23){
			edge_add(elem->edge[1].id, hash_edge_ref );
			edge_add(elem->edge[3].id, hash_edge_ref );
			edge_add(elem->edge[9].id, hash_edge_ref );
			edge_add(elem->edge[11].id, hash_edge_ref );
		}else if(elem->pad==22){
			edge_add(elem->edge[0].id, hash_edge_ref );
			edge_add(elem->edge[2].id, hash_edge_ref );
			edge_add(elem->edge[8].id, hash_edge_ref );
			edge_add(elem->edge[10].id, hash_edge_ref );
		}
		//template 1
		else if(elem->pad==21){
			edge_add(elem->edge[11].id, hash_edge_ref );
		}else if(elem->pad==20){
			edge_add(elem->edge[10].id, hash_edge_ref );
		}else if(elem->pad==19){
			edge_add(elem->edge[9].id, hash_edge_ref );
		}else if(elem->pad==18){
			edge_add(elem->edge[8].id, hash_edge_ref );
		}else if(elem->pad==17){
			edge_add(elem->edge[7].id, hash_edge_ref );
		}else if(elem->pad==16){
			edge_add(elem->edge[6].id, hash_edge_ref );
		}else if(elem->pad==15){
			edge_add(elem->edge[5].id, hash_edge_ref );
		}else if(elem->pad==14){
			edge_add(elem->edge[4].id, hash_edge_ref );
		}else if(elem->pad==13){
			edge_add(elem->edge[3].id, hash_edge_ref );
		}else if(elem->pad==12){
			edge_add(elem->edge[2].id, hash_edge_ref );
		}else if(elem->pad==11){
			edge_add(elem->edge[1].id, hash_edge_ref );
		}else if(elem->pad==10){
			edge_add(elem->edge[0].id, hash_edge_ref );
		}
	}


	for (int iedge = 0; iedge< mesh->shared_edges.elem_count; iedge++)
	{
		size_t position;
		octant_edge_t* shared_ed = (octant_edge_t*) sc_array_index(&mesh->shared_edges, iedge);
		octant_edge_t key;
		key.id = shared_ed->id;
		bool out =  sc_hash_array_lookup(hash_edge_ref, &key, &position);
		shared_ed->ref = out;
	}

#ifdef HEXA_DEBUG_
	if(false)
	{
		fprintf(mesh->fdbg ,"shared edges status\n");
		fprintf(mesh->fdbg ,"shared edges number:%d\n",mesh->shared_edges.elem_count);
		for (int iedge = 0; iedge< mesh->shared_edges.elem_count; iedge++)
		{
			octant_edge_t* shared_ed = (octant_edge_t*) sc_array_index(&mesh->shared_edges, iedge);
			fprintf(mesh->fdbg ,"id: %d  ref:%d  rank:%d\n",shared_ed->id,shared_ed->ref,mesh->mpi_rank);
		}
	}
#endif

}

void Edge_comunicationNew(hexa_tree_t* mesh, sc_hash_array_t* hash_edge_ref, sc_hash_array_t* hash_edge_ref_old)
{

	bool clamped=true;
	sc_hash_array_t* SendTo   = (sc_hash_array_t *) sc_hash_array_new(sizeof(message_t), processors_hash_fn, processors_equal_fn, &clamped);
	sc_hash_array_t* RecvFrom   = (sc_hash_array_t *) sc_hash_array_new(sizeof(message_t), processors_hash_fn, processors_equal_fn, &clamped);

	//rank [n+1] to [n]
	//prepare the buffer
	for(int iedge = 0; iedge < mesh->shared_edges.elem_count; ++iedge)
	{
		shared_edge_t* se = (shared_edge_t*) sc_array_index(&mesh->shared_edges,iedge);
		octant_edge_t key;
		key.id = se->id;
		size_t position;
		bool lse = sc_hash_array_lookup(hash_edge_ref, &key, &position);
		bool lseold = sc_hash_array_lookup(hash_edge_ref_old, &key, &position);

		if(lse)
		{
			for(int j = 0; j < se->listSz; j++)
			{
				if(se->rankList[j] < mesh->mpi_rank)
				{
					message_t* m = (message_t*)sc_hash_array_insert_unique(SendTo,&se->rankList[j],&position);
					if(m!=NULL)
					{
						m->rank  = se->rankList[j];
						sc_array_init(&m->idxs, sizeof(uint64_t));
						uint64_t* p = (uint64_t*) sc_array_push(&m->idxs);
 *p = se->id;
					}
					else
					{
						message_t* m = (message_t*)sc_array_index(&SendTo->a, position);
						uint64_t* p = (uint64_t*) sc_array_push(&m->idxs);
 *p = se->id;
					}
				}
				else if(se->rankList[j] > mesh->mpi_rank)
				{
					message_t* m = (message_t*)sc_hash_array_insert_unique(RecvFrom,&se->rankList[j],&position);
					if(m!=NULL)
					{
						m->rank  = se->rankList[j];
					}
				}
			}
		}
	}

	int max_SendBuffer = 0;
	for(int i = 0; i < SendTo->a.elem_count; i++)
	{
		message_t* m = (message_t*) sc_array_index(&SendTo->a, i);
		max_SendBuffer+= m->idxs.elem_count;
	}
	long long* sendbuf    = (long long*)malloc(max_SendBuffer*sizeof(long long));
	int offset = 0;

	// rank[n+1] send to rank[n]
	for(int i = 0; i < SendTo->a.elem_count; i++)
	{
		//printf("Estou montando o buffer do proc %d\n",mesh->mpi_rank);
		message_t* m = (message_t*) sc_array_index(&SendTo->a, i);
		for(int j = 0; j< m->idxs.elem_count;j++)
		{
			long long *mm = (long long*) sc_array_index(&m->idxs, j);
			sendbuf[offset+j] = (long long) *mm;
		}
		//printf("Estou enviando do proc %d\n",mesh->mpi_rank);
		MPI_Send(&sendbuf[offset], m->idxs.elem_count, MPI_LONG_LONG, m->rank,0,MPI_COMM_WORLD);
		offset += m->idxs.elem_count;
	}

	offset = 0;
	for(int i = 0; i < RecvFrom->a.elem_count; i++){
		message_t* m = (message_t*) sc_array_index(&RecvFrom->a, i);

		//printf("Fazendo a probe do proc %d\n",m->rank);
		MPI_Status status;
		// Probe for an incoming message from process mesh->mpi_rank
		MPI_Probe(m->rank, 0, MPI_COMM_WORLD, &status);
		//printf("Fiz a probe do proc %d\n",m->rank);

		//printf("Descobrindo quantos elementos vem proc %d para o proc %d\n",mesh->mpi_rank,m->rank);
		int number_amount;
		// When probe returns, the status object has the size and other
		// attributes of the incoming message. Get the message size
		MPI_Get_count(&status, MPI_INT, &number_amount);

		// Allocate a buffer to hold the incoming numbers
		long long* number_buf = (long long*)malloc(sizeof(long long) * number_amount);
		//printf("Recebendo os elementos do proc %d para o proc %d\n",mesh->mpi_rank,m->rank);
		// Now receive the message with the allocated buffer
		MPI_Recv(number_buf, number_amount, MPI_LONG_LONG, m->rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		//printf("I dynamically received %d numbers from %d.\n", number_amount,m->rank);
		for(int j = 0; j < number_amount; ++j)
		{
			edge_add(number_buf[j], hash_edge_ref );
		}
		free(number_buf);
	}
	//clean
	sc_hash_array_truncate (SendTo);
	sc_hash_array_truncate (RecvFrom);
	free(sendbuf);

	//rank [n] to [n+1]
	//prepare the buffer
	for(int iedge = 0; iedge < mesh->shared_edges.elem_count; ++iedge)
	{
		shared_edge_t* se = (shared_edge_t*) sc_array_index(&mesh->shared_edges,iedge);
		octant_edge_t key;
		key.id = se->id;
		size_t position;
		bool lse = sc_hash_array_lookup(hash_edge_ref, &key, &position);
		bool lseold = sc_hash_array_lookup(hash_edge_ref_old, &key, &position);

		if(lse)
		{
			for(int j = 0; j < se->listSz; j++)
			{
				if(se->rankList[j] > mesh->mpi_rank)
				{
					message_t* m = (message_t*)sc_hash_array_insert_unique(SendTo,&se->rankList[j],&position);
					if(m!=NULL)
					{
						m->rank  = se->rankList[j];
						sc_array_init(&m->idxs, sizeof(uint64_t));
						uint64_t* p = (uint64_t*) sc_array_push(&m->idxs);
 *p = se->id;
					}
					else
					{
						message_t* m = (message_t*)sc_array_index(&SendTo->a, position);
						uint64_t* p = (uint64_t*) sc_array_push(&m->idxs);
 *p = se->id;
					}
				}
				else if(se->rankList[j] < mesh->mpi_rank)
				{
					message_t* m = (message_t*)sc_hash_array_insert_unique(RecvFrom,&se->rankList[j],&position);
					if(m!=NULL)
					{
						m->rank  = se->rankList[j];
					}
				}
			}
		}
	}

	max_SendBuffer = 0;
	for(int i = 0; i < SendTo->a.elem_count; i++)
	{
		message_t* m = (message_t*) sc_array_index(&SendTo->a, i);
		max_SendBuffer+= m->idxs.elem_count;
	}
	sendbuf    = (long long*)malloc(max_SendBuffer*sizeof(long long));
	offset = 0;

	// rank[n+1] send to rank[n]
	for(int i = 0; i < SendTo->a.elem_count; i++)
	{
		//printf("Estou montando o buffer do proc %d\n",mesh->mpi_rank);
		message_t* m = (message_t*) sc_array_index(&SendTo->a, i);
		for(int j = 0; j< m->idxs.elem_count;j++)
		{
			long long *mm = (long long*) sc_array_index(&m->idxs, j);
			sendbuf[offset+j] = (long long) *mm;
		}
		//printf("Estou enviando do proc %d\n",mesh->mpi_rank);
		MPI_Send(&sendbuf[offset], m->idxs.elem_count, MPI_LONG_LONG, m->rank,0,MPI_COMM_WORLD);
		offset += m->idxs.elem_count;
	}

	offset = 0;
	for(int i = 0; i < RecvFrom->a.elem_count; i++){
		message_t* m = (message_t*) sc_array_index(&RecvFrom->a, i);

		//printf("Fazendo a probe do proc %d\n",m->rank);
		MPI_Status status;
		// Probe for an incoming message from process mesh->mpi_rank
		MPI_Probe(m->rank, 0, MPI_COMM_WORLD, &status);
		//printf("Fiz a probe do proc %d\n",m->rank);

		//printf("Descobrindo quantos elementos vem proc %d para o proc %d\n",mesh->mpi_rank,m->rank);
		int number_amount;
		// When probe returns, the status object has the size and other
		// attributes of the incoming message. Get the message size
		MPI_Get_count(&status, MPI_INT, &number_amount);

		// Allocate a buffer to hold the incoming numbers
		long long* number_buf = (long long*)malloc(sizeof(long long) * number_amount);
		//printf("Recebendo os elementos do proc %d para o proc %d\n",mesh->mpi_rank,m->rank);
		// Now receive the message with the allocated buffer
		MPI_Recv(number_buf, number_amount, MPI_LONG_LONG, m->rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		//printf("1 dynamically received %d numbers from 0.\n", number_amount);
		for(int j = 0; j < number_amount; ++j)
		{
			edge_add(number_buf[j], hash_edge_ref );
		}
		free(number_buf);
	}
	//clean
	sc_hash_array_truncate (SendTo);
	sc_hash_array_truncate (RecvFrom);
	free(sendbuf);

}

void Edge_comunication(hexa_tree_t* mesh, sc_hash_array_t* hash_edge_ref)
{

	size_t position;
	bool clamped=true;
	MPI_Request *requests;
	MPI_Status  *statuses;
	long long   *recvbuf;
	long long   *sendbuf;
	int n_requests;

	n_requests = mesh->comm_map_edge.nrequests;
	recvbuf    = (long long*)malloc(2*mesh->comm_map_edge.max_recvbuf_size*sizeof(long long));
	sendbuf    = (long long*)malloc(2*mesh->comm_map_edge.max_sendbuf_size*sizeof(long long));

	requests = (MPI_Request*) malloc (n_requests*sizeof(MPI_Request));
	statuses = (MPI_Status*)  malloc (n_requests*sizeof(MPI_Status));

	// rank[n+1] send to rank[n]
	int c = 0;

	int max_recvbuf_size = 2*mesh->comm_map_edge.max_recvbuf_size;
	int max_sendbuf_size = 2*mesh->comm_map_edge.max_sendbuf_size;

	printf("max_recvbuf_size:%d max_sendbuf_size:%d n_requests:%d proc:%d\n",max_recvbuf_size,max_sendbuf_size,n_requests,mesh->mpi_rank);
	int offset = 0;

	// post all non-blocking receives
	for(int i = 0; i < mesh->comm_map_edge.RecvFrom.elem_count; ++i)
	{
		message_t *m = (message_t*) sc_array_index(&mesh->comm_map_edge.RecvFrom, i);
		MPI_Irecv(&recvbuf[offset], 2*m->idxs.elem_count, MPI_LONG_LONG, m->rank,0,MPI_COMM_WORLD, &requests[c]);
		offset += 2*m->idxs.elem_count;
		c++;
	}
	assert(offset == max_recvbuf_size);


	offset = 0;
	for(int i=0; i < mesh->comm_map_edge.SendTo.elem_count; i++)
	{
		message_t* m = (message_t*) sc_array_index(&mesh->comm_map_edge.SendTo, i);
		for(int j = 0; j< m->idxs.elem_count;j++)
		{
			long long* mm = (long long*) sc_array_index(&m->idxs, j);

			octant_edge_t key;
			key.id = *mm;
			long long out = sc_hash_array_lookup(hash_edge_ref, &key, &position);
			sendbuf[2*j+offset] = (long long) *mm;
			sendbuf[2*j+offset+1] = (long long) out;
		}
		MPI_Isend(&sendbuf[offset], 2*m->idxs.elem_count, MPI_LONG_LONG, m->rank,0,MPI_COMM_WORLD, &requests[c]);
		offset += 2*m->idxs.elem_count;
		c++;
	}

	assert(offset == max_sendbuf_size);

	assert(c == n_requests);

	MPI_Waitall(n_requests,requests,statuses);


	offset =0;
	//spread the information
	for(int i = 0; i < mesh->comm_map_edge.RecvFrom.elem_count; ++i)
	{
		message_t *m = (message_t*) sc_array_index(&mesh->comm_map_edge.RecvFrom, i);
		for(int j = 0; j < m->idxs.elem_count; ++j)
		{
			if(recvbuf[2*j+offset+1]){
				edge_add(recvbuf[2*j+offset], hash_edge_ref );
			}
		}
		offset += 2*m->idxs.elem_count;
	}

	free(recvbuf);
	free(sendbuf);
	free(requests);
	free(statuses);

	// rank[n] send to rank[n+1]
	recvbuf    = (long long*)malloc(2*mesh->comm_map_edge.max_sendbuf_size*sizeof(long long));
	sendbuf    = (long long*)malloc(2*mesh->comm_map_edge.max_recvbuf_size*sizeof(long long));

	requests = (MPI_Request*) malloc (n_requests*sizeof(MPI_Request));
	statuses = (MPI_Status*)  malloc (n_requests*sizeof(MPI_Status));

	max_recvbuf_size = 2*mesh->comm_map_edge.max_sendbuf_size;
	max_sendbuf_size = 2*mesh->comm_map_edge.max_recvbuf_size;

	c = 0;
	offset = 0;

	// post all non-blocking receives
	for(int i = 0; i < mesh->comm_map_edge.SendTo.elem_count; ++i)
	{
		message_t *m = (message_t*) sc_array_index(&mesh->comm_map_edge.SendTo, i);
		MPI_Irecv(&recvbuf[offset], 2*m->idxs.elem_count, MPI_LONG_LONG, m->rank,0,MPI_COMM_WORLD, &requests[c]);
		offset += 2*m->idxs.elem_count;
		c++;
	}
	assert(offset == max_recvbuf_size);

	offset = 0;
	for(int i=0; i < mesh->comm_map_edge.RecvFrom.elem_count; i++)
	{
		message_t* m = (message_t*) sc_array_index(&mesh->comm_map_edge.RecvFrom, i);
		for(int j = 0; j< m->idxs.elem_count;j++)
		{
			long long* mm = (long long*) sc_array_index(&m->idxs, j);
			octant_edge_t key;
			key.id = *mm;
			int out =  sc_hash_array_lookup(hash_edge_ref, &key, &position);
			sendbuf[2*j+offset] = (long long) *mm;
			sendbuf[2*j+offset+1] = (long long)   out;
		}
		MPI_Isend(&sendbuf[offset], 2*m->idxs.elem_count, MPI_LONG_LONG, m->rank,0,MPI_COMM_WORLD, &requests[c]);
		offset += 2*m->idxs.elem_count;
		c++;
	}

	assert(offset == max_sendbuf_size);

	assert(c == n_requests);

	MPI_Waitall(n_requests,requests,statuses);

	offset =0;
	//add the information in the hash_edge_ref
	for(int i = 0; i < mesh->comm_map_edge.SendTo.elem_count; ++i)
	{
		message_t *m = (message_t*) sc_array_index(&mesh->comm_map_edge.SendTo, i);
		for(int j = 0; j < m->idxs.elem_count; ++j){
			if(recvbuf[2*j+offset+1]){
				edge_add(recvbuf[2*j+offset], hash_edge_ref );
			}
		}
		offset += 2*m->idxs.elem_count;
	}

	free(&recvbuf[0]);
	free(&sendbuf[0]);
	free(requests);
	free(statuses);
}

void Edge_propagation(hexa_tree_t* mesh, std::vector<int>& elements_ids, sc_hash_array_t* hash_edge_ref, sc_hash_array_t* hash_edge, sc_hash_array_t* hash_elem)
{

	for(int iedgeref = 0; iedgeref < hash_edge_ref->a.elem_count; iedgeref++)
	{
		octant_edge_t *edgeref = (octant_edge_t*) sc_array_index(&hash_edge_ref->a, iedgeref);
		size_t position;
		edge_t key;
		key.id = edgeref->id;
		//edgeref->ref = true;
		bool out =  sc_hash_array_lookup(hash_edge, &key, &position);
		if(out)
		{
			edge_t *edge = (edge_t*) sc_array_index(&hash_edge->a, position);
			edge->ref = true;
			for(int iel = 0; iel< edge->list_elem; iel++)
			{
				octant_t *elem = (octant_t*) sc_array_index(&mesh->elements, edge->elem[iel]);

				for(int iedge = 0; iedge < 12; iedge++)
				{
					if(elem->edge[iedge].id == edgeref->id)
					{
						elem->edge[iedge].ref = true;
					}
				}

				elem->pad = -1;
				octant_t *r;
				key.id = elem->id;
				r = (octant_t*) sc_hash_array_insert_unique(hash_elem, &key, &position);
				if(r!=NULL)
				{
					r->id = elem->id;
				}
			}
		}
	}

	//elements_ids.clear();
	for(int iel= elements_ids.size(); iel < hash_elem->a.elem_count; iel++ ){
		octant_t *elem = (octant_t*) sc_array_index(&hash_elem->a, iel);
		elements_ids.push_back(elem->id);
	}



	size_t position;
	octant_t key;

	for(int iel= 0; iel < mesh->elements.elem_count; iel++ )
	{
		octant_t *elem = (octant_t*) sc_array_index(&mesh->elements, iel);

		for (int iedge = 0; iedge < 12; ++iedge)
		{
			bool out = false;
			octant_edge_t key;
			key.id = elem->edge[iedge].id;
			out =  sc_hash_array_lookup(hash_edge_ref, &key, &position);
			if(out)
			{
				octant_t *r;
				key.id = elem->id;
				r = (octant_t*) sc_hash_array_insert_unique(hash_elem, &key, &position);
				if(r!=NULL)
				{
					r->id = elem->id;
				}
				elem->edge[iedge].ref=true;
				elem->pad = -1;
				//elements_ids.push_back(iel);
			}
		}
	}

	//elements_ids.clear();
	for(int iel= elements_ids.size(); iel < hash_elem->a.elem_count; iel++ )
	{
		octant_t *elem = (octant_t*) sc_array_index(&hash_elem->a, iel);
		elements_ids.push_back(elem->id);
	}

}

void CheckOctreeTemplate(hexa_tree_t* mesh, sc_hash_array_t* hash_edge_ref)
{

	bool clamped = true;
	//criando a hash de elementos afetados
	sc_hash_array_t*   hash_elem  = (sc_hash_array_t *)sc_hash_array_new(sizeof(octant_t), el_hash_id, el_equal_id, &clamped);
	size_t position;
	octant_t key;
	std::vector<int> elements_ids;

	for(int ioc = 0; ioc < mesh->oct.elem_count; ioc++)
	{
		octree_t * oc = (octree_t*) sc_array_index (&mesh->oct, ioc);
		for(int iel = 0; iel < 8; iel++)
		{
			elements_ids.push_back(oc->id[iel]);
			octant_t *r;
			key.id = oc->id[iel];
			r = (octant_t*) sc_hash_array_insert_unique(hash_elem, &key, &position);
			if(r!=NULL)
			{
				r->id = elements_ids[iel];
			}
			else
			{
				printf("Verificar o elemento numero %d\n",elements_ids[iel]);
				printf("CheckOctreeTemplate in Pillow Interface\n");
			}
		}
	}

	auto start = std::chrono::steady_clock::now( );
	//TODO trocar o hashid e fn...
	sc_hash_array_t* hash_edge  = (sc_hash_array_t *)sc_hash_array_new(sizeof(edge_t), edget_id_hash, edget_id_equal, &clamped);
	/////////////////
	// create the edge structure
	for (int iel = 0; iel < mesh->elements.elem_count; ++iel)
	{

		octant_t *elem = (octant_t*) sc_array_index(&mesh->elements, iel);

		for (int iedge = 0; iedge < 12; iedge++)
		{
			edge_t key;
			key.id = elem->edge[iedge].id;
			edge_t* edge = (edge_t*) sc_hash_array_insert_unique (hash_edge, &key, &position);
			if(edge != NULL)
			{
				edge->id = key.id;
				edge->list_elem = 1;
				edge->elem[edge->list_elem-1] = elem->id;
				edge->ref = elem->edge[iedge].ref;
			}
			else
			{
				edge = (edge_t*) sc_array_index(&hash_edge->a, position);
				edge->elem[edge->list_elem] = elem->id;
				edge->list_elem++;
			}
		}
	}

	auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>( std::chrono::steady_clock::now( ) - start );
	fprintf(mesh->profile,"        Time in create the edge structure %lld millisecond(s).\n",elapsed.count());

	int iter_count = 0;
	int diff = 50;
	int global_diff = 10;
	while(iter_count < 100 && global_diff != 0)
	{

		sc_hash_array_t* hash_edge_ref_old  = (sc_hash_array_t *) hash_edge_ref;
		//printf("Numero de elementos:%d numero de arestas:%d\n",elements_ids.size(),hash_edge_ref->a.elem_count);
		fprintf(mesh->profile,"        Number of Elements:%d Number of edges:%d\n",elements_ids.size(),hash_edge_ref->a.elem_count);
		start = std::chrono::steady_clock::now( );
		int edgecount = hash_edge_ref->a.elem_count;
		IdentifyTemplate(mesh, elements_ids, hash_edge_ref);
		elapsed = std::chrono::duration_cast<std::chrono::milliseconds>( std::chrono::steady_clock::now( ) - start );
		fprintf(mesh->profile,"        Time in IdentifyTemplate %lld millisecond(s).\n",elapsed.count());
		if(mesh->mpi_rank==0) std::cout << "Time in IdentifyTemplate "<< elapsed.count() <<" millisecond(s)."<< std::endl;

		start = std::chrono::steady_clock::now( );
		Edge_identification( mesh, elements_ids, hash_edge_ref);
		elapsed = std::chrono::duration_cast<std::chrono::milliseconds>( std::chrono::steady_clock::now( ) - start );
		fprintf(mesh->profile,"        Time in Edge_identification %lld millisecond(s).\n",elapsed.count());
		if(mesh->mpi_rank==0) std::cout << "Time in Edge_identification "<< elapsed.count() <<" millisecond(s)."<< std::endl;

		start = std::chrono::steady_clock::now( );
		Edge_propagation (mesh, elements_ids, hash_edge_ref,hash_edge,hash_elem);
		elapsed = std::chrono::duration_cast<std::chrono::milliseconds>( std::chrono::steady_clock::now( ) - start );
		fprintf(mesh->profile,"        Time in Edge_propagation %lld millisecond(s).\n",elapsed.count());
		if(mesh->mpi_rank==0)  std::cout << "Time in Edge_propagation "<< elapsed.count() <<" millisecond(s)."<< std::endl;

		start = std::chrono::steady_clock::now( );
		//TODO add a copy of hash_edge and send only the new
		//edges added to the hash...
		Edge_comunicationNew(mesh, hash_edge_ref,hash_edge_ref_old);
		elapsed = std::chrono::duration_cast<std::chrono::milliseconds>( std::chrono::steady_clock::now( ) - start );
		fprintf(mesh->profile,"        Time in Edge_comunication %lld millisecond(s).\n",elapsed.count());
		if(mesh->mpi_rank==0)  std::cout << "Time in Edge_comunication "<< elapsed.count() <<" millisecond(s)."<< std::endl;
		diff = hash_edge_ref->a.elem_count - edgecount;

		int loc_hash = hash_edge_ref->a.elem_count;
		int glob_hash = 0;
		MPI_Allreduce(&diff, &global_diff, 1, MPI_INT, MPI_SUM,MPI_COMM_WORLD);
		MPI_Allreduce(&loc_hash, &glob_hash, 1, MPI_INT, MPI_SUM,MPI_COMM_WORLD);

		if((iter_count % 1) == 0 && mesh->mpi_rank == 0)
		{
			printf("         Iteration number:%d; diff equals to:%d; number of edges in the hash:%d\n",iter_count,global_diff,glob_hash);
		}

		if(global_diff == 0 && mesh->mpi_rank == 0)
		{
			printf("         The code used %d iteractions to propagate the edge contamination\n",iter_count);
		}
		iter_count++;
	}

	IdentifyTemplate(mesh, elements_ids, hash_edge_ref);

	//sc_hash_array_rip(hash_edge_ref,&mesh->edges_ref);
}
 */

/*
vector<int> RotateHex(int* rot, int* sym)
																																																																																																																																																																																																																																																																																																																																																																																																																																																																																				{

	vector<int> order;
	int aux[8];

	for(int k = 0;k<8;k++) order.push_back(k);

	if(rot[0]==-1)
	{

		for(int i=0;i<8;i++) aux[i]=order[i];

		//order = {3,2,6,7,0,1,5,4};
		order.clear();

		order.push_back(aux[3]);
		order.push_back(aux[2]);
		order.push_back(aux[6]);
		order.push_back(aux[7]);

		order.push_back(aux[0]);
		order.push_back(aux[1]);
		order.push_back(aux[5]);
		order.push_back(aux[4]);

	}

	if(rot[0]==1)
	{

		for(int i=0;i<8;i++) aux[i]=order[i];

		//order = {4,0,7,3,5,1,2,6};
		order.clear();

		order.push_back(aux[4]);
		order.push_back(aux[5]);
		order.push_back(aux[1]);
		order.push_back(aux[0]);

		order.push_back(aux[7]);
		order.push_back(aux[6]);
		order.push_back(aux[2]);
		order.push_back(aux[3]);

	}

	if(rot[1]==1)
	{

		for(int i=0;i<8;i++) aux[i]=order[i];

		//order = {1,5,6,2,0,4,7,3};
		order.clear();

		order.push_back(aux[1]);
		order.push_back(aux[5]);
		order.push_back(aux[6]);
		order.push_back(aux[2]);

		order.push_back(aux[0]);
		order.push_back(aux[4]);
		order.push_back(aux[7]);
		order.push_back(aux[3]);


	}

	if(rot[1]==-1)
	{

		for(int i=0;i<8;i++) aux[i]=order[i];

		//order = {4,5,1,0,7,6,2,3};
		order.clear();

		order.push_back(aux[4]);
		order.push_back(aux[0]);
		order.push_back(aux[3]);
		order.push_back(aux[7]);

		order.push_back(aux[5]);
		order.push_back(aux[1]);
		order.push_back(aux[2]);
		order.push_back(aux[6]);

	}

	if(rot[2]==1)
	{

		for(int i=0;i<8;i++) aux[i]=order[i];

		//order = {3,0,1,2,7,4,5,6};
		order.clear();

		order.push_back(aux[3]);
		order.push_back(aux[0]);
		order.push_back(aux[1]);
		order.push_back(aux[2]);

		order.push_back(aux[7]);
		order.push_back(aux[4]);
		order.push_back(aux[5]);
		order.push_back(aux[6]);

	}

	if(rot[2]==-1)
	{

		for(int i=0;i<8;i++) aux[i]=order[i];

		//order = {1,2,3,0,5,6,7,4};
		order.clear();

		order.push_back(aux[1]);
		order.push_back(aux[2]);
		order.push_back(aux[3]);
		order.push_back(aux[0]);

		order.push_back(aux[5]);
		order.push_back(aux[6]);
		order.push_back(aux[7]);
		order.push_back(aux[4]);

	}

	if(sym[0]==1)
	{

		for(int i=0;i<8;i++) aux[i]=order[i];

		//order = {aux[1],aux[0],aux[3],aux[2],aux[5],aux[4],aux[7],aux[6]};
		order.clear();

		order.push_back(aux[1]);
		order.push_back(aux[0]);
		order.push_back(aux[3]);
		order.push_back(aux[2]);

		order.push_back(aux[5]);
		order.push_back(aux[4]);
		order.push_back(aux[7]);
		order.push_back(aux[6]);
	}

	if(sym[1]==1)
	{

		for(int i=0;i<8;i++) aux[i]=order[i];

		//order = {aux[3],aux[2],aux[6],aux[7],aux[0],aux[1],aux[5],aux[6]};
		order.clear();

		order.push_back(aux[3]);
		order.push_back(aux[2]);
		order.push_back(aux[1]);
		order.push_back(aux[0]);

		order.push_back(aux[7]);
		order.push_back(aux[6]);
		order.push_back(aux[5]);
		order.push_back(aux[4]);
	}

	if(sym[2]==1)
	{

		for(int i=0;i<8;i++) aux[i]=order[i];

		//order = {aux[4],aux[5],aux[6],aux[7],aux[0],aux[1],aux[2],aux[3]};
		order.clear();

		order.push_back(aux[4]);
		order.push_back(aux[5]);
		order.push_back(aux[6]);
		order.push_back(aux[7]);

		order.push_back(aux[0]);
		order.push_back(aux[1]);
		order.push_back(aux[2]);
		order.push_back(aux[3]);
	}

	return order;
																																																																																																																																																																																																																																																																																																																																																																																																																																																																																				}

unsigned node_shared_hash_fn(const void *v, const void *u)
{
	const shared_node_t *q = (const shared_node_t*) v;
	uint32_t a, b, c;

	a = (uint32_t) q->id;
	b = (uint32_t) 0;
	c = (uint32_t) 0;
	sc_hash_mix(a, b, c);
	sc_hash_final(a, b, c);
	return (unsigned) c;
}

int node_shared_equal_fn(const void *v, const void *u, const void *w)
{
	const shared_node_t *e1 = (const shared_node_t*) v;
	const shared_node_t *e2 = (const shared_node_t*) u;

	return (unsigned) (e1->id == e2->id);
}
 */
/*


void ApplyElement(hexa_tree_t* mesh, std::vector<double>& coords, int id, int iel, int* id_node,
		double* local_ref_x, double* local_ref_y, double* local_ref_z, std::vector<int>& ord, sc_hash_array_t* hash_nodes)
{

	int conn_p[8];
	GtsPoint* point[8]={NULL};
	GtsPoint* ref_point[8]={NULL};

	double cord_in_x[8],cord_in_y[8],cord_in_z[8];
	double ref_in_x[8], ref_in_y[8], ref_in_z[8];
	for (int ino = 0; ino < 8; ino++)
	{
		cord_in_x[ino]=coords[3*id_node[ino]+0] ;
		cord_in_y[ino]=coords[3*id_node[ino]+1] ;
		cord_in_z[ino]=coords[3*id_node[ino]+2] ;
		octant_node_t* node = (octant_node_t*) sc_array_index(&mesh->nodes, id_node[ino]);
		ref_in_x[ino] = node->x;
		ref_in_y[ino] = node->y;
		ref_in_z[ino] = node->z;
	}


	double cord_in_ref[3];
	cord_in_ref[0] = 0;
	cord_in_ref[1] = 0;
	cord_in_ref[2] = 0;
	for(int ino = 0; ino < 8; ino++)
	{
		cord_in_ref[0] = local_ref_x[ino];
		cord_in_ref[1] = local_ref_y[ino];
		cord_in_ref[2] = local_ref_z[ino];

		ref_point[ino] = LinearMapHex(cord_in_ref, ref_in_x, ref_in_y, ref_in_z);
		double auxx = ref_point[ino]->x;
		double auxy = ref_point[ino]->y;
		double auxz = ref_point[ino]->z;
		int x = round(auxx);
		int y = round(auxy);
		int z = round(auxz);
		point[ino] = LinearMapHex(cord_in_ref, cord_in_x,cord_in_y,cord_in_z);
		conn_p[ino] = -1;
		conn_p[ino] = AddPoint( mesh, hash_nodes, point[ino], coords, x, y, z);
	}

	int aux[8];
	for(int ino = 0; ino < 8; ino++) aux[ino] = conn_p[ino];
	for(int ino = 0; ino < 8; ino++) conn_p[ord[ino]] = aux[ino];

	if(iel == 0)
	{
		octant_t *elem1 = (octant_t*) sc_array_index(&mesh->elements, id);

		elem1->nodes[0].id = conn_p[0];
		elem1->nodes[1].id = conn_p[1];
		elem1->nodes[2].id = conn_p[2];
		elem1->nodes[3].id = conn_p[3];

		elem1->nodes[4].id = conn_p[4];
		elem1->nodes[5].id = conn_p[5];
		elem1->nodes[6].id = conn_p[6];
		elem1->nodes[7].id = conn_p[7];

		CopyPropEl(mesh,id,elem1);

		for(int ino = 0; ino < 8; ino++)
		{
			octant_node_t *node = (octant_node_t*) sc_array_index(&mesh->nodes, elem1->nodes[ino].id);
			elem1->nodes[ino].x = node->x;
			elem1->nodes[ino].y = node->y;
			elem1->nodes[ino].z = node->z;
		}

	}
	else
	{
		octant_t* elem2 = (octant_t*) sc_array_push(&mesh->elements);

		elem2->id = mesh->elements.elem_count;

		elem2->nodes[0].id = conn_p[0];
		elem2->nodes[1].id = conn_p[1];
		elem2->nodes[2].id = conn_p[2];
		elem2->nodes[3].id = conn_p[3];

		elem2->nodes[4].id = conn_p[4];
		elem2->nodes[5].id = conn_p[5];
		elem2->nodes[6].id = conn_p[6];
		elem2->nodes[7].id = conn_p[7];

		CopyPropEl(mesh,id,elem2);

		for(int ino = 0; ino < 8; ino++)
		{
			octant_node_t *node = (octant_node_t*) sc_array_index(&mesh->nodes, elem2->nodes[ino].id);
			elem2->nodes[ino].x = node->x;
			elem2->nodes[ino].y = node->y;
			elem2->nodes[ino].z = node->z;
		}
	}
}

void ApplyTemplate1(hexa_tree_t* mesh, std::vector<double>& coords, std::vector<int>& elements_ids, int iel, sc_hash_array_t* hash_nodes){

	sc_array_t toto;
	sc_array_init(&toto, sizeof(octant_t));

	octant_t* elemOrig = (octant_t*) sc_array_index(&mesh->elements, elements_ids[iel]);
	octant_t * elem = (octant_t*) sc_array_push(&toto);
	hexa_element_copy(elemOrig,elem);

	double step = double(2)/double(3);
	int id = elements_ids[iel];

	//reference element edge 0
	double local_ref[5][8][3];
	if(true){
		//element 0
		local_ref[0][0][0] = -1;
		local_ref[0][0][1] = -1;
		local_ref[0][0][2] = -1;

		local_ref[0][1][0] = -1+step;
		local_ref[0][1][1] = -1;
		local_ref[0][1][2] = -1;

		local_ref[0][2][0] = -1+step;
		local_ref[0][2][1] = -1+2*step;
		local_ref[0][2][2] = -1;

		local_ref[0][3][0] = -1;
		local_ref[0][3][1] = 1;
		local_ref[0][3][2] = -1;

		local_ref[0][4][0] = -1;
		local_ref[0][4][1] = -1;
		local_ref[0][4][2] = 1;

		local_ref[0][5][0] = -1+step;
		local_ref[0][5][1] = -1;
		local_ref[0][5][2] = -1+2*step;

		local_ref[0][6][0] = -1+step;
		local_ref[0][6][1] = -1+2*step;
		local_ref[0][6][2] = -1+2*step;

		local_ref[0][7][0] = -1;
		local_ref[0][7][1] = 1;
		local_ref[0][7][2] = 1;

		//element 1
		local_ref[1][0][0] = -1+step;
		local_ref[1][0][1] = -1;
		local_ref[1][0][2] = -1;

		local_ref[1][1][0] = -1+2*step;
		local_ref[1][1][1] = -1;
		local_ref[1][1][2] = -1;

		local_ref[1][2][0] = -1+2*step;
		local_ref[1][2][1] = -1+2*step;
		local_ref[1][2][2] = -1;

		local_ref[1][3][0] = -1+step;
		local_ref[1][3][1] = -1+2*step;
		local_ref[1][3][2] = -1;

		local_ref[1][4][0] = -1+step;
		local_ref[1][4][1] = -1;
		local_ref[1][4][2] = -1+2*step;

		local_ref[1][5][0] = -1+2*step;
		local_ref[1][5][1] = -1;
		local_ref[1][5][2] = -1+2*step;

		local_ref[1][6][0] = -1+2*step;
		local_ref[1][6][1] = -1+2*step;
		local_ref[1][6][2] = -1+2*step;

		local_ref[1][7][0] = -1+step;
		local_ref[1][7][1] = -1+2*step;
		local_ref[1][7][2] = -1+2*step;


		//element 2
		local_ref[2][0][0] = -1+2*step;
		local_ref[2][0][1] = -1;
		local_ref[2][0][2] = -1;

		local_ref[2][1][0] =  1;
		local_ref[2][1][1] = -1;
		local_ref[2][1][2] = -1;

		local_ref[2][2][0] =  1;
		local_ref[2][2][1] =  1;
		local_ref[2][2][2] = -1;

		local_ref[2][3][0] = -1+2*step;
		local_ref[2][3][1] = -1+2*step;
		local_ref[2][3][2] = -1;

		local_ref[2][4][0] = -1+2*step;
		local_ref[2][4][1] = -1;
		local_ref[2][4][2] = -1+2*step;

		local_ref[2][5][0] =  1;
		local_ref[2][5][1] = -1;
		local_ref[2][5][2] =  1;

		local_ref[2][6][0] = 1;
		local_ref[2][6][1] = 1;
		local_ref[2][6][2] = 1;

		local_ref[2][7][0] = -1+2*step;
		local_ref[2][7][1] = -1+2*step;
		local_ref[2][7][2] = -1+2*step;


		//element 3
		local_ref[3][0][0] = -1+1*step;
		local_ref[3][0][1] = -1+2*step;
		local_ref[3][0][2] = -1;

		local_ref[3][1][0] = -1+2*step;
		local_ref[3][1][1] = -1+2*step;
		local_ref[3][1][2] = -1;

		local_ref[3][2][0] =  1;
		local_ref[3][2][1] =  1;
		local_ref[3][2][2] = -1;

		local_ref[3][3][0] = -1;
		local_ref[3][3][1] = 1;
		local_ref[3][3][2] = -1;

		local_ref[3][4][0] = -1+1*step;
		local_ref[3][4][1] = -1+2*step;
		local_ref[3][4][2] = -1+2*step;

		local_ref[3][5][0] = -1+2*step;
		local_ref[3][5][1] = -1+2*step;
		local_ref[3][5][2] = -1+2*step;

		local_ref[3][6][0] = 1;
		local_ref[3][6][1] = 1;
		local_ref[3][6][2] = 1;

		local_ref[3][7][0] = -1;
		local_ref[3][7][1] = 1;
		local_ref[3][7][2] = 1;

		//element 4
		local_ref[4][0][0] = -1+step;
		local_ref[4][0][1] = -1;
		local_ref[4][0][2] = -1+2*step;

		local_ref[4][1][0] = -1+2*step;
		local_ref[4][1][1] = -1;
		local_ref[4][1][2] = -1+2*step;

		local_ref[4][2][0] = -1+2*step;
		local_ref[4][2][1] = -1+2*step;
		local_ref[4][2][2] = -1+2*step;

		local_ref[4][3][0] = -1+step;
		local_ref[4][3][1] = -1+2*step;
		local_ref[4][3][2] = -1+2*step;

		local_ref[4][4][0] = -1;
		local_ref[4][4][1] = -1;
		local_ref[4][4][2] = 1;

		local_ref[4][5][0] = 1;
		local_ref[4][5][1] = -1;
		local_ref[4][5][2] = 1;

		local_ref[4][6][0] = 1;
		local_ref[4][6][1] = 1;
		local_ref[4][6][2] = 1;

		local_ref[4][7][0] = -1;
		local_ref[4][7][1] = 1;
		local_ref[4][7][2] = 1;
	}

	//define the rotation of the reference element
	int rot[3];
	int sym[3];
	int id_node[8];
	vector<int> ord;

	if(elem->pad==10){
		//edge 0
		rot[0] = 0;
		rot[1] = 0;
		rot[2] = 0;

		sym[0] = 0;
		sym[1] = 0;
		sym[2] = 0;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}

	}else if(elem->pad==11){
		//edge 1
		rot[0] = 0;
		rot[1] = 0;
		rot[2] = -1;

		sym[0] = 0;
		sym[1] = 0;
		sym[2] = 0;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}

	}else if(elem->pad==12){
		//edge 2
		rot[0] = 0;
		rot[1] = 0;
		rot[2] = 0;

		sym[0] = 0;
		sym[1] = 1;
		sym[2] = 0;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}

	}else if(elem->pad==13){
		//edge 3
		rot[0] = 0;
		rot[1] = 0;
		rot[2] = 1;

		sym[0] = 0;
		sym[1] = 0;
		sym[2] = 0;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}

	}else if(elem->pad==14){
		//edge 4
		rot[0] = 0;
		rot[1] = -1;
		rot[2] = 0;

		sym[0] = 0;
		sym[1] = 0;
		sym[2] = 0;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}

	}else if(elem->pad==15){
		//edge 5
		rot[0] = 0;
		rot[1] = 1;
		rot[2] = 0;

		sym[0] = 0;
		sym[1] = 0;
		sym[2] = 0;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}

	}else if(elem->pad==16){
		//edge 6
		rot[0] = 0;
		rot[1] = 1;
		rot[2] = 0;

		sym[0] = 0;
		sym[1] = 1;
		sym[2] = 0;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}

	}else if(elem->pad==17){
		//edge 7
		rot[0] = 0;
		rot[1] = -1;
		rot[2] = 0;

		sym[0] = 0;
		sym[1] = 1;
		sym[2] = 0;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}

	}else if(elem->pad==18){
		//edge 8
		rot[0] = 0;
		rot[1] = 0;
		rot[2] = 0;

		sym[0] = 0;
		sym[1] = 0;
		sym[2] = 1;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}

	}else if(elem->pad==19){
		//edge 9
		rot[0] = 0;
		rot[1] = 0;
		rot[2] = -1;

		sym[0] = 0;
		sym[1] = 0;
		sym[2] = 1;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}

	}else if(elem->pad==20){
		//edge 10
		rot[0] = -1;
		rot[1] = 0;
		rot[2] = 0;

		sym[0] = 0;
		sym[1] = 1;
		sym[2] = 0;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}

	}else if(elem->pad==21){
		//edge 11
		rot[0] = 0;
		rot[1] = 0;
		rot[2] = 1;

		sym[0] = 0;
		sym[1] = 0;
		sym[2] = 1;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}

	}else{
		printf ("Error in template 1\n");
		exit (EXIT_FAILURE);
	}


	for(int iel = 0; iel < 5; iel++){
		double local_ref_x[8];
		double local_ref_y[8];
		double local_ref_z[8];
		for(int ino = 0; ino < 8; ino++){
			local_ref_x[ino] = local_ref[iel][ino][0];
			local_ref_y[ino] = local_ref[iel][ino][1];
			local_ref_z[ino] = local_ref[iel][ino][2];
		}
		ApplyElement(mesh, coords,  id,  iel,  id_node,  local_ref_x, local_ref_y, local_ref_z,  ord,  hash_nodes);
	}

	sc_array_reset(&toto);
}

void ApplyTemplate2(hexa_tree_t* mesh, std::vector<double>& coords, std::vector<int>& elements_ids, int iel, sc_hash_array_t* hash_nodes){

	sc_array_t toto;
	sc_array_init(&toto, sizeof(octant_t));

	octant_t* elemOrig = (octant_t*) sc_array_index(&mesh->elements, elements_ids[iel]);
	octant_t * elem = (octant_t*) sc_array_push(&toto);
	hexa_element_copy(elemOrig,elem);

	double step = double(2)/double(3);
	int id = elements_ids[iel];

	double cord_in_ref[3];
	cord_in_ref[0] = 0;
	cord_in_ref[1] = 0;
	cord_in_ref[2] = 0;

	//reference element edge 0
	double local_ref[3][8][3];
	if(true){
		//element 0
		local_ref[0][0][0] = -1;
		local_ref[0][0][1] = -1;
		local_ref[0][0][2] = -1;

		local_ref[0][1][0] = 1;
		local_ref[0][1][1] = -1;
		local_ref[0][1][2] = -1;

		local_ref[0][2][0] = 1;
		local_ref[0][2][1] = 1;
		local_ref[0][2][2] = -1;

		local_ref[0][3][0] = -1;
		local_ref[0][3][1] = 1;
		local_ref[0][3][2] = -1;

		local_ref[0][4][0] = -1;
		local_ref[0][4][1] = -1;
		local_ref[0][4][2] = -1+step;

		local_ref[0][5][0] = 1;
		local_ref[0][5][1] = -1;
		local_ref[0][5][2] = -1+step;

		local_ref[0][6][0] = 1;
		local_ref[0][6][1] = 1;
		local_ref[0][6][2] = -1+step;

		local_ref[0][7][0] = -1;
		local_ref[0][7][1] = 1;
		local_ref[0][7][2] = -1+step;

		//element 1
		local_ref[1][0][0] = -1;
		local_ref[1][0][1] = -1;
		local_ref[1][0][2] = -1+step;

		local_ref[1][1][0] = 1;
		local_ref[1][1][1] = -1;
		local_ref[1][1][2] = -1+step;

		local_ref[1][2][0] = 1;
		local_ref[1][2][1] = 1;
		local_ref[1][2][2] = -1+step;

		local_ref[1][3][0] = -1;
		local_ref[1][3][1] = 1;
		local_ref[1][3][2] = -1+step;

		local_ref[1][4][0] = -1;
		local_ref[1][4][1] = -1;
		local_ref[1][4][2] = -1+2*step;

		local_ref[1][5][0] = 1;
		local_ref[1][5][1] = -1;
		local_ref[1][5][2] = -1+2*step;

		local_ref[1][6][0] = 1;
		local_ref[1][6][1] = 1;
		local_ref[1][6][2] = -1+2*step;

		local_ref[1][7][0] = -1;
		local_ref[1][7][1] = 1;
		local_ref[1][7][2] = -1+2*step;


		//element 2
		local_ref[2][0][0] = -1;
		local_ref[2][0][1] = -1;
		local_ref[2][0][2] = -1+2*step;

		local_ref[2][1][0] =  1;
		local_ref[2][1][1] = -1;
		local_ref[2][1][2] = -1+2*step;

		local_ref[2][2][0] =  1;
		local_ref[2][2][1] =  1;
		local_ref[2][2][2] = -1+2*step;

		local_ref[2][3][0] = -1;
		local_ref[2][3][1] = 1;
		local_ref[2][3][2] = -1+2*step;

		local_ref[2][4][0] = -1;
		local_ref[2][4][1] = -1;
		local_ref[2][4][2] = 1;

		local_ref[2][5][0] =  1;
		local_ref[2][5][1] = -1;
		local_ref[2][5][2] =  1;

		local_ref[2][6][0] = 1;
		local_ref[2][6][1] = 1;
		local_ref[2][6][2] = 1;

		local_ref[2][7][0] = -1;
		local_ref[2][7][1] = 1;
		local_ref[2][7][2] = 1;
	}

	//define the rotation of the reference element
	int rot[3];
	int sym[3];
	int id_node[8];
	vector<int> ord;

	if(elem->pad==24){
		//edge 4 5 6 7
		rot[0] = 0;
		rot[1] = 0;
		rot[2] = 0;

		sym[0] = 0;
		sym[1] = 0;
		sym[2] = 0;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}

	}else if(elem->pad==23){
		//edge 1 3 9 11
		rot[0] = 1;
		rot[1] = 0;
		rot[2] = 0;

		sym[0] = 0;
		sym[1] = 0;
		sym[2] = 0;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}

	}else if(elem->pad==22){
		//edge 0 2 8 10
		rot[0] = 0;
		rot[1] = 1;
		rot[2] = 0;

		sym[0] = 0;
		sym[1] = 0;
		sym[2] = 0;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}

	}else{
		printf ("Error in template 2\n");
		exit (EXIT_FAILURE);
	}


	for(int iel = 0; iel < 3; iel++){
		double local_ref_x[8];
		double local_ref_y[8];
		double local_ref_z[8];
		for(int ino = 0; ino < 8; ino++){
			local_ref_x[ino] = local_ref[iel][ino][0];
			local_ref_y[ino] = local_ref[iel][ino][1];
			local_ref_z[ino] = local_ref[iel][ino][2];
		}
		ApplyElement(mesh, coords,  id,  iel,  id_node,  local_ref_x, local_ref_y, local_ref_z,  ord,  hash_nodes);
	}
	sc_array_reset(&toto);
}

void ApplyTemplate3(hexa_tree_t* mesh, std::vector<double>& coords, std::vector<int>& elements_ids, int iel, sc_hash_array_t* hash_nodes){

	sc_array_t toto;
	sc_array_init(&toto, sizeof(octant_t));

	octant_t* elemOrig = (octant_t*) sc_array_index(&mesh->elements, elements_ids[iel]);
	octant_t * elem = (octant_t*) sc_array_push(&toto);
	hexa_element_copy(elemOrig,elem);

	double step = double(2)/double(3);
	int id = elements_ids[iel];

	double cord_in_ref[3];
	cord_in_ref[0] = 0;
	cord_in_ref[1] = 0;
	cord_in_ref[2] = 0;

	//reference element edge 0
	double local_ref[4][8][3];
	if(true){
		//element 0
		local_ref[0][0][0] = -1;
		local_ref[0][0][1] = -1;
		local_ref[0][0][2] = -1;

		local_ref[0][1][0] = 1;
		local_ref[0][1][1] = -1;
		local_ref[0][1][2] = -1;

		local_ref[0][2][0] = 1;
		local_ref[0][2][1] = -1+1*step;
		local_ref[0][2][2] = -1;

		local_ref[0][3][0] = -1;
		local_ref[0][3][1] = -1+1*step;
		local_ref[0][3][2] = -1;

		local_ref[0][4][0] = -1;
		local_ref[0][4][1] = -1;
		local_ref[0][4][2] = 1;

		local_ref[0][5][0] = 1;
		local_ref[0][5][1] = -1;
		local_ref[0][5][2] = 1;

		local_ref[0][6][0] = 1;
		local_ref[0][6][1] = -1+1*step;
		local_ref[0][6][2] = -1+2*step;

		local_ref[0][7][0] = -1;
		local_ref[0][7][1] = -1+1*step;
		local_ref[0][7][2] = -1+2*step;

		//element 1
		local_ref[1][0][0] = -1;
		local_ref[1][0][1] = -1+step;
		local_ref[1][0][2] = -1;

		local_ref[1][1][0] = 1;
		local_ref[1][1][1] = -1+step;
		local_ref[1][1][2] = -1;

		local_ref[1][2][0] = 1;
		local_ref[1][2][1] = -1+2*step;
		local_ref[1][2][2] = -1;

		local_ref[1][3][0] = -1;
		local_ref[1][3][1] = -1+2*step;
		local_ref[1][3][2] = -1;

		local_ref[1][4][0] = -1;
		local_ref[1][4][1] = -1+step;
		local_ref[1][4][2] = -1+2*step;

		local_ref[1][5][0] = 1;
		local_ref[1][5][1] = -1+step;
		local_ref[1][5][2] = -1+2*step;

		local_ref[1][6][0] = 1;
		local_ref[1][6][1] = -1+2*step;
		local_ref[1][6][2] = -1+2*step;

		local_ref[1][7][0] = -1;
		local_ref[1][7][1] = -1+2*step;
		local_ref[1][7][2] = -1+2*step;


		//element 2
		local_ref[2][0][0] = -1;
		local_ref[2][0][1] = -1+2*step;
		local_ref[2][0][2] = -1;

		local_ref[2][1][0] =  1;
		local_ref[2][1][1] = -1+2*step;
		local_ref[2][1][2] = -1;

		local_ref[2][2][0] =  1;
		local_ref[2][2][1] =  1;
		local_ref[2][2][2] = -1;

		local_ref[2][3][0] = -1;
		local_ref[2][3][1] =  1;
		local_ref[2][3][2] = -1;

		local_ref[2][4][0] = -1;
		local_ref[2][4][1] = -1+2*step;
		local_ref[2][4][2] = -1+2*step;

		local_ref[2][5][0] =  1;
		local_ref[2][5][1] = -1+2*step;
		local_ref[2][5][2] =  -1+2*step;

		local_ref[2][6][0] = 1;
		local_ref[2][6][1] = 1;
		local_ref[2][6][2] = 1;

		local_ref[2][7][0] = -1;
		local_ref[2][7][1] = 1;
		local_ref[2][7][2] = 1;

		//element 3
		local_ref[3][0][0] = -1;
		local_ref[3][0][1] = -1+step;
		local_ref[3][0][2] = -1+2*step;

		local_ref[3][1][0] =  1;
		local_ref[3][1][1] = -1+step;
		local_ref[3][1][2] = -1+2*step;

		local_ref[3][2][0] =  1;
		local_ref[3][2][1] = -1+2*step;
		local_ref[3][2][2] = -1+2*step;

		local_ref[3][3][0] = -1;
		local_ref[3][3][1] = -1+2*step;
		local_ref[3][3][2] = -1+2*step;

		local_ref[3][4][0] = -1;
		local_ref[3][4][1] = -1;
		local_ref[3][4][2] = 1;

		local_ref[3][5][0] =  1;
		local_ref[3][5][1] = -1;
		local_ref[3][5][2] =  1;

		local_ref[3][6][0] = 1;
		local_ref[3][6][1] = 1;
		local_ref[3][6][2] = 1;

		local_ref[3][7][0] = -1;
		local_ref[3][7][1] = 1;
		local_ref[3][7][2] = 1;
	}

	//define the rotation of the reference element
	int rot[3];
	int sym[3];
	int id_node[8];
	vector<int> ord;

	if(elem->pad==25){
		//edge 0 2
		rot[0] = 0;
		rot[1] = 0;
		rot[2] = 1;

		sym[0] = 0;
		sym[1] = 0;
		sym[2] = 0;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}

	}else if(elem->pad==26){
		//edge 1 3
		rot[0] = 0;
		rot[1] = 0;
		rot[2] = 0;

		sym[0] = 0;
		sym[1] = 0;
		sym[2] = 0;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}

	}else if(elem->pad==27){
		//edge 8 10
		rot[0] = 0;
		rot[1] = 0;
		rot[2] = 0;

		sym[0] = 0;
		sym[1] = 0;
		sym[2] = 1;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}

	}else if(elem->pad==28){
		//edge 9 11
		rot[0] = 0;
		rot[1] = 0;
		rot[2] = 1;

		sym[0] = 0;
		sym[1] = 0;
		sym[2] = 1;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}

	}else if(elem->pad==29){
		//edge 4 5
		rot[0] = 1;
		rot[1] = 0;
		rot[2] = 0;

		sym[0] = 0;
		sym[1] = 0;
		sym[2] = 0;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}

	}else if(elem->pad==30){
		//edge 5 6
		rot[0] = -1;
		rot[1] = 1;
		rot[2] = 0;

		sym[0] = 0;
		sym[1] = 0;
		sym[2] = 0;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}

	}else if(elem->pad==31){
		//edge 6 7
		rot[0] = -1;
		rot[1] = 0;
		rot[2] = 0;

		sym[0] = 0;
		sym[1] = 0;
		sym[2] = 0;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}

	}else if(elem->pad==32){
		//edge 7 4
		rot[0] = -1;
		rot[1] = -1;
		rot[2] = 0;

		sym[0] = 0;
		sym[1] = 0;
		sym[2] = 0;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}

	}else if(elem->pad==33){
		//edge 0 8
		rot[0] = 1;
		rot[1] = 0;
		rot[2] = 1;

		sym[0] = 0;
		sym[1] = 0;
		sym[2] = 0;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}

	}else if(elem->pad==34){
		//edge 1 9
		rot[0] = 0;
		rot[1] = 1;
		rot[2] = 0;

		sym[0] = 0;
		sym[1] = 0;
		sym[2] = 0;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}

	}else if(elem->pad==35){
		//edge 2 10
		rot[0] = -1;
		rot[1] = 0;
		rot[2] = 1;

		sym[0] = 0;
		sym[1] = 0;
		sym[2] = 0;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}

	}else if(elem->pad==36){
		//edge 3 11
		rot[0] = 0;
		rot[1] = -1;
		rot[2] = 0;

		sym[0] = 0;
		sym[1] = 0;
		sym[2] = 0;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}

	}else{
		printf ("Error in template 3\n");
		exit (EXIT_FAILURE);
	}


	for(int iel = 0; iel < 4; iel++){
		double local_ref_x[8];
		double local_ref_y[8];
		double local_ref_z[8];
		for(int ino = 0; ino < 8; ino++){
			local_ref_x[ino] = local_ref[iel][ino][0];
			local_ref_y[ino] = local_ref[iel][ino][1];
			local_ref_z[ino] = local_ref[iel][ino][2];
		}
		ApplyElement(mesh, coords,  id,  iel,  id_node,  local_ref_x, local_ref_y, local_ref_z,  ord,  hash_nodes);
	}
	sc_array_reset(&toto);
}

void ApplyTemplate4(hexa_tree_t* mesh, std::vector<double>& coords, std::vector<int>& elements_ids, int iel, sc_hash_array_t* hash_nodes){

	sc_array_t toto;
	sc_array_init(&toto, sizeof(octant_t));

	octant_t* elemOrig = (octant_t*) sc_array_index(&mesh->elements, elements_ids[iel]);
	octant_t * elem = (octant_t*) sc_array_push(&toto);
	hexa_element_copy(elemOrig,elem);

	double step = double(2)/double(3);
	int id = elements_ids[iel];

	double cord_in_ref[3];
	cord_in_ref[0] = 0;
	cord_in_ref[1] = 0;
	cord_in_ref[2] = 0;

	//reference element edge 0
	double local_ref[5][8][3];
	if(true){
		//element 0
		local_ref[0][0][0] = -1;
		local_ref[0][0][1] = -1;
		local_ref[0][0][2] = -1;

		local_ref[0][1][0] = -1+step;
		local_ref[0][1][1] = -1;
		local_ref[0][1][2] = -1;

		local_ref[0][2][0] = -1+step;
		local_ref[0][2][1] = -1+step;
		local_ref[0][2][2] = -1;

		local_ref[0][3][0] = -1;
		local_ref[0][3][1] = -1+step;
		local_ref[0][3][2] = -1;

		local_ref[0][4][0] = -1;
		local_ref[0][4][1] = -1;
		local_ref[0][4][2] =  1;

		local_ref[0][5][0] = -1+step;
		local_ref[0][5][1] = -1;
		local_ref[0][5][2] = 1;

		local_ref[0][6][0] = -1+step;
		local_ref[0][6][1] = -1+step;
		local_ref[0][6][2] = 1;

		local_ref[0][7][0] = -1;
		local_ref[0][7][1] = -1+step;
		local_ref[0][7][2] =  1;

		//element 1
		local_ref[1][0][0] = -1+step;
		local_ref[1][0][1] = -1;
		local_ref[1][0][2] = -1;

		local_ref[1][1][0] = -1+2*step;
		local_ref[1][1][1] = -1;
		local_ref[1][1][2] = -1;

		local_ref[1][2][0] = -1+2*step;
		local_ref[1][2][1] = -1+2*step;
		local_ref[1][2][2] = -1;

		local_ref[1][3][0] = -1+step;
		local_ref[1][3][1] = -1+1*step;
		local_ref[1][3][2] = -1;

		local_ref[1][4][0] = -1+step;
		local_ref[1][4][1] = -1;
		local_ref[1][4][2] = 1;

		local_ref[1][5][0] = -1+2*step;
		local_ref[1][5][1] = -1;
		local_ref[1][5][2] = 1;

		local_ref[1][6][0] = -1+2*step;
		local_ref[1][6][1] = -1+2*step;
		local_ref[1][6][2] = 1;

		local_ref[1][7][0] = -1+step;
		local_ref[1][7][1] = -1+1*step;
		local_ref[1][7][2] = 1;


		//element 2
		local_ref[2][0][0] = -1;
		local_ref[2][0][1] = -1+step;
		local_ref[2][0][2] = -1;

		local_ref[2][1][0] = -1+step;
		local_ref[2][1][1] = -1+step;
		local_ref[2][1][2] = -1;

		local_ref[2][2][0] = -1+2*step;
		local_ref[2][2][1] = -1+2*step;
		local_ref[2][2][2] = -1;

		local_ref[2][3][0] = -1;
		local_ref[2][3][1] = -1+2*step;
		local_ref[2][3][2] = -1;

		local_ref[2][4][0] = -1;
		local_ref[2][4][1] = -1+step;
		local_ref[2][4][2] = 1;

		local_ref[2][5][0] =  -1+step;
		local_ref[2][5][1] = -1+step;
		local_ref[2][5][2] =  1;

		local_ref[2][6][0] = -1+2*step;
		local_ref[2][6][1] = -1+2*step;
		local_ref[2][6][2] = 1;

		local_ref[2][7][0] = -1;
		local_ref[2][7][1] = -1+2*step;
		local_ref[2][7][2] = 1;


		//element 3
		local_ref[3][0][0] = -1+2*step;
		local_ref[3][0][1] = -1;
		local_ref[3][0][2] = -1;

		local_ref[3][1][0] =  1;
		local_ref[3][1][1] = -1;
		local_ref[3][1][2] = -1;

		local_ref[3][2][0] =  1;
		local_ref[3][2][1] =  1;
		local_ref[3][2][2] = -1;

		local_ref[3][3][0] = -1+2*step;
		local_ref[3][3][1] = -1+2*step;
		local_ref[3][3][2] = -1;

		local_ref[3][4][0] = -1+2*step;
		local_ref[3][4][1] = -1;
		local_ref[3][4][2] =  1;

		local_ref[3][5][0] = 1;
		local_ref[3][5][1] = -1;
		local_ref[3][5][2] = 1;

		local_ref[3][6][0] = 1;
		local_ref[3][6][1] = 1;
		local_ref[3][6][2] = 1;

		local_ref[3][7][0] = -1+2*step;
		local_ref[3][7][1] = -1+2*step;
		local_ref[3][7][2] = 1;

		//element 4
		local_ref[4][0][0] = -1;
		local_ref[4][0][1] = -1+2*step;
		local_ref[4][0][2] = -1;

		local_ref[4][1][0] = -1+2*step;
		local_ref[4][1][1] = -1+2*step;
		local_ref[4][1][2] = -1;

		local_ref[4][2][0] = 1;
		local_ref[4][2][1] = 1;
		local_ref[4][2][2] = -1;

		local_ref[4][3][0] = -1;
		local_ref[4][3][1] = 1;
		local_ref[4][3][2] = -1;

		local_ref[4][4][0] = -1;
		local_ref[4][4][1] = -1+2*step;
		local_ref[4][4][2] = 1;

		local_ref[4][5][0] = -1+2*step;
		local_ref[4][5][1] = -1+2*step;
		local_ref[4][5][2] = 1;

		local_ref[4][6][0] = 1;
		local_ref[4][6][1] = 1;
		local_ref[4][6][2] = 1;

		local_ref[4][7][0] = -1;
		local_ref[4][7][1] = 1;
		local_ref[4][7][2] = 1;
	}

	//define the rotation of the reference element
	int rot[3];
	int sym[3];
	int id_node[8];
	vector<int> ord;

	if(elem->pad== 37){
		//edge 0 1 8 9
		rot[0] = 0;
		rot[1] = 0;
		rot[2] = 0;

		sym[0] = 1;
		sym[1] = 0;
		sym[2] = 0;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}

	}else if(elem->pad==38){
		//edge 1 2 9 10
		rot[0] = 0;
		rot[1] = 0;
		rot[2] = -1;

		sym[0] = 1;
		sym[1] = 0;
		sym[2] = 0;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}

	}else if(elem->pad==39){
		//edge 2 3 10 11
		rot[0] = 0;
		rot[1] = 0;
		rot[2] = 1;

		sym[0] = 0;
		sym[1] = 0;
		sym[2] = 0;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}

	}else if(elem->pad==40){
		//edge 0 3 8 11
		rot[0] = 0;
		rot[1] = 0;
		rot[2] = 0;

		sym[0] = 0;
		sym[1] = 0;
		sym[2] = 0;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}

	}else if(elem->pad==41){
		//edge 0 4 2 7
		rot[0] = 1;
		rot[1] = 0;
		rot[2] = 0;

		sym[0] = 0;
		sym[1] = 1;
		sym[2] = 0;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}

	}else if(elem->pad==42){
		//edge 0 5 2 6
		rot[0] = 1;
		rot[1] = 0;
		rot[2] = 0;

		sym[0] = 1;
		sym[1] = 1;
		sym[2] = 0;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}

	}else if(elem->pad==43){
		//edge 5 8 6 10
		rot[0] = 1;
		rot[1] = 0;
		rot[2] = 0;

		sym[0] = 1;
		sym[1] = 0;
		sym[2] = 0;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}

	}else if(elem->pad==44){
		//edge 4 8 7 10
		rot[0] = 1;
		rot[1] = 0;
		rot[2] = 0;

		sym[0] = 0;
		sym[1] = 0;
		sym[2] = 0;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}

	}else if(elem->pad==45){
		//edge 3 4  1  5
		rot[0] = 0;
		rot[1] = 1;
		rot[2] = 0;

		sym[0] = 0;
		sym[1] = 0;
		sym[2] = 0;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}

	}else if(elem->pad==46){
		//edge 1 6 3 7
		rot[0] = 0;
		rot[1] = 1;
		rot[2] = 0;

		sym[0] = 0;
		sym[1] = 1;
		sym[2] = 0;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}

	}else if(elem->pad==47){
		//edge 9 6  7 11
		rot[0] = 0;
		rot[1] = 1;
		rot[2] = 0;

		sym[0] = 1;
		sym[1] = 1;
		sym[2] = 0;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}

	}else if(elem->pad==48){
		//edge 4 5 9 11
		rot[0] = 0;
		rot[1] = 1;
		rot[2] = 0;

		sym[0] = 1;
		sym[1] = 0;
		sym[2] = 0;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}

	}else{
		printf ("Error in template 4\n");
		exit (EXIT_FAILURE);
	}


	for(int iel = 0; iel < 5; iel++){
		double local_ref_x[8];
		double local_ref_y[8];
		double local_ref_z[8];
		for(int ino = 0; ino < 8; ino++){
			local_ref_x[ino] = local_ref[iel][ino][0];
			local_ref_y[ino] = local_ref[iel][ino][1];
			local_ref_z[ino] = local_ref[iel][ino][2];
		}
		ApplyElement(mesh, coords,  id,  iel,  id_node,  local_ref_x, local_ref_y, local_ref_z,  ord,  hash_nodes);
	}
	sc_array_reset(&toto);
}

void ApplyTemplate5(hexa_tree_t* mesh, std::vector<double>& coords, std::vector<int>& elements_ids, int iel, sc_hash_array_t* hash_nodes){

	sc_array_t toto;
	sc_array_init(&toto, sizeof(octant_t));

	octant_t* elemOrig = (octant_t*) sc_array_index(&mesh->elements, elements_ids[iel]);
	octant_t * elem = (octant_t*) sc_array_push(&toto);
	hexa_element_copy(elemOrig,elem);

	double step = double(2)/double(3);
	int id = elements_ids[iel];

	double cord_in_ref[3];
	cord_in_ref[0] = 0;
	cord_in_ref[1] = 0;
	cord_in_ref[2] = 0;

	//reference element edge 0
	double local_ref[13][8][3];
	if(true){
		//element 0
		local_ref[0][0][0] = -1;
		local_ref[0][0][1] = -1;
		local_ref[0][0][2] = -1;

		local_ref[0][1][0] = 1;
		local_ref[0][1][1] = -1;
		local_ref[0][1][2] = -1;

		local_ref[0][2][0] = 1;
		local_ref[0][2][1] = 1;
		local_ref[0][2][2] = -1;

		local_ref[0][3][0] = -1;
		local_ref[0][3][1] = 1;
		local_ref[0][3][2] = -1;

		local_ref[0][4][0] = -1;
		local_ref[0][4][1] = -1+step;
		local_ref[0][4][2] = -1+step;

		local_ref[0][5][0] = 1;
		local_ref[0][5][1] = -1+step;
		local_ref[0][5][2] = -1+step;

		local_ref[0][6][0] = 1;
		local_ref[0][6][1] = -1+2*step;
		local_ref[0][6][2] = -1+step;

		local_ref[0][7][0] = -1;
		local_ref[0][7][1] = -1+2*step;
		local_ref[0][7][2] = -1+step;

		//element 1
		local_ref[1][0][0] = -1;
		local_ref[1][0][1] = -1+step;
		local_ref[1][0][2] = -1+step;

		local_ref[1][1][0] = 1;
		local_ref[1][1][1] = -1+step;
		local_ref[1][1][2] = -1+step;

		local_ref[1][2][0] = 1;
		local_ref[1][2][1] = -1+2*step;
		local_ref[1][2][2] = -1+step;

		local_ref[1][3][0] = -1;
		local_ref[1][3][1] = -1+2*step;
		local_ref[1][3][2] = -1+step;

		local_ref[1][4][0] = -1+step;
		local_ref[1][4][1] = -1+step;
		local_ref[1][4][2] = -1+2*step;

		local_ref[1][5][0] = -1+2*step;
		local_ref[1][5][1] = -1+step;
		local_ref[1][5][2] = -1+2*step;

		local_ref[1][6][0] = -1+2*step;
		local_ref[1][6][1] = -1+2*step;
		local_ref[1][6][2] = -1+2*step;

		local_ref[1][7][0] = -1+step;
		local_ref[1][7][1] = -1+2*step;
		local_ref[1][7][2] = -1+2*step;

		//element 2
		local_ref[2][0][0] = -1+step;
		local_ref[2][0][1] = -1+step;
		local_ref[2][0][2] = -1+2*step;

		local_ref[2][1][0] = -1+2*step;
		local_ref[2][1][1] = -1+step;
		local_ref[2][1][2] = -1+2*step;

		local_ref[2][2][0] = -1+2*step;
		local_ref[2][2][1] = -1+2*step;
		local_ref[2][2][2] = -1+2*step;

		local_ref[2][3][0] = -1+step;
		local_ref[2][3][1] = -1+2*step;
		local_ref[2][3][2] = -1+2*step;

		local_ref[2][4][0] = -1+step;
		local_ref[2][4][1] = -1+step;
		local_ref[2][4][2] = 1;

		local_ref[2][5][0] = -1+2*step;
		local_ref[2][5][1] = -1+step;
		local_ref[2][5][2] = 1;

		local_ref[2][6][0] = -1+2*step;
		local_ref[2][6][1] = -1+2*step;
		local_ref[2][6][2] = 1;

		local_ref[2][7][0] = -1+step;
		local_ref[2][7][1] = -1+2*step;
		local_ref[2][7][2] = 1;

		//element 3
		local_ref[3][0][0] = -1;
		local_ref[3][0][1] = -1;
		local_ref[3][0][2] = -1;

		local_ref[3][1][0] =  1;
		local_ref[3][1][1] = -1;
		local_ref[3][1][2] = -1;

		local_ref[3][2][0] =  1;
		local_ref[3][2][1] = -1+step;
		local_ref[3][2][2] = -1+step;

		local_ref[3][3][0] = -1;
		local_ref[3][3][1] = -1+step;
		local_ref[3][3][2] = -1+step;

		local_ref[3][4][0] = -1+step;
		local_ref[3][4][1] = -1;
		local_ref[3][4][2] = -1+step;

		local_ref[3][5][0] = -1+2*step;
		local_ref[3][5][1] = -1;
		local_ref[3][5][2] = -1+step;

		local_ref[3][6][0] = -1+2*step;
		local_ref[3][6][1] = -1+1*step;
		local_ref[3][6][2] = -1+2*step;

		local_ref[3][7][0] = -1+step;
		local_ref[3][7][1] = -1+1*step;
		local_ref[3][7][2] = -1+2*step;

		//element 4
		local_ref[4][0][0] = -1;
		local_ref[4][0][1] = -1+2*step;
		local_ref[4][0][2] = -1+step;

		local_ref[4][1][0] =  1;
		local_ref[4][1][1] = -1+2*step;
		local_ref[4][1][2] = -1+step;

		local_ref[4][2][0] =  1;
		local_ref[4][2][1] =  1;
		local_ref[4][2][2] = -1;

		local_ref[4][3][0] = -1;
		local_ref[4][3][1] =  1;
		local_ref[4][3][2] = -1;

		local_ref[4][4][0] = -1+step;
		local_ref[4][4][1] = -1+2*step;
		local_ref[4][4][2] = -1+2*step;

		local_ref[4][5][0] = -1+2*step;
		local_ref[4][5][1] = -1+2*step;
		local_ref[4][5][2] = -1+2*step;

		local_ref[4][6][0] = -1+2*step;
		local_ref[4][6][1] = 1;
		local_ref[4][6][2] = -1+step;

		local_ref[4][7][0] = -1+step;
		local_ref[4][7][1] = 1;
		local_ref[4][7][2] = -1+step;

		//element 5
		local_ref[5][0][0] = -1;
		local_ref[5][0][1] = -1;
		local_ref[5][0][2] = -1;

		local_ref[5][1][0] = -1+step;
		local_ref[5][1][1] = -1;
		local_ref[5][1][2] = -1+step;

		local_ref[5][2][0] = -1+step;
		local_ref[5][2][1] = -1+step;
		local_ref[5][2][2] = -1+2*step;

		local_ref[5][3][0] = -1;
		local_ref[5][3][1] = -1+step;
		local_ref[5][3][2] = -1+step;

		local_ref[5][4][0] = -1;
		local_ref[5][4][1] = -1;
		local_ref[5][4][2] = 1;

		local_ref[5][5][0] = -1+step;
		local_ref[5][5][1] = -1;
		local_ref[5][5][2] = 1;

		local_ref[5][6][0] = -1+step;
		local_ref[5][6][1] = -1+step;
		local_ref[5][6][2] = 1;

		local_ref[5][7][0] = -1;
		local_ref[5][7][1] = -1+step;
		local_ref[5][7][2] = 1;

		//element 6
		local_ref[6][0][0] = -1+step;
		local_ref[6][0][1] = -1;
		local_ref[6][0][2] = -1+1*step;

		local_ref[6][1][0] = -1+2*step;
		local_ref[6][1][1] = -1;
		local_ref[6][1][2] = -1+step;

		local_ref[6][2][0] = -1+2*step;
		local_ref[6][2][1] = -1+step;
		local_ref[6][2][2] = -1+2*step;

		local_ref[6][3][0] = -1+step;
		local_ref[6][3][1] = -1+step;
		local_ref[6][3][2] = -1+2*step;

		local_ref[6][4][0] = -1+step;
		local_ref[6][4][1] = -1;
		local_ref[6][4][2] = 1;

		local_ref[6][5][0] = -1+2*step;
		local_ref[6][5][1] = -1;
		local_ref[6][5][2] = 1;

		local_ref[6][6][0] = -1+2*step;
		local_ref[6][6][1] = -1+1*step;
		local_ref[6][6][2] = 1;

		local_ref[6][7][0] = -1+step;
		local_ref[6][7][1] = -1+1*step;
		local_ref[6][7][2] = 1;

		//element 7
		local_ref[7][0][0] = -1+2*step;
		local_ref[7][0][1] = -1;
		local_ref[7][0][2] = -1+step;

		local_ref[7][1][0] = 1;
		local_ref[7][1][1] = -1;
		local_ref[7][1][2] = -1;

		local_ref[7][2][0] =  1;
		local_ref[7][2][1] = -1+step;
		local_ref[7][2][2] = -1+1*step;

		local_ref[7][3][0] = -1+2*step;
		local_ref[7][3][1] = -1+1*step;
		local_ref[7][3][2] = -1+2*step;

		local_ref[7][4][0] = -1+2*step;
		local_ref[7][4][1] = -1;
		local_ref[7][4][2] = 1;

		local_ref[7][5][0] = 1;
		local_ref[7][5][1] = -1;
		local_ref[7][5][2] = 1;

		local_ref[7][6][0] = 1;
		local_ref[7][6][1] = -1+1*step;
		local_ref[7][6][2] = 1;

		local_ref[7][7][0] = -1+2*step;
		local_ref[7][7][1] = -1+1*step;
		local_ref[7][7][2] = 1;

		//element 8
		local_ref[8][0][0] = -1+2*step;
		local_ref[8][0][1] = -1+1*step;
		local_ref[8][0][2] = -1+2*step;

		local_ref[8][1][0] = 1;
		local_ref[8][1][1] = -1+1*step;
		local_ref[8][1][2] = -1+1*step;

		local_ref[8][2][0] =  1;
		local_ref[8][2][1] = -1+2*step;
		local_ref[8][2][2] = -1+1*step;

		local_ref[8][3][0] = -1+2*step;
		local_ref[8][3][1] = -1+2*step;
		local_ref[8][3][2] = -1+2*step;

		local_ref[8][4][0] = -1+2*step;
		local_ref[8][4][1] = -1+1*step;
		local_ref[8][4][2] = 1;

		local_ref[8][5][0] = 1;
		local_ref[8][5][1] = -1+1*step;
		local_ref[8][5][2] = 1;

		local_ref[8][6][0] = 1;
		local_ref[8][6][1] = -1+2*step;
		local_ref[8][6][2] = 1;

		local_ref[8][7][0] = -1+2*step;
		local_ref[8][7][1] = -1+2*step;
		local_ref[8][7][2] = 1;

		//element 9
		local_ref[9][0][0] = -1+2*step;
		local_ref[9][0][1] = -1+2*step;
		local_ref[9][0][2] = -1+2*step;

		local_ref[9][1][0] = 1;
		local_ref[9][1][1] = -1+2*step;
		local_ref[9][1][2] = -1+1*step;

		local_ref[9][2][0] =  1;
		local_ref[9][2][1] =  1;
		local_ref[9][2][2] = -1;

		local_ref[9][3][0] = -1+2*step;
		local_ref[9][3][1] =  1;
		local_ref[9][3][2] = -1+1*step;

		local_ref[9][4][0] = -1+2*step;
		local_ref[9][4][1] = -1+2*step;
		local_ref[9][4][2] = 1;

		local_ref[9][5][0] = 1;
		local_ref[9][5][1] = -1+2*step;
		local_ref[9][5][2] = 1;

		local_ref[9][6][0] = 1;
		local_ref[9][6][1] = 1;
		local_ref[9][6][2] = 1;

		local_ref[9][7][0] = -1+2*step;
		local_ref[9][7][1] = 1;
		local_ref[9][7][2] = 1;

		//element 10
		local_ref[10][0][0] = -1+1*step;
		local_ref[10][0][1] = -1+2*step;
		local_ref[10][0][2] = -1+2*step;

		local_ref[10][1][0] = -1+2*step;
		local_ref[10][1][1] = -1+2*step;
		local_ref[10][1][2] = -1+2*step;

		local_ref[10][2][0] =  -1+2*step;
		local_ref[10][2][1] =  1;
		local_ref[10][2][2] = -1+step;

		local_ref[10][3][0] = -1+1*step;
		local_ref[10][3][1] =  1;
		local_ref[10][3][2] = -1+1*step;

		local_ref[10][4][0] = -1+1*step;
		local_ref[10][4][1] = -1+2*step;
		local_ref[10][4][2] = 1;

		local_ref[10][5][0] = -1+2*step;
		local_ref[10][5][1] = -1+2*step;
		local_ref[10][5][2] = 1;

		local_ref[10][6][0] = -1+2*step;
		local_ref[10][6][1] = 1;
		local_ref[10][6][2] = 1;

		local_ref[10][7][0] = -1+1*step;
		local_ref[10][7][1] = 1;
		local_ref[10][7][2] = 1;

		//element 11
		local_ref[11][0][0] = -1;
		local_ref[11][0][1] = -1+2*step;
		local_ref[11][0][2] = -1+1*step;

		local_ref[11][1][0] = -1+1*step;
		local_ref[11][1][1] = -1+2*step;
		local_ref[11][1][2] = -1+2*step;

		local_ref[11][2][0] =  -1+1*step;
		local_ref[11][2][1] =  1;
		local_ref[11][2][2] = -1+step;

		local_ref[11][3][0] = -1;
		local_ref[11][3][1] =  1;
		local_ref[11][3][2] = -1;

		local_ref[11][4][0] = -1;
		local_ref[11][4][1] = -1+2*step;
		local_ref[11][4][2] = 1;

		local_ref[11][5][0] = -1+1*step;
		local_ref[11][5][1] = -1+2*step;
		local_ref[11][5][2] = 1;

		local_ref[11][6][0] = -1+1*step;
		local_ref[11][6][1] = 1;
		local_ref[11][6][2] = 1;

		local_ref[11][7][0] = -1;
		local_ref[11][7][1] = 1;
		local_ref[11][7][2] = 1;

		//element 12
		local_ref[12][0][0] = -1;
		local_ref[12][0][1] = -1+1*step;
		local_ref[12][0][2] = -1+1*step;

		local_ref[12][1][0] = -1+1*step;
		local_ref[12][1][1] = -1+1*step;
		local_ref[12][1][2] = -1+2*step;

		local_ref[12][2][0] =  -1+1*step;
		local_ref[12][2][1] =  -1+2*step;
		local_ref[12][2][2] = -1+2*step;

		local_ref[12][3][0] = -1;
		local_ref[12][3][1] =  -1+2*step;
		local_ref[12][3][2] = -1+step;

		local_ref[12][4][0] = -1;
		local_ref[12][4][1] = -1+1*step;
		local_ref[12][4][2] = 1;

		local_ref[12][5][0] = -1+1*step;
		local_ref[12][5][1] = -1+1*step;
		local_ref[12][5][2] = 1;

		local_ref[12][6][0] = -1+1*step;
		local_ref[12][6][1] = -1+2*step;
		local_ref[12][6][2] = 1;

		local_ref[12][7][0] = -1;
		local_ref[12][7][1] = -1+2*step;
		local_ref[12][7][2] = 1;

	}

	//define the rotation of the reference element
	int rot[3];
	int sym[3];
	int id_node[8];
	vector<int> ord;

	if(elem->pad==49 || elem->pad==55 || elem->pad==56 || elem->pad==57 || elem->pad==58
			|| elem->pad==79 || elem->pad==80 || elem->pad==81 || elem->pad==82){
		//edge 0 1 2 3
		rot[0] = 0;
		rot[1] = 0;
		rot[2] = 0;

		sym[0] = 0;
		sym[1] = 0;
		sym[2] = 1;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}
	}else if(elem->pad==50 || elem->pad==59 || elem->pad==60 || elem->pad==61 || elem->pad==62
			|| elem->pad==83 || elem->pad==84 || elem->pad==85 || elem->pad==86){
		//edge 8 9 10 11
		rot[0] = 0;
		rot[1] = 0;
		rot[2] = 0;

		sym[0] = 0;
		sym[1] = 0;
		sym[2] = 0;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}
	}else if(elem->pad==51 || elem->pad==63 || elem->pad==64 || elem->pad==65 || elem->pad==66
			|| elem->pad==87 || elem->pad==88 || elem->pad==89 || elem->pad==90){
		//edge 0 4 5 8
		rot[0] = -1;
		rot[1] = 0;
		rot[2] = 0;

		sym[0] = 0;
		sym[1] = 0;
		sym[2] = 0;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}
	}else if(elem->pad==52 || elem->pad==67 || elem->pad==68 || elem->pad==69 || elem->pad==70
			|| elem->pad==91 || elem->pad==92 || elem->pad==93 || elem->pad==94){
		//edge 1 5 6 9
		rot[0] = 0;
		rot[1] = -1;
		rot[2] = 0;

		sym[0] = 0;
		sym[1] = 0;
		sym[2] = 0;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}
	}else if(elem->pad==53 || elem->pad==71 || elem->pad==72 || elem->pad==73 || elem->pad==74
			|| elem->pad==95 || elem->pad==96 || elem->pad==97 || elem->pad==98){
		//edge 2 6 7 10
		rot[0] = 1;
		rot[1] = 0;
		rot[2] = 0;

		sym[0] = 0;
		sym[1] = 0;
		sym[2] = 0;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}
	}else if(elem->pad==54 || elem->pad==75 || elem->pad==76 || elem->pad==77 || elem->pad==78
			|| elem->pad==99 || elem->pad==100 || elem->pad==101 || elem->pad==102){
		//edge 3 4 7 11
		rot[0] = 0;
		rot[1] = 1;
		rot[2] = 0;

		sym[0] = 0;
		sym[1] = 0;
		sym[2] = 0;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}
	}else{
		printf ("Error in template 5\n");
		exit (EXIT_FAILURE);
	}



	for(int iel = 0; iel < 13; iel++){
		double local_ref_x[8];
		double local_ref_y[8];
		double local_ref_z[8];
		for(int ino = 0; ino < 8; ino++){
			local_ref_x[ino] = local_ref[iel][ino][0];
			local_ref_y[ino] = local_ref[iel][ino][1];
			local_ref_z[ino] = local_ref[iel][ino][2];
		}
		ApplyElement(mesh, coords,  id,  iel,  id_node,  local_ref_x, local_ref_y, local_ref_z,  ord,  hash_nodes);
	}
	sc_array_reset(&toto);
}

void ApplyTemplate6(hexa_tree_t* mesh, std::vector<double>& coords, std::vector<int>& elements_ids, int iel, sc_hash_array_t* hash_nodes){

	sc_array_t toto;
	sc_array_init(&toto, sizeof(octant_t));

	octant_t* elemOrig = (octant_t*) sc_array_index(&mesh->elements, elements_ids[iel]);
	octant_t * elem = (octant_t*) sc_array_push(&toto);
	hexa_element_copy(elemOrig,elem);

	double step = double(2)/double(3);
	int id = elements_ids[iel];

	double cord_in_ref[3];
	cord_in_ref[0] = 0;
	cord_in_ref[1] = 0;
	cord_in_ref[2] = 0;

	//reference element edge 0
	double local_ref[10][8][3];
	if(true){
		//element 0
		local_ref[0][0][0] = -1;
		local_ref[0][0][1] = -1;
		local_ref[0][0][2] = -1;

		local_ref[0][1][0] = -1+step;
		local_ref[0][1][1] = -1;
		local_ref[0][1][2] = -1;

		local_ref[0][2][0] = -1+2*3*step/double(4);
		local_ref[0][2][1] = -1+step;
		local_ref[0][2][2] = -1;

		local_ref[0][3][0] = -1+3*step/double(4);
		local_ref[0][3][1] = -1+step;
		local_ref[0][3][2] = -1;

		local_ref[0][4][0] = -1;
		local_ref[0][4][1] = -1;
		local_ref[0][4][2] = 1;

		local_ref[0][5][0] = -1+step;
		local_ref[0][5][1] = -1;
		local_ref[0][5][2] = 1;

		local_ref[0][6][0] = -1+2*3*step/double(4);
		local_ref[0][6][1] = -1+step;
		local_ref[0][6][2] = 1;

		local_ref[0][7][0] = -1+3*step/double(4);
		local_ref[0][7][1] = -1+step;
		local_ref[0][7][2] = 1;

		//element 1
		local_ref[1][0][0] = -1+1*step;
		local_ref[1][0][1] = -1;
		local_ref[1][0][2] = -1;

		local_ref[1][1][0] = -1+2*step;
		local_ref[1][1][1] = -1;
		local_ref[1][1][2] = -1;

		local_ref[1][2][0] = -1+3*3*step/double(4);
		local_ref[1][2][1] = -1+step;
		local_ref[1][2][2] = -1;

		local_ref[1][3][0] = -1+2*3*step/double(4);
		local_ref[1][3][1] = -1+step;
		local_ref[1][3][2] = -1;

		local_ref[1][4][0] = -1+1*step;
		local_ref[1][4][1] = -1;
		local_ref[1][4][2] = 1;

		local_ref[1][5][0] = -1+2*step;
		local_ref[1][5][1] = -1;
		local_ref[1][5][2] = 1;

		local_ref[1][6][0] = -1+3*3*step/double(4);
		local_ref[1][6][1] = -1+step;
		local_ref[1][6][2] = 1;

		local_ref[1][7][0] = -1+2*3*step/double(4);
		local_ref[1][7][1] = -1+step;
		local_ref[1][7][2] = 1;


		//element 2
		local_ref[2][0][0] = -1+2*step;
		local_ref[2][0][1] = -1;
		local_ref[2][0][2] = -1;

		local_ref[2][1][0] =  1;
		local_ref[2][1][1] = -1;
		local_ref[2][1][2] = -1;

		local_ref[2][2][0] =  1;
		local_ref[2][2][1] =  -1+step;
		local_ref[2][2][2] = -1;

		local_ref[2][3][0] = -1+3*3*step/double(4);
		local_ref[2][3][1] = -1+step;
		local_ref[2][3][2] = -1;

		local_ref[2][4][0] = -1+2*step;
		local_ref[2][4][1] = -1;
		local_ref[2][4][2] = 1;

		local_ref[2][5][0] =  1;
		local_ref[2][5][1] = -1;
		local_ref[2][5][2] =  1;

		local_ref[2][6][0] = 1;
		local_ref[2][6][1] = -1+step;
		local_ref[2][6][2] = 1;

		local_ref[2][7][0] = -1+3*3*step/double(4);
		local_ref[2][7][1] = -1+step;
		local_ref[2][7][2] = 1;

		//element 3
		local_ref[3][0][0] = -1+1*3*step/double(4);
		local_ref[3][0][1] = -1+step;
		local_ref[3][0][2] = -1;

		local_ref[3][1][0] = -1+2*3*step/double(4);
		local_ref[3][1][1] = -1+step;
		local_ref[3][1][2] = -1;

		local_ref[3][2][0] = -1+2*3*step/double(4);
		local_ref[3][2][1] = -1+2*step;
		local_ref[3][2][2] = -1;

		local_ref[3][3][0] = -1+1*3*step/double(4);
		local_ref[3][3][1] = -1+2*step;
		local_ref[3][3][2] = -1;

		local_ref[3][4][0] = -1+1*3*step/double(4);
		local_ref[3][4][1] = -1+step;
		local_ref[3][4][2] = 1;

		local_ref[3][5][0] = -1+2*3*step/double(4);
		local_ref[3][5][1] = -1+step;
		local_ref[3][5][2] = 1;

		local_ref[3][6][0] = -1+2*3*step/double(4);
		local_ref[3][6][1] = -1+2*step;
		local_ref[3][6][2] = 1;

		local_ref[3][7][0] = -1+1*3*step/double(4);
		local_ref[3][7][1] = -1+2*step;
		local_ref[3][7][2] = 1;

		//element 4
		local_ref[4][0][0] = -1+2*3*step/double(4);
		local_ref[4][0][1] = -1+step;
		local_ref[4][0][2] = -1;

		local_ref[4][1][0] = -1+3*3*step/double(4);
		local_ref[4][1][1] = -1+step;
		local_ref[4][1][2] = -1;

		local_ref[4][2][0] = -1+3*3*step/double(4);
		local_ref[4][2][1] = -1+2*step;
		local_ref[4][2][2] = -1;

		local_ref[4][3][0] = -1+2*3*step/double(4);
		local_ref[4][3][1] = -1+2*step;
		local_ref[4][3][2] = -1;

		local_ref[4][4][0] = -1+2*3*step/double(4);
		local_ref[4][4][1] = -1+step;
		local_ref[4][4][2] = 1;

		local_ref[4][5][0] = -1+3*3*step/double(4);
		local_ref[4][5][1] = -1+step;
		local_ref[4][5][2] = 1;

		local_ref[4][6][0] = -1+3*3*step/double(4);
		local_ref[4][6][1] = -1+2*step;
		local_ref[4][6][2] = 1;

		local_ref[4][7][0] = -1+2*3*step/double(4);
		local_ref[4][7][1] = -1+2*step;
		local_ref[4][7][2] = 1;

		//element 5
		local_ref[5][0][0] = -1+3*3*step/double(4);
		local_ref[5][0][1] = -1+step;
		local_ref[5][0][2] = -1;

		local_ref[5][1][0] = 1;
		local_ref[5][1][1] = -1+step;
		local_ref[5][1][2] = -1;

		local_ref[5][2][0] = 1;
		local_ref[5][2][1] = -1+2*step;
		local_ref[5][2][2] = -1;

		local_ref[5][3][0] = -1+3*3*step/double(4);
		local_ref[5][3][1] = -1+2*step;
		local_ref[5][3][2] = -1;

		local_ref[5][4][0] = -1+3*3*step/double(4);;
		local_ref[5][4][1] = -1+step;
		local_ref[5][4][2] = 1;

		local_ref[5][5][0] = 1;
		local_ref[5][5][1] = -1+step;
		local_ref[5][5][2] = 1;

		local_ref[5][6][0] = 1;
		local_ref[5][6][1] = -1+2*step;
		local_ref[5][6][2] = 1;

		local_ref[5][7][0] = -1+3*3*step/double(4);
		local_ref[5][7][1] = -1+2*step;
		local_ref[5][7][2] = 1;

		//element 6
		local_ref[6][0][0] = -1+1*3*step/double(4);
		local_ref[6][0][1] = -1+2*step;
		local_ref[6][0][2] = -1;

		local_ref[6][1][0] = -1+2*3*step/double(4);
		local_ref[6][1][1] = -1+2*step;
		local_ref[6][1][2] = -1;

		local_ref[6][2][0] = -1+step;
		local_ref[6][2][1] = 1;
		local_ref[6][2][2] = -1;

		local_ref[6][3][0] = -1;
		local_ref[6][3][1] = 1;
		local_ref[6][3][2] = -1;

		local_ref[6][4][0] = -1+1*3*step/double(4);
		local_ref[6][4][1] = -1+2*step;
		local_ref[6][4][2] = 1;

		local_ref[6][5][0] = -1+2*3*step/double(4);
		local_ref[6][5][1] = -1+2*step;
		local_ref[6][5][2] = 1;

		local_ref[6][6][0] = -1+step;
		local_ref[6][6][1] = 1;
		local_ref[6][6][2] = 1;

		local_ref[6][7][0] = -1;
		local_ref[6][7][1] = 1;
		local_ref[6][7][2] = 1;

		//element 7
		local_ref[7][0][0] = -1+2*3*step/double(4);
		local_ref[7][0][1] = -1+2*step;
		local_ref[7][0][2] = -1;

		local_ref[7][1][0] = -1+3*3*step/double(4);
		local_ref[7][1][1] = -1+2*step;
		local_ref[7][1][2] = -1;

		local_ref[7][2][0] = -1+2*step;
		local_ref[7][2][1] = 1;
		local_ref[7][2][2] = -1;

		local_ref[7][3][0] = -1+step;
		local_ref[7][3][1] = 1;
		local_ref[7][3][2] = -1;

		local_ref[7][4][0] = -1+2*3*step/double(4);
		local_ref[7][4][1] = -1+2*step;
		local_ref[7][4][2] = 1;

		local_ref[7][5][0] = -1+3*3*step/double(4);
		local_ref[7][5][1] = -1+2*step;
		local_ref[7][5][2] = 1;

		local_ref[7][6][0] = -1+2*step;
		local_ref[7][6][1] = 1;
		local_ref[7][6][2] = 1;

		local_ref[7][7][0] = -1+step;
		local_ref[7][7][1] = 1;
		local_ref[7][7][2] = 1;


		//element 8
		local_ref[8][0][0] = -1+3*3*step/double(4);
		local_ref[8][0][1] = -1+2*step;
		local_ref[8][0][2] = -1;

		local_ref[8][1][0] = 1;
		local_ref[8][1][1] = -1+2*step;
		local_ref[8][1][2] = -1;

		local_ref[8][2][0] = 1;
		local_ref[8][2][1] = 1;
		local_ref[8][2][2] = -1;

		local_ref[8][3][0] = -1+2*step;
		local_ref[8][3][1] = 1;
		local_ref[8][3][2] = -1;

		local_ref[8][4][0] = -1+3*3*step/double(4);
		local_ref[8][4][1] = -1+2*step;
		local_ref[8][4][2] = 1;

		local_ref[8][5][0] = 1;
		local_ref[8][5][1] = -1+2*step;
		local_ref[8][5][2] = 1;

		local_ref[8][6][0] = 1;
		local_ref[8][6][1] = 1;
		local_ref[8][6][2] = 1;

		local_ref[8][7][0] = -1+2*step;
		local_ref[8][7][1] = 1;
		local_ref[8][7][2] = 1;

		//element 9
		local_ref[9][0][0] = -1;
		local_ref[9][0][1] = -1;
		local_ref[9][0][2] = -1;

		local_ref[9][1][0] = -1+1*3*step/double(4);
		local_ref[9][1][1] = -1+step;
		local_ref[9][1][2] = -1;

		local_ref[9][2][0] = -1+1*3*step/double(4);
		local_ref[9][2][1] = -1+2*step;
		local_ref[9][2][2] = -1;

		local_ref[9][3][0] = -1;
		local_ref[9][3][1] = 1;
		local_ref[9][3][2] = -1;

		local_ref[9][4][0] = -1;
		local_ref[9][4][1] = -1;
		local_ref[9][4][2] = 1;

		local_ref[9][5][0] = -1+1*3*step/double(4);
		local_ref[9][5][1] = -1+step;
		local_ref[9][5][2] = 1;

		local_ref[9][6][0] = -1+1*3*step/double(4);
		local_ref[9][6][1] = -1+2*step;
		local_ref[9][6][2] = 1;

		local_ref[9][7][0] = -1;
		local_ref[9][7][1] = 1;
		local_ref[9][7][2] = 1;
	}

	//define the rotation of the reference element
	int rot[3];
	int sym[3];
	int id_node[8];
	vector<int> ord;

	if(elem->pad==103){
		//edge 0 1 3 8 9 11
		rot[0] = 0;
		rot[1] = 0;
		rot[2] = 1;

		sym[0] = 0;
		sym[1] = 0;
		sym[2] = 0;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}

	}else if(elem->pad==104){
		//edge 0 1 2 8 9 10
		rot[0] = 0;
		rot[1] = 0;
		rot[2] = 0;

		sym[0] = 0;
		sym[1] = 0;
		sym[2] = 0;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}

	}else if(elem->pad==105){
		//edge 1 2 3 9 10 11
		rot[0] = 0;
		rot[1] = 0;
		rot[2] = -1;

		sym[0] = 0;
		sym[1] = 0;
		sym[2] = 0;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}

	}else if(elem->pad==106){
		//edge 0 2 3 8 10 11
		rot[0] = 0;
		rot[1] = 0;
		rot[2] = 0;

		sym[0] = 1;
		sym[1] = 0;
		sym[2] = 0;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}

	}else if(elem->pad==107){
		//edge 0 5 8 2 6 10
		rot[0] = 1;
		rot[1] = 0;
		rot[2] = 0;

		sym[0] = 0;
		sym[1] = 0;
		sym[2] = 0;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}

	}else if(elem->pad==108){
		//edge 5 8 4 6 10 7
		rot[0] = 1;
		rot[1] = 0;
		rot[2] = 1;

		sym[0] = 0;
		sym[1] = 0;
		sym[2] = 0;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}

	}else if(elem->pad==109){
		//edge 0 4 8 2 7 10
		rot[0] = 1;
		rot[1] = 0;
		rot[2] = 0;

		sym[0] = 1;
		sym[1] = 0;
		sym[2] = 0;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}

	}else if(elem->pad==110){
		//edge 5 0 4 6 2 7
		rot[0] = 1;
		rot[1] = 0;
		rot[2] = -1;

		sym[0] = 0;
		sym[1] = 0;
		sym[2] = 0;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}

	}else if(elem->pad==111){
		//edge 4 3 7 5 1 6
		rot[0] = 0;
		rot[1] = 1;
		rot[2] = 0;

		sym[0] = 1;
		sym[1] = 0;
		sym[2] = 0;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}

	}else if(elem->pad==112){
		//edge 1 5 9 3 4 11
		rot[0] = 1;
		rot[1] = 1;
		rot[2] = 0;

		sym[0] = 1;
		sym[1] = 0;
		sym[2] = 0;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}

	}else if(elem->pad==113){
		//edge 5 9 6 4 11 7
		rot[0] = 0;
		rot[1] = 1;
		rot[2] = 0;

		sym[0] = 0;
		sym[1] = 0;
		sym[2] = 0;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}

	}else if(elem->pad==114){
		//edge 1 6 9 3 7 11
		rot[0] = 1;
		rot[1] = 1;
		rot[2] = 0;

		sym[0] = 0;
		sym[1] = 0;
		sym[2] = 0;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}

	}else{
		printf ("Error in template 6\n");
		exit (EXIT_FAILURE);
	}


	for(int iel = 0; iel < 10; iel++){
		double local_ref_x[8];
		double local_ref_y[8];
		double local_ref_z[8];
		for(int ino = 0; ino < 8; ino++){
			local_ref_x[ino] = local_ref[iel][ino][0];
			local_ref_y[ino] = local_ref[iel][ino][1];
			local_ref_z[ino] = local_ref[iel][ino][2];
		}
		ApplyElement(mesh, coords,  id,  iel,  id_node,  local_ref_x, local_ref_y, local_ref_z,  ord,  hash_nodes);
	}
	sc_array_reset(&toto);
}

void ApplyTemplate7(hexa_tree_t* mesh, std::vector<double>& coords, std::vector<int>& elements_ids, int iel, sc_hash_array_t* hash_nodes){

	sc_array_t toto;
	sc_array_init(&toto, sizeof(octant_t));

	octant_t* elemOrig = (octant_t*) sc_array_index(&mesh->elements, elements_ids[iel]);
	octant_t * elem = (octant_t*) sc_array_push(&toto);
	hexa_element_copy(elemOrig,elem);

	double step = double(2)/double(3);
	int id = elements_ids[iel];

	double cord_in_ref[3];
	cord_in_ref[0] = 0;
	cord_in_ref[1] = 0;
	cord_in_ref[2] = 0;

	//reference element edge 0
	double local_ref[7][8][3];
	if(true){
		//element 0
		local_ref[0][0][0] = -1;
		local_ref[0][0][1] = -1;
		local_ref[0][0][2] = -1;

		local_ref[0][1][0] = -1+step;
		local_ref[0][1][1] = -1;
		local_ref[0][1][2] = -1;

		local_ref[0][2][0] = -1+step;
		local_ref[0][2][1] = -1+step;
		local_ref[0][2][2] = -1;

		local_ref[0][3][0] = -1;
		local_ref[0][3][1] = -1+step;
		local_ref[0][3][2] = -1;

		local_ref[0][4][0] = -1;
		local_ref[0][4][1] = -1;
		local_ref[0][4][2] =  -1+step;

		local_ref[0][5][0] = -1+step;
		local_ref[0][5][1] = -1;
		local_ref[0][5][2] = -1+step;

		local_ref[0][6][0] = -1+step;
		local_ref[0][6][1] = -1+step;
		local_ref[0][6][2] = -1+step;

		local_ref[0][7][0] = -1;
		local_ref[0][7][1] = -1+step;
		local_ref[0][7][2] =  -1+step;

		//element 1
		local_ref[1][0][0] = -1+step;
		local_ref[1][0][1] = -1;
		local_ref[1][0][2] = -1;

		local_ref[1][1][0] = -1+2*step;
		local_ref[1][1][1] = -1;
		local_ref[1][1][2] = -1;

		local_ref[1][2][0] = -1+2*step;
		local_ref[1][2][1] = -1+2*step;
		local_ref[1][2][2] = -1;

		local_ref[1][3][0] = -1+step;
		local_ref[1][3][1] = -1+step;
		local_ref[1][3][2] = -1;

		local_ref[1][4][0] = -1+step;
		local_ref[1][4][1] = -1;
		local_ref[1][4][2] = -1+step;

		local_ref[1][5][0] = -1+2*step;
		local_ref[1][5][1] = -1;
		local_ref[1][5][2] = -1+2*step;

		local_ref[1][6][0] = -1+2*step;
		local_ref[1][6][1] = -1+2*step;
		local_ref[1][6][2] = -1+2*step;;

		local_ref[1][7][0] = -1+step;
		local_ref[1][7][1] = -1+step;
		local_ref[1][7][2] = -1+step;


		//element 2
		local_ref[2][0][0] = -1+2*step;
		local_ref[2][0][1] = -1;
		local_ref[2][0][2] = -1;

		local_ref[2][1][0] = 1;
		local_ref[2][1][1] = -1;
		local_ref[2][1][2] = -1;

		local_ref[2][2][0] = 1;
		local_ref[2][2][1] = 1;
		local_ref[2][2][2] = -1;

		local_ref[2][3][0] = -1+2*step;
		local_ref[2][3][1] = -1+2*step;
		local_ref[2][3][2] = -1;

		local_ref[2][4][0] = -1+2*step;
		local_ref[2][4][1] = -1;
		local_ref[2][4][2] = -1+2*step;

		local_ref[2][5][0] =  1;
		local_ref[2][5][1] =  -1;
		local_ref[2][5][2] =  1;

		local_ref[2][6][0] = 1;
		local_ref[2][6][1] = 1;
		local_ref[2][6][2] = 1;

		local_ref[2][7][0] = -1+2*step;
		local_ref[2][7][1] = -1+2*step;
		local_ref[2][7][2] = -1+2*step;


		//element 3
		local_ref[3][0][0] = -1;
		local_ref[3][0][1] = -1+step;
		local_ref[3][0][2] = -1;

		local_ref[3][1][0] = -1+step;
		local_ref[3][1][1] = -1+step;
		local_ref[3][1][2] = -1;

		local_ref[3][2][0] = -1+2*step;
		local_ref[3][2][1] = -1+2*step;
		local_ref[3][2][2] = -1;

		local_ref[3][3][0] = -1;
		local_ref[3][3][1] = -1+2*step;
		local_ref[3][3][2] = -1;

		local_ref[3][4][0] = -1;
		local_ref[3][4][1] = -1+step;
		local_ref[3][4][2] = -1+step;

		local_ref[3][5][0] = -1+step;
		local_ref[3][5][1] = -1+step;
		local_ref[3][5][2] = -1+step;

		local_ref[3][6][0] = -1+2*step;
		local_ref[3][6][1] = -1+2*step;
		local_ref[3][6][2] = -1+2*step;

		local_ref[3][7][0] = -1;
		local_ref[3][7][1] = -1+2*step;
		local_ref[3][7][2] = -1+2*step;

		//element 4
		local_ref[4][0][0] = -1;
		local_ref[4][0][1] = -1;
		local_ref[4][0][2] = -1+1*step;

		local_ref[4][1][0] = -1+1*step;
		local_ref[4][1][1] = -1;
		local_ref[4][1][2] = -1+1*step;

		local_ref[4][2][0] = -1+1*step;
		local_ref[4][2][1] = -1+1*step;
		local_ref[4][2][2] = -1+1*step;

		local_ref[4][3][0] = -1;
		local_ref[4][3][1] = -1+1*step;
		local_ref[4][3][2] = -1+1*step;

		local_ref[4][4][0] = -1;
		local_ref[4][4][1] = -1;
		local_ref[4][4][2] = -1+2*step;

		local_ref[4][5][0] = -1+2*step;
		local_ref[4][5][1] = -1+0*step;
		local_ref[4][5][2] = -1+2*step;

		local_ref[4][6][0] = -1+2*step;
		local_ref[4][6][1] = -1+2*step;
		local_ref[4][6][2] = -1+2*step;

		local_ref[4][7][0] = -1+0*step;
		local_ref[4][7][1] = -1+2*step;
		local_ref[4][7][2] = -1+2*step;

		//element 5
		local_ref[5][0][0] = -1;
		local_ref[5][0][1] = -1;
		local_ref[5][0][2] = -1+2*step;

		local_ref[5][1][0] =  -1+2*step;
		local_ref[5][1][1] = -1;
		local_ref[5][1][2] = -1+2*step;

		local_ref[5][2][0] = -1+2*step;
		local_ref[5][2][1] = -1+2*step;
		local_ref[5][2][2] = -1+2*step;

		local_ref[5][3][0] = -1;
		local_ref[5][3][1] = -1+2*step;
		local_ref[5][3][2] = -1+2*step;

		local_ref[5][4][0] = -1;
		local_ref[5][4][1] = -1;
		local_ref[5][4][2] =  1;

		local_ref[5][5][0] = 1;
		local_ref[5][5][1] = -1;
		local_ref[5][5][2] = 1;

		local_ref[5][6][0] = 1;
		local_ref[5][6][1] = 1;
		local_ref[5][6][2] = 1;

		local_ref[5][7][0] = -1;
		local_ref[5][7][1] = 1;
		local_ref[5][7][2] = 1;

		//element 6
		local_ref[6][0][0] = -1;
		local_ref[6][0][1] = -1+2*step;
		local_ref[6][0][2] = -1;

		local_ref[6][1][0] = -1+2*step;
		local_ref[6][1][1] = -1+2*step;
		local_ref[6][1][2] = -1;

		local_ref[6][2][0] = 1;
		local_ref[6][2][1] = 1;
		local_ref[6][2][2] = -1;

		local_ref[6][3][0] = -1;
		local_ref[6][3][1] = 1;
		local_ref[6][3][2] = -1;

		local_ref[6][4][0] = -1;
		local_ref[6][4][1] = -1+2*step;
		local_ref[6][4][2] =  -1+2*step;;

		local_ref[6][5][0] = -1+2*step;;
		local_ref[6][5][1] = -1+2*step;
		local_ref[6][5][2] = -1+2*step;;

		local_ref[6][6][0] = 1;
		local_ref[6][6][1] = 1;
		local_ref[6][6][2] = 1;

		local_ref[6][7][0] = -1;
		local_ref[6][7][1] = 1;
		local_ref[6][7][2] = 1;
	}

	//define the rotation of the reference element
	int rot[3];
	int sym[3];
	int id_node[8];
	vector<int> ord;

	if(elem->pad==115){
		//edge 0 4 3
		rot[0] = 0;
		rot[1] = 0;
		rot[2] = 0;

		sym[0] = 0;
		sym[1] = 0;
		sym[2] = 0;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}

	}else if(elem->pad==116){
		//edge 0 1 5
		rot[0] = 0;
		rot[1] = 0;
		rot[2] = -1;

		sym[0] = 0;
		sym[1] = 0;
		sym[2] = 0;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}

	}else if(elem->pad==117){
		//edge 1 2 6
		rot[0] = 0;
		rot[1] = 0;
		rot[2] = -1;

		sym[0] = 1;
		sym[1] = 0;
		sym[2] = 0;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}

	}else if(elem->pad==118){
		//edge 2 3 7
		rot[0] = 0;
		rot[1] = 0;
		rot[2] = 1;

		sym[0] = 0;
		sym[1] = 0;
		sym[2] = 0;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}

	}else if(elem->pad==119){
		//edge 4 8 11
		rot[0] = 0;
		rot[1] = 0;
		rot[2] = 0;

		sym[0] = 0;
		sym[1] = 0;
		sym[2] = 1;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}

	}else if(elem->pad==120){
		//edge 8 9 5
		rot[0] = 0;
		rot[1] = 0;
		rot[2] = -1;

		sym[0] = 0;
		sym[1] = 0;
		sym[2] = 1;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}

	}else if(elem->pad==121){
		//edge 9 10 6
		rot[0] = 0;
		rot[1] = 0;
		rot[2] = -1;

		sym[0] = 1;
		sym[1] = 0;
		sym[2] = 1;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}

	}else if(elem->pad==122){
		//edge 11 10 7
		rot[0] = 0;
		rot[1] = 0;
		rot[2] = 1;

		sym[0] = 0;
		sym[1] = 0;
		sym[2] = 1;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}

	}else{
		printf ("Error in template 7\n");
		exit (EXIT_FAILURE);
	}


	for(int iel = 0; iel < 7; iel++){
		double local_ref_x[8];
		double local_ref_y[8];
		double local_ref_z[8];
		for(int ino = 0; ino < 8; ino++){
			local_ref_x[ino] = local_ref[iel][ino][0];
			local_ref_y[ino] = local_ref[iel][ino][1];
			local_ref_z[ino] = local_ref[iel][ino][2];
		}
		ApplyElement(mesh, coords,  id,  iel,  id_node,  local_ref_x, local_ref_y, local_ref_z,  ord,  hash_nodes);
	}
	sc_array_reset(&toto);
}

void ApplyTemplate8(hexa_tree_t* mesh, std::vector<double>& coords, std::vector<int>& elements_ids, int iel, sc_hash_array_t* hash_nodes){

	sc_array_t toto;
	sc_array_init(&toto, sizeof(octant_t));

	octant_t* elemOrig = (octant_t*) sc_array_index(&mesh->elements, elements_ids[iel]);
	octant_t * elem = (octant_t*) sc_array_push(&toto);
	hexa_element_copy(elemOrig,elem);

	double step = double(2)/double(3);
	int id = elements_ids[iel];

	double cord_in_ref[3];
	cord_in_ref[0] = 0;
	cord_in_ref[1] = 0;
	cord_in_ref[2] = 0;

	//reference element edge 0
	double local_ref[9][8][3];
	if(true){
		//element 0
		local_ref[0][0][0] = -1;
		local_ref[0][0][1] = -1;
		local_ref[0][0][2] = -1;

		local_ref[0][1][0] = -1+step;
		local_ref[0][1][1] = -1;
		local_ref[0][1][2] = -1;

		local_ref[0][2][0] = -1+step;
		local_ref[0][2][1] = -1+step;
		local_ref[0][2][2] = -1;

		local_ref[0][3][0] = -1;
		local_ref[0][3][1] = -1+step;
		local_ref[0][3][2] = -1;

		local_ref[0][4][0] = -1;
		local_ref[0][4][1] = -1;
		local_ref[0][4][2] = 1;

		local_ref[0][5][0] = -1+step;
		local_ref[0][5][1] = -1;
		local_ref[0][5][2] = 1;

		local_ref[0][6][0] = -1+step;
		local_ref[0][6][1] = -1+step;
		local_ref[0][6][2] = 1;

		local_ref[0][7][0] = -1;
		local_ref[0][7][1] = -1+step;
		local_ref[0][7][2] = 1;

		//element 1
		local_ref[1][0][0] = -1+step;
		local_ref[1][0][1] = -1;
		local_ref[1][0][2] = -1;

		local_ref[1][1][0] = -1+2*step;
		local_ref[1][1][1] = -1;
		local_ref[1][1][2] = -1;

		local_ref[1][2][0] = -1+2*step;
		local_ref[1][2][1] = -1+step;
		local_ref[1][2][2] = -1;

		local_ref[1][3][0] = -1+step;
		local_ref[1][3][1] = -1+step;
		local_ref[1][3][2] = -1;

		local_ref[1][4][0] = -1+step;
		local_ref[1][4][1] = -1;
		local_ref[1][4][2] = 1;

		local_ref[1][5][0] = -1+2*step;
		local_ref[1][5][1] = -1;
		local_ref[1][5][2] = 1;

		local_ref[1][6][0] = -1+2*step;
		local_ref[1][6][1] = -1+step;
		local_ref[1][6][2] = 1;

		local_ref[1][7][0] = -1+step;
		local_ref[1][7][1] = -1+step;
		local_ref[1][7][2] = 1;

		//element 2
		local_ref[2][0][0] = -1+2*step;
		local_ref[2][0][1] = -1;
		local_ref[2][0][2] = -1;

		local_ref[2][1][0] = 1;
		local_ref[2][1][1] = -1;
		local_ref[2][1][2] = -1;

		local_ref[2][2][0] = 1;
		local_ref[2][2][1] = -1+step;
		local_ref[2][2][2] = -1;

		local_ref[2][3][0] = -1+2*step;
		local_ref[2][3][1] = -1+step;
		local_ref[2][3][2] = -1;

		local_ref[2][4][0] = -1+2*step;
		local_ref[2][4][1] = -1;
		local_ref[2][4][2] = 1;

		local_ref[2][5][0] = 1;
		local_ref[2][5][1] = -1;
		local_ref[2][5][2] = 1;

		local_ref[2][6][0] = 1;
		local_ref[2][6][1] = -1+step;
		local_ref[2][6][2] = 1;

		local_ref[2][7][0] = -1+2*step;
		local_ref[2][7][1] = -1+step;
		local_ref[2][7][2] = 1;

		//element 3
		local_ref[3][0][0] = -1;
		local_ref[3][0][1] = -1+step;
		local_ref[3][0][2] = -1;

		local_ref[3][1][0] = -1+step;
		local_ref[3][1][1] = -1+step;
		local_ref[3][1][2] = -1;

		local_ref[3][2][0] = -1+step;
		local_ref[3][2][1] = -1+2*step;
		local_ref[3][2][2] = -1;

		local_ref[3][3][0] = -1;
		local_ref[3][3][1] = -1+2*step;
		local_ref[3][3][2] = -1;

		local_ref[3][4][0] = -1;
		local_ref[3][4][1] = -1+step;
		local_ref[3][4][2] = 1;

		local_ref[3][5][0] = -1+step;
		local_ref[3][5][1] = -1+step;
		local_ref[3][5][2] = 1;

		local_ref[3][6][0] = -1+step;
		local_ref[3][6][1] = -1+2*step;
		local_ref[3][6][2] = 1;

		local_ref[3][7][0] = -1;
		local_ref[3][7][1] = -1+2*step;
		local_ref[3][7][2] = 1;

		//element 4
		local_ref[4][0][0] = -1+step;
		local_ref[4][0][1] = -1+step;
		local_ref[4][0][2] = -1;

		local_ref[4][1][0] = -1+2*step;
		local_ref[4][1][1] = -1+step;
		local_ref[4][1][2] = -1;

		local_ref[4][2][0] = -1+2*step;
		local_ref[4][2][1] = -1+2*step;
		local_ref[4][2][2] = -1;

		local_ref[4][3][0] = -1+step;
		local_ref[4][3][1] = -1+2*step;
		local_ref[4][3][2] = -1;

		local_ref[4][4][0] = -1+step;
		local_ref[4][4][1] = -1+step;
		local_ref[4][4][2] = 1;

		local_ref[4][5][0] = -1+2*step;
		local_ref[4][5][1] = -1+step;
		local_ref[4][5][2] = 1;

		local_ref[4][6][0] = -1+2*step;
		local_ref[4][6][1] = -1+2*step;
		local_ref[4][6][2] = 1;

		local_ref[4][7][0] = -1+step;
		local_ref[4][7][1] = -1+2*step;
		local_ref[4][7][2] = 1;

		//element 5
		local_ref[5][0][0] = -1+2*step;
		local_ref[5][0][1] = -1+step;
		local_ref[5][0][2] = -1;

		local_ref[5][1][0] = 1;
		local_ref[5][1][1] = -1+step;
		local_ref[5][1][2] = -1;

		local_ref[5][2][0] = 1;
		local_ref[5][2][1] = -1+2*step;
		local_ref[5][2][2] = -1;

		local_ref[5][3][0] = -1+2*step;
		local_ref[5][3][1] = -1+2*step;
		local_ref[5][3][2] = -1;

		local_ref[5][4][0] = -1+2*step;
		local_ref[5][4][1] = -1+step;
		local_ref[5][4][2] = 1;

		local_ref[5][5][0] = 1;
		local_ref[5][5][1] = -1+step;
		local_ref[5][5][2] = 1;

		local_ref[5][6][0] = 1;
		local_ref[5][6][1] = -1+2*step;
		local_ref[5][6][2] = 1;

		local_ref[5][7][0] = -1+2*step;
		local_ref[5][7][1] = -1+2*step;
		local_ref[5][7][2] = 1;

		//element 6
		local_ref[6][0][0] = -1;
		local_ref[6][0][1] = -1+2*step;
		local_ref[6][0][2] = -1;

		local_ref[6][1][0] = -1+step;
		local_ref[6][1][1] = -1+2*step;
		local_ref[6][1][2] = -1;

		local_ref[6][2][0] = -1+step;
		local_ref[6][2][1] = 1;
		local_ref[6][2][2] = -1;

		local_ref[6][3][0] = -1;
		local_ref[6][3][1] = 1;
		local_ref[6][3][2] = -1;

		local_ref[6][4][0] = -1;
		local_ref[6][4][1] = -1+2*step;
		local_ref[6][4][2] = 1;

		local_ref[6][5][0] = -1+step;
		local_ref[6][5][1] = -1+2*step;
		local_ref[6][5][2] = 1;

		local_ref[6][6][0] = -1+step;
		local_ref[6][6][1] = 1;
		local_ref[6][6][2] = 1;

		local_ref[6][7][0] = -1;
		local_ref[6][7][1] = 1;
		local_ref[6][7][2] = 1;

		//element 7
		local_ref[7][0][0] = -1+step;
		local_ref[7][0][1] = -1+2*step;
		local_ref[7][0][2] = -1;

		local_ref[7][1][0] = -1+2*step;
		local_ref[7][1][1] = -1+2*step;
		local_ref[7][1][2] = -1;

		local_ref[7][2][0] = -1+2*step;
		local_ref[7][2][1] = 1;
		local_ref[7][2][2] = -1;

		local_ref[7][3][0] = -1+step;
		local_ref[7][3][1] = 1;
		local_ref[7][3][2] = -1;

		local_ref[7][4][0] = -1+step;
		local_ref[7][4][1] = -1+2*step;
		local_ref[7][4][2] = 1;

		local_ref[7][5][0] = -1+2*step;
		local_ref[7][5][1] = -1+2*step;
		local_ref[7][5][2] = 1;

		local_ref[7][6][0] = -1+2*step;
		local_ref[7][6][1] = 1;
		local_ref[7][6][2] = 1;

		local_ref[7][7][0] = -1+step;
		local_ref[7][7][1] = 1;
		local_ref[7][7][2] = 1;

		//element 8
		local_ref[8][0][0] = -1+2*step;
		local_ref[8][0][1] = -1+2*step;
		local_ref[8][0][2] = -1;

		local_ref[8][1][0] = 1;
		local_ref[8][1][1] = -1+2*step;
		local_ref[8][1][2] = -1;

		local_ref[8][2][0] = 1;
		local_ref[8][2][1] = 1;
		local_ref[8][2][2] = -1;

		local_ref[8][3][0] = -1+2*step;
		local_ref[8][3][1] = 1;
		local_ref[8][3][2] = -1;

		local_ref[8][4][0] = -1+2*step;
		local_ref[8][4][1] = -1+2*step;
		local_ref[8][4][2] = 1;

		local_ref[8][5][0] = 1;
		local_ref[8][5][1] = -1+2*step;
		local_ref[8][5][2] = 1;

		local_ref[8][6][0] = 1;
		local_ref[8][6][1] = 1;
		local_ref[8][6][2] = 1;

		local_ref[8][7][0] = -1+2*step;
		local_ref[8][7][1] = 1;
		local_ref[8][7][2] = 1;
	}

	//define the rotation of the reference element
	int rot[3];
	int sym[3];
	int id_node[8];
	vector<int> ord;

	if(elem->pad==123){
		//edge 0 1 2 3 8 9 10 11
		rot[0] = 0;
		rot[1] = 0;
		rot[2] = 0;

		sym[0] = 0;
		sym[1] = 0;
		sym[2] = 0;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}

	}else if(elem->pad==124){
		//edge 0 5 4 8 2 6 10 7
		rot[0] = 1;
		rot[1] = 0;
		rot[2] = 0;

		sym[0] = 0;
		sym[1] = 0;
		sym[2] = 0;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}

	}else if(elem->pad==125){
		//edge 1 6 8 5 3 7 11 4
		rot[0] = 0;
		rot[1] = 1;
		rot[2] = 0;

		sym[0] = 0;
		sym[1] = 0;
		sym[2] = 0;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}

	}else{
		printf ("Error in template 8\n");
		exit (EXIT_FAILURE);
	}


	for(int iel = 0; iel < 9; iel++){
		double local_ref_x[8];
		double local_ref_y[8];
		double local_ref_z[8];
		for(int ino = 0; ino < 8; ino++){
			local_ref_x[ino] = local_ref[iel][ino][0];
			local_ref_y[ino] = local_ref[iel][ino][1];
			local_ref_z[ino] = local_ref[iel][ino][2];
		}
		ApplyElement(mesh, coords,  id,  iel,  id_node,  local_ref_x, local_ref_y, local_ref_z,  ord,  hash_nodes);
	}
	sc_array_reset(&toto);
}

void ApplyTemplate9(hexa_tree_t* mesh, std::vector<double>& coords, std::vector<int>& elements_ids, int iel, sc_hash_array_t* hash_nodes){

	sc_array_t toto;
	sc_array_init(&toto, sizeof(octant_t));

	octant_t* elemOrig = (octant_t*) sc_array_index(&mesh->elements, elements_ids[iel]);
	octant_t * elem = (octant_t*) sc_array_push(&toto);
	hexa_element_copy(elemOrig,elem);

	double step = double(2)/double(3);
	int id = elements_ids[iel];

	double cord_in_ref[3];
	cord_in_ref[0] = 0;
	cord_in_ref[1] = 0;
	cord_in_ref[2] = 0;

	//reference element edge 0
	double local_ref[31][8][3];
	if(true){
		//element 0
		local_ref[0][0][0] = -1;
		local_ref[0][0][1] = -1;
		local_ref[0][0][2] = -1;

		local_ref[0][1][0] = 1;
		local_ref[0][1][1] = -1;
		local_ref[0][1][2] = -1;

		local_ref[0][2][0] = 1;
		local_ref[0][2][1] = 1;
		local_ref[0][2][2] = -1;

		local_ref[0][3][0] = -1;
		local_ref[0][3][1] = 1;
		local_ref[0][3][2] = -1;

		local_ref[0][4][0] = -1;
		local_ref[0][4][1] = -1+step;
		local_ref[0][4][2] = -1+1*3*step/double(4);

		local_ref[0][5][0] = 1;
		local_ref[0][5][1] = -1+step;
		local_ref[0][5][2] = -1+1*3*step/double(4);

		local_ref[0][6][0] = 1;
		local_ref[0][6][1] = -1+2*step;
		local_ref[0][6][2] = -1+1*3*step/double(4);

		local_ref[0][7][0] = -1;
		local_ref[0][7][1] = -1+2*step;
		local_ref[0][7][2] = -1+1*3*step/double(4);

		//element 1
		local_ref[1][0][0] = -1;
		local_ref[1][0][1] = -1+step;
		local_ref[1][0][2] = -1+1*3*step/double(4);

		local_ref[1][1][0] = 1;
		local_ref[1][1][1] = -1+step;
		local_ref[1][1][2] = -1+1*3*step/double(4);

		local_ref[1][2][0] = 1;
		local_ref[1][2][1] = -1+2*step;
		local_ref[1][2][2] = -1+1*3*step/double(4);

		local_ref[1][3][0] = -1;
		local_ref[1][3][1] = -1+2*step;
		local_ref[1][3][2] = -1+1*3*step/double(4);

		local_ref[1][4][0] = -1+step;
		local_ref[1][4][1] = -1+step;
		local_ref[1][4][2] = -1+3*3*step/double(8);//-1+2*3*step/double(4);

		local_ref[1][5][0] = -1+2*step;
		local_ref[1][5][1] = -1+step;
		local_ref[1][5][2] = -1+3*3*step/double(8);//-1+2*3*step/double(4);

		local_ref[1][6][0] = -1+2*step;
		local_ref[1][6][1] = -1+2*step;
		local_ref[1][6][2] = -1+3*3*step/double(8);//-1+2*3*step/double(4);

		local_ref[1][7][0] = -1+step;
		local_ref[1][7][1] = -1+2*step;
		local_ref[1][7][2] = -1+3*3*step/double(8);//-1+2*3*step/double(4);


		//element 2
		local_ref[2][0][0] = -1;
		local_ref[2][0][1] = -1;
		local_ref[2][0][2] = -1;

		local_ref[2][1][0] =  1;
		local_ref[2][1][1] = -1;
		local_ref[2][1][2] = -1;

		local_ref[2][2][0] =  1;
		local_ref[2][2][1] = -1+step;
		local_ref[2][2][2] = -1+1*3*step/double(4);

		local_ref[2][3][0] = -1;
		local_ref[2][3][1] = -1+step;
		local_ref[2][3][2] = -1+1*3*step/double(4);

		local_ref[2][4][0] = -1+step;
		local_ref[2][4][1] = -1;
		local_ref[2][4][2] = -1+1*3*step/double(4);

		local_ref[2][5][0] = -1+2*step;
		local_ref[2][5][1] = -1;
		local_ref[2][5][2] = -1+1*3*step/double(4);

		local_ref[2][6][0] = -1+2*step;
		local_ref[2][6][1] = -1+1*step;
		local_ref[2][6][2] = -1+3*3*step/double(8);//-1+2*3*step/double(4);

		local_ref[2][7][0] = -1+step;
		local_ref[2][7][1] = -1+1*step;
		local_ref[2][7][2] = -1+3*3*step/double(8);//-1+2*3*step/double(4);

		//element 3
		local_ref[3][0][0] = -1;
		local_ref[3][0][1] = -1+2*step;
		local_ref[3][0][2] = -1+1*3*step/double(4);

		local_ref[3][1][0] =  1;
		local_ref[3][1][1] = -1+2*step;
		local_ref[3][1][2] = -1+1*3*step/double(4);

		local_ref[3][2][0] =  1;
		local_ref[3][2][1] =  1;
		local_ref[3][2][2] = -1;

		local_ref[3][3][0] = -1;
		local_ref[3][3][1] =  1;
		local_ref[3][3][2] = -1;

		local_ref[3][4][0] = -1+step;
		local_ref[3][4][1] = -1+2*step;
		local_ref[3][4][2] = -1+3*3*step/double(8);//-1+2*3*step/double(4);

		local_ref[3][5][0] = -1+2*step;
		local_ref[3][5][1] = -1+2*step;
		local_ref[3][5][2] = -1+3*3*step/double(8);//-1+2*3*step/double(4);

		local_ref[3][6][0] = -1+2*step;
		local_ref[3][6][1] = 1;
		local_ref[3][6][2] = -1+1*3*step/double(4);

		local_ref[3][7][0] = -1+step;
		local_ref[3][7][1] = 1;
		local_ref[3][7][2] = -1+1*3*step/double(4);

		//element 4
		local_ref[4][0][0] = -1+step;
		local_ref[4][0][1] = -1+step;
		local_ref[4][0][2] = -1+3*3*step/double(8);

		local_ref[4][1][0] = -1+2*step;
		local_ref[4][1][1] = -1+step;
		local_ref[4][1][2] = -1+3*3*step/double(8);

		local_ref[4][2][0] = -1+2*step;
		local_ref[4][2][1] = -1+2*step;
		local_ref[4][2][2] = -1+3*3*step/double(8);

		local_ref[4][3][0] = -1+step;
		local_ref[4][3][1] = -1+2*step;
		local_ref[4][3][2] = -1+3*3*step/double(8);

		local_ref[4][4][0] = -1+step;
		local_ref[4][4][1] = -1+step;
		local_ref[4][4][2] = -1+2*3*step/double(4);

		local_ref[4][5][0] = -1+2*step;
		local_ref[4][5][1] = -1+step;
		local_ref[4][5][2] = -1+2*3*step/double(4);

		local_ref[4][6][0] = -1+2*step;
		local_ref[4][6][1] = -1+2*step;
		local_ref[4][6][2] = -1+2*3*step/double(4);

		local_ref[4][7][0] = -1+step;
		local_ref[4][7][1] = -1+2*step;
		local_ref[4][7][2] = -1+2*3*step/double(4);

		//element 5
		local_ref[5][0][0] = -1+step;
		local_ref[5][0][1] = -1+step;
		local_ref[5][0][2] = -1+2*3*step/double(4);

		local_ref[5][1][0] = -1+2*step;
		local_ref[5][1][1] = -1+step;
		local_ref[5][1][2] = -1+2*3*step/double(4);

		local_ref[5][2][0] = -1+2*step;
		local_ref[5][2][1] = -1+2*step;
		local_ref[5][2][2] = -1+2*3*step/double(4);

		local_ref[5][3][0] = -1+step;
		local_ref[5][3][1] = -1+2*step;
		local_ref[5][3][2] = -1+2*3*step/double(4);

		local_ref[5][4][0] = -1+step;
		local_ref[5][4][1] = -1+step;
		local_ref[5][4][2] = -1+3*3*step/double(4);

		local_ref[5][5][0] = -1+2*step;
		local_ref[5][5][1] = -1+step;
		local_ref[5][5][2] = -1+3*3*step/double(4);

		local_ref[5][6][0] = -1+2*step;
		local_ref[5][6][1] = -1+2*step;
		local_ref[5][6][2] = -1+3*3*step/double(4);

		local_ref[5][7][0] = -1+step;
		local_ref[5][7][1] = -1+2*step;
		local_ref[5][7][2] = -1+3*3*step/double(4);

		//element 6
		local_ref[6][0][0] = -1+step;
		local_ref[6][0][1] = -1+step;
		local_ref[6][0][2] = -1+3*3*step/double(4);

		local_ref[6][1][0] = -1+2*step;
		local_ref[6][1][1] = -1+step;
		local_ref[6][1][2] = -1+3*3*step/double(4);

		local_ref[6][2][0] = -1+2*step;
		local_ref[6][2][1] = -1+2*step;
		local_ref[6][2][2] = -1+3*3*step/double(4);

		local_ref[6][3][0] = -1+step;
		local_ref[6][3][1] = -1+2*step;
		local_ref[6][3][2] = -1+3*3*step/double(4);

		local_ref[6][4][0] = -1+step;
		local_ref[6][4][1] = -1+step;
		local_ref[6][4][2] = 1;

		local_ref[6][5][0] = -1+2*step;
		local_ref[6][5][1] = -1+step;
		local_ref[6][5][2] = 1;

		local_ref[6][6][0] = -1+2*step;
		local_ref[6][6][1] = -1+2*step;
		local_ref[6][6][2] = 1;

		local_ref[6][7][0] = -1+step;
		local_ref[6][7][1] = -1+2*step;
		local_ref[6][7][2] = 1;

		//element 7
		local_ref[7][0][0] = -1;
		local_ref[7][0][1] = -1;
		local_ref[7][0][2] = -1;

		local_ref[7][1][0] = -1+step;
		local_ref[7][1][1] = -1;
		local_ref[7][1][2] = -1+1*3*step/double(4);

		local_ref[7][2][0] = -1+step;
		local_ref[7][2][1] = -1+step;
		local_ref[7][2][2] = -1+3*3*step/double(8);

		local_ref[7][3][0] = -1;
		local_ref[7][3][1] = -1+step;
		local_ref[7][3][2] = -1+1*3*step/double(4);

		local_ref[7][4][0] = -1;
		local_ref[7][4][1] = -1;
		local_ref[7][4][2] = -1+step;//-1+2*3*step/double(4);

		local_ref[7][5][0] = -1+step;
		local_ref[7][5][1] = -1;
		local_ref[7][5][2] = -1+2*3*step/double(4);

		local_ref[7][6][0] = -1+step;
		local_ref[7][6][1] = -1+step;
		local_ref[7][6][2] = -1+2*3*step/double(4);

		local_ref[7][7][0] = -1;
		local_ref[7][7][1] = -1+step;
		local_ref[7][7][2] = -1+2*3*step/double(4);


		//element 8
		local_ref[8][0][0] = -1+step;
		local_ref[8][0][1] = -1;
		local_ref[8][0][2] = -1+1*3*step/double(4);;

		local_ref[8][1][0] = -1+2*step;
		local_ref[8][1][1] = -1;
		local_ref[8][1][2] = -1+1*3*step/double(4);;

		local_ref[8][2][0] = -1+2*step;
		local_ref[8][2][1] = -1+step;
		local_ref[8][2][2] = -1+3*3*step/double(8);

		local_ref[8][3][0] = -1+step;
		local_ref[8][3][1] = -1+step;
		local_ref[8][3][2] = -1+3*3*step/double(8);

		local_ref[8][4][0] = -1+step;
		local_ref[8][4][1] = -1;
		local_ref[8][4][2] = -1+2*3*step/double(4);

		local_ref[8][5][0] = -1+2*step;
		local_ref[8][5][1] = -1;
		local_ref[8][5][2] = -1+2*3*step/double(4);

		local_ref[8][6][0] = -1+2*step;
		local_ref[8][6][1] = -1+1*step;
		local_ref[8][6][2] = -1+2*3*step/double(4);

		local_ref[8][7][0] = -1+step;
		local_ref[8][7][1] = -1+1*step;
		local_ref[8][7][2] = -1+2*3*step/double(4);

		//element 9
		local_ref[9][0][0] = -1+2*step;
		local_ref[9][0][1] = -1;
		local_ref[9][0][2] = -1+1*3*step/double(4);

		local_ref[9][1][0] = 1;
		local_ref[9][1][1] = -1;
		local_ref[9][1][2] = -1;

		local_ref[9][2][0] =  1;
		local_ref[9][2][1] = -1+step;
		local_ref[9][2][2] = -1+1*3*step/double(4);

		local_ref[9][3][0] = -1+2*step;
		local_ref[9][3][1] = -1+1*step;
		local_ref[9][3][2] = -1+3*3*step/double(8);

		local_ref[9][4][0] = -1+2*step;
		local_ref[9][4][1] = -1;
		local_ref[9][4][2] = -1+2*3*step/double(4);

		local_ref[9][5][0] = 1;
		local_ref[9][5][1] = -1;
		local_ref[9][5][2] = -1+step;//-1+2*3*step/double(4);

		local_ref[9][6][0] = 1;
		local_ref[9][6][1] = -1+1*step;
		local_ref[9][6][2] = -1+2*3*step/double(4);

		local_ref[9][7][0] = -1+2*step;
		local_ref[9][7][1] = -1+1*step;
		local_ref[9][7][2] = -1+2*3*step/double(4);

		//element 10
		local_ref[10][0][0] = -1+2*step;
		local_ref[10][0][1] = -1+1*step;
		local_ref[10][0][2] = -1+3*3*step/double(8);

		local_ref[10][1][0] = 1;
		local_ref[10][1][1] = -1+1*step;
		local_ref[10][1][2] = -1+1*3*step/double(4);

		local_ref[10][2][0] =  1;
		local_ref[10][2][1] = -1+2*step;
		local_ref[10][2][2] = -1+1*3*step/double(4);

		local_ref[10][3][0] = -1+2*step;
		local_ref[10][3][1] = -1+2*step;
		local_ref[10][3][2] = -1+3*3*step/double(8);

		local_ref[10][4][0] = -1+2*step;
		local_ref[10][4][1] = -1+1*step;
		local_ref[10][4][2] = -1+2*3*step/double(4);

		local_ref[10][5][0] = 1;
		local_ref[10][5][1] = -1+1*step;
		local_ref[10][5][2] = -1+2*3*step/double(4);

		local_ref[10][6][0] = 1;
		local_ref[10][6][1] = -1+2*step;
		local_ref[10][6][2] = -1+2*3*step/double(4);

		local_ref[10][7][0] = -1+2*step;
		local_ref[10][7][1] = -1+2*step;
		local_ref[10][7][2] = -1+2*3*step/double(4);

		//element 11
		local_ref[11][0][0] = -1+2*step;
		local_ref[11][0][1] = -1+2*step;
		local_ref[11][0][2] = -1+3*3*step/double(8);

		local_ref[11][1][0] = 1;
		local_ref[11][1][1] = -1+2*step;
		local_ref[11][1][2] = -1+1*3*step/double(4);

		local_ref[11][2][0] =  1;
		local_ref[11][2][1] =  1;
		local_ref[11][2][2] = -1;

		local_ref[11][3][0] = -1+2*step;
		local_ref[11][3][1] =  1;
		local_ref[11][3][2] = -1+1*3*step/double(4);

		local_ref[11][4][0] = -1+2*step;
		local_ref[11][4][1] = -1+2*step;
		local_ref[11][4][2] = -1+2*3*step/double(4);

		local_ref[11][5][0] = 1;
		local_ref[11][5][1] = -1+2*step;
		local_ref[11][5][2] = -1+2*3*step/double(4);

		local_ref[11][6][0] = 1;
		local_ref[11][6][1] = 1;
		local_ref[11][6][2] = -1+step;//-1+2*3*step/double(4);

		local_ref[11][7][0] = -1+2*step;
		local_ref[11][7][1] = 1;
		local_ref[11][7][2] = -1+2*3*step/double(4);

		//element 12
		local_ref[12][0][0] = -1+1*step;
		local_ref[12][0][1] = -1+2*step;
		local_ref[12][0][2] = -1+3*3*step/double(8);

		local_ref[12][1][0] = -1+2*step;
		local_ref[12][1][1] = -1+2*step;
		local_ref[12][1][2] = -1+3*3*step/double(8);

		local_ref[12][2][0] =  -1+2*step;
		local_ref[12][2][1] =  1;
		local_ref[12][2][2] = -1+1*3*step/double(4);

		local_ref[12][3][0] = -1+1*step;
		local_ref[12][3][1] =  1;
		local_ref[12][3][2] = -1+1*3*step/double(4);

		local_ref[12][4][0] = -1+1*step;
		local_ref[12][4][1] = -1+2*step;
		local_ref[12][4][2] = -1+2*3*step/double(4);

		local_ref[12][5][0] = -1+2*step;
		local_ref[12][5][1] = -1+2*step;
		local_ref[12][5][2] = -1+2*3*step/double(4);

		local_ref[12][6][0] = -1+2*step;
		local_ref[12][6][1] = 1;
		local_ref[12][6][2] = -1+2*3*step/double(4);

		local_ref[12][7][0] = -1+1*step;
		local_ref[12][7][1] = 1;
		local_ref[12][7][2] = -1+2*3*step/double(4);

		//element 13
		local_ref[13][0][0] = -1;
		local_ref[13][0][1] = -1+2*step;
		local_ref[13][0][2] = -1+1*3*step/double(4);

		local_ref[13][1][0] = -1+1*step;
		local_ref[13][1][1] = -1+2*step;
		local_ref[13][1][2] = -1+3*3*step/double(8);

		local_ref[13][2][0] =  -1+1*step;
		local_ref[13][2][1] =  1;
		local_ref[13][2][2] = -1+1*3*step/double(4);

		local_ref[13][3][0] = -1;
		local_ref[13][3][1] =  1;
		local_ref[13][3][2] = -1;

		local_ref[13][4][0] = -1;
		local_ref[13][4][1] = -1+2*step;
		local_ref[13][4][2] = -1+2*3*step/double(4);

		local_ref[13][5][0] = -1+1*step;
		local_ref[13][5][1] = -1+2*step;
		local_ref[13][5][2] = -1+2*3*step/double(4);

		local_ref[13][6][0] = -1+1*step;
		local_ref[13][6][1] = 1;
		local_ref[13][6][2] = -1+2*3*step/double(4);

		local_ref[13][7][0] = -1;
		local_ref[13][7][1] = 1;
		local_ref[13][7][2] = -1+step;//-1+2*3*step/double(4);

		//element 14
		local_ref[14][0][0] = -1;
		local_ref[14][0][1] = -1+1*step;
		local_ref[14][0][2] = -1+1*3*step/double(4);

		local_ref[14][1][0] = -1+1*step;
		local_ref[14][1][1] = -1+1*step;
		local_ref[14][1][2] = -1+3*3*step/double(8);

		local_ref[14][2][0] =  -1+1*step;
		local_ref[14][2][1] =  -1+2*step;
		local_ref[14][2][2] = -1+3*3*step/double(8);

		local_ref[14][3][0] = -1;
		local_ref[14][3][1] =  -1+2*step;
		local_ref[14][3][2] = -1+1*3*step/double(4);

		local_ref[14][4][0] = -1;
		local_ref[14][4][1] = -1+1*step;
		local_ref[14][4][2] = -1+2*3*step/double(4);

		local_ref[14][5][0] = -1+1*step;
		local_ref[14][5][1] = -1+1*step;
		local_ref[14][5][2] = -1+2*3*step/double(4);

		local_ref[14][6][0] = -1+1*step;
		local_ref[14][6][1] = -1+2*step;
		local_ref[14][6][2] = -1+2*3*step/double(4);

		local_ref[14][7][0] = -1;
		local_ref[14][7][1] = -1+2*step;
		local_ref[14][7][2] = -1+2*3*step/double(4);

		//element 15
		local_ref[15][0][0] = -1;
		local_ref[15][0][1] = -1;
		local_ref[15][0][2] = -1+step;

		local_ref[15][1][0] = -1+step;
		local_ref[15][1][1] = -1;
		local_ref[15][1][2] = -1+2*3*step/double(4);

		local_ref[15][2][0] = -1+step;
		local_ref[15][2][1] = -1+step;
		local_ref[15][2][2] = -1+2*3*step/double(4);

		local_ref[15][3][0] = -1;
		local_ref[15][3][1] = -1+step;
		local_ref[15][3][2] = -1+2*3*step/double(4);

		local_ref[15][4][0] = -1;
		local_ref[15][4][1] = -1;
		local_ref[15][4][2] = -1+2*step;//-1+2*3*step/double(4);

		local_ref[15][5][0] = -1+step;
		local_ref[15][5][1] = -1;
		local_ref[15][5][2] = -1+3*3*step/double(4);

		local_ref[15][6][0] = -1+step;
		local_ref[15][6][1] = -1+step;
		local_ref[15][6][2] = -1+3*3*step/double(4);

		local_ref[15][7][0] = -1;
		local_ref[15][7][1] = -1+step;
		local_ref[15][7][2] = -1+3*3*step/double(4);


		//element 16
		local_ref[16][0][0] = -1+step;
		local_ref[16][0][1] = -1;
		local_ref[16][0][2] = -1+2*3*step/double(4);

		local_ref[16][1][0] = -1+2*step;
		local_ref[16][1][1] = -1;
		local_ref[16][1][2] = -1+2*3*step/double(4);

		local_ref[16][2][0] = -1+2*step;
		local_ref[16][2][1] = -1+step;
		local_ref[16][2][2] = -1+2*3*step/double(4);

		local_ref[16][3][0] = -1+step;
		local_ref[16][3][1] = -1+step;
		local_ref[16][3][2] = -1+2*3*step/double(4);

		local_ref[16][4][0] = -1+step;
		local_ref[16][4][1] = -1;
		local_ref[16][4][2] = -1+3*3*step/double(4);

		local_ref[16][5][0] = -1+2*step;
		local_ref[16][5][1] = -1;
		local_ref[16][5][2] = -1+3*3*step/double(4);

		local_ref[16][6][0] = -1+2*step;
		local_ref[16][6][1] = -1+1*step;
		local_ref[16][6][2] = -1+3*3*step/double(4);

		local_ref[16][7][0] = -1+step;
		local_ref[16][7][1] = -1+1*step;
		local_ref[16][7][2] = -1+3*3*step/double(4);

		//element 17
		local_ref[17][0][0] = -1+2*step;
		local_ref[17][0][1] = -1;
		local_ref[17][0][2] = -1+2*3*step/double(4);

		local_ref[17][1][0] = 1;
		local_ref[17][1][1] = -1;
		local_ref[17][1][2] = -1+step;

		local_ref[17][2][0] =  1;
		local_ref[17][2][1] = -1+step;
		local_ref[17][2][2] = -1+2*3*step/double(4);

		local_ref[17][3][0] = -1+2*step;
		local_ref[17][3][1] = -1+1*step;
		local_ref[17][3][2] = -1+2*3*step/double(4);

		local_ref[17][4][0] = -1+2*step;
		local_ref[17][4][1] = -1;
		local_ref[17][4][2] = -1+3*3*step/double(4);

		local_ref[17][5][0] = 1;
		local_ref[17][5][1] = -1;
		local_ref[17][5][2] = -1+2*step;//-1+2*3*step/double(4);

		local_ref[17][6][0] = 1;
		local_ref[17][6][1] = -1+1*step;
		local_ref[17][6][2] = -1+3*3*step/double(4);

		local_ref[17][7][0] = -1+2*step;
		local_ref[17][7][1] = -1+1*step;
		local_ref[17][7][2] = -1+3*3*step/double(4);

		//element 18
		local_ref[18][0][0] = -1+2*step;
		local_ref[18][0][1] = -1+1*step;
		local_ref[18][0][2] = -1+2*3*step/double(4);

		local_ref[18][1][0] = 1;
		local_ref[18][1][1] = -1+1*step;
		local_ref[18][1][2] = -1+2*3*step/double(4);

		local_ref[18][2][0] =  1;
		local_ref[18][2][1] = -1+2*step;
		local_ref[18][2][2] = -1+2*3*step/double(4);

		local_ref[18][3][0] = -1+2*step;
		local_ref[18][3][1] = -1+2*step;
		local_ref[18][3][2] = -1+2*3*step/double(4);

		local_ref[18][4][0] = -1+2*step;
		local_ref[18][4][1] = -1+1*step;
		local_ref[18][4][2] = -1+3*3*step/double(4);

		local_ref[18][5][0] = 1;
		local_ref[18][5][1] = -1+1*step;
		local_ref[18][5][2] = -1+3*3*step/double(4);

		local_ref[18][6][0] = 1;
		local_ref[18][6][1] = -1+2*step;
		local_ref[18][6][2] = -1+3*3*step/double(4);

		local_ref[18][7][0] = -1+2*step;
		local_ref[18][7][1] = -1+2*step;
		local_ref[18][7][2] = -1+3*3*step/double(4);

		//element 19
		local_ref[19][0][0] = -1+2*step;
		local_ref[19][0][1] = -1+2*step;
		local_ref[19][0][2] = -1+2*3*step/double(4);

		local_ref[19][1][0] = 1;
		local_ref[19][1][1] = -1+2*step;
		local_ref[19][1][2] = -1+2*3*step/double(4);

		local_ref[19][2][0] =  1;
		local_ref[19][2][1] =  1;
		local_ref[19][2][2] = -1+step;

		local_ref[19][3][0] = -1+2*step;
		local_ref[19][3][1] =  1;
		local_ref[19][3][2] = -1+2*3*step/double(4);

		local_ref[19][4][0] = -1+2*step;
		local_ref[19][4][1] = -1+2*step;
		local_ref[19][4][2] = -1+3*3*step/double(4);

		local_ref[19][5][0] = 1;
		local_ref[19][5][1] = -1+2*step;
		local_ref[19][5][2] = -1+3*3*step/double(4);

		local_ref[19][6][0] = 1;
		local_ref[19][6][1] = 1;
		local_ref[19][6][2] = -1+2*step;

		local_ref[19][7][0] = -1+2*step;
		local_ref[19][7][1] = 1;
		local_ref[19][7][2] = -1+3*3*step/double(4);

		//element 20
		local_ref[20][0][0] = -1+1*step;
		local_ref[20][0][1] = -1+2*step;
		local_ref[20][0][2] = -1+2*3*step/double(4);

		local_ref[20][1][0] = -1+2*step;
		local_ref[20][1][1] = -1+2*step;
		local_ref[20][1][2] = -1+2*3*step/double(4);

		local_ref[20][2][0] =  -1+2*step;
		local_ref[20][2][1] =  1;
		local_ref[20][2][2] = -1+2*3*step/double(4);

		local_ref[20][3][0] = -1+1*step;
		local_ref[20][3][1] =  1;
		local_ref[20][3][2] = -1+2*3*step/double(4);

		local_ref[20][4][0] = -1+1*step;
		local_ref[20][4][1] = -1+2*step;
		local_ref[20][4][2] = -1+3*3*step/double(4);

		local_ref[20][5][0] = -1+2*step;
		local_ref[20][5][1] = -1+2*step;
		local_ref[20][5][2] = -1+3*3*step/double(4);

		local_ref[20][6][0] = -1+2*step;
		local_ref[20][6][1] = 1;
		local_ref[20][6][2] = -1+3*3*step/double(4);

		local_ref[20][7][0] = -1+1*step;
		local_ref[20][7][1] = 1;
		local_ref[20][7][2] = -1+3*3*step/double(4);

		//element 21
		local_ref[21][0][0] = -1;
		local_ref[21][0][1] = -1+2*step;
		local_ref[21][0][2] = -1+2*3*step/double(4);

		local_ref[21][1][0] = -1+1*step;
		local_ref[21][1][1] = -1+2*step;
		local_ref[21][1][2] = -1+2*3*step/double(4);

		local_ref[21][2][0] =  -1+1*step;
		local_ref[21][2][1] =  1;
		local_ref[21][2][2] = -1+2*3*step/double(4);

		local_ref[21][3][0] = -1;
		local_ref[21][3][1] =  1;
		local_ref[21][3][2] = -1+step;

		local_ref[21][4][0] = -1;
		local_ref[21][4][1] = -1+2*step;
		local_ref[21][4][2] = -1+3*3*step/double(4);

		local_ref[21][5][0] = -1+1*step;
		local_ref[21][5][1] = -1+2*step;
		local_ref[21][5][2] = -1+3*3*step/double(4);

		local_ref[21][6][0] = -1+1*step;
		local_ref[21][6][1] = 1;
		local_ref[21][6][2] = -1+3*3*step/double(4);

		local_ref[21][7][0] = -1;
		local_ref[21][7][1] = 1;
		local_ref[21][7][2] = -1+2*step;//-1+2*3*step/double(4);

		//element 22
		local_ref[22][0][0] = -1;
		local_ref[22][0][1] = -1+1*step;
		local_ref[22][0][2] = -1+2*3*step/double(4);

		local_ref[22][1][0] = -1+1*step;
		local_ref[22][1][1] = -1+1*step;
		local_ref[22][1][2] = -1+2*3*step/double(4);

		local_ref[22][2][0] =  -1+1*step;
		local_ref[22][2][1] =  -1+2*step;
		local_ref[22][2][2] = -1+2*3*step/double(4);

		local_ref[22][3][0] = -1;
		local_ref[22][3][1] =  -1+2*step;
		local_ref[22][3][2] = -1+2*3*step/double(4);

		local_ref[22][4][0] = -1;
		local_ref[22][4][1] = -1+1*step;
		local_ref[22][4][2] = -1+3*3*step/double(4);

		local_ref[22][5][0] = -1+1*step;
		local_ref[22][5][1] = -1+1*step;
		local_ref[22][5][2] = -1+3*3*step/double(4);

		local_ref[22][6][0] = -1+1*step;
		local_ref[22][6][1] = -1+2*step;
		local_ref[22][6][2] = -1+3*3*step/double(4);

		local_ref[22][7][0] = -1;
		local_ref[22][7][1] = -1+2*step;
		local_ref[22][7][2] = -1+3*3*step/double(4);

		//element 23
		local_ref[23][0][0] = -1;
		local_ref[23][0][1] = -1;
		local_ref[23][0][2] = -1+2*step;

		local_ref[23][1][0] = -1+step;
		local_ref[23][1][1] = -1;
		local_ref[23][1][2] = -1+3*3*step/double(4);

		local_ref[23][2][0] = -1+step;
		local_ref[23][2][1] = -1+step;
		local_ref[23][2][2] = -1+3*3*step/double(4);

		local_ref[23][3][0] = -1;
		local_ref[23][3][1] = -1+step;
		local_ref[23][3][2] = -1+3*3*step/double(4);

		local_ref[23][4][0] = -1;
		local_ref[23][4][1] = -1;
		local_ref[23][4][2] = 1;

		local_ref[23][5][0] = -1+step;
		local_ref[23][5][1] = -1;
		local_ref[23][5][2] = 1;

		local_ref[23][6][0] = -1+step;
		local_ref[23][6][1] = -1+step;
		local_ref[23][6][2] = 1;

		local_ref[23][7][0] = -1;
		local_ref[23][7][1] = -1+step;
		local_ref[23][7][2] = 1;


		//element 24
		local_ref[24][0][0] = -1+step;
		local_ref[24][0][1] = -1;
		local_ref[24][0][2] = -1+3*3*step/double(4);

		local_ref[24][1][0] = -1+2*step;
		local_ref[24][1][1] = -1;
		local_ref[24][1][2] = -1+3*3*step/double(4);

		local_ref[24][2][0] = -1+2*step;
		local_ref[24][2][1] = -1+step;
		local_ref[24][2][2] = -1+3*3*step/double(4);

		local_ref[24][3][0] = -1+step;
		local_ref[24][3][1] = -1+step;
		local_ref[24][3][2] = -1+3*3*step/double(4);

		local_ref[24][4][0] = -1+step;
		local_ref[24][4][1] = -1;
		local_ref[24][4][2] = 1;

		local_ref[24][5][0] = -1+2*step;
		local_ref[24][5][1] = -1;
		local_ref[24][5][2] = 1;

		local_ref[24][6][0] = -1+2*step;
		local_ref[24][6][1] = -1+1*step;
		local_ref[24][6][2] = 1;

		local_ref[24][7][0] = -1+step;
		local_ref[24][7][1] = -1+1*step;
		local_ref[24][7][2] = 1;

		//element 25
		local_ref[25][0][0] = -1+2*step;
		local_ref[25][0][1] = -1;
		local_ref[25][0][2] = -1+3*3*step/double(4);

		local_ref[25][1][0] = 1;
		local_ref[25][1][1] = -1;
		local_ref[25][1][2] = -1+2*step;

		local_ref[25][2][0] =  1;
		local_ref[25][2][1] = -1+step;
		local_ref[25][2][2] = -1+3*3*step/double(4);

		local_ref[25][3][0] = -1+2*step;
		local_ref[25][3][1] = -1+1*step;
		local_ref[25][3][2] = -1+3*3*step/double(4);

		local_ref[25][4][0] = -1+2*step;
		local_ref[25][4][1] = -1;
		local_ref[25][4][2] = 1;

		local_ref[25][5][0] = 1;
		local_ref[25][5][1] = -1;
		local_ref[25][5][2] = 1;

		local_ref[25][6][0] = 1;
		local_ref[25][6][1] = -1+1*step;
		local_ref[25][6][2] = 1;

		local_ref[25][7][0] = -1+2*step;
		local_ref[25][7][1] = -1+1*step;
		local_ref[25][7][2] = 1;

		//element 26
		local_ref[26][0][0] = -1+2*step;
		local_ref[26][0][1] = -1+1*step;
		local_ref[26][0][2] = -1+3*3*step/double(4);

		local_ref[26][1][0] = 1;
		local_ref[26][1][1] = -1+1*step;
		local_ref[26][1][2] = -1+3*3*step/double(4);

		local_ref[26][2][0] =  1;
		local_ref[26][2][1] = -1+2*step;
		local_ref[26][2][2] = -1+3*3*step/double(4);

		local_ref[26][3][0] = -1+2*step;
		local_ref[26][3][1] = -1+2*step;
		local_ref[26][3][2] = -1+3*3*step/double(4);

		local_ref[26][4][0] = -1+2*step;
		local_ref[26][4][1] = -1+1*step;
		local_ref[26][4][2] = 1;

		local_ref[26][5][0] = 1;
		local_ref[26][5][1] = -1+1*step;
		local_ref[26][5][2] = 1;

		local_ref[26][6][0] = 1;
		local_ref[26][6][1] = -1+2*step;
		local_ref[26][6][2] = 1;

		local_ref[26][7][0] = -1+2*step;
		local_ref[26][7][1] = -1+2*step;
		local_ref[26][7][2] = 1;

		//element 27
		local_ref[27][0][0] = -1+2*step;
		local_ref[27][0][1] = -1+2*step;
		local_ref[27][0][2] = -1+3*3*step/double(4);

		local_ref[27][1][0] = 1;
		local_ref[27][1][1] = -1+2*step;
		local_ref[27][1][2] = -1+3*3*step/double(4);

		local_ref[27][2][0] =  1;
		local_ref[27][2][1] =  1;
		local_ref[27][2][2] = -1+2*step;

		local_ref[27][3][0] = -1+2*step;
		local_ref[27][3][1] =  1;
		local_ref[27][3][2] = -1+3*3*step/double(4);

		local_ref[27][4][0] = -1+2*step;
		local_ref[27][4][1] = -1+2*step;
		local_ref[27][4][2] = 1;

		local_ref[27][5][0] = 1;
		local_ref[27][5][1] = -1+2*step;
		local_ref[27][5][2] = 1;

		local_ref[27][6][0] = 1;
		local_ref[27][6][1] = 1;
		local_ref[27][6][2] = 1;

		local_ref[27][7][0] = -1+2*step;
		local_ref[27][7][1] = 1;
		local_ref[27][7][2] = 1;

		//element 28
		local_ref[28][0][0] = -1+1*step;
		local_ref[28][0][1] = -1+2*step;
		local_ref[28][0][2] = -1+3*3*step/double(4);

		local_ref[28][1][0] = -1+2*step;
		local_ref[28][1][1] = -1+2*step;
		local_ref[28][1][2] = -1+3*3*step/double(4);

		local_ref[28][2][0] =  -1+2*step;
		local_ref[28][2][1] =  1;
		local_ref[28][2][2] = -1+3*3*step/double(4);

		local_ref[28][3][0] = -1+1*step;
		local_ref[28][3][1] =  1;
		local_ref[28][3][2] = -1+3*3*step/double(4);

		local_ref[28][4][0] = -1+1*step;
		local_ref[28][4][1] = -1+2*step;
		local_ref[28][4][2] = 1;

		local_ref[28][5][0] = -1+2*step;
		local_ref[28][5][1] = -1+2*step;
		local_ref[28][5][2] = 1;

		local_ref[28][6][0] = -1+2*step;
		local_ref[28][6][1] = 1;
		local_ref[28][6][2] = 1;

		local_ref[28][7][0] = -1+1*step;
		local_ref[28][7][1] = 1;
		local_ref[28][7][2] = 1;

		//element 29
		local_ref[29][0][0] = -1;
		local_ref[29][0][1] = -1+2*step;
		local_ref[29][0][2] = -1+3*3*step/double(4);

		local_ref[29][1][0] = -1+1*step;
		local_ref[29][1][1] = -1+2*step;
		local_ref[29][1][2] = -1+3*3*step/double(4);

		local_ref[29][2][0] =  -1+1*step;
		local_ref[29][2][1] =  1;
		local_ref[29][2][2] = -1+3*3*step/double(4);

		local_ref[29][3][0] = -1;
		local_ref[29][3][1] =  1;
		local_ref[29][3][2] = -1+2*step;

		local_ref[29][4][0] = -1;
		local_ref[29][4][1] = -1+2*step;
		local_ref[29][4][2] = 1;

		local_ref[29][5][0] = -1+1*step;
		local_ref[29][5][1] = -1+2*step;
		local_ref[29][5][2] = 1;

		local_ref[29][6][0] = -1+1*step;
		local_ref[29][6][1] = 1;
		local_ref[29][6][2] = 1;

		local_ref[29][7][0] = -1;
		local_ref[29][7][1] = 1;
		local_ref[29][7][2] = 1;

		//element 30
		local_ref[30][0][0] = -1;
		local_ref[30][0][1] = -1+1*step;
		local_ref[30][0][2] = -1+3*3*step/double(4);

		local_ref[30][1][0] = -1+1*step;
		local_ref[30][1][1] = -1+1*step;
		local_ref[30][1][2] = -1+3*3*step/double(4);

		local_ref[30][2][0] =  -1+1*step;
		local_ref[30][2][1] =  -1+2*step;
		local_ref[30][2][2] = -1+3*3*step/double(4);

		local_ref[30][3][0] = -1;
		local_ref[30][3][1] =  -1+2*step;
		local_ref[30][3][2] = -1+3*3*step/double(4);

		local_ref[30][4][0] = -1;
		local_ref[30][4][1] = -1+1*step;
		local_ref[30][4][2] = 1;

		local_ref[30][5][0] = -1+1*step;
		local_ref[30][5][1] = -1+1*step;
		local_ref[30][5][2] = 1;

		local_ref[30][6][0] = -1+1*step;
		local_ref[30][6][1] = -1+2*step;
		local_ref[30][6][2] = 1;

		local_ref[30][7][0] = -1;
		local_ref[30][7][1] = -1+2*step;
		local_ref[30][7][2] = 1;

	}

	//define the rotation of the reference element
	int rot[3];
	int sym[3];
	int id_node[8];
	vector<int> ord;

	if(elem->pad==126){
		//edge 0 1 2 3 4 5 6 7
		rot[0] = 0;
		rot[1] = 0;
		rot[2] = 0;

		sym[0] = 0;
		sym[1] = 0;
		sym[2] = 1;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}
	}else if(elem->pad==127){
		//edge 5 6 7 8 9 10 11
		rot[0] = 0;
		rot[1] = 0;
		rot[2] = 0;

		sym[0] = 0;
		sym[1] = 0;
		sym[2] = 0;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}
	}else if(elem->pad==128){
		//edge 0 5 8 4 1 9 3 11

		rot[0] = -1;
		rot[1] = 0;
		rot[2] = 0;

		sym[0] = 0;
		sym[1] = 0;
		sym[2] = 0;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}
	}else if(elem->pad==129){
		//edge 2 6 10 7 1 9 3 11
		rot[0] = 1;
		rot[1] = 0;
		rot[2] = 0;

		sym[0] = 0;
		sym[1] = 0;
		sym[2] = 0;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}
	}else if(elem->pad==130){
		//edge 1 5 9 6 0 8 10 2
		rot[0] = 0;
		rot[1] = -1;
		rot[2] = 0;

		sym[0] = 0;
		sym[1] = 0;
		sym[2] = 0;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}
	}else if(elem->pad==131){
		//edge 3 7 11 4 0 8 10 2
		rot[0] = 0;
		rot[1] = 1;
		rot[2] = 0;

		sym[0] = 0;
		sym[1] = 0;
		sym[2] = 0;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}
	}else{
		printf ("Error in template 9\n");
		exit (EXIT_FAILURE);
	}


	for(int iel = 0; iel < 31; iel++){
		double local_ref_x[8];
		double local_ref_y[8];
		double local_ref_z[8];
		for(int ino = 0; ino < 8; ino++){
			local_ref_x[ino] = local_ref[iel][ino][0];
			local_ref_y[ino] = local_ref[iel][ino][1];
			local_ref_z[ino] = local_ref[iel][ino][2];
		}
		ApplyElement(mesh, coords,  id,  iel,  id_node,  local_ref_x, local_ref_y, local_ref_z,  ord,  hash_nodes);
	}
	sc_array_reset(&toto);
}

void ApplyTemplate10(hexa_tree_t* mesh, std::vector<double>& coords, std::vector<int>& elements_ids, int iel, sc_hash_array_t* hash_nodes){

	sc_array_t toto;
	sc_array_init(&toto, sizeof(octant_t));

	octant_t* elemOrig = (octant_t*) sc_array_index(&mesh->elements, elements_ids[iel]);
	octant_t * elem = (octant_t*) sc_array_push(&toto);
	hexa_element_copy(elemOrig,elem);

	double step = double(2)/double(3);
	int id = elements_ids[iel];

	double cord_in_ref[3];
	cord_in_ref[0] = 0;
	cord_in_ref[1] = 0;
	cord_in_ref[2] = 0;

	//reference element edge 0
	double local_ref[17][8][3];
	if(true){
		//element 0
		local_ref[0][0][0] = -1;
		local_ref[0][0][1] = -1;
		local_ref[0][0][2] = -1;

		local_ref[0][1][0] = -1+step;
		local_ref[0][1][1] = -1;
		local_ref[0][1][2] = -1;

		local_ref[0][2][0] = -1+step;
		local_ref[0][2][1] = -1+1*3*step/double(4);
		local_ref[0][2][2] = -1;

		local_ref[0][3][0] = -1;
		local_ref[0][3][1] = -1+step;
		local_ref[0][3][2] = -1;

		local_ref[0][4][0] = -1;
		local_ref[0][4][1] = -1;
		local_ref[0][4][2] = -1+step;

		local_ref[0][5][0] = -1+step;
		local_ref[0][5][1] = -1;
		local_ref[0][5][2] = -1+1*3*step/double(4);

		local_ref[0][6][0] = -1+step;
		local_ref[0][6][1] = -1+1*3*step/double(4);
		local_ref[0][6][2] = -1+1*3*step/double(4);

		local_ref[0][7][0] = -1;
		local_ref[0][7][1] = -1+step;
		local_ref[0][7][2] = -1+step;

		//element 1
		local_ref[1][0][0] = -1+1*step;
		local_ref[1][0][1] = -1;
		local_ref[1][0][2] = -1;

		local_ref[1][1][0] = -1+2*step;
		local_ref[1][1][1] = -1;
		local_ref[1][1][2] = -1;

		local_ref[1][2][0] = -1+2*step;
		local_ref[1][2][1] = -1+1*3*step/double(4);
		local_ref[1][2][2] = -1;

		local_ref[1][3][0] = -1+1*step;
		local_ref[1][3][1] = -1+1*3*step/double(4);
		local_ref[1][3][2] = -1;

		local_ref[1][4][0] = -1+1*step;
		local_ref[1][4][1] = -1;
		local_ref[1][4][2] = -1+1*3*step/double(4);

		local_ref[1][5][0] = -1+2*step;
		local_ref[1][5][1] = -1;
		local_ref[1][5][2] = -1+1*3*step/double(4);

		local_ref[1][6][0] = -1+2*step;
		local_ref[1][6][1] = -1+1*3*step/double(4);
		local_ref[1][6][2] = -1+1*3*step/double(4);

		local_ref[1][7][0] = -1+1*step;
		local_ref[1][7][1] = -1+1*3*step/double(4);
		local_ref[1][7][2] = -1+1*3*step/double(4);


		//element 2
		local_ref[2][0][0] = -1+2*step;
		local_ref[2][0][1] = -1;
		local_ref[2][0][2] = -1;

		local_ref[2][1][0] =  1;
		local_ref[2][1][1] = -1;
		local_ref[2][1][2] = -1;

		local_ref[2][2][0] =  1;
		local_ref[2][2][1] = -1+step;
		local_ref[2][2][2] = -1;

		local_ref[2][3][0] = -1+2*step;
		local_ref[2][3][1] = -1+1*3*step/double(4);
		local_ref[2][3][2] = -1;

		local_ref[2][4][0] = -1+2*step;
		local_ref[2][4][1] = -1;
		local_ref[2][4][2] = -1+1*3*step/double(4);

		local_ref[2][5][0] =  1;
		local_ref[2][5][1] = -1;
		local_ref[2][5][2] =  -1+step;

		local_ref[2][6][0] = 1;
		local_ref[2][6][1] = -1+step;
		local_ref[2][6][2] = -1+step;

		local_ref[2][7][0] = -1+2*step;
		local_ref[2][7][1] = -1+1*3*step/double(4);
		local_ref[2][7][2] = -1+1*3*step/double(4);

		//element 3
		local_ref[3][0][0] = -1;
		local_ref[3][0][1] = -1;
		local_ref[3][0][2] = -1+step;

		local_ref[3][1][0] = -1+step;
		local_ref[3][1][1] = -1;
		local_ref[3][1][2] = -1+1*3*step/double(4);

		local_ref[3][2][0] = -1+step;
		local_ref[3][2][1] = -1+1*3*step/double(4);
		local_ref[3][2][2] = -1+1*3*step/double(4);

		local_ref[3][3][0] = -1;
		local_ref[3][3][1] = -1+step;
		local_ref[3][3][2] = -1+step;

		local_ref[3][4][0] = -1;
		local_ref[3][4][1] = -1;
		local_ref[3][4][2] = -1+2*step;

		local_ref[3][5][0] = -1+step;
		local_ref[3][5][1] = -1;
		local_ref[3][5][2] = -1+2*3*step/double(4);

		local_ref[3][6][0] = -1+step;
		local_ref[3][6][1] = -1+2*3*step/double(4);
		local_ref[3][6][2] = -1+2*3*step/double(4);

		local_ref[3][7][0] = -1;
		local_ref[3][7][1] = -1+2*step;
		local_ref[3][7][2] = -1+2*step;

		//element 4
		local_ref[4][0][0] = -1+step;
		local_ref[4][0][1] = -1;
		local_ref[4][0][2] = -1+1*3*step/double(4);

		local_ref[4][1][0] = -1+2*step;
		local_ref[4][1][1] = -1;
		local_ref[4][1][2] = -1+1*3*step/double(4);

		local_ref[4][2][0] = -1+2*step;
		local_ref[4][2][1] = -1+1*3*step/double(4);
		local_ref[4][2][2] = -1+1*3*step/double(4);

		local_ref[4][3][0] = -1+step;
		local_ref[4][3][1] = -1+1*3*step/double(4);
		local_ref[4][3][2] = -1+1*3*step/double(4);

		local_ref[4][4][0] = -1+step;
		local_ref[4][4][1] = -1;
		local_ref[4][4][2] = -1+2*3*step/double(4);

		local_ref[4][5][0] = -1+2*step;
		local_ref[4][5][1] = -1;
		local_ref[4][5][2] = -1+2*3*step/double(4);

		local_ref[4][6][0] = -1+2*step;
		local_ref[4][6][1] = -1+2*3*step/double(4);
		local_ref[4][6][2] = -1+2*3*step/double(4);

		local_ref[4][7][0] = -1+step;
		local_ref[4][7][1] = -1+2*3*step/double(4);
		local_ref[4][7][2] = -1+2*3*step/double(4);

		//element 5
		local_ref[5][0][0] = -1+2*step;
		local_ref[5][0][1] = -1;
		local_ref[5][0][2] = -1+1*3*step/double(4);

		local_ref[5][1][0] =  1;
		local_ref[5][1][1] = -1;
		local_ref[5][1][2] = -1+step;

		local_ref[5][2][0] = 1;
		local_ref[5][2][1] = -1+step;
		local_ref[5][2][2] = -1+step;

		local_ref[5][3][0] = -1+2*step;
		local_ref[5][3][1] = -1+1*3*step/double(4);
		local_ref[5][3][2] = -1+1*3*step/double(4);

		local_ref[5][4][0] = -1+2*step;
		local_ref[5][4][1] = -1;
		local_ref[5][4][2] = -1+2*3*step/double(4);

		local_ref[5][5][0] = 1;
		local_ref[5][5][1] = -1;
		local_ref[5][5][2] = -1+2*step;

		local_ref[5][6][0] = 1;
		local_ref[5][6][1] = -1+2*step;
		local_ref[5][6][2] = -1+2*step;

		local_ref[5][7][0] = -1+2*step;
		local_ref[5][7][1] = -1+2*3*step/double(4);
		local_ref[5][7][2] = -1+2*3*step/double(4);

		//element 6
		local_ref[6][0][0] = -1;
		local_ref[6][0][1] = -1;
		local_ref[6][0][2] = -1+2*step;

		local_ref[6][1][0] = -1+step;
		local_ref[6][1][1] = -1;
		local_ref[6][1][2] = -1+2*3*step/double(4);

		local_ref[6][2][0] = -1+step;
		local_ref[6][2][1] = -1+2*3*step/double(4);
		local_ref[6][2][2] = -1+2*3*step/double(4);

		local_ref[6][3][0] = -1;
		local_ref[6][3][1] = -1+2*step;
		local_ref[6][3][2] = -1+2*step;

		local_ref[6][4][0] = -1;
		local_ref[6][4][1] = -1;
		local_ref[6][4][2] = 1;

		local_ref[6][5][0] = -1+step;
		local_ref[6][5][1] = -1;
		local_ref[6][5][2] = -1+3*3*step/double(4);

		local_ref[6][6][0] = -1+step;
		local_ref[6][6][1] = -1+3*3*step/double(4);
		local_ref[6][6][2] = -1+3*3*step/double(4);

		local_ref[6][7][0] = -1;
		local_ref[6][7][1] = 1;
		local_ref[6][7][2] = 1;

		//element 7
		local_ref[7][0][0] = -1+step;
		local_ref[7][0][1] = -1;
		local_ref[7][0][2] = -1+2*3*step/double(4);

		local_ref[7][1][0] = -1+2*step;
		local_ref[7][1][1] = -1;
		local_ref[7][1][2] = -1+2*3*step/double(4);

		local_ref[7][2][0] = -1+2*step;
		local_ref[7][2][1] = -1+2*3*step/double(4);
		local_ref[7][2][2] = -1+2*3*step/double(4);

		local_ref[7][3][0] = -1+step;
		local_ref[7][3][1] = -1+2*3*step/double(4);
		local_ref[7][3][2] = -1+2*3*step/double(4);

		local_ref[7][4][0] = -1+step;
		local_ref[7][4][1] = -1;
		local_ref[7][4][2] = -1+3*3*step/double(4);

		local_ref[7][5][0] = -1+2*step;
		local_ref[7][5][1] = -1;
		local_ref[7][5][2] = -1+3*3*step/double(4);

		local_ref[7][6][0] = -1+2*step;
		local_ref[7][6][1] = -1+3*3*step/double(4);
		local_ref[7][6][2] = -1+3*3*step/double(4);

		local_ref[7][7][0] = -1+step;
		local_ref[7][7][1] = -1+3*3*step/double(4);
		local_ref[7][7][2] = -1+3*3*step/double(4);

		//element 8
		local_ref[8][0][0] = -1+2*step;
		local_ref[8][0][1] = -1;
		local_ref[8][0][2] = -1+2*3*step/double(4);

		local_ref[8][1][0] =  1;
		local_ref[8][1][1] = -1;
		local_ref[8][1][2] = -1+2*step;

		local_ref[8][2][0] = 1;
		local_ref[8][2][1] = -1+2*step;
		local_ref[8][2][2] = -1+2*step;

		local_ref[8][3][0] = -1+2*step;
		local_ref[8][3][1] = -1+2*3*step/double(4);
		local_ref[8][3][2] = -1+2*3*step/double(4);

		local_ref[8][4][0] = -1+2*step;
		local_ref[8][4][1] = -1;
		local_ref[8][4][2] = -1+3*3*step/double(4);

		local_ref[8][5][0] = 1;
		local_ref[8][5][1] = -1;
		local_ref[8][5][2] = 1;

		local_ref[8][6][0] = 1;
		local_ref[8][6][1] = 1;
		local_ref[8][6][2] = 1;

		local_ref[8][7][0] = -1+2*step;
		local_ref[8][7][1] = -1+3*3*step/double(4);
		local_ref[8][7][2] = -1+3*3*step/double(4);

		//element 9
		local_ref[9][0][0] = -1;
		local_ref[9][0][1] = -1+step;
		local_ref[9][0][2] = -1;

		local_ref[9][1][0] = -1+step;
		local_ref[9][1][1] = -1+1*3*step/double(4);
		local_ref[9][1][2] = -1;

		local_ref[9][2][0] = -1+step;
		local_ref[9][2][1] = -1+2*3*step/double(4);
		local_ref[9][2][2] = -1;

		local_ref[9][3][0] = -1;
		local_ref[9][3][1] = -1+2*step;;
		local_ref[9][3][2] = -1;

		local_ref[9][4][0] = -1;
		local_ref[9][4][1] = -1+step;
		local_ref[9][4][2] = -1+step;

		local_ref[9][5][0] = -1+step;
		local_ref[9][5][1] = -1+1*3*step/double(4);
		local_ref[9][5][2] = -1+1*3*step/double(4);

		local_ref[9][6][0] = -1+step;
		local_ref[9][6][1] = -1+2*3*step/double(4);
		local_ref[9][6][2] = -1+2*3*step/double(4);

		local_ref[9][7][0] = -1;
		local_ref[9][7][1] = -1+2*step;
		local_ref[9][7][2] = -1+2*step;

		//element 10
		local_ref[10][0][0] = -1+step;;
		local_ref[10][0][1] = -1+1*3*step/double(4);
		local_ref[10][0][2] = -1;

		local_ref[10][1][0] = -1+2*step;
		local_ref[10][1][1] = -1+1*3*step/double(4);
		local_ref[10][1][2] = -1;

		local_ref[10][2][0] = -1+2*step;
		local_ref[10][2][1] = -1+2*3*step/double(4);
		local_ref[10][2][2] = -1;

		local_ref[10][3][0] = -1+step;
		local_ref[10][3][1] = -1+2*3*step/double(4);
		local_ref[10][3][2] = -1;

		local_ref[10][4][0] = -1+step;
		local_ref[10][4][1] = -1+1*3*step/double(4);
		local_ref[10][4][2] = -1+1*3*step/double(4);

		local_ref[10][5][0] = -1+2*step;
		local_ref[10][5][1] = -1+1*3*step/double(4);
		local_ref[10][5][2] = -1+1*3*step/double(4);

		local_ref[10][6][0] = -1+2*step;
		local_ref[10][6][1] = -1+2*3*step/double(4);
		local_ref[10][6][2] = -1+2*3*step/double(4);

		local_ref[10][7][0] = -1+step;
		local_ref[10][7][1] = -1+2*3*step/double(4);
		local_ref[10][7][2] = -1+2*3*step/double(4);

		//element 11
		local_ref[11][0][0] = -1+2*step;
		local_ref[11][0][1] = -1+1*3*step/double(4);
		local_ref[11][0][2] = -1;

		local_ref[11][1][0] =  1;
		local_ref[11][1][1] = -1+step;
		local_ref[11][1][2] = -1;

		local_ref[11][2][0] = 1;
		local_ref[11][2][1] = -1+2*step;
		local_ref[11][2][2] = -1;

		local_ref[11][3][0] = -1+2*step;
		local_ref[11][3][1] = -1+2*3*step/double(4);
		local_ref[11][3][2] = -1;

		local_ref[11][4][0] =  -1+2*step;
		local_ref[11][4][1] = -1+1*3*step/double(4);
		local_ref[11][4][2] = -1+1*3*step/double(4);

		local_ref[11][5][0] =  1;
		local_ref[11][5][1] = -1+step;
		local_ref[11][5][2] = -1+step;

		local_ref[11][6][0] = 1;
		local_ref[11][6][1] = -1+2*step;
		local_ref[11][6][2] = -1+2*step;

		local_ref[11][7][0] = -1+2*step;
		local_ref[11][7][1] = -1+2*3*step/double(4);
		local_ref[11][7][2] = -1+2*3*step/double(4);

		//element 12
		local_ref[12][0][0] = -1;
		local_ref[12][0][1] = -1+2*step;
		local_ref[12][0][2] = -1;

		local_ref[12][1][0] = -1+step;
		local_ref[12][1][1] = -1+2*3*step/double(4);
		local_ref[12][1][2] = -1;

		local_ref[12][2][0] = -1+step;
		local_ref[12][2][1] = -1+3*3*step/double(4);
		local_ref[12][2][2] = -1;

		local_ref[12][3][0] = -1;
		local_ref[12][3][1] = 1;
		local_ref[12][3][2] = -1;

		local_ref[12][4][0] = -1;
		local_ref[12][4][1] = -1+2*step;
		local_ref[12][4][2] = -1+2*step;

		local_ref[12][5][0] = -1+step;
		local_ref[12][5][1] = -1+2*3*step/double(4);
		local_ref[12][5][2] = -1+2*3*step/double(4);

		local_ref[12][6][0] = -1+step;
		local_ref[12][6][1] = -1+3*3*step/double(4);
		local_ref[12][6][2] = -1+3*3*step/double(4);

		local_ref[12][7][0] = -1;
		local_ref[12][7][1] = 1;
		local_ref[12][7][2] = 1;

		//element 13
		local_ref[13][0][0] = -1+step;;
		local_ref[13][0][1] = -1+2*3*step/double(4);
		local_ref[13][0][2] = -1;

		local_ref[13][1][0] = -1+2*step;
		local_ref[13][1][1] = -1+2*3*step/double(4);
		local_ref[13][1][2] = -1;

		local_ref[13][2][0] = -1+2*step;
		local_ref[13][2][1] = -1+3*3*step/double(4);
		local_ref[13][2][2] = -1;

		local_ref[13][3][0] = -1+step;
		local_ref[13][3][1] = -1+3*3*step/double(4);
		local_ref[13][3][2] = -1;

		local_ref[13][4][0] = -1+step;
		local_ref[13][4][1] = -1+2*3*step/double(4);
		local_ref[13][4][2] = -1+2*3*step/double(4);

		local_ref[13][5][0] = -1+2*step;
		local_ref[13][5][1] = -1+2*3*step/double(4);
		local_ref[13][5][2] = -1+2*3*step/double(4);

		local_ref[13][6][0] = -1+2*step;
		local_ref[13][6][1] = -1+3*3*step/double(4);
		local_ref[13][6][2] = -1+3*3*step/double(4);

		local_ref[13][7][0] = -1+step;
		local_ref[13][7][1] = -1+3*3*step/double(4);
		local_ref[13][7][2] = -1+3*3*step/double(4);

		//element 14
		local_ref[14][0][0] = -1+2*step;
		local_ref[14][0][1] = -1+2*3*step/double(4);
		local_ref[14][0][2] = -1;

		local_ref[14][1][0] =  1;
		local_ref[14][1][1] = -1+2*step;
		local_ref[14][1][2] = -1;

		local_ref[14][2][0] = 1;
		local_ref[14][2][1] = 1;
		local_ref[14][2][2] = -1;

		local_ref[14][3][0] = -1+2*step;
		local_ref[14][3][1] = -1+3*3*step/double(4);
		local_ref[14][3][2] = -1;

		local_ref[14][4][0] =  -1+2*step;
		local_ref[14][4][1] = -1+2*3*step/double(4);
		local_ref[14][4][2] = -1+2*3*step/double(4);

		local_ref[14][5][0] =  1;
		local_ref[14][5][1] = -1+2*step;
		local_ref[14][5][2] = -1+2*step;

		local_ref[14][6][0] = 1;
		local_ref[14][6][1] = 1;
		local_ref[14][6][2] = 1;

		local_ref[14][7][0] = -1+2*step;
		local_ref[14][7][1] = -1+3*3*step/double(4);
		local_ref[14][7][2] = -1+3*3*step/double(4);

		//element 15
		local_ref[15][0][0] = -1+step;;
		local_ref[15][0][1] = -1;
		local_ref[15][0][2] = -1+3*3*step/double(4);

		local_ref[15][1][0] = -1+2*step;
		local_ref[15][1][1] = -1;
		local_ref[15][1][2] = -1+3*3*step/double(4);

		local_ref[15][2][0] = -1+2*step;
		local_ref[15][2][1] = -1+3*3*step/double(4);
		local_ref[15][2][2] = -1+3*3*step/double(4);

		local_ref[15][3][0] = -1+step;
		local_ref[15][3][1] = -1+3*3*step/double(4);
		local_ref[15][3][2] = -1+3*3*step/double(4);

		local_ref[15][4][0] = -1;
		local_ref[15][4][1] = -1;
		local_ref[15][4][2] = 1;

		local_ref[15][5][0] = 1;
		local_ref[15][5][1] = -1;
		local_ref[15][5][2] = 1;

		local_ref[15][6][0] = 1;
		local_ref[15][6][1] = 1;
		local_ref[15][6][2] = 1;

		local_ref[15][7][0] = -1;
		local_ref[15][7][1] = 1;
		local_ref[15][7][2] = 1;

		//element 16
		local_ref[16][0][0] = -1+1*step;
		local_ref[16][0][1] = -1+3*3*step/double(4);
		local_ref[16][0][2] = -1+3*3*step/double(4);

		local_ref[16][1][0] =  -1+2*step;
		local_ref[16][1][1] = -1+3*3*step/double(4);
		local_ref[16][1][2] = -1+3*3*step/double(4);

		local_ref[16][2][0] = 1;
		local_ref[16][2][1] = 1;
		local_ref[16][2][2] = 1;

		local_ref[16][3][0] = -1;
		local_ref[16][3][1] = 1;
		local_ref[16][3][2] = 1;

		local_ref[16][4][0] = -1+1*step;
		local_ref[16][4][1] = -1+3*3*step/double(4);
		local_ref[16][4][2] = -1;

		local_ref[16][5][0] =  -1+2*step;
		local_ref[16][5][1] = -1+3*3*step/double(4);
		local_ref[16][5][2] = -1;

		local_ref[16][6][0] = 1;
		local_ref[16][6][1] = 1;
		local_ref[16][6][2] = -1;

		local_ref[16][7][0] = -1;
		local_ref[16][7][1] = 1;
		local_ref[16][7][2] = -1;
	}

	//define the rotation of the reference element
	int rot[3];
	int sym[3];
	int id_node[8];
	vector<int> ord;

	if(elem->pad==132){
		//edge 3 0 1 4 5
		rot[0] = 0;
		rot[1] = 0;
		rot[2] = 0;

		sym[0] = 0;
		sym[1] = 0;
		sym[2] = 0;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}

	}else if(elem->pad==133){
		//edge 0 1 2 5 6
		rot[0] = 0;
		rot[1] = 0;
		rot[2] = -1;

		sym[0] = 0;
		sym[1] = 0;
		sym[2] = 0;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}

	}else if(elem->pad==134){
		//edge 1 2 3 6 7
		rot[0] = 0;
		rot[1] = 0;
		rot[2] = 0;

		sym[0] = 0;
		sym[1] = 1;
		sym[2] = 0;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}

	}else if(elem->pad==135){
		//edge 0 2 3 4 7
		rot[0] = 0;
		rot[1] = 0;
		rot[2] = 1;

		sym[0] = 0;
		sym[1] = 0;
		sym[2] = 0;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}

	}else if(elem->pad==136){
		//edge 0 3 4 8 11
		rot[0] = 0;
		rot[1] = -1;
		rot[2] = 0;

		sym[0] = 0;
		sym[1] = 0;
		sym[2] = 0;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}

	}else if(elem->pad==137){
		//edge 0 1 5 8 9
		rot[0] = 0;
		rot[1] = 1;
		rot[2] = 0;

		sym[0] = 0;
		sym[1] = 0;
		sym[2] = 0;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}

	}else if(elem->pad==138){
		//edge 1 2 6 9 10
		rot[0] = 0;
		rot[1] = 1;
		rot[2] = 0;

		sym[0] = 0;
		sym[1] = 1;
		sym[2] = 0;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}

	}else if(elem->pad==139){
		//edge 2 3 7 10 11
		rot[0] = 0;
		rot[1] = -1;
		rot[2] = 0;

		sym[0] = 0;
		sym[1] = 1;
		sym[2] = 0;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}

	}else if(elem->pad==140){
		//edge 8 9 11 4 5
		rot[0] = 0;
		rot[1] = 0;
		rot[2] = 0;

		sym[0] = 0;
		sym[1] = 0;
		sym[2] = 1;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}

	}else if(elem->pad==141){
		//edge 8 9 10 5 6
		rot[0] = 0;
		rot[1] = 0;
		rot[2] = -1;

		sym[0] = 0;
		sym[1] = 0;
		sym[2] = 1;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}

	}else if(elem->pad==142){
		//edge 9 10 11 4 7
		rot[0] = 0;
		rot[1] = 0;
		rot[2] = 0;

		sym[0] = 0;
		sym[1] = 1;
		sym[2] = 1;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}

	}else if(elem->pad==143){
		//edge 8 10 11 4 7
		rot[0] = 0;
		rot[1] = 0;
		rot[2] = 1;

		sym[0] = 0;
		sym[1] = 0;
		sym[2] = 1;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}

	}else{
		printf ("Error in template 10\n");
		exit (EXIT_FAILURE);
	}


	for(int iel = 0; iel < 17; iel++){
		double local_ref_x[8];
		double local_ref_y[8];
		double local_ref_z[8];
		for(int ino = 0; ino < 8; ino++){
			local_ref_x[ino] = local_ref[iel][ino][0];
			local_ref_y[ino] = local_ref[iel][ino][1];
			local_ref_z[ino] = local_ref[iel][ino][2];
		}
		ApplyElement(mesh, coords,  id,  iel,  id_node,  local_ref_x, local_ref_y, local_ref_z,  ord,  hash_nodes);
	}
	sc_array_reset(&toto);
}

void ApplyTemplate11(hexa_tree_t* mesh, std::vector<double>& coords, std::vector<int>& elements_ids, int iel, sc_hash_array_t* hash_nodes){

	sc_array_t toto;
	sc_array_init(&toto, sizeof(octant_t));

	octant_t* elemOrig = (octant_t*) sc_array_index(&mesh->elements, elements_ids[iel]);
	octant_t * elem = (octant_t*) sc_array_push(&toto);
	hexa_element_copy(elemOrig,elem);

	double step = double(2)/double(3);
	int id = elements_ids[iel];
	double cord_in_ref[3];
	cord_in_ref[0] = 0;
	cord_in_ref[1] = 0;
	cord_in_ref[2] = 0;

	//reference element edge 0
	double local_ref[27][8][3];
	if(true){
		//element 0
		local_ref[0][0][0] = -1+0*step;
		local_ref[0][0][1] = -1+0*step;
		local_ref[0][0][2] = -1+0*step;

		local_ref[0][1][0] = -1+1*step;
		local_ref[0][1][1] = -1+0*step;
		local_ref[0][1][2] = -1+0*step;

		local_ref[0][2][0] = -1+1*step;
		local_ref[0][2][1] = -1+1*step;
		local_ref[0][2][2] = -1+0*step;

		local_ref[0][3][0] = -1+0*step;
		local_ref[0][3][1] = -1+1*step;
		local_ref[0][3][2] = -1+0*step;

		local_ref[0][4][0] = -1+0*step;
		local_ref[0][4][1] = -1+0*step;
		local_ref[0][4][2] = -1+1*step;

		local_ref[0][5][0] = -1+1*step;
		local_ref[0][5][1] = -1+0*step;
		local_ref[0][5][2] = -1+1*step;

		local_ref[0][6][0] = -1+1*step;
		local_ref[0][6][1] = -1+1*step;
		local_ref[0][6][2] = -1+1*step;

		local_ref[0][7][0] = -1+0*step;
		local_ref[0][7][1] = -1+1*step;
		local_ref[0][7][2] = -1+1*step;

		//element 1
		local_ref[1][0][0] = -1+1*step;
		local_ref[1][0][1] = -1+0*step;
		local_ref[1][0][2] = -1+0*step;

		local_ref[1][1][0] = -1+2*step;
		local_ref[1][1][1] = -1+0*step;
		local_ref[1][1][2] = -1+0*step;

		local_ref[1][2][0] = -1+2*step;
		local_ref[1][2][1] = -1+1*step;
		local_ref[1][2][2] = -1+0*step;

		local_ref[1][3][0] = -1+1*step;
		local_ref[1][3][1] = -1+1*step;
		local_ref[1][3][2] = -1+0*step;

		local_ref[1][4][0] = -1+1*step;
		local_ref[1][4][1] = -1+0*step;
		local_ref[1][4][2] = -1+1*step;

		local_ref[1][5][0] = -1+2*step;
		local_ref[1][5][1] = -1+0*step;
		local_ref[1][5][2] = -1+1*step;

		local_ref[1][6][0] = -1+2*step;
		local_ref[1][6][1] = -1+1*step;
		local_ref[1][6][2] = -1+1*step;

		local_ref[1][7][0] = -1+1*step;
		local_ref[1][7][1] = -1+1*step;
		local_ref[1][7][2] = -1+1*step;

		//element 2
		local_ref[2][0][0] = -1+2*step;
		local_ref[2][0][1] = -1+0*step;
		local_ref[2][0][2] = -1+0*step;

		local_ref[2][1][0] = -1+3*step;
		local_ref[2][1][1] = -1+0*step;
		local_ref[2][1][2] = -1+0*step;

		local_ref[2][2][0] = -1+3*step;
		local_ref[2][2][1] = -1+1*step;
		local_ref[2][2][2] = -1+0*step;

		local_ref[2][3][0] = -1+2*step;
		local_ref[2][3][1] = -1+1*step;
		local_ref[2][3][2] = -1+0*step;

		local_ref[2][4][0] = -1+2*step;
		local_ref[2][4][1] = -1+0*step;
		local_ref[2][4][2] = -1+1*step;

		local_ref[2][5][0] = -1+3*step;
		local_ref[2][5][1] = -1+0*step;
		local_ref[2][5][2] = -1+1*step;

		local_ref[2][6][0] = -1+3*step;
		local_ref[2][6][1] = -1+1*step;
		local_ref[2][6][2] = -1+1*step;

		local_ref[2][7][0] = -1+2*step;
		local_ref[2][7][1] = -1+1*step;
		local_ref[2][7][2] = -1+1*step;

		//element 3
		local_ref[3][0][0] = -1+0*step;
		local_ref[3][0][1] = -1+1*step;
		local_ref[3][0][2] = -1+0*step;

		local_ref[3][1][0] = -1+1*step;
		local_ref[3][1][1] = -1+1*step;
		local_ref[3][1][2] = -1+0*step;

		local_ref[3][2][0] = -1+1*step;
		local_ref[3][2][1] = -1+2*step;
		local_ref[3][2][2] = -1+0*step;

		local_ref[3][3][0] = -1+0*step;
		local_ref[3][3][1] = -1+2*step;
		local_ref[3][3][2] = -1+0*step;

		local_ref[3][4][0] = -1+0*step;
		local_ref[3][4][1] = -1+1*step;
		local_ref[3][4][2] = -1+1*step;

		local_ref[3][5][0] = -1+1*step;
		local_ref[3][5][1] = -1+1*step;
		local_ref[3][5][2] = -1+1*step;

		local_ref[3][6][0] = -1+1*step;
		local_ref[3][6][1] = -1+2*step;
		local_ref[3][6][2] = -1+1*step;

		local_ref[3][7][0] = -1+0*step;
		local_ref[3][7][1] = -1+2*step;
		local_ref[3][7][2] = -1+1*step;

		//element 4
		local_ref[4][0][0] = -1+1*step;
		local_ref[4][0][1] = -1+1*step;
		local_ref[4][0][2] = -1+0*step;

		local_ref[4][1][0] = -1+2*step;
		local_ref[4][1][1] = -1+1*step;
		local_ref[4][1][2] = -1+0*step;

		local_ref[4][2][0] = -1+2*step;
		local_ref[4][2][1] = -1+2*step;
		local_ref[4][2][2] = -1+0*step;

		local_ref[4][3][0] = -1+1*step;
		local_ref[4][3][1] = -1+2*step;
		local_ref[4][3][2] = -1+0*step;

		local_ref[4][4][0] = -1+1*step;
		local_ref[4][4][1] = -1+1*step;
		local_ref[4][4][2] = -1+1*step;

		local_ref[4][5][0] = -1+2*step;
		local_ref[4][5][1] = -1+1*step;
		local_ref[4][5][2] = -1+1*step;

		local_ref[4][6][0] = -1+2*step;
		local_ref[4][6][1] = -1+2*step;
		local_ref[4][6][2] = -1+1*step;

		local_ref[4][7][0] = -1+1*step;
		local_ref[4][7][1] = -1+2*step;
		local_ref[4][7][2] = -1+1*step;

		//element 5
		local_ref[5][0][0] = -1+2*step;
		local_ref[5][0][1] = -1+1*step;
		local_ref[5][0][2] = -1+0*step;

		local_ref[5][1][0] = -1+3*step;
		local_ref[5][1][1] = -1+1*step;
		local_ref[5][1][2] = -1+0*step;

		local_ref[5][2][0] = -1+3*step;
		local_ref[5][2][1] = -1+2*step;
		local_ref[5][2][2] = -1+0*step;

		local_ref[5][3][0] = -1+2*step;
		local_ref[5][3][1] = -1+2*step;
		local_ref[5][3][2] = -1+0*step;

		local_ref[5][4][0] = -1+2*step;
		local_ref[5][4][1] = -1+1*step;
		local_ref[5][4][2] = -1+1*step;

		local_ref[5][5][0] = -1+3*step;
		local_ref[5][5][1] = -1+1*step;
		local_ref[5][5][2] = -1+1*step;

		local_ref[5][6][0] = -1+3*step;
		local_ref[5][6][1] = -1+2*step;
		local_ref[5][6][2] = -1+1*step;

		local_ref[5][7][0] = -1+2*step;
		local_ref[5][7][1] = -1+2*step;
		local_ref[5][7][2] = -1+1*step;

		//element 6
		local_ref[6][0][0] = -1+0*step;
		local_ref[6][0][1] = -1+2*step;
		local_ref[6][0][2] = -1+0*step;

		local_ref[6][1][0] = -1+1*step;
		local_ref[6][1][1] = -1+2*step;
		local_ref[6][1][2] = -1+0*step;

		local_ref[6][2][0] = -1+1*step;
		local_ref[6][2][1] = -1+3*step;
		local_ref[6][2][2] = -1+0*step;

		local_ref[6][3][0] = -1+0*step;
		local_ref[6][3][1] = -1+3*step;
		local_ref[6][3][2] = -1+0*step;

		local_ref[6][4][0] = -1+0*step;
		local_ref[6][4][1] = -1+2*step;
		local_ref[6][4][2] = -1+1*step;

		local_ref[6][5][0] = -1+1*step;
		local_ref[6][5][1] = -1+2*step;
		local_ref[6][5][2] = -1+1*step;

		local_ref[6][6][0] = -1+1*step;
		local_ref[6][6][1] = -1+3*step;
		local_ref[6][6][2] = -1+1*step;

		local_ref[6][7][0] = -1+0*step;
		local_ref[6][7][1] = -1+3*step;
		local_ref[6][7][2] = -1+1*step;

		//element 7
		local_ref[7][0][0] = -1+1*step;
		local_ref[7][0][1] = -1+2*step;
		local_ref[7][0][2] = -1+0*step;

		local_ref[7][1][0] = -1+2*step;
		local_ref[7][1][1] = -1+2*step;
		local_ref[7][1][2] = -1+0*step;

		local_ref[7][2][0] = -1+2*step;
		local_ref[7][2][1] = -1+3*step;
		local_ref[7][2][2] = -1+0*step;

		local_ref[7][3][0] = -1+1*step;
		local_ref[7][3][1] = -1+3*step;
		local_ref[7][3][2] = -1+0*step;

		local_ref[7][4][0] = -1+1*step;
		local_ref[7][4][1] = -1+2*step;
		local_ref[7][4][2] = -1+1*step;

		local_ref[7][5][0] = -1+2*step;
		local_ref[7][5][1] = -1+2*step;
		local_ref[7][5][2] = -1+1*step;

		local_ref[7][6][0] = -1+2*step;
		local_ref[7][6][1] = -1+3*step;
		local_ref[7][6][2] = -1+1*step;

		local_ref[7][7][0] = -1+1*step;
		local_ref[7][7][1] = -1+3*step;
		local_ref[7][7][2] = -1+1*step;

		//element 8
		local_ref[8][0][0] = -1+2*step;
		local_ref[8][0][1] = -1+2*step;
		local_ref[8][0][2] = -1+0*step;

		local_ref[8][1][0] = -1+3*step;
		local_ref[8][1][1] = -1+2*step;
		local_ref[8][1][2] = -1+0*step;

		local_ref[8][2][0] = -1+3*step;
		local_ref[8][2][1] = -1+3*step;
		local_ref[8][2][2] = -1+0*step;

		local_ref[8][3][0] = -1+2*step;
		local_ref[8][3][1] = -1+3*step;
		local_ref[8][3][2] = -1+0*step;

		local_ref[8][4][0] = -1+2*step;
		local_ref[8][4][1] = -1+2*step;
		local_ref[8][4][2] = -1+1*step;

		local_ref[8][5][0] = -1+3*step;
		local_ref[8][5][1] = -1+2*step;
		local_ref[8][5][2] = -1+1*step;

		local_ref[8][6][0] = -1+3*step;
		local_ref[8][6][1] = -1+3*step;
		local_ref[8][6][2] = -1+1*step;

		local_ref[8][7][0] = -1+2*step;
		local_ref[8][7][1] = -1+3*step;
		local_ref[8][7][2] = -1+1*step;

		//element 9
		local_ref[9][0][0] = -1+0*step;
		local_ref[9][0][1] = -1+0*step;
		local_ref[9][0][2] = -1+1*step;

		local_ref[9][1][0] = -1+1*step;
		local_ref[9][1][1] = -1+0*step;
		local_ref[9][1][2] = -1+1*step;

		local_ref[9][2][0] = -1+1*step;
		local_ref[9][2][1] = -1+1*step;
		local_ref[9][2][2] = -1+1*step;

		local_ref[9][3][0] = -1+0*step;
		local_ref[9][3][1] = -1+1*step;
		local_ref[9][3][2] = -1+1*step;

		local_ref[9][4][0] = -1+0*step;
		local_ref[9][4][1] = -1+0*step;
		local_ref[9][4][2] = -1+2*step;

		local_ref[9][5][0] = -1+1*step;
		local_ref[9][5][1] = -1+0*step;
		local_ref[9][5][2] = -1+2*step;

		local_ref[9][6][0] = -1+1*step;
		local_ref[9][6][1] = -1+1*step;
		local_ref[9][6][2] = -1+2*step;

		local_ref[9][7][0] = -1+0*step;
		local_ref[9][7][1] = -1+1*step;
		local_ref[9][7][2] = -1+2*step;

		//element 10
		local_ref[10][0][0] = -1+1*step;
		local_ref[10][0][1] = -1+0*step;
		local_ref[10][0][2] = -1+1*step;

		local_ref[10][1][0] = -1+2*step;
		local_ref[10][1][1] = -1+0*step;
		local_ref[10][1][2] = -1+1*step;

		local_ref[10][2][0] = -1+2*step;
		local_ref[10][2][1] = -1+1*step;
		local_ref[10][2][2] = -1+1*step;

		local_ref[10][3][0] = -1+1*step;
		local_ref[10][3][1] = -1+1*step;
		local_ref[10][3][2] = -1+1*step;

		local_ref[10][4][0] = -1+1*step;
		local_ref[10][4][1] = -1+0*step;
		local_ref[10][4][2] = -1+2*step;

		local_ref[10][5][0] = -1+2*step;
		local_ref[10][5][1] = -1+0*step;
		local_ref[10][5][2] = -1+2*step;

		local_ref[10][6][0] = -1+2*step;
		local_ref[10][6][1] = -1+1*step;
		local_ref[10][6][2] = -1+2*step;

		local_ref[10][7][0] = -1+1*step;
		local_ref[10][7][1] = -1+1*step;
		local_ref[10][7][2] = -1+2*step;

		//element 11
		local_ref[11][0][0] = -1+2*step;
		local_ref[11][0][1] = -1+0*step;
		local_ref[11][0][2] = -1+1*step;

		local_ref[11][1][0] = -1+3*step;
		local_ref[11][1][1] = -1+0*step;
		local_ref[11][1][2] = -1+1*step;

		local_ref[11][2][0] = -1+3*step;
		local_ref[11][2][1] = -1+1*step;
		local_ref[11][2][2] = -1+1*step;

		local_ref[11][3][0] = -1+2*step;
		local_ref[11][3][1] = -1+1*step;
		local_ref[11][3][2] = -1+1*step;

		local_ref[11][4][0] = -1+2*step;
		local_ref[11][4][1] = -1+0*step;
		local_ref[11][4][2] = -1+2*step;

		local_ref[11][5][0] = -1+3*step;
		local_ref[11][5][1] = -1+0*step;
		local_ref[11][5][2] = -1+2*step;

		local_ref[11][6][0] = -1+3*step;
		local_ref[11][6][1] = -1+1*step;
		local_ref[11][6][2] = -1+2*step;

		local_ref[11][7][0] = -1+2*step;
		local_ref[11][7][1] = -1+1*step;
		local_ref[11][7][2] = -1+2*step;

		//element 12
		local_ref[12][0][0] = -1+0*step;
		local_ref[12][0][1] = -1+1*step;
		local_ref[12][0][2] = -1+1*step;

		local_ref[12][1][0] = -1+1*step;
		local_ref[12][1][1] = -1+1*step;
		local_ref[12][1][2] = -1+1*step;

		local_ref[12][2][0] = -1+1*step;
		local_ref[12][2][1] = -1+2*step;
		local_ref[12][2][2] = -1+1*step;

		local_ref[12][3][0] = -1+0*step;
		local_ref[12][3][1] = -1+2*step;
		local_ref[12][3][2] = -1+1*step;

		local_ref[12][4][0] = -1+0*step;
		local_ref[12][4][1] = -1+1*step;
		local_ref[12][4][2] = -1+2*step;

		local_ref[12][5][0] = -1+1*step;
		local_ref[12][5][1] = -1+1*step;
		local_ref[12][5][2] = -1+2*step;

		local_ref[12][6][0] = -1+1*step;
		local_ref[12][6][1] = -1+2*step;
		local_ref[12][6][2] = -1+2*step;

		local_ref[12][7][0] = -1+0*step;
		local_ref[12][7][1] = -1+2*step;
		local_ref[12][7][2] = -1+2*step;

		//element 13
		local_ref[13][0][0] = -1+1*step;
		local_ref[13][0][1] = -1+1*step;
		local_ref[13][0][2] = -1+1*step;

		local_ref[13][1][0] = -1+2*step;
		local_ref[13][1][1] = -1+1*step;
		local_ref[13][1][2] = -1+1*step;

		local_ref[13][2][0] = -1+2*step;
		local_ref[13][2][1] = -1+2*step;
		local_ref[13][2][2] = -1+1*step;

		local_ref[13][3][0] = -1+1*step;
		local_ref[13][3][1] = -1+2*step;
		local_ref[13][3][2] = -1+1*step;

		local_ref[13][4][0] = -1+1*step;
		local_ref[13][4][1] = -1+1*step;
		local_ref[13][4][2] = -1+2*step;

		local_ref[13][5][0] = -1+2*step;
		local_ref[13][5][1] = -1+1*step;
		local_ref[13][5][2] = -1+2*step;

		local_ref[13][6][0] = -1+2*step;
		local_ref[13][6][1] = -1+2*step;
		local_ref[13][6][2] = -1+2*step;

		local_ref[13][7][0] = -1+1*step;
		local_ref[13][7][1] = -1+2*step;
		local_ref[13][7][2] = -1+2*step;

		//element 14
		local_ref[14][0][0] = -1+2*step;
		local_ref[14][0][1] = -1+1*step;
		local_ref[14][0][2] = -1+1*step;

		local_ref[14][1][0] = -1+3*step;
		local_ref[14][1][1] = -1+1*step;
		local_ref[14][1][2] = -1+1*step;

		local_ref[14][2][0] = -1+3*step;
		local_ref[14][2][1] = -1+2*step;
		local_ref[14][2][2] = -1+1*step;

		local_ref[14][3][0] = -1+2*step;
		local_ref[14][3][1] = -1+2*step;
		local_ref[14][3][2] = -1+1*step;

		local_ref[14][4][0] = -1+2*step;
		local_ref[14][4][1] = -1+1*step;
		local_ref[14][4][2] = -1+2*step;

		local_ref[14][5][0] = -1+3*step;
		local_ref[14][5][1] = -1+1*step;
		local_ref[14][5][2] = -1+2*step;

		local_ref[14][6][0] = -1+3*step;
		local_ref[14][6][1] = -1+2*step;
		local_ref[14][6][2] = -1+2*step;

		local_ref[14][7][0] = -1+2*step;
		local_ref[14][7][1] = -1+2*step;
		local_ref[14][7][2] = -1+2*step;

		//element 15
		local_ref[15][0][0] = -1+0*step;
		local_ref[15][0][1] = -1+2*step;
		local_ref[15][0][2] = -1+1*step;

		local_ref[15][1][0] = -1+1*step;
		local_ref[15][1][1] = -1+2*step;
		local_ref[15][1][2] = -1+1*step;

		local_ref[15][2][0] = -1+1*step;
		local_ref[15][2][1] = -1+3*step;
		local_ref[15][2][2] = -1+1*step;

		local_ref[15][3][0] = -1+0*step;
		local_ref[15][3][1] = -1+3*step;
		local_ref[15][3][2] = -1+1*step;

		local_ref[15][4][0] = -1+0*step;
		local_ref[15][4][1] = -1+2*step;
		local_ref[15][4][2] = -1+2*step;

		local_ref[15][5][0] = -1+1*step;
		local_ref[15][5][1] = -1+2*step;
		local_ref[15][5][2] = -1+2*step;

		local_ref[15][6][0] = -1+1*step;
		local_ref[15][6][1] = -1+3*step;
		local_ref[15][6][2] = -1+2*step;

		local_ref[15][7][0] = -1+0*step;
		local_ref[15][7][1] = -1+3*step;
		local_ref[15][7][2] = -1+2*step;

		//element 16
		local_ref[16][0][0] = -1+1*step;
		local_ref[16][0][1] = -1+2*step;
		local_ref[16][0][2] = -1+1*step;

		local_ref[16][1][0] = -1+2*step;
		local_ref[16][1][1] = -1+2*step;
		local_ref[16][1][2] = -1+1*step;

		local_ref[16][2][0] = -1+2*step;
		local_ref[16][2][1] = -1+3*step;
		local_ref[16][2][2] = -1+1*step;

		local_ref[16][3][0] = -1+1*step;
		local_ref[16][3][1] = -1+3*step;
		local_ref[16][3][2] = -1+1*step;

		local_ref[16][4][0] = -1+1*step;
		local_ref[16][4][1] = -1+2*step;
		local_ref[16][4][2] = -1+2*step;

		local_ref[16][5][0] = -1+2*step;
		local_ref[16][5][1] = -1+2*step;
		local_ref[16][5][2] = -1+2*step;

		local_ref[16][6][0] = -1+2*step;
		local_ref[16][6][1] = -1+3*step;
		local_ref[16][6][2] = -1+2*step;

		local_ref[16][7][0] = -1+1*step;
		local_ref[16][7][1] = -1+3*step;
		local_ref[16][7][2] = -1+2*step;

		//element 17
		local_ref[17][0][0] = -1+2*step;
		local_ref[17][0][1] = -1+2*step;
		local_ref[17][0][2] = -1+1*step;

		local_ref[17][1][0] = -1+3*step;
		local_ref[17][1][1] = -1+2*step;
		local_ref[17][1][2] = -1+1*step;

		local_ref[17][2][0] = -1+3*step;
		local_ref[17][2][1] = -1+3*step;
		local_ref[17][2][2] = -1+1*step;

		local_ref[17][3][0] = -1+2*step;
		local_ref[17][3][1] = -1+3*step;
		local_ref[17][3][2] = -1+1*step;

		local_ref[17][4][0] = -1+2*step;
		local_ref[17][4][1] = -1+2*step;
		local_ref[17][4][2] = -1+2*step;

		local_ref[17][5][0] = -1+3*step;
		local_ref[17][5][1] = -1+2*step;
		local_ref[17][5][2] = -1+2*step;

		local_ref[17][6][0] = -1+3*step;
		local_ref[17][6][1] = -1+3*step;
		local_ref[17][6][2] = -1+2*step;

		local_ref[17][7][0] = -1+2*step;
		local_ref[17][7][1] = -1+3*step;
		local_ref[17][7][2] = -1+2*step;

		//element 18
		local_ref[18][0][0] = -1+0*step;
		local_ref[18][0][1] = -1+0*step;
		local_ref[18][0][2] = -1+2*step;

		local_ref[18][1][0] = -1+1*step;
		local_ref[18][1][1] = -1+0*step;
		local_ref[18][1][2] = -1+2*step;

		local_ref[18][2][0] = -1+1*step;
		local_ref[18][2][1] = -1+1*step;
		local_ref[18][2][2] = -1+2*step;

		local_ref[18][3][0] = -1+0*step;
		local_ref[18][3][1] = -1+1*step;
		local_ref[18][3][2] = -1+2*step;

		local_ref[18][4][0] = -1+0*step;
		local_ref[18][4][1] = -1+0*step;
		local_ref[18][4][2] = -1+3*step;

		local_ref[18][5][0] = -1+1*step;
		local_ref[18][5][1] = -1+0*step;
		local_ref[18][5][2] = -1+3*step;

		local_ref[18][6][0] = -1+1*step;
		local_ref[18][6][1] = -1+1*step;
		local_ref[18][6][2] = -1+3*step;

		local_ref[18][7][0] = -1+0*step;
		local_ref[18][7][1] = -1+1*step;
		local_ref[18][7][2] = -1+3*step;

		//element 19
		local_ref[19][0][0] = -1+1*step;
		local_ref[19][0][1] = -1+0*step;
		local_ref[19][0][2] = -1+2*step;

		local_ref[19][1][0] = -1+2*step;
		local_ref[19][1][1] = -1+0*step;
		local_ref[19][1][2] = -1+2*step;

		local_ref[19][2][0] = -1+2*step;
		local_ref[19][2][1] = -1+1*step;
		local_ref[19][2][2] = -1+2*step;

		local_ref[19][3][0] = -1+1*step;
		local_ref[19][3][1] = -1+1*step;
		local_ref[19][3][2] = -1+2*step;

		local_ref[19][4][0] = -1+1*step;
		local_ref[19][4][1] = -1+0*step;
		local_ref[19][4][2] = -1+3*step;

		local_ref[19][5][0] = -1+2*step;
		local_ref[19][5][1] = -1+0*step;
		local_ref[19][5][2] = -1+3*step;

		local_ref[19][6][0] = -1+2*step;
		local_ref[19][6][1] = -1+1*step;
		local_ref[19][6][2] = -1+3*step;

		local_ref[19][7][0] = -1+1*step;
		local_ref[19][7][1] = -1+1*step;
		local_ref[19][7][2] = -1+3*step;

		//element 20
		local_ref[20][0][0] = -1+2*step;
		local_ref[20][0][1] = -1+0*step;
		local_ref[20][0][2] = -1+2*step;

		local_ref[20][1][0] = -1+3*step;
		local_ref[20][1][1] = -1+0*step;
		local_ref[20][1][2] = -1+2*step;

		local_ref[20][2][0] = -1+3*step;
		local_ref[20][2][1] = -1+1*step;
		local_ref[20][2][2] = -1+2*step;

		local_ref[20][3][0] = -1+2*step;
		local_ref[20][3][1] = -1+1*step;
		local_ref[20][3][2] = -1+2*step;

		local_ref[20][4][0] = -1+2*step;
		local_ref[20][4][1] = -1+0*step;
		local_ref[20][4][2] = -1+3*step;

		local_ref[20][5][0] = -1+3*step;
		local_ref[20][5][1] = -1+0*step;
		local_ref[20][5][2] = -1+3*step;

		local_ref[20][6][0] = -1+3*step;
		local_ref[20][6][1] = -1+1*step;
		local_ref[20][6][2] = -1+3*step;

		local_ref[20][7][0] = -1+2*step;
		local_ref[20][7][1] = -1+1*step;
		local_ref[20][7][2] = -1+3*step;

		//element 21
		local_ref[21][0][0] = -1+0*step;
		local_ref[21][0][1] = -1+1*step;
		local_ref[21][0][2] = -1+2*step;

		local_ref[21][1][0] = -1+1*step;
		local_ref[21][1][1] = -1+1*step;
		local_ref[21][1][2] = -1+2*step;

		local_ref[21][2][0] = -1+1*step;
		local_ref[21][2][1] = -1+2*step;
		local_ref[21][2][2] = -1+2*step;

		local_ref[21][3][0] = -1+0*step;
		local_ref[21][3][1] = -1+2*step;
		local_ref[21][3][2] = -1+2*step;

		local_ref[21][4][0] = -1+0*step;
		local_ref[21][4][1] = -1+1*step;
		local_ref[21][4][2] = -1+3*step;

		local_ref[21][5][0] = -1+1*step;
		local_ref[21][5][1] = -1+1*step;
		local_ref[21][5][2] = -1+3*step;

		local_ref[21][6][0] = -1+1*step;
		local_ref[21][6][1] = -1+2*step;
		local_ref[21][6][2] = -1+3*step;

		local_ref[21][7][0] = -1+0*step;
		local_ref[21][7][1] = -1+2*step;
		local_ref[21][7][2] = -1+3*step;

		//element 22
		local_ref[22][0][0] = -1+1*step;
		local_ref[22][0][1] = -1+1*step;
		local_ref[22][0][2] = -1+2*step;

		local_ref[22][1][0] = -1+2*step;
		local_ref[22][1][1] = -1+1*step;
		local_ref[22][1][2] = -1+2*step;

		local_ref[22][2][0] = -1+2*step;
		local_ref[22][2][1] = -1+2*step;
		local_ref[22][2][2] = -1+2*step;

		local_ref[22][3][0] = -1+1*step;
		local_ref[22][3][1] = -1+2*step;
		local_ref[22][3][2] = -1+2*step;

		local_ref[22][4][0] = -1+1*step;
		local_ref[22][4][1] = -1+1*step;
		local_ref[22][4][2] = -1+3*step;

		local_ref[22][5][0] = -1+2*step;
		local_ref[22][5][1] = -1+1*step;
		local_ref[22][5][2] = -1+3*step;

		local_ref[22][6][0] = -1+2*step;
		local_ref[22][6][1] = -1+2*step;
		local_ref[22][6][2] = -1+3*step;

		local_ref[22][7][0] = -1+1*step;
		local_ref[22][7][1] = -1+2*step;
		local_ref[22][7][2] = -1+3*step;

		//element 23
		local_ref[23][0][0] = -1+2*step;
		local_ref[23][0][1] = -1+1*step;
		local_ref[23][0][2] = -1+2*step;

		local_ref[23][1][0] = -1+3*step;
		local_ref[23][1][1] = -1+1*step;
		local_ref[23][1][2] = -1+2*step;

		local_ref[23][2][0] = -1+3*step;
		local_ref[23][2][1] = -1+2*step;
		local_ref[23][2][2] = -1+2*step;

		local_ref[23][3][0] = -1+2*step;
		local_ref[23][3][1] = -1+2*step;
		local_ref[23][3][2] = -1+2*step;

		local_ref[23][4][0] = -1+2*step;
		local_ref[23][4][1] = -1+1*step;
		local_ref[23][4][2] = -1+3*step;

		local_ref[23][5][0] = -1+3*step;
		local_ref[23][5][1] = -1+1*step;
		local_ref[23][5][2] = -1+3*step;

		local_ref[23][6][0] = -1+3*step;
		local_ref[23][6][1] = -1+2*step;
		local_ref[23][6][2] = -1+3*step;

		local_ref[23][7][0] = -1+2*step;
		local_ref[23][7][1] = -1+2*step;
		local_ref[23][7][2] = -1+3*step;

		//element 24
		local_ref[24][0][0] = -1+0*step;
		local_ref[24][0][1] = -1+2*step;
		local_ref[24][0][2] = -1+2*step;

		local_ref[24][1][0] = -1+1*step;
		local_ref[24][1][1] = -1+2*step;
		local_ref[24][1][2] = -1+2*step;

		local_ref[24][2][0] = -1+1*step;
		local_ref[24][2][1] = -1+3*step;
		local_ref[24][2][2] = -1+2*step;

		local_ref[24][3][0] = -1+0*step;
		local_ref[24][3][1] = -1+3*step;
		local_ref[24][3][2] = -1+2*step;

		local_ref[24][4][0] = -1+0*step;
		local_ref[24][4][1] = -1+2*step;
		local_ref[24][4][2] = -1+3*step;

		local_ref[24][5][0] = -1+1*step;
		local_ref[24][5][1] = -1+2*step;
		local_ref[24][5][2] = -1+3*step;

		local_ref[24][6][0] = -1+1*step;
		local_ref[24][6][1] = -1+3*step;
		local_ref[24][6][2] = -1+3*step;

		local_ref[24][7][0] = -1+0*step;
		local_ref[24][7][1] = -1+3*step;
		local_ref[24][7][2] = -1+3*step;

		//element 25
		local_ref[25][0][0] = -1+1*step;
		local_ref[25][0][1] = -1+2*step;
		local_ref[25][0][2] = -1+2*step;

		local_ref[25][1][0] = -1+2*step;
		local_ref[25][1][1] = -1+2*step;
		local_ref[25][1][2] = -1+2*step;

		local_ref[25][2][0] = -1+2*step;
		local_ref[25][2][1] = -1+3*step;
		local_ref[25][2][2] = -1+2*step;

		local_ref[25][3][0] = -1+1*step;
		local_ref[25][3][1] = -1+3*step;
		local_ref[25][3][2] = -1+2*step;

		local_ref[25][4][0] = -1+1*step;
		local_ref[25][4][1] = -1+2*step;
		local_ref[25][4][2] = -1+3*step;

		local_ref[25][5][0] = -1+2*step;
		local_ref[25][5][1] = -1+2*step;
		local_ref[25][5][2] = -1+3*step;

		local_ref[25][6][0] = -1+2*step;
		local_ref[25][6][1] = -1+3*step;
		local_ref[25][6][2] = -1+3*step;

		local_ref[25][7][0] = -1+1*step;
		local_ref[25][7][1] = -1+3*step;
		local_ref[25][7][2] = -1+3*step;

		//element 26
		local_ref[26][0][0] = -1+2*step;
		local_ref[26][0][1] = -1+2*step;
		local_ref[26][0][2] = -1+2*step;

		local_ref[26][1][0] = -1+3*step;
		local_ref[26][1][1] = -1+2*step;
		local_ref[26][1][2] = -1+2*step;

		local_ref[26][2][0] = -1+3*step;
		local_ref[26][2][1] = -1+3*step;
		local_ref[26][2][2] = -1+2*step;

		local_ref[26][3][0] = -1+2*step;
		local_ref[26][3][1] = -1+3*step;
		local_ref[26][3][2] = -1+2*step;

		local_ref[26][4][0] = -1+2*step;
		local_ref[26][4][1] = -1+2*step;
		local_ref[26][4][2] = -1+3*step;

		local_ref[26][5][0] = -1+3*step;
		local_ref[26][5][1] = -1+2*step;
		local_ref[26][5][2] = -1+3*step;

		local_ref[26][6][0] = -1+3*step;
		local_ref[26][6][1] = -1+3*step;
		local_ref[26][6][2] = -1+3*step;

		local_ref[26][7][0] = -1+2*step;
		local_ref[26][7][1] = -1+3*step;
		local_ref[26][7][2] = -1+3*step;
	}

	int id_node[8];
	std::vector<int> ord;
	for(int ino = 0; ino < 8;ino++){
		id_node[ino] = elem->nodes[ino].id;
		ord.push_back(ino);
	}

	for(int iel = 0; iel < 27; iel++){
		double local_ref_x[8];
		double local_ref_y[8];
		double local_ref_z[8];
		for(int ino = 0; ino < 8; ino++){
			local_ref_x[ino] = local_ref[iel][ino][0];
			local_ref_y[ino] = local_ref[iel][ino][1];
			local_ref_z[ino] = local_ref[iel][ino][2];
		}
		ApplyElement(mesh, coords,  id,  iel,  id_node,  local_ref_x, local_ref_y, local_ref_z,  ord,  hash_nodes);
	}
	sc_array_reset(&toto);
}

void ApplyOctreeTemplate(hexa_tree_t* mesh, std::vector<double>& coords) {

	bool clamped = true;
	sc_hash_array_t*   hash_nodes  = (sc_hash_array_t *)sc_hash_array_new(sizeof(octant_node_t), node_hash_fn , node_equal_fn, &clamped);

	for(int ino = 0; ino < mesh->nodes.elem_count; ino++){
		size_t position;
		octant_node_t *r;
		octant_node_t key;
		octant_node_t* node = (octant_node_t*) sc_array_index (&mesh->nodes, ino);
		key.x = node->x;
		key.y = node->y;
		key.z = node->z;
		r = (octant_node_t*) sc_hash_array_insert_unique(hash_nodes, &key, &position);
		if(r!=NULL){
			r->x = node->x;
			r->y = node->y;
			r->z = node->z;
			r->id = node->id;
		}else{
			printf("Verificar o no numero %d\n",node->id);
			octant_node_t* node_i = (octant_node_t*) sc_array_index (&hash_nodes->a, position);
			printf("Ele foi confundido com o no %d\n", node_i->id);
		}
	}

	assert(hash_nodes->a.elem_count == mesh->nodes.elem_count);

	std::vector<int> elements_ids;
	//for(int ioc = 0; ioc <mesh->oct.elem_count; ioc++){
	//	octree_t * oc = (octree_t*) sc_array_index (&mesh->oct, ioc);
	//	for(int iel = 0; iel<8; iel++) elements_ids.push_back(oc->id[iel]);
	//}
	for(int iel = 0; iel <mesh->elements.elem_count; iel++){
		octant_t * elem = (octant_t*) sc_array_index (&mesh->elements, iel);
		if(elem->pad!=0) elements_ids.push_back(elem->id);
	}

	/////////////////////////////////////////////////////
	//DEBUG
	/////////////////////////////////////////////////////
	//std::vector<int>().swap(elements_ids);
	//elements_ids.push_back(0);
	/////////////////////////////////////////////////////
	//DEBUG
	/////////////////////////////////////////////////////
	for (int iel = 0; iel < elements_ids.size(); ++iel) {

		octant_t *elem = (octant_t*) sc_array_index(&mesh->elements, elements_ids[iel]);
		/////////////////////////////////////////////////////
		//DEBUG
		/////////////////////////////////////////////////////
		//elem->pad = 12;
		//elem->tem = 1;
		/////////////////////////////////////////////////////
		//DEBUG
		/////////////////////////////////////////////////////
		//template 1
		if(elem->tem==1){
			ApplyTemplate1(mesh, coords, elements_ids, iel, hash_nodes);
		}

		//template 2
		else if(elem->tem==2){
			ApplyTemplate2(mesh, coords, elements_ids, iel, hash_nodes);
		}

		//template 3
		else if(elem->tem==3){
			ApplyTemplate3(mesh, coords, elements_ids, iel, hash_nodes);
		}

		//template 4
		else if(elem->tem==4){
			ApplyTemplate4(mesh, coords, elements_ids, iel, hash_nodes);
		}

		//template 5
		else if(elem->tem==5){
			ApplyTemplate5(mesh, coords, elements_ids, iel, hash_nodes);
		}

		//template 6
		else if(elem->tem==6){
			ApplyTemplate6(mesh, coords, elements_ids, iel, hash_nodes);
		}

		//template 7
		else if(elem->tem==7){
			ApplyTemplate7(mesh, coords, elements_ids, iel, hash_nodes);
		}

		//template 8
		else if(elem->tem==8){
			ApplyTemplate8(mesh, coords, elements_ids, iel, hash_nodes);
		}

		//template 9
		else if(elem->tem==9){
			ApplyTemplate9(mesh, coords, elements_ids, iel, hash_nodes);
		}

		//template 10
		else if(elem->tem==10){
			ApplyTemplate10(mesh, coords, elements_ids, iel, hash_nodes);
		}

		//template 11
		else if(elem->tem==11){
			ApplyTemplate11(mesh, coords, elements_ids, iel, hash_nodes);
		}

	}

	//update the vectors
	mesh->local_n_elements = mesh->elements.elem_count;
	mesh->local_n_nodes = mesh->nodes.elem_count;
	MPI_Allreduce(&mesh->local_n_elements, &mesh->total_n_elements, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(&mesh->local_n_nodes, &mesh->total_n_nodes, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);

	if (mesh->mpi_rank == 0 && false) {
		printf("Total number of elements: %lld\n", mesh->local_n_elements);
		printf("Total number of nodes: %lld\n", mesh->local_n_nodes);
	}

	//redefine the connectivity data for write the vtk file
	mesh->part_nodes = NULL;
	mesh->part_nodes = (int*) malloc (mesh->local_n_nodes*sizeof(int));
	for(int i =0; i < mesh->local_n_nodes; i++) mesh->part_nodes[i] = mesh->mpi_rank;


	sc_hash_array_t* shared_nodes    = (sc_hash_array_t *)sc_hash_array_new(sizeof (shared_node_t), node_shared_hash_fn, node_shared_equal_fn, &clamped);
	//insert the shared nodes in the hash_array
	for(int i = 0; i < mesh->nodes.elem_count; i++)
	{
		octant_node_t* node = (octant_node_t*) sc_array_index (&mesh->nodes, i);
		if(node->y == mesh->y_start)
		{
			if(node->x == mesh->x_start)
				hexa_insert_shared_node(shared_nodes,node,mesh->neighbors[0]);
			if(node->x == mesh->x_end)
				hexa_insert_shared_node(shared_nodes,node,mesh->neighbors[2]);

			hexa_insert_shared_node(shared_nodes,node,mesh->neighbors[1]);
			continue;
		}

		if(node->y == mesh->y_end) {
			if(node->x == mesh->x_start)
				hexa_insert_shared_node(shared_nodes,node,mesh->neighbors[6]);
			if(node->x == mesh->x_end)
				hexa_insert_shared_node(shared_nodes,node,mesh->neighbors[8]);

			hexa_insert_shared_node(shared_nodes,node,mesh->neighbors[7]);
			continue;
		}

		if(node->x == mesh->x_start)
			hexa_insert_shared_node(shared_nodes,node,mesh->neighbors[3]);
		if(node->x == mesh->x_end)
			hexa_insert_shared_node(shared_nodes,node,mesh->neighbors[5]);
	}

	//TODO redo this
	/*
	double xs;
	double xe;
	double ys;
	double ye;

	for(int n = 0;n<mesh->nodes.elem_count;n++){
		octant_node_t* node = (octant_node_t*) sc_array_index (&mesh->nodes, n);
		if(node->y == mesh->y_start){ys=coords[3*node->id+1];}
		if(node->y == mesh->y_end)  {ye=coords[3*node->id+1];}
		if(node->x == mesh->x_start)  {xs=coords[3*node->id];}
		if(node->x == mesh->x_end)    {xe=coords[3*node->id];}
	}

	double fy = (ye-ys)/1000;
	double fx = (xe-xs)/1000;

	ys += fy;
	xs += fx;
	ye -= fy;
	xe -= fx;

	for (int iel = 0; iel < mesh->elements.elem_count; ++iel) {
		octant_t *elem = (octant_t*) sc_array_index(&mesh->elements, iel);
		elem->id = iel;
		if(false){
			fprintf(mesh->fdbg,"El:%d\n",elem->id);
			for(int i = 0; i < 8; i++){
				if(ys>=coords[3*elem->nodes[i].id+1]){
					if(xs>=coords[3*elem->nodes[i].id]){
						hexa_insert_shared_node(shared_nodes,&elem->nodes[i],mesh->neighbors[0]);
					}
					if(xe<=coords[3*elem->nodes[i].id]){
						hexa_insert_shared_node(shared_nodes,&elem->nodes[i],mesh->neighbors[2]);
					}
					hexa_insert_shared_node(shared_nodes,&elem->nodes[i],mesh->neighbors[1]);
					//continue;
				}
				if(ye<=coords[3*elem->nodes[i].id+1]){
					if(xs>=coords[3*elem->nodes[i].id]){
						hexa_insert_shared_node(shared_nodes,&elem->nodes[i],mesh->neighbors[6]);
					}
					if(xe<=coords[3*elem->nodes[i].id]){
						hexa_insert_shared_node(shared_nodes,&elem->nodes[i],mesh->neighbors[8]);
					}
					hexa_insert_shared_node(shared_nodes,&elem->nodes[i],mesh->neighbors[7]);
					//continue;
				}

				if(xs>=coords[3*elem->nodes[i].id]){
					hexa_insert_shared_node(shared_nodes,&elem->nodes[i],mesh->neighbors[3]);
				}
				if(xe<=coords[3*elem->nodes[i].id]){
					hexa_insert_shared_node(shared_nodes,&elem->nodes[i],mesh->neighbors[5]);
				}
			}
		}
	}


#ifdef HEXA_DEBUG_
	if(0){
		fprintf(mesh->fdbg, "Shared Nodes: \n");
		fprintf(mesh->fdbg, "Total: %d\n",shared_nodes->a.elem_count);
		for(int i = 0; i < shared_nodes->a.elem_count; ++i)
		{
			shared_node_t* sn = (shared_node_t*) sc_array_index(&shared_nodes->a,i);
			fprintf(mesh->fdbg, "(%d): %d %d %d\n", sn->id, sn->x, sn->y, sn->z);
			fprintf(mesh->fdbg, "     shared with processors: ");
			for(int j = 0; j < sn->listSz; j++){
				fprintf(mesh->fdbg, "%d ", sn->rankList[j]);
			}
			fprintf(mesh->fdbg, "\n");
		}
	}
#endif

	////////////////
	//extract the share nodes from shared_nodes
	//sc_hash_array_truncate(&mesh->shared_nodes);
	sc_hash_array_rip (shared_nodes, &mesh->shared_nodes);
	sc_array_sort(&mesh->shared_nodes,node_comp);

	if(false){
		int local[3];
		size_t              position;

		local[0] = mesh->local_n_nodes    = mesh->nodes.elem_count;
		local[1] = mesh->local_n_elements = mesh->elements.elem_count;

		// node map
		int not_my_nodes    = 0;
		int my_own_nodes    = 0;
		free(mesh->global_id);
		mesh->global_id     = (int64_t*)malloc(sizeof(int64_t)*mesh->local_n_nodes);
		memset(mesh->global_id,-2,mesh->local_n_nodes*sizeof(int64_t));
		sc_hash_array_t* SendTo   = (sc_hash_array_t *) sc_hash_array_new(sizeof(message_t), processors_hash_fn, processors_equal_fn, &clamped);
		sc_hash_array_t* RecvFrom = (sc_hash_array_t *) sc_hash_array_new(sizeof(message_t), processors_hash_fn, processors_equal_fn, &clamped);

		for(int i = 0; i < mesh->shared_nodes.elem_count; ++i)
		{
			shared_node_t* sn = (shared_node_t*) sc_array_index(&mesh->shared_nodes,i);
			for(int j = 0; j < sn->listSz; j++)
			{
				if(sn->rankList[j] < mesh->mpi_rank) {
					message_t* m = (message_t*)sc_hash_array_insert_unique(SendTo,&sn->rankList[j],&position);
					if(m!=NULL)
					{
						m->rank  = sn->rankList[j];
						sc_array_init(&m->idxs, sizeof(uint32_t));
						uint32_t* p = (uint32_t*) sc_array_push(&m->idxs);
 *p = sn->id;
					} else
					{
						message_t* m = (message_t*)sc_array_index(&SendTo->a, position);
						uint32_t* p = (uint32_t*) sc_array_push(&m->idxs);
 *p = sn->id;
					}
					mesh->global_id[sn->id] = -3;
				}
				else if (sn->rankList[j] > mesh->mpi_rank)
				{
					message_t *m = (message_t*)sc_hash_array_insert_unique(RecvFrom,&sn->rankList[j],&position);
					if(m!=NULL)
					{
						m->rank  = sn->rankList[j];
						sc_array_init(&m->idxs, sizeof(uint32_t));
						uint32_t* p = (uint32_t*) sc_array_push(&m->idxs);
 *p = sn->id;
					} else
					{
						message_t* m = (message_t*)sc_array_index(&RecvFrom->a, position);
						uint32_t* p = (uint32_t*) sc_array_push(&m->idxs);
 *p = sn->id;
					}
					mesh->global_id[sn->id] = -1;
				}
			}

			if(mesh->global_id[sn->id] == -1) not_my_nodes++;
			if(mesh->global_id[sn->id] == -3) my_own_nodes++;

		}

		local[0] -= not_my_nodes;

		sc_hash_array_rip(RecvFrom, &mesh->comm_map.RecvFrom );
		sc_hash_array_rip(SendTo  , &mesh->comm_map.SendTo   );

#ifdef HEXA_DEBUG_
		if(0){
			//nodes
			fprintf(mesh->fdbg,"Nodes:\n");
			fprintf(mesh->fdbg, "Recv from %ld processors\n", mesh->comm_map.RecvFrom.elem_count);
			for(int i=0; i < mesh->comm_map.RecvFrom.elem_count; i++)
			{
				message_t* m = (message_t*) sc_array_index(&mesh->comm_map.RecvFrom, i);
				fprintf(mesh->fdbg, "  \n Recv %ld nodes from %d\n", m->idxs.elem_count, m->rank);
				for(int j =0; j < m->idxs.elem_count; j++)
				{
					int32_t *id = (int32_t*) sc_array_index(&m->idxs,j);
					fprintf(mesh->fdbg, "%d ", *id);
					if( (j+1) % 5 == 0 )fprintf(mesh->fdbg, "\n");
				}

			}
			fprintf(mesh->fdbg,"\n");
			fprintf(mesh->fdbg, "Send to %ld processors\n", mesh->comm_map.SendTo.elem_count);
			for(int i=0; i < mesh->comm_map.SendTo.elem_count; i++)
			{
				message_t* m = (message_t*) sc_array_index(&mesh->comm_map.SendTo, i);
				fprintf(mesh->fdbg, "\n Sending %ld nodes from %d\n", m->idxs.elem_count, m->rank);
				for(int j =0; j < m->idxs.elem_count; j++)
				{
					int32_t *id = (int32_t*) sc_array_index(&m->idxs,j);
					fprintf(mesh->fdbg, "%d ", *id);
					if( (j+1) % 5 == 0 )fprintf(mesh->fdbg, "\n");
				}

			}
		}
#endif

		//size for the nodes message
		mesh->comm_map.max_recvbuf_size = 0;
		for(int i = 0; i < mesh->comm_map.RecvFrom.elem_count; i++)
		{
			message_t* m = (message_t*) sc_array_index(&mesh->comm_map.RecvFrom, i);
			mesh->comm_map.max_recvbuf_size +=  m->idxs.elem_count;
		}

		mesh->comm_map.max_sendbuf_size = 0;
		for(int i = 0; i < mesh->comm_map.SendTo.elem_count; i++)
		{
			message_t* m = (message_t*) sc_array_index(&mesh->comm_map.SendTo, i);
			mesh->comm_map.max_sendbuf_size +=  m->idxs.elem_count;
		}

		mesh->comm_map.nrequests = mesh->comm_map.RecvFrom.elem_count +
				mesh->comm_map.SendTo.elem_count;

		int offset = 0;
		MPI_Scan(&local[0], &offset, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

		offset = offset-local[0];
		int count = offset;
		for(int i=0; i < mesh->local_n_nodes; ++i)
			if(mesh->global_id[i] != -1) mesh->global_id[i] = count++;

		int          n_requests;

		MPI_Request *requests;
		MPI_Status  *statuses;
		long long    *recvbuf;
		long long   *sendbuf;

		n_requests = mesh->comm_map.nrequests;
		recvbuf    = (long long*)malloc(mesh->comm_map.max_recvbuf_size*sizeof(long long));
		sendbuf    = (long long*)malloc(mesh->comm_map.max_sendbuf_size*sizeof(long long));

		requests = (MPI_Request*) malloc (n_requests*sizeof(MPI_Request));
		statuses = (MPI_Status*)  malloc (n_requests*sizeof(MPI_Status));
		int c = 0;

		offset = 0;

		// post all non-blocking receives
		for(int i = 0; i < mesh->comm_map.RecvFrom.elem_count; ++i) {
			message_t *m = (message_t*) sc_array_index(&mesh->comm_map.RecvFrom, i);
			MPI_Irecv(&recvbuf[offset], m->idxs.elem_count, MPI_LONG_LONG, m->rank,0,MPI_COMM_WORLD, &requests[c]);
			offset += m->idxs.elem_count;
			c++;
		}

		assert(offset == mesh->comm_map.max_recvbuf_size);

		offset = 0;
		for(int i = 0; i < mesh->comm_map.SendTo.elem_count; ++i) {
			message_t *m = (message_t*) sc_array_index(&mesh->comm_map.SendTo, i);
			for(int j = 0; j < m->idxs.elem_count; ++j)
			{
				int32_t *id = (int32_t*) sc_array_index(&m->idxs,j);
				sendbuf[offset+j] = (long long) mesh->global_id[*id];
			}
			MPI_Isend(&sendbuf[offset], m->idxs.elem_count, MPI_LONG_LONG, m->rank,0,MPI_COMM_WORLD, &requests[c]);
			offset += m->idxs.elem_count;
			c++;
		}
		assert(offset == mesh->comm_map.max_sendbuf_size);

		assert(c == n_requests);

		MPI_Waitall(n_requests,requests,statuses);

		offset = 0;
		for(int i = 0; i < mesh->comm_map.RecvFrom.elem_count; ++i) {
			message_t *m = (message_t*) sc_array_index(&mesh->comm_map.RecvFrom, i);
			for(int j = 0; j < m->idxs.elem_count; ++j)
			{
				int32_t *id = (int32_t*) sc_array_index(&m->idxs,j);
				mesh->global_id[*id] = recvbuf[offset+j];
			}
			offset += m->idxs.elem_count;
		}

		free(recvbuf);
		free(sendbuf);
		free(requests);
		free(statuses);

		mesh->part_nodes = (int*) malloc (mesh->local_n_nodes*sizeof(int));
		for(int i =0; i < mesh->local_n_nodes; i++)
			mesh->part_nodes[i] = mesh->mpi_rank;

		for(int i=0; i < mesh->comm_map.RecvFrom.elem_count; i++)
		{
			message_t* m = (message_t*) sc_array_index(&mesh->comm_map.RecvFrom, i);
			for(int j =0; j < m->idxs.elem_count; j++)
			{
				int32_t *id = (int32_t*) sc_array_index(&m->idxs,j);
				mesh->part_nodes[*id] = m->rank;
			}
		}

	}//do if false
}

 */

/*
void PillowLayer(hexa_tree_t* mesh, std::vector<double>& coords, std::vector<int>& nodes_b_mat){

	bool clamped = true;
	//criando a hash de nos para evitar nos duplicados no pillowing
	sc_hash_array_t*   hash_nodes  = (sc_hash_array_t *)sc_hash_array_new(sizeof(octant_node_t), node_hash_fn , node_equal_fn, &clamped);
	//hash nodes
	for(int ino = 0; ino < mesh->nodes.elem_count; ino++)
	{
		size_t position;
		octant_node_t *r;
		octant_node_t key;
		octant_node_t* node = (octant_node_t*) sc_array_index (&mesh->nodes, ino);
		key.x = node->x;
		key.y = node->y;
		key.z = node->z;
		r = (octant_node_t*) sc_hash_array_insert_unique(hash_nodes, &key, &position);
		if(r!=NULL){
			r->x = node->x;
			r->y = node->y;
			r->z = node->z;
			r->id = node->id;
		}else{
			printf("Verificar o no numero %d\n",node->id);
			octant_node_t* node_i = (octant_node_t*) sc_array_index (&hash_nodes->a, position);
			printf("Ele foi confundido com o no %d\n", node_i->id);
		}
	}

	assert(hash_nodes->a.elem_count == mesh->nodes.elem_count);
	//hash nodes nodes_b_mat
	sc_hash_array_t*   hash_b_mat  = (sc_hash_array_t *)sc_hash_array_new(sizeof(octant_node_t), node_hash_fn, node_equal_fn, &clamped);

	//fazendo hash nos nodes_b_mat
	for(int ino = 0; ino < nodes_b_mat.size(); ino++)
	{
		size_t position;
		octant_node_t *r;
		octant_node_t key;
		octant_node_t* node = (octant_node_t*) sc_array_index (&mesh->nodes, nodes_b_mat[ino]);
		key.x = node->x;
		key.y = node->y;
		key.z = node->z;

		r = (octant_node_t*) sc_hash_array_insert_unique(hash_b_mat, &key, &position);
		if(r!=NULL){
			r->x = key.x;
			r->y = key.y;
			r->z = key.z;
			r->id = node->id;
		}else{
			printf("Verificar o no numero %d\n",node);
		}
	}



	for(int ioc = 0; ioc < mesh->oct.elem_count; ioc++)
	{
		octree_t* oct = (octree_t*) sc_array_index(&mesh->oct,ioc);

		bool surf_oct[8][6];
		int scount[8];
		int ncount[8];
		//set if the octree has a "normal" behaviour or should be
		//treat in a separte case...
		bool normal = true;
		int mat[8];
		for(int iel = 0; iel < 8; iel++)
		{
			octant_t* elem = (octant_t*) sc_array_index(&mesh->elements, oct->id[iel]);
			mat[iel] = elem->n_mat;
			scount[iel] = 0;
			ncount[iel] = 0;
			for(int isurf = 0; isurf < 6; isurf++)
			{
				surf_oct[iel][isurf] = true;
				for(int ino = 0; ino < 4; ino++)
				{
					size_t position;
					octant_node_t key;
					key.x = elem->nodes[FaceNodesMap[isurf][ino]].x;
					key.y = elem->nodes[FaceNodesMap[isurf][ino]].y;
					key.z = elem->nodes[FaceNodesMap[isurf][ino]].z;

					bool lnode = sc_hash_array_lookup(hash_b_mat, &key, &position);
					if(!lnode) surf_oct[iel][isurf] = false;
				}
				if(surf_oct[iel][isurf]) scount[iel]++;
			}
			if(scount[iel] == 0) normal = false;
		}


		if(normal && false)
		{
			for(int iel = 0; iel < 8; iel++)
			{
				sc_array_t toto;
				sc_array_init(&toto, sizeof(octant_t));
				octant_t* elemOrig = (octant_t*) sc_array_index(&mesh->elements, oct->id[iel]);
				octant_t * elem = (octant_t*) sc_array_push(&toto);
				hexa_element_copy(elemOrig,elem);

				bool surf[6];
				int count = scount[iel];
				for(int isurf = 0; isurf < 6; isurf++)
				{
					surf[isurf] = surf_oct[iel][isurf];
				}

				if(count == 0)
				{

				}
				else if(count == 1)
				{
					double step = 1;
					int id = elem->id;

					//reference element z plane cut
					double local_ref[2][8][3];
					if(true)
					{
						//element 0
						local_ref[0][0][0] = -1;
						local_ref[0][0][1] = -1;
						local_ref[0][0][2] = -1;

						local_ref[0][1][0] = 1;
						local_ref[0][1][1] = -1;
						local_ref[0][1][2] = -1;

						local_ref[0][2][0] = 1;
						local_ref[0][2][1] = 1;
						local_ref[0][2][2] = -1;

						local_ref[0][3][0] = -1;
						local_ref[0][3][1] = 1;
						local_ref[0][3][2] = -1;

						local_ref[0][4][0] = -1;
						local_ref[0][4][1] = -1;
						local_ref[0][4][2] = -1+step;

						local_ref[0][5][0] = 1;
						local_ref[0][5][1] = -1;
						local_ref[0][5][2] = -1+step;

						local_ref[0][6][0] = 1;
						local_ref[0][6][1] = 1;
						local_ref[0][6][2] = -1+step;

						local_ref[0][7][0] = -1;
						local_ref[0][7][1] = 1;
						local_ref[0][7][2] = -1+step;

						//element 1
						local_ref[1][0][0] = -1;
						local_ref[1][0][1] = -1;
						local_ref[1][0][2] = -1+step;

						local_ref[1][1][0] = 1;
						local_ref[1][1][1] = -1;
						local_ref[1][1][2] = -1+step;

						local_ref[1][2][0] = 1;
						local_ref[1][2][1] = 1;
						local_ref[1][2][2] = -1+step;

						local_ref[1][3][0] = -1;
						local_ref[1][3][1] = 1;
						local_ref[1][3][2] = -1+step;

						local_ref[1][4][0] = -1;
						local_ref[1][4][1] = -1;
						local_ref[1][4][2] = -1+2*step;

						local_ref[1][5][0] = 1;
						local_ref[1][5][1] = -1;
						local_ref[1][5][2] = -1+2*step;

						local_ref[1][6][0] = 1;
						local_ref[1][6][1] = 1;
						local_ref[1][6][2] = -1+2*step;

						local_ref[1][7][0] = -1;
						local_ref[1][7][1] = 1;
						local_ref[1][7][2] = -1+2*step;
					}

					//define the rotation of the reference element
					int rot[3];
					int sym[3];
					int id_node[8];
					vector<int> ord;

					if(surf[4] || surf[5])
					{
						//edge 4 5 6 7
						rot[0] = 0;
						rot[1] = 0;
						rot[2] = 0;

						sym[0] = 0;
						sym[1] = 0;
						sym[2] = 0;

						ord = RotateHex(rot,sym);
						for(int node_id=0;node_id<8;node_id++) id_node[node_id] = elem->nodes[ord[node_id]].id;
					}
					else if(surf[2] || surf[3])
					{
						//edge 1 3 9 11
						rot[0] = 1;
						rot[1] = 0;
						rot[2] = 0;

						sym[0] = 0;
						sym[1] = 0;
						sym[2] = 0;

						ord = RotateHex(rot,sym);
						for(int node_id=0;node_id<8;node_id++) id_node[node_id] = elem->nodes[ord[node_id]].id;
					}
					else if(surf[0] || surf[1])
					{
						//edge 0 2 8 10
						rot[0] = 0;
						rot[1] = 1;
						rot[2] = 0;

						sym[0] = 0;
						sym[1] = 0;
						sym[2] = 0;

						ord = RotateHex(rot,sym);
						for(int node_id=0;node_id<8;node_id++)id_node[node_id] = elem->nodes[ord[node_id]].id;
					}

					for(int iel = 0; iel < 2; iel++){
						double local_ref_x[8];
						double local_ref_y[8];
						double local_ref_z[8];
						for(int ino = 0; ino < 8; ino++){
							local_ref_x[ino] = local_ref[iel][ino][0];
							local_ref_y[ino] = local_ref[iel][ino][1];
							local_ref_z[ino] = local_ref[iel][ino][2];
						}
						ApplyElement(mesh, coords,  id,  iel,  id_node,  local_ref_x, local_ref_y, local_ref_z,  ord,  hash_nodes);
					}

				}
				else if(count == 2)
				{

					double step = 1;
					int id = elem->id;

					//reference edge 0
					double local_ref[3][8][3];
					if(true)
					{
						//element 0
						local_ref[0][0][0] = -1;
						local_ref[0][0][1] = -1;
						local_ref[0][0][2] = -1;

						local_ref[0][1][0] = 1;
						local_ref[0][1][1] = -1;
						local_ref[0][1][2] = -1;

						local_ref[0][2][0] = 1;
						local_ref[0][2][1] = 1;
						local_ref[0][2][2] = -1;

						local_ref[0][3][0] = -1;
						local_ref[0][3][1] = 1;
						local_ref[0][3][2] = -1;

						local_ref[0][4][0] = -1;
						local_ref[0][4][1] = -1+step;
						local_ref[0][4][2] = -1+step;

						local_ref[0][5][0] = 1;
						local_ref[0][5][1] = -1+step;
						local_ref[0][5][2] = -1+step;

						local_ref[0][6][0] = 1;
						local_ref[0][6][1] = 1;
						local_ref[0][6][2] = -1+step;

						local_ref[0][7][0] = -1;
						local_ref[0][7][1] = 1;
						local_ref[0][7][2] = -1+step;

						//element 1
						local_ref[1][0][0] = -1;
						local_ref[1][0][1] = -1;
						local_ref[1][0][2] = -1;

						local_ref[1][1][0] = 1;
						local_ref[1][1][1] = -1;
						local_ref[1][1][2] = -1;

						local_ref[1][2][0] = 1;
						local_ref[1][2][1] = -1+step;
						local_ref[1][2][2] = -1+step;

						local_ref[1][3][0] = -1;
						local_ref[1][3][1] = -1+step;
						local_ref[1][3][2] = -1+step;

						local_ref[1][4][0] = -1;
						local_ref[1][4][1] = -1;
						local_ref[1][4][2] = -1+2*step;

						local_ref[1][5][0] = 1;
						local_ref[1][5][1] = -1;
						local_ref[1][5][2] = -1+2*step;

						local_ref[1][6][0] = 1;
						local_ref[1][6][1] = -1+step;
						local_ref[1][6][2] = -1+2*step;

						local_ref[1][7][0] = -1;
						local_ref[1][7][1] = -1+step;
						local_ref[1][7][2] = -1+2*step;

						//element 2
						local_ref[2][0][0] = -1;
						local_ref[2][0][1] = -1+1*step;
						local_ref[2][0][2] = -1+1*step;

						local_ref[2][1][0] = 1;
						local_ref[2][1][1] = -1+1*step;
						local_ref[2][1][2] = -1+step;

						local_ref[2][2][0] = 1;
						local_ref[2][2][1] = 1;
						local_ref[2][2][2] = -1+step;

						local_ref[2][3][0] = -1;
						local_ref[2][3][1] = 1;
						local_ref[2][3][2] = -1+step;

						local_ref[2][4][0] = -1;
						local_ref[2][4][1] = -1+1*step;
						local_ref[2][4][2] = -1+2*step;

						local_ref[2][5][0] = 1;
						local_ref[2][5][1] = -1+1*step;
						local_ref[2][5][2] = -1+2*step;

						local_ref[2][6][0] = 1;
						local_ref[2][6][1] = 1;
						local_ref[2][6][2] = -1+2*step;

						local_ref[2][7][0] = -1;
						local_ref[2][7][1] = 1;
						local_ref[2][7][2] = -1+2*step;
					}

					//define the rotation of the reference element
					int rot[3];
					int sym[3];
					int id_node[8];
					vector<int> ord;

					if(surf[5] && surf[2])
					{
						//edge 0
						rot[0] = 0;
						rot[1] = 0;
						rot[2] = 0;

						sym[0] = 0;
						sym[1] = 0;
						sym[2] = 0;

						ord = RotateHex(rot,sym);
						for(int node_id=0;node_id<8;node_id++) id_node[node_id] = elem->nodes[ord[node_id]].id;
					}
					else if(surf[5] && surf[1])
					{
						//edge 1
						rot[0] = 0;
						rot[1] = 0;
						rot[2] = -1;

						sym[0] = 0;
						sym[1] = 0;
						sym[2] = 0;

						ord = RotateHex(rot,sym);
						for(int node_id=0;node_id<8;node_id++)id_node[node_id] = elem->nodes[ord[node_id]].id;
					}
					else if(surf[5] && surf[3])
					{
						//edge 2
						rot[0] = 0;
						rot[1] = 0;
						rot[2] = 0;

						sym[0] = 0;
						sym[1] = 1;
						sym[2] = 0;

						ord = RotateHex(rot,sym);
						for(int node_id=0;node_id<8;node_id++) id_node[node_id] = elem->nodes[ord[node_id]].id;
					}
					else if(surf[5] && surf[0])
					{
						//edge 3
						rot[0] = 0;
						rot[1] = 0;
						rot[2] = 1;

						sym[0] = 0;
						sym[1] = 0;
						sym[2] = 0;

						ord = RotateHex(rot,sym);
						for(int node_id=0;node_id<8;node_id++) id_node[node_id] = elem->nodes[ord[node_id]].id;
					}
					else if(surf[0] && surf[2])
					{
						//edge 4
						rot[0] = 0;
						rot[1] = -1;
						rot[2] = 0;

						sym[0] = 0;
						sym[1] = 0;
						sym[2] = 0;

						ord = RotateHex(rot,sym);
						for(int node_id=0;node_id<8;node_id++) id_node[node_id] = elem->nodes[ord[node_id]].id;
					}
					else if(surf[2] && surf[1])
					{
						//edge 5
						rot[0] = 0;
						rot[1] = 1;
						rot[2] = 0;

						sym[0] = 0;
						sym[1] = 0;
						sym[2] = 0;

						ord = RotateHex(rot,sym);
						for(int node_id=0;node_id<8;node_id++) id_node[node_id] = elem->nodes[ord[node_id]].id;
					}
					else if(surf[1] && surf[3])
					{
						//edge 6
						rot[0] = 0;
						rot[1] = 1;
						rot[2] = 0;

						sym[0] = 0;
						sym[1] = 1;
						sym[2] = 0;

						ord = RotateHex(rot,sym);
						for(int node_id=0;node_id<8;node_id++) id_node[node_id] = elem->nodes[ord[node_id]].id;
					}
					else if(surf[3] && surf[0])
					{
						//edge 7
						rot[0] = 0;
						rot[1] = -1;
						rot[2] = 0;

						sym[0] = 0;
						sym[1] = 1;
						sym[2] = 0;

						ord = RotateHex(rot,sym);
						for(int node_id=0;node_id<8;node_id++)id_node[node_id] = elem->nodes[ord[node_id]].id;
					}
					else if(surf[4] && surf[2])
					{
						//edge 8
						rot[0] = 0;
						rot[1] = 0;
						rot[2] = 0;

						sym[0] = 0;
						sym[1] = 0;
						sym[2] = 1;

						ord = RotateHex(rot,sym);
						for(int node_id=0;node_id<8;node_id++) id_node[node_id] = elem->nodes[ord[node_id]].id;
					}
					else if(surf[4] && surf[1])
					{
						//edge 9
						rot[0] = 0;
						rot[1] = 0;
						rot[2] = -1;

						sym[0] = 0;
						sym[1] = 0;
						sym[2] = 1;

						ord = RotateHex(rot,sym);
						for(int node_id=0;node_id<8;node_id++)id_node[node_id] = elem->nodes[ord[node_id]].id;
					}
					else if(surf[4] && surf[3])
					{
						//edge 10
						rot[0] = -1;
						rot[1] = 0;
						rot[2] = 0;

						sym[0] = 0;
						sym[1] = 1;
						sym[2] = 0;

						ord = RotateHex(rot,sym);
						for(int node_id=0;node_id<8;node_id++) id_node[node_id] = elem->nodes[ord[node_id]].id;
					}
					else if(surf[4] && surf[0])
					{
						//edge 11
						rot[0] = 0;
						rot[1] = 0;
						rot[2] = 1;

						sym[0] = 0;
						sym[1] = 0;
						sym[2] = 1;

						ord = RotateHex(rot,sym);
						for(int node_id=0;node_id<8;node_id++) id_node[node_id] = elem->nodes[ord[node_id]].id;
					}

					for(int iel = 0; iel < 3; iel++){
						double local_ref_x[8];
						double local_ref_y[8];
						double local_ref_z[8];
						for(int ino = 0; ino < 8; ino++){
							local_ref_x[ino] = local_ref[iel][ino][0];
							local_ref_y[ino] = local_ref[iel][ino][1];
							local_ref_z[ino] = local_ref[iel][ino][2];
						}
						ApplyElement(mesh, coords,  id,  iel,  id_node,  local_ref_x, local_ref_y, local_ref_z,  ord,  hash_nodes);
					}
				}
				else if(count == 3)
				{

					double step = 1;
					int id = elem->id;

					//reference node 0
					double local_ref[4][8][3];
					if(true)
					{
						//element 0
						local_ref[0][0][0] = -1+1*step;
						local_ref[0][0][1] = -1+1*step;
						local_ref[0][0][2] = -1;

						local_ref[0][1][0] = -1+2*step;
						local_ref[0][1][1] = -1+1*step;
						local_ref[0][1][2] = -1;

						local_ref[0][2][0] = -1+2*step;
						local_ref[0][2][1] = -1+2*step;
						local_ref[0][2][2] = -1;

						local_ref[0][3][0] = -1+1*step;
						local_ref[0][3][1] = -1+2*step;
						local_ref[0][3][2] = -1;

						local_ref[0][4][0] = -1+1*step;
						local_ref[0][4][1] = -1+1*step;
						local_ref[0][4][2] = -1+1*step;

						local_ref[0][5][0] = -1+2*step;
						local_ref[0][5][1] = -1+1*step;
						local_ref[0][5][2] = -1+1*step;

						local_ref[0][6][0] = -1+2*step;
						local_ref[0][6][1] = -1+2*step;
						local_ref[0][6][2] = -1+1*step;

						local_ref[0][7][0] = -1+1*step;
						local_ref[0][7][1] = -1+2*step;
						local_ref[0][7][2] = -1+1*step;

						//element 1
						local_ref[1][0][0] = -1;
						local_ref[1][0][1] = -1;
						local_ref[1][0][2] = -1;

						local_ref[1][1][0] =  1;
						local_ref[1][1][1] = -1;
						local_ref[1][1][2] = -1;

						local_ref[1][2][0] =  1;
						local_ref[1][2][1] = -1+1*step;
						local_ref[1][2][2] = -1;

						local_ref[1][3][0] = -1+step;
						local_ref[1][3][1] = -1+step;
						local_ref[1][3][2] = -1;

						local_ref[1][4][0] = -1;
						local_ref[1][4][1] = -1;
						local_ref[1][4][2] = 1;

						local_ref[1][5][0] = 1;
						local_ref[1][5][1] = -1;
						local_ref[1][5][2] = 1;

						local_ref[1][6][0] = 1;
						local_ref[1][6][1] = -1+1*step;
						local_ref[1][6][2] = -1+1*step;

						local_ref[1][7][0] = -1+step;
						local_ref[1][7][1] = -1+step;
						local_ref[1][7][2] = -1+step;

						//element 2
						local_ref[2][0][0] = -1;
						local_ref[2][0][1] = -1;
						local_ref[2][0][2] = -1;

						local_ref[2][1][0] = -1+step;
						local_ref[2][1][1] = -1+step;
						local_ref[2][1][2] = -1;

						local_ref[2][2][0] = -1+1*step;
						local_ref[2][2][1] = -1+2*step;
						local_ref[2][2][2] = -1;

						local_ref[2][3][0] = -1;
						local_ref[2][3][1] = -1+2*step;
						local_ref[2][3][2] = -1;

						local_ref[2][4][0] = -1;
						local_ref[2][4][1] = -1;
						local_ref[2][4][2] = -1+2*step;

						local_ref[2][5][0] = -1+step;
						local_ref[2][5][1] = -1+step;
						local_ref[2][5][2] = -1+step;

						local_ref[2][6][0] = -1+1*step;
						local_ref[2][6][1] = -1+2*step;
						local_ref[2][6][2] = -1+1*step;

						local_ref[2][7][0] = -1;
						local_ref[2][7][1] = -1+2*step;
						local_ref[2][7][2] = -1+2*step;

						//element 3
						local_ref[3][0][0] = -1+1*step;
						local_ref[3][0][1] = -1+1*step;
						local_ref[3][0][2] = -1+1*step;

						local_ref[3][1][0] = -1+2*step;
						local_ref[3][1][1] = -1+1*step;
						local_ref[3][1][2] = -1+1*step;

						local_ref[3][2][0] = -1+2*step;
						local_ref[3][2][1] = -1+2*step;
						local_ref[3][2][2] = -1+1*step;

						local_ref[3][3][0] = -1+1*step;
						local_ref[3][3][1] = -1+2*step;
						local_ref[3][3][2] = -1+1*step;

						local_ref[3][4][0] = -1;
						local_ref[3][4][1] = -1;
						local_ref[3][4][2] =  1;

						local_ref[3][5][0] = 1;
						local_ref[3][5][1] = -1;
						local_ref[3][5][2] = 1;

						local_ref[3][6][0] = 1;
						local_ref[3][6][1] = 1;
						local_ref[3][6][2] = 1;

						local_ref[3][7][0] = -1;
						local_ref[3][7][1] = 1;
						local_ref[3][7][2] = 1;
					}
					//define the rotation of the reference element
					int rot[3];
					int sym[3];
					int id_node[8];
					vector<int> ord;

					if(surf[0] && surf[2] && surf[4]){
						//edge 0 4 3
						rot[0] = 0;
						rot[1] = 0;
						rot[2] = 0;

						sym[0] = 0;
						sym[1] = 0;
						sym[2] = 0;

						ord = RotateHex(rot,sym);

						for(int node_id=0;node_id<8;node_id++){
							id_node[node_id] = elem->nodes[ord[node_id]].id;
						}

					}else if(surf[1] && surf[2] && surf[4]){
						//edge 0 1 5
						rot[0] = 0;
						rot[1] = 0;
						rot[2] = -1;

						sym[0] = 0;
						sym[1] = 0;
						sym[2] = 0;

						ord = RotateHex(rot,sym);

						for(int node_id=0;node_id<8;node_id++){
							id_node[node_id] = elem->nodes[ord[node_id]].id;
						}

					}else if(surf[1] && surf[3] && surf[4]){
						//edge 1 2 6
						rot[0] = 0;
						rot[1] = 0;
						rot[2] = -1;

						sym[0] = 1;
						sym[1] = 0;
						sym[2] = 0;

						ord = RotateHex(rot,sym);

						for(int node_id=0;node_id<8;node_id++){
							id_node[node_id] = elem->nodes[ord[node_id]].id;
						}

					}else if(surf[0] && surf[3] && surf[4]){
						//edge 2 3 7
						rot[0] = 0;
						rot[1] = 0;
						rot[2] = 1;

						sym[0] = 0;
						sym[1] = 0;
						sym[2] = 0;

						ord = RotateHex(rot,sym);

						for(int node_id=0;node_id<8;node_id++){
							id_node[node_id] = elem->nodes[ord[node_id]].id;
						}

					}else if(surf[0] && surf[2] && surf[5]){
						//edge 4 8 11
						rot[0] = 0;
						rot[1] = 0;
						rot[2] = 0;

						sym[0] = 0;
						sym[1] = 0;
						sym[2] = 1;

						ord = RotateHex(rot,sym);

						for(int node_id=0;node_id<8;node_id++){
							id_node[node_id] = elem->nodes[ord[node_id]].id;
						}

					}else if(surf[1] && surf[2] && surf[5]){
						//edge 8 9 5
						rot[0] = 0;
						rot[1] = 0;
						rot[2] = -1;

						sym[0] = 0;
						sym[1] = 0;
						sym[2] = 1;

						ord = RotateHex(rot,sym);

						for(int node_id=0;node_id<8;node_id++){
							id_node[node_id] = elem->nodes[ord[node_id]].id;
						}

					}else if(surf[1] && surf[3] && surf[5]){
						//edge 9 10 6
						rot[0] = 0;
						rot[1] = 0;
						rot[2] = -1;

						sym[0] = 1;
						sym[1] = 0;
						sym[2] = 1;

						ord = RotateHex(rot,sym);

						for(int node_id=0;node_id<8;node_id++){
							id_node[node_id] = elem->nodes[ord[node_id]].id;
						}

					}else if(surf[0] && surf[3] && surf[5]){
						//edge 11 10 7
						rot[0] = 0;
						rot[1] = 0;
						rot[2] = 1;

						sym[0] = 0;
						sym[1] = 0;
						sym[2] = 1;

						ord = RotateHex(rot,sym);

						for(int node_id=0;node_id<8;node_id++){
							id_node[node_id] = elem->nodes[ord[node_id]].id;
						}

					}

					for(int iel = 0; iel < 4; iel++){
						double local_ref_x[8];
						double local_ref_y[8];
						double local_ref_z[8];
						for(int ino = 0; ino < 8; ino++){
							local_ref_x[ino] = local_ref[iel][ino][0];
							local_ref_y[ino] = local_ref[iel][ino][1];
							local_ref_z[ino] = local_ref[iel][ino][2];
						}
						ApplyElement(mesh, coords,  id,  iel,  id_node,  local_ref_x, local_ref_y, local_ref_z,  ord,  hash_nodes);
					}

				}
				else
				{
					printf("El:%d nSurf:%d\n",elem->id,count);
					printf("Nao sei o que fazer ainda!");
					for(int isurf = 0; isurf < 6; isurf++) printf("%d ",surf[isurf]);
					printf("\n");
				}
				sc_array_reset(&toto);
			}
		}

		for(int iel = 0; iel < 8; iel++) printf("%d ",mat[iel]);
		printf("\n");
		for(int iel = 0; iel < 8; iel++) printf("%d ",scount[iel]);
		printf("\n");
	}












}

 */

/*
 *
void PillowingInterface(hexa_tree_t* mesh, std::vector<double>& coords, std::vector<int>& nodes_b_mat){
	int elem_old = mesh->elements.elem_count;
	int nodes_old = mesh->nodes.elem_count;
	bool clamped = true;
	fprintf(mesh->profile,"Time inside PillowingInterface\n");

	//redo the mapping in the nodes
	auto start = std::chrono::steady_clock::now( );
	printf("     Redo Node Mapping\n");
	RedoNodeMapping(mesh);
	auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>( std::chrono::steady_clock::now( ) - start );
	fprintf(mesh->profile,"    Redo the mapping in the nodes %lld millisecond(s).\n",elapsed.count());
	//std::cout << "Redo the mapping in the nodes "<< elapsed.count() <<" millisecond(s)."<< std::endl;

	start = std::chrono::steady_clock::now( );
	//criando a hash de nos para evitar nos duplicados no pillowing
	sc_hash_array_t*   hash_nodes  = (sc_hash_array_t *)sc_hash_array_new(sizeof(octant_node_t), node_hash_fn , node_equal_fn, &clamped);
	//vertex hash
	sc_hash_array_t*	  vertex_hash  = (sc_hash_array_t *)sc_hash_array_new(sizeof (octant_vertex_t), vertex_hash_id, vertex_equal_id, &clamped);
	//edge hash
	sc_hash_array_t* hash_edge_ref = sc_hash_array_new(sizeof (octant_edge_t), id_hash, id_equal, &clamped);

	printf("     Building hashs\n");
	BuildHash(mesh, coords, nodes_b_mat, hash_edge_ref);
	elapsed = std::chrono::duration_cast<std::chrono::milliseconds>( std::chrono::steady_clock::now( ) - start );
	fprintf(mesh->profile,"    Time in hash %lld millisecond(s).\n",elapsed.count());
	//std::cout << "Time in hash "<< elapsed.count() <<" millisecond(s)."<< std::endl;

	start = std::chrono::steady_clock::now( );
	printf("     Check octree template\n");
	CheckOctreeTemplate(mesh, hash_edge_ref);
	elapsed = std::chrono::duration_cast<std::chrono::milliseconds>( std::chrono::steady_clock::now( ) - start );
	fprintf(mesh->profile,"    Time in CheckOctreeTemplate %lld millisecond(s).\n",elapsed.count());
	//std::cout << "Time check "<< elapsed.count() <<" millisecond(s)."<< std::endl;

	start = std::chrono::steady_clock::now( );
	printf("     Apply octree template\n");
	ApplyOctreeTemplate(mesh, coords);
	elapsed = std::chrono::duration_cast<std::chrono::milliseconds>( std::chrono::steady_clock::now( ) - start );
	fprintf(mesh->profile,"    Time in ApplyOctreeTemplate %lld millisecond(s).\n",elapsed.count());
	//std::cout << "Time apply "<< elapsed.count() <<" millisecond(s)."<< std::endl;

	//update the vectors
	mesh->local_n_elements = mesh->elements.elem_count;
	mesh->local_n_nodes = mesh->nodes.elem_count;
	MPI_Allreduce(&mesh->local_n_elements, &mesh->total_n_elements, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(&mesh->local_n_nodes, &mesh->total_n_nodes, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
	printf("     %d nodes were created in the pillowing process\n",mesh->nodes.elem_count - nodes_old);
	printf("     %d elements were created in the pillowing process\n",mesh->elements.elem_count - elem_old);

	free(mesh->part_nodes);
	mesh->part_nodes = (int*) malloc (mesh->local_n_nodes*sizeof(int));
	for (int ino = 0; ino < mesh->local_n_nodes; ino++) {
		mesh->part_nodes[ino]=mesh->mpi_rank;
	}

	//for(int ino = 0; ino < nodes_b_mat.size(); ino++){
	//	mesh->part_nodes[nodes_b_mat[ino]] = 1;
	//}

	start = std::chrono::steady_clock::now( );
	printf("     Surface Identification\n");
	SurfaceIdentification(mesh, coords);
	fprintf(mesh->profile,"    Time in SurfaceIdentification %lld millisecond(s).\n",elapsed.count());
	//std::cout << "Time SurfaceIdentification "<< elapsed.count() <<" millisecond(s)."<< std::endl;
}

 */
