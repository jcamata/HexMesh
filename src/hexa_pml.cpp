
#include <gts.h>
#include <gts.h>
#include <glib.h>
#include <cassert>
#include <vector>
#include <iostream>
using namespace std;
#include <sc.h>
#include <sc_io.h>
#include <sc_containers.h>

#include "hexa.h"
#include "pml.h"
#include "hilbert.h"

/*
int8_t SetNodePML(hexa_tree_t* tree, octant_node_t* node) {
	int8_t pml_id = 0;
	if (node->x == 0) pml_id |= PML_X0;
	if (node->x == 3*tree->ncellx) pml_id |= PML_X1;
	if (node->y == 0) pml_id |= PML_Y0;
	if (node->y == 3*tree->ncelly) pml_id |= PML_Y1;
	if (node->z == 0)               pml_id |= PML_Z0;
	if (node->z == 3*tree->max_z) pml_id |= PML_Z1;
	return pml_id;
}

inline int isX0(int8_t pml_id) {
	return ((pml_id & PML_X0) == PML_X0);
}

inline int isX1(int8_t pml_id) {
	return ((pml_id & PML_X1) == PML_X1);
}

inline int isY0(int8_t pml_id) {
	return ((pml_id & PML_Y0) == PML_Y0);
}

inline int isY1(int8_t pml_id) {
	return ((pml_id & PML_Y1) == PML_Y1);
}

inline int isZ0(int8_t pml_id) {
	return ((pml_id & PML_Z0) == PML_Z0);
}

inline int isZ1(int8_t pml_id) {
	return ((pml_id & PML_Z1) == PML_Z1);
}

void SetElemPML(hexa_tree_t* tree, octant_t *elem) {

	int8_t pml_id[8] = {0, 0, 0, 0, 0, 0, 0, 0};
	for (int i = 0; i < 8; ++i) {
		pml_id[i] = SetNodePML(tree, &elem->nodes[i]);
	}

	for (int face = 0; face < 6; ++face) {
		int no1 = face_map[face][0];
		int no2 = face_map[face][1];
		int no3 = face_map[face][2];
		int no4 = face_map[face][3];

		if (isX0(pml_id[no1]) && isX0(pml_id[no2]) && isX0(pml_id[no3]) && isX0(pml_id[no4]))
			elem->pml_id |= PML_X0;

		if (isX1(pml_id[no1]) && isX1(pml_id[no2]) && isX1(pml_id[no3]) && isX1(pml_id[no4]))
			elem->pml_id |= PML_X1;

		if (isY0(pml_id[no1]) && isY0(pml_id[no2]) && isY0(pml_id[no3]) && isY0(pml_id[no4]))
			elem->pml_id |= PML_Y0;

		if (isY1(pml_id[no1]) && isY1(pml_id[no2]) && isY1(pml_id[no3]) && isY1(pml_id[no4]))
			elem->pml_id |= PML_Y1;

		if (isZ0(pml_id[no1]) && isZ0(pml_id[no2]) && isZ0(pml_id[no3]) && isZ0(pml_id[no4]))
			elem->pml_id |= PML_Z0;

		if (isZ1(pml_id[no1]) && isZ1(pml_id[no2]) && isZ1(pml_id[no3]) && isZ1(pml_id[no4]))
			elem->pml_id |= PML_Z1;

	}

}

inline void SetPMLMask(int8_t* mask, int8_t pml_id) {

	mask[PML_CORNER_X0Y0Z0] = (isX0(pml_id) && isY0(pml_id) && isZ0(pml_id));
	mask[PML_CORNER_X1Y0Z0] = (isX1(pml_id) && isY0(pml_id) && isZ0(pml_id));
	mask[PML_CORNER_X0Y1Z0] = (isX0(pml_id) && isY1(pml_id) && isZ0(pml_id));
	mask[PML_CORNER_X1Y1Z0] = (isX1(pml_id) && isY1(pml_id) && isZ0(pml_id));
	mask[PML_CORNER_X0Y0Z1] = (isX0(pml_id) && isY0(pml_id) && isZ1(pml_id));
	mask[PML_CORNER_X1Y0Z1] = (isX1(pml_id) && isY0(pml_id) && isZ1(pml_id));
	mask[PML_CORNER_X0Y1Z1] = (isX0(pml_id) && isY1(pml_id) && isZ1(pml_id));
	mask[PML_CORNER_X1Y1Z1] = (isX1(pml_id) && isY1(pml_id) && isZ1(pml_id));

	mask[PML_EDGE_Z0_X0] = (isX0(pml_id) && isZ0(pml_id));
	mask[PML_EDGE_Z0_X1] = (isX1(pml_id) && isZ0(pml_id));
	mask[PML_EDGE_Z0_Y0] = (isY0(pml_id) && isZ0(pml_id));
	mask[PML_EDGE_Z0_Y1] = (isY1(pml_id) && isZ0(pml_id));

	mask[PML_EDGE_X0_Y0] = (isX0(pml_id) && isY0(pml_id));
	mask[PML_EDGE_X0_Y1] = (isX0(pml_id) && isY1(pml_id));
	mask[PML_EDGE_X1_Y0] = (isX1(pml_id) && isY0(pml_id));
	mask[PML_EDGE_X1_Y1] = (isX1(pml_id) && isY1(pml_id));

	mask[PML_EDGE_Z1_X0] = (isX0(pml_id) && isZ1(pml_id));
	mask[PML_EDGE_Z1_X1] = (isX1(pml_id) && isZ1(pml_id));
	mask[PML_EDGE_Z1_Y0] = (isY0(pml_id) && isZ1(pml_id));
	mask[PML_EDGE_Z1_Y1] = (isY1(pml_id) && isZ1(pml_id));

	mask[PML_FACE_X0] = isX0(pml_id);
	mask[PML_FACE_X1] = isX1(pml_id);
	mask[PML_FACE_Y0] = isY0(pml_id);
	mask[PML_FACE_Y1] = isY1(pml_id);
	mask[PML_FACE_Z0] = isZ0(pml_id);
	mask[PML_FACE_Z1] = isZ1(pml_id);
}

inline void SetPMLMask_corner(int8_t* mask, int8_t pml_id) {

	mask[PML_CORNER_X0Y0Z0] = (isX0(pml_id) && isY0(pml_id) && isZ0(pml_id));
	mask[PML_CORNER_X1Y0Z0] = (isX1(pml_id) && isY0(pml_id) && isZ0(pml_id));
	mask[PML_CORNER_X0Y1Z0] = (isX0(pml_id) && isY1(pml_id) && isZ0(pml_id));
	mask[PML_CORNER_X1Y1Z0] = (isX1(pml_id) && isY1(pml_id) && isZ0(pml_id));
	mask[PML_CORNER_X0Y0Z1] = (isX0(pml_id) && isY0(pml_id) && isZ1(pml_id));
	mask[PML_CORNER_X1Y0Z1] = (isX1(pml_id) && isY0(pml_id) && isZ1(pml_id));
	mask[PML_CORNER_X0Y1Z1] = (isX0(pml_id) && isY1(pml_id) && isZ1(pml_id));
	mask[PML_CORNER_X1Y1Z1] = (isX1(pml_id) && isY1(pml_id) && isZ1(pml_id));

}

inline void SetPMLMask_edge(int8_t* mask, int8_t pml_id) {

	mask[PML_EDGE_Z0_X0] = (isX0(pml_id) && isZ0(pml_id));
	mask[PML_EDGE_Z0_X1] = (isX1(pml_id) && isZ0(pml_id));
	mask[PML_EDGE_Z0_Y0] = (isY0(pml_id) && isZ0(pml_id));
	mask[PML_EDGE_Z0_Y1] = (isY1(pml_id) && isZ0(pml_id));

	mask[PML_EDGE_X0_Y0] = (isX0(pml_id) && isY0(pml_id));
	mask[PML_EDGE_X0_Y1] = (isX0(pml_id) && isY1(pml_id));
	mask[PML_EDGE_X1_Y0] = (isX1(pml_id) && isY0(pml_id));
	mask[PML_EDGE_X1_Y1] = (isX1(pml_id) && isY1(pml_id));

	mask[PML_EDGE_Z1_X0] = (isX0(pml_id) && isZ1(pml_id));
	mask[PML_EDGE_Z1_X1] = (isX1(pml_id) && isZ1(pml_id));
	mask[PML_EDGE_Z1_Y0] = (isY0(pml_id) && isZ1(pml_id));
	mask[PML_EDGE_Z1_Y1] = (isY1(pml_id) && isZ1(pml_id));

}

inline void SetPMLMask_face(int8_t* mask, int8_t pml_id) {

	mask[PML_FACE_X0] = isX0(pml_id);
	mask[PML_FACE_X1] = isX1(pml_id);
	mask[PML_FACE_Y0] = isY0(pml_id);
	mask[PML_FACE_Y1] = isY1(pml_id);
	mask[PML_FACE_Z0] = isZ0(pml_id);
	mask[PML_FACE_Z1] = isZ1(pml_id);
}

 */
unsigned edge_hash_fn(const void *v, const void *u) {
	const node_t *q = (const node_t*) v;
	uint64_t a, b, c;

	a = (double_t) q->coord[0];
	b = (double_t) q->coord[1];
	c = (double_t) q->coord[2];

	sc_hash_mix(a, b, c);
	sc_hash_final(a, b, c);
	return (unsigned) c;
}

int edge_equal_fn(const void *v, const void *u, const void *w) {
	const node_t *e1 = (const node_t*) v;
	const node_t *e2 = (const node_t*) u;

	return (unsigned) ((e1->coord[0] == e2->coord[0]) &&
			(e1->coord[1] == e2->coord[1]) &&
			(e1->coord[2] == e2->coord[2]));
}

void RedoMap(hexa_tree_t* mesh, int layers_x, int layers_y, int layers_z){
	// add the number of PML layers in the mesh
	// it allow us add int points in the mesh
	// keeping a structured mesh
	for(int iel = 0; iel < mesh->elements.elem_count; iel++){
		octant_t * elem = (octant_t*) sc_array_index (&mesh->elements, iel);
		for(int ino = 0; ino < 8; ino++){
			elem->nodes[ino].x = elem->nodes[ino].x + 12*layers_x;
			elem->nodes[ino].y = elem->nodes[ino].y + 12*layers_y;
			elem->nodes[ino].z = elem->nodes[ino].z + 12*layers_z;
		}
		elem->x = 4*elem->x + 4*layers_x;
		elem->y = 4*elem->y + 4*layers_y;
		elem->z = 4*elem->z + 4*layers_z;
	}

	for(int iel = 0; iel < mesh->outsurf.elem_count; iel++){
		octant_t * elem = (octant_t*) sc_array_index (&mesh->outsurf, iel);
		for(int ino = 0; ino < 8; ino++){
			elem->nodes[ino].x = elem->nodes[ino].x + 12*layers_x;
			elem->nodes[ino].y = elem->nodes[ino].y + 12*layers_y;
			elem->nodes[ino].z = elem->nodes[ino].z + 12*layers_z;
		}
		elem->x = 4*elem->x + 4*layers_x;
		elem->y = 4*elem->y + 4*layers_y;
		elem->z = 4*elem->z + 4*layers_z;
	}
	for(int ino = 0; ino < mesh->nodes.elem_count; ino++){
		octant_node_t * node = (octant_node_t*) sc_array_index (&mesh->nodes, ino);
		node->x = node->x + 12*layers_x;
		node->y = node->y + 12*layers_y;;
		node->z = node->z + 12*layers_z;
	}
}

void ExtrudePMLElements(hexa_tree_t* mesh, std::vector<double>& coords) {

	const double X_pml = 8e3;
	const double Y_pml = 8e3;
	const double Z_pml = 8e3;

	const int layers_x = 2;
	const int layers_y = 2;
	const int layers_z = 2;
	int mat_count = 25;
	int n_layers = 2;

	//I should create a toto sc_array
	//it avoid segmentation fault when we perform a
	//push in mesh->elements sc_array due to the
	//realocation of the sc_array
	sc_array_t toto;
	sc_array_init(&toto, sizeof(octant_t));

	RedoMap(mesh,layers_x,layers_y,layers_z);

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
	bool edge, face, point;
	point = true;
	face = true;
	edge = true;

	//mesh->outsurf.elem_count
	for (int i = 0; i < mesh->outsurf.elem_count; ++i) {
		octant_t* elemOrig = (octant_t*) sc_array_index(&mesh->outsurf, i);
		octant_t * elem = (octant_t*) sc_array_push(&toto);

		hexa_element_copy(elemOrig,elem);

		//for(int ino = 0; ino < 8; ino++) printf("%d %d %d %d\n",elem->nodes[ino].id,elem->nodes[ino].x,elem->nodes[ino].y,elem->nodes[ino].z);


		if(face){
			if(elem->surf[0].ext){
				for(int n_l = 0; n_l < layers_x; ++n_l ){

					octant_t* pml_e = (octant_t*) sc_array_push(&mesh->elements);
					pml_e->id = mesh->elements.elem_count+1;

					//nos de referencia
					int aux[4] = {0,3,7,4};
					int node0 = elem->nodes[aux[0]].id;
					int node1 = elem->nodes[aux[1]].id;
					int node2 = elem->nodes[aux[2]].id;
					int node3 = elem->nodes[aux[3]].id;
					double x[8],y[8],z[8];

					x[0] = coords[3*node0+0] - n_l*X_pml/layers_x;
					x[1] = coords[3*node1+0] - (n_l+1)*X_pml/layers_x;
					x[2] = coords[3*node2+0] - (n_l+1)*X_pml/layers_x;
					x[3] = coords[3*node3+0] - n_l*X_pml/layers_x;
					x[4] = coords[3*node0+0] - n_l*X_pml/layers_x;
					x[5] = coords[3*node1+0] - (n_l+1)*X_pml/layers_x;
					x[6] = coords[3*node2+0] - (n_l+1)*X_pml/layers_x;
					x[7] = coords[3*node3+0] - n_l*X_pml/layers_x;

					pml_e->nodes[0].x = elem->nodes[aux[0]].x - 12*(n_l+0);
					pml_e->nodes[1].x = elem->nodes[aux[1]].x - 12*(n_l+1);
					pml_e->nodes[2].x = elem->nodes[aux[2]].x - 12*(n_l+1);
					pml_e->nodes[3].x = elem->nodes[aux[3]].x - 12*(n_l+0);
					pml_e->nodes[4].x = elem->nodes[aux[0]].x - 12*(n_l+0);
					pml_e->nodes[5].x = elem->nodes[aux[1]].x - 12*(n_l+1);
					pml_e->nodes[6].x = elem->nodes[aux[2]].x - 12*(n_l+1);
					pml_e->nodes[7].x = elem->nodes[aux[3]].x - 12*(n_l+0);

					y[0] = coords[3*node0+1];
					y[1] = coords[3*node0+1];
					y[2] = coords[3*node1+1];
					y[3] = coords[3*node1+1];
					y[4] = coords[3*node3+1];
					y[5] = coords[3*node3+1];
					y[6] = coords[3*node2+1];
					y[7] = coords[3*node2+1];

					pml_e->nodes[0].y = elem->nodes[aux[0]].y;
					pml_e->nodes[1].y = elem->nodes[aux[0]].y;
					pml_e->nodes[2].y = elem->nodes[aux[1]].y;
					pml_e->nodes[3].y = elem->nodes[aux[1]].y;
					pml_e->nodes[4].y = elem->nodes[aux[3]].y;
					pml_e->nodes[5].y = elem->nodes[aux[3]].y;
					pml_e->nodes[6].y = elem->nodes[aux[2]].y;
					pml_e->nodes[7].y = elem->nodes[aux[2]].y;

					z[0] = coords[3*node0+2];
					z[1] = coords[3*node0+2];
					z[2] = coords[3*node1+2];
					z[3] = coords[3*node1+2];
					z[4] = coords[3*node3+2];
					z[5] = coords[3*node3+2];
					z[6] = coords[3*node2+2];
					z[7] = coords[3*node2+2];

					pml_e->nodes[0].z = elem->nodes[aux[0]].z;
					pml_e->nodes[1].z = elem->nodes[aux[0]].z;
					pml_e->nodes[2].z = elem->nodes[aux[1]].z;
					pml_e->nodes[3].z = elem->nodes[aux[1]].z;
					pml_e->nodes[4].z = elem->nodes[aux[3]].z;
					pml_e->nodes[5].z = elem->nodes[aux[3]].z;
					pml_e->nodes[6].z = elem->nodes[aux[2]].z;
					pml_e->nodes[7].z = elem->nodes[aux[2]].z;

					for(int ino = 0; ino < 8 ; ino++){
						//definindo ponto p a ser adicionado
						GtsPoint* p = gts_point_new(gts_point_class(),x[ino], y[ino], z[ino]);
						//adicionando ponto p
						int x = pml_e->nodes[ino].x;
						int y = pml_e->nodes[ino].y;
						int z = pml_e->nodes[ino].z;
						pml_e->nodes[ino].id = AddPoint( mesh, hash_nodes, p, coords, x, y, z);
					}
					pml_e->n_mat = elem->n_mat+2;
				}
			}

			if(elem->surf[1].ext){
				for(int n_l = 0; n_l < layers_x; ++n_l ){

					octant_t* pml_e = (octant_t*) sc_array_push(&mesh->elements);
					pml_e->id = mesh->elements.elem_count+1;

					//nos de referencia
					int aux[4] = {1,2,6,5};
					int node0 = elem->nodes[aux[0]].id;
					int node1 = elem->nodes[aux[1]].id;
					int node2 = elem->nodes[aux[2]].id;
					int node3 = elem->nodes[aux[3]].id;
					double x[8],y[8],z[8];

					x[0] = coords[3*node0+0] + n_l*X_pml/layers_x;
					x[1] = coords[3*node1+0] + (n_l+1)*X_pml/layers_x;
					x[2] = coords[3*node2+0] + (n_l+1)*X_pml/layers_x;
					x[3] = coords[3*node3+0] + n_l*X_pml/layers_x;
					x[4] = coords[3*node0+0] + n_l*X_pml/layers_x;
					x[5] = coords[3*node1+0] + (n_l+1)*X_pml/layers_x;
					x[6] = coords[3*node2+0] + (n_l+1)*X_pml/layers_x;
					x[7] = coords[3*node3+0] + n_l*X_pml/layers_x;

					pml_e->nodes[0].x = elem->nodes[aux[0]].x + 12*(n_l+0);
					pml_e->nodes[1].x = elem->nodes[aux[1]].x + 12*(n_l+1);
					pml_e->nodes[2].x = elem->nodes[aux[2]].x + 12*(n_l+1);
					pml_e->nodes[3].x = elem->nodes[aux[3]].x + 12*(n_l+0);
					pml_e->nodes[4].x = elem->nodes[aux[0]].x + 12*(n_l+0);
					pml_e->nodes[5].x = elem->nodes[aux[1]].x + 12*(n_l+1);
					pml_e->nodes[6].x = elem->nodes[aux[2]].x + 12*(n_l+1);
					pml_e->nodes[7].x = elem->nodes[aux[3]].x + 12*(n_l+0);

					y[0] = coords[3*node0+1];
					y[1] = coords[3*node0+1];
					y[2] = coords[3*node1+1];
					y[3] = coords[3*node1+1];
					y[4] = coords[3*node3+1];
					y[5] = coords[3*node3+1];
					y[6] = coords[3*node2+1];
					y[7] = coords[3*node2+1];

					pml_e->nodes[0].y = elem->nodes[aux[0]].y;
					pml_e->nodes[1].y = elem->nodes[aux[0]].y;
					pml_e->nodes[2].y = elem->nodes[aux[1]].y;
					pml_e->nodes[3].y = elem->nodes[aux[1]].y;
					pml_e->nodes[4].y = elem->nodes[aux[3]].y;
					pml_e->nodes[5].y = elem->nodes[aux[3]].y;
					pml_e->nodes[6].y = elem->nodes[aux[2]].y;
					pml_e->nodes[7].y = elem->nodes[aux[2]].y;

					z[0] = coords[3*node0+2];
					z[1] = coords[3*node0+2];
					z[2] = coords[3*node1+2];
					z[3] = coords[3*node1+2];
					z[4] = coords[3*node3+2];
					z[5] = coords[3*node3+2];
					z[6] = coords[3*node2+2];
					z[7] = coords[3*node2+2];

					pml_e->nodes[0].z = elem->nodes[aux[0]].z;
					pml_e->nodes[1].z = elem->nodes[aux[0]].z;
					pml_e->nodes[2].z = elem->nodes[aux[1]].z;
					pml_e->nodes[3].z = elem->nodes[aux[1]].z;
					pml_e->nodes[4].z = elem->nodes[aux[3]].z;
					pml_e->nodes[5].z = elem->nodes[aux[3]].z;
					pml_e->nodes[6].z = elem->nodes[aux[2]].z;
					pml_e->nodes[7].z = elem->nodes[aux[2]].z;

					for(int ino = 0; ino < 8 ; ino++){
						//definindo ponto p a ser adicionado
						GtsPoint* p = gts_point_new(gts_point_class(),x[ino], y[ino], z[ino]);
						//adicionando ponto p
						int x = pml_e->nodes[ino].x;
						int y = pml_e->nodes[ino].y;
						int z = pml_e->nodes[ino].z;
						pml_e->nodes[ino].id = AddPoint( mesh, hash_nodes, p, coords, x, y, z);
					}

					pml_e->n_mat = elem->n_mat+3;

				}
			}

			if(elem->surf[2].ext){
				for(int n_l = 0; n_l < layers_y; ++n_l ){

					octant_t* pml_e = (octant_t*) sc_array_push(&mesh->elements);
					pml_e->id = mesh->elements.elem_count+1;

					//nos de referencia
					int aux[4] = {0,1,4,5};
					int node0 = elem->nodes[aux[0]].id;
					int node1 = elem->nodes[aux[1]].id;
					int node2 = elem->nodes[aux[2]].id;
					int node3 = elem->nodes[aux[3]].id;
					double x[8],y[8],z[8];

					x[0] = coords[3*node0+0];
					x[1] = coords[3*node0+0];
					x[2] = coords[3*node1+0];
					x[3] = coords[3*node1+0];
					x[4] = coords[3*node2+0];
					x[5] = coords[3*node2+0];
					x[6] = coords[3*node3+0];
					x[7] = coords[3*node3+0];

					pml_e->nodes[0].x = elem->nodes[aux[0]].x;
					pml_e->nodes[1].x = elem->nodes[aux[0]].x;
					pml_e->nodes[2].x = elem->nodes[aux[1]].x;
					pml_e->nodes[3].x = elem->nodes[aux[1]].x;
					pml_e->nodes[4].x = elem->nodes[aux[2]].x;
					pml_e->nodes[5].x = elem->nodes[aux[2]].x;
					pml_e->nodes[6].x = elem->nodes[aux[3]].x;
					pml_e->nodes[7].x = elem->nodes[aux[3]].x;

					y[0] = coords[3*node0+1]- n_l*Y_pml/layers_y;
					y[1] = coords[3*node1+1]- (n_l+1)*Y_pml/layers_y;
					y[2] = coords[3*node2+1]- (n_l+1)*Y_pml/layers_y;
					y[3] = coords[3*node3+1]- n_l*Y_pml/layers_y;
					y[4] = coords[3*node0+1] -n_l*Y_pml/layers_y;
					y[5] = coords[3*node1+1]- (n_l+1)*Y_pml/layers_y;
					y[6] = coords[3*node2+1]- (n_l+1)*Y_pml/layers_y;
					y[7] = coords[3*node3+1]- n_l*Y_pml/layers_y;

					pml_e->nodes[0].y = elem->nodes[aux[0]].y - 12*(n_l+0);
					pml_e->nodes[1].y = elem->nodes[aux[1]].y - 12*(n_l+1);
					pml_e->nodes[2].y = elem->nodes[aux[2]].y - 12*(n_l+1);
					pml_e->nodes[3].y = elem->nodes[aux[3]].y - 12*(n_l+0);
					pml_e->nodes[4].y = elem->nodes[aux[0]].y - 12*(n_l+0);
					pml_e->nodes[5].y = elem->nodes[aux[1]].y - 12*(n_l+1);
					pml_e->nodes[6].y = elem->nodes[aux[2]].y - 12*(n_l+1);
					pml_e->nodes[7].y = elem->nodes[aux[3]].y - 12*(n_l+0);

					z[0] = coords[3*node0+2];
					z[1] = coords[3*node0+2];
					z[2] = coords[3*node1+2];
					z[3] = coords[3*node1+2];
					z[4] = coords[3*node2+2];
					z[5] = coords[3*node2+2];
					z[6] = coords[3*node3+2];
					z[7] = coords[3*node3+2];

					pml_e->nodes[0].z = elem->nodes[aux[0]].z;
					pml_e->nodes[1].z = elem->nodes[aux[0]].z;
					pml_e->nodes[2].z = elem->nodes[aux[1]].z;
					pml_e->nodes[3].z = elem->nodes[aux[1]].z;
					pml_e->nodes[4].z = elem->nodes[aux[2]].z;
					pml_e->nodes[5].z = elem->nodes[aux[2]].z;
					pml_e->nodes[6].z = elem->nodes[aux[3]].z;
					pml_e->nodes[7].z = elem->nodes[aux[3]].z;

					for(int ino = 0; ino < 8 ; ino++){
						//definindo ponto p a ser adicionado
						GtsPoint* p = gts_point_new(gts_point_class(),x[ino], y[ino], z[ino]);
						//adicionando ponto p
						int x = pml_e->nodes[ino].x;
						int y = pml_e->nodes[ino].y;
						int z = pml_e->nodes[ino].z;
						pml_e->nodes[ino].id = AddPoint( mesh, hash_nodes, p, coords, x, y, z);
					}
					pml_e->n_mat = elem->n_mat+5;

				}
			}

			if(elem->surf[3].ext){
				for(int n_l = 0; n_l < layers_y; ++n_l ){

					octant_t* pml_e = (octant_t*) sc_array_push(&mesh->elements);
					pml_e->id = mesh->elements.elem_count+1;

					//nos de referencia
					int aux[4] = {3,2,7,6};
					int node0 = elem->nodes[aux[0]].id;
					int node1 = elem->nodes[aux[1]].id;
					int node2 = elem->nodes[aux[2]].id;
					int node3 = elem->nodes[aux[3]].id;
					double x[8],y[8],z[8];

					x[0] = coords[3*node0+0];
					x[1] = coords[3*node0+0];
					x[2] = coords[3*node1+0];
					x[3] = coords[3*node1+0];
					x[4] = coords[3*node2+0];
					x[5] = coords[3*node2+0];
					x[6] = coords[3*node3+0];
					x[7] = coords[3*node3+0];

					pml_e->nodes[0].x = elem->nodes[aux[0]].x;
					pml_e->nodes[1].x = elem->nodes[aux[0]].x;
					pml_e->nodes[2].x = elem->nodes[aux[1]].x;
					pml_e->nodes[3].x = elem->nodes[aux[1]].x;
					pml_e->nodes[4].x = elem->nodes[aux[2]].x;
					pml_e->nodes[5].x = elem->nodes[aux[2]].x;
					pml_e->nodes[6].x = elem->nodes[aux[3]].x;
					pml_e->nodes[7].x = elem->nodes[aux[3]].x;

					y[0] = coords[3*node0+1]+ n_l*Y_pml/layers_y;
					y[1] = coords[3*node1+1]+ (n_l+1)*Y_pml/layers_y;
					y[2] = coords[3*node2+1]+ (n_l+1)*Y_pml/layers_y;
					y[3] = coords[3*node3+1]+ n_l*Y_pml/layers_y;
					y[4] = coords[3*node0+1]+ n_l*Y_pml/layers_y;
					y[5] = coords[3*node1+1]+ (n_l+1)*Y_pml/layers_y;
					y[6] = coords[3*node2+1]+ (n_l+1)*Y_pml/layers_y;
					y[7] = coords[3*node3+1]+ n_l*Y_pml/layers_y;

					pml_e->nodes[0].y = elem->nodes[aux[0]].y + 12*(n_l+0);
					pml_e->nodes[1].y = elem->nodes[aux[1]].y + 12*(n_l+1);
					pml_e->nodes[2].y = elem->nodes[aux[2]].y + 12*(n_l+1);
					pml_e->nodes[3].y = elem->nodes[aux[3]].y + 12*(n_l+0);
					pml_e->nodes[4].y = elem->nodes[aux[0]].y + 12*(n_l+0);
					pml_e->nodes[5].y = elem->nodes[aux[1]].y + 12*(n_l+1);
					pml_e->nodes[6].y = elem->nodes[aux[2]].y + 12*(n_l+1);
					pml_e->nodes[7].y = elem->nodes[aux[3]].y + 12*(n_l+0);

					z[0] = coords[3*node0+2];
					z[1] = coords[3*node0+2];
					z[2] = coords[3*node1+2];
					z[3] = coords[3*node1+2];
					z[4] = coords[3*node2+2];
					z[5] = coords[3*node2+2];
					z[6] = coords[3*node3+2];
					z[7] = coords[3*node3+2];

					pml_e->nodes[0].z = elem->nodes[aux[0]].z;
					pml_e->nodes[1].z = elem->nodes[aux[0]].z;
					pml_e->nodes[2].z = elem->nodes[aux[1]].z;
					pml_e->nodes[3].z = elem->nodes[aux[1]].z;
					pml_e->nodes[4].z = elem->nodes[aux[2]].z;
					pml_e->nodes[5].z = elem->nodes[aux[2]].z;
					pml_e->nodes[6].z = elem->nodes[aux[3]].z;
					pml_e->nodes[7].z = elem->nodes[aux[3]].z;

					for(int ino = 0; ino < 8 ; ino++){
						//definindo ponto p a ser adicionado
						GtsPoint* p = gts_point_new(gts_point_class(),x[ino], y[ino], z[ino]);
						//adicionando ponto p
						int x = pml_e->nodes[ino].x;
						int y = pml_e->nodes[ino].y;
						int z = pml_e->nodes[ino].z;
						pml_e->nodes[ino].id = AddPoint( mesh, hash_nodes, p, coords, x, y, z);
					}

					pml_e->n_mat = elem->n_mat+7;

				}
			}

			if(elem->surf[4].ext){
				for(int n_l = 0; n_l < layers_z; ++n_l ){

					octant_t* pml_e = (octant_t*) sc_array_push(&mesh->elements);
					pml_e->id = mesh->elements.elem_count+1;

					//nos de referencia
					int aux[4] = {4,5,6,7};
					int node0 = elem->nodes[aux[0]].id;
					int node1 = elem->nodes[aux[1]].id;
					int node2 = elem->nodes[aux[2]].id;
					int node3 = elem->nodes[aux[3]].id;
					double x[8],y[8],z[8];

					x[0] = coords[3*node0+0];
					x[1] = coords[3*node1+0];
					x[2] = coords[3*node2+0];
					x[3] = coords[3*node3+0];
					x[4] = coords[3*node0+0];
					x[5] = coords[3*node1+0];
					x[6] = coords[3*node2+0];
					x[7] = coords[3*node3+0];

					pml_e->nodes[0].x = elem->nodes[aux[0]].x;
					pml_e->nodes[1].x = elem->nodes[aux[1]].x;
					pml_e->nodes[2].x = elem->nodes[aux[2]].x;
					pml_e->nodes[3].x = elem->nodes[aux[3]].x;
					pml_e->nodes[4].x = elem->nodes[aux[0]].x;
					pml_e->nodes[5].x = elem->nodes[aux[1]].x;
					pml_e->nodes[6].x = elem->nodes[aux[2]].x;
					pml_e->nodes[7].x = elem->nodes[aux[3]].x;

					y[0] = coords[3*node0+1];
					y[1] = coords[3*node1+1];
					y[2] = coords[3*node2+1];
					y[3] = coords[3*node3+1];
					y[4] = coords[3*node0+1];
					y[5] = coords[3*node1+1];
					y[6] = coords[3*node2+1];
					y[7] = coords[3*node3+1];

					pml_e->nodes[0].y = elem->nodes[aux[0]].y;
					pml_e->nodes[1].y = elem->nodes[aux[1]].y;
					pml_e->nodes[2].y = elem->nodes[aux[2]].y;
					pml_e->nodes[3].y = elem->nodes[aux[3]].y;
					pml_e->nodes[4].y = elem->nodes[aux[0]].y;
					pml_e->nodes[5].y = elem->nodes[aux[1]].y;
					pml_e->nodes[6].y = elem->nodes[aux[2]].y;
					pml_e->nodes[7].y = elem->nodes[aux[3]].y;

					z[0] = coords[3*node0+2] - n_l*Z_pml/layers_z;
					z[1] = coords[3*node1+2] - n_l*Z_pml/layers_z;
					z[2] = coords[3*node2+2] - n_l*Z_pml/layers_z;
					z[3] = coords[3*node3+2] - n_l*Z_pml/layers_z;
					z[4] = coords[3*node0+2] - (n_l+1)*Z_pml/layers_z;
					z[5] = coords[3*node1+2] - (n_l+1)*Z_pml/layers_z;
					z[6] = coords[3*node2+2] - (n_l+1)*Z_pml/layers_z;
					z[7] = coords[3*node3+2] - (n_l+1)*Z_pml/layers_z;

					pml_e->nodes[0].z = elem->nodes[aux[0]].z + 12*(n_l+1);
					pml_e->nodes[1].z = elem->nodes[aux[1]].z + 12*(n_l+1);
					pml_e->nodes[2].z = elem->nodes[aux[2]].z + 12*(n_l+1);
					pml_e->nodes[3].z = elem->nodes[aux[3]].z + 12*(n_l+1);
					pml_e->nodes[4].z = elem->nodes[aux[0]].z + 12*(n_l+0);
					pml_e->nodes[5].z = elem->nodes[aux[1]].z + 12*(n_l+0);
					pml_e->nodes[6].z = elem->nodes[aux[2]].z + 12*(n_l+0);
					pml_e->nodes[7].z = elem->nodes[aux[3]].z + 12*(n_l+0);


					for(int ino = 0; ino < 8 ; ino++){
						//definindo ponto p a ser adicionado
						GtsPoint* p = gts_point_new(gts_point_class(),x[ino], y[ino], z[ino]);
						//adicionando ponto p
						int x = pml_e->nodes[ino].x;
						int y = pml_e->nodes[ino].y;
						int z = pml_e->nodes[ino].z;
						pml_e->nodes[ino].id = AddPoint( mesh, hash_nodes, p, coords, x, y, z);
					}

					pml_e->n_mat = elem->n_mat+11;

				}

			}

			//&& false
			if(elem->surf[5].ext && false){
				for(int n_l = 0; n_l < layers_z; ++n_l ){

					octant_t* pml_e = (octant_t*) sc_array_push(&mesh->elements);
					pml_e->id = mesh->elements.elem_count+1;

					//nos de referencia
					int aux[4] = {0,1,2,3};
					int node0 = elem->nodes[aux[0]].id;
					int node1 = elem->nodes[aux[1]].id;
					int node2 = elem->nodes[aux[2]].id;
					int node3 = elem->nodes[aux[3]].id;
					double x[8],y[8],z[8];


					x[0] = coords[3*node0+0];
					x[1] = coords[3*node1+0];
					x[2] = coords[3*node2+0];
					x[3] = coords[3*node3+0];
					x[4] = coords[3*node0+0];
					x[5] = coords[3*node1+0];
					x[6] = coords[3*node2+0];
					x[7] = coords[3*node3+0];

					pml_e->nodes[0].x = elem->nodes[aux[0]].x;
					pml_e->nodes[1].x = elem->nodes[aux[1]].x;
					pml_e->nodes[2].x = elem->nodes[aux[2]].x;
					pml_e->nodes[3].x = elem->nodes[aux[3]].x;
					pml_e->nodes[4].x = elem->nodes[aux[0]].x;
					pml_e->nodes[5].x = elem->nodes[aux[1]].x;
					pml_e->nodes[6].x = elem->nodes[aux[2]].x;
					pml_e->nodes[7].x = elem->nodes[aux[3]].x;

					y[0] = coords[3*node0+1];
					y[1] = coords[3*node1+1];
					y[2] = coords[3*node2+1];
					y[3] = coords[3*node3+1];
					y[4] = coords[3*node0+1];
					y[5] = coords[3*node1+1];
					y[6] = coords[3*node2+1];
					y[7] = coords[3*node3+1];

					pml_e->nodes[0].y = elem->nodes[aux[0]].y;
					pml_e->nodes[1].y = elem->nodes[aux[1]].y;
					pml_e->nodes[2].y = elem->nodes[aux[2]].y;
					pml_e->nodes[3].y = elem->nodes[aux[3]].y;
					pml_e->nodes[4].y = elem->nodes[aux[0]].y;
					pml_e->nodes[5].y = elem->nodes[aux[1]].y;
					pml_e->nodes[6].y = elem->nodes[aux[2]].y;
					pml_e->nodes[7].y = elem->nodes[aux[3]].y;

					z[0] = coords[3*node0+2] + n_l*Z_pml/layers_z;
					z[1] = coords[3*node1+2] + n_l*Z_pml/layers_z;
					z[2] = coords[3*node2+2] + n_l*Z_pml/layers_z;
					z[3] = coords[3*node3+2] + n_l*Z_pml/layers_z;
					z[4] = coords[3*node0+2] + (n_l+1)*Z_pml/layers_z;
					z[5] = coords[3*node1+2] + (n_l+1)*Z_pml/layers_z;
					z[6] = coords[3*node2+2] + (n_l+1)*Z_pml/layers_z;
					z[7] = coords[3*node3+2] + (n_l+1)*Z_pml/layers_z;

					pml_e->nodes[0].z = elem->nodes[aux[0]].z - 12*(n_l+1);
					pml_e->nodes[1].z = elem->nodes[aux[1]].z - 12*(n_l+1);
					pml_e->nodes[2].z = elem->nodes[aux[2]].z - 12*(n_l+1);
					pml_e->nodes[3].z = elem->nodes[aux[3]].z - 12*(n_l+1);
					pml_e->nodes[4].z = elem->nodes[aux[0]].z - 12*(n_l+0);
					pml_e->nodes[5].z = elem->nodes[aux[1]].z - 12*(n_l+0);
					pml_e->nodes[6].z = elem->nodes[aux[2]].z - 12*(n_l+0);
					pml_e->nodes[7].z = elem->nodes[aux[3]].z - 12*(n_l+0);

					for(int ino = 0; ino < 8 ; ino++){
						//definindo ponto p a ser adicionado
						GtsPoint* p = gts_point_new(gts_point_class(),x[ino], y[ino], z[ino]);
						//adicionando ponto p
						int x = pml_e->nodes[ino].x;
						int y = pml_e->nodes[ino].y;
						int z = pml_e->nodes[ino].z;
						pml_e->nodes[ino].id = AddPoint( mesh, hash_nodes, p, coords, x, y, z);
					}

					pml_e->n_mat = elem->n_mat+9;

				}
			}
		}

		if(edge){
			//&& false
			if(elem->edge[0].ref && false){

				for(int nz = 0; nz < layers_z; nz++){
					for(int ny = 0; ny < layers_y; ny++){

						octant_t* pml_e = (octant_t*) sc_array_push(&mesh->elements);
						pml_e->id = mesh->elements.elem_count+1;

						//nos de referencia
						int aux[2] = {0,1};
						int node0 = elem->nodes[aux[0]].id;
						int node1 = elem->nodes[aux[1]].id;
						double x[8],y[8],z[8];

						x[0] = coords[3*node0+0];
						x[1] = coords[3*node1+0];
						x[2] = coords[3*node1+0];
						x[3] = coords[3*node0+0];
						x[4] = coords[3*node0+0];
						x[5] = coords[3*node1+0];
						x[6] = coords[3*node1+0];
						x[7] = coords[3*node0+0];

						pml_e->nodes[0].x = elem->nodes[aux[0]].x;
						pml_e->nodes[1].x = elem->nodes[aux[1]].x;
						pml_e->nodes[2].x = elem->nodes[aux[1]].x;
						pml_e->nodes[3].x = elem->nodes[aux[0]].x;
						pml_e->nodes[4].x = elem->nodes[aux[0]].x;
						pml_e->nodes[5].x = elem->nodes[aux[1]].x;
						pml_e->nodes[6].x = elem->nodes[aux[1]].x;
						pml_e->nodes[7].x = elem->nodes[aux[0]].x;

						y[0] = coords[3*node0+1] - ny*Y_pml/layers_y;
						y[1] = coords[3*node1+1] - ny*Y_pml/layers_y;
						y[2] = coords[3*node1+1] - (ny+1)*Y_pml/layers_y;
						y[3] = coords[3*node0+1] - (ny+1)*Y_pml/layers_y;
						y[4] = coords[3*node0+1] - ny*Y_pml/layers_y;
						y[5] = coords[3*node1+1] - ny*Y_pml/layers_y;
						y[6] = coords[3*node1+1] - (ny+1)*Y_pml/layers_y;
						y[7] = coords[3*node0+1] - (ny+1)*Y_pml/layers_y;

						pml_e->nodes[0].y = elem->nodes[aux[0]].y - 12*(ny+0);
						pml_e->nodes[1].y = elem->nodes[aux[1]].y - 12*(ny+0);
						pml_e->nodes[2].y = elem->nodes[aux[1]].y - 12*(ny+1);
						pml_e->nodes[3].y = elem->nodes[aux[0]].y - 12*(ny+1);
						pml_e->nodes[4].y = elem->nodes[aux[0]].y - 12*(ny+0);
						pml_e->nodes[5].y = elem->nodes[aux[1]].y - 12*(ny+0);
						pml_e->nodes[6].y = elem->nodes[aux[1]].y - 12*(ny+1);
						pml_e->nodes[7].y = elem->nodes[aux[0]].y - 12*(ny+1);

						z[0] = coords[3*node0+2] + nz*Z_pml/layers_z;
						z[1] = coords[3*node1+2] + nz*Z_pml/layers_z;
						z[2] = coords[3*node1+2] + nz*Z_pml/layers_z;
						z[3] = coords[3*node0+2] + nz*Z_pml/layers_z;
						z[4] = coords[3*node0+2] + (nz+1)*Z_pml/layers_z;
						z[5] = coords[3*node1+2] + (nz+1)*Z_pml/layers_z;
						z[6] = coords[3*node1+2] + (nz+1)*Z_pml/layers_z;
						z[7] = coords[3*node0+2] + (nz+1)*Z_pml/layers_z;

						pml_e->nodes[0].z = elem->nodes[aux[0]].z - 12*(nz+1);
						pml_e->nodes[1].z = elem->nodes[aux[1]].z - 12*(nz+1);
						pml_e->nodes[2].z = elem->nodes[aux[1]].z - 12*(nz+1);
						pml_e->nodes[3].z = elem->nodes[aux[0]].z - 12*(nz+1);
						pml_e->nodes[4].z = elem->nodes[aux[0]].z - 12*(nz+0);
						pml_e->nodes[5].z = elem->nodes[aux[1]].z - 12*(nz+0);
						pml_e->nodes[6].z = elem->nodes[aux[1]].z - 12*(nz+0);
						pml_e->nodes[7].z = elem->nodes[aux[0]].z - 12*(nz+0);

						for(int ino = 0; ino < 8 ; ino++){
							//definindo ponto p a ser adicionado
							GtsPoint* p = gts_point_new(gts_point_class(),x[ino], y[ino], z[ino]);
							//adicionando ponto p
							int x = pml_e->nodes[ino].x;
							int y = pml_e->nodes[ino].y;
							int z = pml_e->nodes[ino].z;
							pml_e->nodes[ino].id = AddPoint( mesh, hash_nodes, p, coords, x, y, z);
						}

						pml_e->n_mat = elem->n_mat+29;
					}
				}

			}

			//&& false
			if(elem->edge[1].ref && false){

				for(int nz = 0; nz < layers_z; nz++){
					for(int nx = 0; nx < layers_x; nx++){

						octant_t* pml_e = (octant_t*) sc_array_push(&mesh->elements);
						pml_e->id = mesh->elements.elem_count+1;

						//nos de referencia
						int aux[2] = {1,2};
						int node0 = elem->nodes[aux[0]].id;
						int node1 = elem->nodes[aux[1]].id;
						double x[8],y[8],z[8];

						x[0] = coords[3*node0+0] + nx*X_pml/layers_x;
						x[1] = coords[3*node0+0] + (nx+1)*X_pml/layers_x;
						x[2] = coords[3*node1+0] + (nx+1)*X_pml/layers_x;
						x[3] = coords[3*node1+0] + nx*X_pml/layers_x;
						x[4] = coords[3*node0+0] + nx*X_pml/layers_x;
						x[5] = coords[3*node0+0] + (nx+1)*X_pml/layers_x;
						x[6] = coords[3*node1+0] + (nx+1)*X_pml/layers_x;
						x[7] = coords[3*node1+0] + nx*X_pml/layers_x;

						pml_e->nodes[0].x = elem->nodes[aux[0]].x + 12*(nx+0);
						pml_e->nodes[1].x = elem->nodes[aux[0]].x + 12*(nx+1);
						pml_e->nodes[2].x = elem->nodes[aux[1]].x + 12*(nx+1);
						pml_e->nodes[3].x = elem->nodes[aux[1]].x + 12*(nx+0);
						pml_e->nodes[4].x = elem->nodes[aux[0]].x + 12*(nx+0);
						pml_e->nodes[5].x = elem->nodes[aux[0]].x + 12*(nx+1);
						pml_e->nodes[6].x = elem->nodes[aux[1]].x + 12*(nx+1);
						pml_e->nodes[7].x = elem->nodes[aux[1]].x + 12*(nx+0);

						y[0] = coords[3*node0+1];
						y[1] = coords[3*node0+1];
						y[2] = coords[3*node1+1];
						y[3] = coords[3*node1+1];
						y[4] = coords[3*node0+1];
						y[5] = coords[3*node0+1];
						y[6] = coords[3*node1+1];
						y[7] = coords[3*node1+1];

						pml_e->nodes[0].y = elem->nodes[aux[0]].y;
						pml_e->nodes[1].y = elem->nodes[aux[0]].y;
						pml_e->nodes[2].y = elem->nodes[aux[1]].y;
						pml_e->nodes[3].y = elem->nodes[aux[1]].y;
						pml_e->nodes[4].y = elem->nodes[aux[0]].y;
						pml_e->nodes[5].y = elem->nodes[aux[0]].y;
						pml_e->nodes[6].y = elem->nodes[aux[1]].y;
						pml_e->nodes[7].y = elem->nodes[aux[1]].y;

						z[0] = coords[3*node0+2] + nz*Z_pml/layers_z;
						z[1] = coords[3*node1+2] + nz*Z_pml/layers_z;
						z[2] = coords[3*node1+2] + nz*Z_pml/layers_z;
						z[3] = coords[3*node0+2] + nz*Z_pml/layers_z;
						z[4] = coords[3*node0+2] + (nz+1)*Z_pml/layers_z;
						z[5] = coords[3*node1+2] + (nz+1)*Z_pml/layers_z;
						z[6] = coords[3*node1+2] + (nz+1)*Z_pml/layers_z;
						z[7] = coords[3*node0+2] + (nz+1)*Z_pml/layers_z;

						pml_e->nodes[0].z = elem->nodes[aux[0]].z - 12*(nz+1);
						pml_e->nodes[1].z = elem->nodes[aux[1]].z - 12*(nz+1);
						pml_e->nodes[2].z = elem->nodes[aux[1]].z - 12*(nz+1);
						pml_e->nodes[3].z = elem->nodes[aux[0]].z - 12*(nz+1);
						pml_e->nodes[4].z = elem->nodes[aux[0]].z - 12*(nz+0);
						pml_e->nodes[5].z = elem->nodes[aux[1]].z - 12*(nz+0);
						pml_e->nodes[6].z = elem->nodes[aux[1]].z - 12*(nz+0);
						pml_e->nodes[7].z = elem->nodes[aux[0]].z - 12*(nz+0);

						for(int ino = 0; ino < 8 ; ino++){
							//definindo ponto p a ser adicionado
							GtsPoint* p = gts_point_new(gts_point_class(),x[ino], y[ino], z[ino]);
							//adicionando ponto p
							int x = pml_e->nodes[ino].x;
							int y = pml_e->nodes[ino].y;
							int z = pml_e->nodes[ino].z;
							pml_e->nodes[ino].id = AddPoint( mesh, hash_nodes, p, coords, x, y, z);
						}

						pml_e->n_mat = elem->n_mat+27;
					}
				}

			}

			//&& false
			if(elem->edge[2].ref && false){

				for(int nz = 0; nz < layers_z; nz++){
					for(int ny = 0; ny < layers_y; ny++){

						octant_t* pml_e = (octant_t*) sc_array_push(&mesh->elements);
						pml_e->id = mesh->elements.elem_count+1;

						//nos de referencia
						int aux[2] = {3,2};
						int node0 = elem->nodes[aux[0]].id;
						int node1 = elem->nodes[aux[1]].id;
						double x[8],y[8],z[8];

						x[0] = coords[3*node0+0];
						x[1] = coords[3*node1+0];
						x[2] = coords[3*node1+0];
						x[3] = coords[3*node0+0];
						x[4] = coords[3*node0+0];
						x[5] = coords[3*node1+0];
						x[6] = coords[3*node1+0];
						x[7] = coords[3*node0+0];

						pml_e->nodes[0].x = elem->nodes[aux[0]].x;
						pml_e->nodes[1].x = elem->nodes[aux[1]].x;
						pml_e->nodes[2].x = elem->nodes[aux[1]].x;
						pml_e->nodes[3].x = elem->nodes[aux[0]].x;
						pml_e->nodes[4].x = elem->nodes[aux[0]].x;
						pml_e->nodes[5].x = elem->nodes[aux[1]].x;
						pml_e->nodes[6].x = elem->nodes[aux[1]].x;
						pml_e->nodes[7].x = elem->nodes[aux[0]].x;

						y[0] = coords[3*node0+1] + ny*Y_pml/layers_y;
						y[1] = coords[3*node1+1] + ny*Y_pml/layers_y;
						y[2] = coords[3*node1+1] + (ny+1)*Y_pml/layers_y;
						y[3] = coords[3*node0+1] + (ny+1)*Y_pml/layers_y;
						y[4] = coords[3*node0+1] + ny*Y_pml/layers_y;
						y[5] = coords[3*node1+1] + ny*Y_pml/layers_y;
						y[6] = coords[3*node1+1] + (ny+1)*Y_pml/layers_y;
						y[7] = coords[3*node0+1] + (ny+1)*Y_pml/layers_y;

						pml_e->nodes[0].y = elem->nodes[aux[0]].y + 12*(ny+0);
						pml_e->nodes[1].y = elem->nodes[aux[1]].y + 12*(ny+0);
						pml_e->nodes[2].y = elem->nodes[aux[1]].y + 12*(ny+1);
						pml_e->nodes[3].y = elem->nodes[aux[0]].y + 12*(ny+1);
						pml_e->nodes[4].y = elem->nodes[aux[0]].y + 12*(ny+0);
						pml_e->nodes[5].y = elem->nodes[aux[1]].y + 12*(ny+0);
						pml_e->nodes[6].y = elem->nodes[aux[1]].y + 12*(ny+1);
						pml_e->nodes[7].y = elem->nodes[aux[0]].y + 12*(ny+1);

						z[0] = coords[3*node0+2] + nz*Z_pml/layers_z;
						z[1] = coords[3*node1+2] + nz*Z_pml/layers_z;
						z[2] = coords[3*node1+2] + nz*Z_pml/layers_z;
						z[3] = coords[3*node0+2] + nz*Z_pml/layers_z;
						z[4] = coords[3*node0+2] + (nz+1)*Z_pml/layers_z;
						z[5] = coords[3*node1+2] + (nz+1)*Z_pml/layers_z;
						z[6] = coords[3*node1+2] + (nz+1)*Z_pml/layers_z;
						z[7] = coords[3*node0+2] + (nz+1)*Z_pml/layers_z;

						pml_e->nodes[0].z = elem->nodes[aux[0]].z - 12*(nz+1);
						pml_e->nodes[1].z = elem->nodes[aux[1]].z - 12*(nz+1);
						pml_e->nodes[2].z = elem->nodes[aux[1]].z - 12*(nz+1);
						pml_e->nodes[3].z = elem->nodes[aux[0]].z - 12*(nz+1);
						pml_e->nodes[4].z = elem->nodes[aux[0]].z - 12*(nz+0);
						pml_e->nodes[5].z = elem->nodes[aux[1]].z - 12*(nz+0);
						pml_e->nodes[6].z = elem->nodes[aux[1]].z - 12*(nz+0);
						pml_e->nodes[7].z = elem->nodes[aux[0]].z - 12*(nz+0);

						for(int ino = 0; ino < 8 ; ino++){
							//definindo ponto p a ser adicionado
							GtsPoint* p = gts_point_new(gts_point_class(),x[ino], y[ino], z[ino]);
							//adicionando ponto p
							int x = pml_e->nodes[ino].x;
							int y = pml_e->nodes[ino].y;
							int z = pml_e->nodes[ino].z;
							pml_e->nodes[ino].id = AddPoint( mesh, hash_nodes, p, coords, x, y, z);
						}

						pml_e->n_mat = elem->n_mat+31;
					}
				}

			}

			//&& false
			if(elem->edge[3].ref && false){

				for(int nz = 0; nz < layers_z; nz++){
					for(int nx = 0; nx < layers_x; nx++){

						octant_t* pml_e = (octant_t*) sc_array_push(&mesh->elements);
						pml_e->id = mesh->elements.elem_count+1;

						//nos de referencia
						int aux[2] = {0,3};
						int node0 = elem->nodes[aux[0]].id;
						int node1 = elem->nodes[aux[1]].id;
						double x[8],y[8],z[8];

						x[0] = coords[3*node0+0] - nx*X_pml/layers_x;
						x[1] = coords[3*node0+0] - (nx+1)*X_pml/layers_x;
						x[2] = coords[3*node1+0] - (nx+1)*X_pml/layers_x;
						x[3] = coords[3*node1+0] - nx*X_pml/layers_x;
						x[4] = coords[3*node0+0] - nx*X_pml/layers_x;
						x[5] = coords[3*node0+0] - (nx+1)*X_pml/layers_x;
						x[6] = coords[3*node1+0] - (nx+1)*X_pml/layers_x;
						x[7] = coords[3*node1+0] - nx*X_pml/layers_x;

						pml_e->nodes[0].x = elem->nodes[aux[0]].x - 12*(nx+0);
						pml_e->nodes[1].x = elem->nodes[aux[0]].x - 12*(nx+1);
						pml_e->nodes[2].x = elem->nodes[aux[1]].x - 12*(nx+1);
						pml_e->nodes[3].x = elem->nodes[aux[1]].x - 12*(nx+0);
						pml_e->nodes[4].x = elem->nodes[aux[0]].x - 12*(nx+0);
						pml_e->nodes[5].x = elem->nodes[aux[0]].x - 12*(nx+1);
						pml_e->nodes[6].x = elem->nodes[aux[1]].x - 12*(nx+1);
						pml_e->nodes[7].x = elem->nodes[aux[1]].x - 12*(nx+0);

						y[0] = coords[3*node0+1];
						y[1] = coords[3*node0+1];
						y[2] = coords[3*node1+1];
						y[3] = coords[3*node1+1];
						y[4] = coords[3*node0+1];
						y[5] = coords[3*node0+1];
						y[6] = coords[3*node1+1];
						y[7] = coords[3*node1+1];

						pml_e->nodes[0].y = elem->nodes[aux[0]].y;
						pml_e->nodes[1].y = elem->nodes[aux[0]].y;
						pml_e->nodes[2].y = elem->nodes[aux[1]].y;
						pml_e->nodes[3].y = elem->nodes[aux[1]].y;
						pml_e->nodes[4].y = elem->nodes[aux[0]].y;
						pml_e->nodes[5].y = elem->nodes[aux[0]].y;
						pml_e->nodes[6].y = elem->nodes[aux[1]].y;
						pml_e->nodes[7].y = elem->nodes[aux[1]].y;

						z[0] = coords[3*node0+2] + nz*Z_pml/layers_z;
						z[1] = coords[3*node1+2] + nz*Z_pml/layers_z;
						z[2] = coords[3*node1+2] + nz*Z_pml/layers_z;
						z[3] = coords[3*node0+2] + nz*Z_pml/layers_z;
						z[4] = coords[3*node0+2] + (nz+1)*Z_pml/layers_z;
						z[5] = coords[3*node1+2] + (nz+1)*Z_pml/layers_z;
						z[6] = coords[3*node1+2] + (nz+1)*Z_pml/layers_z;
						z[7] = coords[3*node0+2] + (nz+1)*Z_pml/layers_z;

						pml_e->nodes[0].z = elem->nodes[aux[0]].z - 12*(nz+1);
						pml_e->nodes[1].z = elem->nodes[aux[1]].z - 12*(nz+1);
						pml_e->nodes[2].z = elem->nodes[aux[1]].z - 12*(nz+1);
						pml_e->nodes[3].z = elem->nodes[aux[0]].z - 12*(nz+1);
						pml_e->nodes[4].z = elem->nodes[aux[0]].z - 12*(nz+0);
						pml_e->nodes[5].z = elem->nodes[aux[1]].z - 12*(nz+0);
						pml_e->nodes[6].z = elem->nodes[aux[1]].z - 12*(nz+0);
						pml_e->nodes[7].z = elem->nodes[aux[0]].z - 12*(nz+0);

						for(int ino = 0; ino < 8 ; ino++){
							//definindo ponto p a ser adicionado
							GtsPoint* p = gts_point_new(gts_point_class(),x[ino], y[ino], z[ino]);
							//adicionando ponto p
							int x = pml_e->nodes[ino].x;
							int y = pml_e->nodes[ino].y;
							int z = pml_e->nodes[ino].z;
							pml_e->nodes[ino].id = AddPoint( mesh, hash_nodes, p, coords, x, y, z);
						}

						pml_e->n_mat = elem->n_mat+25;
					}
				}


			}


			if(elem->edge[4].ref){
				for(int nx = 0; nx < layers_x; nx++){
					for(int ny = 0; ny < layers_y; ny++){

						octant_t* pml_e = (octant_t*) sc_array_push(&mesh->elements);
						pml_e->id = mesh->elements.elem_count+1;

						//nos de referencia
						int aux[2] = {0,4};
						int node0 = elem->nodes[aux[0]].id;
						int node1 = elem->nodes[aux[1]].id;
						double x[8],y[8],z[8];

						x[0] = coords[3*node0+0] - nx*X_pml/layers_x;
						x[1] = coords[3*node0+0] - (nx+1)*X_pml/layers_x;
						x[2] = coords[3*node0+0] - (nx+1)*X_pml/layers_x;
						x[3] = coords[3*node0+0] - nx*X_pml/layers_x;
						x[4] = coords[3*node1+0] - nx*X_pml/layers_x;
						x[5] = coords[3*node1+0] - (nx+1)*X_pml/layers_x;
						x[6] = coords[3*node1+0] - (nx+1)*X_pml/layers_x;
						x[7] = coords[3*node1+0] - nx*X_pml/layers_x;

						pml_e->nodes[0].x = elem->nodes[aux[0]].x - 12*(nx+0);
						pml_e->nodes[1].x = elem->nodes[aux[0]].x - 12*(nx+1);
						pml_e->nodes[2].x = elem->nodes[aux[0]].x - 12*(nx+1);
						pml_e->nodes[3].x = elem->nodes[aux[0]].x - 12*(nx+0);
						pml_e->nodes[4].x = elem->nodes[aux[1]].x - 12*(nx+0);
						pml_e->nodes[5].x = elem->nodes[aux[1]].x - 12*(nx+1);
						pml_e->nodes[6].x = elem->nodes[aux[1]].x - 12*(nx+1);
						pml_e->nodes[7].x = elem->nodes[aux[1]].x - 12*(nx+0);

						y[0] = coords[3*node0+1] - ny*Y_pml/layers_y;
						y[1] = coords[3*node0+1] - ny*Y_pml/layers_y;
						y[2] = coords[3*node0+1] - (ny+1)*Y_pml/layers_y;
						y[3] = coords[3*node0+1] - (ny+1)*Y_pml/layers_y;
						y[4] = coords[3*node1+1] - ny*Y_pml/layers_y;
						y[5] = coords[3*node1+1] - ny*Y_pml/layers_y;
						y[6] = coords[3*node1+1] - (ny+1)*Y_pml/layers_y;
						y[7] = coords[3*node1+1] - (ny+1)*Y_pml/layers_y;

						pml_e->nodes[0].y = elem->nodes[aux[0]].y - 12*(ny+0);
						pml_e->nodes[1].y = elem->nodes[aux[0]].y - 12*(ny+0);
						pml_e->nodes[2].y = elem->nodes[aux[0]].y - 12*(ny+1);
						pml_e->nodes[3].y = elem->nodes[aux[0]].y - 12*(ny+1);
						pml_e->nodes[4].y = elem->nodes[aux[1]].y - 12*(ny+0);
						pml_e->nodes[5].y = elem->nodes[aux[1]].y - 12*(ny+0);
						pml_e->nodes[6].y = elem->nodes[aux[1]].y - 12*(ny+1);
						pml_e->nodes[7].y = elem->nodes[aux[1]].y - 12*(ny+1);

						z[0] = coords[3*node0+2];
						z[1] = coords[3*node0+2];
						z[2] = coords[3*node0+2];
						z[3] = coords[3*node0+2];
						z[4] = coords[3*node1+2];
						z[5] = coords[3*node1+2];
						z[6] = coords[3*node1+2];
						z[7] = coords[3*node1+2];

						pml_e->nodes[0].z = elem->nodes[aux[0]].z;
						pml_e->nodes[1].z = elem->nodes[aux[0]].z;
						pml_e->nodes[2].z = elem->nodes[aux[0]].z;
						pml_e->nodes[3].z = elem->nodes[aux[0]].z;
						pml_e->nodes[4].z = elem->nodes[aux[1]].z;
						pml_e->nodes[5].z = elem->nodes[aux[1]].z;
						pml_e->nodes[6].z = elem->nodes[aux[1]].z;
						pml_e->nodes[7].z = elem->nodes[aux[1]].z;

						for(int ino = 0; ino < 8 ; ino++){
							//definindo ponto p a ser adicionado
							GtsPoint* p = gts_point_new(gts_point_class(),x[ino], y[ino], z[ino]);
							//adicionando ponto p
							int x = pml_e->nodes[ino].x;
							int y = pml_e->nodes[ino].y;
							int z = pml_e->nodes[ino].z;
							pml_e->nodes[ino].id = AddPoint( mesh, hash_nodes, p, coords, x, y, z);
						}

						pml_e->n_mat = elem->n_mat+17;
					}
				}


			}

			if(elem->edge[5].ref){
				for(int nx = 0; nx < layers_x; nx++){
					for(int ny = 0; ny < layers_y; ny++){

						octant_t* pml_e = (octant_t*) sc_array_push(&mesh->elements);
						pml_e->id = mesh->elements.elem_count+1;

						//nos de referencia
						int aux[2] = {1,5};
						int node0 = elem->nodes[aux[0]].id;
						int node1 = elem->nodes[aux[1]].id;
						double x[8],y[8],z[8];

						x[0] = coords[3*node0+0] + nx*X_pml/layers_x;
						x[1] = coords[3*node0+0] + (nx+1)*X_pml/layers_x;
						x[2] = coords[3*node0+0] + (nx+1)*X_pml/layers_x;
						x[3] = coords[3*node0+0] + nx*X_pml/layers_x;
						x[4] = coords[3*node1+0] + nx*X_pml/layers_x;
						x[5] = coords[3*node1+0] + (nx+1)*X_pml/layers_x;
						x[6] = coords[3*node1+0] + (nx+1)*X_pml/layers_x;
						x[7] = coords[3*node1+0] + nx*X_pml/layers_x;

						pml_e->nodes[0].x = elem->nodes[aux[0]].x + 12*(nx+0);
						pml_e->nodes[1].x = elem->nodes[aux[0]].x + 12*(nx+1);
						pml_e->nodes[2].x = elem->nodes[aux[0]].x + 12*(nx+1);
						pml_e->nodes[3].x = elem->nodes[aux[0]].x + 12*(nx+0);
						pml_e->nodes[4].x = elem->nodes[aux[1]].x + 12*(nx+0);
						pml_e->nodes[5].x = elem->nodes[aux[1]].x + 12*(nx+1);
						pml_e->nodes[6].x = elem->nodes[aux[1]].x + 12*(nx+1);
						pml_e->nodes[7].x = elem->nodes[aux[1]].x + 12*(nx+0);

						y[0] = coords[3*node0+1] - ny*Y_pml/layers_y;
						y[1] = coords[3*node0+1] - ny*Y_pml/layers_y;
						y[2] = coords[3*node0+1] - (ny+1)*Y_pml/layers_y;
						y[3] = coords[3*node0+1] - (ny+1)*Y_pml/layers_y;
						y[4] = coords[3*node1+1] - ny*Y_pml/layers_y;
						y[5] = coords[3*node1+1] - ny*Y_pml/layers_y;
						y[6] = coords[3*node1+1] - (ny+1)*Y_pml/layers_y;
						y[7] = coords[3*node1+1] - (ny+1)*Y_pml/layers_y;

						pml_e->nodes[0].y = elem->nodes[aux[0]].y - 12*(ny+0);
						pml_e->nodes[1].y = elem->nodes[aux[0]].y - 12*(ny+0);
						pml_e->nodes[2].y = elem->nodes[aux[0]].y - 12*(ny+1);
						pml_e->nodes[3].y = elem->nodes[aux[0]].y - 12*(ny+1);
						pml_e->nodes[4].y = elem->nodes[aux[1]].y - 12*(ny+0);
						pml_e->nodes[5].y = elem->nodes[aux[1]].y - 12*(ny+0);
						pml_e->nodes[6].y = elem->nodes[aux[1]].y - 12*(ny+1);
						pml_e->nodes[7].y = elem->nodes[aux[1]].y - 12*(ny+1);

						z[0] = coords[3*node0+2];
						z[1] = coords[3*node0+2];
						z[2] = coords[3*node0+2];
						z[3] = coords[3*node0+2];
						z[4] = coords[3*node1+2];
						z[5] = coords[3*node1+2];
						z[6] = coords[3*node1+2];
						z[7] = coords[3*node1+2];

						pml_e->nodes[0].z = elem->nodes[aux[0]].z;
						pml_e->nodes[1].z = elem->nodes[aux[0]].z;
						pml_e->nodes[2].z = elem->nodes[aux[0]].z;
						pml_e->nodes[3].z = elem->nodes[aux[0]].z;
						pml_e->nodes[4].z = elem->nodes[aux[1]].z;
						pml_e->nodes[5].z = elem->nodes[aux[1]].z;
						pml_e->nodes[6].z = elem->nodes[aux[1]].z;
						pml_e->nodes[7].z = elem->nodes[aux[1]].z;

						for(int ino = 0; ino < 8 ; ino++){
							//definindo ponto p a ser adicionado
							GtsPoint* p = gts_point_new(gts_point_class(),x[ino], y[ino], z[ino]);
							//adicionando ponto p
							int x = pml_e->nodes[ino].x;
							int y = pml_e->nodes[ino].y;
							int z = pml_e->nodes[ino].z;
							pml_e->nodes[ino].id = AddPoint( mesh, hash_nodes, p, coords, x, y, z);
						}

						pml_e->n_mat = elem->n_mat+21;
					}
				}
			}

			if(elem->edge[6].ref){
				for(int nx = 0; nx < layers_x; nx++){
					for(int ny = 0; ny < layers_y; ny++){

						octant_t* pml_e = (octant_t*) sc_array_push(&mesh->elements);
						pml_e->id = mesh->elements.elem_count+1;

						//nos de referencia
						int aux[2] = {2,6};
						int node0 = elem->nodes[aux[0]].id;
						int node1 = elem->nodes[aux[1]].id;
						double x[8],y[8],z[8];

						x[0] = coords[3*node0+0] + nx*X_pml/layers_x;
						x[1] = coords[3*node0+0] + (nx+1)*X_pml/layers_x;
						x[2] = coords[3*node0+0] + (nx+1)*X_pml/layers_x;
						x[3] = coords[3*node0+0] + nx*X_pml/layers_x;
						x[4] = coords[3*node1+0] + nx*X_pml/layers_x;
						x[5] = coords[3*node1+0] + (nx+1)*X_pml/layers_x;
						x[6] = coords[3*node1+0] + (nx+1)*X_pml/layers_x;
						x[7] = coords[3*node1+0] + nx*X_pml/layers_x;

						pml_e->nodes[0].x = elem->nodes[aux[0]].x + 12*(nx+0);
						pml_e->nodes[1].x = elem->nodes[aux[0]].x + 12*(nx+1);
						pml_e->nodes[2].x = elem->nodes[aux[0]].x + 12*(nx+1);
						pml_e->nodes[3].x = elem->nodes[aux[0]].x + 12*(nx+0);
						pml_e->nodes[4].x = elem->nodes[aux[1]].x + 12*(nx+0);
						pml_e->nodes[5].x = elem->nodes[aux[1]].x + 12*(nx+1);
						pml_e->nodes[6].x = elem->nodes[aux[1]].x + 12*(nx+1);
						pml_e->nodes[7].x = elem->nodes[aux[1]].x + 12*(nx+0);

						y[0] = coords[3*node0+1] + ny*Y_pml/layers_y;
						y[1] = coords[3*node0+1] + ny*Y_pml/layers_y;
						y[2] = coords[3*node0+1] + (ny+1)*Y_pml/layers_y;
						y[3] = coords[3*node0+1] + (ny+1)*Y_pml/layers_y;
						y[4] = coords[3*node1+1] + ny*Y_pml/layers_y;
						y[5] = coords[3*node1+1] + ny*Y_pml/layers_y;
						y[6] = coords[3*node1+1] + (ny+1)*Y_pml/layers_y;
						y[7] = coords[3*node1+1] + (ny+1)*Y_pml/layers_y;

						pml_e->nodes[0].y = elem->nodes[aux[0]].y + 12*(ny+0);
						pml_e->nodes[1].y = elem->nodes[aux[0]].y + 12*(ny+0);
						pml_e->nodes[2].y = elem->nodes[aux[0]].y + 12*(ny+1);
						pml_e->nodes[3].y = elem->nodes[aux[0]].y + 12*(ny+1);
						pml_e->nodes[4].y = elem->nodes[aux[1]].y + 12*(ny+0);
						pml_e->nodes[5].y = elem->nodes[aux[1]].y + 12*(ny+0);
						pml_e->nodes[6].y = elem->nodes[aux[1]].y + 12*(ny+1);
						pml_e->nodes[7].y = elem->nodes[aux[1]].y + 12*(ny+1);

						z[0] = coords[3*node0+2];
						z[1] = coords[3*node0+2];
						z[2] = coords[3*node0+2];
						z[3] = coords[3*node0+2];
						z[4] = coords[3*node1+2];
						z[5] = coords[3*node1+2];
						z[6] = coords[3*node1+2];
						z[7] = coords[3*node1+2];

						pml_e->nodes[0].z = elem->nodes[aux[0]].z;
						pml_e->nodes[1].z = elem->nodes[aux[0]].z;
						pml_e->nodes[2].z = elem->nodes[aux[0]].z;
						pml_e->nodes[3].z = elem->nodes[aux[0]].z;
						pml_e->nodes[4].z = elem->nodes[aux[1]].z;
						pml_e->nodes[5].z = elem->nodes[aux[1]].z;
						pml_e->nodes[6].z = elem->nodes[aux[1]].z;
						pml_e->nodes[7].z = elem->nodes[aux[1]].z;

						for(int ino = 0; ino < 8 ; ino++){
							//definindo ponto p a ser adicionado
							GtsPoint* p = gts_point_new(gts_point_class(),x[ino], y[ino], z[ino]);
							//adicionando ponto p
							int x = pml_e->nodes[ino].x;
							int y = pml_e->nodes[ino].y;
							int z = pml_e->nodes[ino].z;
							pml_e->nodes[ino].id = AddPoint( mesh, hash_nodes, p, coords, x, y, z);
						}

						pml_e->n_mat = elem->n_mat+23;
					}
				}
			}

			if(elem->edge[7].ref){
				for(int nx = 0; nx < layers_x; nx++){
					for(int ny = 0; ny < layers_y; ny++){

						octant_t* pml_e = (octant_t*) sc_array_push(&mesh->elements);
						pml_e->id = mesh->elements.elem_count+1;

						//nos de referencia
						int aux[2] = {3,7};
						int node0 = elem->nodes[aux[0]].id;
						int node1 = elem->nodes[aux[1]].id;
						double x[8],y[8],z[8];

						x[0] = coords[3*node0+0] - nx*X_pml/layers_x;
						x[1] = coords[3*node0+0] - (nx+1)*X_pml/layers_x;
						x[2] = coords[3*node0+0] - (nx+1)*X_pml/layers_x;
						x[3] = coords[3*node0+0] - nx*X_pml/layers_x;
						x[4] = coords[3*node1+0] - nx*X_pml/layers_x;
						x[5] = coords[3*node1+0] - (nx+1)*X_pml/layers_x;
						x[6] = coords[3*node1+0] - (nx+1)*X_pml/layers_x;
						x[7] = coords[3*node1+0] - nx*X_pml/layers_x;

						pml_e->nodes[0].x = elem->nodes[aux[0]].x - 12*(nx+0);
						pml_e->nodes[1].x = elem->nodes[aux[0]].x - 12*(nx+1);
						pml_e->nodes[2].x = elem->nodes[aux[0]].x - 12*(nx+1);
						pml_e->nodes[3].x = elem->nodes[aux[0]].x - 12*(nx+0);
						pml_e->nodes[4].x = elem->nodes[aux[1]].x - 12*(nx+0);
						pml_e->nodes[5].x = elem->nodes[aux[1]].x - 12*(nx+1);
						pml_e->nodes[6].x = elem->nodes[aux[1]].x - 12*(nx+1);
						pml_e->nodes[7].x = elem->nodes[aux[1]].x - 12*(nx+0);

						y[0] = coords[3*node0+1] + ny*Y_pml/layers_y;
						y[1] = coords[3*node0+1] + ny*Y_pml/layers_y;
						y[2] = coords[3*node0+1] + (ny+1)*Y_pml/layers_y;
						y[3] = coords[3*node0+1] + (ny+1)*Y_pml/layers_y;
						y[4] = coords[3*node1+1] + ny*Y_pml/layers_y;
						y[5] = coords[3*node1+1] + ny*Y_pml/layers_y;
						y[6] = coords[3*node1+1] + (ny+1)*Y_pml/layers_y;
						y[7] = coords[3*node1+1] + (ny+1)*Y_pml/layers_y;

						pml_e->nodes[0].y = elem->nodes[aux[0]].y + 12*(ny+0);
						pml_e->nodes[1].y = elem->nodes[aux[0]].y + 12*(ny+0);
						pml_e->nodes[2].y = elem->nodes[aux[0]].y + 12*(ny+1);
						pml_e->nodes[3].y = elem->nodes[aux[0]].y + 12*(ny+1);
						pml_e->nodes[4].y = elem->nodes[aux[1]].y + 12*(ny+0);
						pml_e->nodes[5].y = elem->nodes[aux[1]].y + 12*(ny+0);
						pml_e->nodes[6].y = elem->nodes[aux[1]].y + 12*(ny+1);
						pml_e->nodes[7].y = elem->nodes[aux[1]].y + 12*(ny+1);

						z[0] = coords[3*node0+2];
						z[1] = coords[3*node0+2];
						z[2] = coords[3*node0+2];
						z[3] = coords[3*node0+2];
						z[4] = coords[3*node1+2];
						z[5] = coords[3*node1+2];
						z[6] = coords[3*node1+2];
						z[7] = coords[3*node1+2];

						pml_e->nodes[0].z = elem->nodes[aux[0]].z;
						pml_e->nodes[1].z = elem->nodes[aux[0]].z;
						pml_e->nodes[2].z = elem->nodes[aux[0]].z;
						pml_e->nodes[3].z = elem->nodes[aux[0]].z;
						pml_e->nodes[4].z = elem->nodes[aux[1]].z;
						pml_e->nodes[5].z = elem->nodes[aux[1]].z;
						pml_e->nodes[6].z = elem->nodes[aux[1]].z;
						pml_e->nodes[7].z = elem->nodes[aux[1]].z;

						for(int ino = 0; ino < 8 ; ino++){
							//definindo ponto p a ser adicionado
							GtsPoint* p = gts_point_new(gts_point_class(),x[ino], y[ino], z[ino]);
							//adicionando ponto p
							int x = pml_e->nodes[ino].x;
							int y = pml_e->nodes[ino].y;
							int z = pml_e->nodes[ino].z;
							pml_e->nodes[ino].id = AddPoint( mesh, hash_nodes, p, coords, x, y, z);
						}

						pml_e->n_mat = elem->n_mat+19;
					}
				}
			}


			if(elem->edge[8].ref){
				for(int nz = 0; nz < layers_z; nz++){
					for(int ny = 0; ny < layers_y; ny++){

						octant_t* pml_e = (octant_t*) sc_array_push(&mesh->elements);
						pml_e->id = mesh->elements.elem_count+1;

						//nos de referencia
						int aux[2] = {4,5};
						int node0 = elem->nodes[aux[0]].id;
						int node1 = elem->nodes[aux[1]].id;
						double x[8],y[8],z[8];

						x[0] = coords[3*node0+0];
						x[1] = coords[3*node1+0];
						x[2] = coords[3*node1+0];
						x[3] = coords[3*node0+0];
						x[4] = coords[3*node0+0];
						x[5] = coords[3*node1+0];
						x[6] = coords[3*node1+0];
						x[7] = coords[3*node0+0];

						pml_e->nodes[0].x = elem->nodes[aux[0]].x;
						pml_e->nodes[1].x = elem->nodes[aux[1]].x;
						pml_e->nodes[2].x = elem->nodes[aux[1]].x;
						pml_e->nodes[3].x = elem->nodes[aux[0]].x;
						pml_e->nodes[4].x = elem->nodes[aux[0]].x;
						pml_e->nodes[5].x = elem->nodes[aux[1]].x;
						pml_e->nodes[6].x = elem->nodes[aux[1]].x;
						pml_e->nodes[7].x = elem->nodes[aux[0]].x;

						y[0] = coords[3*node0+1] - ny*Y_pml/layers_y;
						y[1] = coords[3*node1+1] - ny*Y_pml/layers_y;
						y[2] = coords[3*node1+1] - (ny+1)*Y_pml/layers_y;
						y[3] = coords[3*node0+1] - (ny+1)*Y_pml/layers_y;
						y[4] = coords[3*node0+1] - ny*Y_pml/layers_y;
						y[5] = coords[3*node1+1] - ny*Y_pml/layers_y;
						y[6] = coords[3*node1+1] - (ny+1)*Y_pml/layers_y;
						y[7] = coords[3*node0+1] - (ny+1)*Y_pml/layers_y;

						pml_e->nodes[0].y = elem->nodes[aux[0]].y - 12*(ny+0);
						pml_e->nodes[1].y = elem->nodes[aux[1]].y - 12*(ny+0);
						pml_e->nodes[2].y = elem->nodes[aux[1]].y - 12*(ny+1);
						pml_e->nodes[3].y = elem->nodes[aux[0]].y - 12*(ny+1);
						pml_e->nodes[4].y = elem->nodes[aux[0]].y - 12*(ny+0);
						pml_e->nodes[5].y = elem->nodes[aux[1]].y - 12*(ny+0);
						pml_e->nodes[6].y = elem->nodes[aux[1]].y - 12*(ny+1);
						pml_e->nodes[7].y = elem->nodes[aux[0]].y - 12*(ny+1);

						z[0] = coords[3*node0+2] - nz*Z_pml/layers_z;
						z[1] = coords[3*node1+2] - nz*Z_pml/layers_z;
						z[2] = coords[3*node1+2] - nz*Z_pml/layers_z;
						z[3] = coords[3*node0+2] - nz*Z_pml/layers_z;
						z[4] = coords[3*node0+2] - (nz+1)*Z_pml/layers_z;
						z[5] = coords[3*node1+2] - (nz+1)*Z_pml/layers_z;
						z[6] = coords[3*node1+2] - (nz+1)*Z_pml/layers_z;
						z[7] = coords[3*node0+2] - (nz+1)*Z_pml/layers_z;

						pml_e->nodes[0].z = elem->nodes[aux[0]].z + 12*(nz+1);
						pml_e->nodes[1].z = elem->nodes[aux[1]].z + 12*(nz+1);
						pml_e->nodes[2].z = elem->nodes[aux[1]].z + 12*(nz+1);
						pml_e->nodes[3].z = elem->nodes[aux[0]].z + 12*(nz+1);
						pml_e->nodes[4].z = elem->nodes[aux[0]].z + 12*(nz+0);
						pml_e->nodes[5].z = elem->nodes[aux[1]].z + 12*(nz+0);
						pml_e->nodes[6].z = elem->nodes[aux[1]].z + 12*(nz+0);
						pml_e->nodes[7].z = elem->nodes[aux[0]].z + 12*(nz+0);

						for(int ino = 0; ino < 8 ; ino++){
							//definindo ponto p a ser adicionado
							GtsPoint* p = gts_point_new(gts_point_class(),x[ino], y[ino], z[ino]);
							//adicionando ponto p
							int x = pml_e->nodes[ino].x;
							int y = pml_e->nodes[ino].y;
							int z = pml_e->nodes[ino].z;
							pml_e->nodes[ino].id = AddPoint( mesh, hash_nodes, p, coords, x, y, z);
						}

						pml_e->n_mat = elem->n_mat+29;
					}
				}
			}

			if(elem->edge[9].ref){

				for(int nz = 0; nz < layers_z; nz++){
					for(int nx = 0; nx < layers_x; nx++){

						octant_t* pml_e = (octant_t*) sc_array_push(&mesh->elements);
						pml_e->id = mesh->elements.elem_count+1;

						//nos de referencia
						int aux[2] = {5,6};
						int node0 = elem->nodes[aux[0]].id;
						int node1 = elem->nodes[aux[1]].id;
						double x[8],y[8],z[8];

						x[0] = coords[3*node0+0] + nx*X_pml/layers_x;
						x[1] = coords[3*node0+0] + (nx+1)*X_pml/layers_x;
						x[2] = coords[3*node1+0] + (nx+1)*X_pml/layers_x;
						x[3] = coords[3*node1+0] + nx*X_pml/layers_x;
						x[4] = coords[3*node0+0] + nx*X_pml/layers_x;
						x[5] = coords[3*node0+0] + (nx+1)*X_pml/layers_x;
						x[6] = coords[3*node1+0] + (nx+1)*X_pml/layers_x;
						x[7] = coords[3*node1+0] + nx*X_pml/layers_x;

						pml_e->nodes[0].x = elem->nodes[aux[0]].x + 12*(nx+0);
						pml_e->nodes[1].x = elem->nodes[aux[0]].x + 12*(nx+1);
						pml_e->nodes[2].x = elem->nodes[aux[1]].x + 12*(nx+1);
						pml_e->nodes[3].x = elem->nodes[aux[1]].x + 12*(nx+0);
						pml_e->nodes[4].x = elem->nodes[aux[0]].x + 12*(nx+0);
						pml_e->nodes[5].x = elem->nodes[aux[0]].x + 12*(nx+1);
						pml_e->nodes[6].x = elem->nodes[aux[1]].x + 12*(nx+1);
						pml_e->nodes[7].x = elem->nodes[aux[1]].x + 12*(nx+0);

						y[0] = coords[3*node0+1];
						y[1] = coords[3*node0+1];
						y[2] = coords[3*node1+1];
						y[3] = coords[3*node1+1];
						y[4] = coords[3*node0+1];
						y[5] = coords[3*node0+1];
						y[6] = coords[3*node1+1];
						y[7] = coords[3*node1+1];

						pml_e->nodes[0].y = elem->nodes[aux[0]].y;
						pml_e->nodes[1].y = elem->nodes[aux[0]].y;
						pml_e->nodes[2].y = elem->nodes[aux[1]].y;
						pml_e->nodes[3].y = elem->nodes[aux[1]].y;
						pml_e->nodes[4].y = elem->nodes[aux[0]].y;
						pml_e->nodes[5].y = elem->nodes[aux[0]].y;
						pml_e->nodes[6].y = elem->nodes[aux[1]].y;
						pml_e->nodes[7].y = elem->nodes[aux[1]].y;

						z[0] = coords[3*node0+2] - nz*Z_pml/layers_z;
						z[1] = coords[3*node0+2] - nz*Z_pml/layers_z;
						z[2] = coords[3*node1+2] - nz*Z_pml/layers_z;
						z[3] = coords[3*node1+2] - nz*Z_pml/layers_z;
						z[4] = coords[3*node0+2] - (nz+1)*Z_pml/layers_z;
						z[5] = coords[3*node0+2] - (nz+1)*Z_pml/layers_z;
						z[6] = coords[3*node1+2] - (nz+1)*Z_pml/layers_z;
						z[7] = coords[3*node0+2] - (nz+1)*Z_pml/layers_z;

						pml_e->nodes[0].z = elem->nodes[aux[0]].z + 12*(nz+1);
						pml_e->nodes[1].z = elem->nodes[aux[0]].z + 12*(nz+1);
						pml_e->nodes[2].z = elem->nodes[aux[1]].z + 12*(nz+1);
						pml_e->nodes[3].z = elem->nodes[aux[1]].z + 12*(nz+1);
						pml_e->nodes[4].z = elem->nodes[aux[0]].z + 12*(nz+0);
						pml_e->nodes[5].z = elem->nodes[aux[0]].z + 12*(nz+0);
						pml_e->nodes[6].z = elem->nodes[aux[1]].z + 12*(nz+0);
						pml_e->nodes[7].z = elem->nodes[aux[1]].z + 12*(nz+0);

						for(int ino = 0; ino < 8 ; ino++){
							//definindo ponto p a ser adicionado
							GtsPoint* p = gts_point_new(gts_point_class(),x[ino], y[ino], z[ino]);
							//adicionando ponto p
							int x = pml_e->nodes[ino].x;
							int y = pml_e->nodes[ino].y;
							int z = pml_e->nodes[ino].z;
							pml_e->nodes[ino].id = AddPoint( mesh, hash_nodes, p, coords, x, y, z);
						}

						pml_e->n_mat = elem->n_mat+27;
					}
				}
			}

			if(elem->edge[10].ref){
				for(int nz = 0; nz < layers_z; nz++){
					for(int ny = 0; ny < layers_y; ny++){

						octant_t* pml_e = (octant_t*) sc_array_push(&mesh->elements);
						pml_e->id = mesh->elements.elem_count+1;

						//nos de referencia
						int aux[2] = {7,6};
						int node0 = elem->nodes[aux[0]].id;
						int node1 = elem->nodes[aux[1]].id;
						double x[8],y[8],z[8];

						x[0] = coords[3*node0+0];
						x[1] = coords[3*node1+0];
						x[2] = coords[3*node1+0];
						x[3] = coords[3*node0+0];
						x[4] = coords[3*node0+0];
						x[5] = coords[3*node1+0];
						x[6] = coords[3*node1+0];
						x[7] = coords[3*node0+0];

						pml_e->nodes[0].x = elem->nodes[aux[0]].x;
						pml_e->nodes[1].x = elem->nodes[aux[1]].x;
						pml_e->nodes[2].x = elem->nodes[aux[1]].x;
						pml_e->nodes[3].x = elem->nodes[aux[0]].x;
						pml_e->nodes[4].x = elem->nodes[aux[0]].x;
						pml_e->nodes[5].x = elem->nodes[aux[1]].x;
						pml_e->nodes[6].x = elem->nodes[aux[1]].x;
						pml_e->nodes[7].x = elem->nodes[aux[0]].x;

						y[0] = coords[3*node0+1] + ny*Y_pml/layers_y;
						y[1] = coords[3*node1+1] + ny*Y_pml/layers_y;
						y[2] = coords[3*node1+1] + (ny+1)*Y_pml/layers_y;
						y[3] = coords[3*node0+1] + (ny+1)*Y_pml/layers_y;
						y[4] = coords[3*node0+1] + ny*Y_pml/layers_y;
						y[5] = coords[3*node1+1] + ny*Y_pml/layers_y;
						y[6] = coords[3*node1+1] + (ny+1)*Y_pml/layers_y;
						y[7] = coords[3*node0+1] + (ny+1)*Y_pml/layers_y;

						pml_e->nodes[0].y = elem->nodes[aux[0]].y + 12*(ny+0);
						pml_e->nodes[1].y = elem->nodes[aux[1]].y + 12*(ny+0);
						pml_e->nodes[2].y = elem->nodes[aux[1]].y + 12*(ny+1);
						pml_e->nodes[3].y = elem->nodes[aux[0]].y + 12*(ny+1);
						pml_e->nodes[4].y = elem->nodes[aux[0]].y + 12*(ny+0);
						pml_e->nodes[5].y = elem->nodes[aux[1]].y + 12*(ny+0);
						pml_e->nodes[6].y = elem->nodes[aux[1]].y + 12*(ny+1);
						pml_e->nodes[7].y = elem->nodes[aux[0]].y + 12*(ny+1);

						z[0] = coords[3*node0+2] - nz*Z_pml/layers_z;
						z[1] = coords[3*node1+2] - nz*Z_pml/layers_z;
						z[2] = coords[3*node1+2] - nz*Z_pml/layers_z;
						z[3] = coords[3*node0+2] - nz*Z_pml/layers_z;
						z[4] = coords[3*node0+2] - (nz+1)*Z_pml/layers_z;
						z[5] = coords[3*node1+2] - (nz+1)*Z_pml/layers_z;
						z[6] = coords[3*node1+2] - (nz+1)*Z_pml/layers_z;
						z[7] = coords[3*node0+2] - (nz+1)*Z_pml/layers_z;

						pml_e->nodes[0].z = elem->nodes[aux[0]].z + 12*(nz+1);
						pml_e->nodes[1].z = elem->nodes[aux[1]].z + 12*(nz+1);
						pml_e->nodes[2].z = elem->nodes[aux[1]].z + 12*(nz+1);
						pml_e->nodes[3].z = elem->nodes[aux[0]].z + 12*(nz+1);
						pml_e->nodes[4].z = elem->nodes[aux[0]].z + 12*(nz+0);
						pml_e->nodes[5].z = elem->nodes[aux[1]].z + 12*(nz+0);
						pml_e->nodes[6].z = elem->nodes[aux[1]].z + 12*(nz+0);
						pml_e->nodes[7].z = elem->nodes[aux[0]].z + 12*(nz+0);

						for(int ino = 0; ino < 8 ; ino++){
							//definindo ponto p a ser adicionado
							GtsPoint* p = gts_point_new(gts_point_class(),x[ino], y[ino], z[ino]);
							//adicionando ponto p
							int x = pml_e->nodes[ino].x;
							int y = pml_e->nodes[ino].y;
							int z = pml_e->nodes[ino].z;
							pml_e->nodes[ino].id = AddPoint( mesh, hash_nodes, p, coords, x, y, z);
						}

						pml_e->n_mat = elem->n_mat+31;
					}
				}
			}

			if(elem->edge[11].ref){
				for(int nz = 0; nz < layers_z; nz++){
					for(int nx = 0; nx < layers_x; nx++){

						octant_t* pml_e = (octant_t*) sc_array_push(&mesh->elements);
						pml_e->id = mesh->elements.elem_count+1;

						//nos de referencia
						int aux[2] = {4,7};
						int node0 = elem->nodes[aux[0]].id;
						int node1 = elem->nodes[aux[1]].id;
						double x[8],y[8],z[8];

						x[0] = coords[3*node0+0] - nx*X_pml/layers_x;
						x[1] = coords[3*node0+0] - (nx+1)*X_pml/layers_x;
						x[2] = coords[3*node1+0] - (nx+1)*X_pml/layers_x;
						x[3] = coords[3*node1+0] - nx*X_pml/layers_x;
						x[4] = coords[3*node0+0] - nx*X_pml/layers_x;
						x[5] = coords[3*node0+0] - (nx+1)*X_pml/layers_x;
						x[6] = coords[3*node1+0] - (nx+1)*X_pml/layers_x;
						x[7] = coords[3*node1+0] - nx*X_pml/layers_x;

						pml_e->nodes[0].x = elem->nodes[aux[0]].x - 12*(nx+0);
						pml_e->nodes[1].x = elem->nodes[aux[0]].x - 12*(nx+1);
						pml_e->nodes[2].x = elem->nodes[aux[1]].x - 12*(nx+1);
						pml_e->nodes[3].x = elem->nodes[aux[1]].x - 12*(nx+0);
						pml_e->nodes[4].x = elem->nodes[aux[0]].x - 12*(nx+0);
						pml_e->nodes[5].x = elem->nodes[aux[0]].x - 12*(nx+1);
						pml_e->nodes[6].x = elem->nodes[aux[1]].x - 12*(nx+1);
						pml_e->nodes[7].x = elem->nodes[aux[1]].x - 12*(nx+0);

						y[0] = coords[3*node0+1];
						y[1] = coords[3*node0+1];
						y[2] = coords[3*node1+1];
						y[3] = coords[3*node1+1];
						y[4] = coords[3*node0+1];
						y[5] = coords[3*node0+1];
						y[6] = coords[3*node1+1];
						y[7] = coords[3*node1+1];

						pml_e->nodes[0].y = elem->nodes[aux[0]].y;
						pml_e->nodes[1].y = elem->nodes[aux[0]].y;
						pml_e->nodes[2].y = elem->nodes[aux[1]].y;
						pml_e->nodes[3].y = elem->nodes[aux[1]].y;
						pml_e->nodes[4].y = elem->nodes[aux[0]].y;
						pml_e->nodes[5].y = elem->nodes[aux[0]].y;
						pml_e->nodes[6].y = elem->nodes[aux[1]].y;
						pml_e->nodes[7].y = elem->nodes[aux[1]].y;

						z[0] = coords[3*node0+2] - nz*Z_pml/layers_z;
						z[1] = coords[3*node0+2] - nz*Z_pml/layers_z;
						z[2] = coords[3*node1+2] - nz*Z_pml/layers_z;
						z[3] = coords[3*node1+2] - nz*Z_pml/layers_z;
						z[4] = coords[3*node0+2] - (nz+1)*Z_pml/layers_z;
						z[5] = coords[3*node0+2] - (nz+1)*Z_pml/layers_z;
						z[6] = coords[3*node1+2] - (nz+1)*Z_pml/layers_z;
						z[7] = coords[3*node0+2] - (nz+1)*Z_pml/layers_z;

						pml_e->nodes[0].z = elem->nodes[aux[0]].z + 12*(nz+1);
						pml_e->nodes[1].z = elem->nodes[aux[0]].z + 12*(nz+1);
						pml_e->nodes[2].z = elem->nodes[aux[1]].z + 12*(nz+1);
						pml_e->nodes[3].z = elem->nodes[aux[1]].z + 12*(nz+1);
						pml_e->nodes[4].z = elem->nodes[aux[0]].z + 12*(nz+0);
						pml_e->nodes[5].z = elem->nodes[aux[0]].z + 12*(nz+0);
						pml_e->nodes[6].z = elem->nodes[aux[1]].z + 12*(nz+0);
						pml_e->nodes[7].z = elem->nodes[aux[1]].z + 12*(nz+0);

						for(int ino = 0; ino < 8 ; ino++){
							//definindo ponto p a ser adicionado
							GtsPoint* p = gts_point_new(gts_point_class(),x[ino], y[ino], z[ino]);
							//adicionando ponto p
							int x = pml_e->nodes[ino].x;
							int y = pml_e->nodes[ino].y;
							int z = pml_e->nodes[ino].z;
							pml_e->nodes[ino].id = AddPoint( mesh, hash_nodes, p, coords, x, y, z);
						}

						pml_e->n_mat = elem->n_mat+25;
					}
				}

			}

		}

		if(point){
			if(elem->nodes[0].fixed == -1  && false){

				for(int nx = 0; nx < layers_x; nx++){
					for(int ny = 0; ny < layers_y; ny++){
						for(int nz = 0; nz < layers_z; nz++){

							octant_t* pml_e = (octant_t*) sc_array_push(&mesh->elements);
							pml_e->id = mesh->elements.elem_count+1;

							//nos de referencia
							int aux = 0;
							int node0 = elem->nodes[aux].id;
							double x[8],y[8],z[8];

							x[0] = coords[3*node0+0] - nx*X_pml/layers_x;
							x[1] = coords[3*node0+0] - (nx+1)*X_pml/layers_x;
							x[2] = coords[3*node0+0] - (nx+1)*X_pml/layers_x;
							x[3] = coords[3*node0+0] - nx*X_pml/layers_x;
							x[4] = coords[3*node0+0] - nx*X_pml/layers_x;
							x[5] = coords[3*node0+0] - (nx+1)*X_pml/layers_x;
							x[6] = coords[3*node0+0] - (nx+1)*X_pml/layers_x;
							x[7] = coords[3*node0+0] - nx*X_pml/layers_x;

							pml_e->nodes[0].x = elem->nodes[aux].x - 12*(nx+0);
							pml_e->nodes[1].x = elem->nodes[aux].x - 12*(nx+1);
							pml_e->nodes[2].x = elem->nodes[aux].x - 12*(nx+1);
							pml_e->nodes[3].x = elem->nodes[aux].x - 12*(nx+0);
							pml_e->nodes[4].x = elem->nodes[aux].x - 12*(nx+0);
							pml_e->nodes[5].x = elem->nodes[aux].x - 12*(nx+1);
							pml_e->nodes[6].x = elem->nodes[aux].x - 12*(nx+1);
							pml_e->nodes[7].x = elem->nodes[aux].x - 12*(nx+0);

							y[0] = coords[3*node0+1] - ny*Y_pml/layers_y;
							y[1] = coords[3*node0+1] - ny*Y_pml/layers_y;
							y[2] = coords[3*node0+1] - (ny+1)*Y_pml/layers_y;
							y[3] = coords[3*node0+1] - (ny+1)*Y_pml/layers_y;
							y[4] = coords[3*node0+1] - ny*Y_pml/layers_y;
							y[5] = coords[3*node0+1] - ny*Y_pml/layers_y;
							y[6] = coords[3*node0+1] - (ny+1)*Y_pml/layers_y;
							y[7] = coords[3*node0+1] - (ny+1)*Y_pml/layers_y;

							pml_e->nodes[0].y = elem->nodes[aux].y - 12*(ny+0);
							pml_e->nodes[1].y = elem->nodes[aux].y - 12*(ny+0);
							pml_e->nodes[2].y = elem->nodes[aux].y - 12*(ny+1);
							pml_e->nodes[3].y = elem->nodes[aux].y - 12*(ny+1);
							pml_e->nodes[4].y = elem->nodes[aux].y - 12*(ny+0);
							pml_e->nodes[5].y = elem->nodes[aux].y - 12*(ny+0);
							pml_e->nodes[6].y = elem->nodes[aux].y - 12*(ny+1);
							pml_e->nodes[7].y = elem->nodes[aux].y - 12*(ny+1);

							z[0] = coords[3*node0+2] + nz*Z_pml/layers_z;
							z[1] = coords[3*node0+2] + nz*Z_pml/layers_z;
							z[2] = coords[3*node0+2] + nz*Z_pml/layers_z;
							z[3] = coords[3*node0+2] + nz*Z_pml/layers_z;
							z[4] = coords[3*node0+2] + (nz+1)*Z_pml/layers_z;
							z[5] = coords[3*node0+2] + (nz+1)*Z_pml/layers_z;
							z[6] = coords[3*node0+2] + (nz+1)*Z_pml/layers_z;
							z[7] = coords[3*node0+2] + (nz+1)*Z_pml/layers_z;

							pml_e->nodes[0].z = elem->nodes[aux].z - 12*(nz+1);
							pml_e->nodes[1].z = elem->nodes[aux].z - 12*(nz+1);
							pml_e->nodes[2].z = elem->nodes[aux].z - 12*(nz+1);
							pml_e->nodes[3].z = elem->nodes[aux].z - 12*(nz+1);
							pml_e->nodes[4].z = elem->nodes[aux].z - 12*(nz+0);
							pml_e->nodes[5].z = elem->nodes[aux].z - 12*(nz+0);
							pml_e->nodes[6].z = elem->nodes[aux].z - 12*(nz+0);
							pml_e->nodes[7].z = elem->nodes[aux].z - 12*(nz+0);

							for(int ino = 0; ino < 8 ; ino++){
								//definindo ponto p a ser adicionado
								GtsPoint* p = gts_point_new(gts_point_class(),x[ino], y[ino], z[ino]);
								//adicionando ponto p
								int x = pml_e->nodes[ino].x;
								int y = pml_e->nodes[ino].y;
								int z = pml_e->nodes[ino].z;
								pml_e->nodes[ino].id = AddPoint( mesh, hash_nodes, p, coords, x, y, z);
							}

							pml_e->n_mat = elem->n_mat+12;
						}
					}
				}
			}

			if(elem->nodes[1].fixed == -2 && false){

				for(int nx = 0; nx < layers_x; nx++){
					for(int ny = 0; ny < layers_y; ny++){
						for(int nz = 0; nz < layers_z; nz++){

							octant_t* pml_e = (octant_t*) sc_array_push(&mesh->elements);
							pml_e->id = mesh->elements.elem_count+1;

							//nos de referencia
							int aux = 1;
							int node0 = elem->nodes[5].id;
							double x[8],y[8],z[8];

							x[0] = coords[3*node0+0] + nx*X_pml/layers_x;
							x[1] = coords[3*node0+0] + (nx+1)*X_pml/layers_x;
							x[2] = coords[3*node0+0] + (nx+1)*X_pml/layers_x;
							x[3] = coords[3*node0+0] + nx*X_pml/layers_x;
							x[4] = coords[3*node0+0] + nx*X_pml/layers_x;
							x[5] = coords[3*node0+0] + (nx+1)*X_pml/layers_x;
							x[6] = coords[3*node0+0] + (nx+1)*X_pml/layers_x;
							x[7] = coords[3*node0+0] + nx*X_pml/layers_x;

							pml_e->nodes[0].x = elem->nodes[aux].x + 12*(nx+0);
							pml_e->nodes[1].x = elem->nodes[aux].x + 12*(nx+1);
							pml_e->nodes[2].x = elem->nodes[aux].x + 12*(nx+1);
							pml_e->nodes[3].x = elem->nodes[aux].x + 12*(nx+0);
							pml_e->nodes[4].x = elem->nodes[aux].x + 12*(nx+0);
							pml_e->nodes[5].x = elem->nodes[aux].x + 12*(nx+1);
							pml_e->nodes[6].x = elem->nodes[aux].x + 12*(nx+1);
							pml_e->nodes[7].x = elem->nodes[aux].x + 12*(nx+0);

							y[0] = coords[3*node0+1] - ny*Y_pml/layers_y;
							y[1] = coords[3*node0+1] - ny*Y_pml/layers_y;
							y[2] = coords[3*node0+1] - (ny+1)*Y_pml/layers_y;
							y[3] = coords[3*node0+1] - (ny+1)*Y_pml/layers_y;
							y[4] = coords[3*node0+1] - ny*Y_pml/layers_y;
							y[5] = coords[3*node0+1] - ny*Y_pml/layers_y;
							y[6] = coords[3*node0+1] - (ny+1)*Y_pml/layers_y;
							y[7] = coords[3*node0+1] - (ny+1)*Y_pml/layers_y;

							pml_e->nodes[0].y = elem->nodes[aux].y - 12*(ny+0);
							pml_e->nodes[1].y = elem->nodes[aux].y - 12*(ny+0);
							pml_e->nodes[2].y = elem->nodes[aux].y - 12*(ny+1);
							pml_e->nodes[3].y = elem->nodes[aux].y - 12*(ny+1);
							pml_e->nodes[4].y = elem->nodes[aux].y - 12*(ny+0);
							pml_e->nodes[5].y = elem->nodes[aux].y - 12*(ny+0);
							pml_e->nodes[6].y = elem->nodes[aux].y - 12*(ny+1);
							pml_e->nodes[7].y = elem->nodes[aux].y - 12*(ny+1);

							z[0] = coords[3*node0+2] + nz*Z_pml/layers_z;
							z[1] = coords[3*node0+2] + nz*Z_pml/layers_z;
							z[2] = coords[3*node0+2] + nz*Z_pml/layers_z;
							z[3] = coords[3*node0+2] + nz*Z_pml/layers_z;
							z[4] = coords[3*node0+2] + (nz+1)*Z_pml/layers_z;
							z[5] = coords[3*node0+2] + (nz+1)*Z_pml/layers_z;
							z[6] = coords[3*node0+2] + (nz+1)*Z_pml/layers_z;
							z[7] = coords[3*node0+2] + (nz+1)*Z_pml/layers_z;

							pml_e->nodes[0].z = elem->nodes[aux].z - 12*(nz+1);
							pml_e->nodes[1].z = elem->nodes[aux].z - 12*(nz+1);
							pml_e->nodes[2].z = elem->nodes[aux].z - 12*(nz+1);
							pml_e->nodes[3].z = elem->nodes[aux].z - 12*(nz+1);
							pml_e->nodes[4].z = elem->nodes[aux].z - 12*(nz+0);
							pml_e->nodes[5].z = elem->nodes[aux].z - 12*(nz+0);
							pml_e->nodes[6].z = elem->nodes[aux].z - 12*(nz+0);
							pml_e->nodes[7].z = elem->nodes[aux].z - 12*(nz+0);

							for(int ino = 0; ino < 8 ; ino++){
								//definindo ponto p a ser adicionado
								GtsPoint* p = gts_point_new(gts_point_class(),x[ino], y[ino], z[ino]);
								//adicionando ponto p
								int x = pml_e->nodes[ino].x;
								int y = pml_e->nodes[ino].y;
								int z = pml_e->nodes[ino].z;
								pml_e->nodes[ino].id = AddPoint( mesh, hash_nodes, p, coords, x, y, z);
							}

							pml_e->n_mat = elem->n_mat+13;
						}
					}
				}

			}

			if(elem->nodes[2].fixed == -3  && false){


				for(int nx = 0; nx < layers_x; nx++){
					for(int ny = 0; ny < layers_y; ny++){
						for(int nz = 0; nz < layers_z; nz++){

							octant_t* pml_e = (octant_t*) sc_array_push(&mesh->elements);
							pml_e->id = mesh->elements.elem_count+1;

							//nos de referencia
							int aux = 2;
							int node0 = elem->nodes[6].id;
							double x[8],y[8],z[8];

							x[0] = coords[3*node0+0] + nx*X_pml/layers_x;
							x[1] = coords[3*node0+0] + (nx+1)*X_pml/layers_x;
							x[2] = coords[3*node0+0] + (nx+1)*X_pml/layers_x;
							x[3] = coords[3*node0+0] + nx*X_pml/layers_x;
							x[4] = coords[3*node0+0] + nx*X_pml/layers_x;
							x[5] = coords[3*node0+0] + (nx+1)*X_pml/layers_x;
							x[6] = coords[3*node0+0] + (nx+1)*X_pml/layers_x;
							x[7] = coords[3*node0+0] + nx*X_pml/layers_x;

							pml_e->nodes[0].x = elem->nodes[aux].x + 12*(nx+0);
							pml_e->nodes[1].x = elem->nodes[aux].x + 12*(nx+1);
							pml_e->nodes[2].x = elem->nodes[aux].x + 12*(nx+1);
							pml_e->nodes[3].x = elem->nodes[aux].x + 12*(nx+0);
							pml_e->nodes[4].x = elem->nodes[aux].x + 12*(nx+0);
							pml_e->nodes[5].x = elem->nodes[aux].x + 12*(nx+1);
							pml_e->nodes[6].x = elem->nodes[aux].x + 12*(nx+1);
							pml_e->nodes[7].x = elem->nodes[aux].x + 12*(nx+0);

							y[0] = coords[3*node0+1] + ny*Y_pml/layers_y;
							y[1] = coords[3*node0+1] + ny*Y_pml/layers_y;
							y[2] = coords[3*node0+1] + (ny+1)*Y_pml/layers_y;
							y[3] = coords[3*node0+1] + (ny+1)*Y_pml/layers_y;
							y[4] = coords[3*node0+1] + ny*Y_pml/layers_y;
							y[5] = coords[3*node0+1] + ny*Y_pml/layers_y;
							y[6] = coords[3*node0+1] + (ny+1)*Y_pml/layers_y;
							y[7] = coords[3*node0+1] + (ny+1)*Y_pml/layers_y;

							pml_e->nodes[0].y = elem->nodes[aux].y + 12*(ny+0);
							pml_e->nodes[1].y = elem->nodes[aux].y + 12*(ny+0);
							pml_e->nodes[2].y = elem->nodes[aux].y + 12*(ny+1);
							pml_e->nodes[3].y = elem->nodes[aux].y + 12*(ny+1);
							pml_e->nodes[4].y = elem->nodes[aux].y + 12*(ny+0);
							pml_e->nodes[5].y = elem->nodes[aux].y + 12*(ny+0);
							pml_e->nodes[6].y = elem->nodes[aux].y + 12*(ny+1);
							pml_e->nodes[7].y = elem->nodes[aux].y + 12*(ny+1);

							z[0] = coords[3*node0+2] - nz*Z_pml/layers_z;
							z[1] = coords[3*node0+2] - nz*Z_pml/layers_z;
							z[2] = coords[3*node0+2] - nz*Z_pml/layers_z;
							z[3] = coords[3*node0+2] - nz*Z_pml/layers_z;
							z[4] = coords[3*node0+2] - (nz+1)*Z_pml/layers_z;
							z[5] = coords[3*node0+2] - (nz+1)*Z_pml/layers_z;
							z[6] = coords[3*node0+2] - (nz+1)*Z_pml/layers_z;
							z[7] = coords[3*node0+2] - (nz+1)*Z_pml/layers_z;

							z[0] = coords[3*node0+2] + nz*Z_pml/layers_z;
							z[1] = coords[3*node0+2] + nz*Z_pml/layers_z;
							z[2] = coords[3*node0+2] + nz*Z_pml/layers_z;
							z[3] = coords[3*node0+2] + nz*Z_pml/layers_z;
							z[4] = coords[3*node0+2] + (nz+1)*Z_pml/layers_z;
							z[5] = coords[3*node0+2] + (nz+1)*Z_pml/layers_z;
							z[6] = coords[3*node0+2] + (nz+1)*Z_pml/layers_z;
							z[7] = coords[3*node0+2] + (nz+1)*Z_pml/layers_z;

							pml_e->nodes[0].z = elem->nodes[aux].z - 12*(nz+1);
							pml_e->nodes[1].z = elem->nodes[aux].z - 12*(nz+1);
							pml_e->nodes[2].z = elem->nodes[aux].z - 12*(nz+1);
							pml_e->nodes[3].z = elem->nodes[aux].z - 12*(nz+1);
							pml_e->nodes[4].z = elem->nodes[aux].z - 12*(nz+0);
							pml_e->nodes[5].z = elem->nodes[aux].z - 12*(nz+0);
							pml_e->nodes[6].z = elem->nodes[aux].z - 12*(nz+0);
							pml_e->nodes[7].z = elem->nodes[aux].z - 12*(nz+0);

							for(int ino = 0; ino < 8 ; ino++){
								//definindo ponto p a ser adicionado
								GtsPoint* p = gts_point_new(gts_point_class(),x[ino], y[ino], z[ino]);
								//adicionando ponto p
								int x = pml_e->nodes[ino].x;
								int y = pml_e->nodes[ino].y;
								int z = pml_e->nodes[ino].z;
								pml_e->nodes[ino].id = AddPoint( mesh, hash_nodes, p, coords, x, y, z);
							}

							pml_e->n_mat = elem->n_mat+15;
						}
					}

				}


			}

			if(elem->nodes[3].fixed == -4  && false){

				for(int nx = 0; nx < layers_x; nx++){
					for(int ny = 0; ny < layers_y; ny++){
						for(int nz = 0; nz < layers_z; nz++){

							octant_t* pml_e = (octant_t*) sc_array_push(&mesh->elements);
							pml_e->id = mesh->elements.elem_count+1;

							//nos de referencia
							int aux = 3;
							int node0 = elem->nodes[aux].id;
							double x[8],y[8],z[8];

							x[0] = coords[3*node0+0] - nx*X_pml/layers_x;
							x[1] = coords[3*node0+0] - (nx+1)*X_pml/layers_x;
							x[2] = coords[3*node0+0] - (nx+1)*X_pml/layers_x;
							x[3] = coords[3*node0+0] - nx*X_pml/layers_x;
							x[4] = coords[3*node0+0] - nx*X_pml/layers_x;
							x[5] = coords[3*node0+0] - (nx+1)*X_pml/layers_x;
							x[6] = coords[3*node0+0] - (nx+1)*X_pml/layers_x;
							x[7] = coords[3*node0+0] - nx*X_pml/layers_x;

							pml_e->nodes[0].x = elem->nodes[aux].x - 12*(nx+0);
							pml_e->nodes[1].x = elem->nodes[aux].x - 12*(nx+1);
							pml_e->nodes[2].x = elem->nodes[aux].x - 12*(nx+1);
							pml_e->nodes[3].x = elem->nodes[aux].x - 12*(nx+0);
							pml_e->nodes[4].x = elem->nodes[aux].x - 12*(nx+0);
							pml_e->nodes[5].x = elem->nodes[aux].x - 12*(nx+1);
							pml_e->nodes[6].x = elem->nodes[aux].x - 12*(nx+1);
							pml_e->nodes[7].x = elem->nodes[aux].x - 12*(nx+0);

							y[0] = coords[3*node0+1] + ny*Y_pml/layers_y;
							y[1] = coords[3*node0+1] + ny*Y_pml/layers_y;
							y[2] = coords[3*node0+1] + (ny+1)*Y_pml/layers_y;
							y[3] = coords[3*node0+1] + (ny+1)*Y_pml/layers_y;
							y[4] = coords[3*node0+1] + ny*Y_pml/layers_y;
							y[5] = coords[3*node0+1] + ny*Y_pml/layers_y;
							y[6] = coords[3*node0+1] + (ny+1)*Y_pml/layers_y;
							y[7] = coords[3*node0+1] + (ny+1)*Y_pml/layers_y;

							pml_e->nodes[0].y = elem->nodes[aux].y + 12*(ny+0);
							pml_e->nodes[1].y = elem->nodes[aux].y + 12*(ny+0);
							pml_e->nodes[2].y = elem->nodes[aux].y + 12*(ny+1);
							pml_e->nodes[3].y = elem->nodes[aux].y + 12*(ny+1);
							pml_e->nodes[4].y = elem->nodes[aux].y + 12*(ny+0);
							pml_e->nodes[5].y = elem->nodes[aux].y + 12*(ny+0);
							pml_e->nodes[6].y = elem->nodes[aux].y + 12*(ny+1);
							pml_e->nodes[7].y = elem->nodes[aux].y + 12*(ny+1);

							z[0] = coords[3*node0+2] + nz*Z_pml/layers_z;
							z[1] = coords[3*node0+2] + nz*Z_pml/layers_z;
							z[2] = coords[3*node0+2] + nz*Z_pml/layers_z;
							z[3] = coords[3*node0+2] + nz*Z_pml/layers_z;
							z[4] = coords[3*node0+2] + (nz+1)*Z_pml/layers_z;
							z[5] = coords[3*node0+2] + (nz+1)*Z_pml/layers_z;
							z[6] = coords[3*node0+2] + (nz+1)*Z_pml/layers_z;
							z[7] = coords[3*node0+2] + (nz+1)*Z_pml/layers_z;

							pml_e->nodes[0].z = elem->nodes[aux].z - 12*(nz+1);
							pml_e->nodes[1].z = elem->nodes[aux].z - 12*(nz+1);
							pml_e->nodes[2].z = elem->nodes[aux].z - 12*(nz+1);
							pml_e->nodes[3].z = elem->nodes[aux].z - 12*(nz+1);
							pml_e->nodes[4].z = elem->nodes[aux].z - 12*(nz+0);
							pml_e->nodes[5].z = elem->nodes[aux].z - 12*(nz+0);
							pml_e->nodes[6].z = elem->nodes[aux].z - 12*(nz+0);
							pml_e->nodes[7].z = elem->nodes[aux].z - 12*(nz+0);

							for(int ino = 0; ino < 8 ; ino++){
								//definindo ponto p a ser adicionado
								GtsPoint* p = gts_point_new(gts_point_class(),x[ino], y[ino], z[ino]);
								//adicionando ponto p
								int x = pml_e->nodes[ino].x;
								int y = pml_e->nodes[ino].y;
								int z = pml_e->nodes[ino].z;
								pml_e->nodes[ino].id = AddPoint( mesh, hash_nodes, p, coords, x, y, z);
							}

							pml_e->n_mat = elem->n_mat+14;
						}
					}
				}


			}

			if(elem->nodes[4].fixed == -5){

				for(int nx = 0; nx < layers_x; nx++){
					for(int ny = 0; ny < layers_y; ny++){
						for(int nz = 0; nz < layers_z; nz++){

							octant_t* pml_e = (octant_t*) sc_array_push(&mesh->elements);
							pml_e->id = mesh->elements.elem_count+1;

							//nos de referencia
							int aux = 4;
							int node0 = elem->nodes[aux].id;
							double x[8],y[8],z[8];

							x[0] = coords[3*node0+0] - nx*X_pml/layers_x;
							x[1] = coords[3*node0+0] - (nx+1)*X_pml/layers_x;
							x[2] = coords[3*node0+0] - (nx+1)*X_pml/layers_x;
							x[3] = coords[3*node0+0] - nx*X_pml/layers_x;
							x[4] = coords[3*node0+0] - nx*X_pml/layers_x;
							x[5] = coords[3*node0+0] - (nx+1)*X_pml/layers_x;
							x[6] = coords[3*node0+0] - (nx+1)*X_pml/layers_x;
							x[7] = coords[3*node0+0] - nx*X_pml/layers_x;

							pml_e->nodes[0].x = elem->nodes[aux].x - 12*(nx+0);
							pml_e->nodes[1].x = elem->nodes[aux].x - 12*(nx+1);
							pml_e->nodes[2].x = elem->nodes[aux].x - 12*(nx+1);
							pml_e->nodes[3].x = elem->nodes[aux].x - 12*(nx+0);
							pml_e->nodes[4].x = elem->nodes[aux].x - 12*(nx+0);
							pml_e->nodes[5].x = elem->nodes[aux].x - 12*(nx+1);
							pml_e->nodes[6].x = elem->nodes[aux].x - 12*(nx+1);
							pml_e->nodes[7].x = elem->nodes[aux].x - 12*(nx+0);

							y[0] = coords[3*node0+1] - ny*Y_pml/layers_y;
							y[1] = coords[3*node0+1] - ny*Y_pml/layers_y;
							y[2] = coords[3*node0+1] - (ny+1)*Y_pml/layers_y;
							y[3] = coords[3*node0+1] - (ny+1)*Y_pml/layers_y;
							y[4] = coords[3*node0+1] - ny*Y_pml/layers_y;
							y[5] = coords[3*node0+1] - ny*Y_pml/layers_y;
							y[6] = coords[3*node0+1] - (ny+1)*Y_pml/layers_y;
							y[7] = coords[3*node0+1] - (ny+1)*Y_pml/layers_y;

							pml_e->nodes[0].y = elem->nodes[aux].y - 12*(ny+0);
							pml_e->nodes[1].y = elem->nodes[aux].y - 12*(ny+0);
							pml_e->nodes[2].y = elem->nodes[aux].y - 12*(ny+1);
							pml_e->nodes[3].y = elem->nodes[aux].y - 12*(ny+1);
							pml_e->nodes[4].y = elem->nodes[aux].y - 12*(ny+0);
							pml_e->nodes[5].y = elem->nodes[aux].y - 12*(ny+0);
							pml_e->nodes[6].y = elem->nodes[aux].y - 12*(ny+1);
							pml_e->nodes[7].y = elem->nodes[aux].y - 12*(ny+1);

							z[0] = coords[3*node0+2] - nz*Z_pml/layers_z;
							z[1] = coords[3*node0+2] - nz*Z_pml/layers_z;
							z[2] = coords[3*node0+2] - nz*Z_pml/layers_z;
							z[3] = coords[3*node0+2] - nz*Z_pml/layers_z;
							z[4] = coords[3*node0+2] - (nz+1)*Z_pml/layers_z;
							z[5] = coords[3*node0+2] - (nz+1)*Z_pml/layers_z;
							z[6] = coords[3*node0+2] - (nz+1)*Z_pml/layers_z;
							z[7] = coords[3*node0+2] - (nz+1)*Z_pml/layers_z;

							pml_e->nodes[0].z = elem->nodes[aux].z + 12*(nz+1);
							pml_e->nodes[1].z = elem->nodes[aux].z + 12*(nz+1);
							pml_e->nodes[2].z = elem->nodes[aux].z + 12*(nz+1);
							pml_e->nodes[3].z = elem->nodes[aux].z + 12*(nz+1);
							pml_e->nodes[4].z = elem->nodes[aux].z + 12*(nz+0);
							pml_e->nodes[5].z = elem->nodes[aux].z + 12*(nz+0);
							pml_e->nodes[6].z = elem->nodes[aux].z + 12*(nz+0);
							pml_e->nodes[7].z = elem->nodes[aux].z + 12*(nz+0);

							for(int ino = 0; ino < 8 ; ino++){
								//definindo ponto p a ser adicionado
								GtsPoint* p = gts_point_new(gts_point_class(),x[ino], y[ino], z[ino]);
								//adicionando ponto p
								int x = pml_e->nodes[ino].x;
								int y = pml_e->nodes[ino].y;
								int z = pml_e->nodes[ino].z;
								pml_e->nodes[ino].id = AddPoint( mesh, hash_nodes, p, coords, x, y, z);
							}

							pml_e->n_mat = elem->n_mat+12;
						}
					}
				}

			}

			if(elem->nodes[5].fixed == -6){

				for(int nx = 0; nx < layers_x; nx++){
					for(int ny = 0; ny < layers_y; ny++){
						for(int nz = 0; nz < layers_z; nz++){

							octant_t* pml_e = (octant_t*) sc_array_push(&mesh->elements);
							pml_e->id = mesh->elements.elem_count+1;

							//nos de referencia
							int aux = 5;
							int node0 = elem->nodes[5].id;
							double x[8],y[8],z[8];

							x[0] = coords[3*node0+0] + nx*X_pml/layers_x;
							x[1] = coords[3*node0+0] + (nx+1)*X_pml/layers_x;
							x[2] = coords[3*node0+0] + (nx+1)*X_pml/layers_x;
							x[3] = coords[3*node0+0] + nx*X_pml/layers_x;
							x[4] = coords[3*node0+0] + nx*X_pml/layers_x;
							x[5] = coords[3*node0+0] + (nx+1)*X_pml/layers_x;
							x[6] = coords[3*node0+0] + (nx+1)*X_pml/layers_x;
							x[7] = coords[3*node0+0] + nx*X_pml/layers_x;

							pml_e->nodes[0].x = elem->nodes[aux].x + 12*(nx+0);
							pml_e->nodes[1].x = elem->nodes[aux].x + 12*(nx+1);
							pml_e->nodes[2].x = elem->nodes[aux].x + 12*(nx+1);
							pml_e->nodes[3].x = elem->nodes[aux].x + 12*(nx+0);
							pml_e->nodes[4].x = elem->nodes[aux].x + 12*(nx+0);
							pml_e->nodes[5].x = elem->nodes[aux].x + 12*(nx+1);
							pml_e->nodes[6].x = elem->nodes[aux].x + 12*(nx+1);
							pml_e->nodes[7].x = elem->nodes[aux].x + 12*(nx+0);

							y[0] = coords[3*node0+1] - ny*Y_pml/layers_y;
							y[1] = coords[3*node0+1] - ny*Y_pml/layers_y;
							y[2] = coords[3*node0+1] - (ny+1)*Y_pml/layers_y;
							y[3] = coords[3*node0+1] - (ny+1)*Y_pml/layers_y;
							y[4] = coords[3*node0+1] - ny*Y_pml/layers_y;
							y[5] = coords[3*node0+1] - ny*Y_pml/layers_y;
							y[6] = coords[3*node0+1] - (ny+1)*Y_pml/layers_y;
							y[7] = coords[3*node0+1] - (ny+1)*Y_pml/layers_y;

							pml_e->nodes[0].y = elem->nodes[aux].y - 12*(ny+0);
							pml_e->nodes[1].y = elem->nodes[aux].y - 12*(ny+0);
							pml_e->nodes[2].y = elem->nodes[aux].y - 12*(ny+1);
							pml_e->nodes[3].y = elem->nodes[aux].y - 12*(ny+1);
							pml_e->nodes[4].y = elem->nodes[aux].y - 12*(ny+0);
							pml_e->nodes[5].y = elem->nodes[aux].y - 12*(ny+0);
							pml_e->nodes[6].y = elem->nodes[aux].y - 12*(ny+1);
							pml_e->nodes[7].y = elem->nodes[aux].y - 12*(ny+1);

							z[0] = coords[3*node0+2] - nz*Z_pml/layers_z;
							z[1] = coords[3*node0+2] - nz*Z_pml/layers_z;
							z[2] = coords[3*node0+2] - nz*Z_pml/layers_z;
							z[3] = coords[3*node0+2] - nz*Z_pml/layers_z;
							z[4] = coords[3*node0+2] - (nz+1)*Z_pml/layers_z;
							z[5] = coords[3*node0+2] - (nz+1)*Z_pml/layers_z;
							z[6] = coords[3*node0+2] - (nz+1)*Z_pml/layers_z;
							z[7] = coords[3*node0+2] - (nz+1)*Z_pml/layers_z;

							pml_e->nodes[0].z = elem->nodes[aux].z + 12*(nz+1);
							pml_e->nodes[1].z = elem->nodes[aux].z + 12*(nz+1);
							pml_e->nodes[2].z = elem->nodes[aux].z + 12*(nz+1);
							pml_e->nodes[3].z = elem->nodes[aux].z + 12*(nz+1);
							pml_e->nodes[4].z = elem->nodes[aux].z + 12*(nz+0);
							pml_e->nodes[5].z = elem->nodes[aux].z + 12*(nz+0);
							pml_e->nodes[6].z = elem->nodes[aux].z + 12*(nz+0);
							pml_e->nodes[7].z = elem->nodes[aux].z + 12*(nz+0);

							for(int ino = 0; ino < 8 ; ino++){
								//definindo ponto p a ser adicionado
								GtsPoint* p = gts_point_new(gts_point_class(),x[ino], y[ino], z[ino]);
								//adicionando ponto p
								int x = pml_e->nodes[ino].x;
								int y = pml_e->nodes[ino].y;
								int z = pml_e->nodes[ino].z;
								pml_e->nodes[ino].id = AddPoint( mesh, hash_nodes, p, coords, x, y, z);
							}

							pml_e->n_mat = elem->n_mat+13;
						}
					}
				}
			}

			if(elem->nodes[6].fixed == -7){

				for(int nx = 0; nx < layers_x; nx++){
					for(int ny = 0; ny < layers_y; ny++){
						for(int nz = 0; nz < layers_z; nz++){

							octant_t* pml_e = (octant_t*) sc_array_push(&mesh->elements);
							pml_e->id = mesh->elements.elem_count+1;

							//nos de referencia
							int aux = 6;
							int node0 = elem->nodes[6].id;
							double x[8],y[8],z[8];

							x[0] = coords[3*node0+0] + nx*X_pml/layers_x;
							x[1] = coords[3*node0+0] + (nx+1)*X_pml/layers_x;
							x[2] = coords[3*node0+0] + (nx+1)*X_pml/layers_x;
							x[3] = coords[3*node0+0] + nx*X_pml/layers_x;
							x[4] = coords[3*node0+0] + nx*X_pml/layers_x;
							x[5] = coords[3*node0+0] + (nx+1)*X_pml/layers_x;
							x[6] = coords[3*node0+0] + (nx+1)*X_pml/layers_x;
							x[7] = coords[3*node0+0] + nx*X_pml/layers_x;

							pml_e->nodes[0].x = elem->nodes[aux].x + 12*(nx+0);
							pml_e->nodes[1].x = elem->nodes[aux].x + 12*(nx+1);
							pml_e->nodes[2].x = elem->nodes[aux].x + 12*(nx+1);
							pml_e->nodes[3].x = elem->nodes[aux].x + 12*(nx+0);
							pml_e->nodes[4].x = elem->nodes[aux].x + 12*(nx+0);
							pml_e->nodes[5].x = elem->nodes[aux].x + 12*(nx+1);
							pml_e->nodes[6].x = elem->nodes[aux].x + 12*(nx+1);
							pml_e->nodes[7].x = elem->nodes[aux].x + 12*(nx+0);

							y[0] = coords[3*node0+1] + ny*Y_pml/layers_y;
							y[1] = coords[3*node0+1] + ny*Y_pml/layers_y;
							y[2] = coords[3*node0+1] + (ny+1)*Y_pml/layers_y;
							y[3] = coords[3*node0+1] + (ny+1)*Y_pml/layers_y;
							y[4] = coords[3*node0+1] + ny*Y_pml/layers_y;
							y[5] = coords[3*node0+1] + ny*Y_pml/layers_y;
							y[6] = coords[3*node0+1] + (ny+1)*Y_pml/layers_y;
							y[7] = coords[3*node0+1] + (ny+1)*Y_pml/layers_y;

							pml_e->nodes[0].y = elem->nodes[aux].y + 12*(ny+0);
							pml_e->nodes[1].y = elem->nodes[aux].y + 12*(ny+0);
							pml_e->nodes[2].y = elem->nodes[aux].y + 12*(ny+1);
							pml_e->nodes[3].y = elem->nodes[aux].y + 12*(ny+1);
							pml_e->nodes[4].y = elem->nodes[aux].y + 12*(ny+0);
							pml_e->nodes[5].y = elem->nodes[aux].y + 12*(ny+0);
							pml_e->nodes[6].y = elem->nodes[aux].y + 12*(ny+1);
							pml_e->nodes[7].y = elem->nodes[aux].y + 12*(ny+1);

							z[0] = coords[3*node0+2] - nz*Z_pml/layers_z;
							z[1] = coords[3*node0+2] - nz*Z_pml/layers_z;
							z[2] = coords[3*node0+2] - nz*Z_pml/layers_z;
							z[3] = coords[3*node0+2] - nz*Z_pml/layers_z;
							z[4] = coords[3*node0+2] - (nz+1)*Z_pml/layers_z;
							z[5] = coords[3*node0+2] - (nz+1)*Z_pml/layers_z;
							z[6] = coords[3*node0+2] - (nz+1)*Z_pml/layers_z;
							z[7] = coords[3*node0+2] - (nz+1)*Z_pml/layers_z;

							pml_e->nodes[0].z = elem->nodes[aux].z + 12*(nz+1);
							pml_e->nodes[1].z = elem->nodes[aux].z + 12*(nz+1);
							pml_e->nodes[2].z = elem->nodes[aux].z + 12*(nz+1);
							pml_e->nodes[3].z = elem->nodes[aux].z + 12*(nz+1);
							pml_e->nodes[4].z = elem->nodes[aux].z + 12*(nz+0);
							pml_e->nodes[5].z = elem->nodes[aux].z + 12*(nz+0);
							pml_e->nodes[6].z = elem->nodes[aux].z + 12*(nz+0);
							pml_e->nodes[7].z = elem->nodes[aux].z + 12*(nz+0);

							for(int ino = 0; ino < 8 ; ino++){
								//definindo ponto p a ser adicionado
								GtsPoint* p = gts_point_new(gts_point_class(),x[ino], y[ino], z[ino]);
								//adicionando ponto p
								int x = pml_e->nodes[ino].x;
								int y = pml_e->nodes[ino].y;
								int z = pml_e->nodes[ino].z;
								pml_e->nodes[ino].id = AddPoint( mesh, hash_nodes, p, coords, x, y, z);
							}

							pml_e->n_mat = elem->n_mat+15;
						}
					}

				}

			}

			if(elem->nodes[7].fixed == -8){

				for(int nx = 0; nx < layers_x; nx++){
					for(int ny = 0; ny < layers_y; ny++){
						for(int nz = 0; nz < layers_z; nz++){

							octant_t* pml_e = (octant_t*) sc_array_push(&mesh->elements);
							pml_e->id = mesh->elements.elem_count+1;

							//nos de referencia
							int aux = 7;
							int node0 = elem->nodes[aux].id;
							double x[8],y[8],z[8];

							x[0] = coords[3*node0+0] - nx*X_pml/layers_x;
							x[1] = coords[3*node0+0] - (nx+1)*X_pml/layers_x;
							x[2] = coords[3*node0+0] - (nx+1)*X_pml/layers_x;
							x[3] = coords[3*node0+0] - nx*X_pml/layers_x;
							x[4] = coords[3*node0+0] - nx*X_pml/layers_x;
							x[5] = coords[3*node0+0] - (nx+1)*X_pml/layers_x;
							x[6] = coords[3*node0+0] - (nx+1)*X_pml/layers_x;
							x[7] = coords[3*node0+0] - nx*X_pml/layers_x;

							pml_e->nodes[0].x = elem->nodes[aux].x - 12*(nx+0);
							pml_e->nodes[1].x = elem->nodes[aux].x - 12*(nx+1);
							pml_e->nodes[2].x = elem->nodes[aux].x - 12*(nx+1);
							pml_e->nodes[3].x = elem->nodes[aux].x - 12*(nx+0);
							pml_e->nodes[4].x = elem->nodes[aux].x - 12*(nx+0);
							pml_e->nodes[5].x = elem->nodes[aux].x - 12*(nx+1);
							pml_e->nodes[6].x = elem->nodes[aux].x - 12*(nx+1);
							pml_e->nodes[7].x = elem->nodes[aux].x - 12*(nx+0);

							y[0] = coords[3*node0+1] + ny*Y_pml/layers_y;
							y[1] = coords[3*node0+1] + ny*Y_pml/layers_y;
							y[2] = coords[3*node0+1] + (ny+1)*Y_pml/layers_y;
							y[3] = coords[3*node0+1] + (ny+1)*Y_pml/layers_y;
							y[4] = coords[3*node0+1] + ny*Y_pml/layers_y;
							y[5] = coords[3*node0+1] + ny*Y_pml/layers_y;
							y[6] = coords[3*node0+1] + (ny+1)*Y_pml/layers_y;
							y[7] = coords[3*node0+1] + (ny+1)*Y_pml/layers_y;

							pml_e->nodes[0].y = elem->nodes[aux].y + 12*(ny+0);
							pml_e->nodes[1].y = elem->nodes[aux].y + 12*(ny+0);
							pml_e->nodes[2].y = elem->nodes[aux].y + 12*(ny+1);
							pml_e->nodes[3].y = elem->nodes[aux].y + 12*(ny+1);
							pml_e->nodes[4].y = elem->nodes[aux].y + 12*(ny+0);
							pml_e->nodes[5].y = elem->nodes[aux].y + 12*(ny+0);
							pml_e->nodes[6].y = elem->nodes[aux].y + 12*(ny+1);
							pml_e->nodes[7].y = elem->nodes[aux].y + 12*(ny+1);

							z[0] = coords[3*node0+2] - nz*Z_pml/layers_z;
							z[1] = coords[3*node0+2] - nz*Z_pml/layers_z;
							z[2] = coords[3*node0+2] - nz*Z_pml/layers_z;
							z[3] = coords[3*node0+2] - nz*Z_pml/layers_z;
							z[4] = coords[3*node0+2] - (nz+1)*Z_pml/layers_z;
							z[5] = coords[3*node0+2] - (nz+1)*Z_pml/layers_z;
							z[6] = coords[3*node0+2] - (nz+1)*Z_pml/layers_z;
							z[7] = coords[3*node0+2] - (nz+1)*Z_pml/layers_z;

							pml_e->nodes[0].z = elem->nodes[aux].z + 12*(nz+1);
							pml_e->nodes[1].z = elem->nodes[aux].z + 12*(nz+1);
							pml_e->nodes[2].z = elem->nodes[aux].z + 12*(nz+1);
							pml_e->nodes[3].z = elem->nodes[aux].z + 12*(nz+1);
							pml_e->nodes[4].z = elem->nodes[aux].z + 12*(nz+0);
							pml_e->nodes[5].z = elem->nodes[aux].z + 12*(nz+0);
							pml_e->nodes[6].z = elem->nodes[aux].z + 12*(nz+0);
							pml_e->nodes[7].z = elem->nodes[aux].z + 12*(nz+0);

							for(int ino = 0; ino < 8 ; ino++){
								//definindo ponto p a ser adicionado
								GtsPoint* p = gts_point_new(gts_point_class(),x[ino], y[ino], z[ino]);
								//adicionando ponto p
								int x = pml_e->nodes[ino].x;
								int y = pml_e->nodes[ino].y;
								int z = pml_e->nodes[ino].z;
								pml_e->nodes[ino].id = AddPoint( mesh, hash_nodes, p, coords, x, y, z);
							}

							pml_e->n_mat = elem->n_mat+14;
						}
					}
				}

			}
		}

		sc_array_reset(&toto);
	}

	//update the vectors
	mesh->local_n_elements = mesh->elements.elem_count;
	mesh->local_n_nodes = mesh->nodes.elem_count;
	MPI_Allreduce(&mesh->local_n_elements, &mesh->total_n_elements, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(&mesh->local_n_nodes, &mesh->total_n_nodes, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);

	mesh->ncellx = mesh->ncellx + 2*4*layers_x;
	mesh->ncelly = mesh->ncelly + 2*4*layers_y;
	mesh->max_z = mesh->max_z + 4*layers_y;

	free(mesh->part_nodes);
	mesh->part_nodes = (int*) malloc (mesh->local_n_nodes*sizeof(int));
	for (int ino = 0; ino < mesh->local_n_nodes; ino++) {
		mesh->part_nodes[ino]=mesh->mpi_rank;
	}

	printf(" Ajust material properties\n\n");
	Adjust_material(mesh);

	if(mesh->mpi_rank == 0)
	{
		printf("Total number of elements: %d\n", mesh->total_n_elements);
		printf("Total number of nodes: %d\n", mesh->total_n_nodes);
	}
}
