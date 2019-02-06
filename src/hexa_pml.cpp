
#include <gts.h>
#include <gts.h>
#include <glib.h>
#include <vector>
#include <iostream>
using namespace std;
#include <sc.h>
#include <sc_io.h>
#include <sc_containers.h>

#include "hexa.h"
#include "pml.h"
#include "hilbert.h"

int8_t SetNodePML(hexa_tree_t* tree, octant_node_t* node) {
	int8_t pml_id = 0;
	if (node->x == 0) pml_id |= PML_X0;
	if (node->x == tree->ncellx) pml_id |= PML_X1;
	if (node->y == 0) pml_id |= PML_Y0;
	if (node->y == tree->ncelly) pml_id |= PML_Y1;
	//if(node->z == 0)               pml_id |= PML_Z0;
	if (node->z == tree->max_z) pml_id |= PML_Z1;
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
			(e1->coord[2] == e2->coord[2]) &&
			(e1->node_id == e2->node_id));
}


int AddPoint(hexa_tree_t* mesh, sc_hash_array_t* hash, GtsPoint *p, std::vector<double> &coords) {
	size_t position;
	node_t *r;
	node_t key;
	key.coord[0] = p->x;
	key.coord[1] = p->y;
	key.coord[2] = p->z;

	r = (node_t*) sc_hash_array_insert_unique(hash, &key, &position);

	if (r != NULL) {
		r->coord[0] = key.coord[0];
		r->coord[1] = key.coord[1];
		r->coord[2] = key.coord[2];
		r->node_id = mesh->nodes.elem_count;
		octant_node_t* n = (octant_node_t*) sc_array_push(&mesh->nodes);
		n->id = r->node_id;
		n->x = -1;
		n->y = -1;
		n->z = -1;
		//n->color = -1;
		n->fixed = 0;

		coords.push_back(p->x);
		coords.push_back(p->y);
		coords.push_back(p->z);
		return r->node_id;
	} else {
		r = (node_t*) sc_array_index(&hash->a, position);
		return r->node_id;
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

void ExtrudePMLElements(hexa_tree_t* mesh, std::vector<double>& coords) {

	double X_pml = 30000;
	double Y_pml = 30000;
	double Z_pml = 30000;
	//BUG here, if the number of layers is bigger than X the memory of elem cannot be read.
	int layers_x = 1;
	int layers_y = 1;
	int layers_z = 12;
	int mat_count = 25;
	int n_layers = 2;
	int8_t mask[NPML];
	int32_t npmls[NPML] = {0};

	bool clamped = true;
	sc_hash_array_t* hash_nodes = sc_hash_array_new(sizeof(node_t), edge_hash_fn, edge_equal_fn, &clamped);

	for(int n = 0;n<mesh->nodes.elem_count;n++){
		size_t position;
		node_t *r;
		node_t key;
		octant_node_t* node = (octant_node_t*) sc_array_index (&mesh->nodes, n);
		key.coord[0] = coords[3*node->id+0];
		key.coord[1] = coords[3*node->id+1];
		key.coord[2] = coords[3*node->id+2];
		key.node_id = node->id;

		r = (node_t*) sc_hash_array_insert_unique(hash_nodes, &key, &position);
		if(r!=NULL){
			r->coord[0] = coords[3*node->id+0];
			r->coord[1] = coords[3*node->id+1];
			r->coord[2] = coords[3*node->id+2];
			r->node_id = node->id;
		}else{
			printf("Verificar o no numero %d\n",node->id);
		}
	}

	int n_el=mesh->elements.elem_count;
	//printf("Numero de elementos %d\n",n_el);

	for (int i = 0; i < n_el; ++i) {
		octant_t* elem = (octant_t*) sc_array_index(&mesh->elements, i);
		elem->pml_id = 0;
		bool edge, face, point;
		point = true ;
		face = true;
		edge = true;

		SetElemPML(mesh, elem);
		SetPMLMask(mask, elem->pml_id);

		if(face){

			if(mask[PML_FACE_X0]){

				for(int n_l = 0; n_l < layers_x; ++n_l ){

					octant_t* pml_e = (octant_t*) sc_array_push(&mesh->elements);
					pml_e->id = mesh->elements.elem_count+1;

					//nos de referencia
					int node0 = elem->nodes[0].id;
					int node1 = elem->nodes[3].id;
					int node2 = elem->nodes[7].id;
					int node3 = elem->nodes[4].id;
					double x[8],y[8],z[8];

					x[0] = coords[3*node0+0]- n_l*X_pml/layers_x;
					x[1] = coords[3*node1+0]- (n_l+1)*X_pml/layers_x;
					x[2] = coords[3*node2+0]- (n_l+1)*X_pml/layers_x;
					x[3] = coords[3*node3+0]- n_l*X_pml/layers_x;
					x[4] = coords[3*node0+0] -n_l*X_pml/layers_x;
					x[5] = coords[3*node1+0]- (n_l+1)*X_pml/layers_x;
					x[6] = coords[3*node2+0]- (n_l+1)*X_pml/layers_x;
					x[7] = coords[3*node3+0]- n_l*X_pml/layers_x;

					y[0] = coords[3*node0+1];
					y[1] = coords[3*node0+1];
					y[2] = coords[3*node1+1];
					y[3] = coords[3*node1+1];
					y[4] = coords[3*node3+1];
					y[5] = coords[3*node3+1];
					y[6] = coords[3*node2+1];
					y[7] = coords[3*node2+1];

					z[0] = coords[3*node0+2];
					z[1] = coords[3*node0+2];
					z[2] = coords[3*node1+2];
					z[3] = coords[3*node1+2];
					z[4] = coords[3*node3+2];
					z[5] = coords[3*node3+2];
					z[6] = coords[3*node2+2];
					z[7] = coords[3*node2+2];

					for(int j = 0; j <8 ; j++){
						//definindo ponto p a ser adicionado
						GtsPoint p;
						gts_point_set(&p, x[j], y[j], z[j]);
						//adicionando ponto p
						pml_e->nodes[j].id = AddPoint(mesh, hash_nodes, &p, coords);
					}

					pml_e->n_mat = elem->n_mat+2;

				}

			}
			if(mask[PML_FACE_X1]){

				for(int n_l = 0; n_l < layers_x; ++n_l ){

					octant_t* pml_e = (octant_t*) sc_array_push(&mesh->elements);
					pml_e->id = mesh->elements.elem_count+1;

					//nos de referencia
					int node0 = elem->nodes[1].id;
					int node1 = elem->nodes[2].id;
					int node2 = elem->nodes[6].id;
					int node3 = elem->nodes[5].id;
					double x[8],y[8],z[8];

					x[0] = coords[3*node0+0]+ n_l*X_pml/layers_x;
					x[1] = coords[3*node1+0]+ (n_l+1)*X_pml/layers_x;
					x[2] = coords[3*node2+0]+ (n_l+1)*X_pml/layers_x;
					x[3] = coords[3*node3+0]+ n_l*X_pml/layers_x;
					x[4] = coords[3*node0+0]+ n_l*X_pml/layers_x;
					x[5] = coords[3*node1+0]+ (n_l+1)*X_pml/layers_x;
					x[6] = coords[3*node2+0]+ (n_l+1)*X_pml/layers_x;
					x[7] = coords[3*node3+0]+ n_l*X_pml/layers_x;

					y[0] = coords[3*node0+1];
					y[1] = coords[3*node0+1];
					y[2] = coords[3*node1+1];
					y[3] = coords[3*node1+1];
					y[4] = coords[3*node3+1];
					y[5] = coords[3*node3+1];
					y[6] = coords[3*node2+1];
					y[7] = coords[3*node2+1];

					z[0] = coords[3*node0+2];
					z[1] = coords[3*node0+2];
					z[2] = coords[3*node1+2];
					z[3] = coords[3*node1+2];
					z[4] = coords[3*node3+2];
					z[5] = coords[3*node3+2];
					z[6] = coords[3*node2+2];
					z[7] = coords[3*node2+2];

					for(int j = 0; j <8 ; j++){
						//definindo ponto p a ser adicionado
						GtsPoint p;
						gts_point_set(&p, x[j], y[j], z[j]);
						//adicionando ponto p
						pml_e->nodes[j].id = AddPoint(mesh, hash_nodes, &p, coords);
					}

					pml_e->n_mat = elem->n_mat+3;

				}
			}
			if(mask[PML_FACE_Y0]){

				for(int n_l = 0; n_l < layers_y; ++n_l ){

					octant_t* pml_e = (octant_t*) sc_array_push(&mesh->elements);
					pml_e->id = mesh->elements.elem_count+1;

					//nos de referencia
					int node0 = elem->nodes[0].id;
					int node1 = elem->nodes[1].id;
					int node2 = elem->nodes[4].id;
					int node3 = elem->nodes[5].id;
					double x[8],y[8],z[8];

					x[0] = coords[3*node0+0];
					x[1] = coords[3*node0+0];
					x[2] = coords[3*node1+0];
					x[3] = coords[3*node1+0];
					x[4] = coords[3*node2+0];
					x[5] = coords[3*node2+0];
					x[6] = coords[3*node3+0];
					x[7] = coords[3*node3+0];

					y[0] = coords[3*node0+1]- n_l*Y_pml/layers_y;
					y[1] = coords[3*node1+1]- (n_l+1)*Y_pml/layers_y;
					y[2] = coords[3*node2+1]- (n_l+1)*Y_pml/layers_y;
					y[3] = coords[3*node3+1]- n_l*Y_pml/layers_y;
					y[4] = coords[3*node0+1] -n_l*Y_pml/layers_y;
					y[5] = coords[3*node1+1]- (n_l+1)*Y_pml/layers_y;
					y[6] = coords[3*node2+1]- (n_l+1)*Y_pml/layers_y;
					y[7] = coords[3*node3+1]- n_l*Y_pml/layers_y;

					z[0] = coords[3*node0+2];
					z[1] = coords[3*node0+2];
					z[2] = coords[3*node1+2];
					z[3] = coords[3*node1+2];
					z[4] = coords[3*node2+2];
					z[5] = coords[3*node2+2];
					z[6] = coords[3*node3+2];
					z[7] = coords[3*node3+2];

					for(int j = 0; j <8 ; j++){
						//definindo ponto p a ser adicionado
						GtsPoint p;
						gts_point_set(&p, x[j], y[j], z[j]);
						//adicionando ponto p
						pml_e->nodes[j].id = AddPoint(mesh, hash_nodes, &p, coords);
					}

					pml_e->n_mat = elem->n_mat+5;

				}
			}
			if(mask[PML_FACE_Y1]){

				for(int n_l = 0; n_l < layers_y; ++n_l ){

					octant_t* pml_e = (octant_t*) sc_array_push(&mesh->elements);
					pml_e->id = mesh->elements.elem_count+1;

					//nos de referencia
					int node0 = elem->nodes[3].id;
					int node1 = elem->nodes[2].id;
					int node2 = elem->nodes[7].id;
					int node3 = elem->nodes[6].id;
					double x[8],y[8],z[8];

					x[0] = coords[3*node0+0];
					x[1] = coords[3*node0+0];
					x[2] = coords[3*node1+0];
					x[3] = coords[3*node1+0];
					x[4] = coords[3*node2+0];
					x[5] = coords[3*node2+0];
					x[6] = coords[3*node3+0];
					x[7] = coords[3*node3+0];

					y[0] = coords[3*node0+1]+ n_l*Y_pml/layers_y;
					y[1] = coords[3*node1+1]+ (n_l+1)*Y_pml/layers_y;
					y[2] = coords[3*node2+1]+ (n_l+1)*Y_pml/layers_y;
					y[3] = coords[3*node3+1]+ n_l*Y_pml/layers_y;
					y[4] = coords[3*node0+1]+ n_l*Y_pml/layers_y;
					y[5] = coords[3*node1+1]+ (n_l+1)*Y_pml/layers_y;
					y[6] = coords[3*node2+1]+ (n_l+1)*Y_pml/layers_y;
					y[7] = coords[3*node3+1]+ n_l*Y_pml/layers_y;

					z[0] = coords[3*node0+2];
					z[1] = coords[3*node0+2];
					z[2] = coords[3*node1+2];
					z[3] = coords[3*node1+2];
					z[4] = coords[3*node2+2];
					z[5] = coords[3*node2+2];
					z[6] = coords[3*node3+2];
					z[7] = coords[3*node3+2];

					for(int j = 0; j <8 ; j++){
						//definindo ponto p a ser adicionado
						GtsPoint p;
						gts_point_set(&p, x[j], y[j], z[j]);
						//adicionando ponto p
						pml_e->nodes[j].id = AddPoint(mesh, hash_nodes, &p, coords);
					}

					pml_e->n_mat = elem->n_mat+7;

				}
			}
			if(mask[PML_FACE_Z0]){

				for(int n_l = 0; n_l < layers_z; ++n_l ){

					octant_t* pml_e = (octant_t*) sc_array_push(&mesh->elements);
					pml_e->id = mesh->elements.elem_count+1;

					//int node_id;
					//nos de referencia
					int node0 = elem->nodes[0].id;
					int node1 = elem->nodes[1].id;
					int node2 = elem->nodes[2].id;
					int node3 = elem->nodes[3].id;
					double x[8],y[8],z[8];

					x[0] = coords[3*node0+0];
					x[1] = coords[3*node1+0];
					x[2] = coords[3*node2+0];
					x[3] = coords[3*node3+0];
					x[4] = coords[3*node0+0];
					x[5] = coords[3*node1+0];
					x[6] = coords[3*node2+0];
					x[7] = coords[3*node3+0];

					y[0] = coords[3*node0+1];
					y[1] = coords[3*node1+1];
					y[2] = coords[3*node2+1];
					y[3] = coords[3*node3+1];
					y[4] = coords[3*node0+1];
					y[5] = coords[3*node1+1];
					y[6] = coords[3*node2+1];
					y[7] = coords[3*node3+1];

					z[0] = coords[3*node0+2] + n_l*Z_pml/layers_z;
					z[1] = coords[3*node1+2] + n_l*Z_pml/layers_z;
					z[2] = coords[3*node2+2] + n_l*Z_pml/layers_z;
					z[3] = coords[3*node3+2] + n_l*Z_pml/layers_z;
					z[4] = coords[3*node0+2] + (n_l+1)*Z_pml/layers_z;
					z[5] = coords[3*node1+2] + (n_l+1)*Z_pml/layers_z;
					z[6] = coords[3*node2+2] + (n_l+1)*Z_pml/layers_z;
					z[7] = coords[3*node3+2] + (n_l+1)*Z_pml/layers_z;

					for(int j = 0; j <8 ; j++){
						//definindo ponto p a ser adicionado
						GtsPoint p;
						gts_point_set(&p, x[j], y[j], z[j]);
						//adicionando ponto p
						pml_e->nodes[j].id = AddPoint(mesh, hash_nodes, &p, coords);
					}

					pml_e->n_mat = elem->n_mat+9;

				}



			}
			if(mask[PML_FACE_Z1]){

				for(int n_l = 0; n_l < layers_z; ++n_l ){

					octant_t* pml_e = (octant_t*) sc_array_push(&mesh->elements);
					pml_e->id = mesh->elements.elem_count+1;

					//int node_id;
					//nos de referencia
					int node0 = elem->nodes[4].id;
					int node1 = elem->nodes[5].id;
					int node2 = elem->nodes[6].id;
					int node3 = elem->nodes[7].id;
					double x[8],y[8],z[8];

					x[0] = coords[3*node0+0];
					x[1] = coords[3*node1+0];
					x[2] = coords[3*node2+0];
					x[3] = coords[3*node3+0];
					x[4] = coords[3*node0+0];
					x[5] = coords[3*node1+0];
					x[6] = coords[3*node2+0];
					x[7] = coords[3*node3+0];

					y[0] = coords[3*node0+1];
					y[1] = coords[3*node1+1];
					y[2] = coords[3*node2+1];
					y[3] = coords[3*node3+1];
					y[4] = coords[3*node0+1];
					y[5] = coords[3*node1+1];
					y[6] = coords[3*node2+1];
					y[7] = coords[3*node3+1];

					z[0] = coords[3*node0+2] - n_l*Z_pml/layers_z;
					z[1] = coords[3*node1+2] - n_l*Z_pml/layers_z;
					z[2] = coords[3*node2+2] - n_l*Z_pml/layers_z;
					z[3] = coords[3*node3+2] - n_l*Z_pml/layers_z;
					z[4] = coords[3*node0+2] - (n_l+1)*Z_pml/layers_z;
					z[5] = coords[3*node1+2] - (n_l+1)*Z_pml/layers_z;
					z[6] = coords[3*node2+2] - (n_l+1)*Z_pml/layers_z;
					z[7] = coords[3*node3+2] - (n_l+1)*Z_pml/layers_z;

					for(int j = 0; j <8 ; j++){
						//definindo ponto p a ser adicionado
						GtsPoint p;
						gts_point_set(&p, x[j], y[j], z[j]);
						//adicionando ponto p
						pml_e->nodes[j].id = AddPoint(mesh, hash_nodes, &p, coords);
					}

					pml_e->n_mat = elem->n_mat+11;

				}
			}
		}

		if(point){
			if(mask[PML_CORNER_X0Y0Z0]){

			}
			if(mask[PML_CORNER_X1Y0Z0]){

			}
			if(mask[PML_CORNER_X0Y1Z0]){

			}
			if(mask[PML_CORNER_X1Y1Z0]){

			}
			if(mask[PML_CORNER_X0Y0Z1]){

				for(int nx = 0; nx < layers_x; nx++){
					for(int ny = 0; ny < layers_y; ny++){
						for(int nz = 0; nz < layers_z; nz++){

							octant_t* pml_e = (octant_t*) sc_array_push(&mesh->elements);
							pml_e->id = mesh->elements.elem_count+1;

							//nos de referencia
							int node0 = elem->nodes[4].id;
							double x[8],y[8],z[8];

							x[0] = coords[3*node0+0] - nx*X_pml/layers_x;
							x[1] = coords[3*node0+0] - (nx+1)*X_pml/layers_x;
							x[2] = coords[3*node0+0] - (nx+1)*X_pml/layers_x;
							x[3] = coords[3*node0+0] - nx*X_pml/layers_x;
							x[4] = coords[3*node0+0] - nx*X_pml/layers_x;
							x[5] = coords[3*node0+0] - (nx+1)*X_pml/layers_x;
							x[6] = coords[3*node0+0] - (nx+1)*X_pml/layers_x;
							x[7] = coords[3*node0+0] - nx*X_pml/layers_x;

							y[0] = coords[3*node0+1] - ny*Y_pml/layers_y;
							y[1] = coords[3*node0+1] - ny*Y_pml/layers_y;
							y[2] = coords[3*node0+1] - (ny+1)*Y_pml/layers_y;
							y[3] = coords[3*node0+1] - (ny+1)*Y_pml/layers_y;
							y[4] = coords[3*node0+1] - ny*Y_pml/layers_y;
							y[5] = coords[3*node0+1] - ny*Y_pml/layers_y;
							y[6] = coords[3*node0+1] - (ny+1)*Y_pml/layers_y;
							y[7] = coords[3*node0+1] - (ny+1)*Y_pml/layers_y;

							z[0] = coords[3*node0+2] - nz*Z_pml/layers_z;
							z[1] = coords[3*node0+2] - nz*Z_pml/layers_z;
							z[2] = coords[3*node0+2] - nz*Z_pml/layers_z;
							z[3] = coords[3*node0+2] - nz*Z_pml/layers_z;
							z[4] = coords[3*node0+2] - (nz+1)*Z_pml/layers_z;
							z[5] = coords[3*node0+2] - (nz+1)*Z_pml/layers_z;
							z[6] = coords[3*node0+2] - (nz+1)*Z_pml/layers_z;
							z[7] = coords[3*node0+2] - (nz+1)*Z_pml/layers_z;

							for(int j = 0; j <8 ; j++){
								//definindo ponto p a ser adicionado
								GtsPoint p;
								gts_point_set(&p, x[j], y[j], z[j]);
								//adicionando ponto p
								pml_e->nodes[j].id = AddPoint(mesh, hash_nodes, &p, coords);
							}

							pml_e->n_mat = elem->n_mat+12;
						}
					}
				}

			}
			if(mask[PML_CORNER_X1Y0Z1]){


				for(int nx = 0; nx < layers_x; nx++){
					for(int ny = 0; ny < layers_y; ny++){
						for(int nz = 0; nz < layers_z; nz++){

							octant_t* pml_e = (octant_t*) sc_array_push(&mesh->elements);
							pml_e->id = mesh->elements.elem_count+1;

							//nos de referencia
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

							y[0] = coords[3*node0+1] - ny*Y_pml/layers_y;
							y[1] = coords[3*node0+1] - ny*Y_pml/layers_y;
							y[2] = coords[3*node0+1] - (ny+1)*Y_pml/layers_y;
							y[3] = coords[3*node0+1] - (ny+1)*Y_pml/layers_y;
							y[4] = coords[3*node0+1] - ny*Y_pml/layers_y;
							y[5] = coords[3*node0+1] - ny*Y_pml/layers_y;
							y[6] = coords[3*node0+1] - (ny+1)*Y_pml/layers_y;
							y[7] = coords[3*node0+1] - (ny+1)*Y_pml/layers_y;

							z[0] = coords[3*node0+2] - nz*Z_pml/layers_z;
							z[1] = coords[3*node0+2] - nz*Z_pml/layers_z;
							z[2] = coords[3*node0+2] - nz*Z_pml/layers_z;
							z[3] = coords[3*node0+2] - nz*Z_pml/layers_z;
							z[4] = coords[3*node0+2] - (nz+1)*Z_pml/layers_z;
							z[5] = coords[3*node0+2] - (nz+1)*Z_pml/layers_z;
							z[6] = coords[3*node0+2] - (nz+1)*Z_pml/layers_z;
							z[7] = coords[3*node0+2] - (nz+1)*Z_pml/layers_z;

							for(int j = 0; j <8 ; j++){
								//definindo ponto p a ser adicionado
								GtsPoint p;
								gts_point_set(&p, x[j], y[j], z[j]);
								//adicionando ponto p
								pml_e->nodes[j].id = AddPoint(mesh, hash_nodes, &p, coords);
							}

							pml_e->n_mat = elem->n_mat+13;
						}
					}
				}


			}
			if(mask[PML_CORNER_X0Y1Z1]){

				for(int nx = 0; nx < layers_x; nx++){
					for(int ny = 0; ny < layers_y; ny++){
						for(int nz = 0; nz < layers_z; nz++){

							octant_t* pml_e = (octant_t*) sc_array_push(&mesh->elements);
							pml_e->id = mesh->elements.elem_count+1;

							//nos de referencia
							int node0 = elem->nodes[7].id;
							double x[8],y[8],z[8];

							x[0] = coords[3*node0+0] - nx*X_pml/layers_x;
							x[1] = coords[3*node0+0] - (nx+1)*X_pml/layers_x;
							x[2] = coords[3*node0+0] - (nx+1)*X_pml/layers_x;
							x[3] = coords[3*node0+0] - nx*X_pml/layers_x;
							x[4] = coords[3*node0+0] - nx*X_pml/layers_x;
							x[5] = coords[3*node0+0] - (nx+1)*X_pml/layers_x;
							x[6] = coords[3*node0+0] - (nx+1)*X_pml/layers_x;
							x[7] = coords[3*node0+0] - nx*X_pml/layers_x;

							y[0] = coords[3*node0+1] + ny*Y_pml/layers_y;
							y[1] = coords[3*node0+1] + ny*Y_pml/layers_y;
							y[2] = coords[3*node0+1] + (ny+1)*Y_pml/layers_y;
							y[3] = coords[3*node0+1] + (ny+1)*Y_pml/layers_y;
							y[4] = coords[3*node0+1] + ny*Y_pml/layers_y;
							y[5] = coords[3*node0+1] + ny*Y_pml/layers_y;
							y[6] = coords[3*node0+1] + (ny+1)*Y_pml/layers_y;
							y[7] = coords[3*node0+1] + (ny+1)*Y_pml/layers_y;

							z[0] = coords[3*node0+2] - nz*Z_pml/layers_z;
							z[1] = coords[3*node0+2] - nz*Z_pml/layers_z;
							z[2] = coords[3*node0+2] - nz*Z_pml/layers_z;
							z[3] = coords[3*node0+2] - nz*Z_pml/layers_z;
							z[4] = coords[3*node0+2] - (nz+1)*Z_pml/layers_z;
							z[5] = coords[3*node0+2] - (nz+1)*Z_pml/layers_z;
							z[6] = coords[3*node0+2] - (nz+1)*Z_pml/layers_z;
							z[7] = coords[3*node0+2] - (nz+1)*Z_pml/layers_z;

							for(int j = 0; j <8 ; j++){
								//definindo ponto p a ser adicionado
								GtsPoint p;
								gts_point_set(&p, x[j], y[j], z[j]);
								//adicionando ponto p
								pml_e->nodes[j].id = AddPoint(mesh, hash_nodes, &p, coords);
							}

							pml_e->n_mat = elem->n_mat+14;
						}
					}
				}


			}
			if(mask[PML_CORNER_X1Y1Z1]){

				for(int nx = 0; nx < layers_x; nx++){
					for(int ny = 0; ny < layers_y; ny++){
						for(int nz = 0; nz < layers_z; nz++){

							octant_t* pml_e = (octant_t*) sc_array_push(&mesh->elements);
							pml_e->id = mesh->elements.elem_count+1;

							//nos de referencia
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

							y[0] = coords[3*node0+1] + ny*Y_pml/layers_y;
							y[1] = coords[3*node0+1] + ny*Y_pml/layers_y;
							y[2] = coords[3*node0+1] + (ny+1)*Y_pml/layers_y;
							y[3] = coords[3*node0+1] + (ny+1)*Y_pml/layers_y;
							y[4] = coords[3*node0+1] + ny*Y_pml/layers_y;
							y[5] = coords[3*node0+1] + ny*Y_pml/layers_y;
							y[6] = coords[3*node0+1] + (ny+1)*Y_pml/layers_y;
							y[7] = coords[3*node0+1] + (ny+1)*Y_pml/layers_y;

							z[0] = coords[3*node0+2] - nz*Z_pml/layers_z;
							z[1] = coords[3*node0+2] - nz*Z_pml/layers_z;
							z[2] = coords[3*node0+2] - nz*Z_pml/layers_z;
							z[3] = coords[3*node0+2] - nz*Z_pml/layers_z;
							z[4] = coords[3*node0+2] - (nz+1)*Z_pml/layers_z;
							z[5] = coords[3*node0+2] - (nz+1)*Z_pml/layers_z;
							z[6] = coords[3*node0+2] - (nz+1)*Z_pml/layers_z;
							z[7] = coords[3*node0+2] - (nz+1)*Z_pml/layers_z;

							for(int j = 0; j <8 ; j++){
								//definindo ponto p a ser adicionado
								GtsPoint p;
								gts_point_set(&p, x[j], y[j], z[j]);
								//adicionando ponto p
								pml_e->nodes[j].id = AddPoint(mesh, hash_nodes, &p, coords);
							}

							pml_e->n_mat = elem->n_mat+15;
						}
					}
				}
			}
		}

		if(edge){
			if(mask[PML_EDGE_Z0_X0]){

			}
			if(mask[PML_EDGE_Z0_X1]){

			}
			if(mask[PML_EDGE_Z0_Y0]){

			}
			if(mask[PML_EDGE_Z0_Y1]){

			}
			if(mask[PML_EDGE_X0_Y0]){

				for(int nx = 0; nx < layers_x; nx++){
					for(int ny = 0; ny < layers_y; ny++){

						octant_t* pml_e = (octant_t*) sc_array_push(&mesh->elements);
						pml_e->id = mesh->elements.elem_count+1;

						//nos de referencia
						int node0 = elem->nodes[0].id;
						int node1 = elem->nodes[4].id;
						double x[8],y[8],z[8];

						x[0] = coords[3*node0+0] - nx*X_pml/layers_x;
						x[1] = coords[3*node0+0] - (nx+1)*X_pml/layers_x;
						x[2] = coords[3*node0+0] - (nx+1)*X_pml/layers_x;
						x[3] = coords[3*node0+0] - nx*X_pml/layers_x;
						x[4] = coords[3*node1+0] - nx*X_pml/layers_x;
						x[5] = coords[3*node1+0] - (nx+1)*X_pml/layers_x;
						x[6] = coords[3*node1+0] - (nx+1)*X_pml/layers_x;
						x[7] = coords[3*node1+0] - nx*X_pml/layers_x;

						y[0] = coords[3*node0+1] - ny*Y_pml/layers_y;
						y[1] = coords[3*node0+1] - ny*Y_pml/layers_y;
						y[2] = coords[3*node0+1] - (ny+1)*Y_pml/layers_y;
						y[3] = coords[3*node0+1] - (ny+1)*Y_pml/layers_y;
						y[4] = coords[3*node1+1] - ny*Y_pml/layers_y;
						y[5] = coords[3*node1+1] - ny*Y_pml/layers_y;
						y[6] = coords[3*node1+1] - (ny+1)*Y_pml/layers_y;
						y[7] = coords[3*node1+1] - (ny+1)*Y_pml/layers_y;

						z[0] = coords[3*node0+2];
						z[1] = coords[3*node0+2];
						z[2] = coords[3*node0+2];
						z[3] = coords[3*node0+2];
						z[4] = coords[3*node1+2];
						z[5] = coords[3*node1+2];
						z[6] = coords[3*node1+2];
						z[7] = coords[3*node1+2];

						for(int j = 0; j <8 ; j++){
							//definindo ponto p a ser adicionado
							GtsPoint p;
							gts_point_set(&p, x[j], y[j], z[j]);
							//adicionando ponto p
							pml_e->nodes[j].id = AddPoint(mesh, hash_nodes, &p, coords);
						}

						pml_e->n_mat = elem->n_mat+17;
					}
				}

			}
			if(mask[PML_EDGE_X0_Y1]){


				for(int nx = 0; nx < layers_x; nx++){
					for(int ny = 0; ny < layers_y; ny++){

						octant_t* pml_e = (octant_t*) sc_array_push(&mesh->elements);
						pml_e->id = mesh->elements.elem_count+1;

						//nos de referencia
						int node0 = elem->nodes[3].id;
						int node1 = elem->nodes[7].id;
						double x[8],y[8],z[8];

						x[0] = coords[3*node0+0] - nx*X_pml/layers_x;
						x[1] = coords[3*node0+0] - (nx+1)*X_pml/layers_x;
						x[2] = coords[3*node0+0] - (nx+1)*X_pml/layers_x;
						x[3] = coords[3*node0+0] - nx*X_pml/layers_x;
						x[4] = coords[3*node1+0] - nx*X_pml/layers_x;
						x[5] = coords[3*node1+0] - (nx+1)*X_pml/layers_x;
						x[6] = coords[3*node1+0] - (nx+1)*X_pml/layers_x;
						x[7] = coords[3*node1+0] - nx*X_pml/layers_x;

						y[0] = coords[3*node0+1] + ny*Y_pml/layers_y;
						y[1] = coords[3*node0+1] + ny*Y_pml/layers_y;
						y[2] = coords[3*node0+1] + (ny+1)*Y_pml/layers_y;
						y[3] = coords[3*node0+1] + (ny+1)*Y_pml/layers_y;
						y[4] = coords[3*node1+1] + ny*Y_pml/layers_y;
						y[5] = coords[3*node1+1] + ny*Y_pml/layers_y;
						y[6] = coords[3*node1+1] + (ny+1)*Y_pml/layers_y;
						y[7] = coords[3*node1+1] + (ny+1)*Y_pml/layers_y;

						z[0] = coords[3*node0+2];
						z[1] = coords[3*node0+2];
						z[2] = coords[3*node0+2];
						z[3] = coords[3*node0+2];
						z[4] = coords[3*node1+2];
						z[5] = coords[3*node1+2];
						z[6] = coords[3*node1+2];
						z[7] = coords[3*node1+2];

						for(int j = 0; j <8 ; j++){
							//definindo ponto p a ser adicionado
							GtsPoint p;
							gts_point_set(&p, x[j], y[j], z[j]);
							//adicionando ponto p
							pml_e->nodes[j].id = AddPoint(mesh, hash_nodes, &p, coords);
						}

						pml_e->n_mat = elem->n_mat+19;
					}
				}


			}
			if(mask[PML_EDGE_X1_Y0]){

				for(int nx = 0; nx < layers_x; nx++){
					for(int ny = 0; ny < layers_y; ny++){

						octant_t* pml_e = (octant_t*) sc_array_push(&mesh->elements);
						pml_e->id = mesh->elements.elem_count+1;

						//nos de referencia
						int node0 = elem->nodes[1].id;
						int node1 = elem->nodes[5].id;
						double x[8],y[8],z[8];

						x[0] = coords[3*node0+0] + nx*X_pml/layers_x;
						x[1] = coords[3*node0+0] + (nx+1)*X_pml/layers_x;
						x[2] = coords[3*node0+0] + (nx+1)*X_pml/layers_x;
						x[3] = coords[3*node0+0] + nx*X_pml/layers_x;
						x[4] = coords[3*node1+0] + nx*X_pml/layers_x;
						x[5] = coords[3*node1+0] + (nx+1)*X_pml/layers_x;
						x[6] = coords[3*node1+0] + (nx+1)*X_pml/layers_x;
						x[7] = coords[3*node1+0] + nx*X_pml/layers_x;

						y[0] = coords[3*node0+1] - ny*Y_pml/layers_y;
						y[1] = coords[3*node0+1] - ny*Y_pml/layers_y;
						y[2] = coords[3*node0+1] - (ny+1)*Y_pml/layers_y;
						y[3] = coords[3*node0+1] - (ny+1)*Y_pml/layers_y;
						y[4] = coords[3*node1+1] - ny*Y_pml/layers_y;
						y[5] = coords[3*node1+1] - ny*Y_pml/layers_y;
						y[6] = coords[3*node1+1] - (ny+1)*Y_pml/layers_y;
						y[7] = coords[3*node1+1] - (ny+1)*Y_pml/layers_y;

						z[0] = coords[3*node0+2];
						z[1] = coords[3*node0+2];
						z[2] = coords[3*node0+2];
						z[3] = coords[3*node0+2];
						z[4] = coords[3*node1+2];
						z[5] = coords[3*node1+2];
						z[6] = coords[3*node1+2];
						z[7] = coords[3*node1+2];

						for(int j = 0; j <8 ; j++){
							//definindo ponto p a ser adicionado
							GtsPoint p;
							gts_point_set(&p, x[j], y[j], z[j]);
							//adicionando ponto p
							pml_e->nodes[j].id = AddPoint(mesh, hash_nodes, &p, coords);
						}

						pml_e->n_mat = elem->n_mat+21;
					}
				}


			}
			if(mask[PML_EDGE_X1_Y1]){

				for(int nx = 0; nx < layers_x; nx++){
					for(int ny = 0; ny < layers_y; ny++){

						octant_t* pml_e = (octant_t*) sc_array_push(&mesh->elements);
						pml_e->id = mesh->elements.elem_count+1;

						//nos de referencia
						int node0 = elem->nodes[2].id;
						int node1 = elem->nodes[6].id;
						double x[8],y[8],z[8];

						x[0] = coords[3*node0+0] + nx*X_pml/layers_x;
						x[1] = coords[3*node0+0] + (nx+1)*X_pml/layers_x;
						x[2] = coords[3*node0+0] + (nx+1)*X_pml/layers_x;
						x[3] = coords[3*node0+0] + nx*X_pml/layers_x;
						x[4] = coords[3*node1+0] + nx*X_pml/layers_x;
						x[5] = coords[3*node1+0] + (nx+1)*X_pml/layers_x;
						x[6] = coords[3*node1+0] + (nx+1)*X_pml/layers_x;
						x[7] = coords[3*node1+0] + nx*X_pml/layers_x;

						y[0] = coords[3*node0+1] + ny*Y_pml/layers_y;
						y[1] = coords[3*node0+1] + ny*Y_pml/layers_y;
						y[2] = coords[3*node0+1] + (ny+1)*Y_pml/layers_y;
						y[3] = coords[3*node0+1] + (ny+1)*Y_pml/layers_y;
						y[4] = coords[3*node1+1] + ny*Y_pml/layers_y;
						y[5] = coords[3*node1+1] + ny*Y_pml/layers_y;
						y[6] = coords[3*node1+1] + (ny+1)*Y_pml/layers_y;
						y[7] = coords[3*node1+1] + (ny+1)*Y_pml/layers_y;

						z[0] = coords[3*node0+2];
						z[1] = coords[3*node0+2];
						z[2] = coords[3*node0+2];
						z[3] = coords[3*node0+2];
						z[4] = coords[3*node1+2];
						z[5] = coords[3*node1+2];
						z[6] = coords[3*node1+2];
						z[7] = coords[3*node1+2];

						for(int j = 0; j <8 ; j++){
							//definindo ponto p a ser adicionado
							GtsPoint p;
							gts_point_set(&p, x[j], y[j], z[j]);
							//adicionando ponto p
							pml_e->nodes[j].id = AddPoint(mesh, hash_nodes, &p, coords);
						}

						pml_e->n_mat = elem->n_mat+23;
					}
				}



			}
			if(mask[PML_EDGE_Z1_X0]){

				for(int nz = 0; nz < layers_z; nz++){
					for(int nx = 0; nx < layers_x; nx++){

						octant_t* pml_e = (octant_t*) sc_array_push(&mesh->elements);
						pml_e->id = mesh->elements.elem_count+1;

						//nos de referencia
						int node0 = elem->nodes[4].id;
						int node1 = elem->nodes[7].id;
						double x[8],y[8],z[8];

						x[0] = coords[3*node0+0] - nx*X_pml/layers_x;
						x[1] = coords[3*node0+0] - (nx+1)*X_pml/layers_x;
						x[2] = coords[3*node1+0] - (nx+1)*X_pml/layers_x;
						x[3] = coords[3*node1+0] - nx*X_pml/layers_x;
						x[4] = coords[3*node0+0] - nx*X_pml/layers_x;
						x[5] = coords[3*node0+0] - (nx+1)*X_pml/layers_x;
						x[6] = coords[3*node1+0] - (nx+1)*X_pml/layers_x;
						x[7] = coords[3*node1+0] - nx*X_pml/layers_x;

						y[0] = coords[3*node0+1];
						y[1] = coords[3*node0+1];
						y[2] = coords[3*node1+1];
						y[3] = coords[3*node1+1];
						y[4] = coords[3*node0+1];
						y[5] = coords[3*node0+1];
						y[6] = coords[3*node1+1];
						y[7] = coords[3*node1+1];

						z[0] = coords[3*node0+2] - nz*Z_pml/layers_z;
						z[1] = coords[3*node0+2] - nz*Z_pml/layers_z;
						z[2] = coords[3*node1+2] - nz*Z_pml/layers_z;
						z[3] = coords[3*node1+2] - nz*Z_pml/layers_z;
						z[4] = coords[3*node0+2] - (nz+1)*Z_pml/layers_z;
						z[5] = coords[3*node0+2] - (nz+1)*Z_pml/layers_z;
						z[6] = coords[3*node1+2] - (nz+1)*Z_pml/layers_z;
						z[7] = coords[3*node0+2] - (nz+1)*Z_pml/layers_z;

						for(int j = 0; j <8 ; j++){
							//definindo ponto p a ser adicionado
							GtsPoint p;
							gts_point_set(&p, x[j], y[j], z[j]);
							//adicionando ponto p
							pml_e->nodes[j].id = AddPoint(mesh, hash_nodes, &p, coords);
						}

						pml_e->n_mat = elem->n_mat+25;
					}
				}

			}
			if(mask[PML_EDGE_Z1_X1]){

				for(int nz = 0; nz < layers_z; nz++){
					for(int nx = 0; nx < layers_x; nx++){

						octant_t* pml_e = (octant_t*) sc_array_push(&mesh->elements);
						pml_e->id = mesh->elements.elem_count+1;

						//nos de referencia
						int node0 = elem->nodes[5].id;
						int node1 = elem->nodes[6].id;
						double x[8],y[8],z[8];

						x[0] = coords[3*node0+0] + nx*X_pml/layers_x;
						x[1] = coords[3*node0+0] + (nx+1)*X_pml/layers_x;
						x[2] = coords[3*node1+0] + (nx+1)*X_pml/layers_x;
						x[3] = coords[3*node1+0] + nx*X_pml/layers_x;
						x[4] = coords[3*node0+0] + nx*X_pml/layers_x;
						x[5] = coords[3*node0+0] + (nx+1)*X_pml/layers_x;
						x[6] = coords[3*node1+0] + (nx+1)*X_pml/layers_x;
						x[7] = coords[3*node1+0] + nx*X_pml/layers_x;

						y[0] = coords[3*node0+1];
						y[1] = coords[3*node0+1];
						y[2] = coords[3*node1+1];
						y[3] = coords[3*node1+1];
						y[4] = coords[3*node0+1];
						y[5] = coords[3*node0+1];
						y[6] = coords[3*node1+1];
						y[7] = coords[3*node1+1];

						z[0] = coords[3*node0+2] - nz*Z_pml/layers_z;
						z[1] = coords[3*node0+2] - nz*Z_pml/layers_z;
						z[2] = coords[3*node1+2] - nz*Z_pml/layers_z;
						z[3] = coords[3*node1+2] - nz*Z_pml/layers_z;
						z[4] = coords[3*node0+2] - (nz+1)*Z_pml/layers_z;
						z[5] = coords[3*node0+2] - (nz+1)*Z_pml/layers_z;
						z[6] = coords[3*node1+2] - (nz+1)*Z_pml/layers_z;
						z[7] = coords[3*node0+2] - (nz+1)*Z_pml/layers_z;

						for(int j = 0; j <8 ; j++){
							//definindo ponto p a ser adicionado
							GtsPoint p;
							gts_point_set(&p, x[j], y[j], z[j]);
							//adicionando ponto p
							pml_e->nodes[j].id = AddPoint(mesh, hash_nodes, &p, coords);
						}

						pml_e->n_mat = elem->n_mat+27;
					}
				}


			}
			if(mask[PML_EDGE_Z1_Y0]){

				for(int nz = 0; nz < layers_z; nz++){
					for(int ny = 0; ny < layers_y; ny++){

						octant_t* pml_e = (octant_t*) sc_array_push(&mesh->elements);
						pml_e->id = mesh->elements.elem_count+1;

						//nos de referencia
						int node0 = elem->nodes[4].id;
						int node1 = elem->nodes[5].id;
						double x[8],y[8],z[8];

						x[0] = coords[3*node0+0];
						x[1] = coords[3*node1+0];
						x[2] = coords[3*node1+0];
						x[3] = coords[3*node0+0];
						x[4] = coords[3*node0+0];
						x[5] = coords[3*node1+0];
						x[6] = coords[3*node1+0];
						x[7] = coords[3*node0+0];

						y[0] = coords[3*node0+1] - ny*Y_pml/layers_y;
						y[1] = coords[3*node1+1] - ny*Y_pml/layers_y;
						y[2] = coords[3*node1+1] - (ny+1)*Y_pml/layers_y;
						y[3] = coords[3*node0+1] - (ny+1)*Y_pml/layers_y;
						y[4] = coords[3*node0+1] - ny*Y_pml/layers_y;
						y[5] = coords[3*node1+1] - ny*Y_pml/layers_y;
						y[6] = coords[3*node1+1] - (ny+1)*Y_pml/layers_y;
						y[7] = coords[3*node0+1] - (ny+1)*Y_pml/layers_y;

						z[0] = coords[3*node0+2] - nz*Z_pml/layers_z;
						z[1] = coords[3*node1+2] - nz*Z_pml/layers_z;
						z[2] = coords[3*node1+2] - nz*Z_pml/layers_z;
						z[3] = coords[3*node0+2] - nz*Z_pml/layers_z;
						z[4] = coords[3*node0+2] - (nz+1)*Z_pml/layers_z;
						z[5] = coords[3*node1+2] - (nz+1)*Z_pml/layers_z;
						z[6] = coords[3*node1+2] - (nz+1)*Z_pml/layers_z;
						z[7] = coords[3*node0+2] - (nz+1)*Z_pml/layers_z;

						for(int j = 0; j <8 ; j++){
							//definindo ponto p a ser adicionado
							GtsPoint p;
							gts_point_set(&p, x[j], y[j], z[j]);
							//adicionando ponto p
							pml_e->nodes[j].id = AddPoint(mesh, hash_nodes, &p, coords);
						}

						pml_e->n_mat = elem->n_mat+29;
					}
				}
			}
			if(mask[PML_EDGE_Z1_Y1]){

				for(int nz = 0; nz < layers_z; nz++){
					for(int ny = 0; ny < layers_y; ny++){

						octant_t* pml_e = (octant_t*) sc_array_push(&mesh->elements);
						pml_e->id = mesh->elements.elem_count+1;

						//nos de referencia
						int node0 = elem->nodes[7].id;
						int node1 = elem->nodes[6].id;
						double x[8],y[8],z[8];

						x[0] = coords[3*node0+0];
						x[1] = coords[3*node1+0];
						x[2] = coords[3*node1+0];
						x[3] = coords[3*node0+0];
						x[4] = coords[3*node0+0];
						x[5] = coords[3*node1+0];
						x[6] = coords[3*node1+0];
						x[7] = coords[3*node0+0];

						y[0] = coords[3*node0+1] + ny*Y_pml/layers_y;
						y[1] = coords[3*node1+1] + ny*Y_pml/layers_y;
						y[2] = coords[3*node1+1] + (ny+1)*Y_pml/layers_y;
						y[3] = coords[3*node0+1] + (ny+1)*Y_pml/layers_y;
						y[4] = coords[3*node0+1] + ny*Y_pml/layers_y;
						y[5] = coords[3*node1+1] + ny*Y_pml/layers_y;
						y[6] = coords[3*node1+1] + (ny+1)*Y_pml/layers_y;
						y[7] = coords[3*node0+1] + (ny+1)*Y_pml/layers_y;

						z[0] = coords[3*node0+2] - nz*Z_pml/layers_z;
						z[1] = coords[3*node1+2] - nz*Z_pml/layers_z;
						z[2] = coords[3*node1+2] - nz*Z_pml/layers_z;
						z[3] = coords[3*node0+2] - nz*Z_pml/layers_z;
						z[4] = coords[3*node0+2] - (nz+1)*Z_pml/layers_z;
						z[5] = coords[3*node1+2] - (nz+1)*Z_pml/layers_z;
						z[6] = coords[3*node1+2] - (nz+1)*Z_pml/layers_z;
						z[7] = coords[3*node0+2] - (nz+1)*Z_pml/layers_z;

						for(int j = 0; j <8 ; j++){
							//definindo ponto p a ser adicionado
							GtsPoint p;
							gts_point_set(&p, x[j], y[j], z[j]);
							//adicionando ponto p
							pml_e->nodes[j].id = AddPoint(mesh, hash_nodes, &p, coords);
						}

						pml_e->n_mat = elem->n_mat+31;
					}
				}

			}
		}
	}

	//update the vectors
	mesh->local_n_elements = mesh->elements.elem_count;
	mesh->local_n_nodes = mesh->nodes.elem_count;
	MPI_Allreduce(&mesh->local_n_elements, &mesh->total_n_elements, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(&mesh->local_n_nodes, &mesh->total_n_nodes, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);

	printf(" Ajust material properties\n\n");
	Adjust_material(mesh);

}
