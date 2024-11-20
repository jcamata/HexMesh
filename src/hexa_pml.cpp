
//#include <gts.h>
//#include <gts.h>
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

unsigned edge_hash_fn(const void *v, const void *u)
{
	const node_t *q = (const node_t *)v;
	uint64_t a, b, c;

	a = (double_t)q->coord[0];
	b = (double_t)q->coord[1];
	c = (double_t)q->coord[2];

	sc_hash_mix(a, b, c);
	sc_hash_final(a, b, c);
	return (unsigned)c;
}

int edge_equal_fn(const void *v, const void *u, const void *w)
{
	const node_t *e1 = (const node_t *)v;
	const node_t *e2 = (const node_t *)u;

	return (unsigned)((e1->coord[0] == e2->coord[0]) &&
					  (e1->coord[1] == e2->coord[1]) &&
					  (e1->coord[2] == e2->coord[2]));
}

unsigned pml_hash_fn(const void *v, const void *u)
{
	const pmlmat_t *q = (const pmlmat_t *)v;
	uint32_t a, b, c;

	a = (uint32_t)q->id;
	b = (uint32_t)q->tag;
	c = (uint32_t)q->matref;

	sc_hash_mix(a, b, c);
	sc_hash_final(a, b, c);
	return (unsigned)c;
}

int pml_equal_fn(const void *v, const void *u, const void *w)
{
	const pmlmat_t *e1 = (const pmlmat_t *)v;
	const pmlmat_t *e2 = (const pmlmat_t *)u;

	return (unsigned)((e1->matref == e2->matref) &&
					  (e1->id == e2->id) &&
					  (e1->tag == e2->tag));
}

void RedoMap(hexa_tree_t *mesh, int layers_x, int layers_y, int layers_z)
{
	// add the number of PML layers in the mesh
	// it allow us add int points in the mesh
	// keeping a structured mesh
	for (int iel = 0; iel < mesh->elements.elem_count; iel++)
	{
		octant_t *elem = (octant_t *)sc_array_index(&mesh->elements, iel);
		for (int ino = 0; ino < 8; ino++)
		{
			elem->nodes[ino].x = elem->nodes[ino].x + 12 * layers_x;
			elem->nodes[ino].y = elem->nodes[ino].y + 12 * layers_y;
			elem->nodes[ino].z = elem->nodes[ino].z + 12 * layers_z;
		}
		elem->x = 4 * elem->x + 4 * layers_x;
		elem->y = 4 * elem->y + 4 * layers_y;
		elem->z = 4 * elem->z + 4 * layers_z;
	}

	for (int iel = 0; iel < mesh->outsurf.elem_count; iel++)
	{
		octant_t *elem = (octant_t *)sc_array_index(&mesh->outsurf, iel);
		for (int ino = 0; ino < 8; ino++)
		{
			elem->nodes[ino].x = elem->nodes[ino].x + 12 * layers_x;
			elem->nodes[ino].y = elem->nodes[ino].y + 12 * layers_y;
			elem->nodes[ino].z = elem->nodes[ino].z + 12 * layers_z;
		}
		elem->x = 4 * elem->x + 4 * layers_x;
		elem->y = 4 * elem->y + 4 * layers_y;
		elem->z = 4 * elem->z + 4 * layers_z;
	}
	for (int ino = 0; ino < mesh->nodes.elem_count; ino++)
	{
		octant_node_t *node = (octant_node_t *)sc_array_index(&mesh->nodes, ino);
		node->x = node->x + 12 * layers_x;
		node->y = node->y + 12 * layers_y;
		;
		node->z = node->z + 12 * layers_z;
	}
}

void ExtrudePMLElements(hexa_tree_t *mesh, std::vector<double> &coords)
{

	const double X_pml = 10e3;
	const double Y_pml = 10e3;
	const double Z_pml = 10e3;

	const int layers_x = 2;
	const int layers_y = 2;
	const int layers_z = 2;

	// material.input file 2 SEM3D
	FILE *fp;
	fp = fopen("DebugPML.txt", "w");
	if (fp == NULL)
	{
		printf("Error opening PML file\n");
	}

	// I should create a toto sc_array
	// it avoid segmentation fault when we perform a
	// push in mesh->elements sc_array due to the
	// realocation of the sc_array
	sc_array_t toto;
	sc_array_init(&toto, sizeof(octant_t));

	RedoMap(mesh, layers_x, layers_y, layers_z);

	bool clamped = true;
	sc_hash_array_t *hash_nodes = (sc_hash_array_t *)sc_hash_array_new(sizeof(octant_node_t), node_hash_fn, node_equal_fn, &clamped);

	double xinit = coords[0];
	double yinit = coords[1];
	for (int ino = 0; ino < mesh->nodes.elem_count; ino++)
	{
		size_t position;
		octant_node_t *r;
		octant_node_t key;
		octant_node_t *node = (octant_node_t *)sc_array_index(&mesh->nodes, ino);
		key.x = node->x;
		key.y = node->y;
		key.z = node->z;
		r = (octant_node_t *)sc_hash_array_insert_unique(hash_nodes, &key, &position);
		if (r != NULL)
		{
			r->x = node->x;
			r->y = node->y;
			r->z = node->z;
			r->id = node->id;
		}
	    coords[3*ino + 0] =  coords[3*ino + 0] - xinit;
		coords[3*ino + 1] =  coords[3*ino + 1] - yinit;
	}

	assert(hash_nodes->a.elem_count == mesh->nodes.elem_count);
	bool edge, face, point;
	point = true;
	face = true;
	edge = true;

	sc_hash_array_t *hash_matpml = (sc_hash_array_t *)sc_hash_array_new(sizeof(pmlmat_t), pml_hash_fn, pml_equal_fn, &clamped);

	int min_n_mat = 100;
	int tot_n_mat = 0;
	for (int iel = 0; iel < mesh->elements.elem_count; ++iel) {
		octant_t *elem = (octant_t*) sc_array_index(&mesh->elements, iel);
		if (elem->n_mat>tot_n_mat){
			tot_n_mat = elem->n_mat;
		}
		if (elem->n_mat<min_n_mat){
			min_n_mat = elem->n_mat;
		}
	}

	int MatCount = tot_n_mat+1;
	for (int iel = 0; iel < mesh->outsurf.elem_count; ++iel)
	{
		octant_t *elem = (octant_t *)sc_array_index(&mesh->outsurf, iel);

		std::vector<int> face;
		std::vector<int> edge;
		std::vector<int> corner;

		for (int isurf = 0; isurf < 6; isurf++)
		{
			if (elem->surf[isurf].ext && isurf != 4)
			{
				face.push_back(isurf + 1);
			}
		}

		for (int iedge = 4; iedge < 12; iedge++)
		{
			if (elem->edge[iedge].ref)
			{
				edge.push_back(10 * (iedge + 1));
			}
		}

		for (int ino = 4; ino < 8; ino++)
		{
			if (elem->nodes[ino].color < -29)
			{
				corner.push_back(1000 * (ino + 1));
			}
		}

		size_t position;
		pmlmat_t *r;
		pmlmat_t key;
		int tag = 0;
		int dir = 0;
		for (int isurf = 0; isurf < face.size(); isurf++)
		{
			key.id = face[isurf];
			key.tag = 0;
			key.matref = elem->n_mat;
			r = (pmlmat_t *)sc_hash_array_insert_unique(hash_matpml, &key, &position);
			if (r != NULL)
			{
				r->id = face[isurf];
				r->tag = 0;
				r->matref = elem->n_mat;
				r->mat = MatCount;
				MatCount++;
				r->xmax = -1e10;
				r->xmin =  1e10;
				r->ymax = -1e10;
				r->ymin =  1e10;
				r->zmax = -1e10;
				r->zmin =  1e10;
			}
		}

		for (int iedge = 0; iedge < edge.size(); iedge++)
		{
			key.id = edge[iedge];
			key.tag = 1;
			key.matref = elem->n_mat;
			r = (pmlmat_t *)sc_hash_array_insert_unique(hash_matpml, &key, &position);
			if (r != NULL)
			{
				r->id = edge[iedge];
				r->tag = 1;
				r->matref = elem->n_mat;
				r->mat = MatCount;
				MatCount++;
				r->xmax = -1e10;
				r->xmin =  1e10;
				r->ymax = -1e10;
				r->ymin =  1e10;
				r->zmax = -1e10;
				r->zmin =  1e10;
			}
		}

		for (int ino = 0; ino < corner.size(); ino++)
		{
			key.id = corner[ino];
			key.tag = 2;
			key.matref = elem->n_mat;
			r = (pmlmat_t *)sc_hash_array_insert_unique(hash_matpml, &key, &position);
			if (r != NULL)
			{
				r->id = corner[ino];
				r->tag = 2;
				r->matref = elem->n_mat;
				r->mat = MatCount;
				MatCount++;
				r->xmax = -1e10;
				r->xmin =  1e10;
				r->ymax = -1e10;
				r->ymin =  1e10;
				r->zmax = -1e10;
				r->zmin =  1e10;
			}
		}
	}

	for (int i = 0; i < mesh->outsurf.elem_count; ++i)
	{
		octant_t *elemOrig = (octant_t *)sc_array_index(&mesh->outsurf, i);
		octant_t *elem = (octant_t *)sc_array_push(&toto);

		hexa_element_copy(elemOrig, elem);

		// check the material
		int dir = 0;
		int tag = 0;
		int isurf = 0;
		int iedge = 0;
		int icorner = 0;

		size_t position;
		pmlmat_t key;
		key.tag = tag;
		key.matref = elemOrig->n_mat;

		if (face)
		{
			if (elem->surf[0].ext)
			{
				for (int n_l = 0; n_l < layers_x; ++n_l)
				{

					octant_t *pml_e = (octant_t *)sc_array_push(&mesh->elements);
					pml_e->id = mesh->elements.elem_count + 1;

					// nos de referencia
					int aux[4] = {0, 3, 7, 4};
					int node0 = elem->nodes[aux[0]].id;
					int node1 = elem->nodes[aux[1]].id;
					int node2 = elem->nodes[aux[2]].id;
					int node3 = elem->nodes[aux[3]].id;
					double x[8], y[8], z[8];

					x[0] = coords[3 * node0 + 0] - n_l * X_pml / layers_x;
					x[1] = coords[3 * node1 + 0] - (n_l + 1) * X_pml / layers_x;
					x[2] = coords[3 * node2 + 0] - (n_l + 1) * X_pml / layers_x;
					x[3] = coords[3 * node3 + 0] - n_l * X_pml / layers_x;
					x[4] = coords[3 * node0 + 0] - n_l * X_pml / layers_x;
					x[5] = coords[3 * node1 + 0] - (n_l + 1) * X_pml / layers_x;
					x[6] = coords[3 * node2 + 0] - (n_l + 1) * X_pml / layers_x;
					x[7] = coords[3 * node3 + 0] - n_l * X_pml / layers_x;

					pml_e->nodes[0].x = elem->nodes[aux[0]].x - 12 * (n_l + 0);
					pml_e->nodes[1].x = elem->nodes[aux[1]].x - 12 * (n_l + 1);
					pml_e->nodes[2].x = elem->nodes[aux[2]].x - 12 * (n_l + 1);
					pml_e->nodes[3].x = elem->nodes[aux[3]].x - 12 * (n_l + 0);
					pml_e->nodes[4].x = elem->nodes[aux[0]].x - 12 * (n_l + 0);
					pml_e->nodes[5].x = elem->nodes[aux[1]].x - 12 * (n_l + 1);
					pml_e->nodes[6].x = elem->nodes[aux[2]].x - 12 * (n_l + 1);
					pml_e->nodes[7].x = elem->nodes[aux[3]].x - 12 * (n_l + 0);

					y[0] = coords[3 * node0 + 1];
					y[1] = coords[3 * node0 + 1];
					y[2] = coords[3 * node1 + 1];
					y[3] = coords[3 * node1 + 1];
					y[4] = coords[3 * node3 + 1];
					y[5] = coords[3 * node3 + 1];
					y[6] = coords[3 * node2 + 1];
					y[7] = coords[3 * node2 + 1];

					pml_e->nodes[0].y = elem->nodes[aux[0]].y;
					pml_e->nodes[1].y = elem->nodes[aux[0]].y;
					pml_e->nodes[2].y = elem->nodes[aux[1]].y;
					pml_e->nodes[3].y = elem->nodes[aux[1]].y;
					pml_e->nodes[4].y = elem->nodes[aux[3]].y;
					pml_e->nodes[5].y = elem->nodes[aux[3]].y;
					pml_e->nodes[6].y = elem->nodes[aux[2]].y;
					pml_e->nodes[7].y = elem->nodes[aux[2]].y;

					z[0] = coords[3 * node0 + 2];
					z[1] = coords[3 * node0 + 2];
					z[2] = coords[3 * node1 + 2];
					z[3] = coords[3 * node1 + 2];
					z[4] = coords[3 * node3 + 2];
					z[5] = coords[3 * node3 + 2];
					z[6] = coords[3 * node2 + 2];
					z[7] = coords[3 * node2 + 2];

					pml_e->nodes[0].z = elem->nodes[aux[0]].z;
					pml_e->nodes[1].z = elem->nodes[aux[0]].z;
					pml_e->nodes[2].z = elem->nodes[aux[1]].z;
					pml_e->nodes[3].z = elem->nodes[aux[1]].z;
					pml_e->nodes[4].z = elem->nodes[aux[3]].z;
					pml_e->nodes[5].z = elem->nodes[aux[3]].z;
					pml_e->nodes[6].z = elem->nodes[aux[2]].z;
					pml_e->nodes[7].z = elem->nodes[aux[2]].z;

					for (int ino = 0; ino < 8; ino++)
					{
						// definindo ponto p a ser adicionado
						GtsPoint *p = gts_point_new(gts_point_class(), x[ino], y[ino], z[ino]);
						// adicionando ponto p
						int x = pml_e->nodes[ino].x;
						int y = pml_e->nodes[ino].y;
						int z = pml_e->nodes[ino].z;
						pml_e->nodes[ino].id = AddPoint(mesh, hash_nodes, p, coords, x, y, z);
					}
					isurf = 0;
					key.id = isurf + 1;
					bool lmat = sc_hash_array_lookup(hash_matpml, &key, &position);
					pmlmat_t *pmlT = (pmlmat_t *)sc_array_index(&hash_matpml->a, position);
					pml_e->n_mat = pmlT->mat;
					pmlT->xmin = std::min(pmlT->xmin, x[1]);
					pmlT->xmax = std::max(pmlT->xmax, x[0]);
					pmlT->ymin = 0;
					pmlT->ymax = 0;
					pmlT->zmin = 0;
					pmlT->zmax = 0;
				}
			}

			if (elem->surf[1].ext)
			{
				for (int n_l = 0; n_l < layers_x; ++n_l)
				{

					octant_t *pml_e = (octant_t *)sc_array_push(&mesh->elements);
					pml_e->id = mesh->elements.elem_count + 1;

					// nos de referencia
					int aux[4] = {1, 2, 6, 5};
					int node0 = elem->nodes[aux[0]].id;
					int node1 = elem->nodes[aux[1]].id;
					int node2 = elem->nodes[aux[2]].id;
					int node3 = elem->nodes[aux[3]].id;
					double x[8], y[8], z[8];

					x[0] = coords[3 * node0 + 0] + n_l * X_pml / layers_x;
					x[1] = coords[3 * node1 + 0] + (n_l + 1) * X_pml / layers_x;
					x[2] = coords[3 * node2 + 0] + (n_l + 1) * X_pml / layers_x;
					x[3] = coords[3 * node3 + 0] + n_l * X_pml / layers_x;
					x[4] = coords[3 * node0 + 0] + n_l * X_pml / layers_x;
					x[5] = coords[3 * node1 + 0] + (n_l + 1) * X_pml / layers_x;
					x[6] = coords[3 * node2 + 0] + (n_l + 1) * X_pml / layers_x;
					x[7] = coords[3 * node3 + 0] + n_l * X_pml / layers_x;

					pml_e->nodes[0].x = elem->nodes[aux[0]].x + 12 * (n_l + 0);
					pml_e->nodes[1].x = elem->nodes[aux[1]].x + 12 * (n_l + 1);
					pml_e->nodes[2].x = elem->nodes[aux[2]].x + 12 * (n_l + 1);
					pml_e->nodes[3].x = elem->nodes[aux[3]].x + 12 * (n_l + 0);
					pml_e->nodes[4].x = elem->nodes[aux[0]].x + 12 * (n_l + 0);
					pml_e->nodes[5].x = elem->nodes[aux[1]].x + 12 * (n_l + 1);
					pml_e->nodes[6].x = elem->nodes[aux[2]].x + 12 * (n_l + 1);
					pml_e->nodes[7].x = elem->nodes[aux[3]].x + 12 * (n_l + 0);

					y[0] = coords[3 * node0 + 1];
					y[1] = coords[3 * node0 + 1];
					y[2] = coords[3 * node1 + 1];
					y[3] = coords[3 * node1 + 1];
					y[4] = coords[3 * node3 + 1];
					y[5] = coords[3 * node3 + 1];
					y[6] = coords[3 * node2 + 1];
					y[7] = coords[3 * node2 + 1];

					pml_e->nodes[0].y = elem->nodes[aux[0]].y;
					pml_e->nodes[1].y = elem->nodes[aux[0]].y;
					pml_e->nodes[2].y = elem->nodes[aux[1]].y;
					pml_e->nodes[3].y = elem->nodes[aux[1]].y;
					pml_e->nodes[4].y = elem->nodes[aux[3]].y;
					pml_e->nodes[5].y = elem->nodes[aux[3]].y;
					pml_e->nodes[6].y = elem->nodes[aux[2]].y;
					pml_e->nodes[7].y = elem->nodes[aux[2]].y;

					z[0] = coords[3 * node0 + 2];
					z[1] = coords[3 * node0 + 2];
					z[2] = coords[3 * node1 + 2];
					z[3] = coords[3 * node1 + 2];
					z[4] = coords[3 * node3 + 2];
					z[5] = coords[3 * node3 + 2];
					z[6] = coords[3 * node2 + 2];
					z[7] = coords[3 * node2 + 2];

					pml_e->nodes[0].z = elem->nodes[aux[0]].z;
					pml_e->nodes[1].z = elem->nodes[aux[0]].z;
					pml_e->nodes[2].z = elem->nodes[aux[1]].z;
					pml_e->nodes[3].z = elem->nodes[aux[1]].z;
					pml_e->nodes[4].z = elem->nodes[aux[3]].z;
					pml_e->nodes[5].z = elem->nodes[aux[3]].z;
					pml_e->nodes[6].z = elem->nodes[aux[2]].z;
					pml_e->nodes[7].z = elem->nodes[aux[2]].z;

					for (int ino = 0; ino < 8; ino++)
					{
						// definindo ponto p a ser adicionado
						GtsPoint *p = gts_point_new(gts_point_class(), x[ino], y[ino], z[ino]);
						// adicionando ponto p
						int x = pml_e->nodes[ino].x;
						int y = pml_e->nodes[ino].y;
						int z = pml_e->nodes[ino].z;
						pml_e->nodes[ino].id = AddPoint(mesh, hash_nodes, p, coords, x, y, z);
					}
					isurf = 1;
					key.id = isurf + 1;
					bool lmat = sc_hash_array_lookup(hash_matpml, &key, &position);
					pmlmat_t *pmlT = (pmlmat_t *)sc_array_index(&hash_matpml->a, position);
					pml_e->n_mat = pmlT->mat;
					pmlT->xmin = std::min(pmlT->xmin, x[0]);
					pmlT->xmax = std::max(pmlT->xmax, x[1]);
					pmlT->ymin = 0;
					pmlT->ymax = 0;
					pmlT->zmin = 0;
					pmlT->zmax = 0;
				}
			}

			if (elem->surf[2].ext)
			{
				for (int n_l = 0; n_l < layers_y; ++n_l)
				{

					octant_t *pml_e = (octant_t *)sc_array_push(&mesh->elements);
					pml_e->id = mesh->elements.elem_count + 1;

					// nos de referencia
					int aux[4] = {0, 1, 4, 5};
					int node0 = elem->nodes[aux[0]].id;
					int node1 = elem->nodes[aux[1]].id;
					int node2 = elem->nodes[aux[2]].id;
					int node3 = elem->nodes[aux[3]].id;
					double x[8], y[8], z[8];

					x[0] = coords[3 * node0 + 0];
					x[1] = coords[3 * node0 + 0];
					x[2] = coords[3 * node1 + 0];
					x[3] = coords[3 * node1 + 0];
					x[4] = coords[3 * node2 + 0];
					x[5] = coords[3 * node2 + 0];
					x[6] = coords[3 * node3 + 0];
					x[7] = coords[3 * node3 + 0];

					pml_e->nodes[0].x = elem->nodes[aux[0]].x;
					pml_e->nodes[1].x = elem->nodes[aux[0]].x;
					pml_e->nodes[2].x = elem->nodes[aux[1]].x;
					pml_e->nodes[3].x = elem->nodes[aux[1]].x;
					pml_e->nodes[4].x = elem->nodes[aux[2]].x;
					pml_e->nodes[5].x = elem->nodes[aux[2]].x;
					pml_e->nodes[6].x = elem->nodes[aux[3]].x;
					pml_e->nodes[7].x = elem->nodes[aux[3]].x;

					y[0] = coords[3 * node0 + 1] - n_l * Y_pml / layers_y;
					y[1] = coords[3 * node1 + 1] - (n_l + 1) * Y_pml / layers_y;
					y[2] = coords[3 * node2 + 1] - (n_l + 1) * Y_pml / layers_y;
					y[3] = coords[3 * node3 + 1] - n_l * Y_pml / layers_y;
					y[4] = coords[3 * node0 + 1] - n_l * Y_pml / layers_y;
					y[5] = coords[3 * node1 + 1] - (n_l + 1) * Y_pml / layers_y;
					y[6] = coords[3 * node2 + 1] - (n_l + 1) * Y_pml / layers_y;
					y[7] = coords[3 * node3 + 1] - n_l * Y_pml / layers_y;

					pml_e->nodes[0].y = elem->nodes[aux[0]].y - 12 * (n_l + 0);
					pml_e->nodes[1].y = elem->nodes[aux[1]].y - 12 * (n_l + 1);
					pml_e->nodes[2].y = elem->nodes[aux[2]].y - 12 * (n_l + 1);
					pml_e->nodes[3].y = elem->nodes[aux[3]].y - 12 * (n_l + 0);
					pml_e->nodes[4].y = elem->nodes[aux[0]].y - 12 * (n_l + 0);
					pml_e->nodes[5].y = elem->nodes[aux[1]].y - 12 * (n_l + 1);
					pml_e->nodes[6].y = elem->nodes[aux[2]].y - 12 * (n_l + 1);
					pml_e->nodes[7].y = elem->nodes[aux[3]].y - 12 * (n_l + 0);

					z[0] = coords[3 * node0 + 2];
					z[1] = coords[3 * node0 + 2];
					z[2] = coords[3 * node1 + 2];
					z[3] = coords[3 * node1 + 2];
					z[4] = coords[3 * node2 + 2];
					z[5] = coords[3 * node2 + 2];
					z[6] = coords[3 * node3 + 2];
					z[7] = coords[3 * node3 + 2];

					pml_e->nodes[0].z = elem->nodes[aux[0]].z;
					pml_e->nodes[1].z = elem->nodes[aux[0]].z;
					pml_e->nodes[2].z = elem->nodes[aux[1]].z;
					pml_e->nodes[3].z = elem->nodes[aux[1]].z;
					pml_e->nodes[4].z = elem->nodes[aux[2]].z;
					pml_e->nodes[5].z = elem->nodes[aux[2]].z;
					pml_e->nodes[6].z = elem->nodes[aux[3]].z;
					pml_e->nodes[7].z = elem->nodes[aux[3]].z;

					for (int ino = 0; ino < 8; ino++)
					{
						// definindo ponto p a ser adicionado
						GtsPoint *p = gts_point_new(gts_point_class(), x[ino], y[ino], z[ino]);
						// adicionando ponto p
						int x = pml_e->nodes[ino].x;
						int y = pml_e->nodes[ino].y;
						int z = pml_e->nodes[ino].z;
						pml_e->nodes[ino].id = AddPoint(mesh, hash_nodes, p, coords, x, y, z);
					}
					isurf = 2;
					key.id = isurf + 1;
					bool lmat = sc_hash_array_lookup(hash_matpml, &key, &position);
					pmlmat_t *pmlT = (pmlmat_t *)sc_array_index(&hash_matpml->a, position);
					pml_e->n_mat = pmlT->mat;
					pmlT->xmin = 0;
					pmlT->xmax = 0;
					pmlT->ymin = std::min(pmlT->ymin, y[2]);
					pmlT->ymax = std::max(pmlT->ymax, y[0]);
					pmlT->zmin = 0;
					pmlT->zmax = 0;
				}
			}

			if (elem->surf[3].ext)
			{
				for (int n_l = 0; n_l < layers_y; ++n_l)
				{

					octant_t *pml_e = (octant_t *)sc_array_push(&mesh->elements);
					pml_e->id = mesh->elements.elem_count + 1;

					// nos de referencia
					int aux[4] = {3, 2, 7, 6};
					int node0 = elem->nodes[aux[0]].id;
					int node1 = elem->nodes[aux[1]].id;
					int node2 = elem->nodes[aux[2]].id;
					int node3 = elem->nodes[aux[3]].id;
					double x[8], y[8], z[8];

					x[0] = coords[3 * node0 + 0];
					x[1] = coords[3 * node0 + 0];
					x[2] = coords[3 * node1 + 0];
					x[3] = coords[3 * node1 + 0];
					x[4] = coords[3 * node2 + 0];
					x[5] = coords[3 * node2 + 0];
					x[6] = coords[3 * node3 + 0];
					x[7] = coords[3 * node3 + 0];

					pml_e->nodes[0].x = elem->nodes[aux[0]].x;
					pml_e->nodes[1].x = elem->nodes[aux[0]].x;
					pml_e->nodes[2].x = elem->nodes[aux[1]].x;
					pml_e->nodes[3].x = elem->nodes[aux[1]].x;
					pml_e->nodes[4].x = elem->nodes[aux[2]].x;
					pml_e->nodes[5].x = elem->nodes[aux[2]].x;
					pml_e->nodes[6].x = elem->nodes[aux[3]].x;
					pml_e->nodes[7].x = elem->nodes[aux[3]].x;

					y[0] = coords[3 * node0 + 1] + n_l * Y_pml / layers_y;
					y[1] = coords[3 * node1 + 1] + (n_l + 1) * Y_pml / layers_y;
					y[2] = coords[3 * node2 + 1] + (n_l + 1) * Y_pml / layers_y;
					y[3] = coords[3 * node3 + 1] + n_l * Y_pml / layers_y;
					y[4] = coords[3 * node0 + 1] + n_l * Y_pml / layers_y;
					y[5] = coords[3 * node1 + 1] + (n_l + 1) * Y_pml / layers_y;
					y[6] = coords[3 * node2 + 1] + (n_l + 1) * Y_pml / layers_y;
					y[7] = coords[3 * node3 + 1] + n_l * Y_pml / layers_y;

					pml_e->nodes[0].y = elem->nodes[aux[0]].y + 12 * (n_l + 0);
					pml_e->nodes[1].y = elem->nodes[aux[1]].y + 12 * (n_l + 1);
					pml_e->nodes[2].y = elem->nodes[aux[2]].y + 12 * (n_l + 1);
					pml_e->nodes[3].y = elem->nodes[aux[3]].y + 12 * (n_l + 0);
					pml_e->nodes[4].y = elem->nodes[aux[0]].y + 12 * (n_l + 0);
					pml_e->nodes[5].y = elem->nodes[aux[1]].y + 12 * (n_l + 1);
					pml_e->nodes[6].y = elem->nodes[aux[2]].y + 12 * (n_l + 1);
					pml_e->nodes[7].y = elem->nodes[aux[3]].y + 12 * (n_l + 0);

					z[0] = coords[3 * node0 + 2];
					z[1] = coords[3 * node0 + 2];
					z[2] = coords[3 * node1 + 2];
					z[3] = coords[3 * node1 + 2];
					z[4] = coords[3 * node2 + 2];
					z[5] = coords[3 * node2 + 2];
					z[6] = coords[3 * node3 + 2];
					z[7] = coords[3 * node3 + 2];

					pml_e->nodes[0].z = elem->nodes[aux[0]].z;
					pml_e->nodes[1].z = elem->nodes[aux[0]].z;
					pml_e->nodes[2].z = elem->nodes[aux[1]].z;
					pml_e->nodes[3].z = elem->nodes[aux[1]].z;
					pml_e->nodes[4].z = elem->nodes[aux[2]].z;
					pml_e->nodes[5].z = elem->nodes[aux[2]].z;
					pml_e->nodes[6].z = elem->nodes[aux[3]].z;
					pml_e->nodes[7].z = elem->nodes[aux[3]].z;

					for (int ino = 0; ino < 8; ino++)
					{
						// definindo ponto p a ser adicionado
						GtsPoint *p = gts_point_new(gts_point_class(), x[ino], y[ino], z[ino]);
						// adicionando ponto p
						int x = pml_e->nodes[ino].x;
						int y = pml_e->nodes[ino].y;
						int z = pml_e->nodes[ino].z;
						pml_e->nodes[ino].id = AddPoint(mesh, hash_nodes, p, coords, x, y, z);
					}
					isurf = 3;
					key.id = isurf + 1;
					bool lmat = sc_hash_array_lookup(hash_matpml, &key, &position);
					pmlmat_t *pmlT = (pmlmat_t *)sc_array_index(&hash_matpml->a, position);
					pml_e->n_mat = pmlT->mat;
					pmlT->xmin = 0;
					pmlT->xmax = 0;
					pmlT->ymin = std::min(pmlT->ymin, y[0]);
					pmlT->ymax = std::max(pmlT->ymax, y[2]);
					pmlT->zmin = 0;
					pmlT->zmax = 0;
				}
			}

			if (elem->surf[5].ext)
			{
				for (int n_l = 0; n_l < layers_z; ++n_l)
				{
					octant_t *pml_e = (octant_t *)sc_array_push(&mesh->elements);
					pml_e->id = mesh->elements.elem_count + 1;

					// nos de referencia
					int aux[4] = {4, 5, 6, 7};
					// int aux[4] = {0,1,2,3};
					int node0 = elem->nodes[aux[0]].id;
					int node1 = elem->nodes[aux[1]].id;
					int node2 = elem->nodes[aux[2]].id;
					int node3 = elem->nodes[aux[3]].id;
					double x[8], y[8], z[8];
					double zz[8];

					x[0] = coords[3 * node0 + 0];
					x[1] = coords[3 * node1 + 0];
					x[2] = coords[3 * node2 + 0];
					x[3] = coords[3 * node3 + 0];
					x[4] = coords[3 * node0 + 0];
					x[5] = coords[3 * node1 + 0];
					x[6] = coords[3 * node2 + 0];
					x[7] = coords[3 * node3 + 0];

					pml_e->nodes[0].x = elem->nodes[aux[0]].x;
					pml_e->nodes[1].x = elem->nodes[aux[1]].x;
					pml_e->nodes[2].x = elem->nodes[aux[2]].x;
					pml_e->nodes[3].x = elem->nodes[aux[3]].x;
					pml_e->nodes[4].x = elem->nodes[aux[0]].x;
					pml_e->nodes[5].x = elem->nodes[aux[1]].x;
					pml_e->nodes[6].x = elem->nodes[aux[2]].x;
					pml_e->nodes[7].x = elem->nodes[aux[3]].x;

					y[0] = coords[3 * node0 + 1];
					y[1] = coords[3 * node1 + 1];
					y[2] = coords[3 * node2 + 1];
					y[3] = coords[3 * node3 + 1];
					y[4] = coords[3 * node0 + 1];
					y[5] = coords[3 * node1 + 1];
					y[6] = coords[3 * node2 + 1];
					y[7] = coords[3 * node3 + 1];

					pml_e->nodes[0].y = elem->nodes[aux[0]].y;
					pml_e->nodes[1].y = elem->nodes[aux[1]].y;
					pml_e->nodes[2].y = elem->nodes[aux[2]].y;
					pml_e->nodes[3].y = elem->nodes[aux[3]].y;
					pml_e->nodes[4].y = elem->nodes[aux[0]].y;
					pml_e->nodes[5].y = elem->nodes[aux[1]].y;
					pml_e->nodes[6].y = elem->nodes[aux[2]].y;
					pml_e->nodes[7].y = elem->nodes[aux[3]].y;

					z[0] = coords[3 * node0 + 2] - (n_l + 1) * Z_pml / layers_z;
					z[1] = coords[3 * node1 + 2] - (n_l + 1) * Z_pml / layers_z;
					z[2] = coords[3 * node2 + 2] - (n_l + 1) * Z_pml / layers_z;
					z[3] = coords[3 * node3 + 2] - (n_l + 1) * Z_pml / layers_z;
					z[4] = coords[3 * node0 + 2] - (n_l)*Z_pml / layers_z;
					z[5] = coords[3 * node1 + 2] - (n_l)*Z_pml / layers_z;
					z[6] = coords[3 * node2 + 2] - (n_l)*Z_pml / layers_z;
					z[7] = coords[3 * node3 + 2] - (n_l)*Z_pml / layers_z;

					zz[0] = coords[3 * node0 + 2];
					zz[1] = coords[3 * node1 + 2];
					zz[2] = coords[3 * node2 + 2];
					zz[3] = coords[3 * node3 + 2];
					zz[4] = coords[3 * node0 + 2];
					zz[5] = coords[3 * node1 + 2];
					zz[6] = coords[3 * node2 + 2];
					zz[7] = coords[3 * node3 + 2];

					pml_e->nodes[0].z = elem->nodes[aux[0]].z + 12 * (n_l + 1);
					pml_e->nodes[1].z = elem->nodes[aux[1]].z + 12 * (n_l + 1);
					pml_e->nodes[2].z = elem->nodes[aux[2]].z + 12 * (n_l + 1);
					pml_e->nodes[3].z = elem->nodes[aux[3]].z + 12 * (n_l + 1);
					pml_e->nodes[4].z = elem->nodes[aux[0]].z + 12 * (n_l + 0);
					pml_e->nodes[5].z = elem->nodes[aux[1]].z + 12 * (n_l + 0);
					pml_e->nodes[6].z = elem->nodes[aux[2]].z + 12 * (n_l + 0);
					pml_e->nodes[7].z = elem->nodes[aux[3]].z + 12 * (n_l + 0);

					for (int ino = 0; ino < 8; ino++)
					{
						// definindo ponto p a ser adicionado
						GtsPoint *p = gts_point_new(gts_point_class(), x[ino], y[ino], z[ino]);
						// adicionando ponto p
						int x = pml_e->nodes[ino].x;
						int y = pml_e->nodes[ino].y;
						int z = pml_e->nodes[ino].z;
						pml_e->nodes[ino].id = AddPoint(mesh, hash_nodes, p, coords, x, y, z);
					}
					isurf = 5;
					key.id = isurf + 1;
					bool lmat = sc_hash_array_lookup(hash_matpml, &key, &position);
					pmlmat_t *pmlT = (pmlmat_t *)sc_array_index(&hash_matpml->a, position);
					pml_e->n_mat = pmlT->mat;
					pmlT->xmin = 0;
					pmlT->xmax = 0;
					pmlT->ymin = 0;
					pmlT->ymax = 0;
					pmlT->zmin = std::min(pmlT->zmin, z[0]);
					pmlT->zmax = std::max(pmlT->zmax, z[4]);
				}
			}

			if (elem->surf[4].ext && false)
			{
				for (int n_l = 0; n_l < layers_z; ++n_l)
				{
					// printf("Sou o el %d e entrei no 5\n",elemOrig->id);
					octant_t *pml_e = (octant_t *)sc_array_push(&mesh->elements);
					pml_e->id = mesh->elements.elem_count + 1;

					// nos de referencia
					int aux[4] = {0, 1, 2, 3};
					int node0 = elem->nodes[aux[0]].id;
					int node1 = elem->nodes[aux[1]].id;
					int node2 = elem->nodes[aux[2]].id;
					int node3 = elem->nodes[aux[3]].id;
					double x[8], y[8], z[8];

					x[0] = coords[3 * node0 + 0];
					x[1] = coords[3 * node1 + 0];
					x[2] = coords[3 * node2 + 0];
					x[3] = coords[3 * node3 + 0];
					x[4] = coords[3 * node0 + 0];
					x[5] = coords[3 * node1 + 0];
					x[6] = coords[3 * node2 + 0];
					x[7] = coords[3 * node3 + 0];

					pml_e->nodes[0].x = elem->nodes[aux[0]].x;
					pml_e->nodes[1].x = elem->nodes[aux[1]].x;
					pml_e->nodes[2].x = elem->nodes[aux[2]].x;
					pml_e->nodes[3].x = elem->nodes[aux[3]].x;
					pml_e->nodes[4].x = elem->nodes[aux[0]].x;
					pml_e->nodes[5].x = elem->nodes[aux[1]].x;
					pml_e->nodes[6].x = elem->nodes[aux[2]].x;
					pml_e->nodes[7].x = elem->nodes[aux[3]].x;

					y[0] = coords[3 * node0 + 1];
					y[1] = coords[3 * node1 + 1];
					y[2] = coords[3 * node2 + 1];
					y[3] = coords[3 * node3 + 1];
					y[4] = coords[3 * node0 + 1];
					y[5] = coords[3 * node1 + 1];
					y[6] = coords[3 * node2 + 1];
					y[7] = coords[3 * node3 + 1];

					pml_e->nodes[0].y = elem->nodes[aux[0]].y;
					pml_e->nodes[1].y = elem->nodes[aux[1]].y;
					pml_e->nodes[2].y = elem->nodes[aux[2]].y;
					pml_e->nodes[3].y = elem->nodes[aux[3]].y;
					pml_e->nodes[4].y = elem->nodes[aux[0]].y;
					pml_e->nodes[5].y = elem->nodes[aux[1]].y;
					pml_e->nodes[6].y = elem->nodes[aux[2]].y;
					pml_e->nodes[7].y = elem->nodes[aux[3]].y;

					z[0] = coords[3 * node0 + 2] + n_l * Z_pml / layers_z;
					z[1] = coords[3 * node1 + 2] + n_l * Z_pml / layers_z;
					z[2] = coords[3 * node2 + 2] + n_l * Z_pml / layers_z;
					z[3] = coords[3 * node3 + 2] + n_l * Z_pml / layers_z;
					z[4] = coords[3 * node0 + 2] + (n_l + 1) * Z_pml / layers_z;
					z[5] = coords[3 * node1 + 2] + (n_l + 1) * Z_pml / layers_z;
					z[6] = coords[3 * node2 + 2] + (n_l + 1) * Z_pml / layers_z;
					z[7] = coords[3 * node3 + 2] + (n_l + 1) * Z_pml / layers_z;

					pml_e->nodes[0].z = elem->nodes[aux[0]].z - 12 * (n_l + 1);
					pml_e->nodes[1].z = elem->nodes[aux[1]].z - 12 * (n_l + 1);
					pml_e->nodes[2].z = elem->nodes[aux[2]].z - 12 * (n_l + 1);
					pml_e->nodes[3].z = elem->nodes[aux[3]].z - 12 * (n_l + 1);
					pml_e->nodes[4].z = elem->nodes[aux[0]].z - 12 * (n_l + 0);
					pml_e->nodes[5].z = elem->nodes[aux[1]].z - 12 * (n_l + 0);
					pml_e->nodes[6].z = elem->nodes[aux[2]].z - 12 * (n_l + 0);
					pml_e->nodes[7].z = elem->nodes[aux[3]].z - 12 * (n_l + 0);

					for (int ino = 0; ino < 8; ino++)
					{
						// definindo ponto p a ser adicionado
						GtsPoint *p = gts_point_new(gts_point_class(), x[ino], y[ino], z[ino]);
						// adicionando ponto p
						int x = pml_e->nodes[ino].x;
						int y = pml_e->nodes[ino].y;
						int z = pml_e->nodes[ino].z;
						pml_e->nodes[ino].id = AddPoint(mesh, hash_nodes, p, coords, x, y, z);
					}
					isurf = 4;
					key.id = isurf + 1;
					bool lmat = sc_hash_array_lookup(hash_matpml, &key, &position);
					pmlmat_t *pmlT = (pmlmat_t *)sc_array_index(&hash_matpml->a, position);
					pml_e->n_mat = pmlT->mat;
					pmlT->xmin = 0;
					pmlT->xmax = 0;
					pmlT->ymin = 0;
					pmlT->ymax = 0;
					pmlT->zmin = std::min(pmlT->zmin, z[4]);
					pmlT->zmax = std::max(pmlT->zmax, z[0]);
				}
			}
		}

		if (edge)
		{
			key.tag = 1;
			if (elem->edge[0].ref && false)
			{

				for (int nz = 0; nz < layers_z; nz++)
				{
					for (int ny = 0; ny < layers_y; ny++)
					{

						octant_t *pml_e = (octant_t *)sc_array_push(&mesh->elements);
						pml_e->id = mesh->elements.elem_count + 1;

						// nos de referencia
						int aux[2] = {0, 1};
						int node0 = elem->nodes[aux[0]].id;
						int node1 = elem->nodes[aux[1]].id;
						double x[8], y[8], z[8];

						x[0] = coords[3 * node0 + 0];
						x[1] = coords[3 * node1 + 0];
						x[2] = coords[3 * node1 + 0];
						x[3] = coords[3 * node0 + 0];
						x[4] = coords[3 * node0 + 0];
						x[5] = coords[3 * node1 + 0];
						x[6] = coords[3 * node1 + 0];
						x[7] = coords[3 * node0 + 0];

						pml_e->nodes[0].x = elem->nodes[aux[0]].x;
						pml_e->nodes[1].x = elem->nodes[aux[1]].x;
						pml_e->nodes[2].x = elem->nodes[aux[1]].x;
						pml_e->nodes[3].x = elem->nodes[aux[0]].x;
						pml_e->nodes[4].x = elem->nodes[aux[0]].x;
						pml_e->nodes[5].x = elem->nodes[aux[1]].x;
						pml_e->nodes[6].x = elem->nodes[aux[1]].x;
						pml_e->nodes[7].x = elem->nodes[aux[0]].x;

						y[0] = coords[3 * node0 + 1] - ny * Y_pml / layers_y;
						y[1] = coords[3 * node1 + 1] - ny * Y_pml / layers_y;
						y[2] = coords[3 * node1 + 1] - (ny + 1) * Y_pml / layers_y;
						y[3] = coords[3 * node0 + 1] - (ny + 1) * Y_pml / layers_y;
						y[4] = coords[3 * node0 + 1] - ny * Y_pml / layers_y;
						y[5] = coords[3 * node1 + 1] - ny * Y_pml / layers_y;
						y[6] = coords[3 * node1 + 1] - (ny + 1) * Y_pml / layers_y;
						y[7] = coords[3 * node0 + 1] - (ny + 1) * Y_pml / layers_y;

						pml_e->nodes[0].y = elem->nodes[aux[0]].y - 12 * (ny + 0);
						pml_e->nodes[1].y = elem->nodes[aux[1]].y - 12 * (ny + 0);
						pml_e->nodes[2].y = elem->nodes[aux[1]].y - 12 * (ny + 1);
						pml_e->nodes[3].y = elem->nodes[aux[0]].y - 12 * (ny + 1);
						pml_e->nodes[4].y = elem->nodes[aux[0]].y - 12 * (ny + 0);
						pml_e->nodes[5].y = elem->nodes[aux[1]].y - 12 * (ny + 0);
						pml_e->nodes[6].y = elem->nodes[aux[1]].y - 12 * (ny + 1);
						pml_e->nodes[7].y = elem->nodes[aux[0]].y - 12 * (ny + 1);

						z[0] = coords[3 * node0 + 2] + nz * Z_pml / layers_z;
						z[1] = coords[3 * node1 + 2] + nz * Z_pml / layers_z;
						z[2] = coords[3 * node1 + 2] + nz * Z_pml / layers_z;
						z[3] = coords[3 * node0 + 2] + nz * Z_pml / layers_z;
						z[4] = coords[3 * node0 + 2] + (nz + 1) * Z_pml / layers_z;
						z[5] = coords[3 * node1 + 2] + (nz + 1) * Z_pml / layers_z;
						z[6] = coords[3 * node1 + 2] + (nz + 1) * Z_pml / layers_z;
						z[7] = coords[3 * node0 + 2] + (nz + 1) * Z_pml / layers_z;

						pml_e->nodes[0].z = elem->nodes[aux[0]].z - 12 * (nz + 1);
						pml_e->nodes[1].z = elem->nodes[aux[1]].z - 12 * (nz + 1);
						pml_e->nodes[2].z = elem->nodes[aux[1]].z - 12 * (nz + 1);
						pml_e->nodes[3].z = elem->nodes[aux[0]].z - 12 * (nz + 1);
						pml_e->nodes[4].z = elem->nodes[aux[0]].z - 12 * (nz + 0);
						pml_e->nodes[5].z = elem->nodes[aux[1]].z - 12 * (nz + 0);
						pml_e->nodes[6].z = elem->nodes[aux[1]].z - 12 * (nz + 0);
						pml_e->nodes[7].z = elem->nodes[aux[0]].z - 12 * (nz + 0);

						for (int ino = 0; ino < 8; ino++)
						{
							// definindo ponto p a ser adicionado
							GtsPoint *p = gts_point_new(gts_point_class(), x[ino], y[ino], z[ino]);
							// adicionando ponto p
							int x = pml_e->nodes[ino].x;
							int y = pml_e->nodes[ino].y;
							int z = pml_e->nodes[ino].z;
							pml_e->nodes[ino].id = AddPoint(mesh, hash_nodes, p, coords, x, y, z);
						}
						iedge = 0;
						key.id = 10 * (iedge + 1);
						bool lmat = sc_hash_array_lookup(hash_matpml, &key, &position);
						pmlmat_t *pmlT = (pmlmat_t *)sc_array_index(&hash_matpml->a, position);
						pml_e->n_mat = pmlT->mat;
						//it must be checked
						pmlT->xmin = 0;
						pmlT->xmax = 0;
						pmlT->ymin = std::min(pmlT->zmin, y[0]);
						pmlT->ymax = std::max(pmlT->zmin, y[2]);
						pmlT->zmin = std::min(pmlT->zmin, z[0]);
						pmlT->zmax = std::max(pmlT->zmax, z[4]);
					}
				}
			}

			if (elem->edge[1].ref && false)
			{

				for (int nz = 0; nz < layers_z; nz++)
				{
					for (int nx = 0; nx < layers_x; nx++)
					{

						octant_t *pml_e = (octant_t *)sc_array_push(&mesh->elements);
						pml_e->id = mesh->elements.elem_count + 1;

						// nos de referencia
						int aux[2] = {1, 2};
						int node0 = elem->nodes[aux[0]].id;
						int node1 = elem->nodes[aux[1]].id;
						double x[8], y[8], z[8];

						x[0] = coords[3 * node0 + 0] + nx * X_pml / layers_x;
						x[1] = coords[3 * node0 + 0] + (nx + 1) * X_pml / layers_x;
						x[2] = coords[3 * node1 + 0] + (nx + 1) * X_pml / layers_x;
						x[3] = coords[3 * node1 + 0] + nx * X_pml / layers_x;
						x[4] = coords[3 * node0 + 0] + nx * X_pml / layers_x;
						x[5] = coords[3 * node0 + 0] + (nx + 1) * X_pml / layers_x;
						x[6] = coords[3 * node1 + 0] + (nx + 1) * X_pml / layers_x;
						x[7] = coords[3 * node1 + 0] + nx * X_pml / layers_x;

						pml_e->nodes[0].x = elem->nodes[aux[0]].x + 12 * (nx + 0);
						pml_e->nodes[1].x = elem->nodes[aux[0]].x + 12 * (nx + 1);
						pml_e->nodes[2].x = elem->nodes[aux[1]].x + 12 * (nx + 1);
						pml_e->nodes[3].x = elem->nodes[aux[1]].x + 12 * (nx + 0);
						pml_e->nodes[4].x = elem->nodes[aux[0]].x + 12 * (nx + 0);
						pml_e->nodes[5].x = elem->nodes[aux[0]].x + 12 * (nx + 1);
						pml_e->nodes[6].x = elem->nodes[aux[1]].x + 12 * (nx + 1);
						pml_e->nodes[7].x = elem->nodes[aux[1]].x + 12 * (nx + 0);

						y[0] = coords[3 * node0 + 1];
						y[1] = coords[3 * node0 + 1];
						y[2] = coords[3 * node1 + 1];
						y[3] = coords[3 * node1 + 1];
						y[4] = coords[3 * node0 + 1];
						y[5] = coords[3 * node0 + 1];
						y[6] = coords[3 * node1 + 1];
						y[7] = coords[3 * node1 + 1];

						pml_e->nodes[0].y = elem->nodes[aux[0]].y;
						pml_e->nodes[1].y = elem->nodes[aux[0]].y;
						pml_e->nodes[2].y = elem->nodes[aux[1]].y;
						pml_e->nodes[3].y = elem->nodes[aux[1]].y;
						pml_e->nodes[4].y = elem->nodes[aux[0]].y;
						pml_e->nodes[5].y = elem->nodes[aux[0]].y;
						pml_e->nodes[6].y = elem->nodes[aux[1]].y;
						pml_e->nodes[7].y = elem->nodes[aux[1]].y;

						z[0] = coords[3 * node0 + 2] + nz * Z_pml / layers_z;
						z[1] = coords[3 * node1 + 2] + nz * Z_pml / layers_z;
						z[2] = coords[3 * node1 + 2] + nz * Z_pml / layers_z;
						z[3] = coords[3 * node0 + 2] + nz * Z_pml / layers_z;
						z[4] = coords[3 * node0 + 2] + (nz + 1) * Z_pml / layers_z;
						z[5] = coords[3 * node1 + 2] + (nz + 1) * Z_pml / layers_z;
						z[6] = coords[3 * node1 + 2] + (nz + 1) * Z_pml / layers_z;
						z[7] = coords[3 * node0 + 2] + (nz + 1) * Z_pml / layers_z;

						pml_e->nodes[0].z = elem->nodes[aux[0]].z - 12 * (nz + 1);
						pml_e->nodes[1].z = elem->nodes[aux[1]].z - 12 * (nz + 1);
						pml_e->nodes[2].z = elem->nodes[aux[1]].z - 12 * (nz + 1);
						pml_e->nodes[3].z = elem->nodes[aux[0]].z - 12 * (nz + 1);
						pml_e->nodes[4].z = elem->nodes[aux[0]].z - 12 * (nz + 0);
						pml_e->nodes[5].z = elem->nodes[aux[1]].z - 12 * (nz + 0);
						pml_e->nodes[6].z = elem->nodes[aux[1]].z - 12 * (nz + 0);
						pml_e->nodes[7].z = elem->nodes[aux[0]].z - 12 * (nz + 0);

						for (int ino = 0; ino < 8; ino++)
						{
							// definindo ponto p a ser adicionado
							GtsPoint *p = gts_point_new(gts_point_class(), x[ino], y[ino], z[ino]);
							// adicionando ponto p
							int x = pml_e->nodes[ino].x;
							int y = pml_e->nodes[ino].y;
							int z = pml_e->nodes[ino].z;
							pml_e->nodes[ino].id = AddPoint(mesh, hash_nodes, p, coords, x, y, z);
						}
						iedge = 1;
						key.id = 10 * (iedge + 1);
						bool lmat = sc_hash_array_lookup(hash_matpml, &key, &position);
						pmlmat_t *pmlT = (pmlmat_t *)sc_array_index(&hash_matpml->a, position);
						pml_e->n_mat = pmlT->mat;
						//it must be checked
						pmlT->xmin = std::min(pmlT->zmin, x[0]);
						pmlT->xmax = std::max(pmlT->zmin, x[1]);
						pmlT->ymin = 0;
						pmlT->ymax = 0;
						pmlT->zmin = std::min(pmlT->zmin, z[0]);
						pmlT->zmax = std::max(pmlT->zmax, z[4]);
					}
				}
			}

			if (elem->edge[2].ref && false)
			{

				for (int nz = 0; nz < layers_z; nz++)
				{
					for (int ny = 0; ny < layers_y; ny++)
					{

						octant_t *pml_e = (octant_t *)sc_array_push(&mesh->elements);
						pml_e->id = mesh->elements.elem_count + 1;

						// nos de referencia
						int aux[2] = {3, 2};
						int node0 = elem->nodes[aux[0]].id;
						int node1 = elem->nodes[aux[1]].id;
						double x[8], y[8], z[8];

						x[0] = coords[3 * node0 + 0];
						x[1] = coords[3 * node1 + 0];
						x[2] = coords[3 * node1 + 0];
						x[3] = coords[3 * node0 + 0];
						x[4] = coords[3 * node0 + 0];
						x[5] = coords[3 * node1 + 0];
						x[6] = coords[3 * node1 + 0];
						x[7] = coords[3 * node0 + 0];

						pml_e->nodes[0].x = elem->nodes[aux[0]].x;
						pml_e->nodes[1].x = elem->nodes[aux[1]].x;
						pml_e->nodes[2].x = elem->nodes[aux[1]].x;
						pml_e->nodes[3].x = elem->nodes[aux[0]].x;
						pml_e->nodes[4].x = elem->nodes[aux[0]].x;
						pml_e->nodes[5].x = elem->nodes[aux[1]].x;
						pml_e->nodes[6].x = elem->nodes[aux[1]].x;
						pml_e->nodes[7].x = elem->nodes[aux[0]].x;

						y[0] = coords[3 * node0 + 1] + ny * Y_pml / layers_y;
						y[1] = coords[3 * node1 + 1] + ny * Y_pml / layers_y;
						y[2] = coords[3 * node1 + 1] + (ny + 1) * Y_pml / layers_y;
						y[3] = coords[3 * node0 + 1] + (ny + 1) * Y_pml / layers_y;
						y[4] = coords[3 * node0 + 1] + ny * Y_pml / layers_y;
						y[5] = coords[3 * node1 + 1] + ny * Y_pml / layers_y;
						y[6] = coords[3 * node1 + 1] + (ny + 1) * Y_pml / layers_y;
						y[7] = coords[3 * node0 + 1] + (ny + 1) * Y_pml / layers_y;

						pml_e->nodes[0].y = elem->nodes[aux[0]].y + 12 * (ny + 0);
						pml_e->nodes[1].y = elem->nodes[aux[1]].y + 12 * (ny + 0);
						pml_e->nodes[2].y = elem->nodes[aux[1]].y + 12 * (ny + 1);
						pml_e->nodes[3].y = elem->nodes[aux[0]].y + 12 * (ny + 1);
						pml_e->nodes[4].y = elem->nodes[aux[0]].y + 12 * (ny + 0);
						pml_e->nodes[5].y = elem->nodes[aux[1]].y + 12 * (ny + 0);
						pml_e->nodes[6].y = elem->nodes[aux[1]].y + 12 * (ny + 1);
						pml_e->nodes[7].y = elem->nodes[aux[0]].y + 12 * (ny + 1);

						z[0] = coords[3 * node0 + 2] + nz * Z_pml / layers_z;
						z[1] = coords[3 * node1 + 2] + nz * Z_pml / layers_z;
						z[2] = coords[3 * node1 + 2] + nz * Z_pml / layers_z;
						z[3] = coords[3 * node0 + 2] + nz * Z_pml / layers_z;
						z[4] = coords[3 * node0 + 2] + (nz + 1) * Z_pml / layers_z;
						z[5] = coords[3 * node1 + 2] + (nz + 1) * Z_pml / layers_z;
						z[6] = coords[3 * node1 + 2] + (nz + 1) * Z_pml / layers_z;
						z[7] = coords[3 * node0 + 2] + (nz + 1) * Z_pml / layers_z;

						pml_e->nodes[0].z = elem->nodes[aux[0]].z - 12 * (nz + 1);
						pml_e->nodes[1].z = elem->nodes[aux[1]].z - 12 * (nz + 1);
						pml_e->nodes[2].z = elem->nodes[aux[1]].z - 12 * (nz + 1);
						pml_e->nodes[3].z = elem->nodes[aux[0]].z - 12 * (nz + 1);
						pml_e->nodes[4].z = elem->nodes[aux[0]].z - 12 * (nz + 0);
						pml_e->nodes[5].z = elem->nodes[aux[1]].z - 12 * (nz + 0);
						pml_e->nodes[6].z = elem->nodes[aux[1]].z - 12 * (nz + 0);
						pml_e->nodes[7].z = elem->nodes[aux[0]].z - 12 * (nz + 0);

						for (int ino = 0; ino < 8; ino++)
						{
							// definindo ponto p a ser adicionado
							GtsPoint *p = gts_point_new(gts_point_class(), x[ino], y[ino], z[ino]);
							// adicionando ponto p
							int x = pml_e->nodes[ino].x;
							int y = pml_e->nodes[ino].y;
							int z = pml_e->nodes[ino].z;
							pml_e->nodes[ino].id = AddPoint(mesh, hash_nodes, p, coords, x, y, z);
						}
						iedge = 2;
						key.id = 10 * (iedge + 1);
						bool lmat = sc_hash_array_lookup(hash_matpml, &key, &position);
						pmlmat_t *pmlT = (pmlmat_t *)sc_array_index(&hash_matpml->a, position);
						pml_e->n_mat = pmlT->mat;
						//it must be checked
						pmlT->xmin = 0;
						pmlT->xmax = 0;
						pmlT->ymin = std::min(pmlT->zmin, y[0]);
						pmlT->ymax = std::max(pmlT->zmin, y[2]);
						pmlT->zmin = std::min(pmlT->zmin, z[0]);
						pmlT->zmax = std::max(pmlT->zmax, z[4]);
					}
				}
			}

			if (elem->edge[3].ref && false)
			{

				for (int nz = 0; nz < layers_z; nz++)
				{
					for (int nx = 0; nx < layers_x; nx++)
					{

						octant_t *pml_e = (octant_t *)sc_array_push(&mesh->elements);
						pml_e->id = mesh->elements.elem_count + 1;

						// nos de referencia
						int aux[2] = {0, 3};
						int node0 = elem->nodes[aux[0]].id;
						int node1 = elem->nodes[aux[1]].id;
						double x[8], y[8], z[8];

						x[0] = coords[3 * node0 + 0] - nx * X_pml / layers_x;
						x[1] = coords[3 * node0 + 0] - (nx + 1) * X_pml / layers_x;
						x[2] = coords[3 * node1 + 0] - (nx + 1) * X_pml / layers_x;
						x[3] = coords[3 * node1 + 0] - nx * X_pml / layers_x;
						x[4] = coords[3 * node0 + 0] - nx * X_pml / layers_x;
						x[5] = coords[3 * node0 + 0] - (nx + 1) * X_pml / layers_x;
						x[6] = coords[3 * node1 + 0] - (nx + 1) * X_pml / layers_x;
						x[7] = coords[3 * node1 + 0] - nx * X_pml / layers_x;

						pml_e->nodes[0].x = elem->nodes[aux[0]].x - 12 * (nx + 0);
						pml_e->nodes[1].x = elem->nodes[aux[0]].x - 12 * (nx + 1);
						pml_e->nodes[2].x = elem->nodes[aux[1]].x - 12 * (nx + 1);
						pml_e->nodes[3].x = elem->nodes[aux[1]].x - 12 * (nx + 0);
						pml_e->nodes[4].x = elem->nodes[aux[0]].x - 12 * (nx + 0);
						pml_e->nodes[5].x = elem->nodes[aux[0]].x - 12 * (nx + 1);
						pml_e->nodes[6].x = elem->nodes[aux[1]].x - 12 * (nx + 1);
						pml_e->nodes[7].x = elem->nodes[aux[1]].x - 12 * (nx + 0);

						y[0] = coords[3 * node0 + 1];
						y[1] = coords[3 * node0 + 1];
						y[2] = coords[3 * node1 + 1];
						y[3] = coords[3 * node1 + 1];
						y[4] = coords[3 * node0 + 1];
						y[5] = coords[3 * node0 + 1];
						y[6] = coords[3 * node1 + 1];
						y[7] = coords[3 * node1 + 1];

						pml_e->nodes[0].y = elem->nodes[aux[0]].y;
						pml_e->nodes[1].y = elem->nodes[aux[0]].y;
						pml_e->nodes[2].y = elem->nodes[aux[1]].y;
						pml_e->nodes[3].y = elem->nodes[aux[1]].y;
						pml_e->nodes[4].y = elem->nodes[aux[0]].y;
						pml_e->nodes[5].y = elem->nodes[aux[0]].y;
						pml_e->nodes[6].y = elem->nodes[aux[1]].y;
						pml_e->nodes[7].y = elem->nodes[aux[1]].y;

						z[0] = coords[3 * node0 + 2] + nz * Z_pml / layers_z;
						z[1] = coords[3 * node1 + 2] + nz * Z_pml / layers_z;
						z[2] = coords[3 * node1 + 2] + nz * Z_pml / layers_z;
						z[3] = coords[3 * node0 + 2] + nz * Z_pml / layers_z;
						z[4] = coords[3 * node0 + 2] + (nz + 1) * Z_pml / layers_z;
						z[5] = coords[3 * node1 + 2] + (nz + 1) * Z_pml / layers_z;
						z[6] = coords[3 * node1 + 2] + (nz + 1) * Z_pml / layers_z;
						z[7] = coords[3 * node0 + 2] + (nz + 1) * Z_pml / layers_z;

						pml_e->nodes[0].z = elem->nodes[aux[0]].z - 12 * (nz + 1);
						pml_e->nodes[1].z = elem->nodes[aux[1]].z - 12 * (nz + 1);
						pml_e->nodes[2].z = elem->nodes[aux[1]].z - 12 * (nz + 1);
						pml_e->nodes[3].z = elem->nodes[aux[0]].z - 12 * (nz + 1);
						pml_e->nodes[4].z = elem->nodes[aux[0]].z - 12 * (nz + 0);
						pml_e->nodes[5].z = elem->nodes[aux[1]].z - 12 * (nz + 0);
						pml_e->nodes[6].z = elem->nodes[aux[1]].z - 12 * (nz + 0);
						pml_e->nodes[7].z = elem->nodes[aux[0]].z - 12 * (nz + 0);

						for (int ino = 0; ino < 8; ino++)
						{
							// definindo ponto p a ser adicionado
							GtsPoint *p = gts_point_new(gts_point_class(), x[ino], y[ino], z[ino]);
							// adicionando ponto p
							int x = pml_e->nodes[ino].x;
							int y = pml_e->nodes[ino].y;
							int z = pml_e->nodes[ino].z;
							pml_e->nodes[ino].id = AddPoint(mesh, hash_nodes, p, coords, x, y, z);
						}
						iedge = 3;
						key.id = 10 * (iedge + 1);
						bool lmat = sc_hash_array_lookup(hash_matpml, &key, &position);
						pmlmat_t *pmlT = (pmlmat_t *)sc_array_index(&hash_matpml->a, position);
						pml_e->n_mat = pmlT->mat;
						//it must be checked
						pmlT->xmin = std::min(pmlT->zmin, x[0]);
						pmlT->xmax = std::max(pmlT->zmin, x[1]);
						pmlT->ymin = 0;
						pmlT->ymax = 0;
						pmlT->zmin = std::min(pmlT->zmin, z[0]);
						pmlT->zmax = std::max(pmlT->zmax, z[4]);
					}
				}
			}

			if (elem->edge[4].ref)
			{
				for (int nx = 0; nx < layers_x; nx++)
				{
					for (int ny = 0; ny < layers_y; ny++)
					{

						octant_t *pml_e = (octant_t *)sc_array_push(&mesh->elements);
						pml_e->id = mesh->elements.elem_count + 1;

						// nos de referencia
						int aux[2] = {0, 4};
						int node0 = elem->nodes[aux[0]].id;
						int node1 = elem->nodes[aux[1]].id;
						double x[8], y[8], z[8];

						x[0] = coords[3 * node0 + 0] - nx * X_pml / layers_x;
						x[1] = coords[3 * node0 + 0] - (nx + 1) * X_pml / layers_x;
						x[2] = coords[3 * node0 + 0] - (nx + 1) * X_pml / layers_x;
						x[3] = coords[3 * node0 + 0] - nx * X_pml / layers_x;
						x[4] = coords[3 * node1 + 0] - nx * X_pml / layers_x;
						x[5] = coords[3 * node1 + 0] - (nx + 1) * X_pml / layers_x;
						x[6] = coords[3 * node1 + 0] - (nx + 1) * X_pml / layers_x;
						x[7] = coords[3 * node1 + 0] - nx * X_pml / layers_x;

						pml_e->nodes[0].x = elem->nodes[aux[0]].x - 12 * (nx + 0);
						pml_e->nodes[1].x = elem->nodes[aux[0]].x - 12 * (nx + 1);
						pml_e->nodes[2].x = elem->nodes[aux[0]].x - 12 * (nx + 1);
						pml_e->nodes[3].x = elem->nodes[aux[0]].x - 12 * (nx + 0);
						pml_e->nodes[4].x = elem->nodes[aux[1]].x - 12 * (nx + 0);
						pml_e->nodes[5].x = elem->nodes[aux[1]].x - 12 * (nx + 1);
						pml_e->nodes[6].x = elem->nodes[aux[1]].x - 12 * (nx + 1);
						pml_e->nodes[7].x = elem->nodes[aux[1]].x - 12 * (nx + 0);

						y[0] = coords[3 * node0 + 1] - ny * Y_pml / layers_y;
						y[1] = coords[3 * node0 + 1] - ny * Y_pml / layers_y;
						y[2] = coords[3 * node0 + 1] - (ny + 1) * Y_pml / layers_y;
						y[3] = coords[3 * node0 + 1] - (ny + 1) * Y_pml / layers_y;
						y[4] = coords[3 * node1 + 1] - ny * Y_pml / layers_y;
						y[5] = coords[3 * node1 + 1] - ny * Y_pml / layers_y;
						y[6] = coords[3 * node1 + 1] - (ny + 1) * Y_pml / layers_y;
						y[7] = coords[3 * node1 + 1] - (ny + 1) * Y_pml / layers_y;

						pml_e->nodes[0].y = elem->nodes[aux[0]].y - 12 * (ny + 0);
						pml_e->nodes[1].y = elem->nodes[aux[0]].y - 12 * (ny + 0);
						pml_e->nodes[2].y = elem->nodes[aux[0]].y - 12 * (ny + 1);
						pml_e->nodes[3].y = elem->nodes[aux[0]].y - 12 * (ny + 1);
						pml_e->nodes[4].y = elem->nodes[aux[1]].y - 12 * (ny + 0);
						pml_e->nodes[5].y = elem->nodes[aux[1]].y - 12 * (ny + 0);
						pml_e->nodes[6].y = elem->nodes[aux[1]].y - 12 * (ny + 1);
						pml_e->nodes[7].y = elem->nodes[aux[1]].y - 12 * (ny + 1);

						z[0] = coords[3 * node0 + 2];
						z[1] = coords[3 * node0 + 2];
						z[2] = coords[3 * node0 + 2];
						z[3] = coords[3 * node0 + 2];
						z[4] = coords[3 * node1 + 2];
						z[5] = coords[3 * node1 + 2];
						z[6] = coords[3 * node1 + 2];
						z[7] = coords[3 * node1 + 2];

						pml_e->nodes[0].z = elem->nodes[aux[0]].z;
						pml_e->nodes[1].z = elem->nodes[aux[0]].z;
						pml_e->nodes[2].z = elem->nodes[aux[0]].z;
						pml_e->nodes[3].z = elem->nodes[aux[0]].z;
						pml_e->nodes[4].z = elem->nodes[aux[1]].z;
						pml_e->nodes[5].z = elem->nodes[aux[1]].z;
						pml_e->nodes[6].z = elem->nodes[aux[1]].z;
						pml_e->nodes[7].z = elem->nodes[aux[1]].z;

						for (int ino = 0; ino < 8; ino++)
						{
							// definindo ponto p a ser adicionado
							GtsPoint *p = gts_point_new(gts_point_class(), x[ino], y[ino], z[ino]);
							// adicionando ponto p
							int x = pml_e->nodes[ino].x;
							int y = pml_e->nodes[ino].y;
							int z = pml_e->nodes[ino].z;
							pml_e->nodes[ino].id = AddPoint(mesh, hash_nodes, p, coords, x, y, z);
						}
						iedge = 4;
						key.id = 10 * (iedge + 1);
						bool lmat = sc_hash_array_lookup(hash_matpml, &key, &position);
						pmlmat_t *pmlT = (pmlmat_t *)sc_array_index(&hash_matpml->a, position);
						pml_e->n_mat = pmlT->mat;
						pmlT->xmin = std::min(pmlT->xmin, x[1]);
						pmlT->xmax = std::max(pmlT->xmax, x[0]);
						pmlT->ymin = std::min(pmlT->ymin, y[2]);
						pmlT->ymax = std::max(pmlT->ymax, y[0]);
						pmlT->zmin = 0;
						pmlT->zmax = 0;
					}
				}
			}

			if (elem->edge[5].ref)
			{
				for (int nx = 0; nx < layers_x; nx++)
				{
					for (int ny = 0; ny < layers_y; ny++)
					{

						octant_t *pml_e = (octant_t *)sc_array_push(&mesh->elements);
						pml_e->id = mesh->elements.elem_count + 1;

						// nos de referencia
						int aux[2] = {1, 5};
						int node0 = elem->nodes[aux[0]].id;
						int node1 = elem->nodes[aux[1]].id;
						double x[8], y[8], z[8];

						x[0] = coords[3 * node0 + 0] + nx * X_pml / layers_x;
						x[1] = coords[3 * node0 + 0] + (nx + 1) * X_pml / layers_x;
						x[2] = coords[3 * node0 + 0] + (nx + 1) * X_pml / layers_x;
						x[3] = coords[3 * node0 + 0] + nx * X_pml / layers_x;
						x[4] = coords[3 * node1 + 0] + nx * X_pml / layers_x;
						x[5] = coords[3 * node1 + 0] + (nx + 1) * X_pml / layers_x;
						x[6] = coords[3 * node1 + 0] + (nx + 1) * X_pml / layers_x;
						x[7] = coords[3 * node1 + 0] + nx * X_pml / layers_x;

						pml_e->nodes[0].x = elem->nodes[aux[0]].x + 12 * (nx + 0);
						pml_e->nodes[1].x = elem->nodes[aux[0]].x + 12 * (nx + 1);
						pml_e->nodes[2].x = elem->nodes[aux[0]].x + 12 * (nx + 1);
						pml_e->nodes[3].x = elem->nodes[aux[0]].x + 12 * (nx + 0);
						pml_e->nodes[4].x = elem->nodes[aux[1]].x + 12 * (nx + 0);
						pml_e->nodes[5].x = elem->nodes[aux[1]].x + 12 * (nx + 1);
						pml_e->nodes[6].x = elem->nodes[aux[1]].x + 12 * (nx + 1);
						pml_e->nodes[7].x = elem->nodes[aux[1]].x + 12 * (nx + 0);

						y[0] = coords[3 * node0 + 1] - ny * Y_pml / layers_y;
						y[1] = coords[3 * node0 + 1] - ny * Y_pml / layers_y;
						y[2] = coords[3 * node0 + 1] - (ny + 1) * Y_pml / layers_y;
						y[3] = coords[3 * node0 + 1] - (ny + 1) * Y_pml / layers_y;
						y[4] = coords[3 * node1 + 1] - ny * Y_pml / layers_y;
						y[5] = coords[3 * node1 + 1] - ny * Y_pml / layers_y;
						y[6] = coords[3 * node1 + 1] - (ny + 1) * Y_pml / layers_y;
						y[7] = coords[3 * node1 + 1] - (ny + 1) * Y_pml / layers_y;

						pml_e->nodes[0].y = elem->nodes[aux[0]].y - 12 * (ny + 0);
						pml_e->nodes[1].y = elem->nodes[aux[0]].y - 12 * (ny + 0);
						pml_e->nodes[2].y = elem->nodes[aux[0]].y - 12 * (ny + 1);
						pml_e->nodes[3].y = elem->nodes[aux[0]].y - 12 * (ny + 1);
						pml_e->nodes[4].y = elem->nodes[aux[1]].y - 12 * (ny + 0);
						pml_e->nodes[5].y = elem->nodes[aux[1]].y - 12 * (ny + 0);
						pml_e->nodes[6].y = elem->nodes[aux[1]].y - 12 * (ny + 1);
						pml_e->nodes[7].y = elem->nodes[aux[1]].y - 12 * (ny + 1);

						z[0] = coords[3 * node0 + 2];
						z[1] = coords[3 * node0 + 2];
						z[2] = coords[3 * node0 + 2];
						z[3] = coords[3 * node0 + 2];
						z[4] = coords[3 * node1 + 2];
						z[5] = coords[3 * node1 + 2];
						z[6] = coords[3 * node1 + 2];
						z[7] = coords[3 * node1 + 2];

						pml_e->nodes[0].z = elem->nodes[aux[0]].z;
						pml_e->nodes[1].z = elem->nodes[aux[0]].z;
						pml_e->nodes[2].z = elem->nodes[aux[0]].z;
						pml_e->nodes[3].z = elem->nodes[aux[0]].z;
						pml_e->nodes[4].z = elem->nodes[aux[1]].z;
						pml_e->nodes[5].z = elem->nodes[aux[1]].z;
						pml_e->nodes[6].z = elem->nodes[aux[1]].z;
						pml_e->nodes[7].z = elem->nodes[aux[1]].z;

						for (int ino = 0; ino < 8; ino++)
						{
							// definindo ponto p a ser adicionado
							GtsPoint *p = gts_point_new(gts_point_class(), x[ino], y[ino], z[ino]);
							// adicionando ponto p
							int x = pml_e->nodes[ino].x;
							int y = pml_e->nodes[ino].y;
							int z = pml_e->nodes[ino].z;
							pml_e->nodes[ino].id = AddPoint(mesh, hash_nodes, p, coords, x, y, z);
						}
						iedge = 5;
						key.id = 10 * (iedge + 1);
						bool lmat = sc_hash_array_lookup(hash_matpml, &key, &position);
						pmlmat_t *pmlT = (pmlmat_t *)sc_array_index(&hash_matpml->a, position);
						pml_e->n_mat = pmlT->mat;
						pmlT->xmin = std::min(pmlT->xmin, x[0]);
						pmlT->xmax = std::max(pmlT->xmax, x[1]);
						pmlT->ymin = std::min(pmlT->ymin, y[2]);
						pmlT->ymax = std::max(pmlT->ymax, y[0]);
						pmlT->zmin = 0;
						pmlT->zmax = 0;
					}
				}
			}

			if (elem->edge[6].ref)
			{
				for (int nx = 0; nx < layers_x; nx++)
				{
					for (int ny = 0; ny < layers_y; ny++)
					{

						octant_t *pml_e = (octant_t *)sc_array_push(&mesh->elements);
						pml_e->id = mesh->elements.elem_count + 1;

						// nos de referencia
						int aux[2] = {2, 6};
						int node0 = elem->nodes[aux[0]].id;
						int node1 = elem->nodes[aux[1]].id;
						double x[8], y[8], z[8];

						x[0] = coords[3 * node0 + 0] + nx * X_pml / layers_x;
						x[1] = coords[3 * node0 + 0] + (nx + 1) * X_pml / layers_x;
						x[2] = coords[3 * node0 + 0] + (nx + 1) * X_pml / layers_x;
						x[3] = coords[3 * node0 + 0] + nx * X_pml / layers_x;
						x[4] = coords[3 * node1 + 0] + nx * X_pml / layers_x;
						x[5] = coords[3 * node1 + 0] + (nx + 1) * X_pml / layers_x;
						x[6] = coords[3 * node1 + 0] + (nx + 1) * X_pml / layers_x;
						x[7] = coords[3 * node1 + 0] + nx * X_pml / layers_x;

						pml_e->nodes[0].x = elem->nodes[aux[0]].x + 12 * (nx + 0);
						pml_e->nodes[1].x = elem->nodes[aux[0]].x + 12 * (nx + 1);
						pml_e->nodes[2].x = elem->nodes[aux[0]].x + 12 * (nx + 1);
						pml_e->nodes[3].x = elem->nodes[aux[0]].x + 12 * (nx + 0);
						pml_e->nodes[4].x = elem->nodes[aux[1]].x + 12 * (nx + 0);
						pml_e->nodes[5].x = elem->nodes[aux[1]].x + 12 * (nx + 1);
						pml_e->nodes[6].x = elem->nodes[aux[1]].x + 12 * (nx + 1);
						pml_e->nodes[7].x = elem->nodes[aux[1]].x + 12 * (nx + 0);

						y[0] = coords[3 * node0 + 1] + ny * Y_pml / layers_y;
						y[1] = coords[3 * node0 + 1] + ny * Y_pml / layers_y;
						y[2] = coords[3 * node0 + 1] + (ny + 1) * Y_pml / layers_y;
						y[3] = coords[3 * node0 + 1] + (ny + 1) * Y_pml / layers_y;
						y[4] = coords[3 * node1 + 1] + ny * Y_pml / layers_y;
						y[5] = coords[3 * node1 + 1] + ny * Y_pml / layers_y;
						y[6] = coords[3 * node1 + 1] + (ny + 1) * Y_pml / layers_y;
						y[7] = coords[3 * node1 + 1] + (ny + 1) * Y_pml / layers_y;

						pml_e->nodes[0].y = elem->nodes[aux[0]].y + 12 * (ny + 0);
						pml_e->nodes[1].y = elem->nodes[aux[0]].y + 12 * (ny + 0);
						pml_e->nodes[2].y = elem->nodes[aux[0]].y + 12 * (ny + 1);
						pml_e->nodes[3].y = elem->nodes[aux[0]].y + 12 * (ny + 1);
						pml_e->nodes[4].y = elem->nodes[aux[1]].y + 12 * (ny + 0);
						pml_e->nodes[5].y = elem->nodes[aux[1]].y + 12 * (ny + 0);
						pml_e->nodes[6].y = elem->nodes[aux[1]].y + 12 * (ny + 1);
						pml_e->nodes[7].y = elem->nodes[aux[1]].y + 12 * (ny + 1);

						z[0] = coords[3 * node0 + 2];
						z[1] = coords[3 * node0 + 2];
						z[2] = coords[3 * node0 + 2];
						z[3] = coords[3 * node0 + 2];
						z[4] = coords[3 * node1 + 2];
						z[5] = coords[3 * node1 + 2];
						z[6] = coords[3 * node1 + 2];
						z[7] = coords[3 * node1 + 2];

						pml_e->nodes[0].z = elem->nodes[aux[0]].z;
						pml_e->nodes[1].z = elem->nodes[aux[0]].z;
						pml_e->nodes[2].z = elem->nodes[aux[0]].z;
						pml_e->nodes[3].z = elem->nodes[aux[0]].z;
						pml_e->nodes[4].z = elem->nodes[aux[1]].z;
						pml_e->nodes[5].z = elem->nodes[aux[1]].z;
						pml_e->nodes[6].z = elem->nodes[aux[1]].z;
						pml_e->nodes[7].z = elem->nodes[aux[1]].z;

						for (int ino = 0; ino < 8; ino++)
						{
							// definindo ponto p a ser adicionado
							GtsPoint *p = gts_point_new(gts_point_class(), x[ino], y[ino], z[ino]);
							// adicionando ponto p
							int x = pml_e->nodes[ino].x;
							int y = pml_e->nodes[ino].y;
							int z = pml_e->nodes[ino].z;
							pml_e->nodes[ino].id = AddPoint(mesh, hash_nodes, p, coords, x, y, z);
						}
						iedge = 6;
						key.id = 10 * (iedge + 1);
						bool lmat = sc_hash_array_lookup(hash_matpml, &key, &position);
						pmlmat_t *pmlT = (pmlmat_t *)sc_array_index(&hash_matpml->a, position);
						pml_e->n_mat = pmlT->mat;
						pmlT->xmin = std::min(pmlT->xmin, x[0]);
						pmlT->xmax = std::max(pmlT->xmax, x[1]);
						pmlT->ymin = std::min(pmlT->ymin, y[0]);
						pmlT->ymax = std::max(pmlT->ymax, y[2]);
						pmlT->zmin = 0;
						pmlT->zmax = 0;
					}
				}
			}

			if (elem->edge[7].ref)
			{
				for (int nx = 0; nx < layers_x; nx++)
				{
					for (int ny = 0; ny < layers_y; ny++)
					{

						octant_t *pml_e = (octant_t *)sc_array_push(&mesh->elements);
						pml_e->id = mesh->elements.elem_count + 1;

						// nos de referencia
						int aux[2] = {3, 7};
						int node0 = elem->nodes[aux[0]].id;
						int node1 = elem->nodes[aux[1]].id;
						double x[8], y[8], z[8];

						x[0] = coords[3 * node0 + 0] - nx * X_pml / layers_x;
						x[1] = coords[3 * node0 + 0] - (nx + 1) * X_pml / layers_x;
						x[2] = coords[3 * node0 + 0] - (nx + 1) * X_pml / layers_x;
						x[3] = coords[3 * node0 + 0] - nx * X_pml / layers_x;
						x[4] = coords[3 * node1 + 0] - nx * X_pml / layers_x;
						x[5] = coords[3 * node1 + 0] - (nx + 1) * X_pml / layers_x;
						x[6] = coords[3 * node1 + 0] - (nx + 1) * X_pml / layers_x;
						x[7] = coords[3 * node1 + 0] - nx * X_pml / layers_x;

						pml_e->nodes[0].x = elem->nodes[aux[0]].x - 12 * (nx + 0);
						pml_e->nodes[1].x = elem->nodes[aux[0]].x - 12 * (nx + 1);
						pml_e->nodes[2].x = elem->nodes[aux[0]].x - 12 * (nx + 1);
						pml_e->nodes[3].x = elem->nodes[aux[0]].x - 12 * (nx + 0);
						pml_e->nodes[4].x = elem->nodes[aux[1]].x - 12 * (nx + 0);
						pml_e->nodes[5].x = elem->nodes[aux[1]].x - 12 * (nx + 1);
						pml_e->nodes[6].x = elem->nodes[aux[1]].x - 12 * (nx + 1);
						pml_e->nodes[7].x = elem->nodes[aux[1]].x - 12 * (nx + 0);

						y[0] = coords[3 * node0 + 1] + ny * Y_pml / layers_y;
						y[1] = coords[3 * node0 + 1] + ny * Y_pml / layers_y;
						y[2] = coords[3 * node0 + 1] + (ny + 1) * Y_pml / layers_y;
						y[3] = coords[3 * node0 + 1] + (ny + 1) * Y_pml / layers_y;
						y[4] = coords[3 * node1 + 1] + ny * Y_pml / layers_y;
						y[5] = coords[3 * node1 + 1] + ny * Y_pml / layers_y;
						y[6] = coords[3 * node1 + 1] + (ny + 1) * Y_pml / layers_y;
						y[7] = coords[3 * node1 + 1] + (ny + 1) * Y_pml / layers_y;

						pml_e->nodes[0].y = elem->nodes[aux[0]].y + 12 * (ny + 0);
						pml_e->nodes[1].y = elem->nodes[aux[0]].y + 12 * (ny + 0);
						pml_e->nodes[2].y = elem->nodes[aux[0]].y + 12 * (ny + 1);
						pml_e->nodes[3].y = elem->nodes[aux[0]].y + 12 * (ny + 1);
						pml_e->nodes[4].y = elem->nodes[aux[1]].y + 12 * (ny + 0);
						pml_e->nodes[5].y = elem->nodes[aux[1]].y + 12 * (ny + 0);
						pml_e->nodes[6].y = elem->nodes[aux[1]].y + 12 * (ny + 1);
						pml_e->nodes[7].y = elem->nodes[aux[1]].y + 12 * (ny + 1);

						z[0] = coords[3 * node0 + 2];
						z[1] = coords[3 * node0 + 2];
						z[2] = coords[3 * node0 + 2];
						z[3] = coords[3 * node0 + 2];
						z[4] = coords[3 * node1 + 2];
						z[5] = coords[3 * node1 + 2];
						z[6] = coords[3 * node1 + 2];
						z[7] = coords[3 * node1 + 2];

						pml_e->nodes[0].z = elem->nodes[aux[0]].z;
						pml_e->nodes[1].z = elem->nodes[aux[0]].z;
						pml_e->nodes[2].z = elem->nodes[aux[0]].z;
						pml_e->nodes[3].z = elem->nodes[aux[0]].z;
						pml_e->nodes[4].z = elem->nodes[aux[1]].z;
						pml_e->nodes[5].z = elem->nodes[aux[1]].z;
						pml_e->nodes[6].z = elem->nodes[aux[1]].z;
						pml_e->nodes[7].z = elem->nodes[aux[1]].z;

						for (int ino = 0; ino < 8; ino++)
						{
							// definindo ponto p a ser adicionado
							GtsPoint *p = gts_point_new(gts_point_class(), x[ino], y[ino], z[ino]);
							// adicionando ponto p
							int x = pml_e->nodes[ino].x;
							int y = pml_e->nodes[ino].y;
							int z = pml_e->nodes[ino].z;
							pml_e->nodes[ino].id = AddPoint(mesh, hash_nodes, p, coords, x, y, z);
						}
						iedge = 7;
						key.id = 10 * (iedge + 1);
						bool lmat = sc_hash_array_lookup(hash_matpml, &key, &position);
						pmlmat_t *pmlT = (pmlmat_t *)sc_array_index(&hash_matpml->a, position);
						pml_e->n_mat = pmlT->mat;
						pmlT->xmin = std::min(pmlT->xmin, x[1]);
						pmlT->xmax = std::max(pmlT->xmax, x[0]);
						pmlT->ymin = std::min(pmlT->ymin, y[0]);
						pmlT->ymax = std::max(pmlT->ymax, y[2]);
						pmlT->zmin = 0;
						pmlT->zmax = 0;
					}
				}
			}

			if (elem->edge[8].ref)
			{
				for (int nz = 0; nz < layers_z; nz++)
				{
					for (int ny = 0; ny < layers_y; ny++)
					{

						octant_t *pml_e = (octant_t *)sc_array_push(&mesh->elements);
						pml_e->id = mesh->elements.elem_count + 1;

						// nos de referencia
						int aux[2] = {4, 5};
						int node0 = elem->nodes[aux[0]].id;
						int node1 = elem->nodes[aux[1]].id;
						double x[8], y[8], z[8];

						x[0] = coords[3 * node0 + 0];
						x[1] = coords[3 * node1 + 0];
						x[2] = coords[3 * node1 + 0];
						x[3] = coords[3 * node0 + 0];
						x[4] = coords[3 * node0 + 0];
						x[5] = coords[3 * node1 + 0];
						x[6] = coords[3 * node1 + 0];
						x[7] = coords[3 * node0 + 0];

						pml_e->nodes[0].x = elem->nodes[aux[0]].x;
						pml_e->nodes[1].x = elem->nodes[aux[1]].x;
						pml_e->nodes[2].x = elem->nodes[aux[1]].x;
						pml_e->nodes[3].x = elem->nodes[aux[0]].x;
						pml_e->nodes[4].x = elem->nodes[aux[0]].x;
						pml_e->nodes[5].x = elem->nodes[aux[1]].x;
						pml_e->nodes[6].x = elem->nodes[aux[1]].x;
						pml_e->nodes[7].x = elem->nodes[aux[0]].x;

						y[0] = coords[3 * node0 + 1] - ny * Y_pml / layers_y;
						y[1] = coords[3 * node1 + 1] - ny * Y_pml / layers_y;
						y[2] = coords[3 * node1 + 1] - (ny + 1) * Y_pml / layers_y;
						y[3] = coords[3 * node0 + 1] - (ny + 1) * Y_pml / layers_y;
						y[4] = coords[3 * node0 + 1] - ny * Y_pml / layers_y;
						y[5] = coords[3 * node1 + 1] - ny * Y_pml / layers_y;
						y[6] = coords[3 * node1 + 1] - (ny + 1) * Y_pml / layers_y;
						y[7] = coords[3 * node0 + 1] - (ny + 1) * Y_pml / layers_y;

						pml_e->nodes[0].y = elem->nodes[aux[0]].y - 12 * (ny + 0);
						pml_e->nodes[1].y = elem->nodes[aux[1]].y - 12 * (ny + 0);
						pml_e->nodes[2].y = elem->nodes[aux[1]].y - 12 * (ny + 1);
						pml_e->nodes[3].y = elem->nodes[aux[0]].y - 12 * (ny + 1);
						pml_e->nodes[4].y = elem->nodes[aux[0]].y - 12 * (ny + 0);
						pml_e->nodes[5].y = elem->nodes[aux[1]].y - 12 * (ny + 0);
						pml_e->nodes[6].y = elem->nodes[aux[1]].y - 12 * (ny + 1);
						pml_e->nodes[7].y = elem->nodes[aux[0]].y - 12 * (ny + 1);

						z[0] = coords[3 * node0 + 2] - (nz + 1) * Z_pml / layers_z;
						z[1] = coords[3 * node1 + 2] - (nz + 1) * Z_pml / layers_z;
						z[2] = coords[3 * node1 + 2] - (nz + 1) * Z_pml / layers_z;
						z[3] = coords[3 * node0 + 2] - (nz + 1) * Z_pml / layers_z;
						z[4] = coords[3 * node0 + 2] - (nz)*Z_pml / layers_z;
						z[5] = coords[3 * node1 + 2] - (nz)*Z_pml / layers_z;
						z[6] = coords[3 * node1 + 2] - (nz)*Z_pml / layers_z;
						z[7] = coords[3 * node0 + 2] - (nz)*Z_pml / layers_z;

						pml_e->nodes[0].z = elem->nodes[aux[0]].z + 12 * (nz + 1);
						pml_e->nodes[1].z = elem->nodes[aux[1]].z + 12 * (nz + 1);
						pml_e->nodes[2].z = elem->nodes[aux[1]].z + 12 * (nz + 1);
						pml_e->nodes[3].z = elem->nodes[aux[0]].z + 12 * (nz + 1);
						pml_e->nodes[4].z = elem->nodes[aux[0]].z + 12 * (nz + 0);
						pml_e->nodes[5].z = elem->nodes[aux[1]].z + 12 * (nz + 0);
						pml_e->nodes[6].z = elem->nodes[aux[1]].z + 12 * (nz + 0);
						pml_e->nodes[7].z = elem->nodes[aux[0]].z + 12 * (nz + 0);

						for (int ino = 0; ino < 8; ino++)
						{
							// definindo ponto p a ser adicionado
							GtsPoint *p = gts_point_new(gts_point_class(), x[ino], y[ino], z[ino]);
							// adicionando ponto p
							int x = pml_e->nodes[ino].x;
							int y = pml_e->nodes[ino].y;
							int z = pml_e->nodes[ino].z;
							pml_e->nodes[ino].id = AddPoint(mesh, hash_nodes, p, coords, x, y, z);
						}
						iedge = 8;
						key.id = 10 * (iedge + 1);
						bool lmat = sc_hash_array_lookup(hash_matpml, &key, &position);
						pmlmat_t *pmlT = (pmlmat_t *)sc_array_index(&hash_matpml->a, position);
						pml_e->n_mat = pmlT->mat;
						pmlT->xmin = 0;
						pmlT->xmax = 0;
						pmlT->ymin = std::min(pmlT->ymin, y[2]);
						pmlT->ymax = std::max(pmlT->ymax, y[0]);
						pmlT->zmin = std::min(pmlT->zmin, z[0]);
						pmlT->zmax = std::max(pmlT->zmax, z[4]);
					}
				}
			}

			if (elem->edge[9].ref)
			{

				for (int nz = 0; nz < layers_z; nz++)
				{
					for (int nx = 0; nx < layers_x; nx++)
					{

						octant_t *pml_e = (octant_t *)sc_array_push(&mesh->elements);
						pml_e->id = mesh->elements.elem_count + 1;

						// nos de referencia
						int aux[2] = {5, 6};
						int node0 = elem->nodes[aux[0]].id;
						int node1 = elem->nodes[aux[1]].id;
						double x[8], y[8], z[8];

						x[0] = coords[3 * node0 + 0] + nx * X_pml / layers_x;
						x[1] = coords[3 * node0 + 0] + (nx + 1) * X_pml / layers_x;
						x[2] = coords[3 * node1 + 0] + (nx + 1) * X_pml / layers_x;
						x[3] = coords[3 * node1 + 0] + nx * X_pml / layers_x;
						x[4] = coords[3 * node0 + 0] + nx * X_pml / layers_x;
						x[5] = coords[3 * node0 + 0] + (nx + 1) * X_pml / layers_x;
						x[6] = coords[3 * node1 + 0] + (nx + 1) * X_pml / layers_x;
						x[7] = coords[3 * node1 + 0] + nx * X_pml / layers_x;

						pml_e->nodes[0].x = elem->nodes[aux[0]].x + 12 * (nx + 0);
						pml_e->nodes[1].x = elem->nodes[aux[0]].x + 12 * (nx + 1);
						pml_e->nodes[2].x = elem->nodes[aux[1]].x + 12 * (nx + 1);
						pml_e->nodes[3].x = elem->nodes[aux[1]].x + 12 * (nx + 0);
						pml_e->nodes[4].x = elem->nodes[aux[0]].x + 12 * (nx + 0);
						pml_e->nodes[5].x = elem->nodes[aux[0]].x + 12 * (nx + 1);
						pml_e->nodes[6].x = elem->nodes[aux[1]].x + 12 * (nx + 1);
						pml_e->nodes[7].x = elem->nodes[aux[1]].x + 12 * (nx + 0);

						y[0] = coords[3 * node0 + 1];
						y[1] = coords[3 * node0 + 1];
						y[2] = coords[3 * node1 + 1];
						y[3] = coords[3 * node1 + 1];
						y[4] = coords[3 * node0 + 1];
						y[5] = coords[3 * node0 + 1];
						y[6] = coords[3 * node1 + 1];
						y[7] = coords[3 * node1 + 1];

						pml_e->nodes[0].y = elem->nodes[aux[0]].y;
						pml_e->nodes[1].y = elem->nodes[aux[0]].y;
						pml_e->nodes[2].y = elem->nodes[aux[1]].y;
						pml_e->nodes[3].y = elem->nodes[aux[1]].y;
						pml_e->nodes[4].y = elem->nodes[aux[0]].y;
						pml_e->nodes[5].y = elem->nodes[aux[0]].y;
						pml_e->nodes[6].y = elem->nodes[aux[1]].y;
						pml_e->nodes[7].y = elem->nodes[aux[1]].y;

						z[0] = coords[3 * node0 + 2] - (nz + 1) * Z_pml / layers_z;
						z[1] = coords[3 * node0 + 2] - (nz + 1) * Z_pml / layers_z;
						z[2] = coords[3 * node1 + 2] - (nz + 1) * Z_pml / layers_z;
						z[3] = coords[3 * node1 + 2] - (nz + 1) * Z_pml / layers_z;
						z[4] = coords[3 * node0 + 2] - (nz)*Z_pml / layers_z;
						z[5] = coords[3 * node0 + 2] - (nz)*Z_pml / layers_z;
						z[6] = coords[3 * node1 + 2] - (nz)*Z_pml / layers_z;
						z[7] = coords[3 * node0 + 2] - (nz)*Z_pml / layers_z;

						pml_e->nodes[0].z = elem->nodes[aux[0]].z + 12 * (nz + 1);
						pml_e->nodes[1].z = elem->nodes[aux[0]].z + 12 * (nz + 1);
						pml_e->nodes[2].z = elem->nodes[aux[1]].z + 12 * (nz + 1);
						pml_e->nodes[3].z = elem->nodes[aux[1]].z + 12 * (nz + 1);
						pml_e->nodes[4].z = elem->nodes[aux[0]].z + 12 * (nz + 0);
						pml_e->nodes[5].z = elem->nodes[aux[0]].z + 12 * (nz + 0);
						pml_e->nodes[6].z = elem->nodes[aux[1]].z + 12 * (nz + 0);
						pml_e->nodes[7].z = elem->nodes[aux[1]].z + 12 * (nz + 0);

						for (int ino = 0; ino < 8; ino++)
						{
							// definindo ponto p a ser adicionado
							GtsPoint *p = gts_point_new(gts_point_class(), x[ino], y[ino], z[ino]);
							// adicionando ponto p
							int x = pml_e->nodes[ino].x;
							int y = pml_e->nodes[ino].y;
							int z = pml_e->nodes[ino].z;
							pml_e->nodes[ino].id = AddPoint(mesh, hash_nodes, p, coords, x, y, z);
						}
						iedge = 9;
						key.id = 10 * (iedge + 1);
						bool lmat = sc_hash_array_lookup(hash_matpml, &key, &position);
						pmlmat_t *pmlT = (pmlmat_t *)sc_array_index(&hash_matpml->a, position);
						pml_e->n_mat = pmlT->mat;
						pmlT->xmin = std::min(pmlT->xmin, x[1]);
						pmlT->xmax = std::max(pmlT->xmax, x[0]);
						pmlT->ymin = 0;
						pmlT->ymax = 0;
						pmlT->zmin = std::min(pmlT->zmin, z[0]);
						pmlT->zmax = std::max(pmlT->zmax, z[4]);
					}
				}
			}

			if (elem->edge[10].ref)
			{
				for (int nz = 0; nz < layers_z; nz++)
				{
					for (int ny = 0; ny < layers_y; ny++)
					{

						octant_t *pml_e = (octant_t *)sc_array_push(&mesh->elements);
						pml_e->id = mesh->elements.elem_count + 1;

						// nos de referencia
						int aux[2] = {7, 6};
						int node0 = elem->nodes[aux[0]].id;
						int node1 = elem->nodes[aux[1]].id;
						double x[8], y[8], z[8];

						x[0] = coords[3 * node0 + 0];
						x[1] = coords[3 * node1 + 0];
						x[2] = coords[3 * node1 + 0];
						x[3] = coords[3 * node0 + 0];
						x[4] = coords[3 * node0 + 0];
						x[5] = coords[3 * node1 + 0];
						x[6] = coords[3 * node1 + 0];
						x[7] = coords[3 * node0 + 0];

						pml_e->nodes[0].x = elem->nodes[aux[0]].x;
						pml_e->nodes[1].x = elem->nodes[aux[1]].x;
						pml_e->nodes[2].x = elem->nodes[aux[1]].x;
						pml_e->nodes[3].x = elem->nodes[aux[0]].x;
						pml_e->nodes[4].x = elem->nodes[aux[0]].x;
						pml_e->nodes[5].x = elem->nodes[aux[1]].x;
						pml_e->nodes[6].x = elem->nodes[aux[1]].x;
						pml_e->nodes[7].x = elem->nodes[aux[0]].x;

						y[0] = coords[3 * node0 + 1] + ny * Y_pml / layers_y;
						y[1] = coords[3 * node1 + 1] + ny * Y_pml / layers_y;
						y[2] = coords[3 * node1 + 1] + (ny + 1) * Y_pml / layers_y;
						y[3] = coords[3 * node0 + 1] + (ny + 1) * Y_pml / layers_y;
						y[4] = coords[3 * node0 + 1] + ny * Y_pml / layers_y;
						y[5] = coords[3 * node1 + 1] + ny * Y_pml / layers_y;
						y[6] = coords[3 * node1 + 1] + (ny + 1) * Y_pml / layers_y;
						y[7] = coords[3 * node0 + 1] + (ny + 1) * Y_pml / layers_y;

						pml_e->nodes[0].y = elem->nodes[aux[0]].y + 12 * (ny + 0);
						pml_e->nodes[1].y = elem->nodes[aux[1]].y + 12 * (ny + 0);
						pml_e->nodes[2].y = elem->nodes[aux[1]].y + 12 * (ny + 1);
						pml_e->nodes[3].y = elem->nodes[aux[0]].y + 12 * (ny + 1);
						pml_e->nodes[4].y = elem->nodes[aux[0]].y + 12 * (ny + 0);
						pml_e->nodes[5].y = elem->nodes[aux[1]].y + 12 * (ny + 0);
						pml_e->nodes[6].y = elem->nodes[aux[1]].y + 12 * (ny + 1);
						pml_e->nodes[7].y = elem->nodes[aux[0]].y + 12 * (ny + 1);

						z[0] = coords[3 * node0 + 2] - (nz + 1) * Z_pml / layers_z;
						z[1] = coords[3 * node1 + 2] - (nz + 1) * Z_pml / layers_z;
						z[2] = coords[3 * node1 + 2] - (nz + 1) * Z_pml / layers_z;
						z[3] = coords[3 * node0 + 2] - (nz + 1) * Z_pml / layers_z;
						z[4] = coords[3 * node0 + 2] - (nz)*Z_pml / layers_z;
						z[5] = coords[3 * node1 + 2] - (nz)*Z_pml / layers_z;
						z[6] = coords[3 * node1 + 2] - (nz)*Z_pml / layers_z;
						z[7] = coords[3 * node0 + 2] - (nz)*Z_pml / layers_z;

						pml_e->nodes[0].z = elem->nodes[aux[0]].z + 12 * (nz + 1);
						pml_e->nodes[1].z = elem->nodes[aux[1]].z + 12 * (nz + 1);
						pml_e->nodes[2].z = elem->nodes[aux[1]].z + 12 * (nz + 1);
						pml_e->nodes[3].z = elem->nodes[aux[0]].z + 12 * (nz + 1);
						pml_e->nodes[4].z = elem->nodes[aux[0]].z + 12 * (nz + 0);
						pml_e->nodes[5].z = elem->nodes[aux[1]].z + 12 * (nz + 0);
						pml_e->nodes[6].z = elem->nodes[aux[1]].z + 12 * (nz + 0);
						pml_e->nodes[7].z = elem->nodes[aux[0]].z + 12 * (nz + 0);

						for (int ino = 0; ino < 8; ino++)
						{
							// definindo ponto p a ser adicionado
							GtsPoint *p = gts_point_new(gts_point_class(), x[ino], y[ino], z[ino]);
							// adicionando ponto p
							int x = pml_e->nodes[ino].x;
							int y = pml_e->nodes[ino].y;
							int z = pml_e->nodes[ino].z;
							pml_e->nodes[ino].id = AddPoint(mesh, hash_nodes, p, coords, x, y, z);
						}
						iedge = 10;
						key.id = 10 * (iedge + 1);
						bool lmat = sc_hash_array_lookup(hash_matpml, &key, &position);
						pmlmat_t *pmlT = (pmlmat_t *)sc_array_index(&hash_matpml->a, position);
						pml_e->n_mat = pmlT->mat;
						pmlT->xmin = 0;
						pmlT->xmax = 0;
						pmlT->ymin = std::min(pmlT->ymin, y[0]);
						pmlT->ymax = std::max(pmlT->ymax, y[2]);
						pmlT->zmin = std::min(pmlT->zmin, z[0]);
						pmlT->zmax = std::max(pmlT->zmax, z[4]);
					}
				}
			}

			if (elem->edge[11].ref)
			{
				for (int nz = 0; nz < layers_z; nz++)
				{
					for (int nx = 0; nx < layers_x; nx++)
					{

						octant_t *pml_e = (octant_t *)sc_array_push(&mesh->elements);
						pml_e->id = mesh->elements.elem_count + 1;

						// nos de referencia
						int aux[2] = {4, 7};
						int node0 = elem->nodes[aux[0]].id;
						int node1 = elem->nodes[aux[1]].id;
						double x[8], y[8], z[8];

						x[0] = coords[3 * node0 + 0] - nx * X_pml / layers_x;
						x[1] = coords[3 * node0 + 0] - (nx + 1) * X_pml / layers_x;
						x[2] = coords[3 * node1 + 0] - (nx + 1) * X_pml / layers_x;
						x[3] = coords[3 * node1 + 0] - nx * X_pml / layers_x;
						x[4] = coords[3 * node0 + 0] - nx * X_pml / layers_x;
						x[5] = coords[3 * node0 + 0] - (nx + 1) * X_pml / layers_x;
						x[6] = coords[3 * node1 + 0] - (nx + 1) * X_pml / layers_x;
						x[7] = coords[3 * node1 + 0] - nx * X_pml / layers_x;

						pml_e->nodes[0].x = elem->nodes[aux[0]].x - 12 * (nx + 0);
						pml_e->nodes[1].x = elem->nodes[aux[0]].x - 12 * (nx + 1);
						pml_e->nodes[2].x = elem->nodes[aux[1]].x - 12 * (nx + 1);
						pml_e->nodes[3].x = elem->nodes[aux[1]].x - 12 * (nx + 0);
						pml_e->nodes[4].x = elem->nodes[aux[0]].x - 12 * (nx + 0);
						pml_e->nodes[5].x = elem->nodes[aux[0]].x - 12 * (nx + 1);
						pml_e->nodes[6].x = elem->nodes[aux[1]].x - 12 * (nx + 1);
						pml_e->nodes[7].x = elem->nodes[aux[1]].x - 12 * (nx + 0);

						y[0] = coords[3 * node0 + 1];
						y[1] = coords[3 * node0 + 1];
						y[2] = coords[3 * node1 + 1];
						y[3] = coords[3 * node1 + 1];
						y[4] = coords[3 * node0 + 1];
						y[5] = coords[3 * node0 + 1];
						y[6] = coords[3 * node1 + 1];
						y[7] = coords[3 * node1 + 1];

						pml_e->nodes[0].y = elem->nodes[aux[0]].y;
						pml_e->nodes[1].y = elem->nodes[aux[0]].y;
						pml_e->nodes[2].y = elem->nodes[aux[1]].y;
						pml_e->nodes[3].y = elem->nodes[aux[1]].y;
						pml_e->nodes[4].y = elem->nodes[aux[0]].y;
						pml_e->nodes[5].y = elem->nodes[aux[0]].y;
						pml_e->nodes[6].y = elem->nodes[aux[1]].y;
						pml_e->nodes[7].y = elem->nodes[aux[1]].y;

						z[0] = coords[3 * node0 + 2] - (nz + 1) * Z_pml / layers_z;
						z[1] = coords[3 * node0 + 2] - (nz + 1) * Z_pml / layers_z;
						z[2] = coords[3 * node1 + 2] - (nz + 1) * Z_pml / layers_z;
						z[3] = coords[3 * node1 + 2] - (nz + 1) * Z_pml / layers_z;
						z[4] = coords[3 * node0 + 2] - (nz)*Z_pml / layers_z;
						z[5] = coords[3 * node0 + 2] - (nz)*Z_pml / layers_z;
						z[6] = coords[3 * node1 + 2] - (nz)*Z_pml / layers_z;
						z[7] = coords[3 * node0 + 2] - (nz)*Z_pml / layers_z;

						pml_e->nodes[0].z = elem->nodes[aux[0]].z + 12 * (nz + 1);
						pml_e->nodes[1].z = elem->nodes[aux[0]].z + 12 * (nz + 1);
						pml_e->nodes[2].z = elem->nodes[aux[1]].z + 12 * (nz + 1);
						pml_e->nodes[3].z = elem->nodes[aux[1]].z + 12 * (nz + 1);
						pml_e->nodes[4].z = elem->nodes[aux[0]].z + 12 * (nz + 0);
						pml_e->nodes[5].z = elem->nodes[aux[0]].z + 12 * (nz + 0);
						pml_e->nodes[6].z = elem->nodes[aux[1]].z + 12 * (nz + 0);
						pml_e->nodes[7].z = elem->nodes[aux[1]].z + 12 * (nz + 0);

						for (int ino = 0; ino < 8; ino++)
						{
							// definindo ponto p a ser adicionado
							GtsPoint *p = gts_point_new(gts_point_class(), x[ino], y[ino], z[ino]);
							// adicionando ponto p
							int x = pml_e->nodes[ino].x;
							int y = pml_e->nodes[ino].y;
							int z = pml_e->nodes[ino].z;
							pml_e->nodes[ino].id = AddPoint(mesh, hash_nodes, p, coords, x, y, z);
						}
						iedge = 11;
						key.id = 10 * (iedge + 1);
						bool lmat = sc_hash_array_lookup(hash_matpml, &key, &position);
						pmlmat_t *pmlT = (pmlmat_t *)sc_array_index(&hash_matpml->a, position);
						pml_e->n_mat = pmlT->mat;
						pmlT->xmin = std::min(pmlT->xmin, x[0]);
						pmlT->xmax = std::max(pmlT->xmax, x[1]);
						pmlT->ymin = 0;
						pmlT->ymax = 0;
						pmlT->zmin = std::min(pmlT->zmin, z[0]);
						pmlT->zmax = std::max(pmlT->zmax, z[4]);
					}
				}
			}
		}

		if (point)
		{
			key.tag = 2;
			if (elem->nodes[0].color == -30 && false)
			{

				for (int nx = 0; nx < layers_x; nx++)
				{
					for (int ny = 0; ny < layers_y; ny++)
					{
						for (int nz = 0; nz < layers_z; nz++)
						{

							octant_t *pml_e = (octant_t *)sc_array_push(&mesh->elements);
							pml_e->id = mesh->elements.elem_count + 1;

							// nos de referencia
							int aux = 0;
							int node0 = elem->nodes[aux].id;
							double x[8], y[8], z[8];

							x[0] = coords[3 * node0 + 0] - nx * X_pml / layers_x;
							x[1] = coords[3 * node0 + 0] - (nx + 1) * X_pml / layers_x;
							x[2] = coords[3 * node0 + 0] - (nx + 1) * X_pml / layers_x;
							x[3] = coords[3 * node0 + 0] - nx * X_pml / layers_x;
							x[4] = coords[3 * node0 + 0] - nx * X_pml / layers_x;
							x[5] = coords[3 * node0 + 0] - (nx + 1) * X_pml / layers_x;
							x[6] = coords[3 * node0 + 0] - (nx + 1) * X_pml / layers_x;
							x[7] = coords[3 * node0 + 0] - nx * X_pml / layers_x;

							pml_e->nodes[0].x = elem->nodes[aux].x - 12 * (nx + 0);
							pml_e->nodes[1].x = elem->nodes[aux].x - 12 * (nx + 1);
							pml_e->nodes[2].x = elem->nodes[aux].x - 12 * (nx + 1);
							pml_e->nodes[3].x = elem->nodes[aux].x - 12 * (nx + 0);
							pml_e->nodes[4].x = elem->nodes[aux].x - 12 * (nx + 0);
							pml_e->nodes[5].x = elem->nodes[aux].x - 12 * (nx + 1);
							pml_e->nodes[6].x = elem->nodes[aux].x - 12 * (nx + 1);
							pml_e->nodes[7].x = elem->nodes[aux].x - 12 * (nx + 0);

							y[0] = coords[3 * node0 + 1] - ny * Y_pml / layers_y;
							y[1] = coords[3 * node0 + 1] - ny * Y_pml / layers_y;
							y[2] = coords[3 * node0 + 1] - (ny + 1) * Y_pml / layers_y;
							y[3] = coords[3 * node0 + 1] - (ny + 1) * Y_pml / layers_y;
							y[4] = coords[3 * node0 + 1] - ny * Y_pml / layers_y;
							y[5] = coords[3 * node0 + 1] - ny * Y_pml / layers_y;
							y[6] = coords[3 * node0 + 1] - (ny + 1) * Y_pml / layers_y;
							y[7] = coords[3 * node0 + 1] - (ny + 1) * Y_pml / layers_y;

							pml_e->nodes[0].y = elem->nodes[aux].y - 12 * (ny + 0);
							pml_e->nodes[1].y = elem->nodes[aux].y - 12 * (ny + 0);
							pml_e->nodes[2].y = elem->nodes[aux].y - 12 * (ny + 1);
							pml_e->nodes[3].y = elem->nodes[aux].y - 12 * (ny + 1);
							pml_e->nodes[4].y = elem->nodes[aux].y - 12 * (ny + 0);
							pml_e->nodes[5].y = elem->nodes[aux].y - 12 * (ny + 0);
							pml_e->nodes[6].y = elem->nodes[aux].y - 12 * (ny + 1);
							pml_e->nodes[7].y = elem->nodes[aux].y - 12 * (ny + 1);

							z[0] = coords[3 * node0 + 2] + nz * Z_pml / layers_z;
							z[1] = coords[3 * node0 + 2] + nz * Z_pml / layers_z;
							z[2] = coords[3 * node0 + 2] + nz * Z_pml / layers_z;
							z[3] = coords[3 * node0 + 2] + nz * Z_pml / layers_z;
							z[4] = coords[3 * node0 + 2] + (nz + 1) * Z_pml / layers_z;
							z[5] = coords[3 * node0 + 2] + (nz + 1) * Z_pml / layers_z;
							z[6] = coords[3 * node0 + 2] + (nz + 1) * Z_pml / layers_z;
							z[7] = coords[3 * node0 + 2] + (nz + 1) * Z_pml / layers_z;

							pml_e->nodes[0].z = elem->nodes[aux].z - 12 * (nz + 1);
							pml_e->nodes[1].z = elem->nodes[aux].z - 12 * (nz + 1);
							pml_e->nodes[2].z = elem->nodes[aux].z - 12 * (nz + 1);
							pml_e->nodes[3].z = elem->nodes[aux].z - 12 * (nz + 1);
							pml_e->nodes[4].z = elem->nodes[aux].z - 12 * (nz + 0);
							pml_e->nodes[5].z = elem->nodes[aux].z - 12 * (nz + 0);
							pml_e->nodes[6].z = elem->nodes[aux].z - 12 * (nz + 0);
							pml_e->nodes[7].z = elem->nodes[aux].z - 12 * (nz + 0);

							for (int ino = 0; ino < 8; ino++)
							{
								// definindo ponto p a ser adicionado
								GtsPoint *p = gts_point_new(gts_point_class(), x[ino], y[ino], z[ino]);
								// adicionando ponto p
								int x = pml_e->nodes[ino].x;
								int y = pml_e->nodes[ino].y;
								int z = pml_e->nodes[ino].z;
								pml_e->nodes[ino].id = AddPoint(mesh, hash_nodes, p, coords, x, y, z);
							}
							icorner = 0;
							key.id = 1000 * (icorner + 1);
							bool lmat = sc_hash_array_lookup(hash_matpml, &key, &position);
							pmlmat_t *pmlT = (pmlmat_t *)sc_array_index(&hash_matpml->a, position);
							pml_e->n_mat = pmlT->mat;
							//it should be updated
							pmlT->xmin = std::min(pmlT->xmin, x[0]);
							pmlT->xmax = std::max(pmlT->xmin, x[1]);
							pmlT->ymin = std::min(pmlT->ymin, y[0]);
							pmlT->ymax = std::min(pmlT->ymin, y[2]);
							pmlT->zmin = std::min(pmlT->zmin, z[4]);
							pmlT->zmax = std::max(pmlT->zmax, z[0]);
						}
					}
				}
			}

			if (elem->nodes[1].color == -31 && false)
			{

				for (int nx = 0; nx < layers_x; nx++)
				{
					for (int ny = 0; ny < layers_y; ny++)
					{
						for (int nz = 0; nz < layers_z; nz++)
						{

							octant_t *pml_e = (octant_t *)sc_array_push(&mesh->elements);
							pml_e->id = mesh->elements.elem_count + 1;

							// nos de referencia
							int aux = 1;
							int node0 = elem->nodes[5].id;
							double x[8], y[8], z[8];

							x[0] = coords[3 * node0 + 0] + nx * X_pml / layers_x;
							x[1] = coords[3 * node0 + 0] + (nx + 1) * X_pml / layers_x;
							x[2] = coords[3 * node0 + 0] + (nx + 1) * X_pml / layers_x;
							x[3] = coords[3 * node0 + 0] + nx * X_pml / layers_x;
							x[4] = coords[3 * node0 + 0] + nx * X_pml / layers_x;
							x[5] = coords[3 * node0 + 0] + (nx + 1) * X_pml / layers_x;
							x[6] = coords[3 * node0 + 0] + (nx + 1) * X_pml / layers_x;
							x[7] = coords[3 * node0 + 0] + nx * X_pml / layers_x;

							pml_e->nodes[0].x = elem->nodes[aux].x + 12 * (nx + 0);
							pml_e->nodes[1].x = elem->nodes[aux].x + 12 * (nx + 1);
							pml_e->nodes[2].x = elem->nodes[aux].x + 12 * (nx + 1);
							pml_e->nodes[3].x = elem->nodes[aux].x + 12 * (nx + 0);
							pml_e->nodes[4].x = elem->nodes[aux].x + 12 * (nx + 0);
							pml_e->nodes[5].x = elem->nodes[aux].x + 12 * (nx + 1);
							pml_e->nodes[6].x = elem->nodes[aux].x + 12 * (nx + 1);
							pml_e->nodes[7].x = elem->nodes[aux].x + 12 * (nx + 0);

							y[0] = coords[3 * node0 + 1] - ny * Y_pml / layers_y;
							y[1] = coords[3 * node0 + 1] - ny * Y_pml / layers_y;
							y[2] = coords[3 * node0 + 1] - (ny + 1) * Y_pml / layers_y;
							y[3] = coords[3 * node0 + 1] - (ny + 1) * Y_pml / layers_y;
							y[4] = coords[3 * node0 + 1] - ny * Y_pml / layers_y;
							y[5] = coords[3 * node0 + 1] - ny * Y_pml / layers_y;
							y[6] = coords[3 * node0 + 1] - (ny + 1) * Y_pml / layers_y;
							y[7] = coords[3 * node0 + 1] - (ny + 1) * Y_pml / layers_y;

							pml_e->nodes[0].y = elem->nodes[aux].y - 12 * (ny + 0);
							pml_e->nodes[1].y = elem->nodes[aux].y - 12 * (ny + 0);
							pml_e->nodes[2].y = elem->nodes[aux].y - 12 * (ny + 1);
							pml_e->nodes[3].y = elem->nodes[aux].y - 12 * (ny + 1);
							pml_e->nodes[4].y = elem->nodes[aux].y - 12 * (ny + 0);
							pml_e->nodes[5].y = elem->nodes[aux].y - 12 * (ny + 0);
							pml_e->nodes[6].y = elem->nodes[aux].y - 12 * (ny + 1);
							pml_e->nodes[7].y = elem->nodes[aux].y - 12 * (ny + 1);

							z[0] = coords[3 * node0 + 2] + nz * Z_pml / layers_z;
							z[1] = coords[3 * node0 + 2] + nz * Z_pml / layers_z;
							z[2] = coords[3 * node0 + 2] + nz * Z_pml / layers_z;
							z[3] = coords[3 * node0 + 2] + nz * Z_pml / layers_z;
							z[4] = coords[3 * node0 + 2] + (nz + 1) * Z_pml / layers_z;
							z[5] = coords[3 * node0 + 2] + (nz + 1) * Z_pml / layers_z;
							z[6] = coords[3 * node0 + 2] + (nz + 1) * Z_pml / layers_z;
							z[7] = coords[3 * node0 + 2] + (nz + 1) * Z_pml / layers_z;

							pml_e->nodes[0].z = elem->nodes[aux].z - 12 * (nz + 1);
							pml_e->nodes[1].z = elem->nodes[aux].z - 12 * (nz + 1);
							pml_e->nodes[2].z = elem->nodes[aux].z - 12 * (nz + 1);
							pml_e->nodes[3].z = elem->nodes[aux].z - 12 * (nz + 1);
							pml_e->nodes[4].z = elem->nodes[aux].z - 12 * (nz + 0);
							pml_e->nodes[5].z = elem->nodes[aux].z - 12 * (nz + 0);
							pml_e->nodes[6].z = elem->nodes[aux].z - 12 * (nz + 0);
							pml_e->nodes[7].z = elem->nodes[aux].z - 12 * (nz + 0);

							for (int ino = 0; ino < 8; ino++)
							{
								// definindo ponto p a ser adicionado
								GtsPoint *p = gts_point_new(gts_point_class(), x[ino], y[ino], z[ino]);
								// adicionando ponto p
								int x = pml_e->nodes[ino].x;
								int y = pml_e->nodes[ino].y;
								int z = pml_e->nodes[ino].z;
								pml_e->nodes[ino].id = AddPoint(mesh, hash_nodes, p, coords, x, y, z);
							}
							icorner = 1;
							key.id = 1000 * (icorner + 1);
							bool lmat = sc_hash_array_lookup(hash_matpml, &key, &position);
							pmlmat_t *pmlT = (pmlmat_t *)sc_array_index(&hash_matpml->a, position);
							pml_e->n_mat = pmlT->mat;
							//it should be updated
							pmlT->xmin = std::min(pmlT->xmin, x[0]);
							pmlT->xmax = std::max(pmlT->xmin, x[1]);
							pmlT->ymin = std::min(pmlT->ymin, y[0]);
							pmlT->ymax = std::min(pmlT->ymin, y[2]);
							pmlT->zmin = std::min(pmlT->zmin, z[4]);
							pmlT->zmax = std::max(pmlT->zmax, z[0]);
						}
					}
				}
			}

			if (elem->nodes[2].color == -32 && false)
			{

				for (int nx = 0; nx < layers_x; nx++)
				{
					for (int ny = 0; ny < layers_y; ny++)
					{
						for (int nz = 0; nz < layers_z; nz++)
						{

							octant_t *pml_e = (octant_t *)sc_array_push(&mesh->elements);
							pml_e->id = mesh->elements.elem_count + 1;

							// nos de referencia
							int aux = 2;
							int node0 = elem->nodes[6].id;
							double x[8], y[8], z[8];

							x[0] = coords[3 * node0 + 0] + nx * X_pml / layers_x;
							x[1] = coords[3 * node0 + 0] + (nx + 1) * X_pml / layers_x;
							x[2] = coords[3 * node0 + 0] + (nx + 1) * X_pml / layers_x;
							x[3] = coords[3 * node0 + 0] + nx * X_pml / layers_x;
							x[4] = coords[3 * node0 + 0] + nx * X_pml / layers_x;
							x[5] = coords[3 * node0 + 0] + (nx + 1) * X_pml / layers_x;
							x[6] = coords[3 * node0 + 0] + (nx + 1) * X_pml / layers_x;
							x[7] = coords[3 * node0 + 0] + nx * X_pml / layers_x;

							pml_e->nodes[0].x = elem->nodes[aux].x + 12 * (nx + 0);
							pml_e->nodes[1].x = elem->nodes[aux].x + 12 * (nx + 1);
							pml_e->nodes[2].x = elem->nodes[aux].x + 12 * (nx + 1);
							pml_e->nodes[3].x = elem->nodes[aux].x + 12 * (nx + 0);
							pml_e->nodes[4].x = elem->nodes[aux].x + 12 * (nx + 0);
							pml_e->nodes[5].x = elem->nodes[aux].x + 12 * (nx + 1);
							pml_e->nodes[6].x = elem->nodes[aux].x + 12 * (nx + 1);
							pml_e->nodes[7].x = elem->nodes[aux].x + 12 * (nx + 0);

							y[0] = coords[3 * node0 + 1] + ny * Y_pml / layers_y;
							y[1] = coords[3 * node0 + 1] + ny * Y_pml / layers_y;
							y[2] = coords[3 * node0 + 1] + (ny + 1) * Y_pml / layers_y;
							y[3] = coords[3 * node0 + 1] + (ny + 1) * Y_pml / layers_y;
							y[4] = coords[3 * node0 + 1] + ny * Y_pml / layers_y;
							y[5] = coords[3 * node0 + 1] + ny * Y_pml / layers_y;
							y[6] = coords[3 * node0 + 1] + (ny + 1) * Y_pml / layers_y;
							y[7] = coords[3 * node0 + 1] + (ny + 1) * Y_pml / layers_y;

							pml_e->nodes[0].y = elem->nodes[aux].y + 12 * (ny + 0);
							pml_e->nodes[1].y = elem->nodes[aux].y + 12 * (ny + 0);
							pml_e->nodes[2].y = elem->nodes[aux].y + 12 * (ny + 1);
							pml_e->nodes[3].y = elem->nodes[aux].y + 12 * (ny + 1);
							pml_e->nodes[4].y = elem->nodes[aux].y + 12 * (ny + 0);
							pml_e->nodes[5].y = elem->nodes[aux].y + 12 * (ny + 0);
							pml_e->nodes[6].y = elem->nodes[aux].y + 12 * (ny + 1);
							pml_e->nodes[7].y = elem->nodes[aux].y + 12 * (ny + 1);

							z[0] = coords[3 * node0 + 2] - nz * Z_pml / layers_z;
							z[1] = coords[3 * node0 + 2] - nz * Z_pml / layers_z;
							z[2] = coords[3 * node0 + 2] - nz * Z_pml / layers_z;
							z[3] = coords[3 * node0 + 2] - nz * Z_pml / layers_z;
							z[4] = coords[3 * node0 + 2] - (nz + 1) * Z_pml / layers_z;
							z[5] = coords[3 * node0 + 2] - (nz + 1) * Z_pml / layers_z;
							z[6] = coords[3 * node0 + 2] - (nz + 1) * Z_pml / layers_z;
							z[7] = coords[3 * node0 + 2] - (nz + 1) * Z_pml / layers_z;

							z[0] = coords[3 * node0 + 2] + nz * Z_pml / layers_z;
							z[1] = coords[3 * node0 + 2] + nz * Z_pml / layers_z;
							z[2] = coords[3 * node0 + 2] + nz * Z_pml / layers_z;
							z[3] = coords[3 * node0 + 2] + nz * Z_pml / layers_z;
							z[4] = coords[3 * node0 + 2] + (nz + 1) * Z_pml / layers_z;
							z[5] = coords[3 * node0 + 2] + (nz + 1) * Z_pml / layers_z;
							z[6] = coords[3 * node0 + 2] + (nz + 1) * Z_pml / layers_z;
							z[7] = coords[3 * node0 + 2] + (nz + 1) * Z_pml / layers_z;

							pml_e->nodes[0].z = elem->nodes[aux].z - 12 * (nz + 1);
							pml_e->nodes[1].z = elem->nodes[aux].z - 12 * (nz + 1);
							pml_e->nodes[2].z = elem->nodes[aux].z - 12 * (nz + 1);
							pml_e->nodes[3].z = elem->nodes[aux].z - 12 * (nz + 1);
							pml_e->nodes[4].z = elem->nodes[aux].z - 12 * (nz + 0);
							pml_e->nodes[5].z = elem->nodes[aux].z - 12 * (nz + 0);
							pml_e->nodes[6].z = elem->nodes[aux].z - 12 * (nz + 0);
							pml_e->nodes[7].z = elem->nodes[aux].z - 12 * (nz + 0);

							for (int ino = 0; ino < 8; ino++)
							{
								// definindo ponto p a ser adicionado
								GtsPoint *p = gts_point_new(gts_point_class(), x[ino], y[ino], z[ino]);
								// adicionando ponto p
								int x = pml_e->nodes[ino].x;
								int y = pml_e->nodes[ino].y;
								int z = pml_e->nodes[ino].z;
								pml_e->nodes[ino].id = AddPoint(mesh, hash_nodes, p, coords, x, y, z);
							}
							icorner = 2;
							key.id = 1000 * (icorner + 1);
							bool lmat = sc_hash_array_lookup(hash_matpml, &key, &position);
							pmlmat_t *pmlT = (pmlmat_t *)sc_array_index(&hash_matpml->a, position);
							pml_e->n_mat = pmlT->mat;
							//it should be updated
							pmlT->xmin = std::min(pmlT->xmin, x[0]);
							pmlT->xmax = std::max(pmlT->xmin, x[1]);
							pmlT->ymin = std::min(pmlT->ymin, y[0]);
							pmlT->ymax = std::min(pmlT->ymin, y[2]);
							pmlT->zmin = std::min(pmlT->zmin, z[4]);
							pmlT->zmax = std::max(pmlT->zmax, z[0]);
						}
					}
				}
			}

			if (elem->nodes[3].color == -33 && false)
			{

				for (int nx = 0; nx < layers_x; nx++)
				{
					for (int ny = 0; ny < layers_y; ny++)
					{
						for (int nz = 0; nz < layers_z; nz++)
						{

							octant_t *pml_e = (octant_t *)sc_array_push(&mesh->elements);
							pml_e->id = mesh->elements.elem_count + 1;

							// nos de referencia
							int aux = 3;
							int node0 = elem->nodes[aux].id;
							double x[8], y[8], z[8];

							x[0] = coords[3 * node0 + 0] - nx * X_pml / layers_x;
							x[1] = coords[3 * node0 + 0] - (nx + 1) * X_pml / layers_x;
							x[2] = coords[3 * node0 + 0] - (nx + 1) * X_pml / layers_x;
							x[3] = coords[3 * node0 + 0] - nx * X_pml / layers_x;
							x[4] = coords[3 * node0 + 0] - nx * X_pml / layers_x;
							x[5] = coords[3 * node0 + 0] - (nx + 1) * X_pml / layers_x;
							x[6] = coords[3 * node0 + 0] - (nx + 1) * X_pml / layers_x;
							x[7] = coords[3 * node0 + 0] - nx * X_pml / layers_x;

							pml_e->nodes[0].x = elem->nodes[aux].x - 12 * (nx + 0);
							pml_e->nodes[1].x = elem->nodes[aux].x - 12 * (nx + 1);
							pml_e->nodes[2].x = elem->nodes[aux].x - 12 * (nx + 1);
							pml_e->nodes[3].x = elem->nodes[aux].x - 12 * (nx + 0);
							pml_e->nodes[4].x = elem->nodes[aux].x - 12 * (nx + 0);
							pml_e->nodes[5].x = elem->nodes[aux].x - 12 * (nx + 1);
							pml_e->nodes[6].x = elem->nodes[aux].x - 12 * (nx + 1);
							pml_e->nodes[7].x = elem->nodes[aux].x - 12 * (nx + 0);

							y[0] = coords[3 * node0 + 1] + ny * Y_pml / layers_y;
							y[1] = coords[3 * node0 + 1] + ny * Y_pml / layers_y;
							y[2] = coords[3 * node0 + 1] + (ny + 1) * Y_pml / layers_y;
							y[3] = coords[3 * node0 + 1] + (ny + 1) * Y_pml / layers_y;
							y[4] = coords[3 * node0 + 1] + ny * Y_pml / layers_y;
							y[5] = coords[3 * node0 + 1] + ny * Y_pml / layers_y;
							y[6] = coords[3 * node0 + 1] + (ny + 1) * Y_pml / layers_y;
							y[7] = coords[3 * node0 + 1] + (ny + 1) * Y_pml / layers_y;

							pml_e->nodes[0].y = elem->nodes[aux].y + 12 * (ny + 0);
							pml_e->nodes[1].y = elem->nodes[aux].y + 12 * (ny + 0);
							pml_e->nodes[2].y = elem->nodes[aux].y + 12 * (ny + 1);
							pml_e->nodes[3].y = elem->nodes[aux].y + 12 * (ny + 1);
							pml_e->nodes[4].y = elem->nodes[aux].y + 12 * (ny + 0);
							pml_e->nodes[5].y = elem->nodes[aux].y + 12 * (ny + 0);
							pml_e->nodes[6].y = elem->nodes[aux].y + 12 * (ny + 1);
							pml_e->nodes[7].y = elem->nodes[aux].y + 12 * (ny + 1);

							z[0] = coords[3 * node0 + 2] + nz * Z_pml / layers_z;
							z[1] = coords[3 * node0 + 2] + nz * Z_pml / layers_z;
							z[2] = coords[3 * node0 + 2] + nz * Z_pml / layers_z;
							z[3] = coords[3 * node0 + 2] + nz * Z_pml / layers_z;
							z[4] = coords[3 * node0 + 2] + (nz + 1) * Z_pml / layers_z;
							z[5] = coords[3 * node0 + 2] + (nz + 1) * Z_pml / layers_z;
							z[6] = coords[3 * node0 + 2] + (nz + 1) * Z_pml / layers_z;
							z[7] = coords[3 * node0 + 2] + (nz + 1) * Z_pml / layers_z;

							pml_e->nodes[0].z = elem->nodes[aux].z - 12 * (nz + 1);
							pml_e->nodes[1].z = elem->nodes[aux].z - 12 * (nz + 1);
							pml_e->nodes[2].z = elem->nodes[aux].z - 12 * (nz + 1);
							pml_e->nodes[3].z = elem->nodes[aux].z - 12 * (nz + 1);
							pml_e->nodes[4].z = elem->nodes[aux].z - 12 * (nz + 0);
							pml_e->nodes[5].z = elem->nodes[aux].z - 12 * (nz + 0);
							pml_e->nodes[6].z = elem->nodes[aux].z - 12 * (nz + 0);
							pml_e->nodes[7].z = elem->nodes[aux].z - 12 * (nz + 0);

							for (int ino = 0; ino < 8; ino++)
							{
								// definindo ponto p a ser adicionado
								GtsPoint *p = gts_point_new(gts_point_class(), x[ino], y[ino], z[ino]);
								// adicionando ponto p
								int x = pml_e->nodes[ino].x;
								int y = pml_e->nodes[ino].y;
								int z = pml_e->nodes[ino].z;
								pml_e->nodes[ino].id = AddPoint(mesh, hash_nodes, p, coords, x, y, z);
							}
							icorner = 3;
							key.id = 1000 * (icorner + 1);
							bool lmat = sc_hash_array_lookup(hash_matpml, &key, &position);
							pmlmat_t *pmlT = (pmlmat_t *)sc_array_index(&hash_matpml->a, position);
							pml_e->n_mat = pmlT->mat;
							//it should be updated
							pmlT->xmin = std::min(pmlT->xmin, x[0]);
							pmlT->xmax = std::max(pmlT->xmin, x[1]);
							pmlT->ymin = std::min(pmlT->ymin, y[0]);
							pmlT->ymax = std::min(pmlT->ymin, y[2]);
							pmlT->zmin = std::min(pmlT->zmin, z[4]);
							pmlT->zmax = std::max(pmlT->zmax, z[0]);
						}
					}
				}
			}

			if (elem->nodes[4].color == -34)
			{

				for (int nx = 0; nx < layers_x; nx++)
				{
					for (int ny = 0; ny < layers_y; ny++)
					{
						for (int nz = 0; nz < layers_z; nz++)
						{

							octant_t *pml_e = (octant_t *)sc_array_push(&mesh->elements);
							pml_e->id = mesh->elements.elem_count + 1;

							// nos de referencia
							int aux = 4;
							int node0 = elem->nodes[aux].id;
							double x[8], y[8], z[8];

							x[0] = coords[3 * node0 + 0] - nx * X_pml / layers_x;
							x[1] = coords[3 * node0 + 0] - (nx + 1) * X_pml / layers_x;
							x[2] = coords[3 * node0 + 0] - (nx + 1) * X_pml / layers_x;
							x[3] = coords[3 * node0 + 0] - nx * X_pml / layers_x;
							x[4] = coords[3 * node0 + 0] - nx * X_pml / layers_x;
							x[5] = coords[3 * node0 + 0] - (nx + 1) * X_pml / layers_x;
							x[6] = coords[3 * node0 + 0] - (nx + 1) * X_pml / layers_x;
							x[7] = coords[3 * node0 + 0] - nx * X_pml / layers_x;

							pml_e->nodes[0].x = elem->nodes[aux].x - 12 * (nx + 0);
							pml_e->nodes[1].x = elem->nodes[aux].x - 12 * (nx + 1);
							pml_e->nodes[2].x = elem->nodes[aux].x - 12 * (nx + 1);
							pml_e->nodes[3].x = elem->nodes[aux].x - 12 * (nx + 0);
							pml_e->nodes[4].x = elem->nodes[aux].x - 12 * (nx + 0);
							pml_e->nodes[5].x = elem->nodes[aux].x - 12 * (nx + 1);
							pml_e->nodes[6].x = elem->nodes[aux].x - 12 * (nx + 1);
							pml_e->nodes[7].x = elem->nodes[aux].x - 12 * (nx + 0);

							y[0] = coords[3 * node0 + 1] - ny * Y_pml / layers_y;
							y[1] = coords[3 * node0 + 1] - ny * Y_pml / layers_y;
							y[2] = coords[3 * node0 + 1] - (ny + 1) * Y_pml / layers_y;
							y[3] = coords[3 * node0 + 1] - (ny + 1) * Y_pml / layers_y;
							y[4] = coords[3 * node0 + 1] - ny * Y_pml / layers_y;
							y[5] = coords[3 * node0 + 1] - ny * Y_pml / layers_y;
							y[6] = coords[3 * node0 + 1] - (ny + 1) * Y_pml / layers_y;
							y[7] = coords[3 * node0 + 1] - (ny + 1) * Y_pml / layers_y;

							pml_e->nodes[0].y = elem->nodes[aux].y - 12 * (ny + 0);
							pml_e->nodes[1].y = elem->nodes[aux].y - 12 * (ny + 0);
							pml_e->nodes[2].y = elem->nodes[aux].y - 12 * (ny + 1);
							pml_e->nodes[3].y = elem->nodes[aux].y - 12 * (ny + 1);
							pml_e->nodes[4].y = elem->nodes[aux].y - 12 * (ny + 0);
							pml_e->nodes[5].y = elem->nodes[aux].y - 12 * (ny + 0);
							pml_e->nodes[6].y = elem->nodes[aux].y - 12 * (ny + 1);
							pml_e->nodes[7].y = elem->nodes[aux].y - 12 * (ny + 1);

							z[0] = coords[3 * node0 + 2] - (nz + 1) * Z_pml / layers_z;
							z[1] = coords[3 * node0 + 2] - (nz + 1) * Z_pml / layers_z;
							z[2] = coords[3 * node0 + 2] - (nz + 1) * Z_pml / layers_z;
							z[3] = coords[3 * node0 + 2] - (nz + 1) * Z_pml / layers_z;
							z[4] = coords[3 * node0 + 2] - (nz)*Z_pml / layers_z;
							z[5] = coords[3 * node0 + 2] - (nz)*Z_pml / layers_z;
							z[6] = coords[3 * node0 + 2] - (nz)*Z_pml / layers_z;
							z[7] = coords[3 * node0 + 2] - (nz)*Z_pml / layers_z;

							pml_e->nodes[0].z = elem->nodes[aux].z + 12 * (nz + 1);
							pml_e->nodes[1].z = elem->nodes[aux].z + 12 * (nz + 1);
							pml_e->nodes[2].z = elem->nodes[aux].z + 12 * (nz + 1);
							pml_e->nodes[3].z = elem->nodes[aux].z + 12 * (nz + 1);
							pml_e->nodes[4].z = elem->nodes[aux].z + 12 * (nz + 0);
							pml_e->nodes[5].z = elem->nodes[aux].z + 12 * (nz + 0);
							pml_e->nodes[6].z = elem->nodes[aux].z + 12 * (nz + 0);
							pml_e->nodes[7].z = elem->nodes[aux].z + 12 * (nz + 0);

							for (int ino = 0; ino < 8; ino++)
							{
								// definindo ponto p a ser adicionado
								GtsPoint *p = gts_point_new(gts_point_class(), x[ino], y[ino], z[ino]);
								// adicionando ponto p
								int x = pml_e->nodes[ino].x;
								int y = pml_e->nodes[ino].y;
								int z = pml_e->nodes[ino].z;
								pml_e->nodes[ino].id = AddPoint(mesh, hash_nodes, p, coords, x, y, z);
							}
							icorner = 4;
							key.id = 1000 * (icorner + 1);
							bool lmat = sc_hash_array_lookup(hash_matpml, &key, &position);
							pmlmat_t *pmlT = (pmlmat_t *)sc_array_index(&hash_matpml->a, position);
							pml_e->n_mat = pmlT->mat;
							pmlT->xmin = std::min(pmlT->xmin, x[1]);
							pmlT->xmax = std::max(pmlT->xmax, x[0]);
							pmlT->ymin = std::min(pmlT->ymin, y[2]);
							pmlT->ymax = std::max(pmlT->ymax, y[0]);
							pmlT->zmin = std::min(pmlT->zmin, z[0]);
							pmlT->zmax = std::max(pmlT->zmax, z[4]);
						}
					}
				}
			}

			if (elem->nodes[5].color == -35)
			{

				for (int nx = 0; nx < layers_x; nx++)
				{
					for (int ny = 0; ny < layers_y; ny++)
					{
						for (int nz = 0; nz < layers_z; nz++)
						{

							octant_t *pml_e = (octant_t *)sc_array_push(&mesh->elements);
							pml_e->id = mesh->elements.elem_count + 1;

							// nos de referencia
							int aux = 5;
							int node0 = elem->nodes[5].id;
							double x[8], y[8], z[8];

							x[0] = coords[3 * node0 + 0] + nx * X_pml / layers_x;
							x[1] = coords[3 * node0 + 0] + (nx + 1) * X_pml / layers_x;
							x[2] = coords[3 * node0 + 0] + (nx + 1) * X_pml / layers_x;
							x[3] = coords[3 * node0 + 0] + nx * X_pml / layers_x;
							x[4] = coords[3 * node0 + 0] + nx * X_pml / layers_x;
							x[5] = coords[3 * node0 + 0] + (nx + 1) * X_pml / layers_x;
							x[6] = coords[3 * node0 + 0] + (nx + 1) * X_pml / layers_x;
							x[7] = coords[3 * node0 + 0] + nx * X_pml / layers_x;

							pml_e->nodes[0].x = elem->nodes[aux].x + 12 * (nx + 0);
							pml_e->nodes[1].x = elem->nodes[aux].x + 12 * (nx + 1);
							pml_e->nodes[2].x = elem->nodes[aux].x + 12 * (nx + 1);
							pml_e->nodes[3].x = elem->nodes[aux].x + 12 * (nx + 0);
							pml_e->nodes[4].x = elem->nodes[aux].x + 12 * (nx + 0);
							pml_e->nodes[5].x = elem->nodes[aux].x + 12 * (nx + 1);
							pml_e->nodes[6].x = elem->nodes[aux].x + 12 * (nx + 1);
							pml_e->nodes[7].x = elem->nodes[aux].x + 12 * (nx + 0);

							y[0] = coords[3 * node0 + 1] - ny * Y_pml / layers_y;
							y[1] = coords[3 * node0 + 1] - ny * Y_pml / layers_y;
							y[2] = coords[3 * node0 + 1] - (ny + 1) * Y_pml / layers_y;
							y[3] = coords[3 * node0 + 1] - (ny + 1) * Y_pml / layers_y;
							y[4] = coords[3 * node0 + 1] - ny * Y_pml / layers_y;
							y[5] = coords[3 * node0 + 1] - ny * Y_pml / layers_y;
							y[6] = coords[3 * node0 + 1] - (ny + 1) * Y_pml / layers_y;
							y[7] = coords[3 * node0 + 1] - (ny + 1) * Y_pml / layers_y;

							pml_e->nodes[0].y = elem->nodes[aux].y - 12 * (ny + 0);
							pml_e->nodes[1].y = elem->nodes[aux].y - 12 * (ny + 0);
							pml_e->nodes[2].y = elem->nodes[aux].y - 12 * (ny + 1);
							pml_e->nodes[3].y = elem->nodes[aux].y - 12 * (ny + 1);
							pml_e->nodes[4].y = elem->nodes[aux].y - 12 * (ny + 0);
							pml_e->nodes[5].y = elem->nodes[aux].y - 12 * (ny + 0);
							pml_e->nodes[6].y = elem->nodes[aux].y - 12 * (ny + 1);
							pml_e->nodes[7].y = elem->nodes[aux].y - 12 * (ny + 1);

							z[0] = coords[3 * node0 + 2] - (nz + 1) * Z_pml / layers_z;
							z[1] = coords[3 * node0 + 2] - (nz + 1) * Z_pml / layers_z;
							z[2] = coords[3 * node0 + 2] - (nz + 1) * Z_pml / layers_z;
							z[3] = coords[3 * node0 + 2] - (nz + 1) * Z_pml / layers_z;
							z[4] = coords[3 * node0 + 2] - (nz)*Z_pml / layers_z;
							z[5] = coords[3 * node0 + 2] - (nz)*Z_pml / layers_z;
							z[6] = coords[3 * node0 + 2] - (nz)*Z_pml / layers_z;
							z[7] = coords[3 * node0 + 2] - (nz)*Z_pml / layers_z;

							pml_e->nodes[0].z = elem->nodes[aux].z + 12 * (nz + 1);
							pml_e->nodes[1].z = elem->nodes[aux].z + 12 * (nz + 1);
							pml_e->nodes[2].z = elem->nodes[aux].z + 12 * (nz + 1);
							pml_e->nodes[3].z = elem->nodes[aux].z + 12 * (nz + 1);
							pml_e->nodes[4].z = elem->nodes[aux].z + 12 * (nz + 0);
							pml_e->nodes[5].z = elem->nodes[aux].z + 12 * (nz + 0);
							pml_e->nodes[6].z = elem->nodes[aux].z + 12 * (nz + 0);
							pml_e->nodes[7].z = elem->nodes[aux].z + 12 * (nz + 0);

							for (int ino = 0; ino < 8; ino++)
							{
								// definindo ponto p a ser adicionado
								GtsPoint *p = gts_point_new(gts_point_class(), x[ino], y[ino], z[ino]);
								// adicionando ponto p
								int x = pml_e->nodes[ino].x;
								int y = pml_e->nodes[ino].y;
								int z = pml_e->nodes[ino].z;
								pml_e->nodes[ino].id = AddPoint(mesh, hash_nodes, p, coords, x, y, z);
							}
							icorner = 5;
							key.id = 1000 * (icorner + 1);
							bool lmat = sc_hash_array_lookup(hash_matpml, &key, &position);
							pmlmat_t *pmlT = (pmlmat_t *)sc_array_index(&hash_matpml->a, position);
							pml_e->n_mat = pmlT->mat;
							pmlT->xmin = std::min(pmlT->xmin, x[0]);
							pmlT->xmax = std::max(pmlT->xmax, x[1]);
							pmlT->ymin = std::min(pmlT->ymin, y[2]);
							pmlT->ymax = std::max(pmlT->ymax, y[0]);
							pmlT->zmin = std::min(pmlT->zmin, z[0]);
							pmlT->zmax = std::max(pmlT->zmax, z[4]);
						}
					}
				}
			}

			if (elem->nodes[6].color == -36)
			{

				for (int nx = 0; nx < layers_x; nx++)
				{
					for (int ny = 0; ny < layers_y; ny++)
					{
						for (int nz = 0; nz < layers_z; nz++)
						{

							octant_t *pml_e = (octant_t *)sc_array_push(&mesh->elements);
							pml_e->id = mesh->elements.elem_count + 1;

							// nos de referencia
							int aux = 6;
							int node0 = elem->nodes[6].id;
							double x[8], y[8], z[8];

							x[0] = coords[3 * node0 + 0] + nx * X_pml / layers_x;
							x[1] = coords[3 * node0 + 0] + (nx + 1) * X_pml / layers_x;
							x[2] = coords[3 * node0 + 0] + (nx + 1) * X_pml / layers_x;
							x[3] = coords[3 * node0 + 0] + nx * X_pml / layers_x;
							x[4] = coords[3 * node0 + 0] + nx * X_pml / layers_x;
							x[5] = coords[3 * node0 + 0] + (nx + 1) * X_pml / layers_x;
							x[6] = coords[3 * node0 + 0] + (nx + 1) * X_pml / layers_x;
							x[7] = coords[3 * node0 + 0] + nx * X_pml / layers_x;

							pml_e->nodes[0].x = elem->nodes[aux].x + 12 * (nx + 0);
							pml_e->nodes[1].x = elem->nodes[aux].x + 12 * (nx + 1);
							pml_e->nodes[2].x = elem->nodes[aux].x + 12 * (nx + 1);
							pml_e->nodes[3].x = elem->nodes[aux].x + 12 * (nx + 0);
							pml_e->nodes[4].x = elem->nodes[aux].x + 12 * (nx + 0);
							pml_e->nodes[5].x = elem->nodes[aux].x + 12 * (nx + 1);
							pml_e->nodes[6].x = elem->nodes[aux].x + 12 * (nx + 1);
							pml_e->nodes[7].x = elem->nodes[aux].x + 12 * (nx + 0);

							y[0] = coords[3 * node0 + 1] + ny * Y_pml / layers_y;
							y[1] = coords[3 * node0 + 1] + ny * Y_pml / layers_y;
							y[2] = coords[3 * node0 + 1] + (ny + 1) * Y_pml / layers_y;
							y[3] = coords[3 * node0 + 1] + (ny + 1) * Y_pml / layers_y;
							y[4] = coords[3 * node0 + 1] + ny * Y_pml / layers_y;
							y[5] = coords[3 * node0 + 1] + ny * Y_pml / layers_y;
							y[6] = coords[3 * node0 + 1] + (ny + 1) * Y_pml / layers_y;
							y[7] = coords[3 * node0 + 1] + (ny + 1) * Y_pml / layers_y;

							pml_e->nodes[0].y = elem->nodes[aux].y + 12 * (ny + 0);
							pml_e->nodes[1].y = elem->nodes[aux].y + 12 * (ny + 0);
							pml_e->nodes[2].y = elem->nodes[aux].y + 12 * (ny + 1);
							pml_e->nodes[3].y = elem->nodes[aux].y + 12 * (ny + 1);
							pml_e->nodes[4].y = elem->nodes[aux].y + 12 * (ny + 0);
							pml_e->nodes[5].y = elem->nodes[aux].y + 12 * (ny + 0);
							pml_e->nodes[6].y = elem->nodes[aux].y + 12 * (ny + 1);
							pml_e->nodes[7].y = elem->nodes[aux].y + 12 * (ny + 1);

							z[0] = coords[3 * node0 + 2] - (nz + 1) * Z_pml / layers_z;
							z[1] = coords[3 * node0 + 2] - (nz + 1) * Z_pml / layers_z;
							z[2] = coords[3 * node0 + 2] - (nz + 1) * Z_pml / layers_z;
							z[3] = coords[3 * node0 + 2] - (nz + 1) * Z_pml / layers_z;
							z[4] = coords[3 * node0 + 2] - (nz)*Z_pml / layers_z;
							z[5] = coords[3 * node0 + 2] - (nz)*Z_pml / layers_z;
							z[6] = coords[3 * node0 + 2] - (nz)*Z_pml / layers_z;
							z[7] = coords[3 * node0 + 2] - (nz)*Z_pml / layers_z;

							pml_e->nodes[0].z = elem->nodes[aux].z + 12 * (nz + 1);
							pml_e->nodes[1].z = elem->nodes[aux].z + 12 * (nz + 1);
							pml_e->nodes[2].z = elem->nodes[aux].z + 12 * (nz + 1);
							pml_e->nodes[3].z = elem->nodes[aux].z + 12 * (nz + 1);
							pml_e->nodes[4].z = elem->nodes[aux].z + 12 * (nz + 0);
							pml_e->nodes[5].z = elem->nodes[aux].z + 12 * (nz + 0);
							pml_e->nodes[6].z = elem->nodes[aux].z + 12 * (nz + 0);
							pml_e->nodes[7].z = elem->nodes[aux].z + 12 * (nz + 0);

							for (int ino = 0; ino < 8; ino++)
							{
								// definindo ponto p a ser adicionado
								GtsPoint *p = gts_point_new(gts_point_class(), x[ino], y[ino], z[ino]);
								// adicionando ponto p
								int x = pml_e->nodes[ino].x;
								int y = pml_e->nodes[ino].y;
								int z = pml_e->nodes[ino].z;
								pml_e->nodes[ino].id = AddPoint(mesh, hash_nodes, p, coords, x, y, z);
							}
							icorner = 6;
							key.id = 1000 * (icorner + 1);
							bool lmat = sc_hash_array_lookup(hash_matpml, &key, &position);
							pmlmat_t *pmlT = (pmlmat_t *)sc_array_index(&hash_matpml->a, position);
							pml_e->n_mat = pmlT->mat;
							pmlT->xmin = std::min(pmlT->xmin, x[0]);
							pmlT->xmax = std::max(pmlT->xmax, x[1]);
							pmlT->ymin = std::min(pmlT->ymin, y[0]);
							pmlT->ymax = std::max(pmlT->ymax, y[2]);
							pmlT->zmin = std::min(pmlT->zmin, z[0]);
							pmlT->zmax = std::max(pmlT->zmax, z[4]);
						}
					}
				}
			}

			if (elem->nodes[7].color == -37)
			{

				for (int nx = 0; nx < layers_x; nx++)
				{
					for (int ny = 0; ny < layers_y; ny++)
					{
						for (int nz = 0; nz < layers_z; nz++)
						{

							octant_t *pml_e = (octant_t *)sc_array_push(&mesh->elements);
							pml_e->id = mesh->elements.elem_count + 1;

							// nos de referencia
							int aux = 7;
							int node0 = elem->nodes[aux].id;
							double x[8], y[8], z[8];

							x[0] = coords[3 * node0 + 0] - nx * X_pml / layers_x;
							x[1] = coords[3 * node0 + 0] - (nx + 1) * X_pml / layers_x;
							x[2] = coords[3 * node0 + 0] - (nx + 1) * X_pml / layers_x;
							x[3] = coords[3 * node0 + 0] - nx * X_pml / layers_x;
							x[4] = coords[3 * node0 + 0] - nx * X_pml / layers_x;
							x[5] = coords[3 * node0 + 0] - (nx + 1) * X_pml / layers_x;
							x[6] = coords[3 * node0 + 0] - (nx + 1) * X_pml / layers_x;
							x[7] = coords[3 * node0 + 0] - nx * X_pml / layers_x;

							pml_e->nodes[0].x = elem->nodes[aux].x - 12 * (nx + 0);
							pml_e->nodes[1].x = elem->nodes[aux].x - 12 * (nx + 1);
							pml_e->nodes[2].x = elem->nodes[aux].x - 12 * (nx + 1);
							pml_e->nodes[3].x = elem->nodes[aux].x - 12 * (nx + 0);
							pml_e->nodes[4].x = elem->nodes[aux].x - 12 * (nx + 0);
							pml_e->nodes[5].x = elem->nodes[aux].x - 12 * (nx + 1);
							pml_e->nodes[6].x = elem->nodes[aux].x - 12 * (nx + 1);
							pml_e->nodes[7].x = elem->nodes[aux].x - 12 * (nx + 0);

							y[0] = coords[3 * node0 + 1] + ny * Y_pml / layers_y;
							y[1] = coords[3 * node0 + 1] + ny * Y_pml / layers_y;
							y[2] = coords[3 * node0 + 1] + (ny + 1) * Y_pml / layers_y;
							y[3] = coords[3 * node0 + 1] + (ny + 1) * Y_pml / layers_y;
							y[4] = coords[3 * node0 + 1] + ny * Y_pml / layers_y;
							y[5] = coords[3 * node0 + 1] + ny * Y_pml / layers_y;
							y[6] = coords[3 * node0 + 1] + (ny + 1) * Y_pml / layers_y;
							y[7] = coords[3 * node0 + 1] + (ny + 1) * Y_pml / layers_y;

							pml_e->nodes[0].y = elem->nodes[aux].y + 12 * (ny + 0);
							pml_e->nodes[1].y = elem->nodes[aux].y + 12 * (ny + 0);
							pml_e->nodes[2].y = elem->nodes[aux].y + 12 * (ny + 1);
							pml_e->nodes[3].y = elem->nodes[aux].y + 12 * (ny + 1);
							pml_e->nodes[4].y = elem->nodes[aux].y + 12 * (ny + 0);
							pml_e->nodes[5].y = elem->nodes[aux].y + 12 * (ny + 0);
							pml_e->nodes[6].y = elem->nodes[aux].y + 12 * (ny + 1);
							pml_e->nodes[7].y = elem->nodes[aux].y + 12 * (ny + 1);

							z[0] = coords[3 * node0 + 2] - (nz + 1) * Z_pml / layers_z;
							z[1] = coords[3 * node0 + 2] - (nz + 1) * Z_pml / layers_z;
							z[2] = coords[3 * node0 + 2] - (nz + 1) * Z_pml / layers_z;
							z[3] = coords[3 * node0 + 2] - (nz + 1) * Z_pml / layers_z;
							z[4] = coords[3 * node0 + 2] - (nz)*Z_pml / layers_z;
							z[5] = coords[3 * node0 + 2] - (nz)*Z_pml / layers_z;
							z[6] = coords[3 * node0 + 2] - (nz)*Z_pml / layers_z;
							z[7] = coords[3 * node0 + 2] - (nz)*Z_pml / layers_z;

							pml_e->nodes[0].z = elem->nodes[aux].z + 12 * (nz + 1);
							pml_e->nodes[1].z = elem->nodes[aux].z + 12 * (nz + 1);
							pml_e->nodes[2].z = elem->nodes[aux].z + 12 * (nz + 1);
							pml_e->nodes[3].z = elem->nodes[aux].z + 12 * (nz + 1);
							pml_e->nodes[4].z = elem->nodes[aux].z + 12 * (nz + 0);
							pml_e->nodes[5].z = elem->nodes[aux].z + 12 * (nz + 0);
							pml_e->nodes[6].z = elem->nodes[aux].z + 12 * (nz + 0);
							pml_e->nodes[7].z = elem->nodes[aux].z + 12 * (nz + 0);

							for (int ino = 0; ino < 8; ino++)
							{
								// definindo ponto p a ser adicionado
								GtsPoint *p = gts_point_new(gts_point_class(), x[ino], y[ino], z[ino]);
								// adicionando ponto p
								int x = pml_e->nodes[ino].x;
								int y = pml_e->nodes[ino].y;
								int z = pml_e->nodes[ino].z;
								pml_e->nodes[ino].id = AddPoint(mesh, hash_nodes, p, coords, x, y, z);
							}
							icorner = 7;
							key.id = 1000 * (icorner + 1);
							bool lmat = sc_hash_array_lookup(hash_matpml, &key, &position);
							pmlmat_t *pmlT = (pmlmat_t *)sc_array_index(&hash_matpml->a, position);
							pml_e->n_mat = pmlT->mat;
							pmlT->xmin = std::min(pmlT->xmin, x[1]);
							pmlT->xmax = std::max(pmlT->xmax, x[0]);
							pmlT->ymin = std::min(pmlT->ymin, y[0]);
							pmlT->ymax = std::max(pmlT->ymax, y[2]);
							pmlT->zmin = std::min(pmlT->zmin, z[0]);
							pmlT->zmax = std::max(pmlT->zmax, z[4]);
						}
					}
				}
			}
		}

		sc_array_reset(&toto);
	}

	// update the vectors
	mesh->local_n_elements = mesh->elements.elem_count;
	mesh->local_n_nodes = mesh->nodes.elem_count;
	MPI_Allreduce(&mesh->local_n_elements, &mesh->total_n_elements, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(&mesh->local_n_nodes, &mesh->total_n_nodes, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);

	mesh->ncellx = mesh->ncellx + 2 * 4 * layers_x;
	mesh->ncelly = mesh->ncelly + 2 * 4 * layers_y;
	mesh->max_z = mesh->max_z + 4 * layers_y;

	free(mesh->part_nodes);
	mesh->part_nodes = (int *)malloc(mesh->local_n_nodes * sizeof(int));
	for (int ino = 0; ino < mesh->local_n_nodes; ino++)
	{
		mesh->part_nodes[ino] = mesh->mpi_rank;
	}

	//printf(" Ajust material properties\n\n");
	//Adjust_material(mesh);

	if (mesh->mpi_rank == 0)
	{
		printf("Total number of elements: %d\n", mesh->total_n_elements);
		printf("Total number of nodes: %d\n", mesh->total_n_nodes);
		printf("Total number of materias: %d\n",hash_matpml->a.elem_count);
	}

	// write pml file
	fprintf(fp, "Test\n");
	for(int i = 0; i < hash_matpml->a.elem_count; i++){
		pmlmat_t *pmlT = (pmlmat_t *)sc_array_index(&hash_matpml->a, i);
		if(pmlT->id < 8){
			fprintf(fp, "Hash id: %d, mat: %d, matref: %d\n", pmlT->id, pmlT->mat, pmlT->matref);
			fprintf(fp, "xmin: %f,  xmax: %f\n", pmlT->xmin, pmlT->xmax);
			fprintf(fp, "ymin: %f,  ymax: %f\n", pmlT->ymin, pmlT->ymax);
			fprintf(fp, "zmin: %f,  zmax: %f\n", pmlT->zmin, pmlT->zmax);
		}
	}
	for(int i = 0; i < hash_matpml->a.elem_count; i++){
		pmlmat_t *pmlT = (pmlmat_t *)sc_array_index(&hash_matpml->a, i);
		if(pmlT->id < 200 && pmlT->id > 7){
			fprintf(fp, "Hash id: %d, mat: %d, matref: %d\n", pmlT->id, pmlT->mat, pmlT->matref);
			fprintf(fp, "xmin: %f,  xmax: %f\n", pmlT->xmin, pmlT->xmax);
			fprintf(fp, "ymin: %f,  ymax: %f\n", pmlT->ymin, pmlT->ymax);
			fprintf(fp, "zmin: %f,  zmax: %f\n", pmlT->zmin, pmlT->zmax);
		}
	}
	for(int i = 0; i < hash_matpml->a.elem_count; i++){
		pmlmat_t *pmlT = (pmlmat_t *)sc_array_index(&hash_matpml->a, i);
		if(pmlT->id > 299){
			fprintf(fp, "Hash id: %d, mat: %d, matref: %d\n", pmlT->id, pmlT->mat, pmlT->matref);
			fprintf(fp, "xmin: %f,  xmax: %f\n", pmlT->xmin, pmlT->xmax);
			fprintf(fp, "ymin: %f,  ymax: %f\n", pmlT->ymin, pmlT->ymax);
			fprintf(fp, "zmin: %f,  zmax: %f\n", pmlT->zmin, pmlT->zmax);
		}
	}
	fclose(fp);


	// material.input file 2 SEM3D
	FILE *fp1;
	fp1 = fopen("material.input.FromHexMesh", "w");
	if (fp1 == NULL)
	{
		printf("Error opening material.input file\n");
	}

	double vp, vs, rho;
	vp = 6300;
	vs = 2300;
	rho = 5000;
	tot_n_mat++;
	fprintf(fp1, "%d\n", tot_n_mat+hash_matpml->a.elem_count);
	for(int imat = 0; imat < tot_n_mat; imat++){
		fprintf(fp1, "S %f %f %f %f %f\n",vp,vs,rho,0.0,0.0 );
	}
	for(int imat = 0; imat < hash_matpml->a.elem_count; imat++){
		fprintf(fp1, "P %f %f %f %f %f\n",vp,vs,rho,0.0,0.0);
	}
	fprintf(fp1, "# PML properties\n");
	fprintf(fp1, "# npow,Apow,posX,widthX,posY,widthY,posZ,widthZ,mat\n");

	for(int imat = 0; imat < hash_matpml->a.elem_count; imat++){
		pmlmat_t *pmlT = (pmlmat_t *)sc_array_index(&hash_matpml->a, imat);
		double xx = pmlT->xmin;
		double yy = pmlT->ymin;
		double zz = pmlT->zmin;
		double dx = pmlT->xmax - pmlT->xmin;
		double dy = pmlT->ymax - pmlT->ymin;
		double dz = pmlT->zmax - pmlT->zmin;
		fprintf(fp1, "2 10.000000 %f %f %f %f %f %f %d\n",xx,dx,yy,dy,zz,dz,pmlT->matref+1);
	}

	fclose(fp1);

}
