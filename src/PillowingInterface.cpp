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

#include "hexa.h"
#include "hilbert.h"

#include <ctime>

void PillowingInterface(hexa_tree_t* mesh, std::vector<double>& coords, std::vector<int>& nodes_b_mat){

	//criando a hash de nos para evitar nos duplicados no pillowing
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

	//criando a superficie com as normais para cada ponto do nodes_b_material

	//fazendo vertex hash -> assim consigo acessar os vizinhos num stencil "simples"
	sc_hash_array_t*	indep_vertex    = (sc_hash_array_t *)sc_hash_array_new(sizeof (octant_vertex_t), vertex_hash_id, vertex_equal_id, &clamped);
	//TODO only for the octree structure
	for (int iel = 0; iel < mesh->elements.elem_count; ++iel) {
		size_t  position;
		octant_t *elem = (octant_t*) sc_array_index(&mesh->elements, iel);

		for (int ino = 0; ino < 8; ino++){
			octant_vertex_t key;
			key.id = elem->nodes[ino].id;
			octant_vertex_t* vert = (octant_vertex_t*) sc_hash_array_insert_unique (indep_vertex, &key, &position);
			if(vert != NULL){
				vert->id = elem->nodes[ino].id;
				vert->list_elem = 1;
				vert->elem[vert->list_elem-1] = elem->id;
			}else{
				vert = (octant_vertex_t*) sc_array_index(&indep_vertex->a, position);
				vert->elem[vert->list_elem] = elem->id;
				vert->list_elem++;
			}
		}
	}

	//fazendo hash nos nodes_b_mat
	sc_hash_array_t* hash_b_mat = sc_hash_array_new(sizeof(node_t), edge_hash_fn, edge_equal_fn, &clamped);

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
		}else{
			printf("Verificar o no numero %d\n",node);
		}
	}

	int vert_ord[4][3]= {
			{0, 1, 3},
			{1, 2, 0},
			{2, 3, 1},
			{3, 0, 2}
	};

	// std::vector<double> norm;
	sc_hash_array *hash_normal = sc_hash_array_new(sizeof(normal_t), vertex_hash_id, vertex_equal_id, &clamped);
	double norm_mag = 0.0;
	//entrando em cada vertice de nodes_b_mat e encontrando o vertice,
	// acessando a lista de elementos para calculo da normal
	for(int ino = 0; ino < nodes_b_mat.size(); ino++){

		octant_vertex_t key;
		size_t  position;
		key.id = nodes_b_mat[ino];
		bool tre = (octant_vertex_t*) sc_hash_array_lookup (indep_vertex, &key, &position);

		if(tre){
			octant_vertex_t * vert = (octant_vertex_t*) sc_array_index(&indep_vertex->a, position);
			//printf("achei o no numero:%d\n",vert->id);
			for(int iel = 0; iel < vert->list_elem; iel++){
				octant_t *elem = (octant_t*) sc_array_index(&mesh->elements, vert->elem[iel]);

				//for(int isurf = 0; isurf < 6; isurf += 2){
				for(int isurf = 0; isurf < 6; isurf++){

					if(((elem->nodes[FaceNodesMap[isurf][0]].id == nodes_b_mat[ino]) || (elem->nodes[FaceNodesMap[isurf][1]].id == nodes_b_mat[ino])
							|| (elem->nodes[FaceNodesMap[isurf][2]].id == nodes_b_mat[ino])  || (elem->nodes[FaceNodesMap[isurf][3]].id == nodes_b_mat[ino])) &&
							(elem->nodes[FaceNodesMap[isurf][0]].fixed==1 && elem->nodes[FaceNodesMap[isurf][1]].fixed==1
									&& elem->nodes[FaceNodesMap[isurf][2]].fixed==1 && elem->nodes[FaceNodesMap[isurf][3]].fixed==1)){

						int ive;
						for (int iv = 0 ; iv <4; iv++){
							if(elem->nodes[FaceNodesMap[isurf][iv]].id == nodes_b_mat[ino]){
								ive = 	iv;
							}
						}

						int node0 = elem->nodes[FaceNodesMap[isurf][vert_ord[ive][0]]].id;
						int node1 = elem->nodes[FaceNodesMap[isurf][vert_ord[ive][1]]].id;
						int node2 = elem->nodes[FaceNodesMap[isurf][vert_ord[ive][2]]].id;

						double ax = coords[3*node0 +0] - coords[3*node1 +0];
						double ay = coords[3*node0 +1] - coords[3*node1 +1];
						double az = coords[3*node0 +2] - coords[3*node1 +2];

						double bx = coords[3*node0 +0] - coords[3*node2 +0];
						double by = coords[3*node0 +1] - coords[3*node2 +1];
						double bz = coords[3*node0 +2] - coords[3*node2 +2];

						double nx = -(ay*bz - az*by);
						double ny = -(az*bx - ax*bz);
						double nz = -(ax*by - ay*bx);

						norm_mag = max(norm_mag,sqrt(nx*nx+ny*ny+nz*nz));

						size_t position;
						normal_t key;
						key.id = nodes_b_mat[ino];
						normal_t * normal = (normal_t*) sc_hash_array_insert_unique(hash_normal, &key, &position);

						if(normal!=NULL){
							normal->id = nodes_b_mat[ino];
							normal->list_elem = 1;
							normal->list_face = 1;
							normal->elem[normal->list_face-1] = elem->id;
							normal->face[normal->list_face-1] = isurf;
							normal->n[normal->list_face-1][1] = nx;
							normal->n[normal->list_face-1][2] = ny;
							normal->n[normal->list_face-1][3] = nz;
						}else{
							normal = (normal_t*) sc_array_index(&hash_normal->a, position);
							normal->elem[normal->list_elem] = elem->id;
							normal->face[normal->list_face] = isurf;
							normal->n[normal->list_face][1] = nx;
							normal->n[normal->list_face][2] = ny;
							normal->n[normal->list_face][3] = nz;
							normal->list_elem++;
							normal->list_face++;
						}

						//printf("%d %d %d %f %f %f\n",key.id, vert->elem[iel], isurf,nx,ny,nz );
					}
				}
			}
		}else{

		}
	}

	printf("norma:%f\n",norm_mag);
	double nx = 0;
	double ny = 0;
	double nz = 0;
	//making the mean value of the surface normal
	//TODO check if the normal is good with all the surfaces
	for(int ino = 0; ino > hash_normal->a.elem_count; ino ++){
		normal_t* normal = (normal_t*) sc_array_index (&hash_normal->a, ino);
		for(int isurf = 0; isurf < normal->list_face; isurf++){
			nx += normal->n[isurf][0]/norm_mag;
			ny += normal->n[isurf][1]/norm_mag;
			nz += normal->n[isurf][2]/norm_mag;
		}
		normal->nm[0] = nx/normal->list_face;
		normal->nm[1] = ny/normal->list_face;
		normal->nm[2] = nz/normal->list_face;
	}


	// TODO verify if the normal is the same for all the nodes in the surface


	//update the vectors
	mesh->local_n_elements = mesh->elements.elem_count;
	mesh->local_n_nodes = mesh->nodes.elem_count;
	MPI_Allreduce(&mesh->local_n_elements, &mesh->total_n_elements, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(&mesh->local_n_nodes, &mesh->total_n_nodes, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);

}
