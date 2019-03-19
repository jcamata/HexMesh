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

	GtsPoint p;
	sc_hash_array_t* pillowIds = sc_hash_array_new(sizeof(pillow_t), vertex_hash_id, vertex_equal_id, &clamped);
	double factor = 0.2;
	double factor_xy = 1e-5;
	//creating the nodes for the pillowing
	for(int ino = 0; ino < nodes_b_mat.size(); ino++){
		size_t position;
		size_t position0;
		pillow_t key;
		key.id=nodes_b_mat[ino];

		//find if the node is in the surface (deal with the coastline)
		octant_node_t* node = (octant_node_t*) sc_array_index (&mesh->nodes, nodes_b_mat[ino]);


		pillow_t* pnode = (pillow_t*) sc_hash_array_insert_unique (pillowIds, &key, &position);
		if(pnode != NULL){
			pnode->id = nodes_b_mat[ino];
			double x = coords[3*nodes_b_mat[ino]+0];
			double y = coords[3*nodes_b_mat[ino]+1];

			if(node->z==0){
				printf("Achei node->z==0, no no:%d\n",nodes_b_mat[ino]);
				double z = coords[3*nodes_b_mat[ino]+2];
				x = coords[3*nodes_b_mat[ino]+0]+factor_xy*coords[3*nodes_b_mat[ino]+0];
				//y=coords[3*nodes_b_mat[ino]+1]+factor*coords[3*nodes_b_mat[ino]+1];
				//z = coords[3*nodes_b_mat[ino]+2] + factor;
				gts_point_set(&p, x, y, z);
				int node = AddPoint(mesh, hash_nodes, &p, coords);
				pnode->b = node;
				x = coords[3*nodes_b_mat[ino]+0]-factor_xy*coords[3*nodes_b_mat[ino]+0];
				//y=coords[3*nodes_b_mat[ino]+1]-factor*coords[3*nodes_b_mat[ino]+1];
				//z = coords[3*nodes_b_mat[ino]+2] - factor;
				gts_point_set(&p, x, y, z);
				node = AddPoint(mesh, hash_nodes, &p, coords);
				pnode->a = node;
				//printf("Node:%d, Above:%d, Bellow:%d\n",pnode->id,pnode->a,pnode->b);
			}else{
				if(coords[3*nodes_b_mat[ino]+2] > 0){
					double z = coords[3*nodes_b_mat[ino]+2]+factor*coords[3*nodes_b_mat[ino]+2];
					//z = coords[3*nodes_b_mat[ino]+2] + factor;
					gts_point_set(&p, x, y, z);
					int node = AddPoint(mesh, hash_nodes, &p, coords);
					pnode->a = node;
					z = coords[3*nodes_b_mat[ino]+2]-factor*coords[3*nodes_b_mat[ino]+2];
					//z = coords[3*nodes_b_mat[ino]+2] - factor;
					gts_point_set(&p, x, y, z);
					node = AddPoint(mesh, hash_nodes, &p, coords);
					pnode->b = node;
					//printf("Node:%d, Above:%d, Bellow:%d\n",pnode->id,pnode->a,pnode->b);
				}else{

				}
			}
		}else{
			printf("Error in add pillowing node in PillowingInterface.cpp, node: %d\n",pnode->id);
		}
	}

	//create the elements for the pillowing...
	for(int ioc = 0 ; ioc < mesh->oct.elem_count; ioc++){
		octree_t * oc = (octree_t*) sc_array_index (&mesh->oct, ioc);
		for(int iel = 0; iel < 8 ; iel++){
			if(oc->id[iel]!=-1){
				octant_t *elem = (octant_t*) sc_array_index(&mesh->elements, oc->id[iel]);

				int connec[8];
				for(int i = 0; i <8; i++){
					connec[i] = elem->nodes[i].id;
				}

				int new_nodes[6][4];
				//initialize new_nodes, TODO can be remove
				for(int i = 0; i <6; i++){
					for( int j = 0; j<4;j++){
						new_nodes[i][j]=-1;
					}
				}

				int face_c=0;
				int face_map[6];
				//loop in the faces
				//check if the surface is locked
				//the new nodes are add to new_nodes
				for(int isurf = 0; isurf < 6; isurf++){
					if(elem->nodes[FaceNodesMap[isurf][0]].fixed==1 && elem->nodes[FaceNodesMap[isurf][1]].fixed==1
							&& elem->nodes[FaceNodesMap[isurf][2]].fixed==1 && elem->nodes[FaceNodesMap[isurf][3]].fixed==1){

						for(int ino = 0; ino<4;ino++){
							size_t position;
							pillow_t key;
							key.id = elem->nodes[FaceNodesMap[isurf][ino]].id;
							//printf("Estou procurando elemento %d, no no:%d, na face:%d\n",elem->id,key.id,isurf);
							bool nfound = sc_hash_array_lookup(pillowIds, &key, &position);
							pillow_t* pil = (pillow_t*) sc_array_index(&pillowIds->a,position);
							//printf("Encontrei o no:%d, junto com os nos a:%d e b:%d\n",pil->id,pil->a,pil->b);
							if(iel<4){
								new_nodes[face_c][ino] = pil->a;
							}else{
								new_nodes[face_c][ino] = pil->b;
							}
						}
						face_map[face_c] = isurf;
						face_c++;
					}

				}

				/*
				printf("Para o elemento:%d\n",elem->id);
				for(int isurf = 0; isurf < face_c;isurf++){
					printf("Face:%d\n os nos sao: ",face_map[isurf]);
					for(int i = 0;i<4;i++){
						printf("%d ", new_nodes[isurf][i]);
					}
					printf("\n");
				}
				 */
				bool create_pelem = true;
				//create new element and do the connectivity of elements
				for(int isurf = 0; isurf < face_c;isurf++){
					//printf("Elemento:%d, face:%d\n",elem->id,face_map[isurf]);
					octant_t* pelem;
					if(create_pelem){
						pelem = (octant_t*) sc_array_push(&mesh->elements);
						pelem->id = mesh->elements.elem_count+1;
						pelem->n_mat = elem->n_mat;
						pelem->x = elem->x;
						pelem->y = elem->y;
						pelem->z = elem->z;
						for(int ino = 0; ino<8; ino++){
							pelem->nodes[ino].fixed= 0;
							pelem->nodes[ino].x = elem->nodes[ino].x;
							pelem->nodes[ino].y = elem->nodes[ino].y;
							pelem->nodes[ino].z = elem->nodes[ino].z;
						}
					}
					for(int ino = 0; ino <4; ino++){
						//printf("No do elemento original:%d\n",connec[FaceNodesMap[face_map[isurf]][ino]]);
						//printf("No do elemento modificado:%d\n",elem->nodes[FaceNodesMap[face_map[isurf]][ino]].id);
						//printf("Trocado por:%d\n",new_nodes[isurf][ino]);

						if(create_pelem){
							pelem->nodes[FaceNodesMap[face_map[isurf]][ino]].id = connec[FaceNodesMap[face_map[isurf]][ino]];
							pelem->nodes[FaceNodesMap[face_map[isurf]][ino]].fixed = elem->nodes[FaceNodesMap[face_map[isurf]][ino]].fixed;
							pelem->nodes[FaceNodesMap_inv[face_map[isurf]][ino]].id= new_nodes[isurf][ino];
						}
						//printf("FaceNodesMap:%d\n",FaceNodesMap[face_map[isurf]][ino]);
						elem->nodes[FaceNodesMap[face_map[isurf]][ino]].id = new_nodes[isurf][ino];
					}

					if(create_pelem){
						for(int ino = 0; ino<8; ino++){
							octant_node_t * node = (octant_node_t*) sc_array_index(&mesh->nodes, pelem->nodes[ino].id);
							node->x=elem->nodes[ino].x;
							node->y=elem->nodes[ino].y;
							node->z=elem->nodes[ino].z;
							node->fixed =pelem->nodes[ino].fixed;
						}
					}

				}
				if(create_pelem){
					//free fiwed nodes
					for(int ino = 0; ino<8; ino++){
						elem->nodes[ino].fixed = 0;
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

}
