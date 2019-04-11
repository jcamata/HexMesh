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

int const VerDiag[8] = {6,7,4,5,2,3,0,1};

void Jac(const double* cord_inx, const double* cord_iny, const double* cord_inz, double J_det[8]){

	//cord_in -> matrix de 3 por 8 que tem as coordenadas (x,y,z) dos 8 nos
	double N[3][8];
	double out[3];
	//double J_det[8];

	int ngp = 2;
	double xr[ngp];
	double yr[ngp];
	double zr[ngp];
	double cord_in[3][8];
	int count = 0;
	xr[0] = -sqrt(3)/3;
	xr[1] = sqrt(3)/3;

	yr[0] = xr[0];
	zr[0] = xr[0];
	yr[1] = xr[1];
	zr[1] = xr[1];

	for(int k = 0; k<8;k++){
		cord_in[0][k] = cord_inx[k];
		cord_in[1][k] = cord_iny[k];
		cord_in[2][k] = cord_inz[k];
	}

	double J[3][3];
	for(int i = 0; i<3; i++){
		for(int j = 0; j<3; j++){
			J[i][j] = 0;
		}
	}

	for(int ix = 0; ix < ngp; ix++){
		for(int iy = 0; iy < ngp; iy++){
			for(int iz = 0; iz < ngp; iz++){

				/*
				N[0] = (1-cord_in_ref[0])*(1-cord_in_ref[1])*(1-cord_in_ref[2])/double(8);
				N[1] = (1+cord_in_ref[0])*(1-cord_in_ref[1])*(1-cord_in_ref[2])/double(8);
				N[2] = (1+cord_in_ref[0])*(1+cord_in_ref[1])*(1-cord_in_ref[2])/double(8);
				N[3] = (1-cord_in_ref[0])*(1+cord_in_ref[1])*(1-cord_in_ref[2])/double(8);

				N[4] = (1-cord_in_ref[0])*(1-cord_in_ref[1])*(1+cord_in_ref[2])/double(8);
				N[5] = (1+cord_in_ref[0])*(1-cord_in_ref[1])*(1+cord_in_ref[2])/double(8);
				N[6] = (1+cord_in_ref[0])*(1+cord_in_ref[1])*(1+cord_in_ref[2])/double(8);
				N[7] = (1-cord_in_ref[0])*(1+cord_in_ref[1])*(1+cord_in_ref[2])/double(8);
				 */

				N[0][0] = -(1 - yr[iy])*(1 - zr[iz])/double(8);
				N[1][0] = -(1 - xr[iy])*(1 - zr[iz])/double(8);
				N[2][0] = -(1 - xr[iy])*(1 - yr[iz])/double(8);

				N[0][1] = +(1 - yr[iy])*(1 - zr[iz])/double(8);
				N[1][1] = -(1 - xr[iy])*(1 - zr[iz])/double(8);
				N[2][1] = -(1 - xr[iy])*(1 - yr[iz])/double(8);

				N[0][2] = +(1 - yr[iy])*(1 - zr[iz])/double(8);
				N[1][2] = +(1 - xr[iy])*(1 - zr[iz])/double(8);
				N[2][2] = -(1 - xr[iy])*(1 - yr[iz])/double(8);

				N[0][3] = -(1 - yr[iy])*(1 - zr[iz])/double(8);
				N[1][3] = +(1 - xr[iy])*(1 - zr[iz])/double(8);
				N[2][3] = -(1 - xr[iy])*(1 - yr[iz])/double(8);

				N[0][4] = -(1 - yr[iy])*(1 - zr[iz])/double(8);
				N[1][4] = -(1 - xr[iy])*(1 - zr[iz])/double(8);
				N[2][4] = +(1 - xr[iy])*(1 - yr[iz])/double(8);

				N[0][5] = +(1 - yr[iy])*(1 - zr[iz])/double(8);
				N[1][5] = -(1 - xr[iy])*(1 - zr[iz])/double(8);
				N[2][5] = +(1 - xr[iy])*(1 - yr[iz])/double(8);

				N[0][6] = +(1 - yr[iy])*(1 - zr[iz])/double(8);
				N[1][6] = +(1 - xr[iy])*(1 - zr[iz])/double(8);
				N[2][6] = +(1 - xr[iy])*(1 - yr[iz])/double(8);

				N[0][7] = -(1 - yr[iy])*(1 - zr[iz])/double(8);
				N[1][7] = +(1 - xr[iy])*(1 - zr[iz])/double(8);
				N[2][7] = +(1 - xr[iy])*(1 - yr[iz])/double(8);

				for(int i = 0; i<3; i++){
					for(int j = 0; j<3; j++){
						for(int k = 0; k<8;k++){
							J[i][j] += N[i][k]*cord_in[j][k];
						}
					}
				}

				double aux = J[1][1]*J[2][2]*J[3][3] + J[1][2]*J[2][3]*J[3][1] + J[1][3]*J[2][1]*J[3][2];
				aux = -J[1][3]*J[2][2]*J[3][1] - J[1][1]*J[2][3]*J[3][2] - J[1][2]*J[2][1]*J[3][3];
				J_det[count] = aux;
				count++;
			}
		}
	}

	//return J_det;

}

void PillowingInterface(hexa_tree_t* mesh, std::vector<double>& coords, std::vector<int>& nodes_b_mat){

	int elem_old = mesh->elements.elem_count;
	int nodes_old = mesh->nodes.elem_count;
	bool clamped = true;

	//criando a hash de nos para evitar nos duplicados no pillowing
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
			r->flag = true;
		}else{
			printf("Verificar o no numero %d\n",node);
		}
	}

	//fazendo vertex hash
	sc_hash_array_t*	vertex_hash = (sc_hash_array_t *)sc_hash_array_new(sizeof (octant_vertex_t), vertex_hash_id, vertex_equal_id, &clamped);
	for (int iel = 0; iel < mesh->elements.elem_count; ++iel) {
		size_t  position;
		octant_t *elem = (octant_t*) sc_array_index(&mesh->elements, iel);

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

	sc_hash_array *hash_normal = sc_hash_array_new(sizeof(normal_t), vertex_hash_id, vertex_equal_id, &clamped);
	//entrando em cada vertice de nodes_b_mat e encontrando o vertice na hash de vertices,
	// acessando a lista de elementos para calculo da normal naquele elemento
	//TODO check if the normal is good with all the surfaces
	for(int ino = 0; ino < nodes_b_mat.size(); ino++){
		octant_vertex_t key;
		size_t position;
		key.id = nodes_b_mat[ino];

		bool lvert = (octant_vertex_t*) sc_hash_array_lookup (vertex_hash, &key, &position);

		if(lvert){
			octant_vertex_t * vert = (octant_vertex_t*) sc_array_index(&vertex_hash->a, position);
			//printf("achei o no numero:%d\n",vert->id);
			for(int iel = 0; iel < vert->list_elem; iel++){
				octant_t *elem = (octant_t*) sc_array_index(&mesh->elements, vert->elem[iel]);

				for(int isurf = 0; isurf < 6; isurf += 2){

					bool lockedSurface = true;
					bool lfvinbmat = true;
					for (int iv = 0 ; iv <4; iv++){
						//procurando se os 4 nos da superficie estao na lista de hash_b_mat
						node_t keyN;
						size_t  positionN;
						int snode = FaceNodesMap[isurf][iv];
						keyN.node_id = elem->nodes[snode].id;
						keyN.coord[0] = coords[3*elem->nodes[snode].id+0];
						keyN.coord[1] = coords[3*elem->nodes[snode].id+1];
						keyN.coord[2] = coords[3*elem->nodes[snode].id+2];
						lfvinbmat = (node_t*) sc_hash_array_lookup (hash_b_mat, &keyN, &positionN);
						if(!lfvinbmat) lockedSurface = false;
					}


					if(lockedSurface){

						int ive;
						//find the node to find the normal
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

						double norm = sqrt(nx*nx+ny*ny+nz*nz);

						normal_t keyN;
						keyN.id = nodes_b_mat[ino];
						normal_t * normal = (normal_t*) sc_hash_array_insert_unique(hash_normal, &keyN, &position);

						if(normal!=NULL){
							normal->id = nodes_b_mat[ino];
							normal->list_elem = 1;
							normal->list_face = 1;
							normal->elem[normal->list_face-1] = elem->id;
							normal->face[normal->list_face-1] = isurf;
							normal->n[normal->list_face-1][0] = nx/norm;
							normal->n[normal->list_face-1][1] = ny/norm;
							normal->n[normal->list_face-1][2] = nz/norm;
						}else{
							normal = (normal_t*) sc_array_index(&hash_normal->a, position);
							normal->elem[normal->list_elem] = elem->id;
							normal->face[normal->list_face] = isurf;
							normal->n[normal->list_face][0] = nx/norm;
							normal->n[normal->list_face][1] = ny/norm;
							normal->n[normal->list_face][2] = nz/norm;
							normal->list_elem++;
							normal->list_face++;
						}
						//printf("node:%d element:%d surface:%d nx:%f ny:%f nz:%f\n",key.id, vert->elem[iel], isurf,nx,ny,nz );
					}

				}
			}
		}else{

		}
	}

	//making the mean value of the surface normal
	for(int ino = 0; ino < hash_normal->a.elem_count; ino ++){
		normal_t* normal = (normal_t*) sc_array_index (&hash_normal->a, ino);
		double nx = 0;
		double ny = 0;
		double nz = 0;

		for(int isurf = 0; isurf < normal->list_face; isurf++){
			nx += (normal->n[isurf][0]);
			ny += (normal->n[isurf][1]);
			nz += (normal->n[isurf][2]);
		}
		normal->nm[0] = nx/normal->list_face;
		normal->nm[1] = ny/normal->list_face;
		normal->nm[2] = nz/normal->list_face;

		double norm_mag = sqrt(normal->nm[0]*normal->nm[0] + normal->nm[1]*normal->nm[1] + normal->nm[2]*normal->nm[2]);

		normal->nm[0] = normal->nm[0]/norm_mag;
		normal->nm[1] = normal->nm[1]/norm_mag;
		normal->nm[2] = normal->nm[2]/norm_mag;

	}

	////////////////////////
	GtsPoint p;
	sc_hash_array_t* pillowIds = sc_hash_array_new(sizeof(pillow_t), vertex_hash_id, vertex_equal_id, &clamped);
	double scale = 10;

	//find the scale for the pillowing
	for(int ioc = 0 ; ioc < mesh->oct.elem_count; ioc++){
		octree_t * oct = (octree_t*) sc_array_index(&mesh->oct, ioc);
		for(int i = 0; i<8; i++){
			if(oct->id[i]!=-1){
				octant_t *elem = (octant_t*) sc_array_index(&mesh->elements, oct->id[i]);

				double dx = abs(coords[3*elem->nodes[0].id+0] - coords[3*elem->nodes[1].id+0]);
				double dy = abs(coords[3*elem->nodes[0].id+1] - coords[3*elem->nodes[3].id+1]);
				double dz = abs(coords[3*elem->nodes[0].id+2] - coords[3*elem->nodes[4].id+2]);

				double dmin = min(dx,dy);
				dmin = min(dmin,dz);
				scale = min(scale,dmin);
			}
		}
	}
	if(scale==0) scale = 0.1;
	scale = 0.5*scale;
	//printf("Escala de: %f\n",scale);
	//creating the nodes for the pillowing
	//TODO fazer as bordas do dominio
	for(int ino = 0; ino < nodes_b_mat.size(); ino++){

		//TODO find if the node is in the surface
		octant_node_t* nodeR = (octant_node_t*) sc_array_index (&mesh->nodes, nodes_b_mat[ino]);

		normal_t keyN;
		size_t  positionN;
		keyN.id = nodes_b_mat[ino];
		bool tre = (normal_t*) sc_hash_array_lookup (hash_normal, &keyN, &positionN);
		normal_t* normal = (normal_t*) sc_array_index (&hash_normal->a, positionN);

		size_t positionP;
		pillow_t keyP;
		keyP.id=nodes_b_mat[ino];
		//printf("%d\n",keyP.id);
		pillow_t* pnode = (pillow_t*) sc_hash_array_insert_unique (pillowIds, &keyP, &positionP);
		double x,y,z;
		if(pnode != NULL){
			//node_t keyM;
			//size_t  positionNM;
			//keyM.node_id = nodes_b_mat[ino];
			//keyM.coord[0] = coords[3*nodes_b_mat[ino]+0];
			//keyM.coord[1] = coords[3*nodes_b_mat[ino]+1];
			//keyM.coord[2] = coords[3*nodes_b_mat[ino]+2];
			//bool nnode = (node_t*) sc_hash_array_lookup (hash_b_mat, &keyN, &positionN);
			//node_t* node = (node_t*) sc_array_index (&hash_b_mat->a, positionN);
			if((nodeR->x == mesh->x_start || nodeR->x == mesh->x_end) &&
					(nodeR->y == mesh->y_start || nodeR->y == mesh->y_end)){
				pnode->id = nodes_b_mat[ino];
				double xr = coords[3*nodes_b_mat[ino]+0];
				double yr = coords[3*nodes_b_mat[ino]+1];
				double zr = coords[3*nodes_b_mat[ino]+2];

				x = xr ;
				y = yr ;
				z = zr + scale*normal->nm[2];

				gts_point_set(&p, x, y, z);
				int node = AddPoint(mesh, hash_nodes, &p, coords);
				pnode->a = node;

				x = xr ;
				y = yr ;
				z = zr - scale*normal->nm[2];

				gts_point_set(&p, x, y, z);
				node = AddPoint(mesh, hash_nodes, &p, coords);
				pnode->b = node;
			}else if((nodeR->x == mesh->x_start || nodeR->x == mesh->x_end) &&
					(nodeR->z == 0 || nodeR->z == mesh->max_z)){
				pnode->id = nodes_b_mat[ino];
				double xr = coords[3*nodes_b_mat[ino]+0];
				double yr = coords[3*nodes_b_mat[ino]+1];
				double zr = coords[3*nodes_b_mat[ino]+2];

				x = xr ;
				y = yr + scale*normal->nm[1];
				z = zr ;

				gts_point_set(&p, x, y, z);
				int node = AddPoint(mesh, hash_nodes, &p, coords);
				pnode->a = node;

				x = xr ;
				y = yr - scale*normal->nm[1];
				z = zr ;

				gts_point_set(&p, x, y, z);
				node = AddPoint(mesh, hash_nodes, &p, coords);
				pnode->b = node;
			}else if((nodeR->y == mesh->y_start || nodeR->y == mesh->y_end) &&
					(nodeR->z == 0 || nodeR->z == mesh->max_z)){
				pnode->id = nodes_b_mat[ino];
				double xr = coords[3*nodes_b_mat[ino]+0];
				double yr = coords[3*nodes_b_mat[ino]+1];
				double zr = coords[3*nodes_b_mat[ino]+2];

				x = xr + scale*normal->nm[0];
				y = yr ;
				z = zr ;

				gts_point_set(&p, x, y, z);
				int node = AddPoint(mesh, hash_nodes, &p, coords);
				pnode->a = node;

				x = xr - scale*normal->nm[0];
				y = yr ;
				z = zr ;

				gts_point_set(&p, x, y, z);
				node = AddPoint(mesh, hash_nodes, &p, coords);
				pnode->b = node;
			}else if(nodeR->x == mesh->x_start || nodeR->x == mesh->x_end){
				pnode->id = nodes_b_mat[ino];
				double xr = coords[3*nodes_b_mat[ino]+0];
				double yr = coords[3*nodes_b_mat[ino]+1];
				double zr = coords[3*nodes_b_mat[ino]+2];

				x = xr ;
				y = yr + scale*normal->nm[1];
				z = zr + scale*normal->nm[2];

				gts_point_set(&p, x, y, z);
				int node = AddPoint(mesh, hash_nodes, &p, coords);
				pnode->a = node;

				x = xr ;
				y = yr - scale*normal->nm[1];
				z = zr - scale*normal->nm[2];

				gts_point_set(&p, x, y, z);
				node = AddPoint(mesh, hash_nodes, &p, coords);
				pnode->b = node;
			}else if(nodeR->y == mesh->y_start || nodeR->y == mesh->y_end){
				pnode->id = nodes_b_mat[ino];
				double xr = coords[3*nodes_b_mat[ino]+0];
				double yr = coords[3*nodes_b_mat[ino]+1];
				double zr = coords[3*nodes_b_mat[ino]+2];

				x = xr + scale*normal->nm[0];
				y = yr ;
				z = zr + scale*normal->nm[2];

				gts_point_set(&p, x, y, z);
				int node = AddPoint(mesh, hash_nodes, &p, coords);
				pnode->a = node;

				x = xr - scale*normal->nm[0];
				y = yr ;
				z = zr - scale*normal->nm[2];

				gts_point_set(&p, x, y, z);
				node = AddPoint(mesh, hash_nodes, &p, coords);
				pnode->b = node;
			}else if(nodeR->z == 0 || nodeR->z == mesh->max_z){
				pnode->id = nodes_b_mat[ino];
				double xr = coords[3*nodes_b_mat[ino]+0];
				double yr = coords[3*nodes_b_mat[ino]+1];
				double zr = coords[3*nodes_b_mat[ino]+2];

				x = xr + scale*normal->nm[0];
				y = yr + scale*normal->nm[1];
				z = zr ;

				gts_point_set(&p, x, y, z);
				int node = AddPoint(mesh, hash_nodes, &p, coords);
				pnode->a = node;

				x = xr - scale*normal->nm[0];
				y = yr - scale*normal->nm[1];
				z = zr ;

				gts_point_set(&p, x, y, z);
				node = AddPoint(mesh, hash_nodes, &p, coords);
				pnode->b = node;
			}else{
				pnode->id = nodes_b_mat[ino];
				double xr = coords[3*nodes_b_mat[ino]+0];
				double yr = coords[3*nodes_b_mat[ino]+1];
				double zr = coords[3*nodes_b_mat[ino]+2];

				x = xr + scale*normal->nm[0];
				y = yr + scale*normal->nm[1];
				z = zr + scale*normal->nm[2];

				gts_point_set(&p, x, y, z);
				int node = AddPoint(mesh, hash_nodes, &p, coords);
				pnode->a = node;

				x = xr - scale*normal->nm[0];
				y = yr - scale*normal->nm[1];
				z = zr - scale*normal->nm[2];

				gts_point_set(&p, x, y, z);
				node = AddPoint(mesh, hash_nodes, &p, coords);
				pnode->b = node;
			}
		}else{
			printf("Error in add pillowing node in PillowingInterface.cpp, node: %d\n",keyP.id);
		}
	}

	//create the elements for the pillowing...
	for(int ioc = 0 ; ioc < mesh->oct.elem_count; ioc++){
		octree_t * oc = (octree_t*) sc_array_index (&mesh->oct, ioc);
		for(int iel = 0; iel < 8 ; iel++){
			if(oc->id[iel]!=-1){
				octant_t *elem = (octant_t*) sc_array_index(&mesh->elements, oc->id[iel]);

				int connec[8];
				for(int ino = 0; ino <8; ino++){
					connec[ino] = elem->nodes[ino].id;
				}

				int face_c=0;
				int face_map[6];
				//loop in the faces
				//check if the surface is locked
				//then new nodes are add to new_nodes
				for(int isurf = 0; isurf < 6; isurf++){
					bool lockedSurface = true;
					bool lfvinbmat = true;
					for (int iv = 0 ; iv <4; iv++){
						//procurando se os 4 nos da superficie estao na lista de hash_b_mat
						node_t keyN;
						size_t  positionN;
						int snode = FaceNodesMap[isurf][iv];
						keyN.node_id = elem->nodes[snode].id;
						keyN.coord[0] = coords[3*elem->nodes[snode].id+0];
						keyN.coord[1] = coords[3*elem->nodes[snode].id+1];
						keyN.coord[2] = coords[3*elem->nodes[snode].id+2];
						lfvinbmat = (node_t*) sc_hash_array_lookup (hash_b_mat, &keyN, &positionN);
						if(!lfvinbmat) lockedSurface = false;
					}

					if(lockedSurface){
						/*
						for(int ino = 0; ino<4;ino++){
							size_t position;
							pillow_t key;
							key.id = elem->nodes[FaceNodesMap[isurf][ino]].id;
							//printf("Estou procurando elemento %d, no no:%d, na face:%d\n",elem->id,key.id,isurf);
							bool nfound = sc_hash_array_lookup(pillowIds, &key, &position);
							pillow_t* pil = (pillow_t*) sc_array_index(&pillowIds->a,position);
							//printf("Encontrei o no:%d, junto com os nos a:%d e b:%d\n",pil->id,pil->a,pil->b);
						}
						 */
						face_map[face_c] = isurf;
						face_c++;
					}else{

					}
				}

				int new_nodes[4];
				bool create_pelem = true;
				//create new element and do the connectivity of elements
				for(int isurf = 0; isurf < face_c;isurf++){

					// select the new nodes created in the pillowing process
					// for the surface in the element
					for(int ino = 0; ino <4; ino++){
						//check if the node is in the hash_b_mat
						node_t keyN;
						size_t  positionN;
						int surfnode = FaceNodesMap[face_map[isurf]][ino];
						keyN.node_id = elem->nodes[surfnode].id;
						keyN.coord[0] = coords[3*elem->nodes[surfnode].id+0];
						keyN.coord[1] = coords[3*elem->nodes[surfnode].id+1];
						keyN.coord[2] = coords[3*elem->nodes[surfnode].id+2];
						bool snode = (node_t*) sc_hash_array_lookup (hash_b_mat, &keyN, &positionN);

						if(snode){
							size_t position;
							pillow_t key;
							key.id = elem->nodes[surfnode].id;
							//printf("Estou procurando elemento %d, no no:%d, na face:%d\n",elem->id,key.id,isurf);
							bool nfound = sc_hash_array_lookup(pillowIds, &key, &position);
							pillow_t* pil = (pillow_t*) sc_array_index(&pillowIds->a,position);

							GtsPoint p0;
							GtsPoint pa;
							GtsPoint pb;

							gts_point_set(&pa, coords[3*pil->a+0], coords[3*pil->a+1], coords[3*pil->a+2]);
							gts_point_set(&pb, coords[3*pil->b+0], coords[3*pil->b+1], coords[3*pil->b+2]);

							int nnode = VerDiag[surfnode];
							p0.x = coords[3*elem->nodes[nnode].id+0];
							p0.y = coords[3*elem->nodes[nnode].id+1];
							p0.z = coords[3*elem->nodes[nnode].id+2];

							double dist0a[3];
							double dist0b[3];
							/*
							for(int ive = 0 ; ive < 3; ive++){
								int nnode = VertexVertexMap[surfnode][ive];
								p0.x = coords[3*elem->nodes[nnode].id+0];
								p0.y = coords[3*elem->nodes[nnode].id+1];
								p0.z = coords[3*elem->nodes[nnode].id+2];
								dist0a[ive] = gts_point_distance (&p0,&pa);
								dist0b[ive] = gts_point_distance (&p0,&pb);
							}

							double dista = sqrt(dist0a[0]*dist0a[0]+dist0a[1]*dist0a[1]+dist0a[2]*dist0a[2]);
							double distb = sqrt(dist0b[0]*dist0b[0]+dist0b[1]*dist0b[1]+dist0b[2]*dist0b[2]);
							 */
							dist0a[0] = gts_point_distance (&p0,&pa);
							dist0b[0] = gts_point_distance (&p0,&pb);
							double dista = dist0a[0];
							double distb = dist0b[0];
							if(dista<distb){
								new_nodes[ino] = pil->a;
							}else{
								new_nodes[ino] = pil->b;
							}



						}
					}

					//printf("Elemento:%d, face:%d, new_nodes: %d %d %d %d\n",elem->id,face_map[isurf],new_nodes[0],new_nodes[1],new_nodes[2],new_nodes[3]);

					//printf("Elemento:%d, face:%d\n",elem->id,face_map[isurf]);
					//create the element
					octant_t* pelem;
					if(create_pelem){
						pelem = (octant_t*) sc_array_push(&mesh->elements);
						pelem->id = mesh->elements.elem_count+1;
						pelem->n_mat = elem->n_mat;
						pelem->pad = elem->pad;
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
						//printf("Trocado por:%d\n",new_nodes[ino]);

						if(create_pelem){
							pelem->nodes[FaceNodesMap[face_map[isurf]][ino]].id = connec[FaceNodesMap[face_map[isurf]][ino]];
							pelem->nodes[FaceNodesMap[face_map[isurf]][ino]].fixed = elem->nodes[FaceNodesMap[face_map[isurf]][ino]].fixed;
							pelem->nodes[FaceNodesMap_inv[face_map[isurf]][ino]].id= new_nodes[ino];
							//printf("FaceNodesMap:%d\n",FaceNodesMap[face_map[isurf]][ino]);
							elem->nodes[FaceNodesMap[face_map[isurf]][ino]].id = new_nodes[ino];
						}

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

					if(false){
						printf("El:%d conectividade orignal\n",elem->id);
						for(int ino = 0; ino<8; ino++){
							printf(" %d",connec[ino]);
						}
						printf("\nModificada \n");
						for(int ino = 0; ino<8; ino++){
							printf("%d ",elem->nodes[ino].id);
						}
						printf("\n");
						printf("New El:%d conectividade:\n",pelem->id);
						for(int ino = 0; ino<8; ino++){
							printf("%d ",pelem->nodes[ino].id);
						}
						printf("\n");
					}
				}

				if(face_c==0 && true){
					//printf("Para o elemento %d\n",elem->id);
					int edge_c=0;
					int edge_map[12];
					for(int iedge = 0; iedge < 12; iedge++){
						//procurando se os 2 nos da aresta estao na lista de hash_b_mat
						node_t keyN;
						size_t  positionN;
						bool n0 = false;
						bool n1 = false;
						//printf("edge: %d\n",iedge);
						int enode = EdgeVerticesMap[iedge][0];
						//printf("nos: %d",elem->nodes[enode].id);
						keyN.node_id = elem->nodes[enode].id;
						keyN.coord[0] = coords[3*elem->nodes[enode].id+0];
						keyN.coord[1] = coords[3*elem->nodes[enode].id+1];
						keyN.coord[2] = coords[3*elem->nodes[enode].id+2];
						n0 = (node_t*) sc_hash_array_lookup (hash_b_mat, &keyN, &positionN);

						enode = EdgeVerticesMap[iedge][1];
						//printf(" %d\n",elem->nodes[enode].id);
						keyN.node_id = elem->nodes[enode].id;
						keyN.coord[0] = coords[3*elem->nodes[enode].id+0];
						keyN.coord[1] = coords[3*elem->nodes[enode].id+1];
						keyN.coord[2] = coords[3*elem->nodes[enode].id+2];
						n1 = (node_t*) sc_hash_array_lookup (hash_b_mat, &keyN, &positionN);

						if(n0 && n1){
							edge_map[edge_c] = iedge;
							edge_c++;
						}
					}

					//printf("Para o elemento %d tivemos %d arestas\n",elem->id,edge_c);

					if(edge_c > 0){
						//printf("Elemento:%d\n",elem->id);

						for(int iedge = 0; iedge < edge_c; iedge++){
							for(int ino = 0; ino <2;ino++){
								int enode = EdgeVerticesMap[edge_map[iedge]][ino];
								//printf("No:%d da aresta:%d\n",enode,edge_map[iedge]);
								size_t position;
								pillow_t key;
								key.id = elem->nodes[enode].id;
								//printf("Estou procurando elemento %d, no no:%d, na face:%d\n",elem->id,key.id,isurf);
								bool nfound = sc_hash_array_lookup(pillowIds, &key, &position);
								pillow_t* pil = (pillow_t*) sc_array_index(&pillowIds->a,position);

								GtsPoint p0;
								GtsPoint pa;
								GtsPoint pb;

								gts_point_set(&pa, coords[3*pil->a+0], coords[3*pil->a+1], coords[3*pil->a+2]);
								gts_point_set(&pb, coords[3*pil->b+0], coords[3*pil->b+1], coords[3*pil->b+2]);

								int nnode = VerDiag[enode];
								p0.x = coords[3*elem->nodes[nnode].id+0];
								p0.y = coords[3*elem->nodes[nnode].id+1];
								p0.z = coords[3*elem->nodes[nnode].id+2];

								double dist0a[3];
								double dist0b[3];

								dist0a[0] = gts_point_distance (&p0,&pa);
								dist0b[0] = gts_point_distance (&p0,&pb);
								double dista = dist0a[0];
								double distb = dist0b[0];
								if(dista<distb){
									new_nodes[ino] = pil->a;
								}else{
									new_nodes[ino] = pil->b;
								}
							}
						}

						for(int iedge = 0; iedge < edge_c; iedge++){
							for(int ino = 0; ino <2; ino++){
								int enode = EdgeVerticesMap[edge_map[iedge]][ino];
								//printf("elemento:%d aresta:%d No:%d\n",elem->id,elem->edge[edge_map[iedge]],enode);
								elem->nodes[enode].id = new_nodes[ino];
							}

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


	// for debug
	if(false){
		char fdname[80];
		sprintf(fdname,"normal_%04d.txt", mesh->mpi_rank);
		FILE* treta = fopen(fdname,"w");

		sprintf(fdname,"below_%04d.txt", mesh->mpi_rank);
		FILE* treta1 = fopen(fdname,"w");

		sprintf(fdname,"above_%04d.txt", mesh->mpi_rank);
		FILE* treta2 = fopen(fdname,"w");

		sprintf(fdname,"ids_%04d.txt", mesh->mpi_rank);
		FILE* treta3 = fopen(fdname,"w");

		sprintf(fdname,"nodes_%04d.txt", mesh->mpi_rank);
		FILE* treta4 = fopen(fdname,"w");

		for(int ino = 0; ino < hash_normal->a.elem_count; ino ++){
			normal_t* normal = (normal_t*) sc_array_index (&hash_normal->a, ino);
			int node = 3*normal->id;
			fprintf(treta,"%f %f %f %f %f %f\n",coords[node+0],normal->nm[0],coords[node+1],normal->nm[1],coords[node+2],normal->nm[2]);
		}

		for(int ino = 0; ino < pillowIds->a.elem_count; ino ++){
			pillow_t* pil = (pillow_t*) sc_array_index (&pillowIds->a, ino);
			fprintf(treta3,"%d %d %d\n",pil->id,pil->a,pil->b);

			int node = 3*pil->b;
			fprintf(treta1,"%f %f %f\n",coords[node+0],coords[node+1],coords[node+2]);
			node = 3*pil->a;
			fprintf(treta2,"%f %f %f\n",coords[node+0],coords[node+1],coords[node+2]);
			node = 3*pil->id;
			fprintf(treta4,"%f %f %f\n",coords[node+0],coords[node+1],coords[node+2]);
		}

		fclose(treta);
		fclose(treta1);
		fclose(treta2);
		fclose(treta3);
		fclose(treta4);

		for(int iel = 0 ; iel < mesh->elements.elem_count; iel++){
			octant_t *elem = (octant_t*) sc_array_index(&mesh->elements, iel);
			printf("\n \n");
			printf("Conectividade no final do elemento:%d  eh", elem->id);
			for(int ino = 0; ino<8; ino++){
				printf("%d ",elem->nodes[ino].id);
			}
			printf("\n \n");
		}

	}




	//update the vectors
	mesh->local_n_elements = mesh->elements.elem_count;
	mesh->local_n_nodes = mesh->nodes.elem_count;
	MPI_Allreduce(&mesh->local_n_elements, &mesh->total_n_elements, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(&mesh->local_n_nodes, &mesh->total_n_nodes, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
	printf("Criei %d elementos\n",mesh->elements.elem_count - elem_old);
	printf("Criei %d nos\n",mesh->nodes.elem_count - nodes_old);
}

/*
 *
 * //check if the element is good
					if(false){
						double cord_in_x[8], cord_in_y[8], cord_in_z[8];
						//add the nodes in the coord vector
						for(int ino = 0; ino<8;ino++){
							cord_in_x[ino] = coords[3 * elem->nodes[ino].id + 0];
							cord_in_y[ino] = coords[3 * elem->nodes[ino].id + 1];
							cord_in_z[ino] = coords[3 * elem->nodes[ino].id + 2];
						}

						double test[8];
						double J_det[8];
						Jac(cord_in_x,cord_in_y,cord_in_z, J_det);

						for(int i = 0; i<8 ; i++){
							printf("Elemento:%d, o jacobiano nos pontos de Gauss:%f\n",pelem->id,J_det[i]);
						}
					}

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

	////////////////////////

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

	printf("norma:%f a hash tem:%d elementos\n",norm_mag,hash_normal->a.elem_count);
	double nx = 0;
	double ny = 0;
	double nz = 0;
	//making the mean value of the surface normal
	//TODO check if the normal is good with all the surfaces
	for(int ino = 0; ino < hash_normal->a.elem_count; ino ++){
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

	// for debug
	if(true){
		char fdname[80];
		sprintf(fdname,"normal_%04d.txt", mesh->mpi_rank);
		FILE* treta = fopen(fdname,"w");
		for(int ino = 0; ino < hash_normal->a.elem_count; ino ++){
			normal_t* normal = (normal_t*) sc_array_index (&hash_normal->a, ino);
			int node = 3*normal->id;
			fprintf(treta,"%f %f %f %f %f %f\n",coords[node+0],normal->nm[0],coords[node+1],normal->nm[1],coords[node+2],normal->nm[2]);
		}
		fclose(treta);

	}



	////////////////////////

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
				y=coords[3*nodes_b_mat[ino]+1]+factor_xy*coords[3*nodes_b_mat[ino]+1];
				//z = coords[3*nodes_b_mat[ino]+2] + factor;
				gts_point_set(&p, x, y, z);
				int node = AddPoint(mesh, hash_nodes, &p, coords);
				pnode->b = node;
				x = coords[3*nodes_b_mat[ino]+0]-factor_xy*coords[3*nodes_b_mat[ino]+0];
				y=coords[3*nodes_b_mat[ino]+1]-factor_xy*coords[3*nodes_b_mat[ino]+1];
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
					double z = coords[3*nodes_b_mat[ino]+2]-factor*coords[3*nodes_b_mat[ino]+2];
					//z = coords[3*nodes_b_mat[ino]+2] + factor;
					gts_point_set(&p, x, y, z);
					int node = AddPoint(mesh, hash_nodes, &p, coords);
					pnode->a = node;
					z = coords[3*nodes_b_mat[ino]+2]+factor*coords[3*nodes_b_mat[ino]+2];
					//z = coords[3*nodes_b_mat[ino]+2] - factor;
					gts_point_set(&p, x, y, z);
					node = AddPoint(mesh, hash_nodes, &p, coords);
					pnode->b = node;
					//printf("Node:%d, Above:%d, Bellow:%d\n",pnode->id,pnode->a,pnode->b);
				}
			}
		}else{
			printf("Error in add pillowing node in PillowingInterface.cpp, node: %d\n",pnode->id);
		}
	}

	//create the elements for the pillowing...
	for(int ioc = 0 ; ioc < mesh->oct.elem_count; ioc++){
		octree_t * oc = (octree_t*) sc_array_index (&mesh->oct, ioc);
		for(int iel = 0; iel < 1 ; iel++){
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
/*
				bool create_pelem = false;
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
							printf("FaceNodesMap:%d\n",FaceNodesMap[face_map[isurf]][ino]);
							elem->nodes[FaceNodesMap[face_map[isurf]][ino]].id = new_nodes[isurf][ino];
						}

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

 */

/*
	//criando a hash de nos para evitar nos duplicados no pillowing
	sc_hash_array_t* hash_nodes = sc_hash_array_new(sizeof(node_t), edge_hash_fn, edge_equal_fn, &clamped);
	for(int n = 0; n < mesh->nodes.elem_count; n++){
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

	//fazendo vertex hash
	sc_hash_array_t*	vertex_hash = (sc_hash_array_t *)sc_hash_array_new(sizeof (octant_vertex_t), vertex_hash_id, vertex_equal_id, &clamped);
	for (int iel = 0; iel < mesh->elements.elem_count; ++iel) {
		size_t  position;
		octant_t *elem = (octant_t*) sc_array_index(&mesh->elements, iel);

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

	sc_hash_array *hash_normal = sc_hash_array_new(sizeof(normal_t), vertex_hash_id, vertex_equal_id, &clamped);
	//entrando em cada vertice de nodes_b_mat e encontrando o vertice na hash de vertices,
	// acessando a lista de elementos para calculo da normal naquele elemento
	for(int ino = 0; ino < nodes_b_mat.size(); ino++){

		octant_vertex_t key;
		size_t  position;
		key.id = nodes_b_mat[ino];
		//printf("%d\n",key.id);
		bool tre = (octant_vertex_t*) sc_hash_array_lookup (vertex_hash, &key, &position);

		if(tre){
			octant_vertex_t * vert = (octant_vertex_t*) sc_array_index(&vertex_hash->a, position);
			//printf("achei o no numero:%d\n",vert->id);
			for(int iel = 0; iel < vert->list_elem; iel++){
				octant_t *elem = (octant_t*) sc_array_index(&mesh->elements, vert->elem[iel]);

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

						normal_t keyN;
						keyN.id = nodes_b_mat[ino];
						normal_t * normal = (normal_t*) sc_hash_array_insert_unique(hash_normal, &keyN, &position);

						if(normal!=NULL){
							normal->id = nodes_b_mat[ino];
							normal->list_elem = 1;
							normal->list_face = 1;
							normal->elem[normal->list_face-1] = elem->id;
							normal->face[normal->list_face-1] = isurf;
							normal->n[normal->list_face-1][0] = nx;
							normal->n[normal->list_face-1][1] = ny;
							normal->n[normal->list_face-1][2] = nz;
						}else{
							normal = (normal_t*) sc_array_index(&hash_normal->a, position);
							normal->elem[normal->list_elem] = elem->id;
							normal->face[normal->list_face] = isurf;
							normal->n[normal->list_face][0] = nx;
							normal->n[normal->list_face][1] = ny;
							normal->n[normal->list_face][2] = nz;
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

	double nx = 0;
	double ny = 0;
	double nz = 0;

	//making the mean value of the surface normal
	//TODO check if the normal is good with all the surfaces
	for(int ino = 0; ino < hash_normal->a.elem_count; ino ++){
		normal_t* normal = (normal_t*) sc_array_index (&hash_normal->a, ino);
		for(int isurf = 0; isurf < normal->list_face; isurf++){
			nx += normal->n[isurf][0];
			ny += normal->n[isurf][1];
			nz += normal->n[isurf][2];
		}
		normal->nm[0] = nx/normal->list_face;
		normal->nm[1] = ny/normal->list_face;
		normal->nm[2] = nz/normal->list_face;

		double norm_mag = sqrt(normal->nm[0]*normal->nm[0]+normal->nm[1]*normal->nm[1]+
				normal->nm[2]*normal->nm[2]);

		normal->nm[0] = normal->nm[0]/norm_mag;
		normal->nm[1] = normal->nm[1]/norm_mag;
		normal->nm[2] = normal->nm[2]/norm_mag;
	}

	////////////////////////
	GtsPoint p;
	sc_hash_array_t* pillowIds = sc_hash_array_new(sizeof(pillow_t), vertex_hash_id, vertex_equal_id, &clamped);
	double scale = 10;

	//find the scale
	for(int ioc = 0 ; ioc < mesh->oct.elem_count; ioc++){
		octree_t * oct = (octree_t*) sc_array_index(&mesh->oct, ioc);
		for(int i = 0; i<8; i++){
			if(oct->id[i]!=-1){
				octant_t *elem = (octant_t*) sc_array_index(&mesh->elements, oct->id[i]);

				double dx = abs(coords[3*elem->nodes[0].id+0] - coords[3*elem->nodes[1].id+0]);
				double dy = abs(coords[3*elem->nodes[0].id+1] - coords[3*elem->nodes[3].id+1]);
				double dz = abs(coords[3*elem->nodes[0].id+2] - coords[3*elem->nodes[4].id+2]);

				double dmin = min(dx,dy);
				dmin = min(dmin,dz);
				scale = min(scale,dmin);
			}
		}
	}
	if(scale==0) scale = 0.1;
	scale = 0.9*scale;
	printf("Escala de: %f\n",scale);
	//creating the nodes for the pillowing
	//TODO fazer as bordas do dominio
	for(int ino = 0; ino < nodes_b_mat.size(); ino++){

		//TODO find if the node is in the surface (deal with the coastline)
		//octant_node_t* node = (octant_node_t*) sc_array_index (&mesh->nodes, nodes_b_mat[ino]);

		normal_t keyN;
		size_t  positionN;
		bool tre = (normal_t*) sc_hash_array_lookup (hash_normal, &keyN, &positionN);
		normal_t* normal = (normal_t*) sc_array_index (&hash_normal->a, positionN);

		size_t positionP;
		pillow_t keyP;
		keyP.id=nodes_b_mat[ino];
		//printf("%d\n",keyP.id);
		pillow_t* pnode = (pillow_t*) sc_hash_array_insert_unique (pillowIds, &keyP, &positionP);
		double x,y,z;
		if(pnode != NULL){
			pnode->id = nodes_b_mat[ino];
			double xr = coords[3*nodes_b_mat[ino]+0];
			double yr = coords[3*nodes_b_mat[ino]+1];
			double zr = coords[3*nodes_b_mat[ino]+2];

			x = xr + scale*normal->nm[0];
			y = yr + scale*normal->nm[1];
			z = zr + scale*normal->nm[2];

			gts_point_set(&p, x, y, z);
			int node = AddPoint(mesh, hash_nodes, &p, coords);
			pnode->a = node;

			x = xr - scale*normal->nm[0];
			y = yr - scale*normal->nm[1];
			z = zr - scale*normal->nm[2];

			gts_point_set(&p, x, y, z);
			node = AddPoint(mesh, hash_nodes, &p, coords);
			pnode->b = node;

		}else{
			printf("Error in add pillowing node in PillowingInterface.cpp, node: %d\n",pnode->id);
		}
	}

	// for debug
	if(true){
		char fdname[80];
		sprintf(fdname,"normal_%04d.txt", mesh->mpi_rank);
		FILE* treta = fopen(fdname,"w");

		sprintf(fdname,"below_%04d.txt", mesh->mpi_rank);
		FILE* treta1 = fopen(fdname,"w");

		sprintf(fdname,"above_%04d.txt", mesh->mpi_rank);
		FILE* treta2 = fopen(fdname,"w");

		sprintf(fdname,"ids_%04d.txt", mesh->mpi_rank);
		FILE* treta3 = fopen(fdname,"w");

		sprintf(fdname,"nodes_%04d.txt", mesh->mpi_rank);
		FILE* treta4 = fopen(fdname,"w");

		for(int ino = 0; ino < hash_normal->a.elem_count; ino ++){
			normal_t* normal = (normal_t*) sc_array_index (&hash_normal->a, ino);
			int node = 3*normal->id;
			fprintf(treta,"%f %f %f %f %f %f\n",coords[node+0],normal->nm[0],coords[node+1],normal->nm[1],coords[node+2],normal->nm[2]);
		}

		for(int ino = 0; ino < pillowIds->a.elem_count; ino ++){
			pillow_t* pil = (pillow_t*) sc_array_index (&pillowIds->a, ino);
			fprintf(treta3,"%d %d %d\n",pil->id,pil->a,pil->b);

			int node = 3*pil->b;
			fprintf(treta1,"%f %f %f\n",coords[node+0],coords[node+1],coords[node+2]);
			node = 3*pil->a;
			fprintf(treta2,"%f %f %f\n",coords[node+0],coords[node+1],coords[node+2]);
			node = 3*pil->id;
			fprintf(treta4,"%f %f %f\n",coords[node+0],coords[node+1],coords[node+2]);
		}

		fclose(treta);
		fclose(treta1);
		fclose(treta2);
		fclose(treta3);
		fclose(treta4);
	}

	//mudando a conectividade dos elementos que ja existem
	for(int ioc = 0 ; ioc < mesh->oct.elem_count; ioc++){
		octree_t * oc = (octree_t*) sc_array_index (&mesh->oct, ioc);
		for(int iel = 0; iel < 8 ; iel++){
			if(oc->id[iel]!=-1){
				octant_t *elem = (octant_t*) sc_array_index(&mesh->elements, oc->id[iel]);

				for(int ino = 0; ino < 8; ino++){

					node_t keyN;
					size_t  positionN;
					keyN.node_id = elem->nodes[ino].id;
					keyN.coord[0] = coords[3*elem->nodes[ino].id+0];
					keyN.coord[1] = coords[3*elem->nodes[ino].id+1];
					keyN.coord[2] = coords[3*elem->nodes[ino].id+2];
					bool nnode = (node_t*) sc_hash_array_lookup (hash_b_mat, &keyN, &positionN);

					if(nnode){
						int Nnode = elem->nodes[diagonal_vert[ino]].id;
						int Rnode = elem->nodes[ino].id;
						//printf("Node:%d\n",Rnode);

						double xr = coords[3*Nnode + 0];
						double yr = coords[3*Nnode + 1];
						double zr = coords[3*Nnode + 2];

						size_t positionP;
						pillow_t keyP;
						keyP.id = Rnode;

						bool pil = (pillow_t*) sc_hash_array_lookup (pillowIds, &keyP, &positionP);

						if(pil){
							pillow_t* Npil = (pillow_t*) sc_array_index (&pillowIds->a, positionP);

							double xa = coords[3*Npil->a + 0];
							double ya = coords[3*Npil->a + 1];
							double za = coords[3*Npil->a + 2];

							double xb = coords[3*Npil->b + 0];
							double yb = coords[3*Npil->b + 1];
							double zb = coords[3*Npil->b + 2];

							double dist1 = sqrt((xr-xa)*(xr-xa)+(yr-ya)*(yr-ya)+(zr-za)*(zr-za));
							double dist2 = sqrt((xr-xb)*(xr-xb)+(yr-yb)*(yr-yb)+(zr-zb)*(zr-zb));

							if(dist2 > dist1){
								elem->nodes[ino].id = Npil->a;
								printf("Elemento:%d, no no:%d e sera subistituido por %d\n",elem->id,Rnode,Npil->a);
							}else{
								elem->nodes[ino].id = Npil->b;
								printf("Elemento:%d, no no:%d e sera subistituido por %d\n",elem->id,Rnode,Npil->b);
							}
						}else{

						}
					}else{

					}
				}
			}
		}
	}


	//create the elements for the pillowing...
	//correct the nodes in the elements
	bool create_pelem = false;
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
				//initialize new_nodes
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


	for(int ioc = 0 ; ioc < 0; ioc++){
		//for(int ioc = 0 ; ioc < mesh->oct.elem_count; ioc++){
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


				printf("Para o elemento:%d\n",elem->id);
				for(int isurf = 0; isurf < face_c;isurf++){
					printf("Face:%d\n os nos sao: ",face_map[isurf]);
					for(int i = 0;i<4;i++){
						printf("%d ", new_nodes[isurf][i]);
					}
					printf("\n");
				}


				bool create_pelem = true;
				//create new element and do the connectivity of elements
				for(int isurf = 0; isurf < face_c;isurf++){

					//printf("Elemento:%d, face:%d\n",elem->id,face_map[isurf]);
					octant_t* pelem;
					if(create_pelem){
						pelem = (octant_t*) sc_array_push(&mesh->elements);
						pelem->id = mesh->elements.elem_count+1;
						pelem->n_mat = elem->n_mat;
						pelem->pad = elem->pad;
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
							//printf("FaceNodesMap:%d\n",FaceNodesMap[face_map[isurf]][ino]);
							elem->nodes[FaceNodesMap[face_map[isurf]][ino]].id = new_nodes[isurf][ino];
						}

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
 */
