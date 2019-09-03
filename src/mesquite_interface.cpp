#include <sstream>
#include <vector>
#include <mpi.h>
#include <assert.h>
#include <algorithm>
#include <list>
#include <iostream>
#include <numeric>
#include <random>

//all headers
#include "Mesquite_all_headers.hpp"

#include "hexa.h"

using namespace Mesquite;

void OptLine(hexa_tree_t* mesh, std::vector<double>& coords, sc_hash_array_t* hash_FixedNodes){


	bool deb = true;
	if(deb){
		for(int ino = 0; ino < mesh->nodes.elem_count; ino++) mesh->part_nodes[ino] = -1;
	}

	std::vector<int> aux0;
	std::vector<int> aux1;
	std::vector<int> aux2;
	std::vector<int> aux3;

	std::vector<int> aux4;
	std::vector<int> aux5;
	std::vector<int> aux6;
	std::vector<int> aux7;

	std::vector<int> aux8;
	std::vector<int> aux9;
	std::vector<int> aux10;
	std::vector<int> aux11;

	for(int iel = 0; iel < mesh->outsurf.elem_count; iel++){
		octant_t* elem = (octant_t*) sc_array_index (&mesh->outsurf, iel);

		for(int iedge = 0; iedge < 12; iedge++){
			if(elem->edge[iedge].ref == true){
				if(iedge == 0) {
					aux0.push_back(elem->nodes[EdgeVerticesMap[iedge][0]].id);
					aux0.push_back(elem->nodes[EdgeVerticesMap[iedge][1]].id);
					if(deb) mesh->part_nodes[elem->nodes[EdgeVerticesMap[iedge][0]].id] = 0;
					if(deb) mesh->part_nodes[elem->nodes[EdgeVerticesMap[iedge][1]].id] = 0;
				}
				if(iedge == 1) {
					aux1.push_back(elem->nodes[EdgeVerticesMap[iedge][0]].id);
					aux1.push_back(elem->nodes[EdgeVerticesMap[iedge][1]].id);
					if(deb) mesh->part_nodes[elem->nodes[EdgeVerticesMap[iedge][0]].id] = 1;
					if(deb) mesh->part_nodes[elem->nodes[EdgeVerticesMap[iedge][1]].id] = 1;
				}
				if(iedge == 2) {
					aux2.push_back(elem->nodes[EdgeVerticesMap[iedge][0]].id);
					aux2.push_back(elem->nodes[EdgeVerticesMap[iedge][1]].id);
					if(deb) mesh->part_nodes[elem->nodes[EdgeVerticesMap[iedge][0]].id] = 2;
					if(deb) mesh->part_nodes[elem->nodes[EdgeVerticesMap[iedge][1]].id] = 2;
				}
				if(iedge == 3) {
					aux3.push_back(elem->nodes[EdgeVerticesMap[iedge][0]].id);
					aux3.push_back(elem->nodes[EdgeVerticesMap[iedge][1]].id);
					if(deb) mesh->part_nodes[elem->nodes[EdgeVerticesMap[iedge][0]].id] = 3;
					if(deb) mesh->part_nodes[elem->nodes[EdgeVerticesMap[iedge][1]].id] = 3;
				}

				if(iedge == 4) {
					aux4.push_back(elem->nodes[EdgeVerticesMap[iedge][0]].id);
					aux4.push_back(elem->nodes[EdgeVerticesMap[iedge][1]].id);
					if(deb) mesh->part_nodes[elem->nodes[EdgeVerticesMap[iedge][0]].id] = 4;
					if(deb) mesh->part_nodes[elem->nodes[EdgeVerticesMap[iedge][1]].id] = 4;
				}
				if(iedge == 5) {
					aux5.push_back(elem->nodes[EdgeVerticesMap[iedge][0]].id);
					aux5.push_back(elem->nodes[EdgeVerticesMap[iedge][1]].id);
					if(deb) mesh->part_nodes[elem->nodes[EdgeVerticesMap[iedge][0]].id] = 5;
					if(deb) mesh->part_nodes[elem->nodes[EdgeVerticesMap[iedge][1]].id] = 5;
				}
				if(iedge == 6) {
					aux6.push_back(elem->nodes[EdgeVerticesMap[iedge][0]].id);
					aux6.push_back(elem->nodes[EdgeVerticesMap[iedge][1]].id);
					if(deb) mesh->part_nodes[elem->nodes[EdgeVerticesMap[iedge][0]].id] = 6;
					if(deb) mesh->part_nodes[elem->nodes[EdgeVerticesMap[iedge][1]].id] = 6;
				}
				if(iedge == 7) {
					aux7.push_back(elem->nodes[EdgeVerticesMap[iedge][0]].id);
					aux7.push_back(elem->nodes[EdgeVerticesMap[iedge][1]].id);
					if(deb) mesh->part_nodes[elem->nodes[EdgeVerticesMap[iedge][0]].id] = 7;
					if(deb) mesh->part_nodes[elem->nodes[EdgeVerticesMap[iedge][1]].id] = 7;
				}

				if(iedge == 8) {
					aux8.push_back(elem->nodes[EdgeVerticesMap[iedge][0]].id);
					aux8.push_back(elem->nodes[EdgeVerticesMap[iedge][1]].id);
					if(deb) mesh->part_nodes[elem->nodes[EdgeVerticesMap[iedge][0]].id] = 8;
					if(deb) mesh->part_nodes[elem->nodes[EdgeVerticesMap[iedge][1]].id] = 8;
				}
				if(iedge == 9) {
					aux9.push_back(elem->nodes[EdgeVerticesMap[iedge][0]].id);
					aux9.push_back(elem->nodes[EdgeVerticesMap[iedge][1]].id);
					if(deb) mesh->part_nodes[elem->nodes[EdgeVerticesMap[iedge][0]].id] = 9;
					if(deb) mesh->part_nodes[elem->nodes[EdgeVerticesMap[iedge][1]].id] = 9;
				}
				if(iedge == 10) {
					aux10.push_back(elem->nodes[EdgeVerticesMap[iedge][0]].id);
					aux10.push_back(elem->nodes[EdgeVerticesMap[iedge][1]].id);
					if(deb) mesh->part_nodes[elem->nodes[EdgeVerticesMap[iedge][0]].id] = 10;
					if(deb) mesh->part_nodes[elem->nodes[EdgeVerticesMap[iedge][1]].id] = 10;
				}
				if(iedge == 11) {
					aux11.push_back(elem->nodes[EdgeVerticesMap[iedge][0]].id);
					aux11.push_back(elem->nodes[EdgeVerticesMap[iedge][1]].id);
					if(deb) mesh->part_nodes[elem->nodes[EdgeVerticesMap[iedge][0]].id] = 11;
					if(deb) mesh->part_nodes[elem->nodes[EdgeVerticesMap[iedge][1]].id] = 11;
				}
			}
		}
	}

	for(int iedge = 0; iedge < 12 ; iedge++){
		std::vector<int> aux;

		if(true){
			if(iedge == 0){
				std::vector<double> temp;
				//printf("Tenho %d elementos no vetor\n", aux0.size());
				for(int ino = 0; ino < aux0.size(); ino++){
					temp.push_back(coords[3*aux0[ino]+0]);
					//printf("%f %d\n",temp[ino], aux0[ino]);
				}

				std::vector<int> V(aux0.size());
				int x=0;
				std::iota(V.begin(),V.end(),x++); //Initializing
				sort(V.begin(),V.end(), [&](int i,int j){return temp[i]<temp[j];} );

				//printf("Tenho %d elementos no vetor V\n", V.size());
				//for(int ino = 0; ino < aux0.size(); ino++){
				//	printf("%d\n",V[ino]);
				//}

				for(int ino = 0; ino < aux0.size(); ino++) aux.push_back(aux0[V[ino]]);
			}
			if(iedge == 1){
				std::vector<double> temp;
				for(int ino = 0; ino < aux1.size(); ino++) temp.push_back(coords[3*aux1[ino]+1]);

				std::vector<int> V(aux1.size());
				int x=0;
				std::iota(V.begin(),V.end(),x++); //Initializing
				sort(V.begin(),V.end(), [&](int i,int j){return temp[i]<temp[j];} );

				for(int ino = 0; ino < aux1.size(); ino++) aux.push_back(aux1[V[ino]]);
			}
			if(iedge == 2){
				std::vector<double> temp;
				for(int ino = 0; ino < aux2.size(); ino++) temp.push_back(coords[3*aux2[ino]+0]);

				std::vector<int> V(aux2.size());
				int x=0;
				std::iota(V.begin(),V.end(),x++); //Initializing
				sort(V.begin(),V.end(), [&](int i,int j){return temp[i]<temp[j];} );

				for(int ino = 0; ino < aux2.size(); ino++) aux.push_back(aux2[V[ino]]);
			}
			if(iedge == 3){
				std::vector<double> temp;
				for(int ino = 0; ino < aux3.size(); ino++) temp.push_back(coords[3*aux3[ino]+1]);

				std::vector<int> V(aux3.size());
				int x=0;
				std::iota(V.begin(),V.end(),x++); //Initializing
				sort(V.begin(),V.end(), [&](int i,int j){return temp[i]<temp[j];} );

				for(int ino = 0; ino < aux3.size(); ino++) aux.push_back(aux3[V[ino]]);
			}

			if(iedge == 4){
				std::vector<double> temp;
				for(int ino = 0; ino < aux4.size(); ino++) temp.push_back(coords[3*aux4[ino]+2]);

				std::vector<int> V(aux4.size());
				int x=0;
				std::iota(V.begin(),V.end(),x++); //Initializing
				sort(V.begin(),V.end(), [&](int i,int j){return temp[i]<temp[j];} );

				for(int ino = 0; ino < aux4.size(); ino++) aux.push_back(aux4[V[ino]]);
			}
			if(iedge == 5){
				std::vector<double> temp;
				for(int ino = 0; ino < aux5.size(); ino++) temp.push_back(coords[3*aux5[ino]+2]);

				std::vector<int> V(aux5.size());
				int x=0;
				std::iota(V.begin(),V.end(),x++); //Initializing
				sort(V.begin(),V.end(), [&](int i,int j){return temp[i]<temp[j];} );

				for(int ino = 0; ino < aux5.size(); ino++) aux.push_back(aux5[V[ino]]);
			}
			if(iedge == 6){
				std::vector<double> temp;
				for(int ino = 0; ino < aux6.size(); ino++) temp.push_back(coords[3*aux6[ino]+2]);

				std::vector<int> V(aux6.size());
				int x=0;
				std::iota(V.begin(),V.end(),x++); //Initializing
				sort(V.begin(),V.end(), [&](int i,int j){return temp[i]<temp[j];} );

				for(int ino = 0; ino < aux6.size(); ino++) aux.push_back(aux6[V[ino]]);
			}
			if(iedge == 7){
				std::vector<double> temp;
				for(int ino = 0; ino < aux7.size(); ino++) temp.push_back(coords[3*aux7[ino]+2]);

				std::vector<int> V(aux7.size());
				int x=0;
				std::iota(V.begin(),V.end(),x++); //Initializing
				sort(V.begin(),V.end(), [&](int i,int j){return temp[i]<temp[j];} );

				for(int ino = 0; ino < aux7.size(); ino++) aux.push_back(aux7[V[ino]]);
			}

			if(iedge == 8){
				std::vector<double> temp;
				for(int ino = 0; ino < aux8.size(); ino++) temp.push_back(coords[3*aux8[ino]+0]);

				std::vector<int> V(aux8.size());
				int x=0;
				std::iota(V.begin(),V.end(),x++); //Initializing
				sort(V.begin(),V.end(), [&](int i,int j){return temp[i]<temp[j];} );

				for(int ino = 0; ino < aux8.size(); ino++) aux.push_back(aux8[V[ino]]);
			}
			if(iedge == 9){
				std::vector<double> temp;
				for(int ino = 0; ino < aux9.size(); ino++) temp.push_back(coords[3*aux9[ino]+1]);

				std::vector<int> V(aux9.size());
				int x=0;
				std::iota(V.begin(),V.end(),x++); //Initializing
				sort(V.begin(),V.end(), [&](int i,int j){return temp[i]<temp[j];} );

				for(int ino = 0; ino < aux9.size(); ino++) aux.push_back(aux9[V[ino]]);
			}
			if(iedge == 10){
				std::vector<double> temp;
				for(int ino = 0; ino < aux10.size(); ino++) temp.push_back(coords[3*aux10[ino]+0]);

				std::vector<int> V(aux10.size());
				int x=0;
				std::iota(V.begin(),V.end(),x++); //Initializing
				sort(V.begin(),V.end(), [&](int i,int j){return temp[i]<temp[j];} );

				for(int ino = 0; ino < aux10.size(); ino++) aux.push_back(aux10[V[ino]]);
			}
			if(iedge == 11){
				std::vector<double> temp;
				for(int ino = 0; ino < aux11.size(); ino++) temp.push_back(coords[3*aux11[ino]+1]);

				std::vector<int> V(aux11.size());
				int x=0;
				std::iota(V.begin(),V.end(),x++); //Initializing
				sort(V.begin(),V.end(), [&](int i,int j){return temp[i]<temp[j];} );

				for(int ino = 0; ino < aux11.size(); ino++) aux.push_back(aux11[V[ino]]);
			}
		}

		int niter = 0;
		double tol = 1;
		double factor = 0.01;
		while(niter < 100 && tol >= 0.001){
			for(int ino = 1; ino < (aux.size()-1); ino++){
				node_t key;
				size_t position;
				int node = aux[ino];
				key.coord[0] = coords[3*node+0];
				key.coord[1] = coords[3*node+1];
				key.coord[2] = coords[3*node+2];
				key.node_id = node;

				bool nodel = sc_hash_array_lookup(hash_FixedNodes, &key, &position);
				if(nodel){

				}else{
					int node0 = aux[ino-1];
					int node2 = aux[ino+1];
					double a = coords[3*node+0];
					double b = coords[3*node+1];
					double c = coords[3*node+2];
					coords[3*node+0] = ((coords[3*node0+0] + coords[3*node2+0]) - 2*coords[3*node+0])*factor + coords[3*node+0];
					coords[3*node+1] = ((coords[3*node0+1] + coords[3*node2+1]) - 2*coords[3*node+1])*factor + coords[3*node+1];
					coords[3*node+2] = ((coords[3*node0+2] + coords[3*node2+2]) - 2*coords[3*node+2])*factor + coords[3*node+2];
					a =- coords[3*node+0];
					b =- coords[3*node+1];
					c =- coords[3*node+2];
					tol = sqrt(a*a + b*b + c*c);
					//printf("%f\n",tol);
				}
			}
			niter++;
		}
	}

	//printf("     Please check 1D optimization\n");
	//printf("     WARNING: The 1D Optimization was made only for the vertical lines\n");


}

void OptSurface(hexa_tree_t* mesh, std::vector<double>& coords, sc_hash_array_t* hash_FixedNodes){

	Mesquite::MsqPrintError err(std::cout);

	for(int isurf = 0; isurf < 4 ; isurf ++){

		//hash of the fixed nodes
		bool clamped = true;
		sc_hash_array_t* hash_Fixed = sc_hash_array_new(sizeof(node_t), edge_hash_fn, edge_equal_fn, &clamped);

		for(int ino = 0; ino < hash_FixedNodes->a.elem_count; ino++){
			size_t position;
			node_t *r;
			node_t key;
			node_t* node = (node_t*) sc_array_index (&hash_FixedNodes->a, ino);
			key.coord[0] = coords[3*node->node_id+0];
			key.coord[1] = coords[3*node->node_id+1];
			key.coord[2] = coords[3*node->node_id+2];
			key.node_id = node->node_id;
			r = (node_t*) sc_hash_array_insert_unique(hash_Fixed, &key, &position);
			if (r != NULL) {
				r->coord[0] = key.coord[0];
				r->coord[1] = key.coord[1];
				r->coord[2] = key.coord[2];
				r->node_id = key.node_id;
			} else {

			}
		}

		std::vector<int> elem_aux;
		std::vector<int> aux;
		int nelem=0;
		int nvertices=0;

		//hash for the surface nodes
		sc_hash_array_t* hash_SurfaceNodes = sc_hash_array_new(sizeof(node_t), edge_hash_fn, edge_equal_fn, &clamped);

		if(isurf == 0){
			for(int  iel =0; iel < mesh->outsurf.elem_count; iel++){
				octant_t * elem = (octant_t*) sc_array_index(&mesh->outsurf, iel);

				if(elem->surf[isurf].ext){
					elem_aux.push_back(elem->id);
					for(int ino = 0; ino < 4; ino++){
						node_t key;
						size_t position;
						key.coord[0] = coords[3*elem->nodes[FaceNodesMap[isurf][ino]].id+0];
						key.coord[1] = coords[3*elem->nodes[FaceNodesMap[isurf][ino]].id+1];
						key.coord[2] = coords[3*elem->nodes[FaceNodesMap[isurf][ino]].id+2];
						key.node_id = elem->nodes[FaceNodesMap[isurf][ino]].id;
						node_t* r = (node_t*) sc_hash_array_insert_unique(hash_SurfaceNodes, &key, &position);
						if (r != NULL) {
							r->coord[0] = key.coord[0];
							r->coord[1] = key.coord[1];
							r->coord[2] = key.coord[2];
							r->node_id = key.node_id;
						} else {
						}
					}

					for(int iedge = 0; iedge < 4; iedge++){
						if(elem->edge[FaceEdgesMap[isurf][iedge]].ref == true){
							for(int ino = 0; ino < 2; ino++){
								int node = elem->nodes[EdgeVerticesMap[FaceEdgesMap[isurf][iedge]][ino]].id;
								size_t position;
								node_t key;
								key.coord[0] = coords[3*node+0];
								key.coord[1] = coords[3*node+1];
								key.coord[2] = coords[3*node+2];
								key.node_id = node;
								node_t* r = (node_t*) sc_hash_array_insert_unique(hash_Fixed, &key, &position);
								if (r != NULL) {
									r->coord[0] = key.coord[0];
									r->coord[1] = key.coord[1];
									r->coord[2] = key.coord[2];
									r->node_id = key.node_id;
								}
							}
						}
					}
				}
			}
		}

		if(isurf == 1){
			for(int  iel =0; iel < mesh->outsurf.elem_count; iel++){
				octant_t * elem = (octant_t*) sc_array_index(&mesh->outsurf, iel);

				if(elem->surf[isurf].ext){
					elem_aux.push_back(elem->id);
					for(int ino = 0; ino < 4; ino++){
						node_t key;
						size_t position;
						key.coord[0] = coords[3*elem->nodes[FaceNodesMap[isurf][ino]].id+0];
						key.coord[1] = coords[3*elem->nodes[FaceNodesMap[isurf][ino]].id+1];
						key.coord[2] = coords[3*elem->nodes[FaceNodesMap[isurf][ino]].id+2];
						key.node_id = elem->nodes[FaceNodesMap[isurf][ino]].id;
						node_t* r = (node_t*) sc_hash_array_insert_unique(hash_SurfaceNodes, &key, &position);
						if (r != NULL) {
							r->coord[0] = key.coord[0];
							r->coord[1] = key.coord[1];
							r->coord[2] = key.coord[2];
							r->node_id = key.node_id;
						} else {
						}
					}
					for(int iedge = 0; iedge < 4; iedge++){
						if(elem->edge[FaceEdgesMap[isurf][iedge]].ref == true){
							for(int ino = 0; ino < 2; ino++){
								int node = elem->nodes[EdgeVerticesMap[FaceEdgesMap[isurf][iedge]][ino]].id;
								size_t position;
								node_t key;
								key.coord[0] = coords[3*node+0];
								key.coord[1] = coords[3*node+1];
								key.coord[2] = coords[3*node+2];
								key.node_id = node;
								node_t* r = (node_t*) sc_hash_array_insert_unique(hash_Fixed, &key, &position);
								if (r != NULL) {
									r->coord[0] = key.coord[0];
									r->coord[1] = key.coord[1];
									r->coord[2] = key.coord[2];
									r->node_id = key.node_id;
								}
							}
						}
					}
				}
			}
		}

		if(isurf == 2){
			for(int  iel =0; iel < mesh->outsurf.elem_count; iel++){
				octant_t * elem = (octant_t*) sc_array_index(&mesh->outsurf, iel);

				if(elem->surf[isurf].ext){
					elem_aux.push_back(elem->id);
					for(int ino = 0; ino < 4; ino++){
						node_t key;
						size_t position;
						key.coord[0] = coords[3*elem->nodes[FaceNodesMap[isurf][ino]].id+0];
						key.coord[1] = coords[3*elem->nodes[FaceNodesMap[isurf][ino]].id+1];
						key.coord[2] = coords[3*elem->nodes[FaceNodesMap[isurf][ino]].id+2];
						key.node_id = elem->nodes[FaceNodesMap[isurf][ino]].id;
						node_t* r = (node_t*) sc_hash_array_insert_unique(hash_SurfaceNodes, &key, &position);
						if (r != NULL) {
							r->coord[0] = key.coord[0];
							r->coord[1] = key.coord[1];
							r->coord[2] = key.coord[2];
							r->node_id = key.node_id;
						} else {
						}
					}
					for(int iedge = 0; iedge < 4; iedge++){
						if(elem->edge[FaceEdgesMap[isurf][iedge]].ref == true){
							for(int ino = 0; ino < 2; ino++){
								int node = elem->nodes[EdgeVerticesMap[FaceEdgesMap[isurf][iedge]][ino]].id;
								size_t position;
								node_t key;
								key.coord[0] = coords[3*node+0];
								key.coord[1] = coords[3*node+1];
								key.coord[2] = coords[3*node+2];
								key.node_id = node;
								node_t* r = (node_t*) sc_hash_array_insert_unique(hash_Fixed, &key, &position);
								if (r != NULL) {
									r->coord[0] = key.coord[0];
									r->coord[1] = key.coord[1];
									r->coord[2] = key.coord[2];
									r->node_id = key.node_id;
								}
							}
						}
					}
				}
			}
		}

		if(isurf == 3){
			for(int  iel =0; iel < mesh->outsurf.elem_count; iel++){
				octant_t * elem = (octant_t*) sc_array_index(&mesh->outsurf, iel);

				if(elem->surf[isurf].ext){
					elem_aux.push_back(elem->id);
					for(int ino = 0; ino < 4; ino++){
						node_t key;
						size_t position;
						key.coord[0] = coords[3*elem->nodes[FaceNodesMap[isurf][ino]].id+0];
						key.coord[1] = coords[3*elem->nodes[FaceNodesMap[isurf][ino]].id+1];
						key.coord[2] = coords[3*elem->nodes[FaceNodesMap[isurf][ino]].id+2];
						key.node_id = elem->nodes[FaceNodesMap[isurf][ino]].id;
						node_t* r = (node_t*) sc_hash_array_insert_unique(hash_SurfaceNodes, &key, &position);
						if (r != NULL) {
							r->coord[0] = key.coord[0];
							r->coord[1] = key.coord[1];
							r->coord[2] = key.coord[2];
							r->node_id = key.node_id;
						} else {
						}
					}
					for(int iedge = 0; iedge < 4; iedge++){
						if(elem->edge[FaceEdgesMap[isurf][iedge]].ref == true){
							for(int ino = 0; ino < 2; ino++){
								int node = elem->nodes[EdgeVerticesMap[FaceEdgesMap[isurf][iedge]][ino]].id;
								size_t position;
								node_t key;
								key.coord[0] = coords[3*node+0];
								key.coord[1] = coords[3*node+1];
								key.coord[2] = coords[3*node+2];
								key.node_id = node;
								node_t* r = (node_t*) sc_hash_array_insert_unique(hash_Fixed, &key, &position);
								if (r != NULL) {
									r->coord[0] = key.coord[0];
									r->coord[1] = key.coord[1];
									r->coord[2] = key.coord[2];
									r->node_id = key.node_id;
								}
							}
						}
					}
				}
			}
		}

		if(isurf == 4){
			for(int  iel =0; iel < mesh->outsurf.elem_count; iel++){
				octant_t * elem = (octant_t*) sc_array_index(&mesh->outsurf, iel);

				if(elem->surf[isurf].ext){
					elem_aux.push_back(elem->id);
					for(int ino = 0; ino < 4; ino++){
						node_t key;
						size_t position;
						key.coord[0] = coords[3*elem->nodes[FaceNodesMap[isurf][ino]].id+0];
						key.coord[1] = coords[3*elem->nodes[FaceNodesMap[isurf][ino]].id+1];
						key.coord[2] = coords[3*elem->nodes[FaceNodesMap[isurf][ino]].id+2];
						key.node_id = elem->nodes[FaceNodesMap[isurf][ino]].id;
						node_t* r = (node_t*) sc_hash_array_insert_unique(hash_SurfaceNodes, &key, &position);
						if (r != NULL) {
							r->coord[0] = key.coord[0];
							r->coord[1] = key.coord[1];
							r->coord[2] = key.coord[2];
							r->node_id = key.node_id;
						} else {
						}
					}
					for(int iedge = 0; iedge < 4; iedge++){
						if(elem->edge[FaceEdgesMap[isurf][iedge]].ref == true){
							for(int ino = 0; ino < 2; ino++){
								int node = elem->nodes[EdgeVerticesMap[FaceEdgesMap[isurf][iedge]][ino]].id;
								size_t position;
								node_t key;
								key.coord[0] = coords[3*node+0];
								key.coord[1] = coords[3*node+1];
								key.coord[2] = coords[3*node+2];
								key.node_id = node;
								node_t* r = (node_t*) sc_hash_array_insert_unique(hash_Fixed, &key, &position);
								if (r != NULL) {
									r->coord[0] = key.coord[0];
									r->coord[1] = key.coord[1];
									r->coord[2] = key.coord[2];
									r->node_id = key.node_id;
								}
							}
						}
					}
				}
			}
		}

		if(isurf == 5){
			for(int  iel =0; iel < mesh->outsurf.elem_count; iel++){
				octant_t * elem = (octant_t*) sc_array_index(&mesh->outsurf, iel);

				if(elem->surf[isurf].ext){
					elem_aux.push_back(elem->id);
					for(int ino = 0; ino < 4; ino++){
						node_t key;
						size_t position;
						key.coord[0] = coords[3*elem->nodes[FaceNodesMap[isurf][ino]].id+0];
						key.coord[1] = coords[3*elem->nodes[FaceNodesMap[isurf][ino]].id+1];
						key.coord[2] = coords[3*elem->nodes[FaceNodesMap[isurf][ino]].id+2];
						key.node_id = elem->nodes[FaceNodesMap[isurf][ino]].id;
						node_t* r = (node_t*) sc_hash_array_insert_unique(hash_SurfaceNodes, &key, &position);
						if (r != NULL) {
							r->coord[0] = key.coord[0];
							r->coord[1] = key.coord[1];
							r->coord[2] = key.coord[2];
							r->node_id = key.node_id;
						} else {
						}
					}

					for(int iedge = 0; iedge < 4; iedge++){
						if(elem->edge[FaceEdgesMap[isurf][iedge]].ref == true){
							for(int ino = 0; ino < 2; ino++){
								int node = elem->nodes[EdgeVerticesMap[FaceEdgesMap[isurf][iedge]][ino]].id;
								size_t position;
								node_t key;
								key.coord[0] = coords[3*node+0];
								key.coord[1] = coords[3*node+1];
								key.coord[2] = coords[3*node+2];
								key.node_id = node;
								node_t* r = (node_t*) sc_hash_array_insert_unique(hash_Fixed, &key, &position);
								if (r != NULL) {
									r->coord[0] = key.coord[0];
									r->coord[1] = key.coord[1];
									r->coord[2] = key.coord[2];
									r->node_id = key.node_id;
								}
							}
						}
					}
				}
			}
		}

		nelem = elem_aux.size();
		nvertices = hash_SurfaceNodes->a.elem_count;
		printf("     Mesh normal to %d with: %d elements and %d vertex\n",isurf, nelem,nvertices);

		std::vector<int> conn(4*nelem);
		std::vector<int> conn_aux(4*nelem);
		std::vector<double> Scoors(3*nvertices);
		bool   *fixed_nodes = (bool*)malloc(nvertices*sizeof(bool));

		// montando a malha para passar para o mesquite...
		// conectividade na malha "global"
		int c = 0;
		for(int iel =0; iel < elem_aux.size(); iel++){
			octant_t * elem = (octant_t*) sc_array_index(&mesh->elements, elem_aux[iel]);
			for(int ino = 0; ino < 4; ino++){
				conn_aux[c] = elem->nodes[FaceNodesMapRef[isurf][ino]].id;
				c++;
			}
		}


		//nos fixos e coordenadas dos nos
		for(int ive = 0 ; ive < nvertices ;ive++ ){
			fixed_nodes[ive] = false;
			node_t * node = (node_t*) sc_array_index(&hash_SurfaceNodes->a, ive);
			Scoors[3*ive+0] = coords[3*node->node_id+0];
			Scoors[3*ive+1] = coords[3*node->node_id+1];
			Scoors[3*ive+2] = coords[3*node->node_id+2];
			aux.push_back(node->node_id);

			node_t key;
			size_t position;
			key.coord[0] = coords[3*node->node_id+0];
			key.coord[1] = coords[3*node->node_id+1];
			key.coord[2] = coords[3*node->node_id+2];
			bool tre = sc_hash_array_lookup(hash_Fixed,&key,&position);
			if(tre){
				fixed_nodes[ive] = true;
				mesh->part_nodes[node->node_id] = 1;
			}

			//for(int j = 0; j<hash_FixedNodes->a.elem_count;j++){
			//	node_t* rr= (node_t*) sc_array_index (&hash_FixedNodes->a, j);
			//	if(rr->node_id==node->node_id){
			//		fixed_nodes[ive] = true;
			//		mesh->part_nodes[node->node_id] = 1;
			//	}
			//}


		}

		//refazendo a conectividade com nos de 0 a nvertices...
		//passando uma malha "local"
		for(int i = 0; i < aux.size(); i++){
			for(int j = 0; j < conn_aux.size(); j++){
				if(conn_aux[j]==aux[i]){
					conn[j] = i;
				}
			}
		}

		if(nelem!=0){
			//finalmente mandando a malha pro mesquite...
			Mesquite::MeshImpl mesq_mesh0(nvertices,nelem,Mesquite::QUADRILATERAL, &fixed_nodes[0], &Scoors[0], &conn[0]);

			Vector3D point(Scoors[0],Scoors[1],Scoors[2]);
			Vector3D normal;
			if(isurf == 0){
				normal.set(1,0,0);
				//Mesquite::PlanarDomain plane(normal, point);
			}else if(isurf == 1){
				normal.set(1,0,0);
				//Mesquite::PlanarDomain plane(Vector3D(-1,0,0), point);
			}else if(isurf == 2){
				normal.set(0,-1,0);
				//Mesquite::PlanarDomain plane(Vector3D(0,1,0), point);
			}else if(isurf == 3){
				normal.set(0,1,0);
				//Mesquite::PlanarDomain plane(Vector3D(0,-1,0), point);
			}else if(isurf == 4){
				normal.set(0,0,-1);
				//Mesquite::PlanarDomain plane(Vector3D(0,0,1), point);
			}else if(isurf == 5){
				normal.set(0,0,1);
				//Mesquite::PlanarDomain plane(Vector3D(0,0,-1), point);
			}
			Mesquite::PlanarDomain plane(normal, point);
			MeshDomainAssoc mesh_and_domain0 = MeshDomainAssoc(&mesq_mesh0, &plane);
			//printf("Status da criacao do mesh_domain para normal %d: %d\n",axis, mesh_and_domain0.are_compatible());
			assert(mesh_and_domain0.are_compatible()==1);

			if(false){
				LaplaceWrapper lp_wrapper0;
				//lp_wrapper0.set_vertex_movement_limit_factor(1.e-6);
				lp_wrapper0.set_iteration_limit(20);
				lp_wrapper0.set_cpu_time_limit(120);
				lp_wrapper0.enable_culling(true);
				lp_wrapper0.run_instructions( &mesh_and_domain0, err );
				if (err) std::cout << err << std::endl;
			}

			if(false){
				ShapeImprover smoother;
				IdealWeightInverseMeanRatio extra_metric;
				//smoother.quality_assessor().add_quality_assessment(&extra_metric);
				smoother.set_cpu_time_limit(120);
				smoother.run_instructions( &mesh_and_domain0, err );
				if (err) std::cout << err << std::endl;
			}

			if(false){
				UntangleWrapper::UntangleMetric metric0 = UntangleWrapper::BETA;
				UntangleWrapper un_wrapper0 (metric0);
				//un_wrapper0.set_vertex_movement_limit_factor( 0.005 );
				un_wrapper0.set_outer_iteration_limit(1000);
				un_wrapper0.run_instructions( &mesh_and_domain0, err );
				if (err) std::cout << err << std::endl;
			}

			if(true){
				// creates an intruction queue
				InstructionQueue queue1;

				// creates a mean ratio quality metric ...
				ConditionNumberQualityMetric shape_metric;
				EdgeLengthQualityMetric lapl_met;
				lapl_met.set_averaging_method(QualityMetric::RMS);

				// creates the laplacian smoother  procedures
				LaplacianSmoother lapl1;
				QualityAssessor stop_qa=QualityAssessor(&shape_metric);
				stop_qa.add_quality_assessment(&lapl_met);
				stop_qa.disable_printing_results();

				//**************Set stopping criterion****************
				TerminationCriterion sc2;
				sc2.add_iteration_limit( 10 );
				sc2.add_relative_vertex_movement(1e-5);
				lapl1.set_outer_termination_criterion(&sc2);

				// adds 1 pass of pass1 to mesh_set1
				queue1.add_quality_assessor(&stop_qa,err);
				queue1.set_master_quality_improver(&lapl1, err);
				queue1.add_quality_assessor(&stop_qa,err);
				// adds 1 passes of pass2 to mesh_set1
				//  mesh_set1.add_quality_pass(pass2);

				// launches optimization on mesh_set1
				queue1.run_instructions(&mesh_and_domain0, err);
			}

			if(false){
				ShapeImprover smoother;
				IdealWeightInverseMeanRatio extra_metric;
				smoother.quality_assessor().add_quality_assessment(&extra_metric);
				smoother.set_cpu_time_limit(120);
				smoother.run_instructions( &mesh_and_domain0, err );

				//SizeAdaptShapeWrapper smoother(1e-2);
				//smoother.run_instructions( &mesh_and_domain0, err);

				//SizeAdaptShapeWrapper smoother(1e-2);
				//MeshDomainAssoc mesh_and_domain = MeshDomainAssoc(&mesh, &geom);
				//smoother.run_instructions( &mesh_and_domain, err);
			}

			if(false){
				// creates an intruction queue
				InstructionQueue queue2;

				// creates a mean ratio quality metric ...
				ConditionNumberQualityMetric shape_metric2;
				UntangleBetaQualityMetric untangle(2);
				Randomize pass0(.05);
				// ... and builds an objective function with it
				//LInfTemplate* obj_func = new LInfTemplate(shape_metric);
				LInfTemplate obj_func(&untangle);
				LPtoPTemplate obj_func2(&shape_metric2, 2, err);
				// creates the steepest descent optimization procedures
				ConjugateGradient pass1( &obj_func, err );

				//SteepestDescent* pass2 = new SteepestDescent( obj_func2 );
				ConjugateGradient pass2( &obj_func2, err );
				pass2.use_element_on_vertex_patch();
				pass2.use_global_patch();
				QualityAssessor stop_qa=QualityAssessor(&shape_metric2);
				QualityAssessor stop_qa2=QualityAssessor(&shape_metric2);
				if (true) {
					stop_qa.disable_printing_results();
					stop_qa2.disable_printing_results();
				}

				stop_qa.add_quality_assessment(&untangle);
				// **************Set stopping criterion**************
				//untangle beta should be 0 when untangled
				TerminationCriterion sc11;
				sc11.add_relative_quality_improvement( 0.000001 );
				TerminationCriterion sc31;
				sc31.add_iteration_limit( 10 );
				TerminationCriterion sc_rand;
				sc_rand.add_iteration_limit( 1 );

				//StoppingCriterion sc1(&stop_qa,-1.0,.0000001);
				//StoppingCriterion sc3(&stop_qa2,.9,1.00000001);
				//StoppingCriterion sc2(StoppingCriterion::NUMBER_OF_PASSES,10);
				//StoppingCriterion sc_rand(StoppingCriterion::NUMBER_OF_PASSES,1);
				//either until untangled or 10 iterations
				pass0.set_outer_termination_criterion(&sc_rand);
				pass1.set_outer_termination_criterion(&sc11);
				pass2.set_inner_termination_criterion(&sc31);

				// adds 1 pass of pass1 to mesh_set1
				queue2.add_quality_assessor(&stop_qa,err);
				//queue1.add_preconditioner(pass0,err);MSQ_CHKERR(err);
				//queue1.add_preconditioner(pass1,err);MSQ_CHKERR(err);
				//queue1.set_master_quality_improver(pass2, err); MSQ_CHKERR(err);
				queue2.set_master_quality_improver(&pass1, err);
				queue2.add_quality_assessor(&stop_qa2,err);

				queue2.run_instructions(&mesh_and_domain0, err);
			}

			std::vector<MeshImpl::VertexHandle> vertices;
			mesh_and_domain0.get_mesh()->get_all_vertices(vertices,err);

			if(isurf != 5){
				for (int ino = 0; ino < vertices.size() ; ino++){
					Mesh::VertexHandle vertex = vertices[ino];
					MsqVertex aux_msq;
					mesh_and_domain0.get_mesh()->vertices_get_coordinates( &vertex, &aux_msq, 1, err );
					node_t * node = (node_t*) sc_array_index(&hash_SurfaceNodes->a, ino);
					coords[3*node->node_id+0] = aux_msq[0];
					coords[3*node->node_id+1] = aux_msq[1];
					coords[3*node->node_id+2] = aux_msq[2];
				}
				mesq_mesh0.clear();
			}else{
				for (int ino = 0; ino < vertices.size() ; ino++){
					Mesh::VertexHandle vertex = vertices[ino];
					MsqVertex aux_msq;
					mesh_and_domain0.get_mesh()->vertices_get_coordinates( &vertex, &aux_msq, 1, err );
					node_t * node = (node_t*) sc_array_index(&hash_SurfaceNodes->a, ino);
					coords[3*node->node_id+0] = aux_msq[0];
					coords[3*node->node_id+1] = aux_msq[1];
					//coords[3*node->node_id+2] = aux_msq[2];
				}
				mesq_mesh0.clear();
			}
		}else{
			printf("Please check mesquite_interface! set of elements = 0.\n");
		}

		aux.clear();
		elem_aux.clear();
		conn.clear();
		conn_aux.clear();
		Scoors.clear();
		free(fixed_nodes);
		sc_hash_array_destroy(hash_SurfaceNodes);

	}

}

void RedoCoords(hexa_tree_t* mesh, std::vector<double>& coords){

	GtsPoint *p;
	double dx, dy, dz;
	double d;
	double zmax;

	// Get the surface bounding box
	mesh->tdata.bbox = gts_bbox_surface(gts_bbox_class(), mesh->tdata.s);
	if (mesh->mpi_rank == 0) {
		printf("         Bounding box: \n");
		printf("         x ranges from %f to %f\n", mesh->tdata.bbox->x1, mesh->tdata.bbox->x2);
		printf("         y ranges from %f to %f\n", mesh->tdata.bbox->y1, mesh->tdata.bbox->y2);
		printf("         z ranges from %f to %f\n", mesh->tdata.bbox->z1, mesh->tdata.bbox->z2);
	}

	double Lx = (mesh->tdata.bbox->x2 - mesh->tdata.bbox->x1);
	double Ly = (mesh->tdata.bbox->y2 - mesh->tdata.bbox->y1);
	double zmin = ((Lx < Ly) ? -Lx : -Ly);

	// Get grid-spacing at x and y direction
	dx = 4*(mesh->tdata.bbox->x2 - mesh->tdata.bbox->x1) / (double) mesh->ncellx;
	dy = 4*(mesh->tdata.bbox->y2 - mesh->tdata.bbox->y1) / (double) mesh->ncelly;

	// Build the bounding box tree
	mesh->tdata.bbt = gts_bb_tree_surface(mesh->tdata.s);

	p = gts_point_new(gts_point_class(), 0.0, 0.0, mesh->tdata.bbox->z2);

	for (int i = 0; i < mesh->nodes.elem_count; ++i) {
		octant_node_t* n = (octant_node_t*) sc_array_index(&mesh->nodes, i);
		if(n->z <= 24){
			p->x = mesh->tdata.bbox->x1 + (n->x/12)*dx;
			p->y = mesh->tdata.bbox->y1 + (n->y/12)*dy;

			d = gts_bb_tree_point_distance(mesh->tdata.bbt, p, distance, NULL);
			zmax = mesh->tdata.bbox->z2 - d;

			dz = (zmax - zmin) / (double) mesh->ncellz;
			coords[i * 3 + 2] = zmax - (n->z/12) * dz;
		}
	}
}

void OptVolume(hexa_tree_t* mesh, std::vector<double>& coords, sc_hash_array_t* hash_FixedNodes){
	Mesquite::MsqPrintError err(std::cout);

	int nvertices = mesh->local_n_nodes;
	int nelem     = mesh->local_n_elements;

	std::vector<int> conn(8*nelem);
	bool   *fixed_nodes = (bool*)malloc(nvertices*sizeof(bool));
	size_t *gid         = (size_t*) malloc(nvertices*sizeof(size_t));
	int    *pid         = (int*) malloc (nvertices*sizeof(int));
	int    *mid         = (int*) malloc (nelem*sizeof(int));

	int c = 0;
	int tmp[8];
	for(int  iel = 0; iel < mesh->elements.elem_count; iel++){
		octant_t * elem = (octant_t*) sc_array_index(&mesh->elements, iel);
		mid[iel] = elem->n_mat;
		for(int j = 0; j < 8; j++){
			tmp[j] = elem->nodes[j].id;
		}
		conn[c] = tmp[4]; c++;
		conn[c] = tmp[5]; c++;
		conn[c] = tmp[6]; c++;
		conn[c] = tmp[7]; c++;
		conn[c] = tmp[0]; c++;
		conn[c] = tmp[1]; c++;
		conn[c] = tmp[2]; c++;
		conn[c] = tmp[3]; c++;

		if(elem->tem == -1 && true){
			for(int j = 0; j < 8; j++){
				fixed_nodes[elem->nodes[j].id] = true;
			}
		}
	}

	assert(mesh->nodes.elem_count == mesh->local_n_nodes);

	int bnd_nodes = 0;
	//fixed nodes
	for(int iel = 0; iel < mesh->outsurf.elem_count; iel++){
		octant_t * elem = (octant_t*) sc_array_index(&mesh->outsurf, iel);

		for(int isurf = 0 ; isurf < 6; isurf++){
			if(elem->surf[isurf].ext){
				for(int ino = 0; ino < 4; ino++){
					fixed_nodes[elem->nodes[FaceNodesMap[isurf][ino]].id] = true;
					//mesh->part_nodes[elem->nodes[FaceNodesMap[isurf][ino]].id] = 1;
				}
			}
		}
	}

	for(int ino = 0; ino < hash_FixedNodes->a.elem_count; ino++){
		node_t* node = (node_t*) sc_array_index (&hash_FixedNodes->a, ino);
		fixed_nodes[node->node_id] = true;
	}

	//for(int ino = 0; ino < mesh->shared_nodes.elem_count; ino++){
	//	shared_node_t * node = (shared_node_t*) sc_array_index(&mesh->shared_nodes, ino);
	//	fixed_nodes[node->id] = true;
	//}

	Mesquite::MeshImpl mesq_mesh(nvertices,nelem,Mesquite::HEXAHEDRON, &fixed_nodes[0], &coords[0], &conn[0]);

	if(false){
		/*
		// creates an intruction queue
		InstructionQueue queue;


		// calculate average lambda for mesh
		ReferenceMesh ref_mesh( MeshImpl );
		RefMeshTargetCalculator W_0( &ref_mesh );
		SimpleStats lambda_stats;
		Settings* settings
		MeshUtil tool(mesh, settings);
		tool.lambda_distribution( lambda_stats, err ); MSQ_ERRRTN(err);
		double lambda = lambda_stats.average();



		// create objective function
		IdealShapeTarget W_i;
		LambdaConstant W( lambda, &W_i );
		TShapeSizeB1 tm;
		TQualityMetric mu_0( &W, &tm );
		ElementPMeanP mu( 1.0, &mu_0 );
		PMeanPTemplate of( 1.0, &mu );
		*/
		/*
		// create quality assessor
		EdgeLengthMetric len(0.0);
		ConditionNumberQualityMetric shape_metric;
		QualityAssessor qa=QualityAssessor(&shape_metric);
		qa->add_quality_assessment( &mu );
		qa->add_quality_assessment( &len );
		queue.add_quality_assessor( qa, err );

		// create solver
		TrustRegion solver( &of );
		TerminationCriterion tc, ptc;
		tc.add_absolute_vertex_movement( maxVtxMovement );
		tc.add_iteration_limit( iterationLimit );
		ptc.add_iteration_limit( pmesh ? parallelIterations : 1 );
		solver.set_inner_termination_criterion( &tc );
		solver.set_outer_termination_criterion( &ptc );
		q.set_master_quality_improver( &solver, err ); MSQ_ERRRTN(err);
		q.add_quality_assessor( qa, err );

		// Optimize mesh
		q.run_common( mesh_and_domain, pmesh, settings, err ); MSQ_CHKERR(err);
		 */
	}

	if(false){
		// creates an intruction queue
		InstructionQueue queue1;

		// creates a mean ratio quality metric ...
		ConditionNumberQualityMetric shape_metric;
		UntangleBetaQualityMetric untangle(2);
		//Randomize pass0(.05);
		// ... and builds an objective function with it
		//LInfTemplate* obj_func = new LInfTemplate(shape_metric);
		LInfTemplate obj_func(&untangle);
		//LPtoPTemplate obj_func2(&shape_metric, 2, err);

		// creates the steepest descent optimization procedures
		ConjugateGradient pass1( &obj_func, err );

		//SteepestDescent* pass2 = new SteepestDescent( obj_func2 );
		//ConjugateGradient pass2( &obj_func2, err );

		//pass2.use_element_on_vertex_patch();

		//pass2.use_global_patch();

		QualityAssessor stop_qa=QualityAssessor(&shape_metric);
		//QualityAssessor stop_qa2=QualityAssessor(&shape_metric);

		//stop_qa.disable_printing_results();
		//stop_qa2.disable_printing_results();

		stop_qa.add_quality_assessment(&untangle);
		// **************Set stopping criterion**************
		//untangle beta should be 0 when untangled
		TerminationCriterion sc1;
		sc1.add_iteration_limit( 10 );
		//TerminationCriterion sc_rand;
		//sc_rand.add_iteration_limit( 1 );

		//StoppingCriterion sc1(&stop_qa,-1.0,.0000001);
		//StoppingCriterion sc3(&stop_qa2,.9,1.00000001);
		//StoppingCriterion sc2(StoppingCriterion::NUMBER_OF_PASSES,10);
		//StoppingCriterion sc_rand(StoppingCriterion::NUMBER_OF_PASSES,1);
		//either until untangled or 10 iterations
		//pass0.set_outer_termination_criterion(&sc_rand);
		pass1.set_outer_termination_criterion(&sc1);

		// adds 1 pass of pass1 to mesh_set1
		queue1.add_quality_assessor(&stop_qa,err);
		//queue1.add_preconditioner(pass0,err);MSQ_CHKERR(err);
		//queue1.add_preconditioner(pass1,err);MSQ_CHKERR(err);
		//queue1.set_master_quality_improver(pass2, err); MSQ_CHKERR(err);
		//queue1.set_master_quality_improver(&pass1, err);
		//queue1.add_quality_assessor(&stop_qa2,err);


		// launches optimization on mesh_set1
		queue1.run_instructions(&mesq_mesh, err);


		//print_timing_diagnostics(std::cout);
	}

	if(true){
		// creates an intruction queue
		InstructionQueue queue1;

		// creates a mean ratio quality metric ...
		ConditionNumberQualityMetric shape_metric;
		EdgeLengthQualityMetric lapl_met;
		lapl_met.set_averaging_method(QualityMetric::RMS);

		// creates the laplacian smoother  procedures
		//Here we use SmartLaplacianSmoother
		//it tries to avoid the inversion of the element...
		SmartLaplacianSmoother lapl1;
		//LaplacianSmoother lapl1;
		QualityAssessor stop_qa=QualityAssessor(&shape_metric);
		stop_qa.add_quality_assessment(&lapl_met);
		stop_qa.disable_printing_results();

		//**************Set stopping criterion****************
		TerminationCriterion sc2;
		sc2.add_iteration_limit( 30 );
		sc2.add_cpu_time(120);
		//sc2.add_relative_vertex_movement(1e-5);
		lapl1.set_outer_termination_criterion(&sc2);

		// adds 1 pass of pass1 to mesh_set1
		queue1.add_quality_assessor(&stop_qa,err);
		queue1.set_master_quality_improver(&lapl1, err);
		queue1.add_quality_assessor(&stop_qa,err);
		// adds 1 passes of pass2 to mesh_set1
		//  mesh_set1.add_quality_pass(pass2);

		// launches optimization on mesh_set1
		queue1.run_instructions(&mesq_mesh, err);
	}

	std::vector<MeshImpl::VertexHandle> vertices;
	mesq_mesh.get_all_vertices(vertices, err);

	for (int ino = 0; ino < vertices.size() ; ino++){
		Mesh::VertexHandle vertex = vertices[ino];
		MsqVertex aux;
		mesq_mesh.vertices_get_coordinates( &vertex, &aux, 1, err );
		coords[3*ino+0] = aux[0];
		coords[3*ino+1] = aux[1];
		coords[3*ino+2] = aux[2];
	}

	free(fixed_nodes);
	free(gid);
	free(mid);
	free(pid);

}

void MeshOptimization(hexa_tree_t* mesh, std::vector<double>& coords, std::vector<int> material_fixed_nodes){

	//TODO list
	// create a global id for every id (element, surface, edge, node)
	// generate the mesh with the ghost elements
	// run in serial and local the optimization
	// run again in paralell

	//hash of the fixed nodes
	bool clamped = true;
	sc_hash_array_t* hash_FixedNodes = sc_hash_array_new(sizeof(node_t), edge_hash_fn, edge_equal_fn, &clamped);

	for(int ino = 0; ino < material_fixed_nodes.size(); ino++){
		size_t position;
		node_t *r;
		node_t key;
		octant_node_t* node = (octant_node_t*) sc_array_index (&mesh->nodes, material_fixed_nodes[ino]);
		key.coord[0] = coords[3*node->id+0];
		key.coord[1] = coords[3*node->id+1];
		key.coord[2] = coords[3*node->id+2];
		key.node_id = node->id;
		r = (node_t*) sc_hash_array_insert_unique(hash_FixedNodes, &key, &position);
		if (r != NULL) {
			r->coord[0] = key.coord[0];
			r->coord[1] = key.coord[1];
			r->coord[2] = key.coord[2];
			r->node_id = key.node_id;
		} else {

		}
	}

	////////////////////////////////
	////DEBUG
	////////////////////////////////
	for(int ino = 0; ino < mesh->nodes.elem_count; ino ++){
		octant_node_t* elem = (octant_node_t*) sc_array_index (&mesh->nodes, ino);
		mesh->part_nodes[ino] = -1;
	}
	////////////////////////////////
	////DEBUG
	////////////////////////////////
	printf("     Line Optimization\n");
	OptLine(mesh, coords, hash_FixedNodes);

	//printf("     Surface Optimization\n");
	//OptSurface(mesh, coords, hash_FixedNodes);

	//TODO error here... this procedure "kill" the topo...
	//printf("     Redo coords\n");
	//RedoCoords(mesh, coords);

	printf("     Volume Optimization\n");
	OptVolume(mesh, coords, hash_FixedNodes);
}
