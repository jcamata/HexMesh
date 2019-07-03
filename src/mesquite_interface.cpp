/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

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

void SurfaceOptmimizationMesquite(hexa_tree_t* mesh, std::vector<double>& coords, sc_hash_array_t* hash_FixedNodes,int *snodes,int axis){
	std::vector<int> elem_aux;
	std::vector<int> aux;
	int nelem=0;
	int nvertices=0;
	size_t position;
	node_t *r;
	node_t key;
	Mesquite::MsqPrintError err(std::cout);
	sc_hash_array_truncate(hash_FixedNodes);


	// TODO essa parte deve ser trocada por uma estrutura que tenha a malha da superficie externa
	// que sera usada mais tarde para extrudar CC tipo PML ou similares...
	bool clamped = true;
	sc_hash_array_t* hash_SurfaceNodes = sc_hash_array_new(sizeof(node_t), edge_hash_fn, edge_equal_fn, &clamped);

	if(axis==0){
		for(int  i =0; i < mesh->elements.elem_count; i++){
			octant_t * elem = (octant_t*) sc_array_index(&mesh->elements, i);
			if(elem->x == (mesh->ncellx-1)){
				elem_aux.push_back(elem->id);
				for(int j = 0; j < 4; j++){
					key.coord[0] = coords[3*elem->nodes[snodes[j]].id+0];
					key.coord[1] = coords[3*elem->nodes[snodes[j]].id+1];
					key.coord[2] = coords[3*elem->nodes[snodes[j]].id+2];
					key.node_id = elem->nodes[snodes[j]].id;
					r = (node_t*) sc_hash_array_insert_unique(hash_SurfaceNodes, &key, &position);
					if (r != NULL) {
						r->coord[0] = key.coord[0];
						r->coord[1] = key.coord[1];
						r->coord[2] = key.coord[2];
						r->node_id = key.node_id;
					} else {
					}
				}
			}
		}
	}else if(axis==1){
		for(int  i =0; i < mesh->elements.elem_count; i++){
			octant_t * elem = (octant_t*) sc_array_index(&mesh->elements, i);
			if(elem->x == 0){
				elem_aux.push_back(elem->id);
				for(int j = 0; j < 4; j++){
					key.coord[0] = coords[3*elem->nodes[snodes[j]].id+0];
					key.coord[1] = coords[3*elem->nodes[snodes[j]].id+1];
					key.coord[2] = coords[3*elem->nodes[snodes[j]].id+2];
					key.node_id = elem->nodes[snodes[j]].id;
					r = (node_t*) sc_hash_array_insert_unique(hash_SurfaceNodes, &key, &position);
					if (r != NULL) {
						r->coord[0] = key.coord[0];
						r->coord[1] = key.coord[1];
						r->coord[2] = key.coord[2];
						r->node_id = key.node_id;
					} else {
					}
				}
			}
		}
	}else if(axis==2){
		for(int  i =0; i < mesh->elements.elem_count; i++){
			octant_t * elem = (octant_t*) sc_array_index(&mesh->elements, i);
			if(elem->y == (mesh->ncelly-1)){
				elem_aux.push_back(elem->id);
				for(int j = 0; j < 4; j++){
					key.coord[0] = coords[3*elem->nodes[snodes[j]].id+0];
					key.coord[1] = coords[3*elem->nodes[snodes[j]].id+1];
					key.coord[2] = coords[3*elem->nodes[snodes[j]].id+2];
					key.node_id = elem->nodes[snodes[j]].id;
					r = (node_t*) sc_hash_array_insert_unique(hash_SurfaceNodes, &key, &position);
					if (r != NULL) {
						r->coord[0] = key.coord[0];
						r->coord[1] = key.coord[1];
						r->coord[2] = key.coord[2];
						r->node_id = key.node_id;
					} else {
					}
				}
			}
		}
	}else if(axis==3){

		for(int  i =0; i < mesh->elements.elem_count; i++){
			octant_t * elem = (octant_t*) sc_array_index(&mesh->elements, i);
			if(elem->y == 0){
				elem_aux.push_back(elem->id);
				for(int j = 0; j < 4; j++){
					key.coord[0] = coords[3*elem->nodes[snodes[j]].id+0];
					key.coord[1] = coords[3*elem->nodes[snodes[j]].id+1];
					key.coord[2] = coords[3*elem->nodes[snodes[j]].id+2];
					key.node_id = elem->nodes[snodes[j]].id;
					r = (node_t*) sc_hash_array_insert_unique(hash_SurfaceNodes, &key, &position);
					if (r != NULL) {
						r->coord[0] = key.coord[0];
						r->coord[1] = key.coord[1];
						r->coord[2] = key.coord[2];
						r->node_id = key.node_id;
					} else {
					}
				}
			}
		}
	}

	nelem = elem_aux.size();
	nvertices = hash_SurfaceNodes->a.elem_count;
	printf("Malha de superficie normal a %d tem: %d elementos e %d vertices\n",axis, nelem,nvertices);

	std::vector<int> conn(4*nelem);
	std::vector<int> conn_aux(4*nelem);
	std::vector<double> Scoors(3*nvertices);
	bool   *fixed_nodes = (bool*)malloc(nvertices*sizeof(bool));

	// montando a malha para passar para o mesquite...
	// conectividade na malha "global"
	int c;
	c = 0;
	for(int  i =0; i < elem_aux.size(); i++){
		octant_t * elem = (octant_t*) sc_array_index(&mesh->elements, elem_aux[i]);
		//printf("Elemento:%d. Id do no: ",elem->id);
		for(int j = 0; j < 4; j++){
			//printf("%d ",elem->nodes[snodes[j]].id);
			conn_aux[c] = elem->nodes[snodes[j]].id; c++;
		}
		//printf("\n");
	}

	if(true){
		//colocando os nos das bordas na lista de nos bloqueados...
		for(int  i =0; i < elem_aux.size(); i++){
			octant_t * elem = (octant_t*) sc_array_index(&mesh->elements, elem_aux[i]);
			//printf("Id do elemento:%d\n Id dos nos que ele olha:",elem->id);
			for(int j = 0; j<4;j++){
				int n_node =  snodes[j];

				/*printf(" %d",elem->nodes[n_node].id);

				 */
				if(axis==0 || axis ==1){
					if( elem->nodes[n_node].y == 0 || elem->nodes[n_node].z == 0 || elem->nodes[n_node].y == mesh->ncelly || elem->nodes[n_node].z == (mesh->ncellz) ){
						key.coord[0] = coords[3*elem->nodes[n_node].id+0];
						key.coord[1] = coords[3*elem->nodes[n_node].id+1];
						key.coord[2] = coords[3*elem->nodes[n_node].id+2];
						key.node_id = elem->nodes[n_node].id;
						r = (node_t*) sc_hash_array_insert_unique(hash_FixedNodes, &key, &position);
						if (r != NULL) {
							r->coord[0] = key.coord[0];
							r->coord[1] = key.coord[1];
							r->coord[2] = key.coord[2];
							r->node_id = key.node_id;
						} else {

						}
					}
				}else if(axis==2 || axis ==3){
					if( elem->nodes[n_node].x == 0 || elem->nodes[n_node].z == 0 || elem->nodes[n_node].x == mesh->ncellx || elem->nodes[n_node].z == (mesh->ncellz) ){
						key.coord[0] = coords[3*elem->nodes[n_node].id+0];
						key.coord[1] = coords[3*elem->nodes[n_node].id+1];
						key.coord[2] = coords[3*elem->nodes[n_node].id+2];
						key.node_id = elem->nodes[n_node].id;
						r = (node_t*) sc_hash_array_insert_unique(hash_FixedNodes, &key, &position);
						if (r != NULL) {
							r->coord[0] = key.coord[0];
							r->coord[1] = key.coord[1];
							r->coord[2] = key.coord[2];
							r->node_id = key.node_id;
						} else {

						}
					}
				}
			}
			//printf("\n");
		}
	}

	printf("Numero de nos fixos na face %d é de %d\n",axis,hash_FixedNodes->a.elem_count);
	printf("nos dentro da hash fixa:");
	for(int j = 0; j<hash_FixedNodes->a.elem_count;j++){
		node_t* rr= (node_t*) sc_array_index (&hash_FixedNodes->a, j);
		printf(" %d", rr->node_id);
	}

	//printf("\n nos fixos: ");
	//nos fixos e coordenadas dos nos
	for(int i = 0 ; i<nvertices ;i++ ){
		fixed_nodes[i] = false;
		node_t * node = (node_t*) sc_array_index(&hash_SurfaceNodes->a, i);
		Scoors[3*i+0] = coords[3*node->node_id+0];
		Scoors[3*i+1] = coords[3*node->node_id+1];
		Scoors[3*i+2] = coords[3*node->node_id+2];
		aux.push_back(node->node_id);
		for(int j = 0; j<hash_FixedNodes->a.elem_count;j++){
			node_t* rr= (node_t*) sc_array_index (&hash_FixedNodes->a, j);
			if(rr->node_id==node->node_id){
				fixed_nodes[i] = true;
			}

		}
		//bool tre = sc_hash_array_lookup(hash_FixedNodes,&node,&position);
		//if(tre){
		//	fixed_nodes[i] = true;
		//}
		//printf(" %d", fixed_nodes[i]);
	}
	//printf("\n");

	//refazendo a conectividade com nos de 0 a nvertices...
	//passando uma malha "local"
	for(int i = 0; i <aux.size();i++){
		for(int j = 0; j <conn_aux.size();j++){
			if(conn_aux[j]==aux[i]){
				conn[j] = i;
			}
		}
	}

	//for(int j = 0; j <conn_aux.size();j++){
	//printf("aux:%d novo:%d\n",conn_aux[j], conn[j]);
	//}

	//finalmente mandando a malha pro mesquite...
	Mesquite::MeshImpl mesq_mesh0(nvertices,nelem,Mesquite::QUADRILATERAL, &fixed_nodes[0], &Scoors[0], &conn[0]);

	Mesquite::PlanarDomain plane;
	if(axis==0){
		Mesquite::PlanarDomain plane(Vector3D(1,0,0), (Scoors[0],Scoors[1],Scoors[2]));
	}else if(axis ==1){
		Mesquite::PlanarDomain plane(Vector3D(1,0,0), (Scoors[0],Scoors[1],Scoors[2]));
	}else if(axis ==2){
		Mesquite::PlanarDomain plane(Vector3D(0,1,0), (Scoors[0],Scoors[1],Scoors[2]));
	}else if(axis ==3){
		Mesquite::PlanarDomain plane(Vector3D(0,1,0), (Scoors[0],Scoors[1],Scoors[2]));
	}


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

	//**************Set stopping criterion****************
	TerminationCriterion sc2;
	sc2.add_iteration_limit( 10 );
	lapl1.set_outer_termination_criterion(&sc2);

	// adds 1 pass of pass1 to mesh_set1
	queue1.add_quality_assessor(&stop_qa,err);
	queue1.set_master_quality_improver(&lapl1, err);
	queue1.add_quality_assessor(&stop_qa,err);
	// adds 1 passes of pass2 to mesh_set1
	//  mesh_set1.add_quality_pass(pass2);


	// launches optimization on mesh_set1
	MeshDomainAssoc mesh_and_domain0 = MeshDomainAssoc(&mesq_mesh0, &plane);
	//printf("Status da criacao do mesh_domain para normal %d: %d\n",axis, mesh_and_domain0.are_compatible());

	Timer t;
	queue1.run_instructions(&mesh_and_domain0, err);
	double secs = t.since_birth();

	assert(mesh_and_domain0.are_compatible()==1);

	/*
	UntangleWrapper::UntangleMetric metric0 = UntangleWrapper::BETA;
	UntangleWrapper un_wrapper0 (metric0);
	un_wrapper0.set_vertex_movement_limit_factor( 0.005 );
	un_wrapper0.set_outer_iteration_limit(10);
	un_wrapper0.run_instructions( &mesh_and_domain0, err );
	 */

	/*
	//LaplaceWrapper lp_wrapper0;
	//lp_wrapper0.set_vertex_movement_limit_factor(1.e-5);
	//lp_wrapper0.set_iteration_limit(10);
	//lp_wrapper0.enable_culling(false);
	//lp_wrapper0.run_instructions( &mesh_and_domain0, err );
	 */

	/*
	ShapeImprover smoother;
	IdealWeightInverseMeanRatio extra_metric;
	smoother.quality_assessor().add_quality_assessment(&extra_metric);
	smoother.set_cpu_time_limit(120);
	smoother.run_instructions( &mesh_and_domain0, err );
	 */

	/*
	 * SizeAdaptShapeWrapper smoother(1e-2);
  MeshDomainAssoc mesh_and_domain = MeshDomainAssoc(&mesh, &geom);
  smoother.run_instructions( &mesh_and_domain, err);
	 */

	std::vector<MeshImpl::VertexHandle> vertices;
	mesh_and_domain0.get_mesh()->get_all_vertices(vertices,err);

	for (int ino = 0; ino < vertices.size() ; ino++){
		Mesh::VertexHandle vertex = vertices[ino];
		MsqVertex aux_msq;
		mesh_and_domain0.get_mesh()->vertices_get_coordinates( &vertex, &aux_msq, 1, err );
		node_t * node = (node_t*) sc_array_index(&hash_SurfaceNodes->a, ino);
		coords[3*node->node_id+0] = aux_msq[0];
		coords[3*node->node_id+1] = aux_msq[1];
		coords[3*node->node_id+2] = aux_msq[2];
	}

	aux.clear();
	elem_aux.clear();
	mesq_mesh0.clear();
	conn.clear();
	conn_aux.clear();
	Scoors.clear();
	free(fixed_nodes);
	sc_hash_array_truncate(hash_SurfaceNodes);

}

void OptLine(hexa_tree_t* mesh, std::vector<double>& coords, std::vector<int> material_fixed_nodes){

	//criando a hash para os nos que eu nao posso mover...
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

	double xs;
	double xe;
	double ys;
	double ye;
	double zs;
	double ze;

	for(int ino = 0; ino < mesh->nodes.elem_count; ino++){
		octant_node_t* node = (octant_node_t*) sc_array_index (&mesh->nodes, ino);
		if(node->z == 0){zs=coords[3*node->id+2];}
		if(node->z == (mesh->ncellz))  {ze=coords[3*node->id+2];}
		if(node->y == mesh->y_start){ys=coords[3*node->id+1];}
		if(node->y == mesh->y_end)  {ye=coords[3*node->id+1];}
		if(node->x == mesh->x_start)  {xs=coords[3*node->id];}
		if(node->x == mesh->x_end)    {xe=coords[3*node->id];}
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

	for(int ino = 0; ino < mesh->nodes.elem_count; ino++){
		octant_node_t* node = (octant_node_t*) sc_array_index (&mesh->nodes, ino);

		double x = coords[3*node->id + 0];
		double y = coords[3*node->id + 1];
		double z = coords[3*node->id + 2];

		//line 0
		if(node->z == 0 && y == ys){
			aux0.push_back(node->id);
		}
		//line 1
		if(node->z == 0 && x == xe){
			aux1.push_back(node->id);
		}
		//line 2
		if(node->z == 0 && y == ye){
			aux2.push_back(node->id);
		}
		//line 3
		if(node->z == 0 && x == xs){
			aux3.push_back(node->id);
		}

		//line 4
		if(x == xs && y == ys){
			aux4.push_back(node->id);
		}
		//line 5
		if(x == xe && y == ys){
			aux5.push_back(node->id);
		}
		//line 6
		if(x == xe && y == ye){
			aux6.push_back(node->id);
		}
		//line 7
		if(x == xs && y == ye){
			aux7.push_back(node->id);
		}

		//line 8
		if(z == ze && y == ys){
			aux8.push_back(node->id);
		}
		//line 9
		if(z == ze && x == xe){
			aux9.push_back(node->id);
		}
		//line 10
		if(z == ze && y == ye){
			aux10.push_back(node->id);
		}
		//line 11
		if(z == ze && x == xs){
			aux11.push_back(node->id);
		}


		if(false){
			//line 0
			if(node->z == 0 && node->y == mesh->y_start){
				aux0.push_back(node->id);
			}
			//line 1
			if(node->z == 0 && node->x == mesh->x_end){
				aux1.push_back(node->id);
			}
			//line 2
			if(node->z == 0 && node->y == mesh->y_end){
				aux2.push_back(node->id);
			}
			//line 3
			if(node->z == 0 && node->x == mesh->x_start){
				aux3.push_back(node->id);
			}

			//line 4
			if(node->x == mesh->x_start && node->y == mesh->y_start){
				aux4.push_back(node->id);
			}
			//line 5
			if(node->x == mesh->x_end && node->y == mesh->y_start){
				aux5.push_back(node->id);
			}
			//line 6
			if(node->x == mesh->x_end && node->y == mesh->y_end){
				aux6.push_back(node->id);
			}
			//line 7
			if(node->x == mesh->x_start && node->y == mesh->y_end){
				aux7.push_back(node->id);
			}

			//line 8
			if(node->z == (mesh->ncellz) && node->y == mesh->y_start){
				aux8.push_back(node->id);
			}
			//line 9
			if(node->z == (mesh->ncellz) && node->x == mesh->x_end){
				aux9.push_back(node->id);
			}
			//line 10
			if(node->z == (mesh->ncellz) && node->y == mesh->y_end){
				aux10.push_back(node->id);
			}
			//line 11
			if(node->z == (mesh->ncellz) && node->x == mesh->x_start){
				aux11.push_back(node->id);
			}
		}
	}

	for(int iedge = 0; iedge<12; iedge++){
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

		//printf("edge %d \n",iedge);
		//for(int ino = 0; ino < aux.size(); ino++){
		//	printf("%d ",aux[ino]);
		//	}
		//printf("\n");

		int niter = 0;
		double tol = 1;
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
					coords[3*node+0] = (coords[3*node0+0] + coords[3*node2+0])/2;
					coords[3*node+1] = (coords[3*node0+1] + coords[3*node2+1])/2;
					coords[3*node+2] = (coords[3*node0+2] + coords[3*node2+2])/2;
					double a1 = coords[3*node+0];
					double b1 = coords[3*node+1];
					double c1 = coords[3*node+2];
					tol = sqrt((a-a1)*(a-a1) + (b-b1)*(b-b1) + (c-c1)*(c-c1));
					//printf("%f\n",tol);
				}
			}
			niter++;
		}
	}
	printf("     Please check 1D optimization\n");

}

void OptSurface(hexa_tree_t* mesh, std::vector<double>& coords, std::vector<int> material_fixed_nodes){

	Mesquite::MsqPrintError err(std::cout);

	//criando a hash para os nos que eu nao posso mover...
	bool clamped = true;
	sc_hash_array_t* hash_FixedNodes = sc_hash_array_new(sizeof(node_t), edge_hash_fn, edge_equal_fn, &clamped);
	size_t position;
	node_t *r;
	node_t key;

	for(int n = 0;n<material_fixed_nodes.size();n++){
		octant_node_t* node = (octant_node_t*) sc_array_index (&mesh->nodes, material_fixed_nodes[n]);
		key.coord[0] = coords[3*node->id+0];
		key.coord[1] = coords[3*node->id+1];
		key.coord[2] = coords[3*node->id+2];
		key.node_id = node->id;
		//printf("no fixo, id dele é:%d\n",material_fixed_nodes[n]);
		r = (node_t*) sc_hash_array_insert_unique(hash_FixedNodes, &key, &position);
		if (r != NULL) {
			r->coord[0] = key.coord[0];
			r->coord[1] = key.coord[1];
			r->coord[2] = key.coord[2];
			r->node_id = key.node_id;
		} else {

		}
	}

	double xs;
	double xe;
	double ys;
	double ye;
	double zs;
	double ze;

	for(int n = 0;n<mesh->nodes.elem_count;n++){
		octant_node_t* node = (octant_node_t*) sc_array_index (&mesh->nodes, n);
		if(node->z == 0){zs=coords[3*node->id+2];}
		if(node->z == (mesh->ncellz-1))  {ze=coords[3*node->id+2];}
		if(node->y == mesh->y_start){ys=coords[3*node->id+1];}
		if(node->y == mesh->y_end)  {ye=coords[3*node->id+1];}
		if(node->x == mesh->x_start)  {xs=coords[3*node->id];}
		if(node->x == mesh->x_end)    {xe=coords[3*node->id];}
	}

	for(int treta = 0; treta <1 ; treta ++){

		//x+
		int snodes[4];
		int axis;
		if(treta == 0){
			//x+
			snodes[0]=5;
			snodes[1]=6;
			snodes[2]=2;
			snodes[3]=1;
			axis = 0;
		}else if(treta == 1){
			//x-
			snodes[0]=0;
			snodes[1]=3;
			snodes[2]=7;
			snodes[3]=4;
			axis = 1;
		}else if(treta ==2){
			//y+
			snodes[0]=7;
			snodes[1]=6;
			snodes[2]=2;
			snodes[3]=3;
			axis = 2;
		}else if(treta==3){
			//y-
			snodes[0]=4;
			snodes[1]=5;
			snodes[2]=1;
			snodes[3]=0;
			axis = 3;
		}else if(treta==4){
			//z+
			snodes[0]=0;
			snodes[1]=1;
			snodes[2]=2;
			snodes[3]=3;
			axis = 4;
		}else if(treta==5){
			//z-
			snodes[0]=4;
			snodes[1]=5;
			snodes[2]=6;
			snodes[3]=7;
			axis = 5;
		}

		std::vector<int> elem_aux;
		std::vector<int> aux;
		int nelem=0;
		int nvertices=0;

		// essa parte deve ser trocada por uma estrutura que tenha a malha da superficie externa
		// que sera usada mais tarde para extrudar CC tipo PML ou similares...
		//bool clamped = true;
		sc_hash_array_t* hash_SurfaceNodes = sc_hash_array_new(sizeof(node_t), edge_hash_fn, edge_equal_fn, &clamped);

		if(axis==0){
			for(int  i =0; i < mesh->elements.elem_count; i++){
				octant_t * elem = (octant_t*) sc_array_index(&mesh->elements, i);
				if(elem->x == (mesh->ncellx-1)){
					elem_aux.push_back(elem->id);
					for(int j = 0; j < 4; j++){
						key.coord[0] = coords[3*elem->nodes[snodes[j]].id+0];
						key.coord[1] = coords[3*elem->nodes[snodes[j]].id+1];
						key.coord[2] = coords[3*elem->nodes[snodes[j]].id+2];
						key.node_id = elem->nodes[snodes[j]].id;
						r = (node_t*) sc_hash_array_insert_unique(hash_SurfaceNodes, &key, &position);
						if (r != NULL) {
							r->coord[0] = key.coord[0];
							r->coord[1] = key.coord[1];
							r->coord[2] = key.coord[2];
							r->node_id = key.node_id;
						} else {
						}
					}
				}
			}
		}else if(axis==1){
			for(int  i =0; i < mesh->elements.elem_count; i++){
				octant_t * elem = (octant_t*) sc_array_index(&mesh->elements, i);
				if(elem->x == 0){
					elem_aux.push_back(elem->id);
					for(int j = 0; j < 4; j++){
						key.coord[0] = coords[3*elem->nodes[snodes[j]].id+0];
						key.coord[1] = coords[3*elem->nodes[snodes[j]].id+1];
						key.coord[2] = coords[3*elem->nodes[snodes[j]].id+2];
						key.node_id = elem->nodes[snodes[j]].id;
						r = (node_t*) sc_hash_array_insert_unique(hash_SurfaceNodes, &key, &position);
						if (r != NULL) {
							r->coord[0] = key.coord[0];
							r->coord[1] = key.coord[1];
							r->coord[2] = key.coord[2];
							r->node_id = key.node_id;
						} else {
						}
					}
				}
			}
		}else if(axis==2){
			for(int  i =0; i < mesh->elements.elem_count; i++){
				octant_t * elem = (octant_t*) sc_array_index(&mesh->elements, i);
				if(elem->y == (mesh->ncelly-1)){
					elem_aux.push_back(elem->id);
					for(int j = 0; j < 4; j++){
						key.coord[0] = coords[3*elem->nodes[snodes[j]].id+0];
						key.coord[1] = coords[3*elem->nodes[snodes[j]].id+1];
						key.coord[2] = coords[3*elem->nodes[snodes[j]].id+2];
						key.node_id = elem->nodes[snodes[j]].id;
						r = (node_t*) sc_hash_array_insert_unique(hash_SurfaceNodes, &key, &position);
						if (r != NULL) {
							r->coord[0] = key.coord[0];
							r->coord[1] = key.coord[1];
							r->coord[2] = key.coord[2];
							r->node_id = key.node_id;
						} else {
						}
					}
				}
			}
		}else if(axis==3){

			for(int  i =0; i < mesh->elements.elem_count; i++){
				octant_t * elem = (octant_t*) sc_array_index(&mesh->elements, i);
				if(elem->y == 0){
					elem_aux.push_back(elem->id);
					for(int j = 0; j < 4; j++){
						key.coord[0] = coords[3*elem->nodes[snodes[j]].id+0];
						key.coord[1] = coords[3*elem->nodes[snodes[j]].id+1];
						key.coord[2] = coords[3*elem->nodes[snodes[j]].id+2];
						key.node_id = elem->nodes[snodes[j]].id;
						r = (node_t*) sc_hash_array_insert_unique(hash_SurfaceNodes, &key, &position);
						if (r != NULL) {
							r->coord[0] = key.coord[0];
							r->coord[1] = key.coord[1];
							r->coord[2] = key.coord[2];
							r->node_id = key.node_id;
						} else {
						}
					}
				}
			}
		}else if(axis==4){

			for(int  i =0; i < mesh->elements.elem_count; i++){
				octant_t * elem = (octant_t*) sc_array_index(&mesh->elements, i);
				if(elem->z == 0){
					elem_aux.push_back(elem->id);
					for(int j = 0; j < 4; j++){
						key.coord[0] = coords[3*elem->nodes[snodes[j]].id+0];
						key.coord[1] = coords[3*elem->nodes[snodes[j]].id+1];
						key.coord[2] = coords[3*elem->nodes[snodes[j]].id+2];
						key.node_id = elem->nodes[snodes[j]].id;
						r = (node_t*) sc_hash_array_insert_unique(hash_SurfaceNodes, &key, &position);
						if (r != NULL) {
							r->coord[0] = key.coord[0];
							r->coord[1] = key.coord[1];
							r->coord[2] = key.coord[2];
							r->node_id = key.node_id;
						} else {
						}
					}
				}
			}
		}else if(axis==5){

			int nz_max;
			if(mesh->max_levels<4){
				nz_max = mesh->ncellz-1;
			}else{
				nz_max = mesh->max_z-3*(mesh->max_levels-3);
			}

			for(int  i =0; i < mesh->elements.elem_count; i++){
				octant_t * elem = (octant_t*) sc_array_index(&mesh->elements, i);
				if(elem->z == nz_max){
					elem_aux.push_back(elem->id);
					for(int j = 0; j < 4; j++){
						key.coord[0] = coords[3*elem->nodes[snodes[j]].id+0];
						key.coord[1] = coords[3*elem->nodes[snodes[j]].id+1];
						key.coord[2] = coords[3*elem->nodes[snodes[j]].id+2];
						key.node_id = elem->nodes[snodes[j]].id;
						r = (node_t*) sc_hash_array_insert_unique(hash_SurfaceNodes, &key, &position);
						if (r != NULL) {
							r->coord[0] = key.coord[0];
							r->coord[1] = key.coord[1];
							r->coord[2] = key.coord[2];
							r->node_id = key.node_id;
						} else {
						}
					}
				}
			}
		}

		nelem = elem_aux.size();
		nvertices = hash_SurfaceNodes->a.elem_count;
		printf("Malha de superficie normal a %d tem: %d elementos e %d vertices\n",axis, nelem,nvertices);

		std::vector<int> conn(4*nelem);
		std::vector<int> conn_aux(4*nelem);
		std::vector<double> Scoors(3*nvertices);
		bool   *fixed_nodes = (bool*)malloc(nvertices*sizeof(bool));

		// montando a malha para passar para o mesquite...
		// conectividade na malha "global"
		int c;
		c = 0;
		for(int  i =0; i < elem_aux.size(); i++){
			octant_t * elem = (octant_t*) sc_array_index(&mesh->elements, elem_aux[i]);
			//printf("Elemento:%d. Id do no: ",elem->id);
			for(int j = 0; j < 4; j++){
				//printf("%d ",elem->nodes[snodes[j]].id);
				conn_aux[c] = elem->nodes[snodes[j]].id; c++;
			}
			//printf("\n");
		}

		if(true){
			//colocando os nos das bordas na lista de nos bloqueados...
			for(int  i =0; i < elem_aux.size(); i++){
				octant_t * elem = (octant_t*) sc_array_index(&mesh->elements, elem_aux[i]);
				//printf("Id do elemento:%d\n Id dos nos que ele olha:",elem->id);
				for(int j = 0; j<4;j++){
					int n_node =  snodes[j];

					/*printf(" %d",elem->nodes[n_node].id);

					 */
					if(axis==0 || axis ==1){
						if( elem->nodes[n_node].y == 0 || elem->nodes[n_node].z == 0 || elem->nodes[n_node].y == mesh->ncelly || elem->nodes[n_node].z == (mesh->ncellz) ){
							key.coord[0] = coords[3*elem->nodes[n_node].id+0];
							key.coord[1] = coords[3*elem->nodes[n_node].id+1];
							key.coord[2] = coords[3*elem->nodes[n_node].id+2];
							key.node_id = elem->nodes[n_node].id;
							r = (node_t*) sc_hash_array_insert_unique(hash_FixedNodes, &key, &position);
							if (r != NULL) {
								r->coord[0] = key.coord[0];
								r->coord[1] = key.coord[1];
								r->coord[2] = key.coord[2];
								r->node_id = key.node_id;
							} else {

							}
						}
					}else if(axis==2 || axis ==3){
						if( elem->nodes[n_node].x == 0 || elem->nodes[n_node].z == 0 || elem->nodes[n_node].x == mesh->ncellx || elem->nodes[n_node].z == (mesh->ncellz) ){
							key.coord[0] = coords[3*elem->nodes[n_node].id+0];
							key.coord[1] = coords[3*elem->nodes[n_node].id+1];
							key.coord[2] = coords[3*elem->nodes[n_node].id+2];
							key.node_id = elem->nodes[n_node].id;
							r = (node_t*) sc_hash_array_insert_unique(hash_FixedNodes, &key, &position);
							if (r != NULL) {
								r->coord[0] = key.coord[0];
								r->coord[1] = key.coord[1];
								r->coord[2] = key.coord[2];
								r->node_id = key.node_id;
							} else {

							}
						}
					}else if(axis==4 || axis ==5){
						if( elem->nodes[n_node].x == 0 || elem->nodes[n_node].y == 0 || elem->nodes[n_node].x == mesh->ncellx || elem->nodes[n_node].y == mesh->ncelly
								|| coords[3*elem->nodes[n_node].id+2] > 0){
							key.coord[0] = coords[3*elem->nodes[n_node].id+0];
							key.coord[1] = coords[3*elem->nodes[n_node].id+1];
							key.coord[2] = coords[3*elem->nodes[n_node].id+2];
							key.node_id = elem->nodes[n_node].id;
							r = (node_t*) sc_hash_array_insert_unique(hash_FixedNodes, &key, &position);
							if (r != NULL) {
								r->coord[0] = key.coord[0];
								r->coord[1] = key.coord[1];
								r->coord[2] = key.coord[2];
								r->node_id = key.node_id;
							} else {

							}
						}
					}
				}
				//printf("\n");
			}
		}

		//printf("Numero de nos fixos na face %d é de %d\n",axis,hash_FixedNodes->a.elem_count);
		//printf("nos dentro da hash fixa:");
		for(int j = 0; j<hash_FixedNodes->a.elem_count;j++){
			node_t* rr= (node_t*) sc_array_index (&hash_FixedNodes->a, j);
			//printf(" %d", rr->node_id);
		}
		//printf("\n");
		//printf("nos fixos: ");
		//nos fixos e coordenadas dos nos
		for(int i = 0 ; i<nvertices ;i++ ){
			fixed_nodes[i] = false;
			node_t * node = (node_t*) sc_array_index(&hash_SurfaceNodes->a, i);
			Scoors[3*i+0] = coords[3*node->node_id+0];
			Scoors[3*i+1] = coords[3*node->node_id+1];
			Scoors[3*i+2] = coords[3*node->node_id+2];
			aux.push_back(node->node_id);

			for(int j = 0; j<hash_FixedNodes->a.elem_count;j++){
				node_t* rr= (node_t*) sc_array_index (&hash_FixedNodes->a, j);
				if(rr->node_id==node->node_id){
					fixed_nodes[i] = true;
				}
			}
			//printf(" %d", fixed_nodes[i]);
			/*
			bool tre = sc_hash_array_lookup(hash_FixedNodes,&node->node_id,&position);
			if(tre){
				fixed_nodes[i] = true;
			}

			 */
		}
		//printf("\n");

		//refazendo a conectividade com nos de 0 a nvertices...
		//passando uma malha "local"
		for(int i = 0; i <aux.size();i++){
			for(int j = 0; j <conn_aux.size();j++){
				if(conn_aux[j]==aux[i]){
					conn[j] = i;
				}
			}
		}

		//for(int j = 0; j <conn_aux.size();j++){
		//printf("aux:%d novo:%d\n",conn_aux[j], conn[j]);
		//}

		if(nelem!=0){
			//finalmente mandando a malha pro mesquite...
			Mesquite::MeshImpl mesq_mesh0(nvertices,nelem,Mesquite::QUADRILATERAL, &fixed_nodes[0], &Scoors[0], &conn[0]);


			Vector3D point(Scoors[0],Scoors[1],Scoors[2]);
			Vector3D normal;
			if(axis == 0){
				normal.set(1,0,0);
				//Mesquite::PlanarDomain plane(normal, point);
			}else if(axis == 1){
				normal.set(-1,0,0);
				//Mesquite::PlanarDomain plane(Vector3D(-1,0,0), point);
			}else if(axis == 2){
				normal.set(0,1,0);
				//Mesquite::PlanarDomain plane(Vector3D(0,1,0), point);
			}else if(axis == 3){
				normal.set(0,-1,0);
				//Mesquite::PlanarDomain plane(Vector3D(0,-1,0), point);
			}else if(axis == 4){
				normal.set(0,0,-1);
				//Mesquite::PlanarDomain plane(Vector3D(0,0,1), point);
			}else if(axis == 5){
				normal.set(0,0,1);
				//Mesquite::PlanarDomain plane(Vector3D(0,0,-1), point);
			}
			Mesquite::PlanarDomain plane(normal, point);
			MeshDomainAssoc mesh_and_domain0 = MeshDomainAssoc(&mesq_mesh0, &plane);
			printf("Status da criacao do mesh_domain para normal %d: %d\n",axis, mesh_and_domain0.are_compatible());
			assert(mesh_and_domain0.are_compatible()==1);

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

			UntangleWrapper::UntangleMetric metric0 = UntangleWrapper::BETA;
			UntangleWrapper un_wrapper0 (metric0);
			un_wrapper0.set_vertex_movement_limit_factor( 0.005 );
			un_wrapper0.set_outer_iteration_limit(10);
			un_wrapper0.run_instructions( &mesh_and_domain0, err );

			if (err){
				std::cout << err << std::endl;
			}

			ShapeImprover smoother;
			IdealWeightInverseMeanRatio extra_metric;
			smoother.quality_assessor().add_quality_assessment(&extra_metric);
			smoother.set_cpu_time_limit(120);
			smoother.run_instructions( &mesh_and_domain0, err );
			if (err){
				std::cout << err << std::endl;
			}

			if(false){
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

				//**************Set stopping criterion****************
				TerminationCriterion sc2;
				sc2.add_iteration_limit( 50 );
				sc2.add_relative_vertex_movement(1);
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
			LaplaceWrapper lp_wrapper0;
			lp_wrapper0.set_vertex_movement_limit_factor(1.e-6);
			lp_wrapper0.set_iteration_limit(5);
			lp_wrapper0.enable_culling(true);
			//lp_wrapper0.run_instructions( &mesh_and_domain0, err );
			printf("LAPLACE WRAPPER RETURN INVERTED ELEMENTS\n");
			if (err){
				std::cout << err << std::endl;
			}
			/*
		    ShapeImprover smoother;
		    IdealWeightInverseMeanRatio extra_metric;
		    smoother.quality_assessor().add_quality_assessment(&extra_metric);
		    smoother.set_cpu_time_limit(120);
		    smoother.run_instructions( &mesh_and_domain0, err );

			SizeAdaptShapeWrapper smoother(1e-2);
			smoother.run_instructions( &mesh_and_domain0, err);

	  	    SizeAdaptShapeWrapper smoother(1e-2);
	  	    MeshDomainAssoc mesh_and_domain = MeshDomainAssoc(&mesh, &geom);
	  	    smoother.run_instructions( &mesh_and_domain, err);

			 */
			std::vector<MeshImpl::VertexHandle> vertices;
			mesh_and_domain0.get_mesh()->get_all_vertices(vertices,err);

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

void OptVolume(hexa_tree_t* mesh, std::vector<double>& coords, std::vector<int> material_fixed_nodes){
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
	for(int  i =0; i < mesh->elements.elem_count; i++)
	{
		octant_t * elem = (octant_t*) sc_array_index(&mesh->elements, i);
		mid[i] = elem->n_mat;
		for(int j = 0; j < 8; j++)
		{
			tmp[j] = elem->nodes[j].id;
			//c++;
		}
		conn[c] = tmp[4]; c++;
		conn[c] = tmp[5]; c++;
		conn[c] = tmp[6]; c++;
		conn[c] = tmp[7]; c++;
		conn[c] = tmp[0]; c++;
		conn[c] = tmp[1]; c++;
		conn[c] = tmp[2]; c++;
		conn[c] = tmp[3]; c++;
	}

	assert(mesh->nodes.elem_count == mesh->local_n_nodes);

	int bnd_nodes = 0;
	for(int i = 0; i < mesh->nodes.elem_count; i++)
	{
		fixed_nodes[i] = false;
		gid[i] = (size_t) mesh->global_id[i];
		pid[i] = mesh->part_nodes[i];
		octant_node_t * node = (octant_node_t*) sc_array_index(&mesh->nodes, i);
		assert(node != NULL);

		int node_z = mesh->ncellz -1;
		if(mesh->max_levels>3){
			node_z = mesh->ncellz-2;
		}

		if(node->x == 0 || node->y == 0 || node->z == 0 ||
				node->x == mesh->ncellx || node->y == mesh->ncelly || node->z == node_z )
			//if(node->z == 0 || node->z == (mesh->ncellz-1) )
		{
			fixed_nodes[i] = true;
			bnd_nodes++;
		}

		//if(node->z == mesh->ncellz)
		//std::cout<<"Nodes: " << node->id << std::endl;
	}

	//std::cout<<"BND Nodes: " << bnd_nodes << std::endl;

	for(int i = 0; i < material_fixed_nodes.size(); i++)
	{
		int id = material_fixed_nodes[i];
		fixed_nodes[id] = true;
	}

	Mesquite::MeshImpl mesq_mesh(nvertices,nelem,Mesquite::HEXAHEDRON, &fixed_nodes[0], &coords[0], &conn[0]);

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

	UntangleWrapper::UntangleMetric metric0 = UntangleWrapper::BETA;
	UntangleWrapper un_wrapper0 (metric0);
	un_wrapper0.set_vertex_movement_limit_factor( 0.005 );
	un_wrapper0.set_outer_iteration_limit(10);
	un_wrapper0.run_instructions( &mesq_mesh, err );

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

}

void UntagleMesh(hexa_tree_t* mesh, std::vector<double>& coords, std::vector<int> material_fixed_nodes){

	OptLine(mesh,coords,material_fixed_nodes);
	OptSurface(mesh, coords,material_fixed_nodes);
	//printf("Please check here later... Volumetric UntagleMesh\n");
	//OptVolume(mesh, coords,material_fixed_nodes);

}

void MeshOpt(hexa_tree_t* mesh, std::vector<double>& coords, std::vector<int> material_fixed_nodes){
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
	for(int  i =0; i < mesh->elements.elem_count; i++)
	{
		octant_t * elem = (octant_t*) sc_array_index(&mesh->elements, i);
		mid[i] = elem->n_mat;
		for(int j = 0; j < 8; j++)
		{
			tmp[j] = elem->nodes[j].id;
			//c++;
		}
		conn[c] = tmp[4]; c++;
		conn[c] = tmp[5]; c++;
		conn[c] = tmp[6]; c++;
		conn[c] = tmp[7]; c++;
		conn[c] = tmp[0]; c++;
		conn[c] = tmp[1]; c++;
		conn[c] = tmp[2]; c++;
		conn[c] = tmp[3]; c++;
	}

	assert(mesh->nodes.elem_count == mesh->local_n_nodes);


	int node_z = mesh->ncellz -1;
	if(mesh->max_levels>3){
		node_z = mesh->ncellz-2;
	}

	int bnd_nodes = 0;
	for(int i = 0; i < mesh->nodes.elem_count; i++)
	{
		fixed_nodes[i] = false;
		gid[i] = (size_t) mesh->global_id[i];
		pid[i] = mesh->part_nodes[i];
		octant_node_t * node = (octant_node_t*) sc_array_index(&mesh->nodes, i);
		assert(node != NULL);

		if(node->x == 0 || node->y == 0 || node->z == 0 ||
				node->x == mesh->ncellx || node->y == mesh->ncelly || node->z == node_z )
			//if(node->z == 0 || node->z == (mesh->ncellz-1) )
		{
			fixed_nodes[i] = true;
			bnd_nodes++;
		}

		//if(node->z == mesh->ncellz)
		//std::cout<<"Nodes: " << node->id << std::endl;
	}

	//std::cout<<"BND Nodes: " << bnd_nodes << std::endl;

	for(int i = 0; i < material_fixed_nodes.size(); i++)
	{
		int id = material_fixed_nodes[i];
		fixed_nodes[id] = true;
	}


	Mesquite::MeshImpl mesq_mesh(nvertices,nelem,Mesquite::HEXAHEDRON, &fixed_nodes[0], &coords[0], &conn[0]);

	int default_gid = -1;
	int default_pid = 0;
	int defaut_mat = 0;
	mesq_mesh.tag_create("GLOBAL_ID"   ,Mesquite::Mesh::HANDLE,1,&default_gid, err);
	mesq_mesh.tag_create("PROCESSOR_ID",Mesquite::Mesh::INT   ,1,&default_pid, err);
	mesq_mesh.tag_create("MAT_ID",Mesquite::Mesh::INT   ,1,&defaut_mat, err);

	TagHandle tag_processor_id = mesq_mesh.tag_get("PROCESSOR_ID", err);
	TagHandle tag_global_id    = mesq_mesh.tag_get("GLOBAL_ID", err);
	TagHandle tag_mat_id    = mesq_mesh.tag_get("MAT_ID", err);
	std::vector<MeshImpl::VertexHandle> vertices;
	std::vector<MeshImpl::ElementHandle>	elem_array;

	mesq_mesh.get_all_vertices(vertices, err);
	mesq_mesh.get_all_elements(elem_array,err);
	mesq_mesh.tag_set_vertex_data(tag_global_id   ,vertices.size(),&vertices[0], gid, err);
	mesq_mesh.tag_set_vertex_data(tag_processor_id,vertices.size(),&vertices[0], pid, err);
	mesq_mesh.tag_set_element_data(tag_mat_id,nelem,&elem_array[0],mid,err);

	{
		std::ostringstream out_name;
		out_name << "original_mesh." << mesh->mpi_size << "." << mesh->mpi_rank << ".vtk";
		mesq_mesh.write_vtk(out_name.str().c_str(), err);
	}

	// do Laplacian smooth
	Mesquite::LaplaceWrapper optimizer;
	optimizer.set_vertex_movement_limit_factor(1.e-5);
	optimizer.set_iteration_limit(20);
	optimizer.enable_culling(false);
	optimizer.run_instructions(&mesq_mesh, err);
	if (err) {std::cerr << err << std::endl; }
	{
		std::ostringstream out_name;
		out_name << "smoothed_mesh." << mesh->mpi_size << "." << mesh->mpi_rank << ".vtk";
		mesq_mesh.write_vtk(out_name.str().c_str(), err);
	}

	for (int ino = 0; ino < vertices.size() ; ino++){
		Mesh::VertexHandle vertex = vertices[ino];
		MsqVertex aux;
		mesq_mesh.vertices_get_coordinates( &vertex, &aux, 1, err );
		coords[3*ino+0] = aux[0];
		coords[3*ino+1] = aux[1];
		coords[3*ino+2] = aux[2];
	}


	//TODO procurando os elementos invalidos
	/*
	//procurando os elementos invalidos...
	for (int iel = 0; iel < elem_array.size() ; iel++){
		mesq_mesh->
		inverted_element_tag_name
	}
	 */

	free(fixed_nodes);
	free(pid);
	free(gid);
	free(mid);
}
