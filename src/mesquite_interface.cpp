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

unsigned sel_hash_id(const void *v, const void *u)
{
	const shared_octant_t *q = (const shared_octant_t*) v;
	uint64_t a, b, c;

	a = (uint32_t) 1;
	b = (uint32_t) q->id;
	c = (uint32_t) 0;
	sc_hash_mix(a, b, c);
	sc_hash_final(a, b, c);
	return (unsigned) c;
}

int sel_equal_id(const void *v, const void *u, const void *w)
{
	const shared_octant_t *e1 = (const shared_octant_t*) v;
	const shared_octant_t *e2 = (const shared_octant_t*) u;

	return (unsigned) (e1->id == e2->id);
}

unsigned snode_hash_fn (const void *v, const void *u)
{
	const shared_node_t *q = (const shared_node_t*) v;
	uint32_t            a, b, c;

	a = (uint32_t) q->x;
	b = (uint32_t) q->y;
	c = (uint32_t) q->z;
	sc_hash_mix (a, b, c);
	sc_hash_final (a, b, c);
	return (unsigned) c;
}

int snode_equal_fn (const void *v1, const void *v2, const void *u)
{
	const shared_node_t *q1 = (const shared_node_t*) v1;
	const shared_node_t *q2 = (const shared_node_t*) v2;
	return (q1->x == q2->x && q1->y == q2->y && q1->z == q2->z);
}
void hexa_insert_shared_element(hexa_tree_t* mesh, sc_hash_array_t *shared_element, octant_t* elem, std::vector<double>& coords, int processor)
{
	size_t position;
	shared_octant_t *se;
	shared_octant_t key;
	int i;

	key.id = elem->id;

	if( processor < 0) return;

	se = (shared_octant_t*) sc_hash_array_insert_unique (shared_element, &key, &position);
	if(se != NULL)
	{
		se->id = elem->id;
		se->listSz = 1;
		se->rankList[0] = processor;
		for(int ino = 0; ino < 8; ino++)
		{
			se->coord[ino][0] = coords[3*elem->nodes[ino].id + 0];
			se->coord[ino][1] = coords[3*elem->nodes[ino].id + 1];
			se->coord[ino][2] = coords[3*elem->nodes[ino].id + 2];
			se->nodes[ino].id = mesh->global_id[elem->nodes[ino].id];
			se->nodes[ino].fixed = elem->nodes[ino].fixed;
			se->nodes[ino].x = elem->nodes[ino].x;
			se->nodes[ino].y = elem->nodes[ino].y;
			se->nodes[ino].z = elem->nodes[ino].z;
		}
	} else
	{
		se = (shared_octant_t*) sc_array_index(&shared_element->a, position);
		for(i=0; i < se->listSz; ++i)
			if(se->rankList[i] == processor) break;
		if(i == se->listSz){
			se->rankList[se->listSz] = processor;
			se->listSz++;
		}
	}
}

void OptLine(hexa_tree_t* mesh, std::vector<double>& coords, sc_hash_array_t* hash_FixedNodes)
{

	bool deb = true;
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

	for(int iel = 0; iel < mesh->elements.elem_count; iel++)
	{
		octant_t* elem = (octant_t*) sc_array_index (&mesh->elements, iel);

		for(int ino = 0; ino < 8; ino++)
		{
			//select only the local exterior edges
			//exterior local edges
			//color && fixed = -1;
			if(elem->nodes[ino].fixed == -1 && elem->nodes[ino].color < -10)
			{

				// z-y- = -0-11
				if(elem->nodes[ino].color == -11 || elem->nodes[ino].color == -30  || elem->nodes[ino].color == -31 ) aux0.push_back(elem->nodes[ino].id);
				// z-x+ = -1-11
				if(elem->nodes[ino].color == -12 || elem->nodes[ino].color == -31  || elem->nodes[ino].color == -32 ) aux1.push_back(elem->nodes[ino].id);
				// z-y+ = -2-11
				if(elem->nodes[ino].color == -13 || elem->nodes[ino].color == -32  || elem->nodes[ino].color == -33 ) aux2.push_back(elem->nodes[ino].id);
				// z-x- = -3-11
				if(elem->nodes[ino].color == -14 || elem->nodes[ino].color == -30  || elem->nodes[ino].color == -33 ) aux3.push_back(elem->nodes[ino].id);

				// x-y- = -4-11
				if(elem->nodes[ino].color == -15 || elem->nodes[ino].color == -30  || elem->nodes[ino].color == -34 ) aux4.push_back(elem->nodes[ino].id);
				// x+y- = -5-11
				if(elem->nodes[ino].color == -16 || elem->nodes[ino].color == -31  || elem->nodes[ino].color == -35 ) aux5.push_back(elem->nodes[ino].id);
				// x+y+ = -6-11
				if(elem->nodes[ino].color == -17 || elem->nodes[ino].color == -32  || elem->nodes[ino].color == -36 ) aux6.push_back(elem->nodes[ino].id);
				// x-y+ = -7-11
				if(elem->nodes[ino].color == -18 || elem->nodes[ino].color == -33  || elem->nodes[ino].color == -37 ) aux7.push_back(elem->nodes[ino].id);


				// z+y- = -8-11
				if(elem->nodes[ino].color == -19 || elem->nodes[ino].color == -34  || elem->nodes[ino].color == -35 ) aux8.push_back(elem->nodes[ino].id);
				// z+x+ = -9-11
				if(elem->nodes[ino].color == -20 || elem->nodes[ino].color == -35  || elem->nodes[ino].color == -36 ) aux9.push_back(elem->nodes[ino].id);
				// z+y+ = -10-11
				if(elem->nodes[ino].color == -21 || elem->nodes[ino].color == -36  || elem->nodes[ino].color == -37 ) aux10.push_back(elem->nodes[ino].id);
				// z+x- = -11-11
				if(elem->nodes[ino].color == -22 || elem->nodes[ino].color == -34  || elem->nodes[ino].color == -37 ) aux11.push_back(elem->nodes[ino].id);
			}
		}
	}

	for(int iedge = 0; iedge < 12 ; iedge++)
	{
		std::vector<int> aux;

		if(true)
		{
			if(iedge == 0)
			{
				std::vector<double> temp;
				//printf("Tenho %d elementos no vetor\n", aux0.size());
				for(int ino = 0; ino < aux0.size(); ino++)
				{
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
			if(iedge == 1)
			{
				std::vector<double> temp;
				for(int ino = 0; ino < aux1.size(); ino++) temp.push_back(coords[3*aux1[ino]+1]);

				std::vector<int> V(aux1.size());
				int x=0;
				std::iota(V.begin(),V.end(),x++); //Initializing
				sort(V.begin(),V.end(), [&](int i,int j){return temp[i]<temp[j];} );

				for(int ino = 0; ino < aux1.size(); ino++) aux.push_back(aux1[V[ino]]);
			}
			if(iedge == 2)
			{
				std::vector<double> temp;
				for(int ino = 0; ino < aux2.size(); ino++) temp.push_back(coords[3*aux2[ino]+0]);

				std::vector<int> V(aux2.size());
				int x=0;
				std::iota(V.begin(),V.end(),x++); //Initializing
				sort(V.begin(),V.end(), [&](int i,int j){return temp[i]<temp[j];} );

				for(int ino = 0; ino < aux2.size(); ino++) aux.push_back(aux2[V[ino]]);
			}
			if(iedge == 3)
			{
				std::vector<double> temp;
				for(int ino = 0; ino < aux3.size(); ino++) temp.push_back(coords[3*aux3[ino]+1]);

				std::vector<int> V(aux3.size());
				int x=0;
				std::iota(V.begin(),V.end(),x++); //Initializing
				sort(V.begin(),V.end(), [&](int i,int j){return temp[i]<temp[j];} );

				for(int ino = 0; ino < aux3.size(); ino++) aux.push_back(aux3[V[ino]]);
			}

			if(iedge == 4)
			{
				std::vector<double> temp;
				for(int ino = 0; ino < aux4.size(); ino++) temp.push_back(coords[3*aux4[ino]+2]);

				std::vector<int> V(aux4.size());
				int x=0;
				std::iota(V.begin(),V.end(),x++); //Initializing
				sort(V.begin(),V.end(), [&](int i,int j){return temp[i]<temp[j];} );

				for(int ino = 0; ino < aux4.size(); ino++) aux.push_back(aux4[V[ino]]);
			}
			if(iedge == 5)
			{
				std::vector<double> temp;
				for(int ino = 0; ino < aux5.size(); ino++) temp.push_back(coords[3*aux5[ino]+2]);

				std::vector<int> V(aux5.size());
				int x=0;
				std::iota(V.begin(),V.end(),x++); //Initializing
				sort(V.begin(),V.end(), [&](int i,int j){return temp[i]<temp[j];} );

				for(int ino = 0; ino < aux5.size(); ino++) aux.push_back(aux5[V[ino]]);
			}
			if(iedge == 6)
			{
				std::vector<double> temp;
				for(int ino = 0; ino < aux6.size(); ino++) temp.push_back(coords[3*aux6[ino]+2]);

				std::vector<int> V(aux6.size());
				int x=0;
				std::iota(V.begin(),V.end(),x++); //Initializing
				sort(V.begin(),V.end(), [&](int i,int j){return temp[i]<temp[j];} );

				for(int ino = 0; ino < aux6.size(); ino++) aux.push_back(aux6[V[ino]]);
			}
			if(iedge == 7)
			{
				std::vector<double> temp;
				for(int ino = 0; ino < aux7.size(); ino++) temp.push_back(coords[3*aux7[ino]+2]);

				std::vector<int> V(aux7.size());
				int x=0;
				std::iota(V.begin(),V.end(),x++); //Initializing
				sort(V.begin(),V.end(), [&](int i,int j){return temp[i]<temp[j];} );

				for(int ino = 0; ino < aux7.size(); ino++) aux.push_back(aux7[V[ino]]);
			}

			if(iedge == 8)
			{
				std::vector<double> temp;
				for(int ino = 0; ino < aux8.size(); ino++) temp.push_back(coords[3*aux8[ino]+0]);

				std::vector<int> V(aux8.size());
				int x=0;
				std::iota(V.begin(),V.end(),x++); //Initializing
				sort(V.begin(),V.end(), [&](int i,int j){return temp[i]<temp[j];} );

				for(int ino = 0; ino < aux8.size(); ino++) aux.push_back(aux8[V[ino]]);
			}
			if(iedge == 9)
			{
				std::vector<double> temp;
				for(int ino = 0; ino < aux9.size(); ino++) temp.push_back(coords[3*aux9[ino]+1]);

				std::vector<int> V(aux9.size());
				int x=0;
				std::iota(V.begin(),V.end(),x++); //Initializing
				sort(V.begin(),V.end(), [&](int i,int j){return temp[i]<temp[j];} );

				for(int ino = 0; ino < aux9.size(); ino++) aux.push_back(aux9[V[ino]]);
			}
			if(iedge == 10)
			{
				std::vector<double> temp;
				for(int ino = 0; ino < aux10.size(); ino++) temp.push_back(coords[3*aux10[ino]+0]);

				std::vector<int> V(aux10.size());
				int x=0;
				std::iota(V.begin(),V.end(),x++); //Initializing
				sort(V.begin(),V.end(), [&](int i,int j){return temp[i]<temp[j];} );

				for(int ino = 0; ino < aux10.size(); ino++) aux.push_back(aux10[V[ino]]);
			}
			if(iedge == 11)
			{
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
		double factor = 0.05;
		if(deb)printf("Proc:%d edge:%d aux size:%d\n",mesh->mpi_rank,iedge,aux.size());
		while(niter < 100 && tol >= 0.001 && aux.size() != 0)
		{
			for(int ino = 1; ino < (aux.size()-1); ino++)
			{

				octant_node_t key;
				size_t position;
				int node = aux[ino];
				octant_node_t* n = (octant_node_t*) sc_array_index(&mesh->nodes,aux[ino]);
				key.x = n->x;
				key.y = n->y;
				key.z = n->z;

				bool nodel = sc_hash_array_lookup(hash_FixedNodes, &key, &position);
				if(nodel)
				{

				}
				else
				{
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
				niter++;
			}
		}
	}
}

void OptSurface(hexa_tree_t* mesh, std::vector<double>& coords, sc_hash_array_t* hash_FixedNodes)
{

	Mesquite::MsqPrintError err(std::cout);

	for(int isurf = 0; isurf < 4 ; isurf ++)
	{
		//hash of the fixed nodes
		//I need this...
		//because I made the update for each surface
		bool clamped = true;
		sc_hash_array_t* hash_Fixed = sc_hash_array_new(sizeof(octant_node_t), node_hash_fn , node_equal_fn, &clamped);
		for(int ino = 0; ino < hash_FixedNodes->a.elem_count; ino++)
		{
			size_t position;
			octant_node_t key;
			octant_node_t* node = (octant_node_t*) sc_array_index (&hash_FixedNodes->a, ino);
			key.x = node->x;
			key.y = node->y;
			key.z = node->z;
			key.id = node->id;
			octant_node_t * r = (octant_node_t*) sc_hash_array_insert_unique(hash_Fixed, &key, &position);
			if (r != NULL)
			{
				r->x = key.x;
				r->y = key.y;
				r->z = key.z;
				r->id = key.id;
			}
		}

		std::vector<int> elem_aux;
		//hash for the surface nodes
		sc_hash_array_t* hash_SurfaceNodes = sc_hash_array_new(sizeof(octant_node_t), node_hash_fn , node_equal_fn, &clamped);
		sc_hash_array_t* hash_SurfaceLocal = sc_hash_array_new(sizeof(octant_node_t), node_hash_fn , node_equal_fn, &clamped);

		int nlocal = 0;
		std::vector<int> conn_local;
		for(int  iel = 0; iel < mesh->outsurf.elem_count; iel++)
		{
			octant_t * elem = (octant_t*) sc_array_index(&mesh->outsurf, iel);

			if(elem->surf[isurf].ext)
			{
				elem_aux.push_back(elem->id);
				for(int ino = 0; ino < 4; ino++)
				{
					size_t position;
					octant_node_t key;
					key.x = elem->nodes[FaceNodesMap[isurf][ino]].x;
					key.y = elem->nodes[FaceNodesMap[isurf][ino]].y;
					key.z = elem->nodes[FaceNodesMap[isurf][ino]].z;
					key.id = elem->nodes[FaceNodesMap[isurf][ino]].id;
					octant_node_t* r = (octant_node_t*) sc_hash_array_insert_unique(hash_SurfaceNodes, &key, &position);
					if (r != NULL)
					{
						r->x = key.x;
						r->y = key.y;
						r->z = key.z;
						r->id = key.id;

						r->fixed = 0;
						if(elem->nodes[FaceNodesMap[isurf][ino]].color < -10) r->fixed = -1;
					}

					size_t positionl;
					octant_node_t* rlocal = (octant_node_t*) sc_hash_array_insert_unique(hash_SurfaceLocal, &key, &positionl);
					if (rlocal != NULL)
					{
						rlocal->x = key.x;
						rlocal->y = key.y;
						rlocal->z = key.z;
						rlocal->id = nlocal;
						conn_local.push_back(nlocal);
						nlocal++;
					}else
					{
						octant_node_t* rl = (octant_node_t*) sc_array_index(&hash_SurfaceLocal->a, positionl);
						conn_local.push_back(rl->id);
					}
				}

				for(int ino = 0; ino < 8; ino++) mesh->part_nodes[elem->nodes[ino].id] = 0;

				for(int iedge = 0; iedge < 4; iedge++)
				{
					if(elem->edge[FaceEdgesMap[isurf][iedge]].ref)
					{
						for(int ino = 0; ino < 2; ino++)
						{
							int node = elem->nodes[EdgeVerticesMap[FaceEdgesMap[isurf][iedge]][ino]].id;
							size_t position;
							octant_node_t key;
							key.x = elem->nodes[EdgeVerticesMap[FaceEdgesMap[isurf][iedge]][ino]].x;
							key.y = elem->nodes[EdgeVerticesMap[FaceEdgesMap[isurf][iedge]][ino]].y;
							key.z = elem->nodes[EdgeVerticesMap[FaceEdgesMap[isurf][iedge]][ino]].z;
							octant_node_t* r = (octant_node_t*) sc_hash_array_insert_unique(hash_Fixed, &key, &position);
							if (r != NULL)
							{
								r->x = key.x;
								r->y = key.y;
								r->z = key.z;
								r->id = elem->nodes[EdgeVerticesMap[FaceEdgesMap[isurf][iedge]][ino]].id;;
							}
						}
					}
				}
			}
		}

		int nelem = elem_aux.size();
		int nvertices = hash_SurfaceNodes->a.elem_count;
		printf("     Mesh normal to %d in proc %d with: %d elements and %d vertex\n",isurf, mesh->mpi_rank,nelem,nvertices);

		//nos fixos e coordenadas dos nos
		std::vector<double> Scoors(3*nvertices);
		bool *fixed_nodes = (bool*)malloc(nvertices*sizeof(bool));
		for(int ive = 0 ; ive < nvertices ;ive++ )
		{
			octant_node_t * gnode = (octant_node_t*) sc_array_index(&hash_SurfaceNodes->a, ive);

			octant_node_t key;
			size_t position;
			key.x = gnode->x;
			key.y = gnode->y;
			key.z = gnode->z;

			bool lnode = sc_hash_array_lookup(hash_SurfaceLocal,&key,&position);
			octant_node_t * ln = (octant_node_t*) sc_array_index(&hash_SurfaceLocal->a, position);
			fixed_nodes[ln->id] = false;

			Scoors[3*ln->id+0] = coords[3*gnode->id+0];
			Scoors[3*ln->id+1] = coords[3*gnode->id+1];
			Scoors[3*ln->id+2] = coords[3*gnode->id+2];

			bool tre = sc_hash_array_lookup(hash_Fixed,&key,&position);
			if(tre || gnode->fixed == -1)
			{
				fixed_nodes[ln->id] = true;
				mesh->part_nodes[gnode->id] = 1;
			}
		}

		if(nelem!=0)
		{
			//finalmente mandando a malha pro mesquite...
			Mesquite::MeshImpl mesq_mesh(nvertices,nelem,Mesquite::QUADRILATERAL, &fixed_nodes[0], &Scoors[0], &conn_local[0]);

			int n[3];
			if(isurf == 0)
			{
				n[0] = 1; n[1] = 0; n[2] = 0;
			}
			if(isurf == 1)
			{
				n[0] = 1; n[1] = 0; n[2] = 0;
			}
			if(isurf == 2)
			{
				n[1] = 1; n[0] = 0; n[2] = 0;
			}
			if(isurf == 3)
			{
				n[1] = 1; n[0] = 0; n[2] = 0;
			}
			if(isurf == 4)
			{
				n[2] = 1; n[1] = 0; n[0] = 0;
			}
			if(isurf == 5)
			{
				n[2] = 1; n[1] = 0; n[0] = 0;
			}

			Mesquite::PlanarDomain plane(Vector3D(n[0],n[1],n[2]), Vector3D(Scoors[0],Scoors[1],Scoors[2]));
			MeshDomainAssoc mesh_and_domain = MeshDomainAssoc(&mesq_mesh, &plane);
			assert(mesh_and_domain.are_compatible() == true);
			{
				std::ostringstream out_name;
				out_name << "original_mesh." << isurf << "." << mesh->mpi_size << "."<< mesh->mpi_rank << ".vtk";
				mesq_mesh.write_vtk(out_name.str().c_str(), err);
			}
			//For mesh_and_domain.
			// creates an intruction queue
			InstructionQueue queue1;

			// creates a mean ratio quality metric ...
			// creates an edge length metric ...
			IdealWeightInverseMeanRatio shape_metric(err);
			//LInfTemplate lapl_met(&shape_metric);

			EdgeLengthQualityMetric lapl_met;
			lapl_met.set_averaging_method(QualityMetric::RMS);

			// creates the laplacian smoother  procedures
			//Here we use SmartLaplacianSmoother
			//it tries to avoid the inversion of the element...
			//try to keep this instead of laplacian
			SmartLaplacianSmoother lapl1;
			//LaplacianSmoother lapl1;
			QualityAssessor stop_qa=QualityAssessor(&shape_metric);
			stop_qa.add_quality_assessment(&lapl_met);
			stop_qa.disable_printing_results();

			//**************Set stopping criterion****************
			TerminationCriterion sc2;
			sc2.add_iteration_limit( 30 );
			sc2.add_cpu_time(120);
			sc2.add_relative_vertex_movement(1e-5);
			lapl1.set_outer_termination_criterion(&sc2);
			// adds 1 pass of pass1 to mesh_and_domain
			queue1.add_quality_assessor(&stop_qa,err);
			queue1.set_master_quality_improver(&lapl1, err);
			queue1.add_quality_assessor(&stop_qa,err);
			// launches optimization on mesh_and_domain
			queue1.run_instructions(&mesh_and_domain, err);

			int inverted_elements = 0;
			int inverted_samples = 0;
			bool test = stop_qa.get_inverted_element_count(inverted_elements, inverted_samples ,err);
			if(inverted_elements == 0 && test)
			{
				//get the vertices
				std::vector<MeshImpl::VertexHandle> vertices;
				mesh_and_domain.get_mesh()->get_all_vertices(vertices,err);

				//update the coords vector
				for (int ino = 0; ino < vertices.size() ; ino++)
				{
					Mesh::VertexHandle vertex = vertices[ino];
					MsqVertex aux_msq;
					mesh_and_domain.get_mesh()->vertices_get_coordinates( &vertex, &aux_msq, 1, err );
					octant_node_t * node = (octant_node_t*) sc_array_index(&hash_SurfaceNodes->a, ino);
					coords[3*node->id+0] = aux_msq[0];
					coords[3*node->id+1] = aux_msq[1];
					coords[3*node->id+2] = aux_msq[2];
				}
			}
			else
			{
				printf("     There are inverted elements in the opt mesh. "
						"No modification was made in the mesh\n");
			}
			mesq_mesh.clear();
		}

		elem_aux.clear();
		conn_local.clear();
		Scoors.clear();
		free(fixed_nodes);
		sc_hash_array_destroy(hash_SurfaceNodes);
	}

	//TODO opt surf for isurf four and five
}


void OptVolume(hexa_tree_t* mesh, std::vector<double>& coords, sc_hash_array_t* hash_FixedNodes){

	Mesquite::MsqPrintError err(std::cout);

	int nvertices = mesh->local_n_nodes;
	int nelem     = mesh->local_n_elements;

	std::vector<int> conn(8*nelem);
	bool   *fixed_nodes = (bool*)malloc(nvertices*sizeof(bool));

	int c = 0;
	int tmp[8];
	for(int  iel = 0; iel < mesh->elements.elem_count; iel++)
	{
		octant_t * elem = (octant_t*) sc_array_index(&mesh->elements, iel);
		for(int j = 0; j < 8; j++)
		{
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

		if(elem->tem == -1)
		{
			for(int j = 0; j < 8; j++)
			{
				fixed_nodes[elem->nodes[j].id] = true;
			}
		}

		int isurf;
		isurf = 0;
		elem->surf[isurf].ext = false;
		if(elem->nodes[FaceNodesMap[isurf][0]].x == mesh->x_start &&
				elem->nodes[FaceNodesMap[isurf][1]].x == mesh->x_start &&
				elem->nodes[FaceNodesMap[isurf][2]].x == mesh->x_start &&
				elem->nodes[FaceNodesMap[isurf][3]].x == mesh->x_start)
		{
			for(int ino = 0; ino < 4; ino++) fixed_nodes[elem->nodes[FaceNodesMap[isurf][0]].id] = true;
		}

		isurf = 1;
		elem->surf[isurf].ext = false;
		if(elem->nodes[FaceNodesMap[isurf][0]].x == mesh->x_end &&
				elem->nodes[FaceNodesMap[isurf][1]].x == mesh->x_end &&
				elem->nodes[FaceNodesMap[isurf][2]].x == mesh->x_end &&
				elem->nodes[FaceNodesMap[isurf][3]].x == mesh->x_end)
		{
			for(int ino = 0; ino < 4; ino++) fixed_nodes[elem->nodes[FaceNodesMap[isurf][0]].id] = true;
		}

		isurf = 2;
		elem->surf[isurf].ext = false;
		if(elem->nodes[FaceNodesMap[isurf][0]].y == mesh->y_start &&
				elem->nodes[FaceNodesMap[isurf][1]].y == mesh->y_start &&
				elem->nodes[FaceNodesMap[isurf][2]].y == mesh->y_start &&
				elem->nodes[FaceNodesMap[isurf][3]].y == mesh->y_start)
		{
			for(int ino = 0; ino < 4; ino++) fixed_nodes[elem->nodes[FaceNodesMap[isurf][0]].id] = true;
		}

		isurf = 3;
		elem->surf[isurf].ext = false;
		if(elem->nodes[FaceNodesMap[isurf][0]].y == mesh->y_end &&
				elem->nodes[FaceNodesMap[isurf][1]].y == mesh->y_end &&
				elem->nodes[FaceNodesMap[isurf][2]].y == mesh->y_end &&
				elem->nodes[FaceNodesMap[isurf][3]].y == mesh->y_end)
		{
			for(int ino = 0; ino < 4; ino++) fixed_nodes[elem->nodes[FaceNodesMap[isurf][0]].id] = true;
		}

		isurf = 4;
		elem->surf[isurf].ext = false;
		if(elem->nodes[FaceNodesMap[isurf][0]].z == 0 &&
				elem->nodes[FaceNodesMap[isurf][1]].z == 0 &&
				elem->nodes[FaceNodesMap[isurf][2]].z == 0 &&
				elem->nodes[FaceNodesMap[isurf][3]].z == 0)
		{
			for(int ino = 0; ino < 4; ino++) fixed_nodes[elem->nodes[FaceNodesMap[isurf][0]].id] = true;
		}

		isurf = 5;
		elem->surf[isurf].ext = false;
		if(elem->nodes[FaceNodesMap[isurf][0]].z == 3*mesh->max_z  &&
				elem->nodes[FaceNodesMap[isurf][1]].z == 3*mesh->max_z  &&
				elem->nodes[FaceNodesMap[isurf][2]].z == 3*mesh->max_z  &&
				elem->nodes[FaceNodesMap[isurf][3]].z == 3*mesh->max_z )
		{
			for(int ino = 0; ino < 4; ino++) fixed_nodes[elem->nodes[FaceNodesMap[isurf][0]].id] = true;
		}
	}

	for(int ino = 0; ino < hash_FixedNodes->a.elem_count; ino++)
	{
		node_t* node = (node_t*) sc_array_index (&hash_FixedNodes->a, ino);
		fixed_nodes[node->node_id] = true;
	}

	Mesquite::MeshImpl mesq_mesh(nvertices,nelem,Mesquite::HEXAHEDRON, &fixed_nodes[0], &coords[0], &conn[0]);

	//For mesh_and_domain.
	// creates an intruction queue
	InstructionQueue queue1;

	// creates a mean ratio quality metric ...
	ConditionNumberQualityMetric shape_metric;
	EdgeLengthQualityMetric lapl_met;
	lapl_met.set_averaging_method(QualityMetric::RMS);

	// creates the laplacian smoother  procedures
	//Here we use SmartLaplacianSmoother
	//it tries to avoid the inversion of the element...
	//try to keep this instead of laplacian
	SmartLaplacianSmoother lapl1;
	QualityAssessor stop_qa=QualityAssessor(&shape_metric);
	stop_qa.add_quality_assessment(&lapl_met);
	stop_qa.disable_printing_results();

	//**************Set stopping criterion****************
	TerminationCriterion sc2;
	sc2.add_iteration_limit( 10 );
	sc2.add_cpu_time(120);
	sc2.add_relative_vertex_movement(1e-5);
	lapl1.set_outer_termination_criterion(&sc2);
	// adds 1 pass of pass1 to mesh_and_domain
	queue1.add_quality_assessor(&stop_qa,err);
	queue1.set_master_quality_improver(&lapl1, err);
	queue1.add_quality_assessor(&stop_qa,err);
	// launches optimization on mesh_and_domain
	queue1.run_instructions(&mesq_mesh, err);

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

}

sc_array_t BuildMessage(hexa_tree_t* mesh, std::vector<double>& coords, sc_hash_array_t* shared_element)
{

	bool clamped = true;
	bool deb = false;
	sc_array_t ghostEl;
	sc_array_init(&ghostEl, sizeof(shared_octant_t));

	for(int iel = 0; iel < mesh->elements.elem_count; iel++)
	{
		octant_t* elem = (octant_t*) sc_array_index(&mesh->elements, iel);
		bool el_boundary = false;
		int isurf;
		isurf = 0;
		elem->surf[isurf].ext = false;
		if(elem->nodes[FaceNodesMap[isurf][0]].x == mesh->x_start &&
				elem->nodes[FaceNodesMap[isurf][1]].x == mesh->x_start &&
				elem->nodes[FaceNodesMap[isurf][2]].x == mesh->x_start &&
				elem->nodes[FaceNodesMap[isurf][3]].x == mesh->x_start)
		{
			elem->surf[isurf].ext = true;
			el_boundary = true;
			if(deb) for(int ino = 0; ino < 4; ino++) mesh->part_nodes[elem->nodes[FaceNodesMap[isurf][ino]].id] = isurf+20;
		}

		isurf = 1;
		elem->surf[isurf].ext = false;
		if(elem->nodes[FaceNodesMap[isurf][0]].x == mesh->x_end &&
				elem->nodes[FaceNodesMap[isurf][1]].x == mesh->x_end &&
				elem->nodes[FaceNodesMap[isurf][2]].x == mesh->x_end &&
				elem->nodes[FaceNodesMap[isurf][3]].x == mesh->x_end)
		{
			elem->surf[isurf].ext = true;
			el_boundary = true;
			if(deb) for(int ino = 0; ino < 4; ino++) mesh->part_nodes[elem->nodes[FaceNodesMap[isurf][ino]].id] = isurf+20;
		}

		isurf = 2;
		elem->surf[isurf].ext = false;
		if(elem->nodes[FaceNodesMap[isurf][0]].y == mesh->y_start &&
				elem->nodes[FaceNodesMap[isurf][1]].y == mesh->y_start &&
				elem->nodes[FaceNodesMap[isurf][2]].y == mesh->y_start &&
				elem->nodes[FaceNodesMap[isurf][3]].y == mesh->y_start)
		{
			elem->surf[isurf].ext = true;
			el_boundary = true;
			if(deb) for(int ino = 0; ino < 4; ino++) mesh->part_nodes[elem->nodes[FaceNodesMap[isurf][ino]].id] = isurf+20;
		}

		isurf = 3;
		elem->surf[isurf].ext = false;
		if(elem->nodes[FaceNodesMap[isurf][0]].y == mesh->y_end &&
				elem->nodes[FaceNodesMap[isurf][1]].y == mesh->y_end &&
				elem->nodes[FaceNodesMap[isurf][2]].y == mesh->y_end &&
				elem->nodes[FaceNodesMap[isurf][3]].y == mesh->y_end)
		{
			elem->surf[isurf].ext = true;
			el_boundary = true;
			if(deb) for(int ino = 0; ino < 4; ino++) mesh->part_nodes[elem->nodes[FaceNodesMap[isurf][ino]].id] = isurf+20;
		}

		if(elem->surf[2].ext)
		{
			if(elem->surf[0].ext)
			{
				hexa_insert_shared_element(mesh,shared_element,elem,coords,mesh->neighbors[0]);
			}
			else if(elem->surf[1].ext)
			{
				hexa_insert_shared_element(mesh,shared_element,elem,coords,mesh->neighbors[2]);
			}

			hexa_insert_shared_element(mesh,shared_element,elem,coords,mesh->neighbors[1]);

		}
		if(elem->surf[3].ext)
		{
			if(elem->surf[0].ext)
			{
				hexa_insert_shared_element(mesh,shared_element,elem,coords,mesh->neighbors[6]);
			}
			else if(elem->surf[1].ext)
			{
				hexa_insert_shared_element(mesh,shared_element,elem,coords,mesh->neighbors[8]);
			}

			hexa_insert_shared_element(mesh,shared_element,elem,coords,mesh->neighbors[7]);

		}
		if(elem->surf[0].ext)
		{
			hexa_insert_shared_element(mesh,shared_element,elem,coords,mesh->neighbors[3]);
		}
		if(elem->surf[1].ext)
		{
			hexa_insert_shared_element(mesh,shared_element,elem,coords,mesh->neighbors[5]);
		}
	}

	//if(deb) printf("Proc:%d tenho %d elementos para compartilhar\n",mesh->mpi_rank,shared_element->a.elem_count);

	//mapa de comunicacao
	// comm from proc [n] to proc [n+1]
	if(true)
	{
		sc_hash_array_t* SendTo   = (sc_hash_array_t *) sc_hash_array_new(sizeof(message_el_t), processors_hash_fn, processors_equal_fn, &clamped);
		sc_hash_array_t* RecvFrom = (sc_hash_array_t *) sc_hash_array_new(sizeof(message_el_t), processors_hash_fn, processors_equal_fn, &clamped);

		for(int ino = 0; ino < shared_element->a.elem_count; ++ino)
		{
			size_t position;
			shared_octant_t* se = (shared_octant_t*) sc_array_index(&shared_element->a,ino);
			for(int j = 0; j < se->listSz; j++)
			{
				if(se->rankList[j] > mesh->mpi_rank)
				{
					message_el_t* m = (message_el_t*)sc_hash_array_insert_unique(SendTo,&se->rankList[j],&position);
					if(m!=NULL)
					{
						m->rank  = se->rankList[j];
						sc_array_init(&m->idxs, sizeof(uint32_t));
						sc_array_init(&m->nodes, sizeof(uint32_t));
						sc_array_init(&m->coord, sizeof(double_t));
						uint32_t* p = (uint32_t*) sc_array_push(&m->idxs);
						*p = se->id;
						for(int ino = 0; ino < 8; ino++)
						{
							uint32_t* n = (uint32_t*) sc_array_push(&m->nodes);
							*n = se->nodes[ino].id;
							uint32_t* x = (uint32_t*) sc_array_push(&m->nodes);
							*x = se->nodes[ino].x;
							uint32_t* y = (uint32_t*) sc_array_push(&m->nodes);
							*y = se->nodes[ino].y;
							uint32_t* z = (uint32_t*) sc_array_push(&m->nodes);
							*z = se->nodes[ino].z;
							for(int icord = 0; icord < 3; icord++)
							{
								double_t* c = (double_t*) sc_array_push(&m->coord);
								*c = se->coord[ino][icord];
							}
						}
					}
					else
					{
						message_el_t* m = (message_el_t*)sc_array_index(&SendTo->a, position);
						uint32_t* p = (uint32_t*) sc_array_push(&m->idxs);
						*p = se->id;
						for(int ino = 0; ino < 8; ino++)
						{
							uint32_t* n = (uint32_t*) sc_array_push(&m->nodes);
							*n = se->nodes[ino].id;
							uint32_t* x = (uint32_t*) sc_array_push(&m->nodes);
							*x = se->nodes[ino].x;
							uint32_t* y = (uint32_t*) sc_array_push(&m->nodes);
							*y = se->nodes[ino].y;
							uint32_t* z = (uint32_t*) sc_array_push(&m->nodes);
							*z = se->nodes[ino].z;
							for(int icord = 0; icord < 3; icord++)
							{
								double_t* c = (double_t*) sc_array_push(&m->coord);
								*c = se->coord[ino][icord];
							}
						}
					}
				}
				else if (se->rankList[j] < mesh->mpi_rank)
				{
					message_el_t *m = (message_el_t*)sc_hash_array_insert_unique(RecvFrom,&se->rankList[j],&position);
					if(m!=NULL)
					{
						m->rank  = se->rankList[j];
						sc_array_init(&m->idxs, sizeof(uint32_t));
						sc_array_init(&m->nodes, sizeof(uint32_t));
						sc_array_init(&m->coord, sizeof(double_t));
						uint32_t* p = (uint32_t*) sc_array_push(&m->idxs);
						*p = se->id;
						for(int ino = 0; ino < 8; ino++)
						{
							uint32_t* n = (uint32_t*) sc_array_push(&m->nodes);
							*n = se->nodes[ino].id;
							uint32_t* x = (uint32_t*) sc_array_push(&m->nodes);
							*x = se->nodes[ino].x;
							uint32_t* y = (uint32_t*) sc_array_push(&m->nodes);
							*y = se->nodes[ino].y;
							uint32_t* z = (uint32_t*) sc_array_push(&m->nodes);
							*z = se->nodes[ino].z;
							for(int icord = 0; icord < 3; icord++)
							{
								double_t* c = (double_t*) sc_array_push(&m->coord);
								*c = se->coord[ino][icord];
							}
						}
					}
					else
					{
						message_el_t* m = (message_el_t*)sc_array_index(&RecvFrom->a, position);
						uint32_t* p = (uint32_t*) sc_array_push(&m->idxs);
						*p = se->id;
						for(int ino = 0; ino < 8; ino++)
						{
							uint32_t* n = (uint32_t*) sc_array_push(&m->nodes);
							*n = se->nodes[ino].id;
							uint32_t* x = (uint32_t*) sc_array_push(&m->nodes);
							*x = se->nodes[ino].x;
							uint32_t* y = (uint32_t*) sc_array_push(&m->nodes);
							*y = se->nodes[ino].y;
							uint32_t* z = (uint32_t*) sc_array_push(&m->nodes);
							*z = se->nodes[ino].z;
							for(int icord = 0; icord < 3; icord++)
							{
								double_t* c = (double_t*) sc_array_push(&m->coord);
								*c = se->coord[ino][icord];
							}
						}
					}
				}
			}
		}
		//send the el ids
		std::vector<int> elid;
		if(true)
		{
			//size for the ghost elements message
			int max_recvbuf_size = 0;
			for(int i = 0; i < RecvFrom->a.elem_count; i++)
			{
				message_el_t* m = (message_el_t*) sc_array_index(&RecvFrom->a, i);
				max_recvbuf_size +=  m->idxs.elem_count;
			}

			int max_sendbuf_size = 0;
			for(int i = 0; i < SendTo->a.elem_count; i++)
			{
				message_el_t* m = (message_el_t*) sc_array_index(&SendTo->a, i);
				max_sendbuf_size +=  m->idxs.elem_count;
			}

			int n_requests = RecvFrom->a.elem_count + SendTo->a.elem_count;

			//if(deb)printf("proc:%d send:%d recv:%d nrequest:%d\n",mesh->mpi_rank,max_recvbuf_size,max_sendbuf_size,n_requests);

			//comecando a comunicacao
			long long    * recvbuf    = (long long*)malloc(max_recvbuf_size*sizeof(long long));
			long long    * sendbuf    = (long long*)malloc(max_sendbuf_size*sizeof(long long));

			MPI_Request * requests = (MPI_Request*) malloc (n_requests*sizeof(MPI_Request));
			MPI_Status  * statuses = (MPI_Status*)  malloc (n_requests*sizeof(MPI_Status));
			int c = 0;

			//////////////// comm para o nome dos elementos
			int offset = 0;
			// post all non-blocking receives
			for(int i = 0; i < RecvFrom->a.elem_count; ++i) {
				message_el_t *m = (message_el_t*) sc_array_index(&RecvFrom->a, i);
				MPI_Irecv(&recvbuf[offset], m->idxs.elem_count, MPI_LONG_LONG, m->rank,0,MPI_COMM_WORLD, &requests[c]);
				offset += m->idxs.elem_count;
				c++;
			}

			assert(offset == max_recvbuf_size);

			offset = 0;
			for(int i = 0; i < SendTo->a.elem_count; ++i) {
				message_el_t *m = (message_el_t*) sc_array_index(&SendTo->a, i);
				for(int j = 0; j < m->idxs.elem_count; ++j)
				{
					int32_t *id = (int32_t*) sc_array_index(&m->idxs,j);
					sendbuf[offset+j] = (long long) *id;
				}
				MPI_Isend(&sendbuf[offset], m->idxs.elem_count, MPI_LONG_LONG, m->rank,0,MPI_COMM_WORLD, &requests[c]);
				offset += m->idxs.elem_count;
				c++;
			}
			assert(offset == max_sendbuf_size);

			assert(c == n_requests);

			MPI_Waitall(n_requests,requests,statuses);

			offset = 0;
			for(int i = 0; i < RecvFrom->a.elem_count; ++i) {
				message_el_t *m = (message_el_t*) sc_array_index(&RecvFrom->a, i);
				for(int j = 0; j < m->idxs.elem_count; ++j)
				{
					elid.push_back(recvbuf[offset+j]);
				}
				offset += m->idxs.elem_count;
			}
			free(recvbuf);
			free(sendbuf);
			free(requests);
			free(statuses);
		}

		if(deb)printf("Sou o proc %d e recebi %d elementos\n",mesh->mpi_rank,elid.size());

		//send the global nodes ids
		std::vector<int> nodeid;
		if(true)
		{
			//size for the ghost elements message
			int max_recvbuf_size = 0;
			for(int i = 0; i < RecvFrom->a.elem_count; i++)
			{
				message_el_t* m = (message_el_t*) sc_array_index(&RecvFrom->a, i);
				max_recvbuf_size +=  m->nodes.elem_count;
			}

			int max_sendbuf_size = 0;
			for(int i = 0; i < SendTo->a.elem_count; i++)
			{
				message_el_t* m = (message_el_t*) sc_array_index(&SendTo->a, i);
				max_sendbuf_size +=  m->nodes.elem_count;
			}

			int n_requests = RecvFrom->a.elem_count + SendTo->a.elem_count;

			if(deb)printf("proc:%d send:%d recv:%d nrequest:%d\n",mesh->mpi_rank,max_sendbuf_size,max_recvbuf_size,n_requests);

			//comecando a comunicacao
			long long    * recvbuf    = (long long*)malloc(max_recvbuf_size*sizeof(long long));
			long long    * sendbuf    = (long long*)malloc(max_sendbuf_size*sizeof(long long));

			MPI_Request * requests = (MPI_Request*) malloc (n_requests*sizeof(MPI_Request));
			MPI_Status  * statuses = (MPI_Status*)  malloc (n_requests*sizeof(MPI_Status));
			int c = 0;

			//////////////// comm para o nome dos elementos
			int offset = 0;
			// post all non-blocking receives
			for(int i = 0; i < RecvFrom->a.elem_count; ++i) {
				message_el_t *m = (message_el_t*) sc_array_index(&RecvFrom->a, i);
				MPI_Irecv(&recvbuf[offset], m->nodes.elem_count, MPI_LONG_LONG, m->rank,0,MPI_COMM_WORLD, &requests[c]);
				offset += m->nodes.elem_count;
				c++;
			}

			assert(offset == max_recvbuf_size);

			offset = 0;
			for(int i = 0; i < SendTo->a.elem_count; ++i) {
				message_el_t *m = (message_el_t*) sc_array_index(&SendTo->a, i);
				for(int j = 0; j < m->nodes.elem_count; ++j)
				{
					int32_t *id = (int32_t*) sc_array_index(&m->nodes,j);
					sendbuf[offset+j] = (long long) *id;
				}
				MPI_Isend(&sendbuf[offset], m->nodes.elem_count, MPI_LONG_LONG, m->rank,0,MPI_COMM_WORLD, &requests[c]);
				offset += m->nodes.elem_count;
				c++;
			}
			assert(offset == max_sendbuf_size);

			assert(c == n_requests);

			MPI_Waitall(n_requests,requests,statuses);

			offset = 0;
			for(int i = 0; i < RecvFrom->a.elem_count; ++i) {
				message_el_t *m = (message_el_t*) sc_array_index(&RecvFrom->a, i);
				for(int j = 0; j < m->nodes.elem_count; ++j)
				{
					nodeid.push_back(recvbuf[offset+j]);
				}
				offset += m->nodes.elem_count;
			}
			free(recvbuf);
			free(sendbuf);
			free(requests);
			free(statuses);
		}

		if(deb)printf("Sou o proc %d e recebi %d nos\n",mesh->mpi_rank,nodeid.size());

		//send the cooors
		std::vector<double> coordcom;
		if(true)
		{
			//size for the ghost elements message
			int max_recvbuf_size = 0;
			for(int i = 0; i < RecvFrom->a.elem_count; i++)
			{
				message_el_t* m = (message_el_t*) sc_array_index(&RecvFrom->a, i);
				max_recvbuf_size +=  m->coord.elem_count;
			}

			int max_sendbuf_size = 0;
			for(int i = 0; i < SendTo->a.elem_count; i++)
			{
				message_el_t* m = (message_el_t*) sc_array_index(&SendTo->a, i);
				max_sendbuf_size +=  m->coord.elem_count;
			}

			int n_requests = RecvFrom->a.elem_count + SendTo->a.elem_count;

			if(deb)printf("proc:%d send:%d recv:%d nrequest:%d\n",mesh->mpi_rank,max_sendbuf_size,max_recvbuf_size,n_requests);

			//comecando a comunicacao
			double    * recvbuf    = (double*)malloc(max_recvbuf_size*sizeof(double));
			double    * sendbuf    = (double*)malloc(max_sendbuf_size*sizeof(double));

			MPI_Request * requests = (MPI_Request*) malloc (n_requests*sizeof(MPI_Request));
			MPI_Status  * statuses = (MPI_Status*)  malloc (n_requests*sizeof(MPI_Status));
			int c = 0;

			//////////////// comm para o nome dos elementos
			int offset = 0;
			// post all non-blocking receives
			for(int i = 0; i < RecvFrom->a.elem_count; ++i) {
				message_el_t *m = (message_el_t*) sc_array_index(&RecvFrom->a, i);
				MPI_Irecv(&recvbuf[offset], m->coord.elem_count, MPI_DOUBLE, m->rank,0,MPI_COMM_WORLD, &requests[c]);
				offset += m->coord.elem_count;
				c++;
			}

			assert(offset == max_recvbuf_size);

			offset = 0;
			for(int i = 0; i < SendTo->a.elem_count; ++i) {
				message_el_t *m = (message_el_t*) sc_array_index(&SendTo->a, i);
				for(int j = 0; j < m->coord.elem_count; ++j)
				{
					double_t *id = (double_t*) sc_array_index(&m->coord,j);
					sendbuf[offset+j] = (double) *id;
				}
				MPI_Isend(&sendbuf[offset], m->coord.elem_count, MPI_DOUBLE, m->rank,0,MPI_COMM_WORLD, &requests[c]);
				offset += m->coord.elem_count;
				c++;
			}
			assert(offset == max_sendbuf_size);

			assert(c == n_requests);

			MPI_Waitall(n_requests,requests,statuses);

			offset = 0;
			for(int i = 0; i < RecvFrom->a.elem_count; ++i) {
				message_el_t *m = (message_el_t*) sc_array_index(&RecvFrom->a, i);
				for(int j = 0; j < m->coord.elem_count; ++j)
				{
					coordcom.push_back(recvbuf[offset+j]);
				}
				offset += m->coord.elem_count;
			}
			free(recvbuf);
			free(sendbuf);
			free(requests);
			free(statuses);
		}

		if(deb)printf("Sou o proc %d e recebi %d coords\n",mesh->mpi_rank,coordcom.size());

		if(deb) printf("Sou o processador %d e devo ter mais %d ghost elements\n",mesh->mpi_rank,elid.size());
		int cn = 0;
		int cc = 0;
		int cel = 0;
		for(int i = 0; i < RecvFrom->a.elem_count; i++)
		{
			message_el_t* m = (message_el_t*) sc_array_index(&RecvFrom->a, i);
			for(int iel = 0; iel < m->idxs.elem_count; iel++)
			{
				shared_octant_t* se = (shared_octant_t*) sc_array_push(&ghostEl);
				se->id = elid[cel]; cel++;
				se->rank = m->rank;
				for(int ino = 0; ino < 8; ino++)
				{
					se->nodes[ino].id = nodeid[cn];cn++;
					se->nodes[ino].x = nodeid[cn];cn++;
					se->nodes[ino].y = nodeid[cn];cn++;
					se->nodes[ino].z = nodeid[cn];cn++;
					for(int j = 0; j < 3; j++)
					{
						se->coord[ino][j] = coordcom[cc];
						cc++;
					}
				}
			}
		}
	}

	//mapa de comunicacao
	// comm from proc [n+1] to proc [n]
	if(true)
	{
		sc_hash_array_t* SendTo   = (sc_hash_array_t *) sc_hash_array_new(sizeof(message_el_t), processors_hash_fn, processors_equal_fn, &clamped);
		sc_hash_array_t* RecvFrom = (sc_hash_array_t *) sc_hash_array_new(sizeof(message_el_t), processors_hash_fn, processors_equal_fn, &clamped);

		for(int iel = 0; iel < shared_element->a.elem_count; ++iel)
		{
			size_t position;
			shared_octant_t* se = (shared_octant_t*) sc_array_index(&shared_element->a,iel);
			for(int j = 0; j < se->listSz; j++)
			{
				if(se->rankList[j] < mesh->mpi_rank)
				{
					message_el_t* m = (message_el_t*)sc_hash_array_insert_unique(SendTo,&se->rankList[j],&position);
					if(m!=NULL)
					{
						m->rank  = se->rankList[j];
						sc_array_init(&m->idxs, sizeof(uint32_t));
						sc_array_init(&m->nodes, sizeof(uint32_t));
						sc_array_init(&m->coord, sizeof(double_t));
						uint32_t* p = (uint32_t*) sc_array_push(&m->idxs);
						*p = se->id;
						for(int ino = 0; ino < 8; ino++)
						{
							uint32_t* n = (uint32_t*) sc_array_push(&m->nodes);
							*n = se->nodes[ino].id;
							uint32_t* x = (uint32_t*) sc_array_push(&m->nodes);
							*x = se->nodes[ino].x;
							uint32_t* y = (uint32_t*) sc_array_push(&m->nodes);
							*y = se->nodes[ino].y;
							uint32_t* z = (uint32_t*) sc_array_push(&m->nodes);
							*z = se->nodes[ino].z;
							for(int icord = 0; icord < 3; icord++)
							{
								double_t* c = (double_t*) sc_array_push(&m->coord);
								*c = se->coord[ino][icord];
							}
						}
					}
					else
					{
						message_el_t* m = (message_el_t*)sc_array_index(&SendTo->a, position);
						uint32_t* p = (uint32_t*) sc_array_push(&m->idxs);
						*p = se->id;
						for(int ino = 0; ino < 8; ino++)
						{
							uint32_t* n = (uint32_t*) sc_array_push(&m->nodes);
							*n = se->nodes[ino].id;
							uint32_t* x = (uint32_t*) sc_array_push(&m->nodes);
							*x = se->nodes[ino].x;
							uint32_t* y = (uint32_t*) sc_array_push(&m->nodes);
							*y = se->nodes[ino].y;
							uint32_t* z = (uint32_t*) sc_array_push(&m->nodes);
							*z = se->nodes[ino].z;
							for(int icord = 0; icord < 3; icord++)
							{
								double_t* c = (double_t*) sc_array_push(&m->coord);
								*c = se->coord[ino][icord];
							}
						}
					}
				}
				else if (se->rankList[j] > mesh->mpi_rank)
				{
					message_el_t *m = (message_el_t*)sc_hash_array_insert_unique(RecvFrom,&se->rankList[j],&position);
					if(m!=NULL)
					{
						m->rank  = se->rankList[j];
						sc_array_init(&m->idxs, sizeof(uint32_t));
						sc_array_init(&m->nodes, sizeof(uint32_t));
						sc_array_init(&m->coord, sizeof(double_t));
						uint32_t* p = (uint32_t*) sc_array_push(&m->idxs);
						*p = se->id;
						for(int ino = 0; ino < 8; ino++)
						{
							uint32_t* n = (uint32_t*) sc_array_push(&m->nodes);
							*n = se->nodes[ino].id;
							uint32_t* x = (uint32_t*) sc_array_push(&m->nodes);
							*x = se->nodes[ino].x;
							uint32_t* y = (uint32_t*) sc_array_push(&m->nodes);
							*y = se->nodes[ino].y;
							uint32_t* z = (uint32_t*) sc_array_push(&m->nodes);
							*z = se->nodes[ino].z;
							for(int icord = 0; icord < 3; icord++)
							{
								double_t* c = (double_t*) sc_array_push(&m->coord);
								*c = se->coord[ino][icord];
							}
						}
					}
					else
					{
						message_el_t* m = (message_el_t*)sc_array_index(&RecvFrom->a, position);
						uint32_t* p = (uint32_t*) sc_array_push(&m->idxs);
						*p = se->id;
						for(int ino = 0; ino < 8; ino++)
						{
							uint32_t* n = (uint32_t*) sc_array_push(&m->nodes);
							*n = se->nodes[ino].id;
							uint32_t* x = (uint32_t*) sc_array_push(&m->nodes);
							*x = se->nodes[ino].x;
							uint32_t* y = (uint32_t*) sc_array_push(&m->nodes);
							*y = se->nodes[ino].y;
							uint32_t* z = (uint32_t*) sc_array_push(&m->nodes);
							*z = se->nodes[ino].z;
							for(int icord = 0; icord < 3; icord++)
							{
								double_t* c = (double_t*) sc_array_push(&m->coord);
								*c = se->coord[ino][icord];
							}
						}
					}
				}
			}
		}

		//send the el ids
		std::vector<int> elid;
		if(true)
		{
			//size for the ghost elements message
			int max_recvbuf_size = 0;
			for(int i = 0; i < RecvFrom->a.elem_count; i++)
			{
				message_el_t* m = (message_el_t*) sc_array_index(&RecvFrom->a, i);
				max_recvbuf_size +=  m->idxs.elem_count;
			}

			int max_sendbuf_size = 0;
			for(int i = 0; i < SendTo->a.elem_count; i++)
			{
				message_el_t* m = (message_el_t*) sc_array_index(&SendTo->a, i);
				max_sendbuf_size +=  m->idxs.elem_count;
			}

			int n_requests = RecvFrom->a.elem_count + SendTo->a.elem_count;

			if(deb)printf("proc:%d send:%d recv:%d nrequest:%d\n",mesh->mpi_rank,max_recvbuf_size,max_sendbuf_size,n_requests);

			//comecando a comunicacao
			long long    * recvbuf    = (long long*)malloc(max_recvbuf_size*sizeof(long long));
			long long    * sendbuf    = (long long*)malloc(max_sendbuf_size*sizeof(long long));

			MPI_Request * requests = (MPI_Request*) malloc (n_requests*sizeof(MPI_Request));
			MPI_Status  * statuses = (MPI_Status*)  malloc (n_requests*sizeof(MPI_Status));
			int c = 0;

			//////////////// comm para o nome dos elementos
			int offset = 0;
			// post all non-blocking receives
			for(int i = 0; i < RecvFrom->a.elem_count; ++i) {
				message_el_t *m = (message_el_t*) sc_array_index(&RecvFrom->a, i);
				MPI_Irecv(&recvbuf[offset], m->idxs.elem_count, MPI_LONG_LONG, m->rank,0,MPI_COMM_WORLD, &requests[c]);
				offset += m->idxs.elem_count;
				c++;
			}

			assert(offset == max_recvbuf_size);

			offset = 0;
			for(int i = 0; i < SendTo->a.elem_count; ++i) {
				message_el_t *m = (message_el_t*) sc_array_index(&SendTo->a, i);
				for(int j = 0; j < m->idxs.elem_count; ++j)
				{
					int32_t *id = (int32_t*) sc_array_index(&m->idxs,j);
					sendbuf[offset+j] = (long long) *id;
				}
				MPI_Isend(&sendbuf[offset], m->idxs.elem_count, MPI_LONG_LONG, m->rank,0,MPI_COMM_WORLD, &requests[c]);
				offset += m->idxs.elem_count;
				c++;
			}
			assert(offset == max_sendbuf_size);

			assert(c == n_requests);

			MPI_Waitall(n_requests,requests,statuses);

			offset = 0;
			for(int i = 0; i < RecvFrom->a.elem_count; ++i) {
				message_el_t *m = (message_el_t*) sc_array_index(&RecvFrom->a, i);
				for(int j = 0; j < m->idxs.elem_count; ++j)
				{
					elid.push_back(recvbuf[offset+j]);
				}
				offset += m->idxs.elem_count;
			}
			free(recvbuf);
			free(sendbuf);
			free(requests);
			free(statuses);
		}

		if(deb)printf("Sou o proc %d e recebi %d elementos\n",mesh->mpi_rank,elid.size());

		//send the global nodes ids
		std::vector<int> nodeid;
		if(true)
		{
			//size for the ghost elements message
			int max_recvbuf_size = 0;
			for(int i = 0; i < RecvFrom->a.elem_count; i++)
			{
				message_el_t* m = (message_el_t*) sc_array_index(&RecvFrom->a, i);
				max_recvbuf_size +=  m->nodes.elem_count;
			}

			int max_sendbuf_size = 0;
			for(int i = 0; i < SendTo->a.elem_count; i++)
			{
				message_el_t* m = (message_el_t*) sc_array_index(&SendTo->a, i);
				max_sendbuf_size +=  m->nodes.elem_count;
			}

			int n_requests = RecvFrom->a.elem_count + SendTo->a.elem_count;

			if(deb)printf("proc:%d send:%d recv:%d nrequest:%d\n",mesh->mpi_rank,max_sendbuf_size,max_recvbuf_size,n_requests);

			//comecando a comunicacao
			long long    * recvbuf    = (long long*)malloc(max_recvbuf_size*sizeof(long long));
			long long    * sendbuf    = (long long*)malloc(max_sendbuf_size*sizeof(long long));

			MPI_Request * requests = (MPI_Request*) malloc (n_requests*sizeof(MPI_Request));
			MPI_Status  * statuses = (MPI_Status*)  malloc (n_requests*sizeof(MPI_Status));
			int c = 0;

			//////////////// comm para o nome dos elementos
			int offset = 0;
			// post all non-blocking receives
			for(int i = 0; i < RecvFrom->a.elem_count; ++i) {
				message_el_t *m = (message_el_t*) sc_array_index(&RecvFrom->a, i);
				MPI_Irecv(&recvbuf[offset], m->nodes.elem_count, MPI_LONG_LONG, m->rank,0,MPI_COMM_WORLD, &requests[c]);
				offset += m->nodes.elem_count;
				c++;
			}

			assert(offset == max_recvbuf_size);

			offset = 0;
			for(int i = 0; i < SendTo->a.elem_count; ++i) {
				message_el_t *m = (message_el_t*) sc_array_index(&SendTo->a, i);
				for(int j = 0; j < m->nodes.elem_count; ++j)
				{
					int32_t *id = (int32_t*) sc_array_index(&m->nodes,j);
					sendbuf[offset+j] = (long long) *id;
				}
				MPI_Isend(&sendbuf[offset], m->nodes.elem_count, MPI_LONG_LONG, m->rank,0,MPI_COMM_WORLD, &requests[c]);
				offset += m->nodes.elem_count;
				c++;
			}
			assert(offset == max_sendbuf_size);

			assert(c == n_requests);

			MPI_Waitall(n_requests,requests,statuses);

			offset = 0;
			for(int i = 0; i < RecvFrom->a.elem_count; ++i) {
				message_el_t *m = (message_el_t*) sc_array_index(&RecvFrom->a, i);
				for(int j = 0; j < m->nodes.elem_count; ++j)
				{
					nodeid.push_back(recvbuf[offset+j]);
				}
				offset += m->nodes.elem_count;
			}
			free(recvbuf);
			free(sendbuf);
			free(requests);
			free(statuses);
		}

		if(deb)printf("Sou o proc %d e recebi %d nos\n",mesh->mpi_rank,nodeid.size());

		//send the cooors
		std::vector<double> coordcom;
		if(true)
		{
			//size for the ghost elements message
			int max_recvbuf_size = 0;
			for(int i = 0; i < RecvFrom->a.elem_count; i++)
			{
				message_el_t* m = (message_el_t*) sc_array_index(&RecvFrom->a, i);
				max_recvbuf_size +=  m->coord.elem_count;
			}

			int max_sendbuf_size = 0;
			for(int i = 0; i < SendTo->a.elem_count; i++)
			{
				message_el_t* m = (message_el_t*) sc_array_index(&SendTo->a, i);
				max_sendbuf_size +=  m->coord.elem_count;
			}

			int n_requests = RecvFrom->a.elem_count + SendTo->a.elem_count;

			if(deb)printf("proc:%d send:%d recv:%d nrequest:%d\n",mesh->mpi_rank,max_sendbuf_size,max_recvbuf_size,n_requests);

			//comecando a comunicacao
			double    * recvbuf    = (double*)malloc(max_recvbuf_size*sizeof(double));
			double    * sendbuf    = (double*)malloc(max_sendbuf_size*sizeof(double));

			MPI_Request * requests = (MPI_Request*) malloc (n_requests*sizeof(MPI_Request));
			MPI_Status  * statuses = (MPI_Status*)  malloc (n_requests*sizeof(MPI_Status));
			int c = 0;

			//////////////// comm para o nome dos elementos
			int offset = 0;
			// post all non-blocking receives
			for(int i = 0; i < RecvFrom->a.elem_count; ++i) {
				message_el_t *m = (message_el_t*) sc_array_index(&RecvFrom->a, i);
				MPI_Irecv(&recvbuf[offset], m->coord.elem_count, MPI_DOUBLE, m->rank,0,MPI_COMM_WORLD, &requests[c]);
				offset += m->coord.elem_count;
				c++;
			}

			assert(offset == max_recvbuf_size);

			offset = 0;
			for(int i = 0; i < SendTo->a.elem_count; ++i) {
				message_el_t *m = (message_el_t*) sc_array_index(&SendTo->a, i);
				for(int j = 0; j < m->coord.elem_count; ++j)
				{
					double_t *id = (double_t*) sc_array_index(&m->coord,j);
					sendbuf[offset+j] = (double) *id;
				}
				MPI_Isend(&sendbuf[offset], m->coord.elem_count, MPI_DOUBLE, m->rank,0,MPI_COMM_WORLD, &requests[c]);
				offset += m->coord.elem_count;
				c++;
			}
			assert(offset == max_sendbuf_size);

			assert(c == n_requests);

			MPI_Waitall(n_requests,requests,statuses);

			offset = 0;
			for(int i = 0; i < RecvFrom->a.elem_count; ++i) {
				message_el_t *m = (message_el_t*) sc_array_index(&RecvFrom->a, i);
				for(int j = 0; j < m->coord.elem_count; ++j)
				{
					coordcom.push_back(recvbuf[offset+j]);
				}
				offset += m->coord.elem_count;
			}
			free(recvbuf);
			free(sendbuf);
			free(requests);
			free(statuses);
		}

		if(deb)printf("Sou o proc %d e recebi %d coords\n",mesh->mpi_rank,coordcom.size());

		if(deb) printf("Sou o processador %d e devo ter %d ghost elements\n",mesh->mpi_rank,elid.size());
		int cn = 0;
		int cc = 0;
		int cel = 0;
		for(int i = 0; i < RecvFrom->a.elem_count; i++)
		{
			message_el_t* m = (message_el_t*) sc_array_index(&RecvFrom->a, i);
			for(int iel = 0; iel < m->idxs.elem_count; iel++)
			{
				shared_octant_t* se = (shared_octant_t*) sc_array_push(&ghostEl);
				se->id = elid[cel]; cel++;
				se->rank = m->rank;
				for(int ino = 0; ino < 8; ino++)
				{
					se->nodes[ino].id = nodeid[cn];cn++;
					se->nodes[ino].x = nodeid[cn];cn++;
					se->nodes[ino].y = nodeid[cn];cn++;
					se->nodes[ino].z = nodeid[cn];cn++;
					for(int j = 0; j < 3; j++) {
						se->coord[ino][j] = coordcom[cc];
						cc++;
					}
				}
			}
		}
	}

	if(deb) printf("Sou o processador %d e tenho %d ghost elements\n",mesh->mpi_rank, ghostEl.elem_count);

	return ghostEl;
}

void Mesh2VTK(hexa_tree_t* mesh,std::vector<double>& lcoord,std::vector<int>& conn, bool *fixed_nodes,size_t *gid,int *pid)
{

	char fdname[80];
	sprintf(fdname,"MesqMesh_%04d_%04d.vtk", mesh->mpi_size, mesh->mpi_rank);
	FILE* outfile = fopen(fdname,"w");

	fprintf(outfile,"# vtk DataFile Version 3.0\n");
	fprintf(outfile,"Mesquite Mesh\n");
	fprintf(outfile,"ASCII\n");
	fprintf(outfile,"DATASET UNSTRUCTURED_GRID\n");
	fprintf(outfile,"POINTS %d double\n",lcoord.size()/3);
	for(int ino = 0; ino < lcoord.size()/3; ino++)
	{
		fprintf(outfile,"%f %f %f\n",lcoord[3*ino+0],lcoord[3*ino+1],lcoord[3*ino+2]);
	}

	fprintf(outfile,"CELLS %d %d\n",conn.size()/8, 9*conn.size()/8);
	for(int iel = 0; iel < conn.size()/8; iel++)
	{
		fprintf(outfile,"%d %d %d %d %d %d %d %d %d\n",8,
				conn[8*iel+0],conn[8*iel+1],conn[8*iel+2],conn[8*iel+3],
				conn[8*iel+4],conn[8*iel+5],conn[8*iel+6],conn[8*iel+7]);
	}

	fprintf(outfile,"CELL_TYPES %d\n",conn.size()/8);
	for(int iel = 0; iel < conn.size()/8; iel++)
	{
		fprintf(outfile,"12\n");
	}

	fprintf(outfile,"POINT_DATA %d\n",lcoord.size()/3);
	fprintf(outfile,"SCALARS fixed int\n");
	fprintf(outfile,"LOOKUP_TABLE default\n");
	for(int ino = 0; ino < lcoord.size()/3; ino++)
	{
		fprintf(outfile,"%d\n",fixed_nodes[ino]);
		//fprintf(outfile,"%d\n",0);
	}

	fprintf(outfile,"SCALARS GLOBAL_ID unsigned_long 1\n");
	fprintf(outfile,"LOOKUP_TABLE default\n");
	for(int ino = 0; ino < lcoord.size()/3; ino++)
	{
		fprintf(outfile,"%d\n",gid[ino]);
		//fprintf(outfile,"%d\n",0);
	}

	fprintf(outfile,"SCALARS PROCESSOR_ID unsigned_long 1\n");
	fprintf(outfile,"LOOKUP_TABLE default\n");
	for(int ino = 0; ino < lcoord.size()/3; ino++)
	{
		fprintf(outfile,"%d\n",pid[ino]);
		//fprintf(outfile,"%d\n",1);
	}

	fclose(outfile);
}

void OptVolumeParallel(hexa_tree_t* mesh, std::vector<double>& coords, sc_hash_array_t* hash_FixedNodes){
	Mesquite::MsqPrintError err(std::cout);

	bool deb = false;
	//achando os ids globais para os elementos...
	int local = mesh->elements.elem_count;
	int offset = 0;
	MPI_Scan(&local, &offset, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

	offset = offset-local;
	int count = offset;
	//atualizando os ids globais dos elementos
	for(int iel = 0; iel < mesh->elements.elem_count; iel++)
	{
		octant_t * elem = (octant_t*) sc_array_index(&mesh->elements, iel);
		elem->id = count;
		count++;
	}

	//enviando e recebendo os ghost elements...
	bool clamped = true;
	sc_hash_array_t *shared_element  = (sc_hash_array_t *)sc_hash_array_new(sizeof (shared_octant_t), sel_hash_id ,sel_equal_id , &clamped);
	sc_array_t ghostEl = BuildMessage(mesh, coords, shared_element);

	if(deb) printf("Sou o processador %d e eu tenho %d %d ghost elements\n",mesh->mpi_rank,
			ghostEl.elem_count,shared_element->a.elem_count);

	//pour envoyer la maillage au mesquite:
	//a processor ID of type int for every vertex that determines which processor owns a vertex and is in charge for smoothing this vertex;
	//a global ID of type size_t for every vertex that (at least in combination with the processor ID) is globally unique;
	//all necessary ghost elements and ghost nodes along the partition boundary must be provided.

	//create a hash of shared nodes that we obtain in the shared element
	sc_hash_array_t* nodeswg = (sc_hash_array_t *)sc_hash_array_new(sizeof (shared_node_t), snode_hash_fn, snode_equal_fn, &clamped);
	for(int ino = 0; ino < mesh->nodes.elem_count; ino++)
	{
		octant_node_t* node = (octant_node_t*) sc_array_index(&mesh->nodes,ino);
		size_t position;
		shared_node_t key;
		key.x = node->x;
		key.y = node->y;
		key.z = node->z;
		key.id = node->id;
		shared_node_t* n = (shared_node_t*) sc_hash_array_insert_unique (nodeswg, &key, &position);
		if(n!=NULL)
		{
			n->id = key.id;
			n->x = key.x;
			n->y = key.y;
			n->z = key.z;
			n->listSz = mesh->mpi_rank;
		}
	}

	for(int iel = 0; iel < ghostEl.elem_count; iel++)
	{
		shared_octant_t* elem = (shared_octant_t*) sc_array_index(&ghostEl,iel);

		for(int ino = 0; ino < 8; ino++)
		{
			size_t position;
			shared_node_t key;
			key.x = elem->nodes[ino].x;
			key.y = elem->nodes[ino].y;
			key.z = elem->nodes[ino].z;
			key.id = elem->nodes[ino].id;
			shared_node_t* n = (shared_node_t*) sc_hash_array_insert_unique(nodeswg, &key, &position);
			if(n!=NULL)
			{
				n->id = nodeswg->a.elem_count-1;
				n->x = key.x;
				n->y = key.y;
				n->z = key.z;
				n->listSz = elem->rank;
			}
			else
			{
				shared_node_t* n = (shared_node_t*) sc_array_index(&nodeswg->a,position);
				if(n->listSz < elem->rank) n->listSz = elem->rank;
			}
		}
	}

	//update the number of elements and nodes
	int nvertices = nodeswg->a.elem_count;
	int nelem = mesh->elements.elem_count + ghostEl.elem_count;

	printf("          The mesh that will be send to Mesquite is composed by %d "
			"elements and %d nodes in the proc %d\n", nelem,nvertices,mesh->mpi_rank);

	int c = 0;
	int tmp[8];
	std::vector<int> conn(8*nelem);
	std::vector<double> coords_local(3*nvertices);
	size_t *gid = (size_t*) malloc(nvertices*sizeof(size_t));
	int    *pid = (int*) malloc (nvertices*sizeof(int));
	bool   *fixed_nodes = (bool*)malloc(nvertices*sizeof(bool));

	//TODO pid ta com problema
	for(int ino = 0; ino < nvertices; ino++)
	{
		gid[ino] = -1;
		pid[ino] = -1;
	}

	//for the elements in the proc
	for(int iel = 0; iel < mesh->elements.elem_count; iel++)
	{
		octant_t * elem = (octant_t*) sc_array_index(&mesh->elements, iel);
		for(int ino = 0; ino < 8; ino++)
		{
			tmp[ino] = elem->nodes[ino].id;
			coords_local[3*elem->nodes[ino].id + 0] = coords[3*elem->nodes[ino].id + 0];
			coords_local[3*elem->nodes[ino].id + 1] = coords[3*elem->nodes[ino].id + 1];
			coords_local[3*elem->nodes[ino].id + 2] = coords[3*elem->nodes[ino].id + 2];
		}
		conn[c] = tmp[4]; c++;
		conn[c] = tmp[5]; c++;
		conn[c] = tmp[6]; c++;
		conn[c] = tmp[7]; c++;
		conn[c] = tmp[0]; c++;
		conn[c] = tmp[1]; c++;
		conn[c] = tmp[2]; c++;
		conn[c] = tmp[3]; c++;

		if(elem->tem == -1)
		{
			for(int ino = 0; ino < 8; ino++)
			{
				fixed_nodes[elem->nodes[ino].id] = 1;
			}
		}

		//adding nodes in ghost hash
		for(int ino = 0; ino < 8; ino++)
		{
			size_t position;
			shared_node_t key;
			key.x = elem->nodes[ino].x;
			key.y = elem->nodes[ino].y;
			key.z = elem->nodes[ino].z;
			key.id = elem->nodes[ino].id;
			shared_node_t* n = (shared_node_t*) sc_hash_array_insert_unique (nodeswg, &key, &position);
			if(n!=NULL)
			{
				n->id = key.id;
				n->x = key.x;
				n->y = key.y;
				n->z = key.z;

				//add global id
				gid[n->id] = mesh->global_id[key.id];
				//add proc id
				gid[n->id] = mesh->mpi_rank;
			}
			else
			{
				shared_node_t* n = (shared_node_t*) sc_array_index(&nodeswg->a,position);
				//add global id
				gid[n->id] = mesh->global_id[key.id];
				//add proc id
				pid[n->id] = mesh->mpi_rank;
			}
		}
	}

	//now we add the ghost elements and nodes in the conn,gid,pid
	for(int iel = 0; iel < ghostEl.elem_count; iel++)
	{
		shared_octant_t * elem = (shared_octant_t*) sc_array_index(&ghostEl, iel);
		for(int ino = 0; ino < 8; ino++)
		{
			size_t position;
			shared_node_t key;
			key.x = elem->nodes[ino].x;
			key.y = elem->nodes[ino].y;
			key.z = elem->nodes[ino].z;

			bool lnode = sc_hash_array_lookup(nodeswg,&key,&position);
			shared_node_t * sn = (shared_node_t*) sc_array_index(&nodeswg->a,position);

			//add global id
			//The id in the ghostEl is global
			//The id in nodeswg is local
			gid[sn->id] = elem->nodes[ino].id;
			//add proc id
			if(pid[sn->id] < sn->listSz) pid[sn->id] = sn->listSz;

			tmp[ino] = sn->id;

			coords_local[3*sn->id + 0] = elem->coord[ino][0];
			coords_local[3*sn->id + 1] = elem->coord[ino][1];
			coords_local[3*sn->id + 2] = elem->coord[ino][2];

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

	//TODO refaire: ajouter ghost
	//fixed nodes
	for(int iel = 0; iel < mesh->outsurf.elem_count; iel++)
	{
		octant_t * elem = (octant_t*) sc_array_index(&mesh->outsurf, iel);
		for(int isurf = 0 ; isurf < 6; isurf++)
		{
			if(elem->surf[isurf].ext)
			{
				for(int ino = 0; ino < 4; ino++)
				{
					fixed_nodes[elem->nodes[FaceNodesMap[isurf][ino]].id] = 1;
				}
			}
		}
	}

	for(int ino = 0; ino < hash_FixedNodes->a.elem_count; ino++)
	{
		octant_node_t* node = (octant_node_t*) sc_array_index (&hash_FixedNodes->a, ino);
		fixed_nodes[node->id] = 1;
	}

	//Mesh2VTK(mesh,coords_local,conn,fixed_nodes,gid,pid);

	//char fdname[80];
	//sprintf(fdname,"MesqMesh_%04d_%04d.vtk", mesh->mpi_size, mesh->mpi_rank);
	//Mesquite::MeshImpl parallel_mesh;
	//parallel_mesh.read_vtk(fdname,err);

	Mesquite::MeshImpl parallel_mesh(nvertices,nelem,Mesquite::HEXAHEDRON, &fixed_nodes[0], &coords_local[0], &conn[0]);

	std::vector<MeshImpl::VertexHandle> vertices;
	parallel_mesh.get_all_vertices(vertices, err);

	size_t default_gid = -1;
	int default_pid = 0;
	parallel_mesh.tag_create("GLOBAL_ID"   ,Mesquite::Mesh::HANDLE,1,&default_gid, err);
	parallel_mesh.tag_create("PROCESSOR_ID",Mesquite::Mesh::INT   ,1,&default_pid, err);

	TagHandle tag_processor_id = parallel_mesh.tag_get("PROCESSOR_ID", err);
	TagHandle tag_global_id    = parallel_mesh.tag_get("GLOBAL_ID", err);

	parallel_mesh.tag_set_vertex_data(tag_global_id   ,vertices.size(),&vertices[0], gid, err);
	parallel_mesh.tag_set_vertex_data(tag_processor_id,vertices.size(),&vertices[0], pid, err);

	{
		std::ostringstream out_name;
		out_name << "parallel_mesh." << mesh->mpi_size << "." << mesh->mpi_rank << ".vtk";
		parallel_mesh.write_vtk(out_name.str().c_str(), err);
	}

	//create parallel mesh instance, specifying tags containing parallel data
	Mesquite::ParallelMeshImpl mesq_mesh(&parallel_mesh, "GLOBAL_ID", "PROCESSOR_ID");
	Mesquite::ParallelHelperImpl helper;
	helper.set_communicator(MPI_COMM_WORLD);
	helper.set_parallel_mesh(&mesq_mesh);
	mesq_mesh.set_parallel_helper(&helper);

	//MeshDomain *domain=0;
	//MeshDomainAssoc mesh_and_domain = MeshDomainAssoc(&mesh, domain);

	if(false)
	{
		// creates an intruction queue
		InstructionQueue queue1;

		// creates a mean ratio quality metric ...
		UntangleBetaQualityMetric untangle_metric( 1.e-6 );

		QualityAssessor qa_untangle(&untangle_metric);
		queue1.add_quality_assessor(&qa_untangle, err); MSQ_ERRRTN(err);
		//q1.run_common( mesh_and_domain, pmesh, settings, err );

		LPtoPTemplate untangle_func( 2, &untangle_metric );
		ConjugateGradient untangle_solver( &untangle_func );
		//untangle_solver.set_debugging_level(3);

		//SteepestDescent untangle_solver( &untangle_func );
		TerminationCriterion untangle_inner("<type:untangle_inner>"), untangle_outer("<type:untangle_outer>");
		untangle_solver.use_global_patch();

		untangle_inner.add_absolute_quality_improvement( 0.0 );
		untangle_inner.add_iteration_limit( 20 );

		untangle_outer.add_absolute_quality_improvement( 0.0 );
		untangle_outer.add_iteration_limit( 10 );

		untangle_solver.set_inner_termination_criterion( &untangle_inner );
		untangle_solver.set_outer_termination_criterion( &untangle_outer );

		// adds 1 pass of pass1 to mesh_and_domain
		queue1.add_quality_assessor(&qa_untangle, err); MSQ_ERRRTN(err);
		queue1.set_master_quality_improver( &untangle_solver, err ); MSQ_ERRRTN(err);
		queue1.add_quality_assessor(&qa_untangle, err); MSQ_ERRRTN(err);
		// launches optimization on mesh_and_domain
		queue1.run_instructions(&mesq_mesh, err);
	}

	if(true)
	{
		// creates an intruction queue
		InstructionQueue queue1;

		// creates a mean ratio quality metric ...
		ConditionNumberQualityMetric shape_metric;
		EdgeLengthQualityMetric lapl_met;
		lapl_met.set_averaging_method(QualityMetric::RMS);

		// creates the laplacian smoother  procedures
		//Here we use SmartLaplacianSmoother
		//it tries to avoid the inversion of the element...
		//try to keep this instead of laplacian
		SmartLaplacianSmoother lapl1;
		QualityAssessor stop_qa=QualityAssessor(&shape_metric);
		stop_qa.add_quality_assessment(&lapl_met);
		//queue1.run_common(mesq_mesh,settings,err)
		stop_qa.disable_printing_results();

		//**************Set stopping criterion****************
		// For parallel runs, we generally need to have the inner and outer TerminationCriterion
		// have the same criteria else we can get an infinite loop (see VertexMover::loop_over_mesh)

		TerminationCriterion sc1("<type:SmartLaplacianSmoother_inner>");
		sc1.add_iteration_limit(10);
		sc1.add_cpu_time(10);
		sc1.add_relative_vertex_movement(1e-5);
		sc1.write_iterations("inner.gpt", err);


		TerminationCriterion sc2("<type:SmartLaplacianSmoother_outer>");;
		sc2.add_iteration_limit( 1 );
		sc2.add_cpu_time(120);
		sc2.add_relative_vertex_movement(1e-5);
		sc2.write_iterations("outer.gpt", err);

		lapl1.set_outer_termination_criterion(&sc2);
		lapl1.set_inner_termination_criterion(&sc1);

		// adds 1 pass of pass1 to mesh_and_domain
		queue1.add_quality_assessor(&stop_qa,err);
		queue1.set_master_quality_improver(&lapl1, err);
		queue1.add_quality_assessor(&stop_qa,err);
		// launches optimization on mesh_and_domain
		queue1.run_instructions(&mesq_mesh, err);
	}

	std::vector<MeshImpl::VertexHandle> ver;
	mesq_mesh.get_all_vertices(ver, err);

	/*
	for (int ino = 0; ino < vertices.size() ; ino++)
	{
		Mesh::VertexHandle vertex = vertices[ino];
		MsqVertex aux;
		mesq_mesh.vertices_get_coordinates( &vertex, &aux, 1, err );
		if(pid[ino] == mesh->mpi_rank)
		{
			coords[3*ino+0] = aux[0];
			coords[3*ino+1] = aux[1];
			coords[3*ino+2] = aux[2];
		}
	}
	 */

	{
		std::ostringstream out_name;
		out_name << "parallel_mesh_out." << mesh->mpi_size << "." << mesh->mpi_rank << ".vtk";
		//mesq_mesh.write_vtk(out_name.str().c_str(), err);
	}

	free(fixed_nodes);
	free(gid);
	free(pid);
	sc_hash_array_destroy(nodeswg);

}

void MeshOptimization(hexa_tree_t* mesh, std::vector<double>& coords, std::vector<int> material_fixed_nodes){

	//hash of the fixed nodes
	bool clamped = true;
	sc_hash_array_t* hash_FixedNodes = sc_hash_array_new(sizeof(octant_node_t), node_hash_fn , node_equal_fn, &clamped);

	for(int ino = 0; ino < material_fixed_nodes.size(); ino++){
		size_t position;
		octant_node_t *r;
		octant_node_t key;
		octant_node_t* node = (octant_node_t*) sc_array_index (&mesh->nodes, material_fixed_nodes[ino]);
		key.x = node->x;
		key.y = node->y;
		key.z = node->z;
		key.id = node->id;
		r = (octant_node_t*) sc_hash_array_insert_unique(hash_FixedNodes, &key, &position);
		if (r != NULL) {
			r->x = key.x;
			r->y = key.y;
			r->z = key.z;
			r->id = key.id;
		}
	}

	printf("     Line Optimization\n");
	OptLine(mesh, coords, hash_FixedNodes);

	printf("     Surface Optimization\n");
	OptSurface(mesh, coords, hash_FixedNodes);

	if(mesh->mpi_size == 1 || mesh->mpi_size !=1 && false)
	{
		//Sequential implementation
		//the exterior boundaries were fixed
		//the interior
		printf("     Volume Optimization\n");
		OptVolume(mesh, coords, hash_FixedNodes);
	}
	if(mesh->mpi_size != 1 && false)
	{
		//Parallel implementation
		printf("     Parallel Volume Optimization\n");
		OptVolumeParallel(mesh, coords, hash_FixedNodes);
	}
}
