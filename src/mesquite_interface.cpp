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

void hexa_insert_shared_element(hexa_tree_t* mesh, sc_hash_array_t    *shared_element, octant_t* elem,std::vector<double>& coords, int processor)
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
			se->nodes[ino].id = mesh->global_id[elem->nodes[ino].id];
			se->nodes[ino].fixed = elem->nodes[ino].fixed;
			se->coord[ino][0] = coords[3*elem->nodes[ino].id + 0];
			se->coord[ino][1] = coords[3*elem->nodes[ino].id + 1];
			se->coord[ino][2] = coords[3*elem->nodes[ino].id + 2];
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

	bool deb = false;
	if(deb)
	{
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

	for(int iel = 0; iel < mesh->outsurf.elem_count; iel++)
	{
		octant_t* elem = (octant_t*) sc_array_index (&mesh->outsurf, iel);

		for(int iedge = 0; iedge < 12; iedge++)
		{
			if(elem->edge[iedge].ref == true)
			{
				if(iedge == 0)
				{
					aux0.push_back(elem->nodes[EdgeVerticesMap[iedge][0]].id);
					aux0.push_back(elem->nodes[EdgeVerticesMap[iedge][1]].id);
					if(deb) mesh->part_nodes[elem->nodes[EdgeVerticesMap[iedge][0]].id] = 0;
					if(deb) mesh->part_nodes[elem->nodes[EdgeVerticesMap[iedge][1]].id] = 0;
				}
				if(iedge == 1)
				{
					aux1.push_back(elem->nodes[EdgeVerticesMap[iedge][0]].id);
					aux1.push_back(elem->nodes[EdgeVerticesMap[iedge][1]].id);
					if(deb) mesh->part_nodes[elem->nodes[EdgeVerticesMap[iedge][0]].id] = 1;
					if(deb) mesh->part_nodes[elem->nodes[EdgeVerticesMap[iedge][1]].id] = 1;
				}
				if(iedge == 2)
				{
					aux2.push_back(elem->nodes[EdgeVerticesMap[iedge][0]].id);
					aux2.push_back(elem->nodes[EdgeVerticesMap[iedge][1]].id);
					if(deb) mesh->part_nodes[elem->nodes[EdgeVerticesMap[iedge][0]].id] = 2;
					if(deb) mesh->part_nodes[elem->nodes[EdgeVerticesMap[iedge][1]].id] = 2;
				}
				if(iedge == 3)
				{
					aux3.push_back(elem->nodes[EdgeVerticesMap[iedge][0]].id);
					aux3.push_back(elem->nodes[EdgeVerticesMap[iedge][1]].id);
					if(deb) mesh->part_nodes[elem->nodes[EdgeVerticesMap[iedge][0]].id] = 3;
					if(deb) mesh->part_nodes[elem->nodes[EdgeVerticesMap[iedge][1]].id] = 3;
				}

				if(iedge == 4)
				{
					aux4.push_back(elem->nodes[EdgeVerticesMap[iedge][0]].id);
					aux4.push_back(elem->nodes[EdgeVerticesMap[iedge][1]].id);
					if(deb) mesh->part_nodes[elem->nodes[EdgeVerticesMap[iedge][0]].id] = 4;
					if(deb) mesh->part_nodes[elem->nodes[EdgeVerticesMap[iedge][1]].id] = 4;
				}
				if(iedge == 5)
				{
					aux5.push_back(elem->nodes[EdgeVerticesMap[iedge][0]].id);
					aux5.push_back(elem->nodes[EdgeVerticesMap[iedge][1]].id);
					if(deb) mesh->part_nodes[elem->nodes[EdgeVerticesMap[iedge][0]].id] = 5;
					if(deb) mesh->part_nodes[elem->nodes[EdgeVerticesMap[iedge][1]].id] = 5;
				}
				if(iedge == 6)
				{
					aux6.push_back(elem->nodes[EdgeVerticesMap[iedge][0]].id);
					aux6.push_back(elem->nodes[EdgeVerticesMap[iedge][1]].id);
					if(deb) mesh->part_nodes[elem->nodes[EdgeVerticesMap[iedge][0]].id] = 6;
					if(deb) mesh->part_nodes[elem->nodes[EdgeVerticesMap[iedge][1]].id] = 6;
				}
				if(iedge == 7)
				{
					aux7.push_back(elem->nodes[EdgeVerticesMap[iedge][0]].id);
					aux7.push_back(elem->nodes[EdgeVerticesMap[iedge][1]].id);
					if(deb) mesh->part_nodes[elem->nodes[EdgeVerticesMap[iedge][0]].id] = 7;
					if(deb) mesh->part_nodes[elem->nodes[EdgeVerticesMap[iedge][1]].id] = 7;
				}

				if(iedge == 8)
				{
					aux8.push_back(elem->nodes[EdgeVerticesMap[iedge][0]].id);
					aux8.push_back(elem->nodes[EdgeVerticesMap[iedge][1]].id);
					if(deb) mesh->part_nodes[elem->nodes[EdgeVerticesMap[iedge][0]].id] = 8;
					if(deb) mesh->part_nodes[elem->nodes[EdgeVerticesMap[iedge][1]].id] = 8;
				}
				if(iedge == 9)
				{
					aux9.push_back(elem->nodes[EdgeVerticesMap[iedge][0]].id);
					aux9.push_back(elem->nodes[EdgeVerticesMap[iedge][1]].id);
					if(deb) mesh->part_nodes[elem->nodes[EdgeVerticesMap[iedge][0]].id] = 9;
					if(deb) mesh->part_nodes[elem->nodes[EdgeVerticesMap[iedge][1]].id] = 9;
				}
				if(iedge == 10)
				{
					aux10.push_back(elem->nodes[EdgeVerticesMap[iedge][0]].id);
					aux10.push_back(elem->nodes[EdgeVerticesMap[iedge][1]].id);
					if(deb) mesh->part_nodes[elem->nodes[EdgeVerticesMap[iedge][0]].id] = 10;
					if(deb) mesh->part_nodes[elem->nodes[EdgeVerticesMap[iedge][1]].id] = 10;
				}
				if(iedge == 11)
				{
					aux11.push_back(elem->nodes[EdgeVerticesMap[iedge][0]].id);
					aux11.push_back(elem->nodes[EdgeVerticesMap[iedge][1]].id);
					if(deb) mesh->part_nodes[elem->nodes[EdgeVerticesMap[iedge][0]].id] = 11;
					if(deb) mesh->part_nodes[elem->nodes[EdgeVerticesMap[iedge][1]].id] = 11;
				}
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

				node_t key;
				size_t position;
				int node = aux[ino];
				key.coord[0] = coords[3*node+0];
				key.coord[1] = coords[3*node+1];
				key.coord[2] = coords[3*node+2];
				key.node_id = node;

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
			octant_node_t * r = (octant_node_t*) sc_hash_array_insert_unique(hash_Fixed, &key, &position);
			if (r != NULL)
			{
				r->x = key.x;
				r->y = key.y;
				r->z = key.z;
				r->id = node->id;
			}
		}

		std::vector<int> elem_aux;
		//hash for the surface nodes
		sc_hash_array_t* hash_SurfaceNodes = sc_hash_array_new(sizeof(octant_node_t), node_hash_fn , node_equal_fn, &clamped);
		sc_hash_array_t* hash_NodesLocal = sc_hash_array_new(sizeof(octant_node_t), node_hash_fn , node_equal_fn, &clamped);
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
						r->id = elem->nodes[FaceNodesMap[isurf][ino]].id;
					}

					//refazendo a conectividade com nos de 0 a nvertices...
					//passando uma malha "local"
					octant_node_t* rlocal = (octant_node_t*) sc_hash_array_insert_unique(hash_NodesLocal, &key, &position);
					if (rlocal != NULL)
					{
						rlocal->x = key.x;
						rlocal->y = key.y;
						rlocal->z = key.z;
						rlocal->id = nlocal;
						conn_local.push_back(nlocal);
						nlocal++;
					}
					else
					{
						octant_node_t* rlocal = (octant_node_t*) sc_array_index(&hash_NodesLocal->a, position);
						conn_local.push_back(rlocal->id);
					}
				}

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
		printf("     Mesh normal to %d with: %d elements and %d vertex\n",isurf, nelem,nvertices);

		//nos fixos e coordenadas dos nos
		std::vector<double> Scoors(3*nvertices);
		bool *fixed_nodes = (bool*)malloc(nvertices*sizeof(bool));
		for(int ive = 0 ; ive < nvertices ;ive++ )
		{
			fixed_nodes[ive] = false;
			octant_node_t * node = (octant_node_t*) sc_array_index(&hash_SurfaceNodes->a, ive);
			Scoors[3*ive+0] = coords[3*node->id+0];
			Scoors[3*ive+1] = coords[3*node->id+1];
			Scoors[3*ive+2] = coords[3*node->id+2];

			octant_node_t key;
			size_t position;
			key.x = node->x;
			key.y = node->y;
			key.z = node->z;
			key.id = node->id;
			bool tre = sc_hash_array_lookup(hash_Fixed,&key,&position);
			if(tre)
			{
				fixed_nodes[ive] = true;
				mesh->part_nodes[node->id] = 1;
			}
		}

		if(nelem!=0)
		{
			//finalmente mandando a malha pro mesquite...
			Mesquite::MeshImpl mesq_mesh(nvertices,nelem,Mesquite::QUADRILATERAL, &fixed_nodes[0], &Scoors[0], &conn_local[0]);

			int n[3];
			if(isurf == 0)
			{
				n[0] = -1; n[1] = 0; n[2] = 0;
			}
			if(isurf == 1)
			{
				n[0] = 1; n[1] = 0; n[2] = 0;
			}
			if(isurf == 2)
			{
				n[1] = -1; n[0] = 0; n[2] = 0;
			}
			if(isurf == 3)
			{
				n[1] = 1; n[0] = 0; n[2] = 0;
			}
			if(isurf == 4)
			{
				n[2] = -1; n[1] = 0; n[0] = 0;
			}
			if(isurf == 5)
			{
				n[2] = 1; n[1] = 0; n[0] = 0;
			}
			Mesquite::PlanarDomain plane(Vector3D(n[0],n[1],n[2]), Vector3D(Scoors[0],Scoors[1],Scoors[2]));
			MeshDomainAssoc mesh_and_domain = MeshDomainAssoc(&mesq_mesh, &plane);
			assert(mesh_and_domain.are_compatible() == true);

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
			mesq_mesh.clear();

		}else{
			printf("Please check mesquite_interface! set of elements = 0.\n");
		}

		elem_aux.clear();
		conn_local.clear();
		Scoors.clear();
		free(fixed_nodes);
		sc_hash_array_destroy(hash_SurfaceNodes);
		sc_hash_array_destroy(hash_NodesLocal);
	}

}

void BuildMessage(hexa_tree_t* mesh, std::vector<double>& coords, sc_array_t ghostEl)
{
	//ids, nos se esta fixo ou não e coordenandas
	//além de um saber qual é o processso
	bool clamped = true;
	bool deb = false;
	sc_hash_array_t    *shared_element  = (sc_hash_array_t *)sc_hash_array_new(sizeof (shared_octant_t), sel_hash_id ,sel_equal_id , &clamped);

	for(int iel = 0; iel < mesh->elements.elem_count; iel++)
	{

		octant_t* elem = (octant_t*) sc_array_index(&mesh->elements, iel);
		if(elem->boundary)
		{
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
	}

	//if(deb) printf("Proc:%d tenho %d elementos para compartilhar\n",mesh->mpi_rank,shared_element->a.elem_count);

	//mapa de comunicacao
	// comm from proc [n+1] to proc [n]
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
				for(int j = 0; j < m->idxs.elem_count; ++j)
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

		printf("Sou o processador %d e devo ter %d ghost elements\n",mesh->mpi_rank,elid.size());
		for(int iel = 0; iel < elid.size(); iel++)
		{
			shared_octant_t* se = (shared_octant_t*) sc_array_push(&ghostEl);
			se->id = elid[iel];
			for(int ino = 0; ino < 8; ino++)
			{
				se->nodes[ino].id = nodeid[8*iel+ino];
				for(int j = 0; j < 3; j++) se->coord[ino][j] = coordcom[8*iel+3*ino+j];
			}
		}
	}

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
				for(int j = 0; j < m->idxs.elem_count; ++j)
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

		printf("Sou o processador %d e devo ter mais %d ghost elements\n",mesh->mpi_rank,elid.size());
		for(int iel = 0; iel < elid.size(); iel++)
		{
			shared_octant_t* se = (shared_octant_t*) sc_array_push(&ghostEl);
			se->id = elid[iel];
			for(int ino = 0; ino < 8; ino++)
			{
				se->nodes[ino].id = nodeid[8*iel+ino];
				for(int j = 0; j < 3; j++) se->coord[ino][j] = coordcom[8*iel+3*ino+j];
			}
		}
	}

	printf("Sou o processador %d e tenho %d ghost elements\n",mesh->mpi_rank,ghostEl.elem_count);


}

void OptVolume(hexa_tree_t* mesh, std::vector<double>& coords, sc_hash_array_t* hash_FixedNodes){
	Mesquite::MsqPrintError err(std::cout);

	//achando os ids globais...
	int local = mesh->elements.elem_count;
	int offset = 0;
	MPI_Scan(&local, &offset, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

	offset = offset-local;
	int count = offset;
	for(int iel = 0; iel < mesh->elements.elem_count; iel++)
	{
		octant_t * elem = (octant_t*) sc_array_index(&mesh->elements, iel);
		elem->id = count;
		count++;
	}

	sc_array_t ghostEl;
	sc_array_init(&ghostEl, sizeof(shared_octant_t));
	//montar o as mensagens para serem enviadas...
	BuildMessage(mesh,coords,ghostEl);

	//printf("Sou o processador %d e eu tenho %d ghost elements\n",mesh->mpi_rank,ghostEl.elem_count);



	//
	int nvertices = mesh->local_n_nodes;
	int nelem     = mesh->local_n_elements;
	std::vector<int> conn(8*nelem);
	bool   *fixed_nodes = (bool*)malloc(nvertices*sizeof(bool));
	size_t *gid         = (size_t*) malloc(nvertices*sizeof(size_t));
	int    *pid         = (int*) malloc (nvertices*sizeof(int));
	int    *mid         = (int*) malloc (nelem*sizeof(int));

	int c = 0;
	int tmp[8];
	for(int  iel = 0; iel < mesh->elements.elem_count; iel++)
	{
		octant_t * elem = (octant_t*) sc_array_index(&mesh->elements, iel);
		mid[iel] = elem->n_mat;
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
	}

	//int bnd_nodes = 0;
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
					fixed_nodes[elem->nodes[FaceNodesMap[isurf][ino]].id] = true;
					//mesh->part_nodes[elem->nodes[FaceNodesMap[isurf][ino]].id] = 1;
				}
			}
		}
	}

	for(int ino = 0; ino < hash_FixedNodes->a.elem_count; ino++)
	{
		node_t* node = (node_t*) sc_array_index (&hash_FixedNodes->a, ino);
		fixed_nodes[node->node_id] = true;
	}

	//Mesquite::ParallelMeshImpl mesq_mesh;

	/*
	Mesquite::MeshImpl mesq_mesh(nvertices,nelem,Mesquite::HEXAHEDRON, &fixed_nodes[0], &coords[0], &conn[0]);

//create parallel mesh instance, specifying tags containing parallel data
  Mesquite::ParallelMeshImpl parallel_mesh(&mesh, "GLOBAL_ID", "PROCESSOR_ID");
  Mesquite::ParallelHelperImpl helper;
  helper.set_communicator(MPI_COMM_WORLD);
  helper.set_parallel_mesh(&parallel_mesh);
  parallel_mesh.set_parallel_helper(&helper);

  //do Laplacian smooth
   //SmartLaplacianSmoother
  LaplaceWrapper optimizer;
  optimizer.set_vertex_movement_limit_factor(1.e-10);
  optimizer.set_iteration_limit(2000);
  optimizer.enable_culling(false);
  optimizer.run_instructions(&parallel_mesh, err);

	std::vector<MeshImpl::VertexHandle> vertices;
	mesq_mesh.get_all_vertices(vertices, err);

	for (int ino = 0; ino < vertices.size() ; ino++)
	{
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
	 */

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
	//printf("     Line Optimization\n");
	//OptLine(mesh, coords, hash_FixedNodes);

	//printf("     Surface Optimization\n");
	//OptSurface(mesh, coords, hash_FixedNodes);

	printf("     Volume Optimization\n");
	OptVolume(mesh, coords, hash_FixedNodes);
}

/*
	if(false){

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
		//stop_qa.disable_printing_results();

		//**************Set stopping criterion****************
		TerminationCriterion sc2;
		sc2.add_iteration_limit( 50 );
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
 *
 */
