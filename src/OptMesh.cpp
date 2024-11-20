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
#include "hexa.h"


void OptLine1(hexa_tree_t* mesh, std::vector<double>& coords, sc_hash_array_t* hash_FixedNodes)
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

void MeshOptimization1(hexa_tree_t* mesh, std::vector<double>& coords, std::vector<int> material_fixed_nodes){

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
	OptLine1(mesh, coords, hash_FixedNodes);

	printf("     Surface Optimization\n");
	//OptSurface(mesh, coords, hash_FixedNodes);

	if(mesh->mpi_size == 1 || mesh->mpi_size !=1)
	{
		//Sequential implementation
		//the exterior boundaries were fixed
		//the interior
		printf("     Volume Optimization\n");
		//OptVolume(mesh, coords, hash_FixedNodes);
	}
	if(mesh->mpi_size != 1)
	{
		//Parallel implementation
		printf("     Parallel Volume Optimization\n");
		//OptVolumeParallel(mesh, coords, hash_FixedNodes);
	}
}
