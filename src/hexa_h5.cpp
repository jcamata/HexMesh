#include <iostream>
using std::cout;
using std::endl;
#include <assert.h>

#include <string>
#include "H5Cpp.h"
#include "hdf5.h"
using namespace H5;

#include <vector>
#include <fstream>
#include <iomanip>
using namespace std;

#include "hexa.h"


void hexa_mesh_write_h5(hexa_tree_t mesh, const char* root_name, std::vector<double> coords){

	char filename[80];
	sprintf(filename, "%s_%04d_%04d.h5",root_name , mesh.mpi_size, mesh.mpi_rank);


	H5std_string FILE_NAME(filename);

	/////
	H5std_string DATASET_NAME1("Nodes");
	H5std_string DATASET_NAME2("Sem3D/Hexa8");
	H5std_string DATASET_NAME3("Sem3D/Mat");

	int	RANK = 2;
	hsize_t dims[RANK];               // dataset dimensions
	hsize_t dim[1];
	///////

	// put the data in vectors
	std::vector<int> connect;
	std::vector<int> mat;
	for(int i = 0; i<mesh.elements.elem_count;i++){
		octant_t* h = (octant_t*) sc_array_index(&mesh.elements, i);
		mat.push_back(h->n_mat);
		for(int j=0;j<8;j++){
			connect.push_back(h->nodes[j].id);
		}
	}
	// Create a new file using default property lists.
	H5File file(FILE_NAME, H5F_ACC_TRUNC);

	//write the Nodes:
	dims[0] = mesh.local_n_nodes;
	dims[1] = 3;
	DataSpace *dataspace1 = new DataSpace (RANK, dims);

	// Create the dataset.
	DataSet *dataset1 = new DataSet (file.createDataSet(DATASET_NAME1,
			PredType::IEEE_F64LE, *dataspace1));

	DataSpace mspace1( RANK, dims );

	// Write the data to the dataset
	dataset1->write(&coords[0], PredType::NATIVE_DOUBLE, mspace1,mspace1);

	// Close the current dataset and data space.
	delete dataset1;
	delete dataspace1;

	// Create a group named "/MygGroup" in the file
	Group group(file.createGroup("/Sem3D"));

	//write the Elements:
	//
	dims[0] = mesh.local_n_elements;
	dims[1] = 8;
	dataspace1 = new DataSpace (RANK, dims);

	// Create the dataset in group "SEM3D".
	dataset1 = new DataSet (file.createDataSet(DATASET_NAME2,
			PredType::STD_U64LE, *dataspace1));
	DataSpace mspace2( RANK, dims );

	dataset1->write(&connect[0], PredType::NATIVE_INT,mspace2,mspace2);

	// Close the current dataset and data space.
	delete dataset1;
	delete dataspace1;

	//write the Material:
	//
	dim[0] = mesh.local_n_elements;
	dataspace1 = new DataSpace (1, dim);

	// Create the dataset in group "SEM3D".
	dataset1 = new DataSet (file.createDataSet(DATASET_NAME3,
			PredType::STD_I64LE, *dataspace1));
	DataSpace mspace3( RANK, dim );

	dataset1->write(&mat[0], PredType::NATIVE_INT,mspace3,mspace3);

	// Close the current dataset and data space.
	delete dataset1;
	delete dataspace1;

	coords.clear();
	connect.clear();
	mat.clear();

	sprintf(filename, "%s_%04d_%04d.h5.xmf",root_name , mesh.mpi_size, mesh.mpi_rank);

	FILE *fid = fopen (filename,"w");
	sprintf(filename, "%s_%04d_%04d.h5",root_name , mesh.mpi_size, mesh.mpi_rank);


	fprintf(fid,"<?xml version=\"1.0\" ?>\n");
	fprintf(fid,"<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\">\n");
	fprintf(fid,"<Xdmf Version=\"2.0\" xmlns:xi=\"http://www.w3.org/2001/XInclude\">\n");
	fprintf(fid,"<Domain>\n");
	fprintf(fid,"<Grid GridType=\"Uniform\" Name=\"main\"><Geometry Type=\"XYZ\">\n");
	fprintf(fid,"<DataItem Dimensions=\"%d 3\" Format=\"HDF\" NumberType=\"Float\" Precision=\"8\">%s:/Nodes</DataItem>\n",mesh.local_n_nodes,filename);
	fprintf(fid,"</Geometry>\n");
	fprintf(fid,"<Topology NumberOfElements=\"%d\" Type=\"Hexahedron\">\n",mesh.local_n_elements);
	fprintf(fid,"<DataItem Dimensions=\"%d 8\" Format=\"HDF\" NumberType=\"UInt\" Precision=\"8\">%s:/Sem3D/Hexa8</DataItem>\n",mesh.local_n_elements,filename);
	fprintf(fid,"</Topology>\n");
	fprintf(fid,"<Attribute AttributeType=\"Scalar\" Center=\"Cell\" Dimensions=\"%d\" Name=\"Mat\">\n",mesh.local_n_elements);
	fprintf(fid,"<DataItem Dimensions=\"%d\" Format=\"HDF\" NumberType=\"Int\" Precision=\"8\">%s:/Sem3D/Mat</DataItem>\n",mesh.local_n_elements,filename);
	fprintf(fid,"</Attribute>\n");
	fprintf(fid,"</Grid></Domain></Xdmf>");


	fclose (fid);
}

