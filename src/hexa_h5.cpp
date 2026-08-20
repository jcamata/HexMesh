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

void hexa_mesh_write_h5(hexa_tree_t *mesh, const char* root_name, std::vector<double> coords, const std::vector<double> *dtcrit)
{

	char filename[80];
	sprintf(filename, "%s_%04d_%04d.h5",root_name , mesh->mpi_size, mesh->mpi_rank);


	H5std_string FILE_NAME(filename);

	/////
	H5std_string DATASET_NAME1("Nodes");
	H5std_string DATASET_NAME2("Sem3D/Hexa8");
	H5std_string DATASET_NAME3("Sem3D/Mat");
	H5std_string DATASET_NAME4("Sem3D/Pad");
	H5std_string DATASET_NAME5("Sem3D/PillowType");

	int	RANK = 2;
	hsize_t dims[RANK];               // dataset dimensions
	hsize_t dim[1];
	///////
	unsigned int assign_elem_nodes[8];
	assign_elem_nodes[0] = 4;
	assign_elem_nodes[1] = 5;
	assign_elem_nodes[2] = 6;
	assign_elem_nodes[3] = 7;
	assign_elem_nodes[4] = 0;
	assign_elem_nodes[5] = 1;
	assign_elem_nodes[6] = 2;
	assign_elem_nodes[7] = 3;
	
	// put the data in vectors
	std::vector<int> connect;
	std::vector<int> mat;
	std::vector<int> pad;
	std::vector<int> pillow_type;
	for(int i = 0; i<mesh->elements.elem_count;i++){
		octant_t* h = (octant_t*) sc_array_index(&mesh->elements, i);
		mat.push_back(h->n_mat);
		pad.push_back(h->pad);
		for(int j=0;j<8;j++){
			connect.push_back(h->nodes[assign_elem_nodes[j]].id);
		}
		// PillowType:
		//   0 = mat-1, not at interface
		//   1 = mat-1, interface layer (has surface node fixed==1)
		//   2 = mat-0, original element (nodes remapped to pillow positions)
		//   3 = pillow element (new element created by Pillowing, level==-1)
		int pt;
		if (h->n_mat == 0 && h->level == -1) {
			pt = 3;
		} else if (h->n_mat == 0) {
			pt = 2;
		} else {
			bool on_interface = false;
			for (int j = 0; j < 8; j++)
				if (h->nodes[j].fixed == 1) { on_interface = true; break; }
			pt = on_interface ? 1 : 0;
		}
		pillow_type.push_back(pt);
	}

	// Create a new file using default property lists.
	H5File file(filename, H5F_ACC_TRUNC);

	//write the Nodes:
	dims[0] = mesh->local_n_nodes;
	dims[1] = 3;
	DataSpace *dataspace1 = new DataSpace (RANK, dims);

	// Create the dataset.
	DataSet *dataset1 = new DataSet (file.createDataSet("Nodes",
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
	dims[0] = mesh->local_n_elements;
	dims[1] = 8;
	dataspace1 = new DataSpace (RANK, dims);

	// Create the dataset in group "SEM3D".
	dataset1 = new DataSet (file.createDataSet("Sem3D/Hexa8",
			PredType::STD_U64LE, *dataspace1));
	DataSpace mspace2( RANK, dims );

	dataset1->write(&connect[0], PredType::NATIVE_INT,mspace2,mspace2);

	// Close the current dataset and data space.
	delete dataset1;
	delete dataspace1;

	//write the Material:
	//
	dim[0] = mesh->local_n_elements;
	dataspace1 = new DataSpace (1, dim);
	RANK = 1;
	// Create the dataset in group "SEM3D".
	dataset1 = new DataSet (file.createDataSet("Sem3D/Mat",
			PredType::STD_I64LE, *dataspace1));
	DataSpace mspace3( RANK, dim );

	dataset1->write(&mat[0], PredType::NATIVE_INT,mspace3,mspace3);

	// Close the current dataset and data space.
	delete dataset1;
	delete dataspace1;

	//write the Pad:
	//
	dim[0] = mesh->local_n_elements;
	dataspace1 = new DataSpace (1, dim);
	RANK = 1;
	// Create the dataset in group "SEM3D".
	dataset1 = new DataSet (file.createDataSet("Sem3D/Pad",
			PredType::STD_I64LE, *dataspace1));
	DataSpace mspace4( RANK, dim );

	dataset1->write(&pad[0], PredType::NATIVE_INT,mspace4,mspace4);

	// Close the current dataset and data space.
	delete dataset1;
	delete dataspace1;

	//write the PillowType:
	//
	dim[0] = mesh->local_n_elements;
	dataspace1 = new DataSpace (1, dim);
	RANK = 1;
	dataset1 = new DataSet (file.createDataSet("Sem3D/PillowType",
			PredType::STD_I64LE, *dataspace1));
	DataSpace mspace5( RANK, dim );
	dataset1->write(&pillow_type[0], PredType::NATIVE_INT, mspace5, mspace5);
	delete dataset1;
	delete dataspace1;

	// write DtCrit if provided
	if (dtcrit != NULL && (int32_t)dtcrit->size() >= mesh->local_n_elements) {
		dim[0] = mesh->local_n_elements;
		dataspace1 = new DataSpace(1, dim);
		RANK = 1;
		dataset1 = new DataSet(file.createDataSet("Sem3D/DtCrit", PredType::IEEE_F64LE, *dataspace1));
		DataSpace mspace_dt(RANK, dim);
		dataset1->write(&(*dtcrit)[0], PredType::NATIVE_DOUBLE, mspace_dt, mspace_dt);
		delete dataset1;
		delete dataspace1;
	}

	// write the Centroids (one point per element) — diagnostic point cloud so
	// degenerate / zero-volume elements (which do not render as cells in
	// ParaView) are still visible as points: confirms the element EXISTS rather
	// than being truly missing. Colour the point cloud by Mat / PillowType.
	std::vector<double> centroids(3 * (size_t)mesh->local_n_elements);
	for (int i = 0; i < mesh->local_n_elements; i++) {
		double cx = 0, cy = 0, cz = 0;
		for (int j = 0; j < 8; j++) {
			int nid = connect[8*i + j];
			cx += coords[3*nid+0]; cy += coords[3*nid+1]; cz += coords[3*nid+2];
		}
		centroids[3*i+0] = cx/8.0; centroids[3*i+1] = cy/8.0; centroids[3*i+2] = cz/8.0;
	}
	RANK = 2;
	dims[0] = mesh->local_n_elements; dims[1] = 3;
	dataspace1 = new DataSpace(RANK, dims);
	dataset1 = new DataSet(file.createDataSet("Sem3D/Centroids", PredType::IEEE_F64LE, *dataspace1));
	{
		DataSpace mspace6(RANK, dims);
		dataset1->write(&centroids[0], PredType::NATIVE_DOUBLE, mspace6, mspace6);
	}
	delete dataset1;
	delete dataspace1;
	centroids.clear();

	//coords.clear();
	connect.clear();
	mat.clear();
	pad.clear();
	pillow_type.clear();

	sprintf(filename, "%s_%04d_%04d.h5.xmf",root_name , mesh->mpi_size, mesh->mpi_rank);

	FILE *fid = fopen (filename,"w");
	sprintf(filename, "%s_%04d_%04d.h5",root_name , mesh->mpi_size, mesh->mpi_rank);


	fprintf(fid,"<?xml version=\"1.0\" ?>\n");
	fprintf(fid,"<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\">\n");
	fprintf(fid,"<Xdmf Version=\"2.0\" xmlns:xi=\"http://www.w3.org/2001/XInclude\">\n");
	fprintf(fid,"<Domain>\n");
	fprintf(fid,"<Grid GridType=\"Uniform\" Name=\"main\"><Geometry Type=\"XYZ\">\n");
	fprintf(fid,"<DataItem Dimensions=\"%d 3\" Format=\"HDF\" NumberType=\"Float\" Precision=\"8\">%s:/Nodes</DataItem>\n",mesh->local_n_nodes,filename);
	fprintf(fid,"</Geometry>\n");
	fprintf(fid,"<Topology NumberOfElements=\"%d\" Type=\"Hexahedron\">\n",mesh->local_n_elements);
	fprintf(fid,"<DataItem Dimensions=\"%d 8\" Format=\"HDF\" NumberType=\"UInt\" Precision=\"8\">%s:/Sem3D/Hexa8</DataItem>\n",mesh->local_n_elements,filename);
	fprintf(fid,"</Topology>\n");

	fprintf(fid,"<Attribute AttributeType=\"Scalar\" Center=\"Cell\" Dimensions=\"%d\" Name=\"Mat\">\n",mesh->local_n_elements);
	fprintf(fid,"<DataItem Dimensions=\"%d\" Format=\"HDF\" NumberType=\"Int\" Precision=\"8\">%s:/Sem3D/Mat</DataItem>\n",mesh->local_n_elements,filename);
	fprintf(fid,"</Attribute>\n");

	fprintf(fid,"<Attribute AttributeType=\"Scalar\" Center=\"Cell\" Dimensions=\"%d\" Name=\"Pad\">\n",mesh->local_n_elements);
	fprintf(fid,"<DataItem Dimensions=\"%d\" Format=\"HDF\" NumberType=\"Int\" Precision=\"8\">%s:/Sem3D/Pad</DataItem>\n",mesh->local_n_elements,filename);
	fprintf(fid,"</Attribute>\n");

	fprintf(fid,"<Attribute AttributeType=\"Scalar\" Center=\"Cell\" Dimensions=\"%d\" Name=\"PillowType\">\n",mesh->local_n_elements);
	fprintf(fid,"<DataItem Dimensions=\"%d\" Format=\"HDF\" NumberType=\"Int\" Precision=\"8\">%s:/Sem3D/PillowType</DataItem>\n",mesh->local_n_elements,filename);
	fprintf(fid,"</Attribute>\n");

	if (dtcrit != NULL && (int32_t)dtcrit->size() >= mesh->local_n_elements) {
		fprintf(fid,"<Attribute AttributeType=\"Scalar\" Center=\"Cell\" Dimensions=\"%d\" Name=\"DtCrit\">\n",mesh->local_n_elements);
		fprintf(fid,"<DataItem Dimensions=\"%d\" Format=\"HDF\" NumberType=\"Float\" Precision=\"8\">%s:/Sem3D/DtCrit</DataItem>\n",mesh->local_n_elements,filename);
		fprintf(fid,"</Attribute>\n");
	}

	fprintf(fid,"</Grid>\n");

	// second grid: one point per element at its centroid (diagnostic point
	// cloud). A point with no rendered cell around it => the element exists but
	// is degenerate (zero-volume/inverted); a gap with no point => truly missing.
	fprintf(fid,"<Grid GridType=\"Uniform\" Name=\"centroids\">\n");
	fprintf(fid,"<Topology TopologyType=\"Polyvertex\" NumberOfElements=\"%d\" NodesPerElement=\"1\"/>\n",mesh->local_n_elements);
	fprintf(fid,"<Geometry Type=\"XYZ\">\n");
	fprintf(fid,"<DataItem Dimensions=\"%d 3\" Format=\"HDF\" NumberType=\"Float\" Precision=\"8\">%s:/Sem3D/Centroids</DataItem>\n",mesh->local_n_elements,filename);
	fprintf(fid,"</Geometry>\n");
	fprintf(fid,"<Attribute AttributeType=\"Scalar\" Center=\"Node\" Name=\"Mat\">\n");
	fprintf(fid,"<DataItem Dimensions=\"%d\" Format=\"HDF\" NumberType=\"Int\" Precision=\"8\">%s:/Sem3D/Mat</DataItem>\n",mesh->local_n_elements,filename);
	fprintf(fid,"</Attribute>\n");
	fprintf(fid,"<Attribute AttributeType=\"Scalar\" Center=\"Node\" Name=\"PillowType\">\n");
	fprintf(fid,"<DataItem Dimensions=\"%d\" Format=\"HDF\" NumberType=\"Int\" Precision=\"8\">%s:/Sem3D/PillowType</DataItem>\n",mesh->local_n_elements,filename);
	fprintf(fid,"</Attribute>\n");
	fprintf(fid,"</Grid>\n");

	fprintf(fid,"</Domain></Xdmf>");


	fclose (fid);
}

void hexa_mesh_write_quality_h5(hexa_tree_t *mesh, const char *root_name,
                                const std::vector<double> &coords,
                                const std::vector<hex_quality_t> &qualities,
                                const std::vector<double> *dtcrit) {
	if (!mesh || mesh->elements.elem_count == 0) return;

	int n_elem = mesh->elements.elem_count;
	int n_nodes = coords.size() / 3;

	char filename[128];
	snprintf(filename, sizeof(filename), "%s_%04d_%04d.h5", root_name, mesh->mpi_size, mesh->mpi_rank);

	H5File file(filename, H5F_ACC_TRUNC);

	// Remap to VTK/h5 corner order
	static const int assign_elem_nodes[8] = {4, 5, 6, 7, 0, 1, 2, 3};

	std::vector<int> connect;
	connect.reserve(n_elem * 8);
	std::vector<int> mat, pad, pillow_type;
	mat.reserve(n_elem);
	pad.reserve(n_elem);
	pillow_type.reserve(n_elem);

	std::vector<double> v_scaledJac(n_elem), v_volume(n_elem), v_condNum(n_elem);
	std::vector<double> v_edgeRatio(n_elem), v_skew(n_elem), v_shape(n_elem);
	std::vector<double> v_oddy(n_elem), v_diagRatio(n_elem), v_taper(n_elem);
	std::vector<double> v_stretch(n_elem), v_minAngle(n_elem);

	for (int i = 0; i < n_elem; i++) {
		octant_t *h = (octant_t *) sc_array_index(&mesh->elements, i);
		mat.push_back(h->n_mat);
		pad.push_back(h->pad);
		for (int j = 0; j < 8; j++)
			connect.push_back(h->nodes[assign_elem_nodes[j]].id);

		int pt;
		if (h->n_mat == 0 && h->level == -1) pt = 3;
		else if (h->n_mat == 0) pt = 2;
		else {
			bool on_interface = false;
			for (int j = 0; j < 8; j++)
				if (h->nodes[j].fixed == 1) { on_interface = true; break; }
			pt = on_interface ? 1 : 0;
		}
		pillow_type.push_back(pt);

		if (i < (int)qualities.size()) {
			const auto &q = qualities[i];
			v_scaledJac[i]  = q.scaledJacobian;
			v_volume[i]     = q.volume;
			v_condNum[i]    = q.conditionNumber;
			v_edgeRatio[i]  = q.edgeRatio;
			v_skew[i]       = q.skew;
			v_shape[i]      = q.shape;
			v_oddy[i]       = q.oddy;
			v_diagRatio[i]  = q.diagonalRatio;
			v_taper[i]      = q.taper;
			v_stretch[i]    = q.stretch;
			v_minAngle[i]   = q.minFaceAngle;
		}
	}

	// Nodes Dataset
	hsize_t dims2D[2] = {(hsize_t)n_nodes, 3};
	DataSpace dspace_nodes(2, dims2D);
	DataSet dataset_nodes(file.createDataSet("Nodes", PredType::IEEE_F64LE, dspace_nodes));
	dataset_nodes.write(coords.data(), PredType::NATIVE_DOUBLE, dspace_nodes, dspace_nodes);

	// Sem3D Group
	Group group(file.createGroup("/Sem3D"));

	// Connectivity
	dims2D[0] = n_elem; dims2D[1] = 8;
	DataSpace dspace_conn(2, dims2D);
	DataSet dataset_conn(file.createDataSet("Sem3D/Hexa8", PredType::STD_U64LE, dspace_conn));
	dataset_conn.write(connect.data(), PredType::NATIVE_INT, dspace_conn, dspace_conn);

	// Mat, Pad, PillowType
	hsize_t dim1D[1] = {(hsize_t)n_elem};
	DataSpace dspace_1D(1, dim1D);

	DataSet dataset_mat(file.createDataSet("Sem3D/Mat", PredType::STD_I64LE, dspace_1D));
	dataset_mat.write(mat.data(), PredType::NATIVE_INT, dspace_1D, dspace_1D);

	DataSet dataset_pad(file.createDataSet("Sem3D/Pad", PredType::STD_I64LE, dspace_1D));
	dataset_pad.write(pad.data(), PredType::NATIVE_INT, dspace_1D, dspace_1D);

	DataSet dataset_pt(file.createDataSet("Sem3D/PillowType", PredType::STD_I64LE, dspace_1D));
	dataset_pt.write(pillow_type.data(), PredType::NATIVE_INT, dspace_1D, dspace_1D);

	// Helper macro for double datasets
	auto write_dbl_dataset = [&](const char *dset_name, const std::vector<double> &data) {
		DataSet ds(file.createDataSet(dset_name, PredType::IEEE_F64LE, dspace_1D));
		ds.write(data.data(), PredType::NATIVE_DOUBLE, dspace_1D, dspace_1D);
	};

	write_dbl_dataset("Sem3D/ScaledJacobian", v_scaledJac);
	write_dbl_dataset("Sem3D/Volume", v_volume);
	write_dbl_dataset("Sem3D/ConditionNumber", v_condNum);
	write_dbl_dataset("Sem3D/EdgeRatio", v_edgeRatio);
	write_dbl_dataset("Sem3D/Skew", v_skew);
	write_dbl_dataset("Sem3D/Shape", v_shape);
	write_dbl_dataset("Sem3D/Oddy", v_oddy);
	write_dbl_dataset("Sem3D/DiagonalRatio", v_diagRatio);
	write_dbl_dataset("Sem3D/Taper", v_taper);
	write_dbl_dataset("Sem3D/Stretch", v_stretch);
	write_dbl_dataset("Sem3D/MinFaceAngle", v_minAngle);

	if (dtcrit != NULL && (int32_t)dtcrit->size() >= n_elem) {
		write_dbl_dataset("Sem3D/DtCrit", *dtcrit);
	}

	// Centroids
	std::vector<double> centroids(3 * (size_t)n_elem);
	for (int i = 0; i < n_elem; i++) {
		double cx = 0, cy = 0, cz = 0;
		for (int j = 0; j < 8; j++) {
			int nid = connect[8 * i + j];
			cx += coords[3 * nid + 0];
			cy += coords[3 * nid + 1];
			cz += coords[3 * nid + 2];
		}
		centroids[3 * i + 0] = cx / 8.0;
		centroids[3 * i + 1] = cy / 8.0;
		centroids[3 * i + 2] = cz / 8.0;
	}
	dims2D[0] = n_elem; dims2D[1] = 3;
	DataSpace dspace_cent(2, dims2D);
	DataSet dataset_cent(file.createDataSet("Sem3D/Centroids", PredType::IEEE_F64LE, dspace_cent));
	dataset_cent.write(centroids.data(), PredType::NATIVE_DOUBLE, dspace_cent, dspace_cent);

	// XDMF Wrapper
	sprintf(filename, "%s_%04d_%04d.h5.xmf", root_name, mesh->mpi_size, mesh->mpi_rank);
	FILE *fid = fopen(filename, "w");
	sprintf(filename, "%s_%04d_%04d.h5", root_name, mesh->mpi_size, mesh->mpi_rank);

	fprintf(fid, "<?xml version=\"1.0\" ?>\n");
	fprintf(fid, "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\">\n");
	fprintf(fid, "<Xdmf Version=\"2.0\" xmlns:xi=\"http://www.w3.org/2001/XInclude\">\n");
	fprintf(fid, "<Domain>\n");
	fprintf(fid, "<Grid GridType=\"Uniform\" Name=\"main\"><Geometry Type=\"XYZ\">\n");
	fprintf(fid, "<DataItem Dimensions=\"%d 3\" Format=\"HDF\" NumberType=\"Float\" Precision=\"8\">%s:/Nodes</DataItem>\n", n_nodes, filename);
	fprintf(fid, "</Geometry>\n");
	fprintf(fid, "<Topology NumberOfElements=\"%d\" Type=\"Hexahedron\">\n", n_elem);
	fprintf(fid, "<DataItem Dimensions=\"%d 8\" Format=\"HDF\" NumberType=\"UInt\" Precision=\"8\">%s:/Sem3D/Hexa8</DataItem>\n", n_elem, filename);
	fprintf(fid, "</Topology>\n");

	auto write_xmf_attr = [&](const char *name, const char *path, const char *type = "Float", const char *number_type = "Float") {
		fprintf(fid, "<Attribute AttributeType=\"%s\" Center=\"Cell\" Dimensions=\"%d\" Name=\"%s\">\n", type, n_elem, name);
		fprintf(fid, "<DataItem Dimensions=\"%d\" Format=\"HDF\" NumberType=\"%s\" Precision=\"8\">%s:%s</DataItem>\n", n_elem, number_type, filename, path);
		fprintf(fid, "</Attribute>\n");
	};

	write_xmf_attr("Mat", "/Sem3D/Mat", "Scalar", "Int");
	write_xmf_attr("Pad", "/Sem3D/Pad", "Scalar", "Int");
	write_xmf_attr("PillowType", "/Sem3D/PillowType", "Scalar", "Int");
	write_xmf_attr("ScaledJacobian", "/Sem3D/ScaledJacobian");
	write_xmf_attr("Volume", "/Sem3D/Volume");
	write_xmf_attr("ConditionNumber", "/Sem3D/ConditionNumber");
	write_xmf_attr("EdgeRatio", "/Sem3D/EdgeRatio");
	write_xmf_attr("Skew", "/Sem3D/Skew");
	write_xmf_attr("Shape", "/Sem3D/Shape");
	write_xmf_attr("Oddy", "/Sem3D/Oddy");
	write_xmf_attr("DiagonalRatio", "/Sem3D/DiagonalRatio");
	write_xmf_attr("Taper", "/Sem3D/Taper");
	write_xmf_attr("Stretch", "/Sem3D/Stretch");
	write_xmf_attr("MinFaceAngle", "/Sem3D/MinFaceAngle");

	if (dtcrit != NULL && (int32_t)dtcrit->size() >= n_elem) {
		write_xmf_attr("DtCrit", "/Sem3D/DtCrit");
	}

	fprintf(fid, "</Grid>\n");

	// Diagnostic centroids point cloud
	fprintf(fid, "<Grid GridType=\"Uniform\" Name=\"centroids\">\n");
	fprintf(fid, "<Topology TopologyType=\"Polyvertex\" NumberOfElements=\"%d\" NodesPerElement=\"1\"/>\n", n_elem);
	fprintf(fid, "<Geometry Type=\"XYZ\">\n");
	fprintf(fid, "<DataItem Dimensions=\"%d 3\" Format=\"HDF\" NumberType=\"Float\" Precision=\"8\">%s:/Sem3D/Centroids</DataItem>\n", n_elem, filename);
	fprintf(fid, "</Geometry>\n");
	fprintf(fid, "<Attribute AttributeType=\"Scalar\" Center=\"Node\" Name=\"Mat\">\n");
	fprintf(fid, "<DataItem Dimensions=\"%d\" Format=\"HDF\" NumberType=\"Int\" Precision=\"8\">%s:/Sem3D/Mat</DataItem>\n", n_elem, filename);
	fprintf(fid, "</Attribute>\n");
	fprintf(fid, "<Attribute AttributeType=\"Scalar\" Center=\"Node\" Name=\"PillowType\">\n");
	fprintf(fid, "<DataItem Dimensions=\"%d\" Format=\"HDF\" NumberType=\"Int\" Precision=\"8\">%s:/Sem3D/PillowType</DataItem>\n", n_elem, filename);
	fprintf(fid, "</Attribute>\n");
	fprintf(fid, "<Attribute AttributeType=\"Scalar\" Center=\"Node\" Name=\"ScaledJacobian\">\n");
	fprintf(fid, "<DataItem Dimensions=\"%d\" Format=\"HDF\" NumberType=\"Float\" Precision=\"8\">%s:/Sem3D/ScaledJacobian</DataItem>\n", n_elem, filename);
	fprintf(fid, "</Attribute>\n");
	fprintf(fid, "</Grid>\n");

	fprintf(fid, "</Domain></Xdmf>");

	fclose(fid);
}
