#include <vector>
#include <cstdio>
#include <sc.h>
#include <sc_containers.h>

#include "hexa.h"

// ponytail: no real node smoothing yet (Mesquite's was disabled, see
// mesquite_interface.cpp) -- this pass only computes and dumps hex quality
// metrics, one line per element, for offline inspection. Wire in an actual
// optimizer here when one is ready; MeshOptimization's signature/call site
// in main.cpp does not need to change.
void MeshOptimization(hexa_tree_t *mesh, std::vector<double> &coords, std::vector<int> material_fixed_nodes) {
	hexQualitySelfTest();

	// h5-output node order (see assign_elem_nodes = {4,5,6,7,0,1,2,3} in
	// hexa_h5.cpp); hexQualityMetrics assumes nodes already in this order.
	static const int ord[8] = {4,5,6,7,0,1,2,3};

	char fname[64];
	sprintf(fname, "HexQuality_%04d_%04d.txt", mesh->mpi_size, mesh->mpi_rank);
	FILE *fout = fopen(fname, "w");
	fprintf(fout, "# elem_id volume edge_ratio jacobian condition_number skew scaled_jacobian shape oddy\n");

	int ne = mesh->elements.elem_count;
	hex_quality_t worst{}, sum{};
	worst.edgeRatio = worst.conditionNumber = worst.oddy = -1e300;
	worst.volume = worst.jacobianDet = worst.skew = worst.scaledJacobian = worst.shape = 1e300;

	for (int iel = 0; iel < ne; iel++) {
		octant_t *e = (octant_t*) sc_array_index(&mesh->elements, iel);

		double nodes[8][3];
		for (int i = 0; i < 8; i++) {
			int id = e->nodes[ord[i]].id;
			nodes[i][0] = coords[3*id];
			nodes[i][1] = coords[3*id+1];
			nodes[i][2] = coords[3*id+2];
		}

		hex_quality_t q;
		hexQualityMetrics(nodes, &q);

		fprintf(fout, "%lld %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e\n",
			(long long)e->id, q.volume, q.edgeRatio, q.jacobianDet, q.conditionNumber,
			q.skew, q.scaledJacobian, q.shape, q.oddy);

		if (q.volume < worst.volume) worst.volume = q.volume;
		if (q.edgeRatio > worst.edgeRatio) worst.edgeRatio = q.edgeRatio;
		if (q.jacobianDet < worst.jacobianDet) worst.jacobianDet = q.jacobianDet;
		if (q.conditionNumber > worst.conditionNumber) worst.conditionNumber = q.conditionNumber;
		if (q.skew > worst.skew) worst.skew = q.skew;
		if (q.scaledJacobian < worst.scaledJacobian) worst.scaledJacobian = q.scaledJacobian;
		if (q.shape < worst.shape) worst.shape = q.shape;
		if (q.oddy > worst.oddy) worst.oddy = q.oddy;

		sum.volume += q.volume; sum.edgeRatio += q.edgeRatio; sum.jacobianDet += q.jacobianDet;
		sum.conditionNumber += q.conditionNumber; sum.skew += q.skew;
		sum.scaledJacobian += q.scaledJacobian; sum.shape += q.shape; sum.oddy += q.oddy;
	}
	fclose(fout);

	if (ne > 0) {
		printf("     Hex quality (%d elements, worst / mean): volume %.3e/%.3e  edge_ratio %.3f/%.3f  "
			"jacobian %.3e/%.3e  cond# %.3f/%.3f  skew %.3f/%.3f  scaled_jac %.3f/%.3f  shape %.3f/%.3f  oddy %.3f/%.3f\n",
			ne, worst.volume, sum.volume/ne, worst.edgeRatio, sum.edgeRatio/ne,
			worst.jacobianDet, sum.jacobianDet/ne, worst.conditionNumber, sum.conditionNumber/ne,
			worst.skew, sum.skew/ne, worst.scaledJacobian, sum.scaledJacobian/ne,
			worst.shape, sum.shape/ne, worst.oddy, sum.oddy/ne);
		fprintf(mesh->profile, "Hex quality report written to %s (%d elements)\n", fname, ne);
	}
}
