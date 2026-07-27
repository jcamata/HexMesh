#include <cstdio>
#include <cmath>
#include <vector>
#include <sc.h>
#include <sc_containers.h>
#include "hexa.h"
#include "verify_mesh.h"
#include "mesh_geom.h"

void load_elem_xyz(hexa_tree_t *mesh, const std::vector<double> &coords,
                   int iel, double X[8], double Y[8], double Z[8]) {
	octant_t *e = (octant_t *) sc_array_index(&mesh->elements, iel);
	for (int i = 0; i < 8; i++) {
		int id = e->nodes[mgeom::H5_ORD[i]].id;
		X[i] = coords[3*id + 0];
		Y[i] = coords[3*id + 1];
		Z[i] = coords[3*id + 2];
	}
}

static double global_min_edge(hexa_tree_t *mesh, const std::vector<double> &coords) {
	static const int E[12][2] = {
		{0,1},{1,2},{2,3},{3,0},{4,5},{5,6},{6,7},{7,4},{0,4},{1,5},{2,6},{3,7}
	};
	double h = 1e300;
	for (int iel = 0; iel < mesh->elements.elem_count; iel++) {
		octant_t *e = (octant_t *) sc_array_index(&mesh->elements, iel);
		for (int k = 0; k < 12; k++) {
			int a = e->nodes[E[k][0]].id, b = e->nodes[E[k][1]].id;
			if (a == b) continue;
			double dx=coords[3*a]-coords[3*b], dy=coords[3*a+1]-coords[3*b+1], dz=coords[3*a+2]-coords[3*b+2];
			double len = std::sqrt(dx*dx+dy*dy+dz*dz);
			if (len > 1e-4 && len < h) h = len;
		}
	}
	return (h < 1e299) ? h : 1.0;
}

MeshAnalysis analyze_mesh(hexa_tree_t *mesh, const std::vector<double> &coords) {
	MeshAnalysis a;
	a.reference_sign = 1;
	a.n_inverted = 0;
	a.h_min = 0.0;
	a.global_min_sj = 1e300;
	a.global_max_sj = -1e300;

	int n = mesh->elements.elem_count;

	// Pass 1: dominant orientation = sign of the summed signed volume.
	double vol_sum = 0.0;
	std::vector<double> vols(n), sjs(n);
	for (int iel = 0; iel < n; iel++) {
		double X[8],Y[8],Z[8];
		load_elem_xyz(mesh, coords, iel, X, Y, Z);
		vols[iel] = mgeom::hex_signed_volume(X, Y, Z);
		sjs[iel]  = mgeom::hex_min_corner_sj(X, Y, Z);
		vol_sum += vols[iel];
	}
	a.reference_sign = mgeom::reference_sign(vol_sum);

	// Pass 2: classify against the reference sign.
	for (int iel = 0; iel < n; iel++) {
		double signed_sj = sjs[iel] * a.reference_sign;
		if (signed_sj < a.global_min_sj) a.global_min_sj = signed_sj;
		if (signed_sj > a.global_max_sj) a.global_max_sj = signed_sj;
		if (mgeom::is_inverted(vols[iel], sjs[iel], a.reference_sign)) {
			a.n_inverted++;
			a.inverted_ids.push_back(iel);
		}
	}
	a.h_min = global_min_edge(mesh, coords);
	return a;
}

int VerifyFacePlanarity(hexa_tree_t *mesh, const std::vector<double> &coords, double tol) {
	if (!mesh || mesh->elements.elem_count == 0 || coords.empty()) return 0;

	int n = mesh->elements.elem_count;
	double max_warp = 0.0, max_gap = 0.0, sum_warp = 0.0;
	int worst_id = -1, n_bad = 0;
	int buckets[4] = {0,0,0,0};   // <=1e-6, <=1e-4, <=1e-2, >1e-2

	std::vector<int> bad_ids;
	for (int iel = 0; iel < n; iel++) {
		double X[8],Y[8],Z[8], gap = 0.0;
		load_elem_xyz(mesh, coords, iel, X, Y, Z);
		double w = mgeom::hex_max_face_warp(X, Y, Z, &gap);
		sum_warp += w;
		if (w > max_warp) { max_warp = w; max_gap = gap; worst_id = iel; }
		if (w <= 1e-6)      buckets[0]++;
		else if (w <= 1e-4) buckets[1]++;
		else if (w <= 1e-2) buckets[2]++;
		else                buckets[3]++;
		if (w > tol) { n_bad++; if ((int)bad_ids.size() < 20) bad_ids.push_back(iel); }
	}

	printf(" =========================================================\n");
	printf("   FACE PLANARITY CHECK (tol = %.3e)\n", tol);
	printf(" =========================================================\n");
	printf("    Elements scanned:            %d\n", n);
	printf("    Max warp (gap/mean edge):    %.6e  (elem %d, gap %.6e m)\n", max_warp, worst_id, max_gap);
	printf("    Mean warp:                   %.6e\n", sum_warp / n);
	printf("    Distribution: planar(<=1e-6) %d | <=1e-4 %d | <=1e-2 %d | >1e-2 %d\n",
	       buckets[0], buckets[1], buckets[2], buckets[3]);
	printf("    Elements above tol:          %d\n", n_bad);
	if (n_bad == 0) {
		printf("    SUCCESS: every face is planar within tolerance.\n");
	} else {
		printf("    WARNING: %d elements with warped faces. First few ids:\n", n_bad);
		for (size_t i = 0; i < bad_ids.size(); i++) printf("      elem id=%d\n", bad_ids[i]);
	}
	printf(" =========================================================\n\n");
	return n_bad;
}

int VerifyMeshInversion(hexa_tree_t *mesh, const std::vector<double> *coords) {
	if (!mesh || mesh->elements.elem_count == 0 || !coords || coords->empty()) return 0;

	MeshAnalysis a = analyze_mesh(mesh, *coords);

	printf(" =========================================================\n");
	printf("   MESH INVERSION CHECK (reference sign = %+d)\n", a.reference_sign);
	printf(" =========================================================\n");
	printf("    Elements scanned:            %d\n", (int)mesh->elements.elem_count);
	printf("    Signed Scaled-Jacobian range: [%.6f, %.6f]\n", a.global_min_sj, a.global_max_sj);
	printf("    Minimum edge length (h_min):  %.6e m\n", a.h_min);
	printf("    Inverted elements:            %d\n", a.n_inverted);
	if (a.n_inverted == 0) {
		printf("    SUCCESS: 0 inverted elements (all consistent with reference orientation).\n");
	} else {
		printf("    WARNING: %d inverted elements. First few ids:\n", a.n_inverted);
		int lim = a.n_inverted < 20 ? a.n_inverted : 20;
		for (int i = 0; i < lim; i++) printf("      elem id=%d\n", a.inverted_ids[i]);
	}
	printf(" =========================================================\n\n");
	return a.n_inverted;
}
