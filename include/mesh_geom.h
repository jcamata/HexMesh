#ifndef MESH_GEOM_H
#define MESH_GEOM_H

#include <cstdint>
#include <cmath>

// Pure geometry / orientation / lock helpers. NO libsc / p4est / gts includes,
// so this header is unit-testable standalone. All hex helpers assume the caller
// has ALREADY placed the 8 corner coords into X/Y/Z indexed by the reordered
// slot; callers that read elem->nodes[] must index through H5_ORD (see below).
namespace mgeom {

// h5 writer node order (assign_elem_nodes). Volume/Jacobian sign is wrong
// without this reorder applied to elem->nodes[].
static const int H5_ORD[8] = {4, 5, 6, 7, 0, 1, 2, 3};

// corner-adjacency for the 8 scaled-Jacobian corners
static const int CORNER_NB[8][3] = {
	{1,3,4},{2,0,5},{3,1,6},{0,2,7},{7,5,0},{4,6,1},{5,7,2},{6,4,3}
};

// 6-tet decomposition for exact signed volume
static const int T6[6][4] = {
	{0,1,2,6},{0,2,3,6},{0,3,7,6},{0,7,4,6},{0,4,5,6},{0,5,1,6}
};

enum { LOCK_X = 1, LOCK_Y = 2, LOCK_Z = 4 };

inline double hex_signed_volume(const double X[8], const double Y[8], const double Z[8]) {
	double v = 0.0;
	for (int t = 0; t < 6; t++) {
		const int *T = T6[t];
		double ax=X[T[1]]-X[T[0]], ay=Y[T[1]]-Y[T[0]], az=Z[T[1]]-Z[T[0]];
		double bx=X[T[2]]-X[T[0]], by=Y[T[2]]-Y[T[0]], bz=Z[T[2]]-Z[T[0]];
		double cx=X[T[3]]-X[T[0]], cy=Y[T[3]]-Y[T[0]], cz=Z[T[3]]-Z[T[0]];
		v += (ax*(by*cz-bz*cy) - ay*(bx*cz-bz*cx) + az*(bx*cy-by*cx)) / 6.0;
	}
	return v;
}

inline double hex_min_corner_sj(const double X[8], const double Y[8], const double Z[8]) {
	double min_sj = 1e300;
	for (int k = 0; k < 8; k++) {
		int i0=k, i1=CORNER_NB[k][0], i2=CORNER_NB[k][1], i3=CORNER_NB[k][2];
		double ax=X[i1]-X[i0], ay=Y[i1]-Y[i0], az=Z[i1]-Z[i0];
		double bx=X[i2]-X[i0], by=Y[i2]-Y[i0], bz=Z[i2]-Z[i0];
		double cx=X[i3]-X[i0], cy=Y[i3]-Y[i0], cz=Z[i3]-Z[i0];
		double det = ax*(by*cz-bz*cy) - ay*(bx*cz-bz*cx) + az*(bx*cy-by*cx);
		double na=std::sqrt(ax*ax+ay*ay+az*az);
		double nb=std::sqrt(bx*bx+by*by+bz*bz);
		double nc=std::sqrt(cx*cx+cy*cy+cz*cz);
		double denom = na*nb*nc;
		double sj = (denom > 1e-15) ? det/denom : 0.0;
		if (sj < min_sj) min_sj = sj;
	}
	return min_sj;
}

// the 6 faces as node rings (reordered slot order): bottom, top, 4 sides
static const int FACE[6][4] = {
	{0,1,2,3},{4,5,6,7},{0,1,5,4},{1,2,6,5},{2,3,7,6},{3,0,4,7}
};

// Out-of-plane gap of a quad, in length units: the distance between its two
// diagonal lines. Exactly 0 iff the 4 corners are coplanar. Symmetric in the
// corners, unlike "distance of p3 to plane(p0,p1,p2)".
inline double quad_planarity_gap(const double q[4][3]) {
	double d1[3] = {q[2][0]-q[0][0], q[2][1]-q[0][1], q[2][2]-q[0][2]};
	double d2[3] = {q[3][0]-q[1][0], q[3][1]-q[1][1], q[3][2]-q[1][2]};
	double n[3] = {d1[1]*d2[2]-d1[2]*d2[1], d1[2]*d2[0]-d1[0]*d2[2], d1[0]*d2[1]-d1[1]*d2[0]};
	double nn = std::sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]);
	if (nn < 1e-15) return 0.0;   // degenerate/parallel diagonals: nothing to measure
	double w[3] = {q[1][0]-q[0][0], q[1][1]-q[0][1], q[1][2]-q[0][2]};
	return std::fabs(w[0]*n[0]+w[1]*n[1]+w[2]*n[2]) / nn;
}

inline double quad_mean_edge(const double q[4][3]) {
	double s = 0.0;
	for (int k = 0; k < 4; k++) {
		int a = k, b = (k+1)%4;
		double dx=q[b][0]-q[a][0], dy=q[b][1]-q[a][1], dz=q[b][2]-q[a][2];
		s += std::sqrt(dx*dx+dy*dy+dz*dz);
	}
	return s / 4.0;
}

// Worst face non-planarity of a hex: returns the gap/mean-edge ratio
// (dimensionless); *gap_out, if given, gets the same face's gap in meters.
inline double hex_max_face_warp(const double X[8], const double Y[8], const double Z[8],
                                double *gap_out = nullptr) {
	double worst = 0.0, worst_gap = 0.0;
	for (int f = 0; f < 6; f++) {
		double q[4][3];
		for (int k = 0; k < 4; k++) {
			int i = FACE[f][k];
			q[k][0] = X[i]; q[k][1] = Y[i]; q[k][2] = Z[i];
		}
		double gap = quad_planarity_gap(q);
		double h = quad_mean_edge(q);
		double r = (h > 1e-12) ? gap / h : 0.0;
		if (r > worst) { worst = r; worst_gap = gap; }
	}
	if (gap_out) *gap_out = worst_gap;
	return worst;
}

// Two hexes sharing a quad face must extend to OPPOSITE sides of it. Project each
// element's own OTHER 4 corners' centroid onto the shared face's normal: a healthy
// pair lands on opposite sides, the SAME side means one has folded back through the
// face into the other's volume. This is real 3D interpenetration and the corner
// Jacobian does NOT see it -- both elements can be locally valid and non-inverted
// while overlapping. q = the 4 shared corners in either element's traversal order,
// ca / cb = the two elements' complement centroids.
inline bool faces_folded(const double q[4][3], const double ca[3], const double cb[3]) {
	double fc[3], d02[3], d13[3];
	for (int d = 0; d < 3; d++) {
		fc[d]  = 0.25 * (q[0][d] + q[1][d] + q[2][d] + q[3][d]);
		d02[d] = q[2][d] - q[0][d];
		d13[d] = q[3][d] - q[1][d];
	}
	double n[3] = { d02[1]*d13[2]-d02[2]*d13[1],
	                d02[2]*d13[0]-d02[0]*d13[2],
	                d02[0]*d13[1]-d02[1]*d13[0] };
	double da = 0.0, db = 0.0;
	for (int d = 0; d < 3; d++) { da += (ca[d]-fc[d])*n[d]; db += (cb[d]-fc[d])*n[d]; }
	return da * db > 0.0;
}

inline int reference_sign(double signed_vol_sum) {
	return (signed_vol_sum >= 0.0) ? 1 : -1;
}

inline bool is_inverted(double vol, double min_sj, int ref) {
	return (vol * ref <= 0.0) || (min_sj * ref <= 0.0);
}

inline uint8_t constant_axes_mask(const int cx[4], const int cy[4], const int cz[4]) {
	uint8_t m = 0;
	if (cx[0]==cx[1] && cx[1]==cx[2] && cx[2]==cx[3]) m |= LOCK_X;
	if (cy[0]==cy[1] && cy[1]==cy[2] && cy[2]==cy[3]) m |= LOCK_Y;
	if (cz[0]==cz[1] && cz[1]==cz[2] && cz[2]==cz[3]) m |= LOCK_Z;
	return m;
}

inline void apply_lock(uint8_t mask, double& dx, double& dy, double& dz) {
	if (mask & LOCK_X) dx = 0.0;
	if (mask & LOCK_Y) dy = 0.0;
	if (mask & LOCK_Z) dz = 0.0;
}

} // namespace mgeom
#endif // MESH_GEOM_H
