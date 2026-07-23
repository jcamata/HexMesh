#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include "hexa.h"

// Node order for every function below: bottom face 0,1,2,3 (CCW), top face
// 4,5,6,7 (CCW), verticals 0-4,1-5,2-6,3-7. Callers MUST pass nodes already
// permuted into this order (see hexa.h comment above hex_quality_t) -- this
// file has no dependency on hexa_tree_t/octant_t on purpose, so it stays
// trivially unit-testable.

static const int HEX_EDGES[12][2] = {
	{0,1},{1,2},{2,3},{3,0}, // bottom
	{4,5},{5,6},{6,7},{7,4}, // top
	{0,4},{1,5},{2,6},{3,7}  // verticals
};

// nb[k] = the 3 neighbors of corner k along its outgoing edges
static const int CORNER_NB[8][3] = {
	{1,3,4},{2,0,5},{3,1,6},{0,2,7},{7,5,0},{4,6,1},{5,7,2},{6,4,3}
};

static void sub3(const double a[3], const double b[3], double out[3]) {
	out[0] = a[0]-b[0]; out[1] = a[1]-b[1]; out[2] = a[2]-b[2];
}
static void add3(double a[3], const double b[3]) {
	a[0] += b[0]; a[1] += b[1]; a[2] += b[2];
}
static double dot3(const double a[3], const double b[3]) {
	return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}
static double norm3(const double a[3]) {
	return sqrt(dot3(a,a));
}
static void cross3(const double a[3], const double b[3], double out[3]) {
	out[0] = a[1]*b[2] - a[2]*b[1];
	out[1] = a[2]*b[0] - a[0]*b[2];
	out[2] = a[0]*b[1] - a[1]*b[0];
}

static double det3x3(double J[3][3]) {
	return J[0][0] * (J[1][1] * J[2][2] - J[1][2] * J[2][1]) -
	       J[0][1] * (J[1][0] * J[2][2] - J[1][2] * J[2][0]) +
	       J[0][2] * (J[1][0] * J[2][1] - J[1][1] * J[2][0]);
}

// Volume via 5-tet decomposition.
static double computeHexVolume(double nodes[8][3]) {
	int tetrahedra[5][4] = {
		{0, 1, 3, 4},
		{1, 2, 3, 6},
		{1, 4, 5, 6},
		{3, 4, 6, 7},
		{1, 3, 4, 6}
	};

	double volume = 0.0;
	for (int i = 0; i < 5; i++) {
		double *n0 = nodes[tetrahedra[i][0]];
		double *n1 = nodes[tetrahedra[i][1]];
		double *n2 = nodes[tetrahedra[i][2]];
		double *n3 = nodes[tetrahedra[i][3]];

		double v1[3] = {n1[0]-n0[0], n1[1]-n0[1], n1[2]-n0[2]};
		double v2[3] = {n2[0]-n0[0], n2[1]-n0[1], n2[2]-n0[2]};
		double v3[3] = {n3[0]-n0[0], n3[1]-n0[1], n3[2]-n0[2]};

		double det = v1[0]*(v2[1]*v3[2]-v2[2]*v3[1]) -
		             v1[1]*(v2[0]*v3[2]-v2[2]*v3[0]) +
		             v1[2]*(v2[0]*v3[1]-v2[1]*v3[0]);

		volume += fabs(det) / 6.0;
	}
	return volume;
}

// Longest edge / shortest edge, over all 12 edges (1 = ideal cube).
static double hexEdgeRatio(double nodes[8][3]) {
	double lmin = 1e300, lmax = 0.0;
	for (int e = 0; e < 12; e++) {
		double d[3];
		sub3(nodes[HEX_EDGES[e][1]], nodes[HEX_EDGES[e][0]], d);
		double l = norm3(d);
		if (l < lmin) lmin = l;
		if (l > lmax) lmax = l;
	}
	if (lmin < 1e-14) return 1e10; // degenerate edge
	return lmax / lmin;
}

// Jacobian determinant and Frobenius condition number at the element
// centroid (natural coords xi=eta=zeta=0). Condition number is
// ||J||*||J^-1||/3, normalized to 1 for an ideal cube.
static void computeCentroidJacobian(double nodes[8][3], double *jacobianDet, double *condNumber) {
	double xi = 0.0, eta = 0.0, zeta = 0.0;

	double dNdxi[8][3] = {
		{-0.125*(1-eta)*(1-zeta), -0.125*(1-xi)*(1-zeta), -0.125*(1-xi)*(1-eta)},
		{ 0.125*(1-eta)*(1-zeta), -0.125*(1+xi)*(1-zeta), -0.125*(1+xi)*(1-eta)},
		{ 0.125*(1+eta)*(1-zeta),  0.125*(1+xi)*(1-zeta), -0.125*(1+xi)*(1+eta)},
		{-0.125*(1+eta)*(1-zeta),  0.125*(1-xi)*(1-zeta), -0.125*(1-xi)*(1+eta)},
		{-0.125*(1-eta)*(1+zeta), -0.125*(1-xi)*(1+zeta),  0.125*(1-xi)*(1-eta)},
		{ 0.125*(1-eta)*(1+zeta), -0.125*(1+xi)*(1+zeta),  0.125*(1+xi)*(1-eta)},
		{ 0.125*(1+eta)*(1+zeta),  0.125*(1+xi)*(1+zeta),  0.125*(1+xi)*(1+eta)},
		{-0.125*(1+eta)*(1+zeta),  0.125*(1-xi)*(1+zeta),  0.125*(1-xi)*(1+eta)}
	};

	double J[3][3] = {{0.0}};
	for (int i = 0; i < 8; i++) {
		J[0][0] += dNdxi[i][0]*nodes[i][0];
		J[0][1] += dNdxi[i][1]*nodes[i][0];
		J[0][2] += dNdxi[i][2]*nodes[i][0];
		J[1][0] += dNdxi[i][0]*nodes[i][1];
		J[1][1] += dNdxi[i][1]*nodes[i][1];
		J[1][2] += dNdxi[i][2]*nodes[i][1];
		J[2][0] += dNdxi[i][0]*nodes[i][2];
		J[2][1] += dNdxi[i][1]*nodes[i][2];
		J[2][2] += dNdxi[i][2]*nodes[i][2];
	}

	*jacobianDet = det3x3(J);

	double frobenius = 0.0;
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			frobenius += J[i][j]*J[i][j];
	frobenius = sqrt(frobenius);

	double det = *jacobianDet;
	if (fabs(det) < 1e-10) {
		*condNumber = 1e10; // degenerate
		return;
	}

	double invJ[3][3];
	invJ[0][0] = (J[1][1]*J[2][2]-J[1][2]*J[2][1]) / det;
	invJ[0][1] = -(J[0][1]*J[2][2]-J[0][2]*J[2][1]) / det;
	invJ[0][2] = (J[0][1]*J[1][2]-J[0][2]*J[1][1]) / det;
	invJ[1][0] = -(J[1][0]*J[2][2]-J[1][2]*J[2][0]) / det;
	invJ[1][1] = (J[0][0]*J[2][2]-J[0][2]*J[2][0]) / det;
	invJ[1][2] = -(J[0][0]*J[1][2]-J[0][2]*J[1][0]) / det;
	invJ[2][0] = (J[1][0]*J[2][1]-J[1][1]*J[2][0]) / det;
	invJ[2][1] = -(J[0][0]*J[2][1]-J[0][1]*J[2][0]) / det;
	invJ[2][2] = (J[0][0]*J[1][1]-J[0][1]*J[1][0]) / det;

	double frobeniusInv = 0.0;
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			frobeniusInv += invJ[i][j]*invJ[i][j];
	frobeniusInv = sqrt(frobeniusInv);

	*condNumber = frobenius * frobeniusInv / 3.0;
}

// Verdict-style skew: build the 3 principal axes by summing the 4 parallel
// edges in each logical direction, then take the largest absolute cosine
// between any pair of axes. 0 = ideal cube (orthogonal axes), 1 = worst.
static double hexSkew(double nodes[8][3]) {
	static const int xiEdges[4][2]   = {{0,1},{3,2},{4,5},{7,6}};
	static const int etaEdges[4][2]  = {{0,3},{1,2},{4,7},{5,6}};
	static const int zetaEdges[4][2] = {{0,4},{1,5},{2,6},{3,7}};

	double X1[3] = {0,0,0}, X2[3] = {0,0,0}, X3[3] = {0,0,0}, d[3];
	for (int i = 0; i < 4; i++) { sub3(nodes[xiEdges[i][1]],   nodes[xiEdges[i][0]],   d); add3(X1, d); }
	for (int i = 0; i < 4; i++) { sub3(nodes[etaEdges[i][1]],  nodes[etaEdges[i][0]],  d); add3(X2, d); }
	for (int i = 0; i < 4; i++) { sub3(nodes[zetaEdges[i][1]], nodes[zetaEdges[i][0]], d); add3(X3, d); }

	double n1 = norm3(X1), n2 = norm3(X2), n3 = norm3(X3);
	if (n1 < 1e-14 || n2 < 1e-14 || n3 < 1e-14) return 1.0; // degenerate -> worst skew

	double c12 = fabs(dot3(X1,X2)) / (n1*n2);
	double c13 = fabs(dot3(X1,X3)) / (n1*n3);
	double c23 = fabs(dot3(X2,X3)) / (n2*n3);

	double m = c12;
	if (c13 > m) m = c13;
	if (c23 > m) m = c23;
	return m;
}

// Per-corner metrics computed in a single pass over the 8 corners, each from
// the 3 edge vectors a,b,c leaving that corner (A = [a b c]):
//  - scaled Jacobian: det(A) normalized by |a||b||c|, in [-1,1] (1 = ideal)
//  - shape (Verdict):  3*det(A)^(2/3) / (|a|^2+|b|^2+|c|^2), in (0,1] (1 = ideal)
//  - Oddy: departure from orthogonality/equal edge length, >=0 (0 = ideal)
// Element-level value is the worst corner: min for scaledJacobian/shape, max for Oddy.
static void hexCornerMetrics(double nodes[8][3], double *scaledJacMin, double *shapeMin, double *oddyMax) {
	double sjMin = 1e300, shMin = 1e300, odMax = -1e300;

	for (int k = 0; k < 8; k++) {
		double a[3], b[3], c[3];
		sub3(nodes[CORNER_NB[k][0]], nodes[k], a);
		sub3(nodes[CORNER_NB[k][1]], nodes[k], b);
		sub3(nodes[CORNER_NB[k][2]], nodes[k], c);

		double bxc[3];
		cross3(b, c, bxc);
		double detA = dot3(a, bxc);

		double na = norm3(a), nb = norm3(b), nc = norm3(c);
		double prod = na*nb*nc;
		double sj = (prod > 1e-14) ? detA/prod : -1.0;
		if (sj < sjMin) sjMin = sj;

		double sumSq = na*na + nb*nb + nc*nc;
		double shp = (detA > 0.0 && sumSq > 1e-14) ? 3.0*pow(detA, 2.0/3.0)/sumSq : 0.0;
		if (shp < shMin) shMin = shp;

		double daa = dot3(a,a), dbb = dot3(b,b), dcc = dot3(c,c);
		double dab = dot3(a,b), dac = dot3(a,c), dbc = dot3(b,c);
		double frob2 = daa*daa + dbb*dbb + dcc*dcc + 2.0*(dab*dab + dac*dac + dbc*dbc);
		double tr = daa + dbb + dcc;
		double oddy = (detA > 1e-14) ? (frob2 - tr*tr/3.0) / pow(detA, 4.0/3.0) : 1e10;
		if (oddy > odMax) odMax = oddy;
	}

	*scaledJacMin = sjMin;
	*shapeMin = shMin;
	*oddyMax = odMax;
}

void hexQualityMetrics(double nodes[8][3], hex_quality_t *q) {
	if (!nodes || !q) {
		printf("Error: Null pointer passed to hexQualityMetrics.\n");
		return;
	}

	q->volume = computeHexVolume(nodes);
	q->edgeRatio = hexEdgeRatio(nodes);
	computeCentroidJacobian(nodes, &q->jacobianDet, &q->conditionNumber);
	q->skew = hexSkew(nodes);
	hexCornerMetrics(nodes, &q->scaledJacobian, &q->shape, &q->oddy);

	if (q->jacobianDet <= 0.0) {
		printf("Warning: centroid Jacobian determinant is non-positive (%.4f), indicating an invalid element.\n", q->jacobianDet);
	}
}

// ponytail: cheap assert-based self-check on a unit cube, run once from
// MeshOptimization (src/optimize_mesh.cpp) before the per-element loop.
void hexQualitySelfTest() {
	double cube[8][3] = {
		{0,0,0}, {1,0,0}, {1,1,0}, {0,1,0},
		{0,0,1}, {1,0,1}, {1,1,1}, {0,1,1}
	};
	hex_quality_t q;
	hexQualityMetrics(cube, &q);

	assert(fabs(q.volume - 1.0) < 1e-9);
	assert(fabs(q.edgeRatio - 1.0) < 1e-9);
	assert(fabs(q.jacobianDet - 1.0) < 1e-9);
	assert(fabs(q.conditionNumber - 1.0) < 1e-6);
	assert(fabs(q.skew) < 1e-9);
	assert(fabs(q.scaledJacobian - 1.0) < 1e-9);
	assert(fabs(q.shape - 1.0) < 1e-9);
	assert(fabs(q.oddy) < 1e-9);
}
