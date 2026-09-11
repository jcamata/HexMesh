// Patch-local elliptic untangling of the pillow layer.
//
// The pattern search in pillow.cpp moves one buffer node at a time, so it cannot fix a node that
// is shared by two hexes which spoil each other: every single-node move that helps one breaks the
// other. Measured after the node-pinch despeckle, that coupling is what the whole residual is
// made of (Argostoli ref3: 10 of the 42 buffer nodes involved carry two inverted hexes each;
// mauna_loa_small: all six inverted elements).
//
// The fix is to move a whole cluster at once, with an energy that tolerates det <= 0 so it can
// start from the tangled state -- the regularised elliptic energy of Garanzha et al., the same
// formulation used by hexsmoothing / robustPolycube, written here directly against our own mesh
// structures so nothing external is pulled in:
//
//     chi(eps,d) = (d + sqrt(eps^2 + d^2)) / 2                        (d  > 0)
//                = eps^2 / (2 (sqrt(eps^2 + d^2) - d))                (d <= 0)
//     E = sum_t w_t [ (1-theta) ||J||_F^2 / chi^(2/3) + theta (1+d^2) / chi ]
//
// chi is a smooth positive surrogate for det J, so the energy is finite and differentiable even
// where the element is inverted; shrinking eps over a few outer rounds pulls the determinants
// back to positive. Minimised by nonlinear conjugate gradient with an Armijo backtracking line
// search -- patches hold a few dozen free nodes, so there is no reason for anything heavier.
//
// It stays LOCAL: one small solve per connected cluster of inverted elements, the rest of the
// mesh frozen. Nothing global is assembled, so an MPI rank only ever touches its own patches.
//
// Every hex is split into its 8 corner tets, each compared against a right-angled reference of
// the SAME corner size -- the target is "make this corner orthogonal", which is exactly the
// scaled-Jacobian criterion used everywhere else here, and it is scale-correct: a unit reference
// cube would ask 500 m elements to become 1 m ones.
//
// Self-check: UNTANGLE_SELFCHECK=1 compares the analytic gradient against finite differences on
// the first patch and prints the worst relative error.

#include <vector>
#include <array>
#include <unordered_map>
#include <unordered_set>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <algorithm>

#include "hexa.h"
#include "mesh_geom.h"
#include "verify_mesh.h"
#include "pillow_untangle.h"

namespace {

// One corner tet: the four patch-local node indices, the inverse reference edge lengths
// (the reference frame is right-angled, so R is diagonal and R^-1 is just three numbers),
// and the reference volume used as the integration weight.
struct Tet {
	int v[4];
	double inv[3];
	double w;
};

inline double chi(double eps, double d) {
	double e2 = eps*eps;
	if (d > 0.0) return 0.5 * (d + std::sqrt(e2 + d*d));
	return 0.5 * e2 / (std::sqrt(e2 + d*d) - d);
}
inline double chi_prime(double eps, double d) {
	return 0.5 + 0.5 * d / std::sqrt(eps*eps + d*d);
}

// Energy and gradient of the whole patch. `pos` holds every patch node (locked ones included);
// the gradient is accumulated for all of them and masked by the caller.
double energy_grad(const std::vector<Tet> &tets, const std::vector<double> &pos,
                   double eps, double theta, std::vector<double> *grad, double *detmin_out)
{
	double E = 0.0, detmin = 1e300;
	if (grad) std::fill(grad->begin(), grad->end(), 0.0);

	for (size_t it = 0; it < tets.size(); it++) {
		const Tet &t = tets[it];
		const double *p0 = &pos[3*t.v[0]];
		// J = P * R^-1, with P the physical edge matrix at this corner (columns) and R the
		// right-angled reference; R diagonal makes the product a per-column scaling.
		double J[3][3];
		for (int c = 0; c < 3; c++) {
			const double *pc = &pos[3*t.v[c+1]];
			for (int r = 0; r < 3; r++) J[r][c] = (pc[r] - p0[r]) * t.inv[c];
		}
		double d = J[0][0]*(J[1][1]*J[2][2] - J[1][2]*J[2][1])
		         - J[0][1]*(J[1][0]*J[2][2] - J[1][2]*J[2][0])
		         + J[0][2]*(J[1][0]*J[2][1] - J[1][1]*J[2][0]);
		if (d < detmin) detmin = d;

		double S = 0.0;
		for (int r = 0; r < 3; r++) for (int c = 0; c < 3; c++) S += J[r][c]*J[r][c];

		double c1 = chi(eps, d), c2 = std::pow(c1, 2.0/3.0), c3 = chi_prime(eps, d);
		double f = S / c2;
		double g = (1.0 + d*d) / c1;
		E += t.w * ((1.0 - theta)*f + theta*g);
		if (!grad) continue;

		// cofactor matrix: d(det J)/dJ
		double K[3][3];
		K[0][0] =  (J[1][1]*J[2][2] - J[1][2]*J[2][1]);
		K[0][1] = -(J[1][0]*J[2][2] - J[1][2]*J[2][0]);
		K[0][2] =  (J[1][0]*J[2][1] - J[1][1]*J[2][0]);
		K[1][0] = -(J[0][1]*J[2][2] - J[0][2]*J[2][1]);
		K[1][1] =  (J[0][0]*J[2][2] - J[0][2]*J[2][0]);
		K[1][2] = -(J[0][0]*J[2][1] - J[0][1]*J[2][0]);
		K[2][0] =  (J[0][1]*J[1][2] - J[0][2]*J[1][1]);
		K[2][1] = -(J[0][0]*J[1][2] - J[0][2]*J[1][0]);
		K[2][2] =  (J[0][0]*J[1][1] - J[0][1]*J[1][0]);

		// dE/dJ = w [ (1-theta) (2J/c2 - (2/3) f c3 / c1 * K) + theta ((2d - g c3)/c1) K ]
		double kf = -(2.0/3.0) * f * c3 / c1;
		double kg = (2.0*d - g*c3) / c1;
		double dEdJ[3][3];
		for (int r = 0; r < 3; r++) for (int c = 0; c < 3; c++)
			dEdJ[r][c] = t.w * ((1.0 - theta) * (2.0*J[r][c]/c2 + kf*K[r][c]) + theta * kg*K[r][c]);

		// chain through J = P R^-1, then through P's columns to the four corner positions
		double *g0 = &(*grad)[3*t.v[0]];
		for (int c = 0; c < 3; c++) {
			double *gc = &(*grad)[3*t.v[c+1]];
			for (int r = 0; r < 3; r++) {
				double v = dEdJ[r][c] * t.inv[c];
				gc[r] += v;
				g0[r] -= v;
			}
		}
	}
	if (detmin_out) *detmin_out = detmin;
	return E;
}

// Nonlinear conjugate gradient (Polak-Ribiere with restarts) plus Armijo backtracking, over the
// free coordinates only. Small problems, so the plain method is enough and has no state to tune.
void minimize(const std::vector<Tet> &tets, std::vector<double> &pos,
              const std::vector<int> &free_dof, double eps, double theta, int maxit)
{
	const int n = (int) free_dof.size();
	if (n == 0) return;
	std::vector<double> grad(pos.size()), gprev(n, 0.0), dir(n, 0.0), pos0(pos);

	double E = energy_grad(tets, pos, eps, theta, &grad, NULL);
	for (int it = 0; it < maxit; it++) {
		double gg = 0.0, gy = 0.0, gg_prev = 0.0;
		for (int i = 0; i < n; i++) {
			double gi = grad[free_dof[i]];
			gg += gi*gi;
			gy += gi * (gi - gprev[i]);
			gg_prev += gprev[i]*gprev[i];
		}
		if (gg < 1e-24) break;
		double beta = (it == 0 || gg_prev < 1e-24) ? 0.0 : std::max(0.0, gy / gg_prev);
		double slope = 0.0;
		for (int i = 0; i < n; i++) {
			dir[i] = -grad[free_dof[i]] + beta * dir[i];
			slope += dir[i] * grad[free_dof[i]];
		}
		if (slope >= 0.0) {                       // not a descent direction: restart on -grad
			slope = 0.0;
			for (int i = 0; i < n; i++) { dir[i] = -grad[free_dof[i]]; slope += dir[i]*grad[free_dof[i]]; }
		}
		for (int i = 0; i < n; i++) gprev[i] = grad[free_dof[i]];

		// step scaled so the first trial moves the worst node by a fraction of the local size
		double dmax = 0.0;
		for (int i = 0; i < n; i++) dmax = std::max(dmax, std::fabs(dir[i]));
		if (dmax < 1e-30) break;
		double step = 1.0 / dmax;                 // first trial moves the fastest coordinate by 1
		pos0 = pos;
		bool moved = false;
		for (int ls = 0; ls < 20; ls++) {
			for (int i = 0; i < n; i++) pos[free_dof[i]] = pos0[free_dof[i]] + step * dir[i];
			double Enew = energy_grad(tets, pos, eps, theta, &grad, NULL);
			if (Enew < E + 1e-4 * step * slope) { E = Enew; moved = true; break; }
			step *= 0.5;
		}
		if (!moved) { pos = pos0; break; }
		energy_grad(tets, pos, eps, theta, &grad, NULL);
	}
}

// Put a buffer node back inside its feasible set: on the interior side of its interface node,
// thickness within [alpha_min, 1] of the placed one, direction within the cone. Same set the
// pattern search uses -- the optimiser may propose anything, the layer stays visible.
void clamp_to_cone(const BufferAnchor &a, double alpha_min, double angle_cos, double c[3])
{
	double w[3] = { c[0]-a.orig_xyz[0], c[1]-a.orig_xyz[1], c[2]-a.orig_xyz[2] };
	double l = std::sqrt(w[0]*w[0] + w[1]*w[1] + w[2]*w[2]);
	if (l < 1e-12) {
		for (int k = 0; k < 3; k++) c[k] = a.orig_xyz[k] + alpha_min*a.t_nom*a.u[k];
		return;
	}
	double wh[3] = { w[0]/l, w[1]/l, w[2]/l };
	double cs = wh[0]*a.u[0] + wh[1]*a.u[1] + wh[2]*a.u[2];
	if (cs < angle_cos) {
		double t[3] = { wh[0]-cs*a.u[0], wh[1]-cs*a.u[1], wh[2]-cs*a.u[2] };
		double tl = std::sqrt(t[0]*t[0] + t[1]*t[1] + t[2]*t[2]);
		double sn = std::sqrt(std::max(0.0, 1.0 - angle_cos*angle_cos));
		for (int k = 0; k < 3; k++)
			wh[k] = angle_cos*a.u[k] + (tl > 1e-12 ? sn*t[k]/tl : 0.0);
	}
	if (l > a.t_nom)            l = a.t_nom;
	if (l < alpha_min*a.t_nom)  l = alpha_min*a.t_nom;
	for (int k = 0; k < 3; k++) c[k] = a.orig_xyz[k] + l*wh[k];
}

bool elem_inverted(hexa_tree_t *mesh, const std::vector<double> &coords, int ie, int ref)
{
	double X[8], Y[8], Z[8];
	load_elem_xyz(mesh, coords, ie, X, Y, Z);
	return mgeom::is_inverted(mgeom::hex_signed_volume(X, Y, Z),
	                          mgeom::hex_min_corner_sj(X, Y, Z), ref);
}

// Analytic gradient against central differences, on the patch at hand. Prints the worst relative
// error; anything above ~1e-5 means the chain rule above is wrong.
void selfcheck(const std::vector<Tet> &tets, std::vector<double> pos,
               const std::vector<int> &free_dof, double eps, double theta)
{
	std::vector<double> grad(pos.size());
	energy_grad(tets, pos, eps, theta, &grad, NULL);
	double worst = 0.0;
	int n = (int) std::min<size_t>(free_dof.size(), 30);
	for (int i = 0; i < n; i++) {
		int k = free_dof[i];
		double h = std::max(1e-4, 1e-6 * std::fabs(pos[k]));
		double keep = pos[k];
		pos[k] = keep + h; double Ep = energy_grad(tets, pos, eps, theta, NULL, NULL);
		pos[k] = keep - h; double Em = energy_grad(tets, pos, eps, theta, NULL, NULL);
		pos[k] = keep;
		double fd = (Ep - Em) / (2*h);
		double den = std::max(1.0, std::max(std::fabs(fd), std::fabs(grad[k])));
		worst = std::max(worst, std::fabs(fd - grad[k]) / den);
	}
	printf("    [untangle selfcheck] worst relative gradient error over %d dofs: %.3e\n", n, worst);
}

} // namespace

int EllipticPatchUntangle(hexa_tree_t *mesh, std::vector<double> &coords, int ref,
                          const std::unordered_map<int, BufferAnchor> &anchors,
                          double alpha_min, double angle_cos, int *n_patches_out)
{
	const int ne = (int) mesh->elements.elem_count;
	std::vector<char> bad(ne, 0);
	int nbad = 0;
	for (int ie = 0; ie < ne; ie++)
		if (elem_inverted(mesh, coords, ie, ref)) { bad[ie] = 1; nbad++; }
	if (nbad == 0) { if (n_patches_out) *n_patches_out = 0; return 0; }

	// free node -> the elements that use it (only buffer nodes may move)
	std::unordered_map<int, std::vector<int> > n2e;
	for (int ie = 0; ie < ne; ie++) {
		octant_t *e = (octant_t *) sc_array_index(&mesh->elements, ie);
		for (int k = 0; k < 8; k++) {
			int id = e->nodes[k].id;
			if (anchors.count(id)) n2e[id].push_back(ie);
		}
	}

	// Patches: inverted elements joined transitively through a shared free node. The element list
	// then takes in everything incident to those nodes -- the neighbours are what keep the solve
	// honest; they are present but locked.
	std::vector<char> seen(ne, 0);
	int n_patches = 0, n_fixed = 0;
	bool do_selfcheck = (getenv("UNTANGLE_SELFCHECK") != NULL);

	for (int seed = 0; seed < ne; seed++) {
		if (!bad[seed] || seen[seed]) continue;

		std::unordered_set<int> pfree, pelem;
		std::vector<int> stack(1, seed);
		seen[seed] = 1;
		while (!stack.empty()) {
			int ie = stack.back(); stack.pop_back();
			pelem.insert(ie);
			octant_t *e = (octant_t *) sc_array_index(&mesh->elements, ie);
			for (int k = 0; k < 8; k++) {
				int id = e->nodes[k].id;
				if (!anchors.count(id)) continue;
				if (!pfree.insert(id).second) continue;
				for (int je : n2e[id]) {
					pelem.insert(je);
					if (bad[je] && !seen[je]) { seen[je] = 1; stack.push_back(je); }
				}
			}
		}
		if (pfree.empty()) continue;

		// --- local numbering ----------------------------------------------------------------
		std::unordered_map<int,int> g2l;
		std::vector<int> l2g;
		std::vector<int> elems(pelem.begin(), pelem.end());
		for (int ie : elems) {
			octant_t *e = (octant_t *) sc_array_index(&mesh->elements, ie);
			for (int k = 0; k < 8; k++) {
				int id = e->nodes[k].id;
				if (g2l.emplace(id, (int)l2g.size()).second) l2g.push_back(id);
			}
		}
		const int nv = (int) l2g.size();
		std::vector<double> pos(3*nv);
		for (int i = 0; i < nv; i++)
			for (int d = 0; d < 3; d++) pos[3*i+d] = coords[3*l2g[i]+d];
		std::vector<int> free_dof;
		free_dof.reserve(3*pfree.size());
		for (int id : pfree) {
			int li = g2l[id];
			for (int d = 0; d < 3; d++) free_dof.push_back(3*li+d);
		}

		// --- corner tets with a right-angled reference of the same size ----------------------
		std::vector<Tet> tets;
		tets.reserve(8*elems.size());
		double X[8], Y[8], Z[8];
		for (int ie : elems) {
			octant_t *e = (octant_t *) sc_array_index(&mesh->elements, ie);
			load_elem_xyz(mesh, coords, ie, X, Y, Z);
			for (int k = 0; k < 8; k++) {
				const int a = mgeom::CORNER_NB[k][0];
				const int b = mgeom::CORNER_NB[k][1];
				const int c = mgeom::CORNER_NB[k][2];
				double la = std::sqrt((X[a]-X[k])*(X[a]-X[k]) + (Y[a]-Y[k])*(Y[a]-Y[k]) + (Z[a]-Z[k])*(Z[a]-Z[k]));
				double lb = std::sqrt((X[b]-X[k])*(X[b]-X[k]) + (Y[b]-Y[k])*(Y[b]-Y[k]) + (Z[b]-Z[k])*(Z[b]-Z[k]));
				double lc = std::sqrt((X[c]-X[k])*(X[c]-X[k]) + (Y[c]-Y[k])*(Y[c]-Y[k]) + (Z[c]-Z[k])*(Z[c]-Z[k]));
				if (la < 1e-9) la = 1e-9;
				if (lb < 1e-9) lb = 1e-9;
				if (lc < 1e-9) lc = 1e-9;
				Tet t;
				// a mirrored mesh (ref = -1) is made right-handed by swapping two corners, so
				// the reference frame stays positively oriented in either case
				t.v[0] = g2l[e->nodes[mgeom::H5_ORD[k]].id];
				t.v[1] = g2l[e->nodes[mgeom::H5_ORD[a]].id];
				t.v[2] = g2l[e->nodes[mgeom::H5_ORD[ref > 0 ? b : c]].id];
				t.v[3] = g2l[e->nodes[mgeom::H5_ORD[ref > 0 ? c : b]].id];
				double l1 = ref > 0 ? lb : lc, l2 = ref > 0 ? lc : lb;
				t.inv[0] = 1.0/la; t.inv[1] = 1.0/l1; t.inv[2] = 1.0/l2;
				t.w = la * l1 * l2 / 6.0;
				tets.push_back(t);
			}
		}

		// --- solve, shrinking the regularisation ---------------------------------------------
		std::vector<double> saved(pos);
		const double theta = 1e-3;
		double detmin = 0.0;
		energy_grad(tets, pos, 1.0, theta, NULL, &detmin);
		if (do_selfcheck) {
			double e0 = detmin > 0 ? 0.1 : std::sqrt(0.01 + 0.04*detmin*detmin);
			selfcheck(tets, pos, free_dof, e0, theta);
			do_selfcheck = false;
		}
		for (int outer = 0; outer < 8; outer++) {
			energy_grad(tets, pos, 1.0, theta, NULL, &detmin);
			// eps large enough to smooth over the inverted tets, then decayed toward zero
			double e0 = std::pow(0.5, outer);
			double eps = detmin > 0.0 ? 0.5*e0 : std::sqrt(e0*e0 + 0.04*detmin*detmin);
			minimize(tets, pos, free_dof, eps, theta, 60);
		}

		// --- accept only if it actually removes inverted elements -----------------------------
		int before = 0;
		for (int ie : elems) if (elem_inverted(mesh, coords, ie, ref)) before++;
		for (int id : pfree) {
			int li = g2l[id];
			double c3[3] = { pos[3*li+0], pos[3*li+1], pos[3*li+2] };
			clamp_to_cone(anchors.at(id), alpha_min, angle_cos, c3);
			for (int d = 0; d < 3; d++) coords[3*id+d] = c3[d];
		}
		int after = 0;
		for (int ie : elems) if (elem_inverted(mesh, coords, ie, ref)) after++;
		if (after >= before) {                       // no gain: put the patch back untouched
			for (int id : pfree) {
				int li = g2l[id];
				for (int d = 0; d < 3; d++) coords[3*id+d] = saved[3*li+d];
			}
		} else {
			n_fixed += before - after;
		}
		n_patches++;
	}

	if (n_patches_out) *n_patches_out = n_patches;
	return n_fixed;
}
