#include <cstdio>
#include <cmath>
#include <vector>
#include <algorithm>
#include <iostream>
#include <array>
#include "stability.h"
#include "verify_mesh.h"

// 8-node Hexahedron reference coordinates
static const double XI_NODE[8]   = {-1.0,  1.0,  1.0, -1.0, -1.0,  1.0,  1.0, -1.0};
static const double ETA_NODE[8]  = {-1.0, -1.0,  1.0,  1.0, -1.0, -1.0,  1.0,  1.0};
static const double ZETA_NODE[8] = {-1.0, -1.0, -1.0, -1.0,  1.0,  1.0,  1.0,  1.0};

// Legendre polynomial and its derivatives via Bonnet recurrence
static void legendre_poly(int n, double x, double &L, double &dL, double &d2L) {
	if (n == 0) { L = 1.0; dL = 0.0; d2L = 0.0; return; }
	if (n == 1) { L = x; dL = 1.0; d2L = 0.0; return; }

	double L_prev2 = 1.0;
	double L_prev1 = x;
	double L_curr = 0.0;

	for (int k = 2; k <= n; k++) {
		L_curr = ((2.0 * k - 1.0) * x * L_prev1 - (k - 1.0) * L_prev2) / (double)k;
		L_prev2 = L_prev1;
		L_prev1 = L_curr;
	}
	L = L_curr;

	// dL/dx = n/(1-x^2) * (L_{n-1}(x) - x * L_n(x))
	if (std::abs(1.0 - x*x) > 1e-14) {
		dL = (double)n / (1.0 - x*x) * (L_prev2 - x * L_curr);
		d2L = (2.0 * x * dL - (double)(n * (n + 1)) * L_curr) / (1.0 - x*x);
	} else {
		// Limit as x -> +/- 1
		double sgn = (x > 0) ? 1.0 : ((n % 2 == 0) ? 1.0 : -1.0);
		dL = 0.5 * (double)(n * (n + 1)) * sgn;
		d2L = 0.125 * (double)(n * (n + 1) * (n * (n + 1) - 2)) * ((x > 0) ? 1.0 : ((n % 2 == 0) ? -1.0 : 1.0));
	}
}

// Generate Gauss-Lobatto-Legendre nodes and weights of arbitrary order N in [1, 20]
void compute_gll_nodes_and_weights(int order, std::vector<double> &nodes, std::vector<double> &weights) {
	int N = std::max(1, std::min(order, 20));
	int n_points = N + 1;
	nodes.resize(n_points);
	weights.resize(n_points);

	if (N == 1) {
		nodes[0] = -1.0; nodes[1] = 1.0;
		weights[0] = 1.0; weights[1] = 1.0;
		return;
	}

	nodes[0] = -1.0;
	nodes[N] =  1.0;

	// Interior nodes: roots of L'_N(x) = 0 via Newton-Raphson
	for (int j = 1; j < N; j++) {
		// Initial guess from Chebyshev-Lobatto points
		double x = -std::cos((double)j * M_PI / (double)N);

		for (int it = 0; it < 100; it++) {
			double L, dL, d2L;
			legendre_poly(N, x, L, dL, d2L);
			if (std::abs(dL) < 1e-15 || std::abs(d2L) < 1e-15) break;
			double dx = dL / d2L;
			x -= dx;
			if (std::abs(dx) < 1e-15) break;
		}
		nodes[j] = x;
	}

	// Compute weights: w_i = 2 / (N*(N+1) * [L_N(x_i)]^2)
	double factor = 2.0 / (double)(N * (N + 1));
	for (int i = 0; i <= N; i++) {
		double L, dL, d2L;
		legendre_poly(N, nodes[i], L, dL, d2L);
		weights[i] = factor / (L * L);
	}
}

// D[a*n+b] = l'_b(x_a): derivada do polinomio de Lagrange do no b avaliada no no a.
// E' o operador de derivacao 1D da base espectral; o produto tensorial dele nas tres
// direcoes da' o gradiente sem montar matriz nenhuma.
static void gll_derivative_matrix(const std::vector<double> &x, std::vector<double> &D) {
	const int n = (int)x.size();
	D.assign((size_t)n * n, 0.0);
	for (int a = 0; a < n; a++) {
		for (int b = 0; b < n; b++) {
			if (a == b) {
				double sm = 0.0;
				for (int m = 0; m < n; m++) if (m != a) sm += 1.0 / (x[a] - x[m]);
				D[(size_t)a * n + b] = sm;
			} else {
				double num = 1.0, den = 1.0;
				for (int m = 0; m < n; m++) if (m != b && m != a) num *= (x[a] - x[m]);
				for (int m = 0; m < n; m++) if (m != b) den *= (x[b] - x[m]);
				D[(size_t)a * n + b] = num / den;
			}
		}
	}
}

// Compute element stability for a single 8-node hex using GLL order N
ElementStability compute_hex_gll_stability(int elem_id, int mat_id, const Material &mat,
                                           const double X[8], const double Y[8], const double Z[8],
                                           int gll_order,
                                           const std::vector<double> &gll_nodes,
                                           const std::vector<double> &gll_weights) {
	ElementStability res;
	res.elem_id = elem_id;
	res.mat_id = mat_id;
	res.dt_crit = 0.0;
	res.lambda_max = 0.0;
	res.omega_max = 0.0;
	res.dt_cfl = 0.0;
	res.min_h_gll = 0.0;
	res.jac_ratio = 0.0;
	res.valid = false;

	double vp = mat.vp > 0.0 ? mat.vp : 6000.0;
	int N = std::max(1, std::min(gll_order, 20));
	int n1d = N + 1;

	// Evaluate 3D GLL grid coordinates inside the hexahedron
	std::vector<std::vector<std::vector<std::array<double, 3>>>> gll_pos(
		n1d, std::vector<std::vector<std::array<double, 3>>>(
			n1d, std::vector<std::array<double, 3>>(n1d)));

	double min_det_J = 1e300;
	double max_det_J = -1e300;
	double min_h_gll = 1e300;
	// Metrica por ponto GLL, guardada para o limite de Irons. Ate' 2026-09 o
	// jacobiano era calculado e deitado fora, e dt_crit vinha da CFL classica.
	const size_t npt = (size_t)n1d * n1d * n1d;
	std::vector<double> detv(npt, 0.0), Jinv_all(9 * npt, 0.0);

	for (int i = 0; i < n1d; i++) {
		double xi = gll_nodes[i];
		for (int j = 0; j < n1d; j++) {
			double eta = gll_nodes[j];
			for (int k = 0; k < n1d; k++) {
				double zeta = gll_nodes[k];

				double px = 0.0, py = 0.0, pz = 0.0;
				double dxdxi = 0.0, dxdeta = 0.0, dxdzeta = 0.0;
				double dydxi = 0.0, dydeta = 0.0, dydzeta = 0.0;
				double dzdxi = 0.0, dzdeta = 0.0, dzdzeta = 0.0;

				for (int a = 0; a < 8; a++) {
					double Na = 0.125 * (1.0 + XI_NODE[a] * xi) * (1.0 + ETA_NODE[a] * eta) * (1.0 + ZETA_NODE[a] * zeta);
					px += Na * X[a];
					py += Na * Y[a];
					pz += Na * Z[a];

					double dNdxi   = 0.125 * XI_NODE[a]   * (1.0 + ETA_NODE[a] * eta) * (1.0 + ZETA_NODE[a] * zeta);
					double dNdeta  = 0.125 * ETA_NODE[a]  * (1.0 + XI_NODE[a] * xi)   * (1.0 + ZETA_NODE[a] * zeta);
					double dNdzeta = 0.125 * ZETA_NODE[a] * (1.0 + XI_NODE[a] * xi)   * (1.0 + ETA_NODE[a] * eta);

					dxdxi += dNdxi * X[a]; dxdeta += dNdeta * X[a]; dxdzeta += dNdzeta * X[a];
					dydxi += dNdxi * Y[a]; dydeta += dNdeta * Y[a]; dydzeta += dNdzeta * Y[a];
					dzdxi += dNdxi * Z[a]; dzdeta += dNdeta * Z[a]; dzdzeta += dNdzeta * Z[a];
				}

				gll_pos[i][j][k] = {px, py, pz};

				double detJ = dxdxi * (dydeta * dzdzeta - dydzeta * dzdeta)
				            - dxdeta * (dydxi * dzdzeta - dydzeta * dzdxi)
				            + dxdzeta * (dydxi * dzdeta - dydeta * dzdxi);

				if (detJ < min_det_J) min_det_J = detJ;
				if (detJ > max_det_J) max_det_J = detJ;

				const size_t p = ((size_t)i * n1d + j) * n1d + k;
				detv[p] = detJ;
				if (detJ > 1e-12) {
					const double id = 1.0 / detJ;
					double *Ji = &Jinv_all[9 * p];
					Ji[0] =  (dydeta * dzdzeta - dydzeta * dzdeta) * id;
					Ji[1] = -(dxdeta * dzdzeta - dxdzeta * dzdeta) * id;
					Ji[2] =  (dxdeta * dydzeta - dxdzeta * dydeta) * id;
					Ji[3] = -(dydxi * dzdzeta - dydzeta * dzdxi) * id;
					Ji[4] =  (dxdxi * dzdzeta - dxdzeta * dzdxi) * id;
					Ji[5] = -(dxdxi * dydzeta - dxdzeta * dydxi) * id;
					Ji[6] =  (dydxi * dzdeta - dydeta * dzdxi) * id;
					Ji[7] = -(dxdxi * dzdeta - dxdeta * dzdxi) * id;
					Ji[8] =  (dxdxi * dydeta - dxdeta * dydxi) * id;
				}
			}
		}
	}

	// Criterio de degenerescencia RELATIVO, nao absoluto.
	//
	// Ate' 2026-09 isto era `min_det_J <= 1e-12`. Mas detJ escala com o volume do
	// elemento: num hexaedro de 200 m vale ~1e6, portanto um limiar fixo de 1e-12
	// so' apanha jacobianos exatamente nulos e deixa passar elementos quase
	// degenerados. Com a CFL classica isso nunca se via -- ela devolve um valor
	// sensato para qualquer geometria porque so' mede distancias. Com o limite de
	// Irons aparecem dt de 1e-20 e contrastes de 1e17 nesta malha.
	//
	// O limiar aqui e' apenas uma GUARDA NUMERICA: abaixo dele o jacobiano nao
	// e' invertivel em dupla precisao e nao ha' operador que calcular. NAO e' um
	// filtro de qualidade -- o gerador nao decide o que e' um elemento bom.
	// Quem consome a malha filtra a jusante pelo campo JacRatio, que vai
	// exportado; medido no domo (ref=3), a razao tem mediana 0.92 e P1 = 5e-2,
	// e descartar abaixo de 1e-2 custa 0.2% dos elementos, abaixo de 0.1 custa
	// 2.2%. Essa escolha e' de quem usa, nao de quem gera.
	const double JAC_REL_TOL = 1e-14;
	if (min_det_J <= 0.0 || max_det_J <= 0.0 ||
	    min_det_J <= JAC_REL_TOL * max_det_J) {
		// Invertido ou degenerado em relacao ao proprio tamanho
		res.jac_ratio = (max_det_J > 0.0) ? (min_det_J / max_det_J) : 0.0;
		return res;
	}
	res.jac_ratio = min_det_J / max_det_J;

	// Compute minimum GLL grid spacing along coordinate lines inside the element
	for (int i = 0; i < n1d; i++) {
		for (int j = 0; j < n1d; j++) {
			for (int k = 0; k < n1d; k++) {
				if (i + 1 < n1d) {
					double dx = gll_pos[i+1][j][k][0] - gll_pos[i][j][k][0];
					double dy = gll_pos[i+1][j][k][1] - gll_pos[i][j][k][1];
					double dz = gll_pos[i+1][j][k][2] - gll_pos[i][j][k][2];
					double ds = std::sqrt(dx*dx + dy*dy + dz*dz);
					if (ds < min_h_gll) min_h_gll = ds;
				}
				if (j + 1 < n1d) {
					double dx = gll_pos[i][j+1][k][0] - gll_pos[i][j][k][0];
					double dy = gll_pos[i][j+1][k][1] - gll_pos[i][j][k][1];
					double dz = gll_pos[i][j+1][k][2] - gll_pos[i][j][k][2];
					double ds = std::sqrt(dx*dx + dy*dy + dz*dz);
					if (ds < min_h_gll) min_h_gll = ds;
				}
				if (k + 1 < n1d) {
					double dx = gll_pos[i][j][k+1][0] - gll_pos[i][j][k][0];
					double dy = gll_pos[i][j][k+1][1] - gll_pos[i][j][k][1];
					double dz = gll_pos[i][j][k+1][2] - gll_pos[i][j][k][2];
					double ds = std::sqrt(dx*dx + dy*dy + dz*dz);
					if (ds < min_h_gll) min_h_gll = ds;
				}
			}
		}
	}

	if (min_h_gll <= 1e-12) return res;

	// Spectral Element CFL time step formula: dt_crit = C_CFL * h_min_GLL / Vp
	// Standard SEM CFL limit for 3D elasticity with GLL integration
	// ---- tempo de transito, exportado ao lado do limite de Irons -------
	// h_min/vp: a menor distancia entre pontos GLL adjacentes dividida pela
	// velocidade. E' a grandeza CRUA, SEM constante CFL -- quem consome escolhe
	// a sua (0.6/sqrt(3), 0.5, o que for), e a razao DtCFL/DtCrit fica a medir
	// so' a diferenca entre geometria e operador, sem uma convencao no meio.
	//
	// E' cega ao cisalhamento: mede distancias, nao o operador. Num hexaedro
	// regular acompanha o limite de Irons a menos de um fator constante; num
	// hexaedro deformado por pillowing diverge por ordens de grandeza.
	res.valid = true;
	res.min_h_gll = min_h_gll;
	res.dt_cfl = min_h_gll / vp;

	// ---- limite de Irons: lambda_max(M_e^-1 K_e) por iteracao de potencia ---
	// Operador acustico local, aplicado matrix-free por produto tensorial:
	//   g_d = D_d u                              (gradiente de referencia)
	//   f_i = sum_k G_ik g_k,  G = w*detJ*c^2*Jinv*Jinv^T
	//   K u = sum_d D_d^T f_d
	// A massa concentrada e' M = w*detJ*rho; rho cancela em M^-1 K. Itera-se na
	// forma simetrizada M^-1/2 K M^-1/2, cuja norma e' o proprio lambda_max.
	{
		std::vector<double> Dm;
		gll_derivative_matrix(gll_nodes, Dm);
		const double c2 = vp * vp;
		std::vector<double> G(6 * npt, 0.0), sq(npt, 0.0);
		for (size_t p = 0; p < npt; p++) {
			const int ii = (int)(p / ((size_t)n1d * n1d));
			const int jj = (int)((p / n1d) % n1d);
			const int kk = (int)(p % n1d);
			const double w = gll_weights[ii] * gll_weights[jj] * gll_weights[kk];
			const double wd = w * detv[p];
			if (wd <= 0.0) return res;          // ja' filtrado por min_det_J, defensivo
			sq[p] = 1.0 / std::sqrt(wd);
			const double *Ji = &Jinv_all[9 * p];
			// G = wd * c2 * Jinv^T * Jinv  (simetrica: 00,01,02,11,12,22)
			//
			// A ordem importa e nao e' obvia. Com J[i][j] = dx_i/dxi_j, tem-se
			// grad_xi u = J^T grad_x u, logo grad_x u = J^-T grad_xi u e
			//   (grad_x u).(grad_x v) = grad_xi u^T (Jinv^T Jinv) grad_xi v.
			// Escrever Jinv*Jinv^T em vez disto da' o mesmo resultado em
			// hexaedros alinhados (onde J e' quase diagonal) e erra ate' 2x nos
			// deformados -- verificado contra eig() densa de K_e 125x125.
			double M3[3][3];
			for (int a2 = 0; a2 < 3; a2++)
				for (int b2 = 0; b2 < 3; b2++) {
					double t = 0.0;
					for (int k2 = 0; k2 < 3; k2++) t += Ji[3 * k2 + a2] * Ji[3 * k2 + b2];
					M3[a2][b2] = t * wd * c2;
				}
			G[0 * npt + p] = M3[0][0]; G[1 * npt + p] = M3[0][1]; G[2 * npt + p] = M3[0][2];
			G[3 * npt + p] = M3[1][1]; G[4 * npt + p] = M3[1][2]; G[5 * npt + p] = M3[2][2];
		}

		std::vector<double> v(npt), y(npt), t1(npt), g0(npt), g1(npt), g2(npt);
		for (size_t p = 0; p < npt; p++) v[p] = std::sin(1.0 + 0.7 * (double)p); // deterministico
		double nv = 0.0; for (double a2 : v) nv += a2 * a2;
		nv = std::sqrt(nv); for (double &a2 : v) a2 /= nv;

		double lam = 0.0;
		const int MAXIT = 200;
		for (int it = 0; it < MAXIT; it++) {
			for (size_t p = 0; p < npt; p++) t1[p] = sq[p] * v[p];
			// gradientes de referencia
			for (int i2 = 0; i2 < n1d; i2++)
			for (int j2 = 0; j2 < n1d; j2++)
			for (int k2 = 0; k2 < n1d; k2++) {
				const size_t p = ((size_t)i2 * n1d + j2) * n1d + k2;
				double a0 = 0.0, a1 = 0.0, a2 = 0.0;
				for (int x = 0; x < n1d; x++) {
					a0 += Dm[(size_t)i2 * n1d + x] * t1[((size_t)x * n1d + j2) * n1d + k2];
					a1 += Dm[(size_t)j2 * n1d + x] * t1[((size_t)i2 * n1d + x) * n1d + k2];
					a2 += Dm[(size_t)k2 * n1d + x] * t1[((size_t)i2 * n1d + j2) * n1d + x];
				}
				g0[p] = a0; g1[p] = a1; g2[p] = a2;
			}
			// f = G g
			for (size_t p = 0; p < npt; p++) {
				const double x0 = g0[p], x1 = g1[p], x2 = g2[p];
				const double f0 = G[0 * npt + p] * x0 + G[1 * npt + p] * x1 + G[2 * npt + p] * x2;
				const double f1 = G[1 * npt + p] * x0 + G[3 * npt + p] * x1 + G[4 * npt + p] * x2;
				const double f2 = G[2 * npt + p] * x0 + G[4 * npt + p] * x1 + G[5 * npt + p] * x2;
				g0[p] = f0; g1[p] = f1; g2[p] = f2;
			}
			// K u = sum_d D_d^T f_d
			for (int i2 = 0; i2 < n1d; i2++)
			for (int j2 = 0; j2 < n1d; j2++)
			for (int k2 = 0; k2 < n1d; k2++) {
				const size_t p = ((size_t)i2 * n1d + j2) * n1d + k2;
				double a0 = 0.0;
				for (int x = 0; x < n1d; x++) {
					a0 += Dm[(size_t)x * n1d + i2] * g0[((size_t)x * n1d + j2) * n1d + k2];
					a0 += Dm[(size_t)x * n1d + j2] * g1[((size_t)i2 * n1d + x) * n1d + k2];
					a0 += Dm[(size_t)x * n1d + k2] * g2[((size_t)i2 * n1d + j2) * n1d + x];
				}
				y[p] = sq[p] * a0;
			}
			double ny = 0.0; for (double a2 : y) ny += a2 * a2;
			ny = std::sqrt(ny);
			if (ny <= 0.0) return res;
			for (size_t p = 0; p < npt; p++) v[p] = y[p] / ny;
			if (it > 20 && std::abs(ny - lam) <= 1e-10 * ny) { lam = ny; break; }
			lam = ny;
		}
		if (!(lam > 0.0) || !std::isfinite(lam)) return res;
		res.lambda_max = lam;
		res.omega_max = std::sqrt(lam);
		res.dt_crit = 2.0 / res.omega_max;
	}

	return res;
}

MeshStabilityAnalysis analyze_mesh_stability(hexa_tree_t *mesh, const std::vector<double> &coords, int gll_order) {
	MeshStabilityAnalysis out;
	out.dt_crit_min = 1e300;
	out.dt_crit_mean = 0.0;
	out.dt_crit_p01 = 0.0;
	out.dt_crit_p05 = 0.0;
	out.dt_crit_p50 = 0.0;
	out.dt_crit_max = 0.0;
	out.max_freq = 0.0;
	out.worst_elem_id = -1;
	out.worst_mat_id = -1;
	out.n_unstable_elems = 0;

	if (!mesh || mesh->elements.elem_count == 0 || coords.empty()) return out;

	int ne = mesh->elements.elem_count;
	out.elem_stability.resize(ne);

	int N = std::max(1, std::min(gll_order, 20));
	std::vector<double> gll_nodes, gll_weights;
	compute_gll_nodes_and_weights(N, gll_nodes, gll_weights);

	std::vector<double> valid_dts;
	valid_dts.reserve(ne);

	double sum_dt = 0.0;

	for (int iel = 0; iel < ne; iel++) {
		octant_t *elem = (octant_t *) sc_array_index(&mesh->elements, iel);
		int mat_id = elem->n_mat;
		Material mat;
		if (mat_id >= 0 && mat_id < (int)mesh->input.materials.size()) {
			mat = mesh->input.materials[mat_id];
		} else if (!mesh->input.materials.empty()) {
			mat = mesh->input.materials[0];
		} else {
			mat.type = "S"; mat.vp = 6000.0; mat.vs = 3400.0; mat.rho = 2700.0;
		}

		double X[8], Y[8], Z[8];
		load_elem_xyz(mesh, coords, iel, X, Y, Z);

		ElementStability es = compute_hex_gll_stability(iel, mat_id, mat, X, Y, Z, N, gll_nodes, gll_weights);
		out.elem_stability[iel] = es;

		if (es.valid && es.dt_crit > 0.0) {
			valid_dts.push_back(es.dt_crit);
			sum_dt += es.dt_crit;
			if (es.dt_crit < out.dt_crit_min) {
				out.dt_crit_min = es.dt_crit;
				out.worst_elem_id = iel;
				out.worst_mat_id = mat_id;
			}
			if (es.dt_crit > out.dt_crit_max) {
				out.dt_crit_max = es.dt_crit;
			}
			if (es.omega_max > out.max_freq) {
				out.max_freq = es.omega_max;
			}
		} else {
			out.n_unstable_elems++;
		}
	}

	if (!valid_dts.empty()) {
		out.dt_crit_mean = sum_dt / valid_dts.size();
		std::sort(valid_dts.begin(), valid_dts.end());
		size_t nv = valid_dts.size();
		out.dt_crit_p01 = valid_dts[std::max((size_t)0, (size_t)(nv * 0.01))];
		out.dt_crit_p05 = valid_dts[std::max((size_t)0, (size_t)(nv * 0.05))];
		out.dt_crit_p50 = valid_dts[std::max((size_t)0, (size_t)(nv * 0.50))];
	} else {
		out.dt_crit_min = 0.0;
	}

	return out;
}

void print_stability_report(const char *label, const MeshStabilityAnalysis &analysis, int gll_order) {
	int N = std::max(1, std::min(gll_order, 20));
	printf("\n =========================================================\n");
	printf("   SPECTRAL ELEMENT GLL (N=%d) STABILITY REPORT [%s]\n", N, label ? label : "MESH");
	printf(" =========================================================\n");
	printf("    Total elements evaluated:    %zu\n", analysis.elem_stability.size());
	printf("    GLL Quadrature Order:        N = %d ((N+1)^3 = %d points/elem)\n", N, (N+1)*(N+1)*(N+1));
	printf("    Unstable/inverted elements:  %d\n", analysis.n_unstable_elems);
	printf("    Critical Time Step (SEM CFL): %.6e s  (dt_min)\n", analysis.dt_crit_min);
	printf("    Maximum Frequency (omega):   %.6e rad/s (%.2f Hz)\n", analysis.max_freq, analysis.max_freq / (2.0 * M_PI));
	printf("    Worst Element ID:            %d (Material %d)\n", analysis.worst_elem_id, analysis.worst_mat_id);
	printf("    Distribution of dt_crit:\n");
	printf("      1st  Percentile (P01):    %.6e s\n", analysis.dt_crit_p01);
	printf("      5th  Percentile (P05):    %.6e s\n", analysis.dt_crit_p05);
	printf("      50th Percentile (Median): %.6e s\n", analysis.dt_crit_p50);
	printf("      Mean dt_crit:             %.6e s\n", analysis.dt_crit_mean);
	printf("      Max  dt_crit:             %.6e s\n", analysis.dt_crit_max);
	printf(" =========================================================\n\n");
}
