#ifndef PILLOW_GEOMETRY_H
#define PILLOW_GEOMETRY_H

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include "mesh_geom.h"

// Geometry in VTK/H5 corner order, independent of the mesh containers.
namespace pillowgeom {
using Point = std::array<double, 3>;
using Hex = std::array<Point, 8>;
inline Point sub(const Point &a, const Point &b) {
    return {a[0]-b[0], a[1]-b[1], a[2]-b[2]};
}
inline Point cross(const Point &a, const Point &b) {
    return {a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]};
}
inline double dot(const Point &a, const Point &b) {
    return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
}

inline double jacobian(const Hex &p, double u, double v, double w) {
    static constexpr int bits[8][3] = {
        {0,0,0},{1,0,0},{1,1,0},{0,1,0},{0,0,1},{1,0,1},{1,1,1},{0,1,1}
    };
    Point a{}, b{}, c{};
    for (int i=0; i<8; ++i) {
        double x=bits[i][0] ? u : 1-u;
        double y=bits[i][1] ? v : 1-v;
        double z=bits[i][2] ? w : 1-w;
        // Translation-free evaluation avoids cancellation at large map coordinates.
        Point q=sub(p[i],p[0]);
        for (int d=0; d<3; ++d) {
            a[d] += (bits[i][0] ? 1 : -1)*y*z*q[d];
            b[d] += (bits[i][1] ? 1 : -1)*x*z*q[d];
            c[d] += (bits[i][2] ? 1 : -1)*x*y*q[d];
        }
    }
    return dot(cross(a,b),c);
}

// Columns of J = dX/d(u,v,w) on the [0,1]^3 reference cube.
inline void jacobian_matrix(const Hex &p, double u, double v, double w, double J[3][3]) {
    static constexpr int bits[8][3] = {
        {0,0,0},{1,0,0},{1,1,0},{0,1,0},{0,0,1},{1,0,1},{1,1,1},{0,1,1}
    };
    for (int r=0; r<3; ++r) for (int c=0; c<3; ++c) J[r][c]=0;
    for (int i=0; i<8; ++i) {
        const double x=bits[i][0] ? u : 1-u, y=bits[i][1] ? v : 1-v, z=bits[i][2] ? w : 1-w;
        const double g[3]={(bits[i][0] ? 1 : -1)*y*z, (bits[i][1] ? 1 : -1)*x*z, (bits[i][2] ? 1 : -1)*x*y};
        Point q=sub(p[i],p[0]);
        for (int r=0; r<3; ++r) for (int c=0; c<3; ++c) J[r][c]+=q[r]*g[c];
    }
}

// Smallest singular value of a 3x3 matrix: sqrt of the smallest eigenvalue of J^T J,
// closed form (trigonometric solution of the symmetric cubic).
inline double sigma_min3(const double J[3][3]) {
    double A[3][3];
    for (int r=0; r<3; ++r) for (int c=0; c<3; ++c) {
        A[r][c]=0; for (int k=0; k<3; ++k) A[r][c]+=J[k][r]*J[k][c];
    }
    const double p1=A[0][1]*A[0][1]+A[0][2]*A[0][2]+A[1][2]*A[1][2];
    const double q=(A[0][0]+A[1][1]+A[2][2])/3;
    const double p2=(A[0][0]-q)*(A[0][0]-q)+(A[1][1]-q)*(A[1][1]-q)+(A[2][2]-q)*(A[2][2]-q)+2*p1;
    if (p2<=1e-300) return std::sqrt(std::max(0.0,q));
    const double p=std::sqrt(p2/6);
    double B[3][3];
    for (int r=0; r<3; ++r) for (int c=0; c<3; ++c) B[r][c]=(A[r][c]-(r==c ? q : 0))/p;
    const double detB=B[0][0]*(B[1][1]*B[2][2]-B[1][2]*B[2][1])
                     -B[0][1]*(B[1][0]*B[2][2]-B[1][2]*B[2][0])
                     +B[0][2]*(B[1][0]*B[2][1]-B[1][1]*B[2][0]);
    const double phi=std::acos(std::max(-1.0,std::min(1.0,detB/2)))/3;
    const double lmax=q+2*p*std::cos(phi), lmin=q+2*p*std::cos(phi+2.0943951023931957);
    const double lmid=3*q-lmin-lmax;
    // lmin itself cancels on very anisotropic J; s1*s2*s3 = |det J| does not
    const double det=std::fabs(J[0][0]*(J[1][1]*J[2][2]-J[1][2]*J[2][1])
                              -J[0][1]*(J[1][0]*J[2][2]-J[1][2]*J[2][0])
                              +J[0][2]*(J[1][0]*J[2][1]-J[1][1]*J[2][0]));
    return (lmax>0 && lmid>0) ? det/std::sqrt(lmax*lmid) : std::sqrt(std::max(0.0,lmin));
}

// dt proxy of a trilinear hex, in metres: min of sigma_min(J) over the order-4 GLL points
// (the default gll_order). The GLL eigen-dt of a spectral element scales as sigma_min(J)/vp
// up to a per-order constant (Spearman 0.99 against DtCrit on kefalonia_small1 and hawaii
// ref3, tools/dt_report.py), at a few thousand flops instead of a power iteration.
// The GLL points and not the 27 vertex/mid points: on hawaii the 3^3 grid missed the thinnest
// spot by 2-3x on some elements and the optimizer then thinned exactly those (dt/sigma min
// ratio 8.1e-06 with 3^3, 1.8e-05 with the GLL points, median 2.3e-05 both).
// gll=false: the 3^3 grid, 5x cheaper, for the pillow inset search (same dt as the GLL
// points on kefalonia_small1 and hawaii ref3, at a fifth of the pillow time).
inline double size_metric(const Hex &p, bool gll=true) {
    static const double a=0.5*(1-std::sqrt(3.0/7.0));
    const int n=gll ? 5 : 3;
    const double g5[5]={0,a,0.5,1-a,1}, g3[3]={0,0.5,1};
    const double *g=gll ? g5 : g3;
    double s=std::numeric_limits<double>::infinity(), J[3][3];
    for (int k=0; k<n; ++k) for (int j=0; j<n; ++j) for (int i=0; i<n; ++i) {
        jacobian_matrix(p,g[i],g[j],g[k],J);
        s=std::min(s,sigma_min3(J));
    }
    return s;
}

inline double mesh_volume(const Hex &p) {
    double x[8], y[8], z[8];
    for (int k=0; k<8; ++k) { x[k]=p[k][0]; y[k]=p[k][1]; z[k]=p[k][2]; }
    return mgeom::hex_signed_volume(x,y,z);
}

// det J of a trilinear hex is tri-quadratic. A positive Bernstein hull certifies
// the ENTIRE reference cube, including its interior, not just sampled corners.
using Coefficients = std::array<double,27>;
inline bool positive_hull(const Coefficients &b, double tolerance, int depth) {
    if (*std::min_element(b.begin(),b.end()) > tolerance) return true;
    for (int k : {0,2}) for (int j : {0,2}) for (int i : {0,2})
        if (!(b[i+3*j+9*k] > tolerance)) return false;
    if (depth == 0) return false; // uncertain is never accepted as valid
    // Subdivide one axis at a time with de Casteljau; cycle through all axes.
    const int stride = depth%3 == 0 ? 1 : (depth%3 == 1 ? 3 : 9);
    Coefficients left=b, right=b;
    for (int n=0; n<27; ++n) if ((n/stride)%3 == 0) {
        double a=b[n], c=b[n+stride], e=b[n+2*stride];
        double ac=(a+c)*0.5, ce=(c+e)*0.5, mid=(ac+ce)*0.5;
        left[n+stride]=ac; left[n+2*stride]=mid;
        right[n]=mid; right[n+stride]=ce;
    }
    return positive_hull(left,tolerance,depth-1) && positive_hull(right,tolerance,depth-1);
}
inline bool valid(const Hex &p, int orientation) {
    Coefficients b;
    double scale=0;
    for (const auto &q : p) for (double t : q) if (!std::isfinite(t)) return false;
    if (!(orientation*mesh_volume(p)>0)) return false;
    for (int k=0; k<3; ++k) for (int j=0; j<3; ++j) for (int i=0; i<3; ++i) {
        double d=orientation*jacobian(p,0.5*i,0.5*j,0.5*k);
        if (!(d > 0) || !std::isfinite(d)) return false;
        b[i+3*j+9*k]=d;
        scale=std::max(scale,d);
    }
    for (int stride : {1,3,9})
        for (int n=0; n<27; ++n) if ((n/stride)%3 == 0)
            b[n+stride]=2*b[n+stride]-0.5*(b[n]+b[n+2*stride]);
    return positive_hull(b,128*std::numeric_limits<double>::epsilon()*scale,18);
}

// All faces use an outward ring in VTK/H5 order.
static constexpr int faces[6][4] = {
    {0,3,2,1},{4,5,6,7},{0,1,5,4},{1,2,6,5},{2,3,7,6},{3,0,4,7}
};

inline std::array<int,8> shell_corners(int face, int shift) {
    // 0..7 denote inner nodes, 8..15 outer nodes. Preserve the host's local
    // boundary-face index: PML extrusion uses that index to choose its axis.
    static constexpr int flip[6]={4,4,3,1,3,1};
    std::array<int,8> corners;
    for (int k=0; k<4; ++k) {
        int slot=faces[face][k], source=faces[face][(k+shift)%4];
        corners[slot]=8+source;
        corners[slot^flip[face]]=source;
    }
    return corners;
}

inline Hex shell_geometry(const Hex &outer, const Hex &inner, int face, int shift) {
    Hex shell;
    auto corners=shell_corners(face,shift);
    for (int k=0; k<8; ++k) shell[k]=corners[k]<8 ? inner[corners[k]] : outer[corners[k]-8];
    return shell;
}

inline int shell_rotation(const Hex &outer, const Hex &inner, int face, int orientation) {
    // Six-tet volume uses a chosen body diagonal and differs from integrated
    // trilinear volume on warped quads. A cyclic reference-cube rotation keeps
    // the exact same geometry, faces and Jacobian while choosing a compatible
    // diagonal for the existing mesh verifier/exporters.
    for (int shift=0; shift<4; ++shift) {
        Hex shell=shell_geometry(outer,inner,face,shift);
        if (valid(shell,orientation)) return shift;
    }
    return -1;
}

inline Point interpolate(const Hex &p, double u, double v, double w) {
    static constexpr int bits[8][3] = {
        {0,0,0},{1,0,0},{1,1,0},{0,1,0},{0,0,1},{1,0,1},{1,1,1},{0,1,1}
    };
    Point q=p[0];
    for (int k=0; k<8; ++k) {
        double weight=(bits[k][0] ? u : 1-u)*(bits[k][1] ? v : 1-v)*(bits[k][2] ? w : 1-w);
        for (int d=0; d<3; ++d) q[d]+=weight*(p[k][d]-p[0][d]);
    }
    return q;
}

// A singleton shrink set: one inner core and six closed pillow cells. Original
// outer quads are untouched, so neighbouring shrink sets cannot leave hanging
// nodes or cracks, even at material pinches or triple junctions.
// interface_faces: bit f set = host face f (faces[] order) lies on the material interface.
// Insets of the six faces (default_insets):
//  half_layer (pillowFaceInset = 1, the default): per axis, a face on the interface 0.5 and the
//    opposite one 0.25; both on the interface 0.375 each; neither 1/3 each. The shell on the
//    interface is half the host, as the half-element layer of the old extruded pillow (lattice x2,
//    buffer node at the midpoint), and the core and the opposite shell share the other half.
//  otherwise: 0.175 on interface faces, 1/3 on the others (all faces = the old uniform 0.175).
inline void default_insets(unsigned interface_faces, bool half_layer, double t[6]) {
    static constexpr int opposite[6]={1,0,4,5,2,3};   // faces[] order: w0 w1 v0 u1 v1 u0
    for (int f=0; f<6; ++f) {
        const bool me=interface_faces>>f & 1u, op=interface_faces>>opposite[f] & 1u;
        t[f]=!half_layer ? (me ? 0.175 : 1.0/3.0)
            : me && op ? 0.375 : me ? 0.5 : op ? 0.25 : 1.0/3.0;
    }
}

inline bool build(const Hex &outer, int orientation, Hex &inner, double *inset_out=nullptr,
                  unsigned interface_faces=0x3F, bool half_layer=false, double *t_out=nullptr) {
    if (!valid(outer,orientation)) return false;
    static constexpr int bits[8][3] = {
        {0,0,0},{1,0,0},{1,1,0},{0,1,0},{0,0,1},{1,0,1},{1,1,1},{0,1,1}
    };
    // face f of faces[] is the reference plane (axis, side): w=0, w=1, v=0, u=1, v=1, u=0
    static constexpr int face_axis[6]={2,2,1,0,1,0}, face_side[6]={0,1,0,1,1,0};
    double t0[6];
    default_insets(interface_faces,half_layer,t0);
    // Inset in the HOST'S reference cube, rather than toward its physical
    // centroid: a valid warped hex need not have a star-shaped physical kernel.
    // The core is an exact restriction of its trilinear map. Validate all six
    // connecting cells as well, backing off only this host if they do not fit.
    double scale=1.0;
    for (int attempt=0; attempt<20; ++attempt, scale*=0.5) {
        double lo[3], hi[3];
        for (int f=0; f<6; ++f) {
            const double t=scale*t0[f];
            if (face_side[f]) hi[face_axis[f]]=1-t; else lo[face_axis[f]]=t;
        }
        for (int k=0; k<8; ++k)
            inner[k]=interpolate(outer,bits[k][0] ? hi[0] : lo[0],
                                       bits[k][1] ? hi[1] : lo[1],
                                       bits[k][2] ? hi[2] : lo[2]);
        bool accepted=valid(inner,orientation);
        for (int f=0; f<6; ++f) {
            if (shell_rotation(outer,inner,f,orientation)<0) { accepted=false; break; }
        }
        if (accepted) {
            if (inset_out) *inset_out=scale*(*std::min_element(t0,t0+6));
            if (t_out) for (int f=0; f<6; ++f) t_out[f]=scale*t0[f];
            return true;
        }
    }
    return false;
}

// Core of the host restricted to [lo,hi]^3, from the six face insets t[f] (faces[] order).
inline void core_from_insets(const Hex &outer, const double t[6], Hex &inner) {
    static constexpr int bits[8][3] = {
        {0,0,0},{1,0,0},{1,1,0},{0,1,0},{0,0,1},{1,0,1},{1,1,1},{0,1,1}
    };
    static constexpr int face_axis[6]={2,2,1,0,1,0}, face_side[6]={0,1,0,1,1,0};
    double lo[3], hi[3];
    for (int f=0; f<6; ++f) { if (face_side[f]) hi[face_axis[f]]=1-t[f]; else lo[face_axis[f]]=t[f]; }
    for (int k=0; k<8; ++k)
        inner[k]=interpolate(outer,bits[k][0] ? hi[0] : lo[0],bits[k][1] ? hi[1] : lo[1],bits[k][2] ? hi[2] : lo[2]);
}

// Smallest size_metric over the core and the six shells; -1 if any of them is not
// certified valid (or has no compatible diagonal, see shell_rotation).
inline double children_size(const Hex &outer, const Hex &inner, int orientation, double floor_to_beat) {
    double s=size_metric(inner,false);
    Hex shells[6];
    for (int f=0; f<6; ++f) { shells[f]=shell_geometry(outer,inner,f,0); s=std::min(s,size_metric(shells[f],false)); }
    if (!(s>floor_to_beat)) return s;            // cannot win: skip the expensive certification
    if (!valid(inner,orientation)) return -1;
    for (int f=0; f<6; ++f) if (shell_rotation(outer,inner,f,orientation)<0) return -1;
    return s;
}

// Per-face insets chosen to maximise the thinnest child (the dt proxy), all seven certified.
// Hosts are independent -- the inner nodes are private to the host and the outer faces are
// untouched -- so each host can pick its own insets without affecting conformity.
// The uniform backoff of build() halves all six insets when any one shell fails, which on
// anisotropic coastline hosts left shells 1e-4 of the host's own size_metric.
// Starts from build()'s answer (so it is never worse than it), then coordinate ascent per face.
inline bool build_optimized(const Hex &outer, int orientation, Hex &inner, double *inset_out=nullptr,
                            unsigned interface_faces=0x3F, bool half_layer=false) {
    double first=0, t[6];
    if (!build(outer,orientation,inner,&first,interface_faces,half_layer,t)) return false;
    static constexpr int face_axis[6]={2,2,1,0,1,0};
    double best=children_size(outer,inner,orientation,-1);
    // A regular host at the uniform 0.175 inset gives children at 0.175 of its own metric:
    // only hosts well below that have a thin child worth searching for.
    const double host=size_metric(outer,false);
    if (best>=0.15*host) { if (inset_out) *inset_out=first; return true; }
    static constexpr double grid[]={0.5,0.375,1.0/3,0.29,0.25,0.21,0.175,0.14,0.11,0.085,0.065,0.05,
                                    0.035,0.025,0.015,0.01,0.006,0.003};
    Hex trial;
    for (int sweep=0; sweep<3; ++sweep) {
        bool improved=false;
        for (int f=0; f<6; ++f) {
            for (double g : grid) {
                if (g==t[f]) continue;
                double tt[6]; std::copy(t,t+6,tt); tt[f]=g;
                double other=0;   // the opposite face on the same axis: keep the core open
                for (int k=0; k<6; ++k) if (k!=f && face_axis[k]==face_axis[f]) other=tt[k];
                if (g+other>0.75) continue;
                core_from_insets(outer,tt,trial);
                const double s=children_size(outer,trial,orientation,best*(1+1e-9));
                if (s>best*(1+1e-9)) { best=s; std::copy(tt,tt+6,t); inner=trial; improved=true; }
            }
        }
        if (!improved) break;
    }
    // Still thin: the core as a reference-space box cannot follow a sheared host. Free its 8
    // nodes (pattern search, 6 axis moves, step from the node's distance to its host corner).
    if (best<0.15*host) {
        for (int sweep=0; sweep<4; ++sweep) {
            bool improved=false;
            for (int k=0; k<8; ++k) {
                const double r=std::sqrt(dot(sub(inner[k],outer[k]),sub(inner[k],outer[k])));
                for (double step=0.2*r; step>0.2*r/16; step*=0.5) {
                    bool moved=false;
                    for (int m=0; m<6; ++m) {
                        trial=inner;
                        trial[k][m/2]+=(m%2 ? -step : step);
                        const double s=children_size(outer,trial,orientation,best*(1+1e-9));
                        if (s>best*(1+1e-9)) { best=s; inner=trial; moved=improved=true; }
                    }
                    if (moved) break;
                }
            }
            if (!improved) break;
        }
    }
    if (inset_out) *inset_out=*std::min_element(t,t+6);
    return true;
}
} // namespace pillowgeom
#endif
