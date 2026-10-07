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
inline bool build(const Hex &outer, int orientation, Hex &inner, double *inset_out=nullptr) {
    if (!valid(outer,orientation)) return false;
    static constexpr int bits[8][3] = {
        {0,0,0},{1,0,0},{1,1,0},{0,1,0},{0,0,1},{1,0,1},{1,1,1},{0,1,1}
    };
    // Inset in the HOST'S reference cube, rather than toward its physical
    // centroid: a valid warped hex need not have a star-shaped physical kernel.
    // The core is an exact restriction of its trilinear map. Validate all six
    // connecting cells as well, backing off only this host if they do not fit.
    double inset=0.175;
    for (int attempt=0; attempt<20; ++attempt, inset*=0.5) {
        for (int k=0; k<8; ++k)
            inner[k]=interpolate(outer,bits[k][0] ? 1-inset : inset,
                                       bits[k][1] ? 1-inset : inset,
                                       bits[k][2] ? 1-inset : inset);
        bool accepted=valid(inner,orientation);
        for (int f=0; f<6; ++f) {
            if (shell_rotation(outer,inner,f,orientation)<0) { accepted=false; break; }
        }
        if (accepted) {
            if (inset_out) *inset_out=inset;
            return true;
        }
    }
    return false;
}
} // namespace pillowgeom
#endif
