// Auto-generated reference template definitions and orientation data for HexMesh.
// Pinned source commit: 1f146f1e7ddd85ff5f63700dc10f197387cba6d1 (origin/template_ref).
// Gate A geometry library: 11 canonical template families and 135 supported pad orientations.

#ifndef REFINEMENT_TEMPLATES_H
#define REFINEMENT_TEMPLATES_H

#include <array>
#include <vector>
#include <cstdint>
#include <string>

namespace refinement {

// Exact rational coordinates with base 24 in [0, 1]^3
struct RationalCoord3 {
    int8_t x, y, z; // numerator in 0..24, denominator is 24

    constexpr double x_double() const { return static_cast<double>(x) / 24.0; }
    constexpr double y_double() const { return static_cast<double>(y) / 24.0; }
    constexpr double z_double() const { return static_cast<double>(z) / 24.0; }

    constexpr bool operator==(const RationalCoord3 &o) const {
        return x == o.x && y == o.y && z == o.z;
    }
    constexpr bool operator!=(const RationalCoord3 &o) const { return !(*this == o); }
    constexpr bool operator<(const RationalCoord3 &o) const {
        if (x != o.x) return x < o.x;
        if (y != o.y) return y < o.y;
        return z < o.z;
    }
};

// Exact rational coordinates with base 24 in [0, 1]^2 (for boundary face signatures)
struct RationalCoord2 {
    int8_t u, v; // numerator in 0..24, denominator is 24

    constexpr double u_double() const { return static_cast<double>(u) / 24.0; }
    constexpr double v_double() const { return static_cast<double>(v) / 24.0; }

    constexpr bool operator==(const RationalCoord2 &o) const {
        return u == o.u && v == o.v;
    }
    constexpr bool operator!=(const RationalCoord2 &o) const { return !(*this == o); }
    constexpr bool operator<(const RationalCoord2 &o) const {
        if (u != o.u) return u < o.u;
        return v < o.v;
    }
};

// 3D Point
struct Point3 {
    double x = 0.0, y = 0.0, z = 0.0;
};

// Hex corners: 8 vertices in standard VTK/HexMesh order:
// 0: (0,0,0), 1: (1,0,0), 2: (1,1,0), 3: (0,1,0),
// 4: (0,0,1), 5: (1,0,1), 6: (1,1,1), 7: (0,1,1)
using HexCorners = std::array<Point3, 8>;

// Child hex connectivity: 8 indices into canonical reference vertex array
using HexConnectivity = std::array<int, 8>;

// Canonical template definition
struct CanonicalTemplate {
    int template_id; // 1..11
    int num_children;
    int num_unique_vertices;
    const RationalCoord3* vertices;     // array of unique vertices (indices 0..7 are parent corners)
    const HexConnectivity* children;    // array of child element connectivities
};

// Orientation mapping: maps a legacy pad code to (template_id, parent_order, flip_parity, edge_mask)
struct PadOrientation {
    int pad;
    int template_id;
    int rot[3];
    int sym[3];
    std::array<int, 8> parent_order; // permutation of parent corners
    bool flip_parity;                // true if reflection parity requires child corner reversal
    uint16_t edge_mask;              // 12-bit edge mask recognized by this pad
};

// 2D Face Signature representing a face's quadrilateral tiling
struct FaceSignature {
    int face_idx; // 0..5
    std::vector<std::array<RationalCoord2, 4>> quads; // sorted canonical 2D quads

    bool operator==(const FaceSignature &o) const {
        return face_idx == o.face_idx && quads == o.quads;
    }
    bool operator!=(const FaceSignature &o) const { return !(*this == o); }
};

// Validation report for Gate A inspection
struct TemplateValidationReport {
    bool valid = false;
    int child_count = 0;
    double total_volume = 0.0;
    double min_jacobian = 1e9;
    int internal_faces = 0;
    int boundary_faces = 0;
    int incompatible_faces = 0;
    std::string error_message;
};

// Template retrieval API
const CanonicalTemplate& get_canonical_template(int template_id);
const PadOrientation* get_pad_orientation(int pad);
const std::vector<int>& get_all_pad_codes();

// Permute parent corners according to RotateHex logic
std::array<int, 8> rotate_hex(const int rot[3], const int sym[3]);

// Geometry interpolation and Jacobian
Point3 trilinear_map(const Point3 &ref_pt, const HexCorners &parent);
Point3 trilinear_map(const RationalCoord3 &ref_pt, const HexCorners &parent);
double hex_jacobian_at(const HexCorners &hex, double u, double v, double w);
void hex_volume_and_jacobian_bounds(const HexCorners &hex, double &volume, double &min_jac);

// Child generation for a specific pad and parent hex
std::vector<HexCorners> generate_children_for_pad(int pad, const HexCorners &parent);

// Face signature computation
FaceSignature compute_face_signature(int template_id, int face_idx);
FaceSignature compute_face_signature_for_pad(int pad, int face_idx);

// Gate A validation function
TemplateValidationReport validate_template_pad(int pad, const HexCorners &parent);

} // namespace refinement

#endif // REFINEMENT_TEMPLATES_H
