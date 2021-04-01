// =============================================================================
// PROJECT CHRONO - http://projectchrono.org
//
// Copyright (c) 2016 projectchrono.org
// All rights reserved.
//
// Use of this source code is governed by a BSD-style license that can be found
// in the LICENSE file at the top level of the distribution and at
// http://projectchrono.org/license-chrono.txt.
//
// =============================================================================
// Authors: Radu Serban
// =============================================================================
//
// Header file for ChCNarrowphaseRUtils.
// This file defines various low-level utility functions used in specialized
// pair-wise collision detection (e.g., finding the closest point on a shape to
// a specified point).
//
// =============================================================================

#pragma once

#include "chrono_multicore/math/ChMulticoreMath.h"
#include "chrono_multicore/math/matrix.h"

namespace chrono {
namespace collision {

/// @addtogroup multicore_collision
/// @{

/// This utility function returns the normal to the triangular face defined by
/// the vertices A, B, and C. The face is assumed to be non-degenerate.
/// Note that order of vertices is important!
real3 face_normal(const real3& A, const real3& B, const real3& C) {
    real3 v1 = B - A;
    real3 v2 = C - A;
    real3 n = Cross(v1, v2);
    real len = Length(n);

    return n / len;
}

/// This utility function takes the location 'P' and snaps it to the closest
/// point on the triangular face with given vertices (A, B, and C). The result
/// is returned in 'res'. Both 'P' and 'res' are assumed to be specified in
/// the same frame as the face vertices. This function returns 'true' if the
/// result is on an edge of this face and 'false' if the result is inside the
/// triangle.
/// Code from Ericson, "Real-time collision detection", 2005, pp. 141
bool snap_to_face(const real3& A, const real3& B, const real3& C, const real3& P, real3& res) {
    real3 AB = B - A;
    real3 AC = C - A;

    // Check if P in vertex region outside A
    real3 AP = P - A;
    real d1 = Dot(AB, AP);
    real d2 = Dot(AC, AP);
    if (d1 <= 0 && d2 <= 0) {
        res = A;  // barycentric coordinates (1,0,0)
        return true;
    }

    // Check if P in vertex region outside B
    real3 BP = P - B;
    real d3 = Dot(AB, BP);
    real d4 = Dot(AC, BP);
    if (d3 >= 0 && d4 <= d3) {
        res = B;  // barycentric coordinates (0,1,0)
        return true;
    }

    // Check if P in edge region of AB
    real vc = d1 * d4 - d3 * d2;
    if (vc <= 0 && d1 >= 0 && d3 <= 0) {
        // Return projection of P onto AB
        real v = d1 / (d1 - d3);
        res = A + v * AB;  // barycentric coordinates (1-v,v,0)
        return true;
    }

    // Check if P in vertex region outside C
    real3 CP = P - C;
    real d5 = Dot(AB, CP);
    real d6 = Dot(AC, CP);
    if (d6 >= 0 && d5 <= d6) {
        res = C;  // barycentric coordinates (0,0,1)
        return true;
    }

    // Check if P in edge region of AC
    real vb = d5 * d2 - d1 * d6;
    if (vb <= 0 && d2 >= 0 && d6 <= 0) {
        // Return projection of P onto AC
        real w = d2 / (d2 - d6);
        res = A + w * AC;  // barycentric coordinates (1-w,0,w)
        return true;
    }

    // Check if P in edge region of BC
    real va = d3 * d6 - d5 * d4;
    if (va <= 0 && (d4 - d3) >= 0 && (d5 - d6) >= 0) {
        // Return projection of P onto BC
        real w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        res = B + w * (C - B);  // barycentric coordinates (0,1-w,w)
        return true;
    }

    // P inside face region. Return projection of P onto face
    // barycentric coordinates (u,v,w)
    real denom = 1 / (va + vb + vc);
    real v = vb * denom;
    real w = vc * denom;
    res = A + v * AB + w * AC;  // = u*A + v*B + w*C  where  (u = 1 - v - w)
    return false;
}

/// This utility function snaps the specified location to a point on a cylinder
/// with given radius and half-length. The in/out location is assumed to be
/// specified in the frame of the cylinder (in this frame the cylinder is assumed
/// to be centered at the origin and aligned with the Y axis).  The return code
/// indicates the feature of the cylinder that caused snapping.
///   code = 0 indicates and interior point
///   code = 1 indicates snapping to one of the cylinder caps
///   code = 2 indicates snapping to the cylinder side
///   code = 3 indicates snapping to one of the cylinder edges
uint snap_to_cylinder(const real& rad, const real& hlen, real3& loc) {
    uint code = 0;

    if (loc.y > hlen) {
        code |= 1;
        loc.y = hlen;
    } else if (loc.y < -hlen) {
        code |= 1;
        loc.y = -hlen;
    }

    real d2 = loc.x * loc.x + loc.z * loc.z;

    if (d2 > rad * rad) {
        code |= 2;
        real d = Sqrt(d2);
        loc.x *= (rad / d);
        loc.z *= (rad / d);
    }

    return code;
}

/// This utility function snaps the specified location to a point on a box with
/// given half-dimensions. The in/out location is assumed to be specified in
/// the frame of the box (which is therefore assumed to be an AABB centered at
/// the origin).  The return code indicates the box axes that caused snapping.
///   - first bit (least significant) corresponds to x-axis
///   - second bit corresponds to y-axis
///   - third bit corresponds to z-axis
///
/// Therefore:
///   code = 0 indicates an interior point
///   code = 1 or code = 2 or code = 4  indicates snapping to a face
///   code = 3 or code = 5 or code = 6  indicates snapping to an edge
///   code = 7 indicates snapping to a corner
uint snap_to_box(const real3& hdims, real3& loc) {
    uint code = 0;

    if (Abs(loc.x) > hdims.x) {
        code |= 1;
        loc.x = (loc.x > 0) ? hdims.x : -hdims.x;
    }
    if (Abs(loc.y) > hdims.y) {
        code |= 2;
        loc.y = (loc.y > 0) ? hdims.y : -hdims.y;
    }
    if (Abs(loc.z) > hdims.z) {
        code |= 4;
        loc.z = (loc.z > 0) ? hdims.z : -hdims.z;
    }

    return code;
}

// ================TEST=========================
void clamp_value(real& val, const real& upper_bd, const real& lower_bd) {
    if (val < lower_bd) {
        val = lower_bd;
    }

    if (val > upper_bd) {
        val = upper_bd;
    }
}

/// TEST ============================
real3 snap_to_face_box(const real3& pt_on_box, uint code, const real3& pt_to_snap, const real3& hdims) {
    uint side = code >> 1;
    real3 A;
    real3 B;
    real3 C;
    real3 result = pt_to_snap;

    if (side == 0) {
        result.x = pt_on_box.x;
        clamp_value(result.y, hdims.y, -hdims.y);
        clamp_value(result.z, hdims.z, -hdims.z);
    } else if (side == 1) {
        result.y = pt_on_box.y;
        clamp_value(result.x, hdims.x, -hdims.x);
        clamp_value(result.z, hdims.z, -hdims.z);
    } else {
        result.z = pt_on_box.z;
        clamp_value(result.x, hdims.x, -hdims.x);
        clamp_value(result.y, hdims.y, -hdims.y);
    }

    return result;
}
/// TEST ====================
real3 snap_to_edge(const real3& pt_on_box, uint code, const real3& pt_to_snap, const real3& hdims) {
    uint axis = (~code & 7) >> 1;

    real3 pt1 = pt_on_box;
    real3 pt2 = pt_on_box;

    if (axis == 0) {
        pt1.x = -hdims.x;
        pt2.x = hdims.x;
    } else if (axis == 1) {
        pt1.y = -hdims.y;
        pt2.y = hdims.y;
    } else {
        pt1.z = -hdims.z;
        pt2.z = hdims.z;
    }

    real3 edge = pt2 - pt1;
    real t = Dot((pt_to_snap - pt1), edge) / (Length(edge) * Length(edge));

    return pt1 + t * edge;
}

// ======================================TEST===================
void get_edge_corners(const real3& pt_on_edge, uint code, real3* corners) {
    code = ((~code & 7) >> 1);

    corners[0] = pt_on_edge;
    corners[1] = pt_on_edge;

    if (code == 0) {
        corners[1].x = -corners[1].x;
    } else if (code == 1) {
        corners[1].y = -corners[1].y;
    } else {
        corners[1].z = -corners[1].z;
    }
}

// ================TEST================================
void get_face_corners(const real3& pt_on_face, uint code, real3* corners) {
    code = code >> 1;

    corners[0] = pt_on_face;
    corners[1] = pt_on_face;

    if (code == 0) {
        corners[1].y = -corners[1].y;
        corners[2] = corners[1];
        corners[2].z = -corners[1].z;
        corners[3] = pt_on_face;
        corners[3].z = corners[2].z;
    } else if (code == 1) {
        corners[1].z = -corners[1].z;
        corners[2] = corners[1];
        corners[2].x = -corners[1].x;
        corners[3] = pt_on_face;
        corners[3].x = corners[2].x;
    } else {
        corners[1].x = -corners[1].x;
        corners[2] = corners[1];
        corners[2].y = -corners[1].y;
        corners[3] = pt_on_face;
        corners[3].y = corners[2].y;
    }
}

bool point_contact_face(const real3& pt_on_face,
                        uint code,
                        const real3& pt_to_snap,
                        real3& result,
                        real& dist,
                        const real3& hdims) {
    code = code >> 1;
    if (abs(pt_to_snap.x) > hdims.x)
        return false;
    if (abs(pt_to_snap.y) > hdims.y)
        return false;
    if (abs(pt_to_snap.z) > hdims.z)
        return false;

    result = pt_to_snap;
    if (code == 0) {
        result.x = pt_on_face.x;
        dist = abs(pt_to_snap.x - pt_on_face.x);
    } else if (code == 1) {
        result.y = pt_on_face.y;
        dist = abs(pt_to_snap.y - pt_on_face.y);
    } else {
        result.z = pt_on_face.z;
        dist = abs(pt_to_snap.z - pt_on_face.z);
    }
    return dist > 1e-6;
}

// =====================TEST=======================
bool edge_contact_edge(const real3& pt_on_edge,
                       uint code,
                       const real3& pt1,
                       const real3& pt2,
                       real3& loc1,
                       real3& loc2,
                       const real3& hdims) {
    // Calculate the denominator of the solution. If it is zero, the edges are
    // parallel and there is no contact.
    code = (~code & 7) >> 1;
    real3 seg = pt2 - pt1;
    real segLen2 = Length(seg) * Length(seg);
    real denom;
    if (code == 0) {
        denom = segLen2 - seg.x * seg.x;
    } else if (code == 1) {
        denom = segLen2 - seg.y * seg.y;
    } else {
        denom = segLen2 - seg.z * seg.z;
    }

    if (denom < segLen2 * 0.025)
        return false;

    // Solve for the closest point on the edge. If this point is not on the edge
    // there is no contact.
    real3 delta = pt_on_edge - pt1;
    real tC;

    if (code == 0) {
        tC = Dot(seg, delta) - seg.x * delta.x;
    } else if (code == 1) {
        tC = Dot(seg, delta) - seg.y * delta.y;
    } else {
        tC = Dot(seg, delta) - seg.z * delta.z;
    }
    if (tC <= 0 || tC >= denom)
        return false;

    tC = tC / denom;
    loc2 = pt1 + seg * tC;

    // Calculate the closest point on this box's edge. There is contact only if
    // the point is on the edge.
    real sC;

    if (code == 0) {
        sC = seg.x * tC - delta.x;
    } else if (code == 1) {
        sC = seg.y * tC - delta.y;
    } else {
        sC = seg.z * tC - delta.z;
    }

    loc1 = pt_on_edge;
    if (code == 0) {
        loc1.x = loc1.x + sC;
        return abs(loc1.x) <= hdims.x;
    } else if (code == 1) {
        loc1.y = loc1.y + sC;
        return abs(loc1.y) <= hdims.y;
    } else {
        loc1.z = loc1.z + sC;
        return abs(loc1.z) <= hdims.z;
    }
}

/// This utility function returns the corner of a box of given dimensions that
/// if farthest in the direction 'dir', which is assumed to be given in the frame of the box.
void box_farthest_corner(const real3& hdims, const real3& dir, real3& corner) {
    corner.x = (dir.x < 0) ? hdims.x : -hdims.x;
    corner.y = (dir.y < 0) ? hdims.y : -hdims.y;
    corner.z = (dir.z < 0) ? hdims.z : -hdims.z;
}

/// This utility function returns the corner of a box of given dimensions that
/// if closest in the direction 'dir', which is assumed to be given in the frame of the box.
void box_closest_corner(const real3& hdims, const real3& dir, real3& corner) {
    corner.x = (dir.x > 0) ? hdims.x : -hdims.x;
    corner.y = (dir.y > 0) ? hdims.y : -hdims.y;
    corner.z = (dir.z > 0) ? hdims.z : -hdims.z;
}

/// This utility function returns a code that indicates the closest feature of
/// a box in the specified direction. The direction 'dir' is assumed to be
/// given in the frame of the box. The return code encodes the box axes that
/// define the closest feature:
///   - first bit (least significant) corresponds to x-axis
///   - second bit corresponds to y-axis
///   - third bit corresponds to z-axis
///
/// Therefore:
///   code = 0 indicates a degenerate direction (within a threshold)
///   code = 1 or code = 2 or code = 4  indicates a face
///   code = 3 or code = 5 or code = 6  indicates an edge
///   code = 7 indicates a corner
uint box_closest_feature(const real3& dir, const real3& hdims) {
    real threshold = 0.05;
    return ((Abs(dir.x) > threshold) << 0) | ((Abs(dir.y) > threshold) << 1) | ((Abs(dir.z) > threshold) << 2);
    /*
    real3 box_standard_vec;
    box_standard_vec.x = hdims.x;
    box_standard_vec.y = hdims.y;
    box_standard_vec.z = hdims.z;
    box_standard_vec = Normalize(box_standard_vec);

    return ((Abs(dir.x) > box_standard_vec.x) << 0) | ((Abs(dir.y) > box_standard_vec.y) << 1) |
           ((Abs(dir.z) > box_standard_vec.z) << 2);
           */
}

/// This function returns a boolean indicating whether or not a box1 with
/// dimensions hdims1 intersects a second box with the dimensions hdims2.
/// The check is performed in the local frame of box1. The transform from the
/// other box is given through 'pos' and 'rot'. If an intersection exists, the
/// direction of smallest intersection is returned in 'dir'.
///
/// This check is performed by testing 15 possible separating planes between the
/// two boxes (Gottschalk, Lin, Manocha - Siggraph96).
bool box_intersects_box(const real3& hdims1,
                        const real3& hdims2,
                        const real3& pos,
                        const quaternion& rot,
                        real3& dir,
                        real& ret_overlap) {
    Mat33 R(rot);
    Mat33 Rabs = Abs(R);
    real minOverlap = FLT_MAX;
    real overlap;
    real r1, r2;

    // 1. Test the axes of box1 (3 cases)
    // x-axis
    r2 = Dot(Rabs.row(0), hdims2);
    overlap = hdims1.x + r2 - Abs(pos.x);
    if (overlap < 0) {
        return false;
    }
    if (overlap < minOverlap) {
        dir = real3(1, 0, 0);
        minOverlap = overlap;
    }
    // y-axis
    r2 = Dot(Rabs.row(1), hdims2);
    overlap = hdims1.y + r2 - Abs(pos.y);
    if (overlap < 0) {
        return false;
    }
    if (overlap < minOverlap) {
        dir = real3(0, 1, 0);
        minOverlap = overlap;
    }
    // z-axis
    r2 = Dot(Rabs.row(2), hdims2);
    overlap = hdims1.z + r2 - Abs(pos.z);
    if (overlap < 0) {
        return false;
    }
    if (overlap < minOverlap) {
        dir = real3(0, 0, 1);
        minOverlap = overlap;
    }

    // 2. Test the axes of box2 (3 cases)
    // x-axis
    r1 = Dot(Rabs.col(0), hdims1);
    overlap = r1 + hdims2.x - Abs(Dot(R.col(0), pos));
    if (overlap < 0) {
        return false;
    }
    if (overlap < minOverlap) {
        dir = R.col(0);
        minOverlap = overlap;
    }
    // y-axis
    r1 = Dot(Rabs.col(1), hdims1);
    overlap = r1 + hdims2.y - Abs(Dot(R.col(1), pos));
    if (overlap < 0) {
        return false;
    }
    if (overlap < minOverlap) {
        dir = R.col(1);
        minOverlap = overlap;
    }
    // z-axis
    r1 = Dot(Rabs.col(2), hdims1);
    overlap = r1 + hdims2.z - Abs(Dot(R.col(2), pos));
    if (overlap < 0) {
        return false;
    }
    if (overlap < minOverlap) {
        dir = R.col(2);
        minOverlap = overlap;
    }

    // 3. Test the planes that are orthogonal (the cross-product) to pairs of axes
    // of the two boxes (9 cases)

    // case 1  A0 x B0

    real axisLen;

    r1 = hdims1.y * Rabs[2] + hdims1.z * Rabs[1];
    r2 = hdims2.y * Rabs[8] + hdims2.z * Rabs[4];
    overlap = r1 + r2 - (Abs(pos[2] * R[1] - pos[1] * R[2]));
    axisLen = Length(Cross(real3(1, 0, 0), R.col(0)));
    if (axisLen > 0.1) {
        overlap = overlap * axisLen;
        if (overlap < 0) {
            return false;
        }
        if (overlap < minOverlap) {
            dir = Normalize(Cross(real3(1, 0, 0), R.col(0)));
            minOverlap = overlap;
        }
    }

    // case 2  A0 x B1
    r1 = hdims1.y * Rabs[6] + hdims1.z * Rabs[5];
    r2 = hdims2.x * Rabs[8] + hdims2.z * Rabs[0];
    overlap = r1 + r2 - (Abs(pos[2] * R[5] - pos[1] * R[6]));
    axisLen = Length(Cross(real3(1, 0, 0), R.col(1)));
    if (axisLen > 0.1) {
        overlap = overlap * axisLen;
        if (overlap < 0) {
            return false;
        }
        if (overlap < minOverlap) {
            dir = Normalize(Cross(real3(1, 0, 0), R.col(1)));
            minOverlap = overlap;
        }
    }

    // case 3  Test axis L = A0 x B2
    r1 = hdims1.y * Rabs[10] + hdims1.z * Rabs[9];
    r2 = hdims2.x * Rabs[4] + hdims2.y * Rabs[0];
    overlap = r1 + r2 - (Abs(pos[2] * R[9] - pos[1] * R[10]));
    axisLen = Length(Cross(real3(1, 0, 0), R.col(2)));
    if (axisLen > 0.1) {
        overlap = overlap * axisLen;
        if (overlap < 0) {
            return false;
        }
        if (overlap < minOverlap) {
            dir = Normalize(Cross(real3(1, 0, 0), R.col(2)));
            minOverlap = overlap;
        }
    }

    // case 4 Test axis L = A1 x B0
    r1 = hdims1.x * Rabs[2] + hdims1.z * Rabs[0];
    r2 = hdims2.y * Rabs[9] + hdims2.z * Rabs[5];
    overlap = r1 + r2 - (Abs(pos[0] * R[2] - pos[2] * R[0]));
    axisLen = Length(Cross(real3(0, 1, 0), R.col(0)));
    if (axisLen > 0.1) {
        overlap = overlap * axisLen;
        if (overlap < 0) {
            return false;
        }
        if (overlap < minOverlap) {
            dir = Normalize(Cross(real3(0, 1, 0), R.col(0)));
            minOverlap = overlap;
        }
    }

    // case 5 Test axis L = A1 x B1
    r1 = hdims1.x * Rabs[6] + hdims1.z * Rabs[4];
    r2 = hdims2.x * Rabs[9] + hdims2.z * Rabs[1];
    overlap = r1 + r2 - (Abs(pos[0] * R[6] - pos[2] * R[4]));
    axisLen = Length(Cross(real3(0, 1, 0), R.col(1)));
    if (axisLen > 0.1) {
        overlap = overlap * axisLen;
        if (overlap < 0) {
            return false;
        }
        if (overlap < minOverlap) {
            dir = Normalize(Cross(real3(0, 1, 0), R.col(1)));
            minOverlap = overlap;
        }
    }

    // case 6 Test axis L = A1 x B2
    r1 = hdims1.x * Rabs[10] + hdims1.z * Rabs[8];
    r2 = hdims2.x * Rabs[5] + hdims2.y * Rabs[1];
    overlap = r1 + r2 - (Abs(pos[0] * R[10] - pos[2] * R[8]));
    axisLen = Length(Cross(real3(0, 1, 0), R.col(2)));
    if (axisLen > 0.1) {
        overlap = overlap * axisLen;
        if (overlap < 0) {
            return false;
        }
        if (overlap < minOverlap) {
            dir = Normalize(Cross(real3(0, 1, 0), R.col(2)));
            minOverlap = overlap;
        }
    }

    // case 7 Test axis L = A2 x B0
    r1 = hdims1.x * Rabs[1] + hdims1.y * Rabs[0];
    r2 = hdims2.y * Rabs[10] + hdims2.z * Rabs[6];
    overlap = r1 + r2 - (Abs(pos[1] * R[0] - pos[0] * R[1]));
    axisLen = Length(Cross(real3(0, 0, 1), R.col(0)));
    if (axisLen > 0.1) {
        overlap = overlap * axisLen;
        if (overlap < 0) {
            return false;
        }
        if (overlap < minOverlap) {
            dir = Normalize(Cross(real3(0, 0, 1), R.col(0)));
            minOverlap = overlap;
        }
    }

    // case 8 Test axis L = A2 x B1
    r1 = hdims1.x * Rabs[5] + hdims1.y * Rabs[4];
    r2 = hdims2.x * Rabs[10] + hdims2.z * Rabs[2];
    overlap = r1 + r2 - (Abs(pos[1] * R[4] - pos[0] * R[5]));
    axisLen = Length(Cross(real3(0, 0, 1), R.col(1)));
    if (axisLen > 0.1) {
        overlap = overlap * axisLen;
        if (overlap < 0) {
            return false;
        }
        if (overlap < minOverlap) {
            dir = Normalize(Cross(real3(0, 0, 1), R.col(1)));
            minOverlap = overlap;
        }
    }

    // case 9 Test axis L = A2 x B2
    r1 = hdims1.x * Rabs[9] + hdims1.y * Rabs[8];
    r2 = hdims2.x * Rabs[6] + hdims2.y * Rabs[2];
    overlap = r1 + r2 - (Abs(pos[1] * R[8] - pos[0] * R[9]));
    axisLen = Length(Cross(real3(0, 0, 1), R.col(2)));
    if (axisLen > 0.1) {
        overlap = overlap * axisLen;
        if (overlap < 0)
            return false;
        if (overlap < minOverlap) {
            dir = Normalize(Cross(real3(0, 0, 1), R.col(2)));
            minOverlap = overlap;
        }
    }

    ret_overlap = minOverlap;

    return true;
}

/// @} multicore_colision

}  // end namespace collision
}  // end namespace chrono
