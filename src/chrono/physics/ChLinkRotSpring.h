// =============================================================================
// PROJECT CHRONO - http://projectchrono.org
//
// Copyright (c) 2014 projectchrono.org
// All rights reserved.
//
// Use of this source code is governed by a BSD-style license that can be found
// in the LICENSE file at the top level of the distribution and at
// http://projectchrono.org/license-chrono.txt.
//
// =============================================================================
// Authors: Radu Serban
// =============================================================================

#ifndef CH_LINK_ROTSPRING_H
#define CH_LINK_ROTSPRING_H

#include "chrono/physics/ChLinkMarkers.h"

namespace chrono {

/// Class for rotational spring-damper elements with the torque specified through a callback object.
/// The torque is applied in the current direction of the relative axis of rotation.
/// While a rotational spring-damper can be associated with any pair of bodies in the system, it is
/// typically used in cases where the mechanism kinematics are such that the two bodies have a single
/// relative rotational degree of freedom (e.g, between two bodies connected through a revolute,
/// cylindrical, or screw joint).
class ChApi ChLinkRotSpring : public ChLinkMarkers {

  public:
    ChLinkRotSpring();
    ChLinkRotSpring(const ChLinkRotSpring& other);
    virtual ~ChLinkRotSpring() {}

    /// "Virtual" copy constructor (covariant return type).
    virtual ChLinkRotSpring* Clone() const override { return new ChLinkRotSpring(*this); }

    /// Get the current rotation angle about the relative rotation axis.
    double GetRotSpringAngle() const { return relAngle; }

    /// Get the current relative axis of rotation.
    const ChVector<>& GetRotSpringAxis() const { return relAxis; }

    /// Get the current relative angular speed about the common rotation axis.
    double GetRotSpringSpeed() const { return Vdot(relWvel, relAxis); }

    /// Get the current generated torque.
    double GetRotSpringTorque() const { return m_torque; }

    double GetRotSpringCoef() const { return m_spring_coef; }
    double GetRotSpringRestAngle() const { return m_rest_angle; }
    double GetRotSpringDampingCoef() const { return m_damping_coef; }

    void SetRotSpringDampingCoef(double m_d) { m_damping_coef = m_d; }
    void SetRotSpringRestAngle(double m_a) { m_rest_angle = m_a; }
    void SetRotSpringCoef(double m_s) { m_spring_coef = m_s; }

    /// Method to allow serialization of transient data to archives.
    virtual void ArchiveOUT(ChArchiveOut& marchive) override;

    /// Method to allow deserialization of transient data from archives.
    virtual void ArchiveIN(ChArchiveIn& marchive) override;

  protected:
    /// Include the rotational spring-damper custom torque.
    virtual void UpdateForces(double time) override;

    double m_torque;                              ///< resulting torque along relative axis of rotation
    double m_spring_coef;
    double m_rest_angle;
    double m_damping_coef;
};

CH_CLASS_VERSION(ChLinkRotSpring, 0)

}  // end namespace chrono

#endif
