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

#include "chrono/physics/ChLinkRotSpring.h"

namespace chrono {

// Register into the object factory, to enable run-time dynamic creation and persistence
CH_FACTORY_REGISTER(ChLinkRotSpring)

ChLinkRotSpring::ChLinkRotSpring() : m_torque(0), m_damping_coef(0), m_rest_angle(0), m_spring_coef(0)  {}

ChLinkRotSpring::ChLinkRotSpring(const ChLinkRotSpring& other) : ChLinkMarkers(other) {
    m_torque = other.m_torque;
    m_damping_coef = other.m_damping_coef;
    m_rest_angle = other.m_rest_angle;
    m_spring_coef = other.m_spring_coef;

}

void ChLinkRotSpring::UpdateForces(double time) {
    // Allow the base class to update itself (possibly adding its own forces).
    ChLinkMarkers::UpdateForces(time);

    // Invoke the provided functor to evaluate torque.
    double relAngle_dt = Vdot(relWvel, relAxis);

    m_torque = -m_spring_coef * (relAngle - m_rest_angle) - m_damping_coef * relAngle_dt;

    // Add to existing torque.
    C_torque += Vmul(relAxis, m_torque);
}

void ChLinkRotSpring::ArchiveOUT(ChArchiveOut& marchive) {
    // version number
    marchive.VersionWrite<ChLinkRotSpring>();

    // serialize parent class
    ChLinkMarkers::ArchiveOUT(marchive);
}

/// Method to allow de serialization of transient data from archives.
void ChLinkRotSpring::ArchiveIN(ChArchiveIn& marchive) {
    // version number
    /*int version =*/ marchive.VersionRead<ChLinkRotSpring>();

    // deserialize parent class
    ChLinkMarkers::ArchiveIN(marchive);
}

}  // end namespace chrono
