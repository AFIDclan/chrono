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

#include <cstdio>
#include <vector>
#include <cmath>

#include "chrono_multicore/physics/ChSystemMulticore.h"

#include "chrono/ChConfig.h"
#include "chrono/utils/ChUtilsCreators.h"
#include "chrono/utils/ChUtilsInputOutput.h"

#ifdef CHRONO_OPENGL
    #include "chrono_opengl/ChOpenGLWindow.h"
#endif

using namespace chrono;
using namespace chrono::collision;

// Tilt angle (about global Y axis) of the container.
double tilt_angle = 0;

// Number of balls: (2 * count_X + 1) * (2 * count_Y + 1)
int count_X = 2;
int count_Y = 2;

// Material properties (same on bin and balls)
float Y = 2e6f;
float mu = 0.4f;
float cr = 0.4f;

// -----------------------------------------------------------------------------
// Create a bin consisting of five boxes attached to the ground.
// -----------------------------------------------------------------------------
void AddContainer(ChSystemMulticoreSMC* sys) {
    // IDs for the two bodies
    int binId = -200;

    // Create a common material
    auto mat = chrono_types::make_shared<ChMaterialSurfaceSMC>();
    mat->SetYoungModulus(Y);
    mat->SetFriction(mu);
    mat->SetRestitution(cr);

    // Create the containing bin (4 x 4 x 1)
    auto bin = chrono_types::make_shared<ChBody>(chrono_types::make_shared<ChCollisionModelMulticore>());
    bin->SetIdentifier(binId);
    bin->SetMass(1);
    bin->SetPos(ChVector<>(0, 0, 0));
    bin->SetRot(Q_from_AngY(tilt_angle));
    bin->SetCollide(true);
    bin->SetBodyFixed(true);

    ChVector<> hdim(2, 2, 0.5);
    double hthick = 0.1;

    bin->GetCollisionModel()->ClearModel();
    utils::AddBoxGeometry(bin.get(), mat, ChVector<>(hdim.x(), hdim.y(), hthick), ChVector<>(0, 0, -hthick));
    // utils::AddBoxGeometry(bin.get(), mat, ChVector<>(hthick, hdim.y(), hdim.z()),
    //                      ChVector<>(-hdim.x() - hthick, 0, hdim.z()));
    // utils::AddBoxGeometry(bin.get(), mat, ChVector<>(hthick, hdim.y(), hdim.z()),
    //                      ChVector<>(hdim.x() + hthick, 0, hdim.z()));
    // utils::AddBoxGeometry(bin.get(), mat, ChVector<>(hdim.x(), hthick, hdim.z()),
    //                      ChVector<>(0, -hdim.y() - hthick, hdim.z()));
    // utils::AddBoxGeometry(bin.get(), mat, ChVector<>(hdim.x(), hthick, hdim.z()),
    //                      ChVector<>(0, hdim.y() + hthick, hdim.z()));
    bin->GetCollisionModel()->BuildModel();

    sys->AddBody(bin);
}

// -----------------------------------------------------------------------------
// Create the falling spherical objects in a uniform rectangular grid.
// -----------------------------------------------------------------------------
void AddCollisionBox(ChSystemMulticore* sys) {
    // Common material
    auto boxMat = chrono_types::make_shared<ChMaterialSurfaceSMC>();
    boxMat->SetYoungModulus(Y);
    boxMat->SetFriction(mu);
    boxMat->SetRestitution(cr);
    boxMat->SetAdhesion(0);  // Magnitude of the adhesion in Constant adhesion model

    int boxId = 0;
    double mass = 10;
    ChVector<> size = ChVector<>(0.15, 0.1, 0.1);
    ChVector<> inertia = 1;

    for (int i = 0; i < 1; i++) {
        ChVector<> pos(0.0 + 0.4 * i, 1.6 + 0.4 * i, 0.3 + 0.3 * i);
        ChQuaternion<> rot = Q_from_Euler123(ChVector<>(0, 0, 0));

        auto box = chrono_types::make_shared<ChBody>(chrono_types::make_shared<ChCollisionModelMulticore>());
        box->SetIdentifier(boxId++);
        box->SetMass(mass);
        box->SetInertiaXX(inertia);
        box->SetPos(pos);
        // box->SetRot(rot);
        // ChQuaternion<> rot = Q_from_Euler123(ChVector<>(0, CH_C_PI / 6, 0));
        box->SetRot(rot);
        if (i == 0) {
            box->SetBodyFixed(false);
        }
        if (i == 1) {
            box->SetBodyFixed(false);
        }

        box->SetCollide(true);

        box->GetCollisionModel()->ClearModel();
        utils::AddBoxGeometry(box.get(), boxMat, size);
        box->GetCollisionModel()->BuildModel();

        sys->AddBody(box);
    }
}

// -----------------------------------------------------------------------------
// Create the system, specify simulation parameters, and run simulation loop.
// -----------------------------------------------------------------------------
int main(int argc, char* argv[]) {
    GetLog() << "Copyright (c) 2017 projectchrono.org\nChrono version: " << CHRONO_VERSION << "\n\n";

    // Simulation parameters
    // ---------------------

    double gravity = 9.81;
    double time_step = 1e-4;
    double time_end = 20;

    double out_fps = 50;

    uint max_iteration = 500;
    real tolerance = 1e-3;

    // Create system
    // -------------

    ChSystemMulticoreSMC msystem;

    // Set number of threads
    msystem.SetNumThreads(8);

    // Set gravitational acceleration
    msystem.Set_G_acc(ChVector<>(0, 0, -gravity));

    // Set solver parameters
    msystem.GetSettings()->solver.max_iteration_bilateral = max_iteration;
    msystem.GetSettings()->solver.tolerance = tolerance;

    msystem.GetSettings()->collision.narrowphase_algorithm = NarrowPhaseType::NARROWPHASE_HYBRID_MPR;
    msystem.GetSettings()->collision.bins_per_axis = vec3(10, 10, 10);

    // The following two lines are optional, since they are the default options. They are added for future reference,
    // i.e. when needed to change those models.
    msystem.GetSettings()->solver.contact_force_model = ChSystemSMC::ContactForceModel::Hertz;
    msystem.GetSettings()->solver.adhesion_force_model = ChSystemSMC::AdhesionForceModel::Constant;

    // Create the fixed and moving bodies
    // ----------------------------------
    AddContainer(&msystem);
    AddCollisionBox(&msystem);

    // Perform the simulation
    // ----------------------

#ifdef CHRONO_OPENGL
    opengl::ChOpenGLWindow& gl_window = opengl::ChOpenGLWindow::getInstance();
    gl_window.Initialize(1280, 720, "boxSMC", &msystem);
    gl_window.SetCamera(ChVector<>(4, -5, 3), ChVector<>(0, 0, 0), ChVector<>(0, 0, 1));
    gl_window.SetRenderMode(opengl::WIREFRAME);

    // Uncomment the following two lines for the OpenGL manager to automatically
    // run the simulation in an infinite loop.
    // gl_window.StartDrawLoop(time_step);
    // return 0;

    while (true) {
        if (gl_window.Active()) {
            gl_window.DoStepDynamics(time_step);
            gl_window.Render();
            if (gl_window.Running()) {
                // Print cumulative contact force on container bin.
                // real3 frc = msystem.GetBodyContactForce(0);
                // std::cout << frc.x << "  " << frc.y << "  " << frc.z << std::endl;
            }
        } else {
            break;
        }
    }
#else
    // Run simulation for specified time
    int num_steps = (int)std::ceil(time_end / time_step);
    double time = 0;

    for (int i = 0; i < num_steps; i++) {
        msystem.DoStepDynamics(time_step);
        time += time_step;
    }
#endif

    return 0;
}
