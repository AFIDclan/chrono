// =============================================================================
// PROJECT CHRONO - http://projectchrono.org
//
// Copyright (c) 2021 projectchrono.org
// All rights reserved.
//
// Use of this source code is governed by a BSD-style license that can be found
// in the LICENSE file at the top level of the distribution and at
// http://projectchrono.org/license-chrono.txt.
//
// =============================================================================
// Authors: Jason Zhou
// =============================================================================
//
// NASA Curiosity Mars Rover Model Class.
// This class contains model for NASA's 6-wheel mars rover curiosity
//
// =============================================================================

#ifndef CURIOSITY_H
#define CURIOSITY_H

#include <array>
#include <fstream>
#include <string>
#include <unordered_map>
#include <vector>

#include "chrono/assets/ChColor.h"
#include "chrono/physics/ChLinkMotorRotationSpeed.h"
#include "chrono/physics/ChLinkDistance.h"
#include "chrono/physics/ChSystem.h"

#include "chrono_models/ChApiModels.h"

#include "chrono/physics/ChInertiaUtils.h"

namespace chrono {

/// Namespace with classes for the Curiosity model.
namespace curiosity {

/// @addtogroup robot_models_curiosity
/// @{

/// Curiosity wheel/suspension identifiers.
enum WheelID {
    LF,  ///< left front
    RF,  ///< right front
    LM,  ///< left middle
    RM,  ///< right middle
    LB,  ///< left back
    RB   ///< right back
};

/// Curiosity chassis type
/// Full Curiosity mass: 840.2 kg
/// Scarecrow barebone mass: 310.2 kg
enum Chassis_Type{
  FullRover,
  Scarecrow
};

/// Base class definition of the Curiosity Rover Part.
class CH_MODELS_API Curiosity_Part {
  public:
    Curiosity_Part(const std::string& name,
               bool fixed,
               std::shared_ptr<ChMaterialSurface> mat,
               ChSystem* system,
               const ChVector<>& body_pos,
               const ChQuaternion<>& body_rot,
               std::shared_ptr<ChBodyAuxRef> chassis_body,
               bool collide);
    virtual ~Curiosity_Part() {}

    /// Return the name of the part.
    const std::string& GetName() const { return m_name; }

    /// Set the name of the part.
    void SetName(const std::string& name) { m_name = name; }

    /// Return the ChBody of the corresponding Curiosity part.
    std::shared_ptr<ChBodyAuxRef> GetBody() const { return m_body; }

    /// Return the ChBody of the chassis wrt the Curiosity part.
    std::shared_ptr<ChBodyAuxRef> GetChassis() const { return m_chassis; }

    /// Return the Position of the Curiosity part.
    const ChVector<>& GetPos() const { return m_body->GetFrame_REF_to_abs().GetPos(); }

    /// Return the Rotation of the Curiosity part.
    const ChQuaternion<>& GetRot() const { return m_body->GetFrame_REF_to_abs().GetRot(); }

  protected:

    /// Initialize the visulization mesh of the Curiosity part.
    void AddVisualizationAssets();

    /// Initialize the collision mesh of the Curiosity part.
    void AddCollisionShapes();

    /// Enable/disable collision.
    void SetCollide(bool state);

    std::string m_name;                        ///< subsystem name
    std::shared_ptr<ChBodyAuxRef> m_body;      ///< rigid body
    std::shared_ptr<ChMaterialSurface> m_mat;  ///< contact material (shared among all shapes)

    std::string m_mesh_name;                  ///< visualization mesh name
    ChVector<> m_offset;                      ///< offset for visualization mesh
    ChColor m_color;                          ///< visualization asset color
    ChSystem* m_system;                       ///< system which Curiosity Part belongs to
    std::shared_ptr<ChBodyAuxRef> m_chassis;  ///< the chassis body for the rover

    ChVector<> m_pos;      ///< Curiosity part's relative position wrt the chassis
    ChQuaternion<> m_rot;  ///< Curiosity part's relative rotation wrt the chassis
    double m_density;      ///< Curiosity part's density

    bool m_collide; ///< Curiosity part's collision indicator
    bool m_fixed; ///< Curiosity part's fixed indication
};

/// Curiosity rover Chassis.
class CH_MODELS_API Curiosity_Chassis : public Curiosity_Part {
  public:
    Curiosity_Chassis(const std::string& name,
                  bool fixed,
                  std::shared_ptr<ChMaterialSurface> mat,
                  ChSystem* system,
                  const ChVector<>& body_pos,
                  const ChQuaternion<>& body_rot,
                  bool collide, 
                  Chassis_Type chassis_type);
    ~Curiosity_Chassis() {}

    /// Initialize the chassis at the specified (absolute) position.
    void Initialize();

    /// Enable/disable collision for the rover chassis.
    void SetCollide(bool state);

  private:

    Chassis_Type m_chassis_type = Chassis_Type::FullRover;

    /// Translate the chassis by the specified value.
    void Translate(const ChVector<>& shift);
    friend class Curiosity_Rover;
};

/// Curiosity rover Wheel.
class CH_MODELS_API Curiosity_Wheel : public Curiosity_Part {
  public:
    Curiosity_Wheel(const std::string& name,
                bool fixed,
                std::shared_ptr<ChMaterialSurface> mat,
                ChSystem* system,
                const ChVector<>& body_pos,
                const ChQuaternion<>& body_rot,
                std::shared_ptr<ChBodyAuxRef> chassis,
                bool collide);
    ~Curiosity_Wheel() {}

    /// Initialize the wheel at the specified (absolute) position.
    void Initialize();

    /// Enable/disable collision for the wheel.
    void SetCollide(bool state);

  private:

    /// Translate the chassis by the specified value.
    void Translate(const ChVector<>& shift);
    friend class Curiosity_Rover;
};


/// Curiosity rover Connecting Arm.
class CH_MODELS_API Curiosity_Arm : public Curiosity_Part {
  public:
    Curiosity_Arm(const std::string& name,
                bool fixed,
                std::shared_ptr<ChMaterialSurface> mat,
                ChSystem* system,
                const ChVector<>& body_pos,
                const ChQuaternion<>& body_rot,
                std::shared_ptr<ChBodyAuxRef> chassis,
                bool collide,
                const int& side);   // 0 == FL, 1 == FR, 2 == BL, 3 == BR
    ~Curiosity_Arm() {}

    /// Initialize the wheel at the specified (absolute) position.
    void Initialize();

    /// Enable/disable collision for the wheel.
    void SetCollide(bool state);

  private:

    /// Translate the chassis by the specified value.
    void Translate(const ChVector<>& shift);
    friend class Curiosity_Rover;
};

/// Curiosity rover steering rod.
class CH_MODELS_API Curiosity_Steer : public Curiosity_Part {
  public:
    Curiosity_Steer(const std::string& name,
                bool fixed,
                std::shared_ptr<ChMaterialSurface> mat,
                ChSystem* system,
                const ChVector<>& body_pos,
                const ChQuaternion<>& body_rot,
                std::shared_ptr<ChBodyAuxRef> chassis,
                bool collide);
    ~Curiosity_Steer() {}

    /// Initialize the wheel at the specified (absolute) position.
    void Initialize();

    /// Enable/disable collision for the wheel.
    void SetCollide(bool state);

  private:

    /// Translate the chassis by the specified value.
    void Translate(const ChVector<>& shift);
    friend class Curiosity_Rover;
};

/// Curiosity rover steering rod.
class CH_MODELS_API Curiosity_Balancer : public Curiosity_Part {
  public:
    Curiosity_Balancer(const std::string& name,
                bool fixed,
                std::shared_ptr<ChMaterialSurface> mat,
                ChSystem* system,
                const ChVector<>& body_pos,
                const ChQuaternion<>& body_rot,
                std::shared_ptr<ChBodyAuxRef> chassis,
                bool collide,
                const int& side);   // 0 - > L, 1 - > R, 2 - > M
    ~Curiosity_Balancer() {}

    /// Initialize the wheel at the specified (absolute) position.
    void Initialize();

    /// Enable/disable collision for the wheel.
    void SetCollide(bool state);

  private:

    /// Translate the chassis by the specified value.
    void Translate(const ChVector<>& shift);
    friend class Curiosity_Rover;
};

/// Curiosity rover class.
/// This class encapsulates the location and rotation information of all Curiosity parts wrt the chassis.
/// This class should be the entry point to create a complete rover.
class CH_MODELS_API CuriosityRover {
  public:
    CuriosityRover(ChSystem* system, 
              const ChVector<>& rover_pos, 
              const ChQuaternion<>& rover_rot, 
              std::shared_ptr<ChMaterialSurface> wheel_mat,
              Chassis_Type chassis_type=Chassis_Type::FullRover);
    CuriosityRover(ChSystem* system, 
              const ChVector<>& rover_pos, 
              const ChQuaternion<>& rover_rot,
              Chassis_Type chassis_type=Chassis_Type::FullRover);
    ~CuriosityRover();

    /// Initialize the Curiosity rover using current parameters.
    void Initialize();

    /// Get the ChSystem
    ChSystem* GetSystem() { return m_system; }

    /// Get Chassis Type
    Chassis_Type GetChassisType();

    /// Set Motor Speed
    void SetMotorSpeed(double rad_speed, WheelID id);

    /// Get wheel speed
    ChVector<> GetWheelSpeed(WheelID id);

    /// Get wheel angular velocity
    ChQuaternion<> GetWheelAngVel(WheelID id);

    /// Get wheel contact force
    ChVector<> GetWheelContactForce(WheelID id);

    /// Get wheel contact torque
    ChVector<> GetWheelContactTorque(WheelID id);

    /// Get wheel total applied force
    ChVector<> GetWheelAppliedForce(WheelID id);

    /// Get wheel total applied torque
    ChVector<> GetWheelAppliedTorque(WheelID id);

    /// Get the chassis body
    std::shared_ptr<ChBodyAuxRef> GetChassisBody();

    /// Get the wheel body
    std::shared_ptr<ChBodyAuxRef> GetWheelBody(WheelID id);

    /// Get total rover mass
    double GetRoverMass();

    /// Get total wheel mass
    double GetWheelMass();

  private:

    /// This function initializes all parameters for the rover
    /// Note: The rover will not be constructed in the ChSystem until Initialize() is called
    void Create();

    ChSystem* m_system;  ///< pointer to the Chrono system

    Chassis_Type m_chassis_type = Chassis_Type::FullRover;  ///< curiosity chassis type

    bool m_custom_wheel_mat;  ///< bool indicating whether the wheel material is customized

    std::shared_ptr<Curiosity_Chassis> m_chassis;               ///< rover chassis
    std::vector<std::shared_ptr<Curiosity_Wheel>> m_wheels;     ///< rover wheels - 1:LF, 2:RF, 3:LM, 4:RM, 5:LB, 6:RB
    std::vector<std::shared_ptr<Curiosity_Arm>> m_arms;         ///< rover arms
    std::vector<std::shared_ptr<Curiosity_Steer>> m_steers;     ///< rover steering rods
    std::vector<std::shared_ptr<Curiosity_Balancer>> m_balancers;   ///< rover balancers parts

    ChQuaternion<> m_rover_rot; ///< rover placement rotation
    ChVector<> m_rover_pos;     ///< rover placement position

    std::vector<std::shared_ptr<ChLinkMotorRotationSpeed>> m_motors;  ///< vector to store motors
                                                                      ///< 0-LF,1-RF,2-LM,3-RM,4-LB,5-RB

    std::vector<std::shared_ptr<ChFunction_Const>> m_motors_func;  ///< constant motor angular speed func

    // model parts material
    std::shared_ptr<ChMaterialSurface> m_chassis_material;  ///< chassis contact material
    std::shared_ptr<ChMaterialSurface> m_wheel_material;    ///< wheel contact material (shared across limbs)
    std::shared_ptr<ChMaterialSurface> m_arm_material;    ///< arm contact material (shared across suspension arms)
    std::shared_ptr<ChMaterialSurface> m_steer_material;    ///< steer contact material
    std::shared_ptr<ChMaterialSurface> m_balancer_material;    ///< balancer contact material

};

/// @} robot_models_curiosity

}  // namespace curiosity
}  // namespace chrono
#endif