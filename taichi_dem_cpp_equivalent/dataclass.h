#pragma once
#include "pch.h"

/*
//======================================
// data class definition
//======================================
*/
// Particle in DEM
// Denver Pilphis : keep spherical shape first, added particle attributes to make the particle kinematically complete
class Grain
{
public:
    // Shape attributes
    Real radius; // Radius, double
    Real mass; // Mass, double
    DEMMatrix inertia; // Moment of inertia tensor in LOCAL coordinate, 3 * 3 matrix with double
    // Material attributes
    // At this stage, the material attributes are written here
    Real density; // Density, double
    Real elasticModulus; // Elastic modulus, double
    Real poissonRatio; // Poisson's ratio, double
    // Translational attributes, all in GLOBAL coordinates
    Vector3 position; // Position, Vector3
    Vector3 velocity; // Velocity, Vector3
    Vector3 acceleration; // Acceleration, Vector3
    Vector3 force; // Force, Vector3

    // Rotational attributes, all in GLOBAL coordinates
    Vector4 quaternion; // Quaternion, Vector4, order in[w, x, y, z]
    Vector3 omega; // Angular velocity, Vector3
    Vector3 omega_dot; // Angular acceleration, Vector3
    Vector3 moment; // Total moment(including torque), Vector3
};

// Wall in DEM
// Only premitive wall is implemented
class Wall
{
public:
    // Wall equation : Ax + By + Cz - D = 0
    // Reference : Peng and Hanley(2019) Contact detection between convex polyhedra and superquadrics in discrete element codes.
    // https ://doi.org/10.1016/j.powtec.2019.07.082
    // Eq. (8)
    Vector3 normal; // Outer normal vector of the wall, [A, B, C]
    Real distance; // Distance between origin and the wall, D
    // Material properties
    Real density; // Density of the wall
    Real elasticModulus; // Elastic modulus of the wall
    Real poissonRatio; // Poisson's ratio of the wall
};

// Contact in DEM
// In this example, the Edinburgh Bond Particle Model(EBPM), along with Hertz - Mindlin model, is implemented
// Reference: Brown et al. (2014) A bond model for DEM simulation of cementitious materialsand deformable structures.
// https ://doi.org/10.1007/s10035-014-0494-4
// Reference : Mindlin and Deresiewicz(1953) Elastic spheres in contact under varying oblique forces.
// https ://doi.org/10.1115/1.4010702
class Contact
{
public:
    Contact();
    // Contact status
    // Integer isActive; // Deprecated // Contact exists : 1 - exist 0 - not exist
    Integer isBonded; // Contact is bonded : 1 - bonded, use EBPM 0 - unbonded, use Hertz - Mindlin
    // Common Parameters
    DEMMatrix rotationMatrix; // Rotation matrix from global to local system of the contact
    Vector3 position; // Position of contact point in GLOBAL coordinate
    // EBPM parameters
    Real radius_ratio; // Section radius ratio
    // radius : Real // Section radius : r = rratio * min(r1, r2), temporarily calculated in evaluation
    Real length; // Length of the bond
    Real elasticModulus; // Elastic modulus of the bond
    Real poissonRatio; // Possion's ratio of the bond
    Real compressiveStrength; // Compressive strength of the bond
    Real tensileStrength; // Tensile strength of the bond
    Real shearStrength; // Shear strength of the bond
    Vector3 force_a; // Contact force at side a in LOCAL coordinate
    Vector3 moment_a; // Contact moment / torque at side a in LOCAL coordinate
    Vector3 force_b; // Contact force at side b in LOCAL coordinate
    Vector3 moment_b; // Contact moment / torque at side b in LOCAL coordinate
    // Hertz - Mindlin parameters
    // normalStiffness: Real // Normal stiffness, middleware parameter only
    // shearStiffness: Real // Shear stiffness, middleware parameter only
    Real coefficientFriction; // Friction coefficient, double
    Real coefficientRestitution; // Coefficient of resitution, double
    Real coefficientRollingResistance; // Coefficient of rolling resistance, double
    Vector3 shear_displacement; // Shear displacement stored in the contact
};

class DEMSolverConfig
{
public:
    DEMSolverConfig();

    Vector3 gravity;
    Real dt;
    Real target_time;
    Integer nsteps;
    Real saving_interval_time;
    Integer saving_interval_steps;
};

class DEMSolver
{
public:
    DEMSolver(const DEMSolverConfig& input_config);
    void save(const std::string& file_name, const Real& time);
    void init_particle_fields(const std::string& file_name = "input.p4p");
    void init_grid_fields(const Integer& grid_n);
    void clear_state();
    void apply_body_force();
    void update();
    void evaluate(const Integer& i, const Integer& j);
    void resolve(const Integer& i, const Integer& j);
    void evaluate_wall(const Integer& i, const Integer& j);
    void resolve_wall();
    void bond_detect(const Integer& i, const Integer& j);
    void bond();
    void contact();
    void run_simulation();
    void init_simulation();

    DEMSolverConfig config;
    // element fields
    StructFieldObj<Grain> gf;
    StructFieldObj2<Contact> cf;
    StructFieldObj<Wall> wf;
    StructFieldObj2<Contact> wcf;
    StructField<Integer> particle_id;
    // grid fields
    StructField2<Integer> list_head;
    StructField2<Integer> list_cur;
    StructField2<Integer> list_tail;
    StructField2<Integer> grain_count;
    StructField<Integer> column_sum;
    StructField2<Integer> prefix_sum;
};

