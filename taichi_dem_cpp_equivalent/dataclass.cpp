#include "dataclass.h"

DEMSolverConfig::DEMSolverConfig()
{
    // Gravity, a global parameter
    // Denver Pilphis : in this example, we assign no gravity
    gravity = Vector3(gravity_x, gravity_y, gravity_z);
    // Time step, a global parameter
    dt = time_increment; // Larger dt might lead to unstable results.
    target_time = total_time;
    // No. of steps for run, a global parameter
    nsteps = (int)(target_time / dt);
    saving_interval_time = save_interval_time;
    saving_interval_steps = (int)(saving_interval_time / dt);
}

DEMSolver::DEMSolver(const DEMSolverConfig& input_config)
{
    config = DEMSolverConfig(input_config);
}

void DEMSolver::save(const std::string& file_name, const Real& etime)
{
    // P4P file for particles
    std::ofstream fp(file_name + ".p4p", std::ios::out|std::ios::app);
    Integer n = gf.rows();
    fp << "TIMESTEP  PARTICLES" << std::endl;
    fp << etime << " " << n << std::endl;

    fp << "ID  GROUP  VOLUME  MASS  PX  PY  PZ  VX  VY  VZ" << std::endl;
    for (int i = 0; i < n; ++i)
    {
        // GROUP omitted
        Integer group = 0;
        Real volume = gf[i]->mass / gf[i]->density;
        Real mass = gf[i]->mass;
        Vector3 p = gf[i]->position;
        Vector3 v = gf[i]->velocity;
        Real px = p[0];
        Real py = p[1];
        Real pz = p[2];
        Real vx = v[0];
        Real vy = v[1];
        Real vz = v[2];
        fp << (i + 1) << " " << group << " " << volume << " " << mass << " "
           << px << " " << py << " " << pz << " "
           << vx << " " << vy << " " << vz << std::endl;
    }
    fp << std::endl;
    fp.close();

    // P4C file for contacts
    std::vector<std::string> ccache = { "P1  P2  CX  CY  CZ  FX  FY  FZ  CONTACT_IS_BONDED\n" };
    Integer ncontact = 0;
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
        {
            if (i >= j) continue;
            if (!cf(i, j)) continue; // if (!cf(i, j)->isActive) continue;
            // GROUP omitted
            Integer p1 = i + 1;
            Integer p2 = j + 1;
            Real cx = cf(i, j)->position[0];
            Real cy = cf(i, j)->position[1];
            Real cz = cf(i, j)->position[2];
            Real fx = cf(i, j)->force_a[0]; // Denver Pilphis : the force is not stored in present version
            Real fy = cf(i, j)->force_a[1]; // Denver Pilphis : the force is not stored in present version
            Real fz = cf(i, j)->force_a[2]; // Denver Pilphis : the force is not stored in present version
            Integer bonded = cf(i, j)->isBonded;
            ccache.emplace_back(
                std::to_string(p1) + " " +
                std::to_string(p2) + " " +
                std::to_string(cx) + " " +
                std::to_string(cy) + " " +
                std::to_string(cz) + " " +
                std::to_string(fx) + " " +
                std::to_string(fy) + " " +
                std::to_string(fz) + " " +
                std::to_string(bonded) + "\n"
            );
            ++ncontact;
        }

    fp.open(file_name + ".p4c", std::ios::out|std::ios::app);
    fp << "TIMESTEP  CONTACTS" << std::endl;
    fp << etime << " " << ncontact << std::endl;
    for (const std::string& str : ccache)
        fp << str;
    fp << std::endl;
    fp.close();
}

void DEMSolver::init_particle_fields(const std::string& file_name)
{
    std::ifstream fp(file_name, std::ios::in);
    std::string line;
    std::getline(fp, line);
    Integer n; Real etime;
    fp >> etime >> n;
    std::getline(fp, line); // Finish off this line
    Integer nwall = 1; // Denver Pilphis: hard coding

    // Initialize particles
    gf.resize(n, 1); gf.fill(nullptr);
    cf.resize(n, n); cf.fill(nullptr);
    wf.resize(nwall, 1); wf.fill(nullptr);
    wcf.resize(n, nwall); wcf.fill(nullptr);
    particle_id.resize(n, 1);

    std::getline(fp, line); // "ID  GROUP  VOLUME  MASS  PX  PY  PZ  VX  VY  VZ" line

    // Processing particles
    Integer i, group; Real volume, mass, px, py, pz, vx, vy, vz;
    while (fp >> i >> group >> volume >> mass >> px >> py >> pz >> vx >> vy >> vz)
    {
        --i; // Particle id starts from 1 in P4P
        const Real density = mass / volume;
        const Real radius = cbrt(volume * 3.0 / 4.0 / M_PI);
        const Real inertia = 2.0 / 5.0 * mass * radius * radius;
        gf[i] = new Grain();
        gf[i]->density = density;
        gf[i]->mass = mass;
        gf[i]->radius = radius;
        gf[i]->elasticModulus = 7e10; // Denver Pilphis : hard coding need to be modified in the future
        gf[i]->poissonRatio = 0.25; // Denver Pilphis : hard coding need to be modified in the future
        gf[i]->position = Vector3(px, py, pz);
        gf[i]->velocity = Vector3(vx, vy, vz);
        gf[i]->acceleration = Vector3(0.0, 0.0, 0.0);
        gf[i]->force = Vector3(0.0, 0.0, 0.0);
        gf[i]->quaternion = Vector4(1.0, 0.0, 0.0, 0.0);
        gf[i]->omega = Vector3(0.0, 0.0, 0.0);
        gf[i]->omega_dot = Vector3(0.0, 0.0, 0.0);
        gf[i]->moment = Vector3(0.0, 0.0, 0.0);
        gf[i]->inertia = inertia * DEMMatrix::Identity();
    }

    fp.close();
    // Input wall
    // Denver Pilphis : hard coding : need to be modified in the future
    for (int j = 0; j < wf.rows(); ++j)
    {
        wf[j] = new Wall();
        wf[j]->normal = Vector3(wall_normal_x, wall_normal_y, wall_notmal_z); // Outer normal vector of the wall, [A, B, C]
        wf[j]->distance = wall_distance; // Distance between origin and the wall, D
            // Material properties
        wf[j]->density = 7800.0; // Density of the wall
        wf[j]->elasticModulus = 2e11; // Elastic modulus of the wall
        wf[j]->poissonRatio = 0.25; // Poisson's ratio of the wall
    }
}

void DEMSolver::init_grid_fields(const Integer& grid_n)
{
    // GPU grid parameters
    list_head.resize(grid_n, grid_n);
    list_cur.resize(grid_n, grid_n);
    list_tail.resize(grid_n, grid_n);
    grain_count.resize(grid_n, grid_n);
    column_sum.resize(grid_n, 1);
    prefix_sum.resize(grid_n, grid_n);
}

void DEMSolver::clear_state()
{
    for (Grain * const p : gf)
    {
        p->force = Vector3(0.0, 0.0, 0.0);
        p->moment = Vector3(0.0, 0.0, 0.0);
    }
}

void DEMSolver::apply_body_force()
{
    Vector3 g = config.gravity;
    for (Grain* const p : gf)
    {
        // Gravity
        p->force += p->mass * g;

        // GLOBAL damping
        // Add GLOBAL damping for EBPM, GLOBAL damping is assigned to particles
        Real t_d = 0.0; // Denver Pilphis : hard coding - should be modified in the future
        Vector3 damp_force = Vector3(0.0, 0.0, 0.0);
        Vector3 damp_moment = Vector3(0.0, 0.0, 0.0);
        for (int j = 0; j < 3; ++j)
        { 
            damp_force[j] = -t_d * abs(p->force[j]) * sgn(p->velocity[j]);
            damp_moment[j] = -t_d * abs(p->moment[j]) * sgn(p->omega[j]);
        }
        p->force += damp_force;
        p->moment += damp_moment;
    }

}

// NVE integrator
void DEMSolver::update()
{
    // For debug only
    Real total_kinematic_energy = 0.0;
    Real dt = config.dt;
    for (Grain * const p : gf)
    {
        // Translational
        // Velocity Verlet integrator is adopted
        // Reference: https://www.algorithm-archive.org/contents/verlet_integration/verlet_integration.html
        p->acceleration = p->force / p->mass;
        p->position += p->velocity * dt + 0.5 * p->acceleration * pow(dt, 2);
        p->velocity += p->acceleration * dt;
        // Rotational
        // Angular acceleration should be calculated via Euler's equation for rigid body
        // Reference: Peng et al. (2021) Critical time step for discrete element method simulations of convex particles with central symmetry.
        // https ://doi.org/10.1002/nme.6568
        // Eqs. (5) - (16)
        DEMMatrix rotational_matrix = quat2RotMatrix(p->quaternion);
        Vector3 moment_local = rotational_matrix * p->moment;
        Vector3 omega_local = rotational_matrix * p->omega;
        Vector3 omega_dot_local = p->inertia.inverse() * (moment_local - omega_local.cross(p->inertia * omega_local));
        p->omega_dot = rotational_matrix.inverse() * omega_dot_local;
        // Update particle orientation
        // Reference: Lu et al. (2015) Discrete element models for non - spherical particle systems : From theoretical developments to applications.
        // http ://dx.doi.org/10.1016/j.ces.2014.11.050
        // Eq. (6)
        // Originally from Langston et al. (2004) Distinct element modelling of non - spherical frictionless particle flow.
        // https://doi.org/10.1016/j.ces.2003.10.008
        Real dq0 = -0.5 * (p->quaternion[1] * p->omega[0] + p->quaternion[2] * p->omega[1] + p->quaternion[3] * p->omega[2]);
        Real dq1 = +0.5 * (p->quaternion[0] * p->omega[0] - p->quaternion[3] * p->omega[1] + p->quaternion[2] * p->omega[2]);
        Real dq2 = +0.5 * (p->quaternion[3] * p->omega[0] + p->quaternion[0] * p->omega[1] + p->quaternion[1] * p->omega[2]);
        Real dq3 = +0.5 * (-p->quaternion[2] * p->omega[0] + p->quaternion[1] * p->omega[1] + p->quaternion[0] * p->omega[2]);
        p->quaternion[0] += dq0;
        p->quaternion[1] += dq1;
        p->quaternion[2] += dq2;
        p->quaternion[3] += dq3;
        p->quaternion.normalize();
        // Update angular velocity
        p->omega += p->omega_dot * dt;

        // For debug only
        total_kinematic_energy += p->mass / 2.0 * p->velocity.dot(p->velocity);
    }
    std::cout << total_kinematic_energy << std::endl;
}

// Contact resolution and evaluation
void DEMSolver::evaluate(const Integer& i, const Integer& j)
{
    Real dt = config.dt;
    // Contact resolution
    // Find out rotation matrix
    // https://math.stackexchange.com/questions/180418/calculate-rotation-matrix-to-align-vector-a-to-vector-b-in-3d
    Vector3 a = (gf[j]->position - gf[i]->position).normalized();
    Vector3 b = Vector3(1.0, 0.0, 0.0); // Local x coordinate
    Vector3 v = a.cross(b);
    Real s = v.norm();
    Real c = a.dot(b);
    if (s < DoublePrecisionTolerance)
    {
        if (c > 0.0)
            cf(i, j)->rotationMatrix = DEMMatrix::Identity();
        else
            cf(i, j)->rotationMatrix << -1.0, 0.0, 0.0,
                                        0.0, 1.0, 0.0,
                                        0.0, 0.0, -1.0;
    }
    else
    {
        DEMMatrix vx;
        vx << 0.0, -v[2], v[1],
              v[2], 0.0, -v[0],
              -v[1], v[0], 0.0;
        cf(i, j)->rotationMatrix = DEMMatrix::Identity() + vx + ((1.0 - c) / pow(s, 2)) * vx * vx;
    }
    cf(i, j)->length = (gf[j]->position - gf[i]->position).norm();

    // Contact evaluation(with contact model)
    if (cf(i, j)->isBonded) // Bonded, use EBPM
    {
        cf(i, j)->position = 0.5 * (gf[i]->position + gf[j]->position);
        Vector3 disp_a = cf(i, j)->rotationMatrix * gf[i]->velocity * dt;
        Vector3 disp_b = cf(i, j)->rotationMatrix * gf[j]->velocity * dt;
        Vector3 rot_a = cf(i, j)->rotationMatrix * gf[i]->omega * dt;
        Vector3 rot_b = cf(i, j)->rotationMatrix * gf[j]->omega * dt;
        EBPMForceDisplacementVector dispVector;
        for (int k = 0; k < 3; ++k)
        {
            // Offsets: 0, 3, 6, 9
            dispVector[0 + k] = disp_a[k];
            dispVector[3 + k] = rot_a[k];
            dispVector[6 + k] = disp_b[k];
            dispVector[9 + k] = rot_b[k];
        }
        Real r_b = cf(i, j)->radius_ratio * std::min(gf[i]->radius, gf[j]->radius);
        Real L_b = cf(i, j)->length;
        Real E_b = cf(i, j)->elasticModulus;
        Real nu = cf(i, j)->poissonRatio;
        Real I_b = pow(r_b, 4) * M_PI / 4.0;
        Real phi = 20.0 / 3.0 * pow(r_b, 2) / pow(L_b, 2) * (1.0 + nu);
        Real A_b = M_PI * pow(r_b, 2);
        Real k1 = E_b * A_b / L_b;
        Real k2 = 12.0 * E_b * I_b / pow(L_b, 3) / (1.0 + phi);
        Real k3 = 6.0 * E_b * I_b / pow(L_b, 2) / (1.0 + phi);
        Real k4 = E_b * I_b / L_b / (1.0 + nu);
        Real k5 = E_b * I_b * (4.0 + phi) / L_b / (1.0 + phi);
        Real k6 = E_b * I_b * (2.0 - phi) / L_b / (1.0 + phi);
        EBPMStiffnessMatrix K;
        K << k1, 0, 0, 0, 0, 0, -k1, 0, 0, 0, 0, 0,
            0, k2, 0, 0, 0, k3, 0, -k2, 0, 0, 0, k3,
            0, 0, k2, 0, -k3, 0, 0, 0, -k2, 0, -k3, 0,
            0, 0, 0, k4, 0, 0, 0, 0, 0, -k4, 0, 0,
            0, 0, -k3, 0, k5, 0, 0, 0, k3, 0, k6, 0,
            0, k3, 0, 0, 0, k5, 0, -k3, 0, 0, 0, k6,
            -k1, 0, 0, 0, 0, 0, k1, 0, 0, 0, 0, 0,
            // K(7, 5) is WRONG in original EBPM document
            // ¦¤Fay + ¦¤Fby is nonzero
            // which does not satisfy the equilibrium
            // Acknowledgement to Dr. Xizhong Chen in 
            // Department of Chemical and Biological Engineering,
            // The University of Sheffield
            // Reference: Chen et al. (2022) A comparative assessment and unification of bond models in DEM simulations.
            // https://doi.org/10.1007/s10035-021-01187-2
            0, -k2, 0, 0, 0, -k3, 0, k2, 0, 0, 0, -k3,
            0, 0, -k2, 0, k3, 0, 0, 0, k2, 0, k3, 0,
            0, 0, 0, -k4, 0, 0, 0, 0, 0, k4, 0, 0,
            0, 0, -k3, 0, k6, 0, 0, 0, k3, 0, k5, 0,
            0, k3, 0, 0, 0, k6, 0, -k3, 0, 0, 0, k5;
        EBPMForceDisplacementVector forceVector = K * dispVector;
        cf(i, j)->force_a += Vector3(forceVector[0], forceVector[1], forceVector[2]);
        cf(i, j)->moment_a += Vector3(forceVector[3], forceVector[4], forceVector[5]);
        cf(i, j)->force_b += Vector3(forceVector[6], forceVector[7], forceVector[8]);
        cf(i, j)->moment_b += Vector3(forceVector[9], forceVector[10], forceVector[11]);

        if ((cf(i, j)->force_a + cf(i, j)->force_b).norm() > DoublePrecisionTolerance)
            throw "Force equillibrium error!";
        if (abs(cf(i, j)->moment_a[0] + cf(i, j)->moment_b[0]) > DoublePrecisionTolerance)
            throw "Torque equillibrium error!";

        // Check whether the bond fails
        Real sigma_c_a = cf(i, j)->force_b[0] / A_b - r_b / I_b * sqrt(pow(cf(i, j)->moment_a[1], 2) + pow(cf(i, j)->moment_a[2], 2));
        Real sigma_c_b = cf(i, j)->force_b[0] / A_b - r_b / I_b * sqrt(pow(cf(i, j)->moment_b[1], 2) + pow(cf(i, j)->moment_b[2], 2));
        Real sigma_c_max = -std::min(sigma_c_a, sigma_c_b);
        Real sigma_t_a = sigma_c_a;
        Real sigma_t_b = sigma_c_b;
        Real sigma_t_max = std::max(sigma_t_a, sigma_t_b);
        Real tau_max = abs(cf(i, j)->moment_a[0]) * r_b / 2.0 / I_b + 4.0 / 3.0 / A_b * sqrt(pow(cf(i, j)->force_a[1], 2) + pow(cf(i, j)->force_a[2], 2));
        if (sigma_c_max >= cf(i, j)->compressiveStrength) // Compressive failure
        {
            // cf(i, j)->isBonded = 0
            // cf(i, j)->isActive = 0
            delete cf(i, j);
            cf(i, j) = nullptr;
            // For debug
            std::cout << "Bond compressive failure at: " << i << ", " << j << std::endl;
        }
        else if (sigma_t_max >= cf(i, j)->tensileStrength) // Tensile failure
        {
            // cf(i, j)->isBonded = 0
            // cf(i, j)->isActive = 0
            delete cf(i, j);
            cf(i, j) = nullptr;
            // For debug
            std::cout << "Bond tensile failure at: " << i << ", " << j << std::endl;
        }
        else if (tau_max >= cf(i, j)->shearStrength) // Shear failure
        {
            // cf(i, j)->isBonded = 0
            // cf(i, j)->isActive = 0
            delete cf(i, j);
            cf(i, j) = nullptr;
            // For debug
            std::cout << "Bond shear failure at: " << i << ", " << j << std::endl;
        }
        else // Intact bond, need to conduct force to particles
        {
            // Notice the inverse of signs due to Newton's third law
            // and LOCAL to GLOBAL coordinates
            gf[i]->force += cf(i, j)->rotationMatrix.inverse() * (-cf(i, j)->force_a);
            gf[j]->force += cf(i, j)->rotationMatrix.inverse() * (-cf(i, j)->force_b);
            gf[i]->moment += cf(i, j)->rotationMatrix.inverse() * (-cf(i, j)->moment_a);
            gf[j]->moment += cf(i, j)->rotationMatrix.inverse() * (-cf(i, j)->moment_b);
        }
    }
    else // Non - bonded, use Hertz - Mindlin
    {
        // Calculation relative translational and rotational displacements
        // Need to include the impact of particle rotation in contact relative translational displacement
        // Reference: Peng et al. (2021) Critical time step for discrete element method simulations of convex particles with central symmetry.
        // https ://doi.org/10.1002/nme.6568
        // Eqs. (1) - (2)
        // Implementation reference : https://github.com/CFDEMproject/LIGGGHTS-PUBLIC/blob/master/src/surface_model_default.h
        // Lines 140 - 189
        Real gap = cf(i, j)->length - gf[i]->radius - gf[j]->radius; // gap must be negative to ensure an intact contact
        Real delta_n = abs(gap); // For parameter calculation only

        // For debug only
        Real delta_n_ratio = delta_n / std::min(gf[i]->radius, gf[j]->radius);
        if (delta_n_ratio > 0.05)
            std::cout << "Overlap limit exceeded: " << delta_n_ratio << std::endl;

        cf(i, j)->position = gf[i]->position + (gf[j]->position - gf[i]->position).normalized() * (gf[i]->radius - delta_n / 2.0);
        Vector3 r_i = cf(i, j)->position - gf[i]->position;
        Vector3 r_j = cf(i, j)->position - gf[j]->position;
        // Velocity of a point on the surface of a rigid body
        
        /*
        Vector3 temp0 = gf[i]->velocity;
        Vector3 temp1 = gf[j]->velocity;
        Vector3 temp2 = gf[i]->omega;
        Vector3 temp3 = gf[j]->omega;
        */

        Vector3 v_c_i = gf[i]->omega.cross(r_i) + gf[i]->velocity;
        Vector3 v_c_j = gf[j]->omega.cross(r_j) + gf[j]->velocity;
        Vector3 v_c = cf(i, j)->rotationMatrix * (v_c_j - v_c_i); // LOCAL coordinate
        // Parameter calculation
        // Reference: https://www.cfdem.com/media/DEM/docu/gran_model_hertz.html
        Real Y_star = 1.0 / ((1.0 - pow(gf[i]->poissonRatio, 2)) / gf[i]->elasticModulus + (1.0 - pow(gf[j]->poissonRatio, 2)) / gf[j]->elasticModulus);
        Real G_star = 1.0 / (2.0 * (2.0 - gf[i]->poissonRatio) * (1.0 + gf[i]->poissonRatio) / gf[i]->elasticModulus + 2.0 * (2.0 - gf[j]->poissonRatio) * (1.0 + gf[j]->poissonRatio) / gf[j]->elasticModulus);
        Real R_star = 1.0 / (1.0 / gf[i]->radius + 1.0 / gf[j]->radius);
        Real m_star = 1.0 / (1.0 / gf[i]->mass + 1.0 / gf[j]->mass);
        Real beta = log(cf(i, j)->coefficientRestitution) / sqrt(pow(log(cf(i, j)->coefficientRestitution), 2) + pow(M_PI, 2));
        Real S_n = 2.0 * Y_star * sqrt(R_star * delta_n);
        Real S_t = 8.0 * G_star * sqrt(R_star * delta_n);
        Real k_n = 4.0 / 3.0 * Y_star * sqrt(R_star * delta_n);
        Real gamma_n = -2.0 * beta * sqrt(5.0 / 6.0 * S_n * m_star); // Check whether gamma_n >= 0
        Real k_t = 8.0 * G_star * sqrt(R_star * delta_n);
        Real gamma_t = -2.0 * beta * sqrt(5.0 / 6.0 * S_t * m_star); // Check whether gamma_t >= 0

        // Shear displacement increments
        Vector3 shear_increment = v_c * dt;
        shear_increment[0] = 0.0; // Remove the normal direction
        cf(i, j)->shear_displacement += shear_increment;
        // Normal direction - LOCAL - the force towards particle j
        Vector3 F = Vector3(0.0, 0.0, 0.0);
        // Reference: Peng et al. (2021) Critical time step for discrete element method simulations of convex particles with central symmetry.
        // https://doi.org/10.1002/nme.6568
        // Eq. (29)
        // Be aware of signs
        F[0] = -k_n * gap - gamma_n * v_c[0];
        // if (F[0] < 0.0)
        //    std::cout << "Normal force problem: " << F[0] << std::endl;
        // Shear direction - LOCAL - the force towards particle j
        Vector3 try_shear_force = -k_t * cf(i, j)->shear_displacement;
        if (try_shear_force.norm() >= cf(i, j)->coefficientFriction * F[0]) // Sliding
        {
            Real ratio = cf(i, j)->coefficientFriction * F[0] / try_shear_force.norm();
            F[1] = try_shear_force[1] * ratio;
            F[2] = try_shear_force[2] * ratio;
            cf(i, j)->shear_displacement[1] = F[1] / k_t;
            cf(i, j)->shear_displacement[2] = F[2] / k_t;
        }
        else // No sliding
        {
            F[1] = try_shear_force[1] - gamma_t * v_c[1];
            F[2] = try_shear_force[2] - gamma_t * v_c[2];
        }
        
        // No moment is conducted in Hertz - Mindlin model

        // Assigning contact force to particles
        // Notice the inverse of signs due to Newton's third law
        // and LOCAL to GLOBAL coordinates
        Vector3 F_i_global = cf(i, j)->rotationMatrix.inverse() * (-F);
        Vector3 F_j_global = cf(i, j)->rotationMatrix.inverse() * F;
        gf[i]->force += F_i_global;
        gf[j]->force += F_j_global;
        // As the force is at contact position
        // additional moments will be assigned to particles
        // Reference: Peng et al. (2021) Critical time step for discrete element method simulations of convex particles with central symmetry.
        // https ://doi.org/10.1002/nme.6568
        // Eqs. (3) - (4)
        gf[i]->moment += r_i.cross(F_i_global);
        gf[j]->moment += r_j.cross(F_j_global);
    }
}

// Particle-particle contact detection
void DEMSolver::resolve(const Integer& i, const Integer& j)
{
    // Particle - particle contacts
    if (cf(i, j)) // if (cf(i, j)->isActive) // Existing contact
        if (cf(i, j)->isBonded) // Bonded contact
            evaluate(i, j); // Bonded contact must exist.Go to evaluationand if bond fails, the contact state will change thereby.
        else // Non - bonded contact, should check whether two particles are still in contact
        {
            if (-gf[i]->radius - gf[j]->radius + (gf[j]->position - gf[i]->position).norm() < 0) // Use PFC's gap < 0 criterion
                evaluate(i, j);
            else
            {
                // cf(i, j).isActive = 0
                delete cf(i, j);
                cf(i, j) = nullptr;
            }
                
        }
    else
    {
        if (-gf[i]->radius - gf[j]->radius + (gf[j]->position - gf[i]->position).norm() < 0) // Use PFC's gap < 0 criterion
        {
            cf(i, j) = new Contact(); // Hertz - Mindlin model
            // cf(i, j)->isActive = 1;
            cf(i, j)->isBonded = 0;
            cf(i, j)->coefficientFriction = 0.3; // Denver Pilphis : hard coding need to be modified in the future
            cf(i, j)->coefficientRestitution = 0.9; // Denver Pilphis : hard coding need to be modified in the future
            cf(i, j)->coefficientRollingResistance = 0.01; // Denver Pilphis : hard coding need to be modified in the future
            cf(i, j)->shear_displacement = Vector3(0.0, 0.0, 0.0);
            
            evaluate(i, j); // Send to evaluation using Hertz - Mindlin contact model
        }
    }
}

// Particle - wall contact evaluation
// Contact model is Hertz - Mindlin
void DEMSolver::evaluate_wall(const Integer& i, const Integer& j)
{
    Real dt = config.dt;
    // Contact resolution
    // Find out rotation matrix
    // https://math.stackexchange.com/questions/180418/calculate-rotation-matrix-to-align-vector-a-to-vector-b-in-3d
    Vector3 a = wf[j]->normal;
    Vector3 b = Vector3(1.0, 0.0, 0.0); // Local x coordinate
    Vector3 v = a.cross(b);
    Real s = v.norm();
    Real c = a.dot(b);
    if (s < DoublePrecisionTolerance)
    {
        if (c > 0.0)
            wcf(i, j)->rotationMatrix = DEMMatrix::Identity();
        else
            wcf(i, j)->rotationMatrix << -1.0, 0.0, 0.0,
                                         0.0, 1.0, 0.0,
                                         0.0, 0.0, -1.0;
    }
    else
    {
        DEMMatrix vx;
        vx << 0.0, -v[2], v[1],
            v[2], 0.0, -v[0],
            -v[1], v[0], 0.0;
        wcf(i, j)->rotationMatrix = DEMMatrix::Identity() + vx + ((1.0 - c) / pow(s, 2)) * vx * vx;
    }

    // Calculate relative translational and rotational displacements
    Real distance = gf[i]->position.dot(wf[j]->normal) - wf[j]->distance; // Distance < 0 means that particle is beneath the plane
    Real gap = abs(distance) - gf[i]->radius; // gap must be negative
    Real delta_n = abs(gap); // For parameter calculation only

    // For debug only
    Real delta_n_ratio = delta_n / gf[i]->radius;
    if (delta_n_ratio > 0.05)
        std::cout << "Overlap limit exceeded: " << delta_n_ratio << std::endl;

    Vector3 r_i = -distance * wf[j]->normal / abs(distance) * (abs(distance) + delta_n / 2.0);
    wcf(i, j)->position = gf[i]->position + r_i;
    // Velocity of a point on the surface of a rigid body
    Vector3 v_c_i = gf[i]->omega.cross(r_i) + gf[i]->velocity;
    Vector3 v_c = wcf(i, j)->rotationMatrix * (-v_c_i); // LOCAL coordinate
    // Parameter calculation
    // Reference: https://www.cfdem.com/media/DEM/docu/gran_model_hertz.html
    Real Y_star = 1.0 / ((1.0 - pow(gf[i]->poissonRatio, 2)) / gf[i]->elasticModulus + (1.0 - pow(wf[j]->poissonRatio, 2)) / wf[j]->elasticModulus);
    Real G_star = 1.0 / (2.0 * (2.0 - gf[i]->poissonRatio) * (1.0 + gf[i]->poissonRatio) / gf[i]->elasticModulus + 2.0 * (2.0 - wf[j]->poissonRatio) * (1.0 + wf[j]->poissonRatio) / wf[j]->elasticModulus);
    Real R_star = gf[i]->radius;
    Real m_star = gf[i]->mass;
    Real beta = log(wcf(i, j)->coefficientRestitution) / sqrt(pow(log(wcf(i, j)->coefficientRestitution), 2) + pow(M_PI, 2));
    Real S_n = 2.0 * Y_star * sqrt(R_star * delta_n);
    Real S_t = 8.0 * G_star * sqrt(R_star * delta_n);
    Real k_n = 4.0 / 3.0 * Y_star * sqrt(R_star * delta_n);
    Real gamma_n = -2.0 * beta * sqrt(5.0 / 6.0 * S_n * m_star); // Check whether gamma_n >= 0
    Real k_t = 8.0 * G_star * sqrt(R_star * delta_n);
    Real gamma_t = -2.0 * beta * sqrt(5.0 / 6.0 * S_t * m_star); // Check whether gamma_t >= 0

    // Shear displacement increments
    Vector3 shear_increment = v_c * dt;
    shear_increment[0] = 0.0; // Remove the normal direction
    wcf(i, j)->shear_displacement += shear_increment;
    // Normal direction - LOCAL - the force towards the wall
    Vector3 F = Vector3(0.0, 0.0, 0.0);
    // Reference: Peng et al. (2021) Critical time step for discrete element method simulations of convex particles with central symmetry.
    // https ://doi.org/10.1002/nme.6568
    // Eq. (29)
    // Be aware of signs
    F[0] = -k_n * gap - gamma_n * v_c[0];
    // Shear direction - LOCAL - the force towards the wall
    Vector3 try_shear_force = -k_t * wcf(i, j)->shear_displacement;
    if (try_shear_force.norm() >= wcf(i, j)->coefficientFriction * F[0]) // Sliding
    {
        Real ratio = wcf(i, j)->coefficientFriction * F[0] / try_shear_force.norm();
        F[1] = try_shear_force[1] * ratio;
        F[2] = try_shear_force[2] * ratio;
        wcf(i, j)->shear_displacement[1] = F[1] / k_t;
        wcf(i, j)->shear_displacement[2] = F[2] / k_t;
    }
    else // No sliding
    {
        F[1] = try_shear_force[1] - gamma_t * v_c[1];
        F[2] = try_shear_force[2] - gamma_t * v_c[2];
    }
       
    // No moment is conducted in Hertz - Mindlin model

    // Assigning contact force to particles
    // Notice the inverse of signs due to Newton's third law
    // and LOCAL to GLOBAL coordinates
    // As the force is at contact position
    // additional moments will be assigned to particles
    Vector3 F_i_global = wcf(i, j)->rotationMatrix.inverse() * (-F);
    gf[i]->force += F_i_global;
    // Reference: Peng et al. (2021) Critical time step for discrete element method simulations of convex particles with central symmetry.
    // https ://doi.org/10.1002/nme.6568
    // Eqs. (3) - (4)
    gf[i]->moment += r_i.cross(F_i_global);
}

// Particle-wall contact detection
void DEMSolver::resolve_wall()
{
    // Denver Pilphis : particle - wall neighboring search has not been implemented
    // and thus all particles will be contact detection with the wall
    for (int i = 0; i < gf.rows(); ++i)
        for (int j = 0; j < wf.rows(); ++j)
        {
            // Particle - wall contacts
            if (wcf(i, j)) // if (wcf(i, j).isActive) : // Existing contact
            {
                if (abs(gf[i]->position.dot(wf[j]->normal) - wf[j]->distance) >= gf[i]->radius) // Non - contact
                {
                    // wcf(i, j).isActive = 0
                    delete wcf(i, j);
                    wcf(i, j) = nullptr;
                }
                else // Contact
                    evaluate_wall(i, j);
            }
            else
            {
                if (abs(gf[i]->position.dot(wf[j]->normal) - wf[j]->distance) < gf[i]->radius) // Contact
                {
                    wcf(i, j) = new Contact(); // Hertz - Mindlin model
                    // wcf(i, j)->isActive = 1;
                    wcf(i, j)->isBonded = 0;
                    wcf(i, j)->coefficientFriction = 0.35; // Denver Pilphis : hard coding need to be modified in the future
                    wcf(i, j)->coefficientRestitution = 0.7; // Denver Pilphis : hard coding need to be modified in the future
                    wcf(i, j)->coefficientRollingResistance = 0.01; // Denver Pilphis : hard coding need to be modified in the future
                    wcf(i, j)->shear_displacement = Vector3(0.0, 0.0, 0.0);
                    evaluate_wall(i, j);
                }
            }
        }
}

// Using CONTACT RADIUS of the spheres
// To determine whether a bond is assigned between two particles
void DEMSolver::bond_detect(const Integer& i, const Integer& j)
{
    Real contact_radius_i = gf[i]->radius * 1.1; // Denver Pilphis : hard coding need to be modified in the future
    Real contact_radius_j = gf[j]->radius * 1.1; // Denver Pilphis : hard coding need to be modified in the future
    if ((gf[j]->position - gf[i]->position).norm() - contact_radius_i - contact_radius_j < 0.0)
    {
        cf(i, j) = new Contact(); // Forced to bond contact
        // cf(i, j)->isActive = 1;
        cf(i, j)->isBonded = 1;
        // EBPM parameters
        // Denver Pilphis : hard coding need to be modified in the future
        cf(i, j)->radius_ratio = 0.5;
        cf(i, j)->elasticModulus = 28e9;
        cf(i, j)->poissonRatio = 0.2;
        cf(i, j)->compressiveStrength = 3e8;
        cf(i, j)->tensileStrength = 6e7;
        cf(i, j)->shearStrength = 6e7;
        cf(i, j)->force_a = Vector3(0.0, 0.0, 0.0);
        cf(i, j)->moment_a = Vector3(0.0, 0.0, 0.0);
        cf(i, j)->force_b = Vector3(0.0, 0.0, 0.0);
        cf(i, j)->moment_b = Vector3(0.0, 0.0, 0.0);
    }
}

// Similar to contact, but runs only once at the beginning
void DEMSolver::bond()
{
    Integer n = gf.rows();
    for (int i = 0; i < n; ++i)
        for (int j = i + 1; j < n; ++j)
            bond_detect(i, j);
}

// Neighboring search
void DEMSolver::contact()
{
    Integer n = gf.rows();
    for (int i = 0; i < n; ++i)
        for (int j = i + 1; j < n; ++j)
            resolve(i, j);
}

// Everything that needs to be done within a time step
void DEMSolver::run_simulation()
{
    clear_state();
    contact();
    resolve_wall();
    apply_body_force();
    update();
}

void DEMSolver::init_simulation()
{
    bond();
}

// Must control the contact status
Contact::Contact()
{
    // isActive = 0;
    isBonded = 0;
}
