/**
 * Self-gravitating disc
 *
 * A self-gravitating disc is integrated using
 * the leap frog integrator. Collisions are not resolved.
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "rebound.h"


void heartbeat(struct reb_simulation* const r);

double experiment_inclination_deg;
double experiment_satellite_central_mass;
double experiment_satellite_disk_mass;
int experiment_N_satellite;
double experiment_rp;

/*
 * Additional morphology-analysis state.
 *
 * initial_R_by_hash stores the initial cylindrical radius of each
 * primary-disk particle.  Hashes 1..10000 are used as stable particle IDs.
 */
#define N_PRIMARY_DISK 10000 
#define N_RADIAL_BINS 4
#define N_SATELLITE_DISK (N_PRIMARY_DISK/10)

static double initial_R_by_hash[N_PRIMARY_DISK + 1];
static int initial_R_initialized = 0;

void add_satellite(
    struct reb_simulation* const r,
    double inclination,
    double satellite_central_mass,
    double satellite_disk_mass,
    int N_satellite,
    double rp
){

    // Physical size of satellite disk
    double a_min = 0.5;     // Minimum semi-major axis of satellite disk particles 
    double a_max = 1.8;     // Maximum semi-major axis of satellite disk particles

    // Initial orbital position
    double x0 = 8.0;
    double y0 = 0.0;
    double z0 = 0.0;

    // Total masses of primary and satellite
    double primary_mass = 1.0 + 0.2;
    double satellite_mass = satellite_central_mass + satellite_disk_mass;

    // Initial orbital distance
    double r0 = x0;

    // Total mass of the two galaxies
    double total_mass = primary_mass + satellite_mass;

    // Total speed for a parabolic orbit at r0
    double v_orbit = sqrt(
    2.0 * r->G * total_mass / r0
    );

    // Angular momentum for a parabolic orbit
    // with periapsis rp
    double h = sqrt(
        2.0 * r->G * total_mass * rp
    );

    // Tangential component of velocity
    double v_t = h / r0;

    // Radial component of velocity
    double v_r = sqrt(
        v_orbit*v_orbit - v_t*v_t
    );

    // Satellite initially moves inward toward the primary
    double vx_orbit = -v_r;

    // Incline the tangential velocity
    double vy_orbit = v_t * cos(inclination);
    double vz_orbit = v_t * sin(inclination);

    // Add satellite central mass
    struct reb_particle satellite_star = {0};

    satellite_star.m = satellite_central_mass;
    satellite_star.hash = 10001;

    satellite_star.x = x0;
    satellite_star.y = y0;
    satellite_star.z = z0;

    satellite_star.vx = vx_orbit;
    satellite_star.vy = vy_orbit;
    satellite_star.vz = vz_orbit;

    reb_simulation_add(r, satellite_star);

    // Add satellite disk particles
    for (int i = 0; i < N_satellite; i++){
        struct reb_particle pt = {0};
        pt.hash = 10002 + i;

        // Internal position within satellite disk
        double a = reb_random_powerlaw(
            r,
            a_min,
            a_max,
            -1.5
        );

        double phi = reb_random_uniform(
            r,
            0,
            2.0 * M_PI
        );

        pt.x = a * cos(phi) + x0;
        pt.y = a * sin(phi) + y0;
        pt.z = a * reb_random_normal(r, 0.001) + z0;

        // Satellite internal gravitational mass
        double mu = satellite_star.m
            + satellite_disk_mass *
            (
                pow(a, -3.0/2.0) - pow(a_min, -3.0/2.0)
            )
            /
            (
                pow(a_max, -3.0/2.0) - pow(a_min, -3.0/2.0)
            );

        // Internal circular velocity
        double vkep = sqrt(r->G * mu / a);

        // Internal disk velocity
        pt.vx = vkep * sin(phi) + vx_orbit;
        pt.vy = -vkep * cos(phi) + vy_orbit;
        pt.vz = vz_orbit;

        // Particle mass
        pt.m = satellite_disk_mass / (double)N_satellite;

        reb_simulation_add(r, pt);
    }
}

int main(int argc, char* argv[]){
    struct reb_simulation* const r = reb_simulation_create();
    
    r->rand_seed = 12345;   // Set the random seed for reproducibility

    // Start the REBOUND visualization server. This
    // allows you to visualize the simulation by pointing 
    // your web browser to http://localhost:1234
    //reb_simulation_start_server(r, 1234);

    const double t_end = 100;

    // Setup constants
    r->integrator       = REB_INTEGRATOR_LEAPFROG; // Leapfrog integrator
    r->gravity          = REB_GRAVITY_TREE;        // Tree code gravity
    r->boundary         = REB_BOUNDARY_OPEN;       // Open boundary conditions
    r->opening_angle2   = 1.5;          // This constant determines the accuracy of the tree code gravity estimate.
    r->G                = 1;            // Gravitational constant
    r->softening        = 0.02;         // Gravitational softening length
    r->dt               = 1.5e-2;         // Timestep
    const double boxsize = 80.0;        // Size of the simulation box
    reb_simulation_configure_box(r,boxsize,1,1,1);  // Configure the simulation box for the tree code

    // Setup particles
    double disc_mass = 2e-1;    // Total disc mass
    int N = N_PRIMARY_DISK;            // Number of particles
    // Initial conditions
    struct reb_particle star = {0};
    star.m         = 1;
    reb_simulation_add(r, star);
    for (int i=0;i<N;i++){
        struct reb_particle pt = {0};
        pt.hash = i + 1;
        double a_min = 1.02;
        double a_max = 4.25;
        double a = reb_random_powerlaw(r, a_min, a_max, -1.5);
        double phi     = reb_random_uniform(r, 0,2.*M_PI);
        pt.x         = a*cos(phi);
        pt.y         = a*sin(phi);
        pt.z         = a*reb_random_normal(r, 0.001);
        double mu = star.m + disc_mass *(pow(a,-3./2.)-pow(a_min,-3./2.)) /(pow(a_max,-3./2.)-pow(a_min,-3./2.));
        double vkep     = sqrt(r->G*mu/a);
        pt.vx         =  vkep * sin(phi);
        pt.vy         = -vkep * cos(phi);
        pt.vz         = 0;
        pt.m         = disc_mass/(double)N;
        reb_simulation_add(r, pt);
        
    }
    // Satellite Galaxy / experiment parameters
    // Change these values when setting up a new test run.

    double inclination_deg = 45.0;        // Orbital inclination in degrees
    double satellite_central_mass = 0.1; // Satellite central mass
    double satellite_disk_mass = 0.02;   // Satellite disk mass
    int N_satellite = N_SATELLITE_DISK;               // Satellite disk particles
    double rp = 2.5;                      // Periapsis / closest approach

    double inclination = inclination_deg * M_PI / 180.0;

    // Store experiment parameters for the heartbeat output.
    experiment_inclination_deg = inclination_deg;
    experiment_satellite_central_mass = satellite_central_mass;
    experiment_satellite_disk_mass = satellite_disk_mass;
    experiment_N_satellite = N_satellite;
    experiment_rp = rp;

    add_satellite(
        r,
        inclination,
        satellite_central_mass,
        satellite_disk_mass,
        N_satellite,
        rp
    );

    r->heartbeat = heartbeat;
    reb_simulation_integrate(r, t_end);
}

void heartbeat(struct reb_simulation* const r){
    if (reb_simulation_output_check(r, 1.0)){

        /*
         * Print the experiment parameters once at the start.
         */
        if (r->t < 1.0e-12){
            double satellite_total_mass =
                experiment_satellite_central_mass +
                experiment_satellite_disk_mass;

            printf("# primary_mass = 1.200000\n");
            printf("# satellite_central_mass = %.6f\n",
                   experiment_satellite_central_mass);
            printf("# satellite_disk_mass = %.6f\n",
                   experiment_satellite_disk_mass);
            printf("# satellite_total_mass = %.6f\n",
                   satellite_total_mass);
            printf("# satellite_disk_particles = %d\n",
                   experiment_N_satellite);
            printf("# mass_ratio = %.6f\n",
                   1.2 / satellite_total_mass);
            printf("# inclination_deg = %.1f\n",
                   experiment_inclination_deg);
            printf("# rp = %.2f\n",
                   experiment_rp);

            printf("# radial_bins = 0-1, 1-2, 2-3, 3-4.25\n");

            printf("# columns:\n");
            printf("# 1 t\n");
            printf("# 2 N_primary_present\n");
            printf("# 3 f_primary_present\n");
            printf("# 4 N_core\n");
            printf("# 5 f_core\n");
            printf("# 6 R_rms_all\n");
            printf("# 7 z_rms_all\n");
            printf("# 8 R_rms_core\n");
            printf("# 9 z_rms_core\n");
            printf("# 10 z_core_cm\n");
            printf("# 11 E\n");
            printf("# 12 Lz_all\n");
            printf("# 13 x_cm\n");
            printf("# 14 z_cm\n");
            printf("# 15 deltaR_mean_abs\n");
            printf("# 16 deltaR_rms\n");
            printf("# 17 A1_all\n");
            printf("# 18 z_rms_0_1\n");
            printf("# 19 z_rms_1_2\n");
            printf("# 20 z_rms_2_3\n");
            printf("# 21 z_rms_3_4.25\n");
            printf("# 22 A1_0_1\n");
            printf("# 23 A1_1_2\n");
            printf("# 24 A1_2_3\n");
            printf("# 25 A1_3_4.25\n");
        }

        const int N_primary_disk = N_PRIMARY_DISK;
        const double R_core = 2.5;
        const int N_core_min = 500;

        /*
         * Radial bins for spatially resolved morphology.
         *
         * These cover the original primary disk.  Particles scattered
         * beyond 4.25 are still included in the global diagnostics and
         * radial-migration diagnostics.
         */
        const double R_bin_edges[N_RADIAL_BINS + 1] =
            {0.0, 1.0, 2.0, 3.0, 4.25};

        /*
         * ---------------------------------------------------------
         * 1. Find primary-disk particles and their centre of mass.
         * ---------------------------------------------------------
         */
        int N_primary_present = 0;

        double total_mass = 0.0;
        double x_cm = 0.0;
        double y_cm = 0.0;
        double z_cm = 0.0;

        for (int i = 0; i < r->N; i++){

            struct reb_particle p = r->particles[i];

            if (p.hash >= 1 && p.hash <= N_primary_disk){

                N_primary_present++;

                total_mass += p.m;

                x_cm += p.m * p.x;
                y_cm += p.m * p.y;
                z_cm += p.m * p.z;
            }
        }

        /*
         * If somehow no primary particles remain, stop the
         * diagnostic calculation.
         */
        if (N_primary_present == 0){
            printf("%f 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0\n",
                   r->t);
            return;
        }

        x_cm /= total_mass;
        y_cm /= total_mass;
        z_cm /= total_mass;

        /*
         * ---------------------------------------------------------
         * 2. Store the initial cylindrical radius of each primary
         *    particle the first time the heartbeat runs.
         *
         *    We use the primary-disk COM at t=0 as the reference.
         * ---------------------------------------------------------
         */
        if (!initial_R_initialized){

            for (int i = 0; i < r->N; i++){

                struct reb_particle p = r->particles[i];

                if (p.hash >= 1 && p.hash <= N_primary_disk){

                    double dx = p.x - x_cm;
                    double dy = p.y - y_cm;

                    initial_R_by_hash[p.hash] =
                        sqrt(dx*dx + dy*dy);
                }
            }

            initial_R_initialized = 1;
        }

        /*
         * ---------------------------------------------------------
         * 3. Diagnostics for the entire primary disk.
         * ---------------------------------------------------------
         */
        double R2_sum_all = 0.0;
        double z2_sum_all = 0.0;
        double Lz_all = 0.0;

        /*
         * Radial-migration diagnostics:
         *
         * deltaR = R_current - R_initial.
         *
         * Mean signed deltaR can cancel between inward and outward
         * migration, so we record both mean absolute displacement
         * and RMS displacement.
         */
        double deltaR_abs_sum = 0.0;
        double deltaR2_sum = 0.0;

        /*
         * Global m=1 lopsidedness:
         *
         * A1 = |sum exp(i phi)| / N.
         *
         * Because all primary disk particles have equal mass,
         * particle counting is sufficient here.
         */
        double A1_cos_all = 0.0;
        double A1_sin_all = 0.0;

        /*
         * Radially resolved z_rms and A1.
         */
        int N_bin[N_RADIAL_BINS] = {0, 0, 0, 0};
        double z2_sum_bin[N_RADIAL_BINS] = {0.0, 0.0, 0.0, 0.0};
        double A1_cos_bin[N_RADIAL_BINS] = {0.0, 0.0, 0.0, 0.0};
        double A1_sin_bin[N_RADIAL_BINS] = {0.0, 0.0, 0.0, 0.0};

        for (int i = 0; i < r->N; i++){

            struct reb_particle p = r->particles[i];

            if (p.hash >= 1 && p.hash <= N_primary_disk){

                double dx = p.x - x_cm;
                double dy = p.y - y_cm;
                double dz = p.z - z_cm;

                double R2 = dx*dx + dy*dy;
                double R = sqrt(R2);

                R2_sum_all += R2;
                z2_sum_all += dz*dz;

                Lz_all += p.m * (
                    dx*p.vy - dy*p.vx
                );

                /*
                 * Radial displacement from the particle's own
                 * initial radius.
                 */
                double deltaR = R - initial_R_by_hash[p.hash];

                deltaR_abs_sum += fabs(deltaR);
                deltaR2_sum += deltaR*deltaR;

                /*
                 * Global m=1 mode.
                 */
                if (R > 0.0){
                    double cos_phi = dx / R;
                    double sin_phi = dy / R;

                    A1_cos_all += cos_phi;
                    A1_sin_all += sin_phi;
                }

                /*
                 * Radially resolved measurements.
                 */
                for (int b = 0; b < N_RADIAL_BINS; b++){

                    if (R >= R_bin_edges[b] &&
                        R < R_bin_edges[b+1]){

                        N_bin[b]++;
                        z2_sum_bin[b] += dz*dz;

                        if (R > 0.0){
                            A1_cos_bin[b] += dx / R;
                            A1_sin_bin[b] += dy / R;
                        }

                        break;
                    }
                }
            }
        }

        double R_rms_all = sqrt(
            R2_sum_all / N_primary_present
        );

        double z_rms_all = sqrt(
            z2_sum_all / N_primary_present
        );

        double deltaR_mean_abs =
            deltaR_abs_sum / N_primary_present;

        double deltaR_rms =
            sqrt(deltaR2_sum / N_primary_present);

        double A1_all =
            sqrt(
                A1_cos_all*A1_cos_all +
                A1_sin_all*A1_sin_all
            ) / N_primary_present;

        /*
         * Radial-bin z_rms and A1.
         *
         * If a bin has no particles, output -1.
         */
        double z_rms_bin[N_RADIAL_BINS];
        double A1_bin[N_RADIAL_BINS];

        for (int b = 0; b < N_RADIAL_BINS; b++){

            z_rms_bin[b] = -1.0;
            A1_bin[b] = -1.0;

            if (N_bin[b] > 0){

                z_rms_bin[b] =
                    sqrt(z2_sum_bin[b] / N_bin[b]);

                A1_bin[b] =
                    sqrt(
                        A1_cos_bin[b]*A1_cos_bin[b] +
                        A1_sin_bin[b]*A1_sin_bin[b]
                    ) / N_bin[b];
            }
        }

        /*
         * ---------------------------------------------------------
         * 4. Find primary particles inside the core.
         *
         * Core definition:
         *
         *                 R < 2.5
         *
         * R is measured relative to the primary disk COM.
         * ---------------------------------------------------------
         */
        int N_core = 0;

        double core_mass = 0.0;

        double x_core_cm = 0.0;
        double y_core_cm = 0.0;
        double z_core_cm = 0.0;

        for (int i = 0; i < r->N; i++){

            struct reb_particle p = r->particles[i];

            if (p.hash >= 1 && p.hash <= N_primary_disk){

                double dx = p.x - x_cm;
                double dy = p.y - y_cm;

                double R2 = dx*dx + dy*dy;

                if (R2 < R_core*R_core){

                    N_core++;

                    core_mass += p.m;

                    x_core_cm += p.m*p.x;
                    y_core_cm += p.m*p.y;
                    z_core_cm += p.m*p.z;
                }
            }
        }

        /*
         * ---------------------------------------------------------
         * 5. Core centre of mass.
         * ---------------------------------------------------------
         */
        if (N_core > 0){

            x_core_cm /= core_mass;
            y_core_cm /= core_mass;
            z_core_cm /= core_mass;
        }

        /*
         * ---------------------------------------------------------
         * 6. Core RMS quantities.
         * ---------------------------------------------------------
         */
        double R2_sum_core = 0.0;
        double z2_sum_core = 0.0;

        if (N_core >= N_core_min){

            for (int i = 0; i < r->N; i++){

                struct reb_particle p = r->particles[i];

                if (p.hash >= 1 && p.hash <= N_primary_disk){

                    double dx_global = p.x - x_cm;
                    double dy_global = p.y - y_cm;

                    double R2_global =
                        dx_global*dx_global +
                        dy_global*dy_global;

                    if (R2_global < R_core*R_core){

                        double dx_core =
                            p.x - x_core_cm;

                        double dy_core =
                            p.y - y_core_cm;

                        double dz_core =
                            p.z - z_core_cm;

                        double R2_core =
                            dx_core*dx_core +
                            dy_core*dy_core;

                        R2_sum_core += R2_core;
                        z2_sum_core += dz_core*dz_core;
                    }
                }
            }
        }

        /*
         * ---------------------------------------------------------
         * 7. Core RMS values.
         *
         * If fewer than 500 primary particles are in the core,
         * mark the measurement as invalid with -1.
         * ---------------------------------------------------------
         */
        double R_rms_core = -1.0;
        double z_rms_core = -1.0;

        if (N_core >= N_core_min){

            R_rms_core = sqrt(
                R2_sum_core / N_core
            );

            z_rms_core = sqrt(
                z2_sum_core / N_core
            );
        }

        /*
         * ---------------------------------------------------------
         * 8. Fraction of original primary disk inside core.
         * ---------------------------------------------------------
         */
        double f_core =
            (double)N_core /
            (double)N_primary_disk;

        /*
         * ---------------------------------------------------------
         * 9. Fraction of original primary disk still present.
         * ---------------------------------------------------------
         */
        double f_primary_present =
            (double)N_primary_present /
            (double)N_primary_disk;

        /*
         * ---------------------------------------------------------
         * 10. Total energy.
         * ---------------------------------------------------------
         */
        double E = reb_simulation_energy(r);

        /*
         * ---------------------------------------------------------
         * 11. Output.
         *
         * The new diagnostics are:
         *
         * deltaR_mean_abs:
         *     Mean absolute change in each particle's cylindrical
         *     radius from its own initial radius.
         *
         * deltaR_rms:
         *     RMS radial displacement from the initial radius.
         *
         * A1_all:
         *     Global m=1 lopsidedness.
         *
         * z_rms_0_1 ... z_rms_3_4.25:
         *     Vertical thickness in four radial bins.
         *
         * A1_0_1 ... A1_3_4.25:
         *     m=1 lopsidedness in the same radial bins.
         * ---------------------------------------------------------
         */
        printf(
            "%f %d %e %d %e %e %e %e %e %e "
            "%e %e %e %e %e %e %e "
            "%e %e %e %e %e %e %e %e\n",
            r->t,
            N_primary_present,
            f_primary_present,
            N_core,
            f_core,
            R_rms_all,
            z_rms_all,
            R_rms_core,
            z_rms_core,
            z_core_cm,
            E,
            Lz_all,
            x_cm,
            z_cm,
            deltaR_mean_abs,
            deltaR_rms,
            A1_all,
            z_rms_bin[0],
            z_rms_bin[1],
            z_rms_bin[2],
            z_rms_bin[3],
            A1_bin[0],
            A1_bin[1],
            A1_bin[2],
            A1_bin[3]
        );
    }
}
