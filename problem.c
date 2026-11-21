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

static double initial_R_by_hash[N_PRIMARY_DISK + 1];
static int initial_R_initialized = 0;

/* Initial values for whole-system conservation diagnostics. */
static double E_initial = 0.0;
static struct reb_vec3d L_initial;
static double Px_initial = 0.0;
static double Py_initial = 0.0;
static double Pz_initial = 0.0;
static double P_scale_initial = 0.0;
static FILE *conservation_file = NULL;
static FILE *orbit_file = NULL;

/* Running minima of the two galaxy-separation measurements. */
static double min_central_separation = HUGE_VAL;
static double min_com_separation = HUGE_VAL;
static double time_min_central_separation = 0.0;
static double time_min_com_separation = 0.0;

/* Calculate total linear momentum and a characteristic momentum scale. */
void calculate_total_momentum(
    struct reb_simulation* const r,
    double *Px,
    double *Py,
    double *Pz,
    double *Pscale
){
    *Px = 0.0;
    *Py = 0.0;
    *Pz = 0.0;
    *Pscale = 0.0;

    for (int i = 0; i < r->N; i++){
        struct reb_particle p = r->particles[i];

        *Px += p.m * p.vx;
        *Py += p.m * p.vy;
        *Pz += p.m * p.vz;

        double speed = sqrt(
            p.vx*p.vx + p.vy*p.vy + p.vz*p.vz
        );
        *Pscale += p.m * speed;
    }
}


/*
 * Measure galaxy separation in two ways:
 * 1. Between the two central particles (hashes 0 and 10001).
 * 2. Between the mass-weighted COMs of all currently present particles
 *    belonging to each galaxy.
 */
int calculate_galaxy_separations(
    struct reb_simulation* const r,
    double *d_central,
    double *d_com
){
    double xp = 0.0, yp = 0.0, zp = 0.0, mp = 0.0;
    double xs = 0.0, ys = 0.0, zs = 0.0, ms = 0.0;

    double xpc = 0.0, ypc = 0.0, zpc = 0.0;
    double xsc = 0.0, ysc = 0.0, zsc = 0.0;
    int found_pc = 0, found_sc = 0;

    for (int i = 0; i < r->N; i++){
        struct reb_particle p = r->particles[i];

        if (p.hash == 0 ||
            (p.hash >= 1 && p.hash <= N_PRIMARY_DISK)){
            mp += p.m;
            xp += p.m * p.x;
            yp += p.m * p.y;
            zp += p.m * p.z;
        }

        if (p.hash == 10001 ||
            (p.hash >= 10002 &&
             p.hash < 10002 + experiment_N_satellite)){
            ms += p.m;
            xs += p.m * p.x;
            ys += p.m * p.y;
            zs += p.m * p.z;
        }

        if (p.hash == 0){
            xpc = p.x; ypc = p.y; zpc = p.z;
            found_pc = 1;
        }
        if (p.hash == 10001){
            xsc = p.x; ysc = p.y; zsc = p.z;
            found_sc = 1;
        }
    }

    if (mp <= 0.0 || ms <= 0.0 || !found_pc || !found_sc){
        return 0;
    }

    xp /= mp; yp /= mp; zp /= mp;
    xs /= ms; ys /= ms; zs /= ms;

    double dxc = xsc - xpc;
    double dyc = ysc - ypc;
    double dzc = zsc - zpc;

    double dxm = xs - xp;
    double dym = ys - yp;
    double dzm = zs - zp;

    *d_central = sqrt(dxc*dxc + dyc*dyc + dzc*dzc);
    *d_com = sqrt(dxm*dxm + dym*dym + dzm*dzm);

    return 1;
}

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

    /*
     * Nominal two-body (Keplerian) pericentre.
     * The velocity is initialized using the total galaxy masses and
     * the initial central-position separation r0. In the extended,
     * evolving N-body system, the measured minimum central-particle
     * or galaxy-COM separation need not equal this nominal rp.
     */
    // Angular momentum for a parabolic orbit with nominal pericentre rp
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
    r->opening_angle2   = 0.25;          // This constant determines the accuracy of the tree code gravity estimate.(1.5)
    r->G                = 1;            // Gravitational constant
    r->softening        = 0.02;         // Gravitational softening length
    r->dt               = 0.015;         // Timestep
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

    double inclination_deg = 45.0;        // Orbital inclination in degrees 45, 0, 90
    double satellite_central_mass = 0.1; // Satellite central mass 1:10, 1:5, 1:20
    double satellite_disk_mass = 0.02;   // Satellite disk mass 1:10, 1:5, 1:20
    const double satellite_particle_mass = 2.0e-5;
    int N_satellite = (int)llround(
        satellite_disk_mass / satellite_particle_mass
    );  // Keeps satellite disk particle mass consistent across mass ratios

    double rp = 2.5;  // Nominal Keplerian pericentre used to set initial velocity

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

    /* Establish whole-system conservation baselines after all particles
       have been added and before the first integration step. */
    E_initial = reb_simulation_energy(r);
    L_initial = reb_simulation_angular_momentum(r);

    calculate_total_momentum(
        r,
        &Px_initial,
        &Py_initial,
        &Pz_initial,
        &P_scale_initial
    );

    conservation_file = fopen("conservation.txt", "w");
    if (conservation_file == NULL){
        fprintf(stderr, "Error opening conservation.txt\n");
        exit(EXIT_FAILURE);
    }

    fprintf(conservation_file,
        "# t E dE_frac Lx Ly Lz dL_abs dL_frac Px Py Pz dP_abs dP_frac N\n"
    );

    orbit_file = fopen("orbit_separation.txt", "w");
    if (orbit_file == NULL){
        fprintf(stderr, "Error opening orbit_separation.txt\n");
        exit(EXIT_FAILURE);
    }
    fprintf(orbit_file,
        "# t d_central d_com min_d_central t_min_central min_d_com t_min_com\n"
    );

    r->heartbeat = heartbeat;
    reb_simulation_integrate(r, t_end);

    if (conservation_file != NULL){
        fclose(conservation_file);
        conservation_file = NULL;
    }

    if (orbit_file != NULL){
        fclose(orbit_file);
        orbit_file = NULL;
    }

    reb_simulation_free(r);
}

void heartbeat(struct reb_simulation* const r){

    /*
     * Measure orbit separations every integration step. This gives a
     * better estimate of the true minimum than sampling only every 1.0
     * time unit. Morphology output remains at its original 1.0 cadence.
     */
    double d_central = -1.0;
    double d_com = -1.0;

    if (calculate_galaxy_separations(r, &d_central, &d_com)){
        if (d_central < min_central_separation){
            min_central_separation = d_central;
            time_min_central_separation = r->t;
        }
        if (d_com < min_com_separation){
            min_com_separation = d_com;
            time_min_com_separation = r->t;
        }

        if (orbit_file != NULL){
            fprintf(orbit_file,
                "%.6f %.12e %.12e %.12e %.6f %.12e %.6f\n",
                r->t, d_central, d_com,
                min_central_separation, time_min_central_separation,
                min_com_separation, time_min_com_separation
            );
            fflush(orbit_file);
        }
    }

    if (reb_simulation_output_check(r, 1.0)){

        /* Whole-system conservation diagnostics.  These include both
           galaxies, unlike the primary-only morphology quantities below. */
        double E = reb_simulation_energy(r);
        double energy_error =
            (E - E_initial) / fmax(fabs(E_initial), 1e-12);

        struct reb_vec3d L = reb_simulation_angular_momentum(r);

        double dLx = L.x - L_initial.x;
        double dLy = L.y - L_initial.y;
        double dLz = L.z - L_initial.z;
        double dL = sqrt(dLx*dLx + dLy*dLy + dLz*dLz);

        double L_initial_mag = sqrt(
            L_initial.x*L_initial.x +
            L_initial.y*L_initial.y +
            L_initial.z*L_initial.z
        );
        double L_error = dL / fmax(L_initial_mag, 1e-12);

        double Px, Py, Pz, Pscale;
        calculate_total_momentum(r, &Px, &Py, &Pz, &Pscale);
        (void)Pscale;  /* The normalization uses the initial momentum scale. */

        double dPx = Px - Px_initial;
        double dPy = Py - Py_initial;
        double dPz = Pz - Pz_initial;
        double dP = sqrt(dPx*dPx + dPy*dPy + dPz*dPz);
        double P_error = dP / fmax(P_scale_initial, 1e-12);

        if (conservation_file != NULL){
            fprintf(conservation_file,
                "%.6f %.12e %.12e "
                "%.12e %.12e %.12e %.12e %.12e "
                "%.12e %.12e %.12e %.12e %.12e %d\n",
                r->t, E, energy_error,
                L.x, L.y, L.z, dL, L_error,
                Px, Py, Pz, dP, P_error, r->N
            );
            fflush(conservation_file);
        }

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
            printf("# rp_nominal_keplerian = %.2f\n",
                   experiment_rp);
            printf("# rp_definition = nominal two-body pericentre used for initial velocity\n");

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
            printf("# 12 Lz_primary_internal\n");
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
            printf("# 26 dx_cm_core\n");
            printf("# 27 dy_cm_core\n");
            printf("# 28 dz_cm_core\n");
            printf("# 29 dcm_core\n");
            printf("# 30 A2_all\n");
            printf("# 31 A2_0_1\n");
            printf("# 32 A2_1_2\n");
            printf("# 33 A2_2_3\n");
            printf("# 34 A2_3_4.25\n");
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
        double vx_cm = 0.0;
        double vy_cm = 0.0;
        double vz_cm = 0.0;

        for (int i = 0; i < r->N; i++){

            struct reb_particle p = r->particles[i];

            if (p.hash >= 1 && p.hash <= N_primary_disk){

                N_primary_present++;

                total_mass += p.m;

                x_cm += p.m * p.x;
                y_cm += p.m * p.y;
                z_cm += p.m * p.z;
                vx_cm += p.m * p.vx;
                vy_cm += p.m * p.vy;
                vz_cm += p.m * p.vz;
            }
        }

        /*
         * If somehow no primary particles remain, stop the
         * diagnostic calculation.
         */
        if (N_primary_present == 0){
            printf("%f "
                   "0 0 0 0 0 0 0 0 0 0 "
                   "0 0 0 0 0 0 0 0 0 0 "
                   "0 0 0 0 0 0 0 0 0 0 "
                   "0 0 0\n",
                   r->t);
            return;
        }

        x_cm /= total_mass;
        y_cm /= total_mass;
        z_cm /= total_mass;
        vx_cm /= total_mass;
        vy_cm /= total_mass;
        vz_cm /= total_mass;

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
         * Global m=2 two-fold asymmetry:
         * A1 = |sum exp(i phi)| / N.
         *
         * Because all primary disk particles have equal mass,
         * particle counting is sufficient here.
         */
        double A1_cos_all = 0.0;
        double A1_sin_all = 0.0;
        double A2_cos_all = 0.0;
        double A2_sin_all = 0.0;
        /*
         * Radially resolved z_rms, A1, and A2.
         */
        int N_bin[N_RADIAL_BINS] = {0, 0, 0, 0};
        double z2_sum_bin[N_RADIAL_BINS] = {0.0, 0.0, 0.0, 0.0};
        double A1_cos_bin[N_RADIAL_BINS] = {0.0, 0.0, 0.0, 0.0};
        double A1_sin_bin[N_RADIAL_BINS] = {0.0, 0.0, 0.0, 0.0};
        double A2_cos_bin[N_RADIAL_BINS] = {0.0, 0.0, 0.0, 0.0};
        double A2_sin_bin[N_RADIAL_BINS] = {0.0, 0.0, 0.0, 0.0};

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
                /* Internal angular momentum about the primary COM. */
                Lz_all += p.m * (
                    dx*(p.vy - vy_cm) - dy*(p.vx - vx_cm)
                );

                double deltaR = R - initial_R_by_hash[p.hash];
                deltaR_abs_sum += fabs(deltaR);
                deltaR2_sum += deltaR*deltaR;

                if (R > 0.0){
                    double cos_phi = dx / R;
                    double sin_phi = dy / R;
                    A1_cos_all += cos_phi;
                    A1_sin_all += sin_phi;
                    A2_cos_all += cos_phi*cos_phi - sin_phi*sin_phi;
                    A2_sin_all += 2.0*cos_phi*sin_phi;
                }

                /* Radial-bin calculations must remain inside this
                   primary-particle block, where R, dx, dy, and dz exist. */
                for (int b = 0; b < N_RADIAL_BINS; b++){
                    if (R >= R_bin_edges[b] &&
                        R < R_bin_edges[b+1]){

                        N_bin[b]++;
                        z2_sum_bin[b] += dz*dz;

                        if (R > 0.0){
                            double cos_phi = dx / R;
                            double sin_phi = dy / R;
                            A1_cos_bin[b] += cos_phi;
                            A1_sin_bin[b] += sin_phi;
                            A2_cos_bin[b] += cos_phi*cos_phi - sin_phi*sin_phi;
                            A2_sin_bin[b] += 2.0*cos_phi*sin_phi;
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
        
        double A2_all =
            sqrt(
                A2_cos_all*A2_cos_all +
                A2_sin_all*A2_sin_all
            ) / N_primary_present;
        /*
         * Radial-bin z_rms, A1, and A2.
         *
         * If a bin has no particles, output -1.
         */
        double z_rms_bin[N_RADIAL_BINS];
        double A1_bin[N_RADIAL_BINS];
        double A2_bin[N_RADIAL_BINS];

        for (int b = 0; b < N_RADIAL_BINS; b++){

            z_rms_bin[b] = -1.0;
            A1_bin[b] = -1.0;
            A2_bin[b] = -1.0;
            
            if (N_bin[b] > 0){

                z_rms_bin[b] =
                    sqrt(z2_sum_bin[b] / N_bin[b]);

                A1_bin[b] =
                    sqrt(
                        A1_cos_bin[b]*A1_cos_bin[b] +
                        A1_sin_bin[b]*A1_sin_bin[b]
                    ) / N_bin[b];

                A2_bin[b] =
                    sqrt(
                        A2_cos_bin[b]*A2_cos_bin[b] +
                        A2_sin_bin[b]*A2_sin_bin[b]
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
         * 5b. Compare whole-primary COM with core COM.
         *
         * These measure how far the centre of the whole primary
         * disk has moved relative to its central core.
         * ---------------------------------------------------------
         */
        double dx_cm_core = -1.0;
        double dy_cm_core = -1.0;
        double dz_cm_core = -1.0;
        double dcm_core   = -1.0;

        if (N_core >= N_core_min){

            dx_cm_core = x_cm - x_core_cm;
            dy_cm_core = y_cm - y_core_cm;
            dz_cm_core = z_cm - z_core_cm;

            dcm_core = sqrt(
            dx_cm_core*dx_cm_core +
            dy_cm_core*dy_cm_core +
            dz_cm_core*dz_cm_core
    );
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
         * 10. Output.
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
         * A2_all:
         *     Global m=2 lopsidedness.
         *
         * z_rms_0_1 ... z_rms_3_4.25:
         *     Vertical thickness in four radial bins.
         *
         * A1_0_1 ... A1_3_4.25:
         *     m=1 lopsidedness in the same radial bins.
         * 
         * A2_0_1 ... A2_3_4.25:
         *     m=2 lopsidedness in the same radial bins.
         * ---------------------------------------------------------
         */
        printf(
            "%f %d %e %d "
            "%e %e %e %e %e %e %e %e %e %e "
            "%e %e %e %e %e %e %e %e %e %e "
            "%e %e %e %e %e %e %e %e %e %e\n",
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
            A1_bin[3],
            dx_cm_core,
            dy_cm_core,
            dz_cm_core,
            dcm_core,
            A2_all,
            A2_bin[0],
            A2_bin[1],
            A2_bin[2],
            A2_bin[3]
        );
    }
}
