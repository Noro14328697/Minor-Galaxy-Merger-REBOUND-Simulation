# Galaxy Merger Simulations with REBOUND

This repository contains the modified `problem.c` source file and its
`Makefile` for running collisionless galaxy-merger simulations using
[REBOUND](https://rebound.readthedocs.io/).

The code is intended to be used inside REBOUND's existing
**`examples/selfgravitating_disk/`** example directory. It is not a
standalone REBOUND installation.

## Installation and placement

1.  Download or clone the REBOUND source code and make sure its
    dependencies are installed.

2.  Locate the `examples/selfgravitating_disk/` directory in your
    REBOUND checkout.

3.  Place the repository's `problem.c` and `Makefile` in that directory,
    replacing the example files if prompted.

4.  From inside `examples/selfgravitating_disk/`, build the program
    using:

    ``` bash
    make
    ```

5.  Run the executable produced by the Makefile. For the standard
    REBOUND example target, this is typically:

    ``` bash
    ./problem
    ```

    If your Makefile uses a different executable name, run that target
    instead.

Run the program from the example directory so that its output files are
created there.

## Simulation overview

The code models an encounter between a primary galaxy and a satellite
galaxy using a collisionless, self-gravitating N-body simulation.

The simulation uses:

-   Newtonian gravity with (G=1).
-   REBOUND's leapfrog integrator.
-   REBOUND's tree-gravity calculation.
-   Open boundary conditions.
-   Gravitational softening length of `0.02`.
-   Timestep of `0.015` simulation time units.
-   Tree opening-angle parameter `opening_angle2 = 0.25`.
-   A fixed random seed (`12345`) for reproducible initial particle
    realizations.
-   A total integration time of `100` simulation time units.

The primary galaxy consists of a central particle of mass `1.0` and a
disk of total mass `0.2`, represented by `10,000` particles.

The satellite contains a central particle and a disk. Its disk particle
count is calculated from the selected disk mass so that the satellite
disk particle mass remains `2.0e-5` across the mass-ratio experiments.

The satellite begins at a separation of `8.0` simulation length units.
Its initial velocity is set using a nominal two-body parabolic orbit.
The input `rp` is the **nominal Keplerian pericentre used to calculate
the initial velocity**; the actual closest approach in the extended,
evolving N-body system is measured during the simulation and need not
equal this input value.

The satellite disk is initially in the (xy)-plane. The inclination
parameter changes the direction of the satellite's bulk orbital velocity
relative to the primary disk plane.

## Setting the experiment parameters

The main experiment parameters are set near the beginning of `main()` in
`problem.c`:

``` c
double inclination_deg = 0.0;
double satellite_central_mass = 0.1;
double satellite_disk_mass = 0.02;
double rp = 2.5;
```

Change these values before compiling to run a different experiment.

The parameter sweeps used in the associated project are:

  -----------------------------------------------------------------------
  Experiment              Values varied           Other parameters held
                                                  fixed
  ----------------------- ----------------------- -----------------------
  Orbital inclination     `0°`, `45°`, `90°`      Mass ratio `1:10`,
                                                  nominal `rp = 2.5`

  Satellite-to-primary    `1:20`, `1:10`, `1:5`   Inclination `45°`,
  mass ratio                                      nominal `rp = 2.5`

  Nominal pericentre      `1.5`, `2.5`, `3.5`     Inclination `45°`, mass
                                                  ratio `1:10`
  -----------------------------------------------------------------------

For the mass-ratio experiments, change both satellite mass parameters
consistently:

  -----------------------------------------------------------------------
  Total mass      Satellite central     Satellite disk     Satellite disk
  ratio                        mass               mass          particles
  -------------- ------------------ ------------------ ------------------
  `1:20`                     `0.05`             `0.01`              `500`

  `1:10`                     `0.10`             `0.02`             `1000`

  `1:5`                      `0.20`             `0.04`             `2000`
  -----------------------------------------------------------------------

The primary galaxy's total mass is `1.2`. The satellite total masses in
the table are `0.06`, `0.12`, and `0.24`, respectively.

After changing parameters, rebuild with `make` before running the
simulation again.

## Requirements

-   A working REBOUND source checkout.
-   A C compiler supported by the REBOUND build.
-   The REBOUND dependencies required by the parent project.
-   The `problem.c` and `Makefile` from this repository, placed in
    `examples/selfgravitating_disk/`.

For REBOUND installation and build details, consult the [official
REBOUND documentation](https://rebound.readthedocs.io/).
