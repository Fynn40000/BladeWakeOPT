#=##############################################################################
# DESCRIPTION
    Simulation of a single wind turbine rotor.

# AUTHORSHIP
  * Author          : Fynn Gerhardy
  * Email           : fygerh@gmail.com
  * Created         : Dec 2023
  * Last updated    : Dec 2023
  * License         : -
=###############################################################################
#include(joinpath("/home/fynn/Repositories/BladeWakeOPT/01_single_turbine_simulation/single_turbine_simulation.jl"))

import FLOWUnsteady as uns
import FLOWVLM as vlm
import FLOWVPM as vpm

start_simulation_path = splitdir(@__FILE__)[1]
include(joinpath(start_simulation_path, "..", "functions", "OwnFunctions.jl"))  # Read file with module definition first
using .OwnFunctions                                                             # Include all self defined functions

turbine_name    = "NREL-5MW"

# ----------------- Fidelity Options --------------------------------------------------------

fidelity        = "mid"                     # options: "low", "mid", "high"
run_length      = 2                        # number of revolutions to run => defines the length of the simulation

cut_wake_mode   = "plane"                   # cut all particles away 
x_loc           = 12                        # x location from wich on the particles will be cut away (x coordinate = x_loc*2*R in meters)

# ----------------- Postprocessing and Visualization ----------------------------------------
postprocessing  = true                      # perform postprocessing in general???
paraview        = false                     # Whether to visualize with Paraview
plot_bladeloads = true                      # postprocess the blade loads and plot the radial distribution
postprocess_fdom= true                      # postprocess the fluid domains last timestep and calculate velocity field etc.
debug           = true                      # enables calculation of coefficients such as cn, ct, cl, cd
show_bladeload_plots = false                # show the bladeload plots on display after postprocessing?
postprocess_all_tsteps_fdom = false          # postprocess all timesteps of the fluid domain?
                                            # => NOTE: if true, you need to set up the script that gets called 
                                            #          via this file with your desired inputs (see "6) POSTPROCESSING")

all_tstep_method = "last_Rev"               # method to use when postprocessing several timesteps of the fluiddomain
all_tstep_stepwidth = 72                     # take each ...th timestep
all_tstep_n_lastRevs = 20                    # number of (last) revolutions that will be postprocessed

# ----------------- GEOMETRY PARAMETERS ----------------------------------------

# Rotor geometry
rotor_file      = "NREL5MW.csv"                                             # Rotor geometry
data_path       = joinpath(start_simulation_path, "..", "00_database")      # Path to rotor database
pitch           = 0.0                                                       # (deg) collective pitch of blades
CW              = true                                                      # Clock-wise rotation
xfoil           = false                                                     # Whether to run XFOIL
ncrit           = 9                                                         # Turbulence criterion for XFOIL
turbine_flag    = true                                                      # This rotor is a turbine

# NOTE: If `xfoil=true`, XFOIL will be run to generate the airfoil polars used
#       by blade elements before starting the simulation. XFOIL is run
#       on the airfoil contours found in `rotor_file` at the corresponding
#       local Reynolds and Mach numbers along the blade.
#       Alternatively, the user can provide pre-computer airfoil polars using
#       `xfoil=false` and pointing to polar files through `rotor_file`.

# Discretization
if fidelity == "low"
    n               = 50#20                         # Number of blade elements per blade
elseif fidelity == "mid"
    n               = 50                                                
elseif fidelity == "high"
    n               = 50                         
end
r               = 1/10                       # Geometric expansion of elements               

# NOTE: Here a geometric expansion of 1/5 means that the spacing between the
#       tip elements is 1/5 of the spacing between the hub elements. Refine the
#       discretization towards the blade tip like this in order to better
#       resolve the tip vortex.

# Read radius of this rotor and number of blades
R, B            = uns.read_rotor(rotor_file; data_path=data_path)[[1,3]]

# Set angle of attack boundaries manually here 
# (this is an optional input when generating a rotor and will only calculate the 
# specified boundaries for the (rediscretized) control point airfoil locations/positions)
aoa_bounds = uns.DataFrames.DataFrame(uns.CSV.File(joinpath(data_path, "rotors", "NREL5MW_aoa_bounds.csv")))
aoa_bounds = Array{Float64, 2}(aoa_bounds)

# ----------------- SIMULATION PARAMETERS --------------------------------------

# Operating conditions
RPM             = 12.1*8/11.4#12.1                      # RPM
J               = 8/((RPM/60)*2*R)#11.4/((RPM/60)*2*R)       # Advance ratio Vinf/(nD)
AOA             = 0.0                       # (deg) Angle of attack (incidence angle)

rho             = 1.225                     # (kg/m^3) air density
mu              = 1.81e-5                   # (kg/ms) air dynamic viscosity
speedofsound    = 342.35                    # (m/s) speed of sound

magVinf         = J*RPM/60*(2*R)
magVinfx        = magVinf*cosd(AOA)         # wind velocity in x direction
Vinf(X, t)      = magVinf*[cosd(AOA), sind(AOA), 0] # (m/s) freestream velocity vector

ReD             = 2*pi*RPM/60*R * rho/mu * 2*R      # Diameter-based Reynolds number
Matip           = 2*pi*RPM/60*R / speedofsound      # Tip Mach number

println("""
    RPM:    $(RPM)
    Vinf:   $(Vinf(zeros(3), 0)) m/s
    tsr:    $((2*pi*RPM/60*R)/(J*((RPM/60)*2*R)))
    Matip:  $(round(Matip, digits=3))
    ReD:    $(round(ReD, digits=0))
""")

# ----------------- SOLVER PARAMETERS ------------------------------------------

# Aerodynamic solver
VehicleType     = uns.UVLMVehicle           # Unsteady solver
# VehicleType   = uns.QVLMVehicle           # Quasi-steady solver => for low fidelity simulation (uses blade-element momentum theory?!?!)
const_solution  = VehicleType==uns.QVLMVehicle  # Whether to assume that the
                                                #  solution is constant or not
# Time parameters
nrevs           = run_length                # Number of revolutions in simulation

if fidelity == "low"
    nsteps_per_rev  = 36                    # Time steps per revolution
elseif fidelity == "mid"
    nsteps_per_rev  = 72
elseif fidelity == "high"
    nsteps_per_rev  = 360
end

nsteps          = const_solution ? 2 : nrevs*nsteps_per_rev # Number of time steps
ttot            = nsteps/nsteps_per_rev / (RPM/60)       # (s) total simulation time

# VPM particle shedding
if fidelity == "low"
    p_per_step      = 2                      # Sheds per time step
    shed_starting   = false                  # Whether to shed starting vortex
elseif fidelity == "mid"
    p_per_step      = 2
    shed_starting   = false
elseif fidelity == "high"
    p_per_step      = 2
    shed_starting   = true
end

shed_unsteady   = true                        # Whether to shed vorticity from unsteady loading
unsteady_shedcrit = 0.001                     # Shed unsteady loading whenever circulation
                                              #  fluctuates by more than this ratio (default = 0.01)
max_particles   = ((2*n+1)*B)*nsteps*p_per_step + 1 # Maximum number of particles

# Regularization
if fidelity == "low"
    sigma_rotor_surf= R/10                      # Rotor-on-VPM smoothing radius
    sigmafactor_vpmonvlm= 1                     # Shrink particles by this factor when
                                                #  calculating VPM-on-VLM/Rotor induced velocities
elseif fidelity == "mid"
    sigma_rotor_surf= R/10
    sigmafactor_vpmonvlm= 1                                   
elseif fidelity == "high"
    sigma_rotor_surf= R/80
    sigmafactor_vpmonvlm= 5.5                   
end

lambda_vpm      = 2.125                     # VPM core overlap
                                            # VPM smoothing radius
sigma_vpm_overwrite = lambda_vpm * 2*pi*R/(nsteps_per_rev*p_per_step)

# Rotor solver
vlm_rlx         = 0.5                       # VLM relaxation <-- this also applied to rotors
hubtiploss_correction = vlm.hubtiploss_correction_prandtl#vlm.hubtiploss_nocorrection# Hub and tip loss correction (options: vlm.hubtiploss_nocorrection, 
                                                                                           #vlm.hubtiploss_correction_prandtl
                                                                                           #vlm.hubtiploss_correction_modprandtl)

# VPM solver
if fidelity == "low"
    vpm_integration = vpm.euler                 # VPM temporal integration scheme (=> vpm.rungekutta3 = default)
    vpm_SFS         = vpm.SFS_none              # VPM LES subfilter-scale model
elseif fidelity == "mid"
    vpm_integration = vpm.rungekutta3  
    vpm_SFS         = vpm.SFS_none
elseif fidelity == "high"
    vpm_integration = vpm.rungekutta3
    vpm_SFS = vpm.SFS_Cd_twolevel_nobackscatter
end
# vpm_SFS       = vpm.SFS_Cd_twolevel_nobackscatter
# vpm_SFS       = vpm.SFS_Cd_threelevel_nobackscatter
# vpm_SFS       = vpm.DynamicSFS(vpm.Estr_fmm, vpm.pseudo3level_positive;
#                                   alpha=0.999, maxC=1.0,
#                                   clippings=[vpm.clipping_backscatter])
# vpm_SFS       = vpm.DynamicSFS(vpm.Estr_fmm, vpm.pseudo3level_positive;
#                                   alpha=0.999, rlxf=0.005, minC=0, maxC=1
#                                   clippings=[vpm.clipping_backscatter],
#                                   controls=[vpm.control_sigmasensor],
#                                   )

vpm_viscous     = vpm.Inviscid()            # VPM viscous diffusion scheme

# NOTE: In most practical situations, open rotors operate at a Reynolds number
#       high enough that viscous diffusion in the wake is negligible.
#       Hence, it does not make much of a difference whether we run the
#       simulation with viscous diffusion enabled or not.


if VehicleType == uns.QVLMVehicle
    # NOTE: If the quasi-steady solver is used, this mutes warnings regarding
    #       potential colinear vortex filaments. This is needed since the
    #       quasi-steady solver will probe induced velocities at the lifting
    #       line of the blade
    uns.vlm.VLMSolver._mute_warning(true)
end


# ----------------- WAKE TREATMENT ---------------------------------------------
# NOTE: It is known in the CFD community that rotor simulations with an
#       impulsive RPM start (*i.e.*, 0 to RPM in the first time step, as opposed
#       to gradually ramping up the RPM) leads to the hub "fountain effect",
#       with the root wake reversing the flow near the hub.
#       The fountain eventually goes away as the wake develops, but this happens
#       very slowly, which delays the convergence of the simulation to a steady
#       state. To accelerate convergence, here we define a wake treatment
#       procedure that suppresses the hub wake for the first three revolutions,
#       avoiding the fountain effect altogether.
#       This is especially helpful in low and mid-fidelity simulations.

if fidelity == "low"
    suppress_fountain   = false                  # Toggle
elseif fidelity == "mid"
    suppress_fountain   = false
elseif fidelity == "high"
    suppress_fountain   = false
end

# Supress wake shedding on blade elements inboard of this r/R radial station
no_shedding_Rthreshold = suppress_fountain ? 0.35 : 0.0


# Supress wake shedding for this many time steps
no_shedding_nstepsthreshold = 3*nsteps_per_rev

omit_shedding = []          # Index of blade elements to supress wake shedding

# Function to suppress or activate wake shedding
function wake_treatment_supress(sim, args...; optargs...)

    # Case: start of simulation -> suppress shedding
    if sim.nt == 1

        # Identify blade elements on which to suppress shedding
        for i in 1:vlm.get_m(rotor)
            HS = vlm.getHorseshoe(rotor, i)
            CP = HS[5]

            if uns.vlm.norm(CP - vlm._get_O(rotor)) <= no_shedding_Rthreshold*R
                push!(omit_shedding, i)
            end
        end
    end

    # Case: sufficient time steps -> enable shedding
    if sim.nt == no_shedding_nstepsthreshold

        # Flag to stop suppressing
        omit_shedding .= -1

    end

    return false
end


# ----------------- 0) SET THE SIMULATION NAME AND SAVE PATH --------------------------------------
run_name        = turbine_name*"_$(nrevs)Revs_$(n)BE_$(nsteps_per_rev)steps"    # Name of this simulation
save_path       = joinpath(start_simulation_path, "data_out", run_name)         # Where to save this simulation           
save_path_post  = joinpath(save_path, "postprocessing")                         # Where to save postprocessing plots


# ----------------- 1) VEHICLE DEFINITION --------------------------------------
println("Generating geometry...")

# Generate rotor
rotor = uns.generate_rotor(rotor_file; pitch=pitch,
                                        n=n, CW=CW, blade_r=r,
                                        altReD=[RPM, J, mu/rho],
                                        xfoil=xfoil,
                                        ncrit=ncrit,
                                        turbine_flag=turbine_flag,
                                        aoa_bounds=aoa_bounds,
                                        data_path=data_path,
                                        verbose=true,
                                        verbose_xfoil=false,
                                        plot_disc=true,
                                        save_polars=save_path
                                        );


println("Generating vehicle...")

# Generate vehicle
system = vlm.WingSystem()                   # System of all FLOWVLM objects
vlm.addwing(system, "Rotor", rotor)

rotors = [rotor];                           # Defining this rotor as its own system
rotor_systems = (rotors, );                 # All systems of rotors

wake_system = vlm.WingSystem()              # System that will shed a VPM wake
                                            # NOTE: Do NOT include rotor when using the quasi-steady solver
if VehicleType != uns.QVLMVehicle
    vlm.addwing(wake_system, "Rotor", rotor)
end

vehicle = VehicleType(   system;
                            rotor_systems=rotor_systems,
                            wake_system=wake_system
                         );

# NOTE: Through the `rotor_systems` keyword argument to `uns.VLMVehicle` we
#       have declared any systems (groups) of rotors that share a common RPM.
#       We will later declare the control inputs to each rotor system when we
#       define the `uns.KinematicManeuver`.


# ------------- 2) MANEUVER DEFINITION -----------------------------------------
# Non-dimensional translational velocity of vehicle over time
Vvehicle(t) = zeros(3)

# Angle of the vehicle over time
anglevehicle(t) = zeros(3)

# RPM control input over time (RPM over `RPMref`)
RPMcontrol(t) = 1.0

angles = ()                                 # Angle of each tilting system (none)
RPMs = (RPMcontrol, )                       # RPM of each rotor system

maneuver = uns.KinematicManeuver(angles, RPMs, Vvehicle, anglevehicle)

# NOTE: `FLOWUnsteady.KinematicManeuver` defines a maneuver with prescribed
#       kinematics. `Vvehicle` defines the velocity of the vehicle (a vector)
#       over time. `anglevehicle` defines the attitude of the vehicle over time.
#       `angle` defines the tilting angle of each tilting system over time.
#       `RPM` defines the RPM of each rotor system over time.
#       Each of these functions receives a nondimensional time `t`, which is the
#       simulation time normalized by the total time `ttot`, from 0 to
#       1, beginning to end of simulation. They all return a nondimensional
#       output that is then scaled by either a reference velocity (`Vref`) or
#       a reference RPM (`RPMref`). Defining the kinematics and controls of the
#       maneuver in this way allows the user to have more control over how fast
#       to perform the maneuver, since the total time, reference velocity and
#       RPM are then defined in the simulation parameters shown below.


# ------------- 3) SIMULATION DEFINITION ---------------------------------------

Vref = 0.0                                  # Reference velocity to scale maneuver by
RPMref = RPM                                # Reference RPM to scale maneuver by

Vinit = Vref*Vvehicle(0)                    # Initial vehicle velocity
Winit = pi/180*(anglevehicle(1e-6) - anglevehicle(0))/(1e-6*ttot)  # Initial angular velocity

simulation = uns.Simulation(vehicle, maneuver, Vref, RPMref, ttot;
                                                    Vinit=Vinit, Winit=Winit);


# ------------- 4) MONITORS DEFINITIONS ----------------------------------------
figs, figaxs = [], []                       # Figures generated by monitor

#=
# Generate turbine monitor
println("Generating monitor...")
monitor_rotor = uns.generate_monitor_turbines(rotors, J, rho, RPM, nsteps, magVinfx, turbine_flag;
                                            t_scale=RPM/60,        # Scaling factor for time in plots
                                            t_lbl="Revolutions",   # Label for time axis
                                            out_figs=figs,
                                            out_figaxs=figaxs,
                                            save_path=save_path,
                                            run_name=run_name,
                                            figname="turbine monitor",
                                            )
=#

monitor_rotor = OwnFunctions.generate_monitor_turbines(rotors, J, rho, RPM, nsteps, magVinfx, turbine_flag;
                                            t_scale=RPM/60,        # Scaling factor for time in plots
                                            t_lbl="Revolutions",   # Label for time axis
                                            out_figs=figs,
                                            out_figaxs=figaxs,
                                            save_path=save_path,
                                            run_name=run_name,
                                            figname="turbine monitor",
                                            cut_wake_mode=cut_wake_mode,
                                            cut_wake_loc=[1/(2*R*x_loc), 0, 0]
                                            )


# ------------- 5) RUN SIMULATION ----------------------------------------------
println("Running simulation...")

# Concatenate monitors and wake treatment procedure into one runtime function
runtime_function = uns.concatenate(monitor_rotor, wake_treatment_supress)

uns.run_simulation(simulation, nsteps;
                    # ----- SIMULATION OPTIONS -------------
                    Vinf=Vinf,
                    rho=rho, mu=mu, sound_spd=speedofsound,
                    # ----- SOLVERS OPTIONS ----------------
                    p_per_step=p_per_step,
                    max_particles=max_particles,
                    vpm_integration=vpm_integration,
                    vpm_viscous=vpm_viscous,
                    vpm_SFS=vpm_SFS,
                    sigma_vlm_surf=sigma_rotor_surf,
                    sigma_rotor_surf=sigma_rotor_surf,
                    sigma_vpm_overwrite=sigma_vpm_overwrite,
                    sigmafactor_vpmonvlm=sigmafactor_vpmonvlm, # (experimental) shrinks the particles by this factor when calculating VPM-on-VLM/Rotor induced velocities
                    vlm_vortexsheet = false, # Whether to spread surface circulation as a vortex sheet in the VPM (turns ASM on; ALM if false)
                    vlm_rlx=vlm_rlx, # VLM relaxation (>0.9 can cause divergence, <0.2 slows simulation too much, deactivated with <0)
                    hubtiploss_correction=hubtiploss_correction, # Hub and tip loss correction of rotors (ignored in quasi-steady solver)
                    wake_coupled = true,          # true => VLM is used; false => BEM is used
                    shed_starting=shed_starting,
                    shed_unsteady=shed_unsteady,
                    unsteady_shedcrit=unsteady_shedcrit,
                    omit_shedding=omit_shedding,
                    extra_runtime_function=runtime_function,#monitor_rotor,
                    # ----- OUTPUT OPTIONS ------------------
                    save_path=save_path,
                    run_name=run_name,
                    debug=debug
                    );

# Save the calculated aoa boundaries as csv files if desired
OwnFunctions._create_csv(rotor._r, rotor._aoa_bound_min, save_path, "aoa_min_bound";
                         x_name="r(m)",
                         y_name="AOA_min(°)"
                         )
OwnFunctions._create_csv(rotor._r, rotor._aoa_bound_max, save_path, "aoa_max_bound";
                         x_name="r(m)",
                         y_name="AOA_max(°)"
                         )

# ------------- 6) POSTPROCESSING ----------------------------------------------
if postprocessing
    println("\nPostprocessing...")
    println("\n     => Blade loads and last time step fluiddomain...\n")
    #Post-process monitor plots
    fdom_suffixes = OwnFunctions.single_turbine_simulation_postprocessing(
                                                            save_path, save_path_post, run_name, R, AOA;
                                                            # ----- POSTPROCESSING EXECUTION -------------
                                                            plot_bladeloads=plot_bladeloads,       # postprocessing the bladeloads?
                                                            postprocess_fdom=postprocess_fdom,     # postprocessing the fluiddomain?
                                                            # ----- SETTINGS FOR POSTPROCESSING -------------
                                                            Vinf = Vinf,                           # Freestream Velocity
                                                            magVinfx = magVinfx,                   # Freestream Velocity in x direction
                                                            sim_time = ttot,                       # Overall (real) simulation time in seconds
                                                            rev_to_average_idx=nrevs,              # Revolution to wich the postprocessing should be applied on
                                                            nrevs_to_average=1,                    # number of Revolutions to average for postprocessing the bladeloads
                                                            num_elements=n,                        # number of blade elements per blade
                                                            tsteps = [nsteps-1],                   # time steps to be postprocessed
                                                            debug=debug,                           # postprocess dimensionless coefficients too? => NOTE: debug statement must be set to true for uns.run_simulation. Otherwise the simulation files will not contain the coefficient data.
                                                            suppress_plots=!show_bladeload_plots,  # suppresses the plots to show up on the display
                                                            gridsize_x_y=0.25,                     # grid size of x-y fluid domain plane in meters
                                                            gridsize_y_z=0.25,                     # grid size of y-z fluid domain plane in meters
                                                            cylindrical_grid = true,               # if true, the y-z plane will be calculated as a cylindrical grid and the wake velocity profiles will be saved within a .csv file
                                                                                                   # this grid will be set automatically with the turbine diameter as its diameter
                                                            verbose = false,
                                                            cylinder_radius = 2                    # set the cylinder radius (factor*R)
                                                            )

    # postprocess all timesteps via ths script and call paraview to visualize them...
    if postprocess_all_tsteps_fdom
        println("\n     => All time steps fluiddomain...\n")
        include(joinpath("..", "04_fluiddomain_visualization", "postprocess_fluiddomain.jl"))
    end

    # ----------------- 7) VISUALIZATION -------------------------------------------
    if paraview
        println("Calling Paraview...")

        # Files to open in Paraview
        files = joinpath(save_path, run_name*"_pfield...xmf;")
        for bi in 1:B
            global files
            files *= run_name*"_Rotor_Blade$(bi)_loft...vtk;"
            files *= run_name*"_Rotor_Blade$(bi)_vlm...vtk;"
        end

        if postprocess_fdom
            for suffix in fdom_suffixes
                global files
                files *= run_name*suffix
            end
        end


        # Call Paraview
        run(`paraview --data=$(files)`)

    end

    
end

println("\n\nSIMULATION DONE!\nPROCESS EXITED NORMALLY!")
