#=##############################################################################
# DESCRIPTION
    Evaluate several manually defined timesteps of the fluid domain. by evaluating a single simulation

# AUTHORSHIP
  * Author          : Fynn Gerhardy
  * Email           : fygerh@gmail.com
  * Created         : Jan 2023
  * Last updated    : Jan 2023
  * License         : -
=###############################################################################
#include(joinpath("/home/fynn/Repositories/BladeWakeOPT/04_fluiddomain_visualization/postprocess_fluiddomain.jl"))

import FLOWUnsteady as uns
import FLOWVLM as vlm
import FLOWVPM as vpm

this_file_path = splitdir(@__FILE__)[1]
include(joinpath(this_file_path, "..", "functions", "OwnFunctions.jl"))  # Read file with module definition first
using .OwnFunctions                                                             # Include all self defined functions




# ----------------- 0) SPECIFY SIMULATION TO BE POSTPROCESSED AND FOLDER TO STORE DATA IN ------

# Folders and paths
sim_name = "NREL-5MW_10Revs_50BE_360steps"                                                                                      # Name of simulation to be evaluated (typically the last folder name of "simulation_path")
simulation_path = joinpath(this_file_path, "..", "01_single_turbine_simulation", "data_out", sim_name)   # Folder of simulation to be evaluated


call_paraview = true                                                              # call paraview after postprocessing

# turbine tip radius in (m)
Radius = 63.0

# Freestream reference wind speed
magVinf_fdom    = 11.4                                                             # (m/s) Magnitude of free stream wind speed
AOA_fdom        = 0.0                                                              # (deg) Angle of attack (incidence angle)
Vinf_fdom(X, t)      = magVinf_fdom*[cosd(AOA_fdom), sind(AOA_fdom), 0]            # wind speed in global coordinatesystem

# Time steps to evaluate
tstep_method    = "last_Rev"                                                       # (Options: "manual", "all", "last_Rev") tstep_method defines the timesteps to be calculated

# => when "manual", these timesteps will be calculated
#tsteps_manual = [2379]                                                         # set by specific timesteps manually chosen by the user
tsteps_manual = collect(2020:10:2380)                                                   # set by own start, step and end time (start:step:end)

# => following variables are necessary when using tstep_method = "all" or "last_Rev"
nrevs_fdom      = 40                                                               # number of revolutions the simulation was simulated with
nsteps_per_rev_fdom  = 72                                                          # number of steps per revolution the simulation was simulated with
stepwidth       = 6                                                                # set this to e.g. 2 if you want to calculate each second timestep, to 3 if you want to calculate each third timestep, ... and so on
n_lastRevs      = 1                                                               # number of last revolutions to be calculated (used if tstep_method = "last_Rev")

# grid to be calculated
calc_grid_x_y   = true
calc_grid_y_z   = true

gridsize_x_y    = 0.25      # grid size of x-y fluid domain plane in meters
gridsize_y_z    = 0.25      # grid size of y-z fluid domain plane in meters

z_locs          = [0]                                   # z coordinate location of x-y-plane = z_loc*2*Radius in meters
x_locs          = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]       # x coordinate location of y-z-plane = x_loc*2*Radius in meters

# x-y grid boundaries of x-y-plane (=> factor*2*Radius in meters) => Hub = coordinate origin
x_b_min_for_x_y = -0.2
y_b_min_for_x_y = -0.7
x_b_max_for_x_y = 13                 
y_b_max_for_x_y = 0.7

# y-z grid boundaries of y-z-plane (=> factor*2*Radius in meters) => Hub = coordinate origin
cylindrical_grid = true     # if true, a cylindrical grid with the size of the turbine radius is calculated (disregards the following variables)
y_b_min_for_y_z = -0.7
z_b_min_for_y_z = -0.7
y_b_max_for_y_z = 0.7
z_b_max_for_y_z = 0.7



# ----------------- 0) OVERWRITE VARIABLES IF THIS FILE WAS STARTED VIA THE "single_turbine_simulation" -------------------------------------------
# overwrite some variables, if this scrips gets executed via "01_single_turbine_simulation/single_turbine_simulation.jl"
if @isdefined run_name
    sim_name = run_name
end

if @isdefined save_path
    simulation_path = save_path
end

if @isdefined paraview                                                            
    call_paraview = paraview
end

if @isdefined R
    Radius = R
end

if @isdefined AOA
    AOA_fdom = AOA
end

if @isdefined magVinf
    magVinf_fdom = magVinf
    Vinf_fdom(X, t)      = magVinf_fdom*[cosd(AOA_fdom), sind(AOA_fdom), 0]
end

if @isdefined nrevs
    nrevs_fdom = nrevs
end

if @isdefined nsteps_per_rev
    nsteps_per_rev_fdom = nsteps_per_rev
end

if @isdefined all_tstep_method
    tstep_method = all_tstep_method
end

if @isdefined all_tstep_stepwidth
    stepwidth = all_tstep_stepwidth
end

if @isdefined all_tstep_n_lastRevs
    n_lastRevs = all_tstep_n_lastRevs
end

# set the foldername and save path automatically
    # Give the postprocessed fluid domain a name. 
    # => This will be the name of the folder all fluiddomain files will be stored in.
if tstep_method == "manual"
    folder_name_fluiddomain = sim_name*"-manual_tsteps"
elseif tstep_method == "last_Rev"
    folder_name_fluiddomain = sim_name*"-every_$(stepwidth)_tstep-last_$(n_lastRevs)_Revs"
elseif tstep_method == "all"
    folder_name_fluiddomain = sim_name*"-every_$(stepwidth)_tstep"                                                                                                                                                                           
end

save_path_fdom = joinpath(this_file_path, "data_out", folder_name_fluiddomain)

# ----------------- 1) SET TIMESTEPS TO BE EVALUATED -------------------------------------------
if tstep_method == "manual"
    tsteps = tsteps_manual

elseif tstep_method == "last_Rev"
    if n_lastRevs >= nrevs_fdom
        error("n_lastRevs is bigger then overall number of revolutions calculated for that Simulation!!!")
    end

    nsteps = nrevs_fdom*nsteps_per_rev_fdom # Number of time steps
    nsteps_start = nsteps - n_lastRevs*nsteps_per_rev_fdom
    tsteps = collect(nsteps_start-1:stepwidth:nsteps-1)

elseif tstep_method == "all"
    nsteps = nrevs_fdom*nsteps_per_rev_fdom # Number of time steps
    tsteps = collect(1:stepwidth:nsteps-1)
end

# ----------------- 2)POSTPROCESS PLANES -------------------------------------------------------
file_suffixes = []
if calc_grid_x_y
    for z in z_locs
        z_loc = z   # z coordinate location of plane = z_loc*2*Radius in meters
        file_suffix = "_x_y_atz$(z_loc)D"
        push!(file_suffixes, file_suffix)
        vol_thickness = 0                         # plane thickness in meters => set this value to 0 to get a 2D solution (plane)
        gridsize_m = gridsize_x_y                 # grid size in meters => 0.5 equals 0.5m

        # grid resolution
        x_res = Radius*(1/gridsize_m)
        y_res = Radius*(1/gridsize_m)
        z_res = 1
        # grid minimum boundaries (bound_factor*2*Radius in meters)
        z_b_min = (-(vol_thickness/2)/(2*Radius))+z_loc # choose volume that has a thickness of vol_thickness meter in the z-dimension (-0.5m under z_loc and +0.5m over z_loc)
        # grid maximum boundaries (bound_factor*2*Radius in meters)
        z_b_max = ((vol_thickness/2)/(2*Radius))+z_loc # choose volume that has a thickness of vol_thickness meter in the z-dimension (-0.5m under z_loc and +0.5m over z_loc)

        # calculate the plane
        if length(file_suffixes) == 1
            first_fdom = true
        else
            first_fdom = false
        end
        OwnFunctions.postprocess_fluiddomain(simulation_path, sim_name, file_suffix, Radius, AOA_fdom, tsteps;
                                             # ----- OPTIONAL ARGUMENTS ----------
                                             save_path       = save_path_fdom,               # folder to store data in (if nothing, data will be stored under simulation_path)
                                             # ----- FREESTREAM VELOCITY ---------
                                             Vinf            = Vinf_fdom,
                                             # ----- GRID OPTIONS ----------------  
                                             x_resolution    = x_res,                   # discretization of grid in x-direction
                                             y_resolution    = y_res,                   # discretization of grid in y-direction
                                             z_resolution    = z_res,                   # discretization of grid in z-direction
                                             x_bound_min     = x_b_min_for_x_y,         # minimum bounds in x-direction (bound_factor*2*Radius in meters)
                                             y_bound_min     = y_b_min_for_x_y,         # minimum bounds in y-direction (bound_factor*2*Radius in meters)
                                             z_bound_min     = z_b_min,                 # minimum bounds in z-direction (bound_factor*2*Radius in meters)
                                             x_bound_max     = x_b_max_for_x_y,         # maximum bounds in x-direction (bound_factor*2*Radius in meters)
                                             y_bound_max     = y_b_max_for_x_y,         # maximum bounds in y-direction (bound_factor*2*Radius in meters)
                                             z_bound_max     = z_b_max,                 # maximum bounds in z-direction (bound_factor*2*Radius in meters)
                                             prompt          = false,
                                             verbose         = true,
                                             debug           = false,                   # saves fdom grid as a file
                                             # ----- OTHER ------------------------- 
                                             first_fdom      = first_fdom               # prevents question if fluiddomain postprocessing folder should be removed or not
                                             )
    end
end

if calc_grid_y_z
    for x in x_locs # x represents the x-coordinate stations a y-z-plane will be calculated for
        x_loc = x   # x coordinate location of plane = x_loc*2*Radius in meters
        file_suffix = "_y_z_atx$(x_loc)D"
        push!(file_suffixes, file_suffix)
        vol_thickness = 0                         # plane thickness in meters => set this value to 0 to get a 2D solution (plane)
        gridsize_m = gridsize_y_z                 # grid size in meters => 0.5 equals 0.5m



        if cylindrical_grid
            # NOTE: The plane is created in the x_y plane first, then rotated into the y-z plane
                  # => therefore the x_loc is used to define the z boundaries (location in the z plane) and then rotated to -90 deg around the y axis
            
            # grid resolution
            x_res = Radius*(1/gridsize_m)
            y_res = 49                                                # Angular Step = 360/y_res in (°) !!!For some reason an error occurs when y_res >= 50!!! Nice would be: y_res = 2*pi*Radius*(1/gridsize_m) to have the minimum gridsize at the outer cylindrical plane section
            z_res = Radius*(1/gridsize_m)                                  # THIS WILL DO NOTHING IF vol_thickness == 0 (meaning z_b_min==z_b_max)
            # grid minimum boundaries (bound_factor*2*Radius in meters)
            x_b_min = 0.0
            #y_b_min_for_y_z = 0.0
            z_b_min_for_y_z = (-(vol_thickness/2)/(2*Radius))+x_loc
            # grid maximum boundaries (bound_factor*2*Radius in meters)
            x_b_max = 1.0
            #y_b_max_for_y_z = 2*pi
            z_b_max_for_y_z = ((vol_thickness/2)/(2*Radius))+x_loc
            
          else
            # grid resolution
            x_res = 1
            y_res = Radius*(1/gridsize_m)
            z_res = Radius*(1/gridsize_m)
            # grid minimum boundaries (bound_factor*2*Radius in meters)
            x_b_min = (-(vol_thickness/2)/(2*Radius))+x_loc # choose volume that has a thickness of vol_thickness meter in the x-dimension (-0.5m before x_loc and +0.5m behind x_loc)
            # grid maximum boundaries (bound_factor*2*Radius in meters)
            x_b_max = ((vol_thickness/2)/(2*Radius))+x_loc # choose volume that has a thickness of vol_thickness meter in the x-dimension (-0.5m before x_loc and +0.5m behind x_loc)

          end

        # calculate the plane
        if length(file_suffixes) == 1
            first_fdom = true
        else
            first_fdom = false
        end
        OwnFunctions.postprocess_fluiddomain(simulation_path, sim_name, file_suffix, Radius, AOA_fdom, tsteps;
                                             # ----- OPTIONAL ARGUMENTS ----------
                                             save_path       = save_path_fdom,               # folder to store data in (if nothing, data will be stored under simulation_path)
                                             # ----- FREESTREAM VELOCITY ---------
                                             Vinf            = Vinf_fdom,
                                             # ----- GRID OPTIONS ----------------  
                                             x_resolution    = x_res,                   # discretization of grid in x-direction
                                             y_resolution    = y_res,                   # discretization of grid in y-direction
                                             z_resolution    = z_res,                   # discretization of grid in z-direction
                                             x_bound_min     = x_b_min,                 # minimum bounds in x-direction (bound_factor*2*Radius in meters)
                                             y_bound_min     = y_b_min_for_y_z,         # minimum bounds in y-direction (bound_factor*2*Radius in meters)
                                             z_bound_min     = z_b_min_for_y_z,         # minimum bounds in z-direction (bound_factor*2*Radius in meters)
                                             x_bound_max     = x_b_max,                 # maximum bounds in x-direction (bound_factor*2*Radius in meters)
                                             y_bound_max     = y_b_max_for_y_z,         # maximum bounds in y-direction (bound_factor*2*Radius in meters)
                                             z_bound_max     = z_b_max_for_y_z,         # maximum bounds in z-direction (bound_factor*2*Radius in meters)
                                             cylindrical_grid=cylindrical_grid,         # a cylindrical grid will be calculated if true
                                             prompt          = false,
                                             verbose         = true,
                                             debug           = false,                   # saves fdom grid as a file
                                             # ----- OTHER ------------------------- 
                                             first_fdom      = first_fdom               # prevents question if fluiddomain postprocessing folder should be removed or not
                                             )
    end
end



# ----------------- 3) PARAVIEW ----------------------------------------------------

if call_paraview
    # --------------- save the suffixes of the _fdom files that will be opened in paraview -------------------
    if length(tsteps) == 1
        for i in 1:length(file_suffixes)
        file_suffixes[i] = file_suffixes[i]*"_fdom.$(tsteps[1]).xmf;"
        end
    elseif length(tsteps) > 1
    for i in 1:length(file_suffixes)
        file_suffixes[i] = file_suffixes[i]*"_fdom...xmf;"
    end
    end


    println("Calling Paraview...")
    files_fdom = joinpath(save_path_fdom, sim_name*file_suffixes[1])
    if length(file_suffixes)>1
        for suffix in file_suffixes[2:end]
            global files_fdom
            files_fdom *= sim_name*suffix
        end
    end
    println(files_fdom)

    # Call Paraview
    run(`paraview --data=$(files_fdom)`)

end
