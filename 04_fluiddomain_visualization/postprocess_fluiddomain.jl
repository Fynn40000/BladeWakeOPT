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
simulation_path = joinpath(this_file_path, "..", "01_single_turbine_simulation", "data_out", "NREL5MW_turbine_simulation")   # Folder of simulation to be evaluated
sim_name = "NREL5MW_turbine_simulation"                                                                                      # Name of simulation to be evaluated (typically the last folder name of "simulation_path")
folder_name_fluiddomain = sim_name*"-lowfid_15Revs_ev2nd_tstep"                                                              # Give the postprocessed fluid domain a name. 
                                                                                                                             # => This will be the name of the folder all fluiddomain files will be stored in.
save_path_fdom = joinpath(this_file_path, "data_out", folder_name_fluiddomain)                                                    # Folder to store postprocessed fluiddomain files in


call_paraview = true                                                              # call paraview after postprocessing
if @isdefined paraview                                                            # overwrite the call_paraview variable if this scrips gets executed via the single turbine simulation script
    call_paraview = paraview
end

# turbine tip radius in (m)
R = 63.0

# Freestream reference wind speed
magVinf         = 11.4                                                             # (m/s) Magnitude of free stream wind speed
AOA             = 0.0                                                              # (deg) Angle of attack (incidence angle)
Vinf(X, t)      = magVinf*[cosd(AOA), sind(AOA), 0]                                # wind speed in global coordinatesystem

# Time steps to evaluate
tstep_method    = "manual"                                                            # tstep_method defines the timesteps to be calculated (Options: "manual", "all")
                                                                                   # => when "manual", jump to next section and specify the timesteps to be evaluated manually
# => following variables are necessary when using tstep_method = "all"
nrevs           = 15                                                               # number of revolutions the simulation was simulated with
nsteps_per_rev  = 36                                                               # number of steps per revolution the simulation was simulated with
stepwidth       = 2                                                                # set this to e.g. 2 if you want to calculate each second timestep, to 3 if you want to calculate each third timestep, ... and so on

# grid to be calculated
calc_grid_x_y   = true
calc_grid_y_z   = false

gridsize_x_y    = 0.25      # grid size of x-y fluid domain plane in meters
gridsize_y_z    = 0.25      # grid size of y-z fluid domain plane in meters

z_locs          = [0]       # z coordinate location of plane = z_loc*2*R in meters
x_locs          = [1]       # x coordinate location of plane = x_loc*2*R in meters

# x-y grid boundaries of x-y-plane (=> factor*2*R in meters) => Hub = coordinate origin
x_b_min_for_x_y = -0.2
y_b_min_for_x_y = -0.7
x_b_max_for_x_y = 10                 
y_b_max_for_x_y = 0.7

# x-y grid boundaries of y-z-plane (=> factor*2*R in meters) => Hub = coordinate origin
cylindrical_grid = true     # if true, a cylindrical grid with the size of the turbine radius is calculated
y_b_min_for_y_z = -0.7
z_b_min_for_y_z = -0.7
y_b_max_for_y_z = 0.7
z_b_max_for_y_z = 0.7



# ----------------- 1) SET TIMESTEPS TO BE EVALUATED -------------------------------------------
if tstep_method == "manual"

    tsteps = [10,15,20] # set by specific timesteps
    #tsteps = collect(1:2:10)      # set by own start, step and end time step (start:step:end)

elseif tstep_method == "all"
    nsteps = nrevs*nsteps_per_rev # Number of time steps
    tsteps = collect(1:stepwidth:nsteps-1)
end

# ----------------- 2)POSTPROCESS PLANES -------------------------------------------------------
file_suffixes = []
if calc_grid_x_y
    for z in z_locs
        z_loc = z   # z coordinate location of plane = z_loc*2*R in meters
        file_suffix = "_x_y_atz$(z_loc)D"
        push!(file_suffixes, file_suffix)
        vol_thickness = 0                         # plane thickness in meters => set this value to 0 to get a 2D solution (plane)
        gridsize_m = gridsize_x_y                 # grid size in meters => 0.5 equals 0.5m

        # grid resolution
        x_res = R*(1/gridsize_m)
        y_res = R*(1/gridsize_m)
        z_res = 1
        # grid minimum boundaries (bound_factor*2*R in meters)
        z_b_min = (-(vol_thickness/2)/(2*R))+z_loc # choose volume that has a thickness of vol_thickness meter in the z-dimension (-0.5m under z_loc and +0.5m over z_loc)
        # grid maximum boundaries (bound_factor*2*R in meters)
        z_b_max = ((vol_thickness/2)/(2*R))+z_loc # choose volume that has a thickness of vol_thickness meter in the z-dimension (-0.5m under z_loc and +0.5m over z_loc)

        # calculate the plane
        if length(file_suffixes) == 1
            first_fdom = true
        else
            first_fdom = false
        end
        OwnFunctions.postprocess_fluiddomain(simulation_path, sim_name, file_suffix, R, AOA, tsteps;
                                             # ----- OPTIONAL ARGUMENTS ----------
                                             save_path       = save_path_fdom,               # folder to store data in (if nothing, data will be stored under simulation_path)
                                             # ----- FREESTREAM VELOCITY ---------
                                             Vinf            = Vinf,
                                             # ----- GRID OPTIONS ----------------  
                                             x_resolution    = x_res,                   # discretization of grid in x-direction
                                             y_resolution    = y_res,                   # discretization of grid in y-direction
                                             z_resolution    = z_res,                   # discretization of grid in z-direction
                                             x_bound_min     = x_b_min_for_x_y,         # minimum bounds in x-direction (bound_factor*2*R in meters)
                                             y_bound_min     = y_b_min_for_x_y,         # minimum bounds in y-direction (bound_factor*2*R in meters)
                                             z_bound_min     = z_b_min,                 # minimum bounds in z-direction (bound_factor*2*R in meters)
                                             x_bound_max     = x_b_max_for_x_y,         # maximum bounds in x-direction (bound_factor*2*R in meters)
                                             y_bound_max     = y_b_max_for_x_y,         # maximum bounds in y-direction (bound_factor*2*R in meters)
                                             z_bound_max     = z_b_max,                 # maximum bounds in z-direction (bound_factor*2*R in meters)
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
        x_loc = x   # x coordinate location of plane = x_loc*2*R in meters
        file_suffix = "_y_z_atx$(x_loc)D"
        push!(file_suffixes, file_suffix)
        vol_thickness = 0                         # plane thickness in meters => set this value to 0 to get a 2D solution (plane)
        gridsize_m = gridsize_y_z                 # grid size in meters => 0.5 equals 0.5m



        if cylindrical_grid
            # NOTE: The plane is created in the x_y plane first, then rotated into the y-z plane
                  # => therefore the x_loc is used to define the z boundaries (location in the z plane) and then rotated to -90 deg around the y axis
            
            # grid resolution
            x_res = R*(1/gridsize_m)
            y_res = 49                                                # Angular Step = 360/y_res in (°) !!!For some reason an error occurs when y_res >= 50!!! Nice would be: y_res = 2*pi*R*(1/gridsize_m) to have the minimum gridsize at the outer cylindrical plane section
            z_res = R*(1/gridsize_m)                                  # THIS WILL DO NOTHING IF vol_thickness == 0 (meaning z_b_min==z_b_max)
            # grid minimum boundaries (bound_factor*2*R in meters)
            x_b_min = 0.0
            #y_b_min_for_y_z = 0.0
            z_b_min_for_y_z = (-(vol_thickness/2)/(2*R))+x_loc
            # grid maximum boundaries (bound_factor*2*R in meters)
            x_b_max = 1.0
            #y_b_max_for_y_z = 2*pi
            z_b_max_for_y_z = ((vol_thickness/2)/(2*R))+x_loc
            

    
          else
            # grid resolution
            x_res = 1
            y_res = R*(1/gridsize_m)
            z_res = R*(1/gridsize_m)
            # grid minimum boundaries (bound_factor*2*R in meters)
            x_b_min = (-(vol_thickness/2)/(2*R))+x_loc # choose volume that has a thickness of vol_thickness meter in the x-dimension (-0.5m before x_loc and +0.5m behind x_loc)
            # grid maximum boundaries (bound_factor*2*R in meters)
            x_b_max = ((vol_thickness/2)/(2*R))+x_loc # choose volume that has a thickness of vol_thickness meter in the x-dimension (-0.5m before x_loc and +0.5m behind x_loc)

          end

        # calculate the plane
        if length(file_suffixes) == 1
            first_fdom = true
        else
            first_fdom = false
        end
        OwnFunctions.postprocess_fluiddomain(simulation_path, sim_name, file_suffix, R, AOA, tsteps;
                                             # ----- OPTIONAL ARGUMENTS ----------
                                             save_path       = save_path_fdom,               # folder to store data in (if nothing, data will be stored under simulation_path)
                                             # ----- FREESTREAM VELOCITY ---------
                                             Vinf            = Vinf,
                                             # ----- GRID OPTIONS ----------------  
                                             x_resolution    = x_res,                   # discretization of grid in x-direction
                                             y_resolution    = y_res,                   # discretization of grid in y-direction
                                             z_resolution    = z_res,                   # discretization of grid in z-direction
                                             x_bound_min     = x_b_min,                 # minimum bounds in x-direction (bound_factor*2*R in meters)
                                             y_bound_min     = y_b_min_for_y_z,         # minimum bounds in y-direction (bound_factor*2*R in meters)
                                             z_bound_min     = z_b_min_for_y_z,         # minimum bounds in z-direction (bound_factor*2*R in meters)
                                             x_bound_max     = x_b_max,                 # maximum bounds in x-direction (bound_factor*2*R in meters)
                                             y_bound_max     = y_b_max_for_y_z,         # maximum bounds in y-direction (bound_factor*2*R in meters)
                                             z_bound_max     = z_b_max_for_y_z,         # maximum bounds in z-direction (bound_factor*2*R in meters)
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
