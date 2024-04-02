#=##############################################################################
# DESCRIPTION
    This file performs a parameter study and loops over the number of blade elements.

# AUTHORSHIP
  * Author          : Fynn Gerhardy
  * Email           : fygerh@gmail.com
  * Created         : Jan 2024
  * Last updated    : Jan 2024
  * License         : -
=###############################################################################
#include(joinpath("/home/fynn/Repositories/BladeWakeOPT/02_parameter_study/ParameterSpace/ParameterStudy-Twist_Angle.jl"))

# ----------------- IMPORT PACKAGES --------------------------------------------------------
start_simulation_path = splitdir(@__FILE__)[1]
include(joinpath(start_simulation_path, "..", "..", "functions", "OwnFunctions.jl"))        # Read file with module definition first
using .OwnFunctions                                                                         # Include all self defined functions
import DataFrames


# ----------------- SETTINGS ---------------------------------------------------------------
# => FOLDER SETTINGS:
save_folder                 = "ParameterStudy-Twist_Angle"                       # Folder to store all simulations in
param_study_folderprefix    = "NREL5MW_TwistAngle"                                # Prefix of the folders each simulation will be stored in
save_path                   = joinpath(start_simulation_path, save_folder)          # Where to save this parameter study

# => ROTOR SPECIFICATIONS:
rotor_file      = "NREL5MW.csv"                                                     # Rotor geometry
data_path       = joinpath(start_simulation_path, "..", "..", "00_database")        # Path to rotor database

# read in the different twist distributions to be varied
twists = OwnFunctions.DataFrame(OwnFunctions.CSV.File(joinpath(start_simulation_path, "scaled_twists.csv")))

# => SIMULATION LENGTH SETTINGS:
nrevs           = 33                                                                # Number of revolutions to run
n               = 75                                                                # Number of blade elements per blade
nsteps_per_rev = 120                                                                # Number of steps per revolution


# => OPERATING CONDITIONS:
RPM             = 12.1                                                              # Rotational speed (1/min)
magVinf         = 11.4                                                                 # Free stream wind speed (m/s)
AOA             = 0.0                                                               # Angle of Attack (deg)

rho             = 1.225                                                             # Air density (kg/m^3)
mu              = 1.789e-5                                                           # Air dynamic viscosity (kg/ms)
speedofsound    = 342.35                                                            # Speed of sound (m/s)

# => POSTPROCESSING AND VISUALIZATION
postprocessing          = true                                                      # Perform postprocessing in general???
debug                   = true                                                      # Enables calculation of coefficients such as cn, ct, cl, cd
plot_bladeloads         = true                                                      # Postprocess the blade loads and plot the radial distribution (plots will be saved in "postprocessing folder")
show_bladeload_plots    = false                                                     # Show the bladeload plots on display after each simulation?
postprocess_fdom        = true                                                      # Postprocess the fluid domain and calculate velocity field of the last timestep?
paraview                = false                                                     # Whether to visualize with Paraview 
                                                                                    # (if true, the parameter study will stop after each turbine simulation...)
cylindrical_grid        = true                                                      # if true, the y-z plane will be calculated as a cylindrical grid and the wake velocity profiles will be saved within a .csv file
                                                                                    # this grid will be set automatically with the turbine diameter as its diameter

# => OTHER:
p_per_step      = 2                                                                 # Particles shed per step



# ----------------- START PARAMETERSTUDY ---------------------------------------------------

println("\n################################################################################")
println("START PARAMETERSTUDY $(param_study_folderprefix)")
println("################################################################################")
start_time_overall = time()

for i in 1:2:DataFrames.ncol(twists)
    start_time = time()
    twist_name = chop(DataFrames.names(twists)[i], tail=2)
    twist = Array{Float64, 2}(twists[:,i:i+1])

    println("\n--------------------------------------------------------------------------------")
    println("--------------------------------------------------------------------------------")
    println("--------------------------------------------------------------------------------")
    println("--------------------------------------------------------------------------------")
    println("EVALUATING "*param_study_folderprefix*"_$(twist_name) SIMULATION...")

    save_path_temp = joinpath(save_path, param_study_folderprefix*"_$(twist_name)")
    run_name = param_study_folderprefix*"_$(twist_name)"

    OwnFunctions.start_single_turbine_simulation(
                                    # ---- ESSENTIAL ARGUMENTS ---------
                                    start_simulation_path,
                                    save_path_temp,  
                                    run_name,
                                    rotor_file, 
                                    data_path,
                                    nrevs,
                                    nsteps_per_rev,
                                    n,
                                    RPM,
                                    magVinf,
                                    AOA;
                                    # ---- OPTIONAL ARGUMENTS ---------
                                    twist_overwrite =   twist,                            # scaled new twist distribution
                                    rho             =   rho,
                                    mu              =   mu,
                                    speedofsound    =   speedofsound,
                                    p_per_step      =   p_per_step,
                                    fidelity_extension      =   "mid",                     # options: "low", "mid", "high"
                                    x_loc           = 12,                                  # x location from wich on the particles will be cut away (x coordinate = x_loc*2*R in meters)
                                    # => POSTPROCESSING AND VISUALIZATION
                                    dt_fdom                 = 2,                           # every ...th timestep of the last revolution is postprocessed (fluiddomain)
                                    postprocessing          = postprocessing,
                                    debug                   = debug,
                                    plot_bladeloads         = plot_bladeloads,
                                    show_bladeload_plots    = show_bladeload_plots,
                                    postprocess_fdom        = postprocess_fdom,
                                    paraview                = paraview,
                                    cylindrical_grid        = cylindrical_grid,              # if true, the y-z plane will be calculated as a cylindrical grid and the wake velocity profiles will be saved within a .csv file
                                                                                            # this grid will be set automatically with the turbine diameter as its diameter
                                    cylinder_radius         = 2                          # set the cylinder radius (factor*R)
                                    )

    elapsed_time_s = time() - start_time
    # save time
    open(joinpath(save_path_temp, "00_elapsed_time_overall.txt"), "w") do file
      write(file, string(floor(elapsed_time_s)))
    end

    elapsed_time_min = floor(elapsed_time_s/60)
    elapsed_time_s = elapsed_time_s - (60 * elapsed_time_min)
    elapsed_time_h = floor(elapsed_time_min/60)
    elapsed_time_min = elapsed_time_min - (60 * elapsed_time_h)
    println("SIMULATION TOOK $(elapsed_time_h) HOURS, $(elapsed_time_min) MINUTES AND $(elapsed_time_s) SECONDS!")

end

elapsed_time_overall_s = time() - start_time_overall
elapsed_time_overall_min = floor(elapsed_time_overall_s/60)
elapsed_time_overall_s = elapsed_time_overall_s - (60 * elapsed_time_overall_min)
elapsed_time_overall_h = floor(elapsed_time_overall_min/60)
elapsed_time_overall_min = elapsed_time_overall_min - (60 * elapsed_time_overall_h)
println("\n################################################################################")
println("END PARAMETERSTUDY $(param_study_folderprefix)")
println("PARAMETER STUDY TOOK $(elapsed_time_overall_h) HOURS, $(elapsed_time_overall_min) MINUTES AND $(elapsed_time_overall_s) SECONDS!")
println("################################################################################")
