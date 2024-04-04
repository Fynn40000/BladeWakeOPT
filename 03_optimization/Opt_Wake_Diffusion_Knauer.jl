#=##############################################################################
# DESCRIPTION
    Peform optimization of wake diffusion rotor inspired by Knauer

# AUTHORSHIP
  * Author          : Fynn Gerhardy
  * Email           : fygerh@gmail.com
  * Created         : Mar 2024
  * Last updated    : Mar 2024
  * License         : -
=###############################################################################
#include(joinpath("/home/fynn/Repositories/BladeWakeOPT/03_optimization/Opt_Wake_Diffusion_Knauer.jl"))


# ----------------- IMPORT PACKAGES --------------------------------------------------------                                                                        # Include all self defined functions
using Surrogates
using Plots

start_simulation_path = splitdir(@__FILE__)[1]
include(joinpath(start_simulation_path, "..", "functions", "OwnFunctions.jl"))        # Read file with module definition first
using .OwnFunctions

# ----------------- TURBINE PROPERTIES --------------------------------------------------------   
R_hub = 0.024         # relative hub radius
R_lcaf = 0.111110569  # relative radius of last circular airfoil
R = 1.0               # relative tip radius

# ----------------- SEARCHSPACE BOUNDARIES -------------------------------------------------------- 
LB_ventilation=R_lcaf
UB_ventilation=0.8
LB_SF=-1
UB_SF=1

# ----------------- SAMPLING PLAN SETTINGS -------------------------------------------------------- 
use_freezed_points = true
lower_bound = [LB_ventilation, LB_SF]
upper_bound = [UB_ventilation, UB_SF]
  # !!! the following variables only apply when "use_freezed_points = false" !!!
method = "optimized"                                                                  # sampling method (Options: "optimized", "random")
Nsamples = 10
Ndimensions = 2
dim_ranges = [(lower_bound[1], upper_bound[1]), (lower_bound[2], upper_bound[2])]
generations = 2000                                                                     # used when method == "optimization"

# ----------------- SIMULATION SETTINGS -------------------------------------------------------- 
# => FOLDER SETTINGS:
save_folder                 = "solutions-Knauer_opt"                                  # Folder to store all simulations in
param_study_folderprefix    = "NREL5MW_solution"                                      # Prefix of the folders each simulation will be stored in
save_path                   = joinpath(start_simulation_path, save_folder)            # Where to save this parameter study

# => ROTOR SPECIFICATIONS:
rotor_file      = "NREL5MW.csv"                                                       # Rotor geometry
data_path       = joinpath(start_simulation_path, "..", "00_database")          # Path to rotor database

# => SIMULATION LENGTH SETTINGS:
nrevs           = 33                                                                  # Number of revolutions to run
n               = 75                                                                  # Number of blade elements per blade
nsteps_per_rev = 120                                                              # Number of steps per revolution

# => OPERATING CONDITIONS:
RPM             = 12.1                                                                # Rotational speed (1/min)
magVinf         = 11.4                                                                # Free stream wind speed (m/s)
AOA             = 0.0                                                                 # Angle of Attack (deg)

rho             = 1.225                                                               # Air density (kg/m^3)
mu              = 1.789e-5                                                            # Air dynamic viscosity (kg/ms)
speedofsound    = 342.35                                                              # Speed of sound (m/s)

# => POSTPROCESSING AND VISUALIZATION
postprocessing          = true                                                        # Perform postprocessing in general???
debug                   = true                                                        # Enables calculation of coefficients such as cn, ct, cl, cd
plot_bladeloads         = true                                                        # Postprocess the blade loads and plot the radial distribution (plots will be saved in "postprocessing folder")
show_bladeload_plots    = false                                                       # Show the bladeload plots on display after each simulation?
postprocess_fdom        = true                                                        # Postprocess the fluid domain and calculate velocity field of the last timestep?
paraview                = false                                                       # Whether to visualize with Paraview 
                                                                                      # (if true, the parameter study will stop after each turbine simulation...)
cylindrical_grid        = true                                                        # if true, the y-z plane will be calculated as a cylindrical grid and the wake velocity profiles will be saved within a .csv file
                                                                                      # this grid will be set automatically with the turbine diameter as its diameter

# => OTHER:
p_per_step      = 2                                                                   # Particles shed per step


# ----------------- GET THE SOLUTIONS TO BE CALCULATED -------------------------------------------------------- 
if !use_freezed_points
  solution_points = OwnFunctions.get_solution_points(method, Nsamples, Ndimensions, dim_ranges, generations)
else
  # solution_points = [0.34074037933333334 1.0; 
  #                    0.1876538391111111 -0.11111111111111116; 
  #                    0.41728364944444446 0.33333333333333326; 
  #                    0.5703701896666667 -1.0; 
  #                    0.6469134597777778 0.7777777777777777; 
  #                    0.111110569 0.5555555555555556; 
  #                    0.723456729888889 0.11111111111111116; 
  #                    0.4938269195555555 -0.33333333333333337; 
  #                    0.2641971092222222 -0.7777777777777778; 
  #                    0.8 -0.5555555555555556]
  solution_points = [0.32 -0.3;
                     0.39 -0.7;
                     0.47 0.9;
                     0.55 0.25;
                     0.62 -0.45]
                    
                    # Points 11-15
                    #  0.111110569 -1.0;
                    #  0.8 -1.0;
                    #  0.69 -0.9;
                    #  0.21 -0.55;
                    #  0.25 0.4;
                    # Points 16-20
                    #  0.32 -0.3;
                    #  0.39 -0.7;
                    #  0.47 0.9;
                    #  0.55 0.25;
                    #  0.62 -0.45

                     #0.111110569 1.0;   => End points represent "low lift rotor" (scaled to twist angle UB) => use 0.34074037933333334 1.0 vaules
                     #0.8 1.0;           => End points represent "low lift rotor" (scaled to twist angle UB) => use 0.34074037933333334 1.0 vaules
                     
end

# save the solution points in 2D array again
sol_points = []
for i in 1:size(solution_points, 1)
  push!(sol_points, [solution_points[i, 1], solution_points[i, 2]])
end

# println(solution_points)
# scatter(solution_points[:, 1], solution_points[:, 2], label="", xlabel="ventilation area", ylabel="scaling amount", legend=:topright)#, aspect_ratio=:equal


# ----------------- GET REFERENCE TWIST DISTRIBUTION TO BE SCALED -------------------------------------------------------- 
ref_aoa_min_df = OwnFunctions.DataFrame(OwnFunctions.CSV.File(joinpath(start_simulation_path, "aoa_min_bound.csv")))
ref_aoa_min_arr = Array{Float64, 2}(ref_aoa_min_df)
ref_aoa_max_df = OwnFunctions.DataFrame(OwnFunctions.CSV.File(joinpath(start_simulation_path, "aoa_max_bound.csv")))
ref_aoa_max_arr = Array{Float64, 2}(ref_aoa_max_df)
ref_aoa_df = OwnFunctions.DataFrame(OwnFunctions.CSV.File(joinpath(start_simulation_path, "ref_aoa.csv")))
ref_aoa_arr = Array{Float64, 2}(ref_aoa_df)
ref_beta_df = OwnFunctions.DataFrame(OwnFunctions.CSV.File(joinpath(start_simulation_path, "ref_beta.csv")))
ref_beta_arr = Array{Float64, 2}(ref_beta_df)

ref_x = ref_aoa_arr[:,1]
ref_aoa_min = ref_aoa_min_arr[:,2]
ref_aoa_max = ref_aoa_max_arr[:,2]
ref_aoa = ref_aoa_arr[:,2]
ref_beta = ref_beta_arr[:,2]

    # adjust aoa max and aoa min if neccessary
ref_aoa_min_new = deepcopy(ref_aoa_min)
ref_aoa_max_new = deepcopy(ref_aoa_max)
Offset = zeros(length(ref_aoa))
for i in 1:length(ref_aoa)
    if ref_aoa_max[i] < (-1 * ref_aoa[i])
      ref_aoa_max_new[i] = -1 * ref_aoa[i]
    end
    
    Offset[i] = ref_aoa_max_new[i] - ref_aoa_max[i]
    ref_aoa_min_new[i] = ref_aoa_min[i] + Offset[i]
end

ref_phi_min = (ref_beta .* -1) .- ref_aoa_max_new
ref_phi_max = (ref_beta .* -1) .- ref_aoa_min
ref_phi_max[1] = ref_phi_max[3] 
ref_phi_max[2] = ref_phi_max[3]
ref_phi = (ref_beta .* -1) .- (ref_aoa .* -1)


# ----------------- LOOP THROUGH DIFFERENT TWIST DISTRIBUTIONS -------------------------------------------------------- 
# 1. calculate controlpoints and bezier curves
# 2. create scaling curve by appending the bezier curves to each other
# 3. calculate scaled twist distribution
# 4. start simulation and save data using the scaled twist distribution

println("\n################################################################################")
println("START OPTIMIZATION $(param_study_folderprefix)")
println("################################################################################")
start_time_overall = time()
# LOOP - Simulations
for (sol_idx, solution) in enumerate(sol_points)

  start_time = time()
  r_ventilation = solution[1]
  Psi_at_R = solution[2]
  
  #name = "_rvent$(round(r_ventilation, digits=3))_Psi$(round(Psi_at_R, digits=3))"
  name = "_$(sol_idx)"

  # 1. calculate controlpoints of this solution to create the bezier curves
  if r_ventilation <= (R_lcaf + ((R - R_lcaf) / 2))
    chi11_x = (r_ventilation - ((r_ventilation - R_lcaf)/2))
    chi12_x = (r_ventilation + ((r_ventilation - R_lcaf)/2))
    chi22_x = (chi12_x + ((R - UB_ventilation) / 2))

  elseif r_ventilation > (R_lcaf + ((R - R_lcaf) / 2))
    chi11_x = (r_ventilation - ((R - r_ventilation)/2))
    chi12_x = (r_ventilation + ((R - r_ventilation)/2))
    chi22_x = (chi12_x + ((R - UB_ventilation) / 2))

  end

  CPs_B1 = [[R_lcaf, 1.0],                              # chi01
            [chi11_x, 1.0],                             # chi11
            [r_ventilation, 1.0]]                       # chi21

  CPs_B2 = [[r_ventilation, 1.0],                       # chi02
            [chi12_x, 1.0],                             # chi12
            [chi22_x, Psi_at_R],                        # chi22
            [R, Psi_at_R]]                              # chi32

  ################################################################
  # CPs_B1 = [[R_lcaf, 1.0],                              # chi01
  #           [0.23022195, 1.0],                             # chi11
  #           [0.34933333, 1.0]]                       # chi21

  # CPs_B2 = [[0.34933333, 1.0],                       # chi02
  #           [0.46844472, 1.0],                             # chi12
  #           [0.46844472+0.1, -1.0],                        # chi22
  #           [R, -1.0]]                              # chi32
  ################################################################

  # calculate the bezier curves
  t_values = collect(range(0, stop=1, length=100))
  SC1 = [OwnFunctions.bezier_point(t, CPs_B1) for t in t_values]  # first bezier curve
  SC2 = [OwnFunctions.bezier_point(t, CPs_B2) for t in t_values]  # second bezier curve
  SC0 = [[R_hub, SC1[1][2]], [R_lcaf, SC1[1][2]]]                 # part of scaling curve within hub area

  # 2. create scaling curve by appending the bezier curves to each other
  SC = vcat([SC0[1]], SC1[1:end-1], SC2)                          # scaling curve
  SC_x = [point[1] for point in SC]                               # x-values of scaling curve (from x=R_hub to x=R)
  SC_y = [point[2] for point in SC]                               # y-values of scaling curve (scaling factors)
  
  # 3. calculate scaled twist distribution
  twist_scaled = OwnFunctions.get_scaled_twist(ref_x, SC_x, SC_y, ref_phi, ref_phi_min, ref_phi_max)
  col1 = reshape(ref_x, :, 1)
  col2 = reshape(twist_scaled, :, 1)
  twist = hcat(Float64.(col1), Float64.(col2))

  # 4. start simulation and save data using the scaled twist distribution
  println("\n--------------------------------------------------------------------------------")
  println("--------------------------------------------------------------------------------")
  println("--------------------------------------------------------------------------------")
  println("--------------------------------------------------------------------------------")
  println("EVALUATING " * param_study_folderprefix * name *" SIMULATION...")

  save_path_temp = joinpath(save_path, param_study_folderprefix * name)
  run_name = param_study_folderprefix * name

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
  # save parameter settings
  open(joinpath(save_path_temp, "00_solution_parameters.txt"), "w") do file
    write(file, string(r_ventilation))
    write(file, "\n" * string(Psi_at_R))
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
println("END KNAUER OPTIMIZATION")
println("PARAMETER STUDY TOOK $(elapsed_time_overall_h) HOURS, $(elapsed_time_overall_min) MINUTES AND $(elapsed_time_overall_s) SECONDS!")
println("################################################################################")

#plot(SC[:, 1], SC[:, 2], label="", xlabel="r/R", ylabel="Psi", legend=:topright)#, aspect_ratio=:equal
#plot(x, y, label="", xlabel="r/R", ylabel="Psi", legend=:topright)#, aspect_ratio=:equal