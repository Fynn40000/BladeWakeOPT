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

# turbine properties
R_hub = 0.024         # relative hub radius
R_lcaf = 0.111110569  # relative radius of last circular airfoil
R = 1.0               # relative tip radius



# SAMPLING PLAN SETTINGS
method = "optimized" # sampling method (Options: "optimized", "random")
Nsamples = 15
Ndimensions = 2
lower_bound = [0, -1]
upper_bound = [1, 1]
dim_ranges = [(lower_bound[1], upper_bound[1]), (lower_bound[2], upper_bound[2])]
generations = 200 # used when method == "optimization"



solution_points = OwnFunctions.get_solution_points(method, Nsamples, Ndimensions, dim_ranges, generations)
println(solution_points)
scatter(solution_points[:, 1], solution_points[:, 2], label="", xlabel="ventilation area", ylabel="scaling amount", legend=:topright)#, aspect_ratio=:equal




# no_lift_point = 