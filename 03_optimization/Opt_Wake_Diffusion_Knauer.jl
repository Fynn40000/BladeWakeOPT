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

LB_ventilation=R_lcaf
UB_ventilation=0.8
LB_SF=-1
UB_SF=1


# SAMPLING PLAN SETTINGS
method = "optimized" # sampling method (Options: "optimized", "random")
Nsamples = 15
Ndimensions = 2
lower_bound = [LB_ventilation, LB_SF]
upper_bound = [UB_ventilation, UB_SF]
dim_ranges = [(lower_bound[1], upper_bound[1]), (lower_bound[2], upper_bound[2])]
generations = 200 # used when method == "optimization"



#solution_points = OwnFunctions.get_solution_points(method, Nsamples, Ndimensions, dim_ranges, generations)
# solution_points = [0.20952334485714286 0.2857142857142858; 
#                   0.3079361207142857 0.7142857142857142; 
#                   0.7015872241428572 0.1428571428571428; 
#                   0.3571425086428571 -0.8571428571428572; 
#                   0.16031695692857142 -0.7142857142857143; 
#                   0.8 -0.2857142857142857; 
#                   0.5047616724285714 0.4285714285714286; 
#                   0.6031744482857143 -1.0; 
#                   0.4555552845 -0.5714285714285714; 
#                   0.5539680603571429 1.0; 
#                   0.4063488965714286 0.0; 
#                   0.7507936120714287 0.5714285714285714; 
#                   0.6523808362142858 -0.4285714285714286; 
#                   0.2587297327857143 -0.1428571428571429; 
#                   0.111110569 0.8571428571428572]

solution_points = [0.4063488965714286 0.1428571428571428;
                   0.3079361207142857 -0.2857142857142857;
                   0.3571425086428571 0.8571428571428572;
                   0.7015872241428572 1.0;
                   0.4555552845 -1.0;
                   0.5047616724285714 -0.4285714285714286; 
                   0.7507936120714287 0.4285714285714286; 
                   0.8 -0.5714285714285714; 
                   0.20952334485714286 -0.7142857142857143; 
                   0.111110569 -0.1428571428571429; 
                   0.16031695692857142 0.7142857142857142; 
                   0.6031744482857143 0.0; 
                   0.2587297327857143 0.2857142857142858; 
                   0.5539680603571429 0.5714285714285714; 
                   0.6523808362142858 -0.8571428571428572]

println(solution_points)
scatter(solution_points[:, 1], solution_points[:, 2], label="", xlabel="ventilation area", ylabel="scaling amount", legend=:topright)#, aspect_ratio=:equal



# no_lift_point = 